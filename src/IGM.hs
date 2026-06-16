{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

module IGM where

{-
Module      : HASKy.IGM
Description : Inter-Galactic Medium evolution
Copyright   : (c) Oleksii Sokoliuk, 20256
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

This module is the most important part of the HASKy package. It
defines two evolutionary equations for the IGM and ISM and enables us to
observe the chemical evolution of various metals in the IGM/ISM at various
redshifts.
-}

import Control.Lens.Combinators (each, over)
import Control.Parallel.Strategies
import Cosmology
import Data.Bifunctor
import Data.List (transpose)
import qualified Data.Map as M
import qualified Data.Vector as V
import HMF
import Helper
import Lookup
import Math.GaussianQuadratureIntegration
import SMF
import StarFormation

data IMFKind
  = Salpeter
  | Kroupa
  | Chabrier
  | SN_Ia
  deriving (Eq, Show, Read)

data RemnantKind
  = Pereira
  | WoosleyWeaver
  deriving (Eq, Show, Read)

data HNeKind
  = Constant_HNe
  | Grimmett
  deriving (Eq, Show, Read)

data PNSKind
  = Arcones
  deriving (Eq, Show, Read)

data EjectaOutflow a
  = MkEjectaOutflow
  { e_AGB :: a,
    e_AGB_Element :: a,
    e_CCSN :: a,
    e_CCSN_Element :: a,
    e_HNe_Element :: a,
    e_MRSNe_Element :: a,
    e_SNe_Ia :: a,
    e_SNe_Ia_Element :: a,
    e_Novae_Element :: a,
    e_NSM :: a,
    e_NSM_Element :: a,
    o_Wind :: a
  }
  deriving (Eq, Show, Read, Functor)

-- | Models for the Initial Mass Function:
--  * Single power-law [Salpeter et al. 1995]
--  * Broken power-law [Kroupa et al. 2001]
--  * Log-normal [Chabrier et al. 2003] converted to [Mpc^-3]
initialMassFunction :: IMFKind -> Mstar -> Double
initialMassFunction iKind m =
  let (alpha0, alpha1, alpha2, m1, m2, k0, k1, k2) =
        (-0.3, -1.3, -2.3, 0.08, 0.5, 0.58, k0 * m1 ** (alpha0 - alpha1), k1 * m2 ** (alpha1 - alpha2))
      (a_Ch, b_Ch, center_Ch, sigma_Ch) =
        (0.85, 0.24, 0.079, 0.69)
   in case iKind of
        Salpeter -> m ** (-2.35)
        Kroupa
          | m < m1 -> k0 * m ** alpha0
          | m >= m1 && m < m2 -> k1 * m ** alpha1
          | otherwise -> k2 * m ** alpha2
        Chabrier
          | m < 1 -> a_Ch * exp (-(log m - log center_Ch) ** 2 / (2 * sigma_Ch ** 2))
          | otherwise -> b_Ch * m ** (-1.3)
        SN_Ia -> m ** (-1.35)

-- | Next two functions normalise IMF between m_inf = 0.1 Msol and m_sup = 100 Msol
-- so that it can acts as a PDF in that range
imfNormalisation :: ReferenceCosmology -> IMFKind -> Mstar -> Mstar -> Double
imfNormalisation MkCosmology {prec} iKind m_down m_up =
  let m_arr = [m_down, m_down + 0.1 .. m_up]
      integrand m = m * initialMassFunction iKind m
      integrator = makeIntegrator prec
   in integrator integrand m_down m_up

normalisedInitialMassFunction :: ReferenceCosmology -> IMFKind -> Mstar -> Mstar -> Mstar -> Double
normalisedInitialMassFunction cosmology iKind m_down m_up m =
  let norm = imfNormalisation cosmology iKind m_down m_up
   in (initialMassFunction iKind m) / norm

-- Define a normalised IMF for a specific IMF kind and Cosmology kind
normImf :: ReferenceCosmology -> IMFKind -> Double -> Double
normImf cosmology iKind =
  let IGMParams {..} = defaultIGMParams
      MassLimits {..} = masses
   in normalisedInitialMassFunction cosmology iKind 0.01 mUp

-- More so for the sake of simplicity, define a separate normalised SN Ia IMF
normImfSN :: ReferenceCosmology -> Double -> Double -> Double -> Double
normImfSN cosmology =
  normalisedInitialMassFunction
    cosmology
    SN_Ia

-- | Mass of a remnant produced by the supernova,
-- calculated according to various works in the units of [Msol].
-- Note that below a certain mass, stellar lifetimes are higher than the
-- t_0 (age of the universe), so there will be no remnant
massRemnant :: Mstar -> RemnantKind -> Double -> Double
massRemnant m rKind metalFraction =
  case rKind of
    -- Remnant masses taken from the [Woosley & Weaver 95]
    -- and [Hoek & Groenewegen 1996]
    WoosleyWeaver
      | m < 0.9 -> 0
      | m <= 8 ->
          let (mi, mr) = (unzip . M.toList) (remnantMediumMass metalFraction)
           in makeInterp mi mr $ m
      | m <= 40 ->
          let (mi, mr) = (unzip . M.toList) (remnantHighMass metalFraction)
           in makeInterp mi mr $ m
      | otherwise ->
          extrapolate
            m
            ( (\x y -> [x, y])
                <$> M.findWithDefault 0.0 35
                <*> M.findWithDefault 0.0 40
                $ remnantHighMass
                  metalFraction
            )
            [35, 40]
    -- Remnant masses from the [Pereira & Miranda 2010]
    -- Note that unlike WW95 work, these are metallicity-independent
    Pereira
      | m < 1 -> 0
      | m <= 8 -> 0.1156 * m + 0.4551
      | m <= 10 -> 1.35
      | m < 25 -> 1.4
      | m < 140 -> 13 / 24 * (m - 20)
      | m <= 260 -> 0
      | m <= 500 -> m

-- | Lifetime of a main sequence star in relation to it's mass,
-- taken from the work of [Maeder & Meynet 1989] and the extrapolation to
-- m > 60 Msol is taken from the [Romano et al. 2005]
tauMS :: Mstar -> CosmicTime
tauMS m
  | m <= 1.3 = 10 ** (-0.6545 * log10 m + 1)
  | m <= 3 = 10 ** (-3.7 * log10 m + 1.35)
  | m <= 7 = 10 ** (-2.51 * log10 m + 0.77)
  | m <= 15 = 10 ** (-1.78 * log10 m + 0.17)
  | m <= 60 = 10 ** (-0.86 * log10 m - 0.94)
  | otherwise = 1.2 * m ** (-1.85) + 0.003

-- | Mass of the main sequence star that dies at age t [Gyr].
-- This is practically an inverse of tauMS function defined above
mDynamicalTime :: CosmicTime -> Mstar -> Mstar
mDynamicalTime t m_up =
  let mass_range = [0.1, 0.1 + 0.5 .. m_up]
      interp_tau = makeInterp (tauMS <$> mass_range) mass_range
   in interp_tau t

-- | Dynamical mass of a star at redshift z
mDynamicalRedshift :: ReferenceCosmology -> Double -> Double
mDynamicalRedshift cosmology z =
  let IGMParams {..} = defaultIGMParams
      MassLimits {..} = masses
   in mDynamicalTime (interpT cosmology z) mUp

-- | Lower limit for a mass of a CCSN
mDown :: ReferenceCosmology -> Double -> Double
mDown cosmology z = maximum [8, mDynamicalRedshift cosmology z]

-- | Progenitor mass - proto NS mass relation from []
massPNS :: PNSKind -> Mstar -> Mstar
massPNS kind_PNS =
  case kind_PNS of
    Arcones ->
      let masses_proto = [1.4, 1.6, 1.8, 2.0]
          masses_prog = [13, 15, 20, 40]
       in makeInterp masses_prog masses_proto

-- | Fraction of CCSNe that are HNe for higher masses (fiducially, for M >= 20 Msol)
-- There are two models:
--  * Constant fraction, a toy model
--  * More complicated, metallicity dependent [Grimmett et al. 2020] model
fractionHNe :: HNeKind -> Double -> Double -> Double
fractionHNe hneKind epsHNe0 metalFraction =
  case hneKind of
    Constant_HNe -> epsHNe0
    Grimmett -> maximum (epsHNe0 * exp (-metalFraction / 0.001), 0.001)

-- | Construct yield of SNe type II by interpolating over AGB, SAGB, ECSNe and CCSNe yields
-- If neutrino-driven yields are turned on, add yields to SN II yields following the
-- initial mass - NS mass relation
yieldII :: ReferenceStarFormationModel -> Double -> Double
yieldII yields m =
  y1 m + y2 m
  where
    (mSAGB, ySAGB) = yield_sagb yields
    yECSN = yield_ecsn yields
    (mII, yII)
      | m < 8.0 = yield_agb yields
      | m >= 8.8 && m < 9.0 && yECSN /= 0 = ([8.8, 9.0], [yECSN, yECSN])
      | m < 11.0 && null ySAGB = yield_ccsn yields
      | m < 11.0 = (mSAGB, ySAGB)
      | otherwise = yield_ccsn yields
    y1 = uncurry makeInterp $ yield_psn yields
    y2 = makeInterp mII yII

-- | Calculate total mass ejected by an AGB star from the total remnant mass
massEjectedAGB :: RemnantKind -> Double -> Double -> Double
massEjectedAGB rKind metalFraction m = m - massRemnant m rKind metalFraction

-- | Calculate total mass of a particular isotope ejected by an AGB star
massEjectedAGB_Element :: ReferenceStarFormationModel -> RemnantKind -> Double -> Double -> Double
massEjectedAGB_Element yields rKind metalFraction m =
  massEjectedAGB rKind metalFraction m - yieldII yields m

-- | Calculate metal fraction at redshift z
metalFractionAtZ :: ReferenceCosmology -> Double -> Metallicity -> Double
metalFractionAtZ cosmology z metalFraction = metalFraction . interpT cosmology $ z

-- | Calculate a yield from Novae star
yield_Novae :: ReferenceStarFormationModel -> Double
yield_Novae yields =
  let IGMParams {..} = defaultIGMParams
      Efficiencies {..} = effs
   in epsCO * (yield_co yields) + (1 - epsCO) * (yield_one yields)

-- | Calculate a yield from an HNe star as a function of mass
yield_HNe :: ReferenceStarFormationModel -> Double -> Double
yield_HNe yields = yieldInterp
  where
    (mHNe, yHNe) = yield_hne yields
    yieldInterp = makeInterp mHNe yHNe

-- | Helper function to construct an IMF normalised integral
mkIntegrand ::
  ReferenceCosmology -> (Double -> Double) -> (Double -> Double) -> Double -> (Double -> Double) -> Double -> Double -> Double
mkIntegrand cosmology norm sfrd offset yield z m =
  norm m
    * sfrd (interpT cosmology z - offset)
    * yield m

-- | Similarly, a helper function to construct an IMF normalised integral with normalised IMF and no time delay
integrand0 :: ReferenceCosmology -> IMFKind -> (Double -> Double) -> (Double -> Double) -> Double -> Double -> Double
integrand0 cosmology iKind sfrd yield = mkIntegrand cosmology (normImf cosmology iKind) sfrd 0 yield

-- | Some of the terms (IGM/ISM outflows for all mass and a specific element yield, SFRD),
-- to be used in the next function
interGalacticMediumTerms ::
  ReferenceCosmology ->
  ReferenceStarFormationModel ->
  PowerSpectrum ->
  RemnantKind ->
  IMFKind ->
  SMFKind ->
  HMFKind ->
  WKind ->
  HNeKind ->
  Metallicity ->
  Mhalo ->
  SFRD ->
  Redshift ->
  EjectaOutflow Double
interGalacticMediumTerms cosmology@MkCosmology {prec} yields pk rKind iKind sKind hKind wKind hneKind metalFraction mh_min sfrd z =
  let -- Unwrap all of the user-defined parameters
      IGMParams {..} = defaultIGMParams
      PhysicalConstants {..} = phys
      MassLimits {..} = masses
      Efficiencies {..} = effs
      DelayTimes {..} = delays
      NovaeParams {..} = novae
      ECSNParams {..} = ecsn

      -- Integrands
      integrand_AGB m = integrand0 cosmology iKind sfrd (massEjectedAGB rKind (metalFractionAtZ cosmology z metalFraction)) z m
      integrand_AGB_Element m = integrand0 cosmology iKind sfrd (massEjectedAGB_Element yields rKind (metalFractionAtZ cosmology z metalFraction)) z m
      integrand_ECSN_Element m = mkIntegrand cosmology (normImf cosmology iKind) sfrd (tauMS m) (yieldII yields) z m
      integrand_Wind = mkIntegrand cosmology (normImf cosmology iKind) sfrd 0 (\m -> 2 * snEnergy / (kmsErgMsol * escapeVelocitySq cosmology pk hKind wKind mh_min z)) z
      integrand_CCSN m = (1 - fractionHNe hneKind epsHNe0 (metalFractionAtZ cosmology z metalFraction) - epsMRSNe) * mkIntegrand cosmology (normImf cosmology iKind) sfrd (tauMS m) (const 1) z m
      integrand_CCSN_Element m = (1 - fractionHNe hneKind epsHNe0 (metalFractionAtZ cosmology z metalFraction) - epsMRSNe) * mkIntegrand cosmology (normImf cosmology iKind) sfrd (tauMS m) (yieldII yields) z m
      integrand_NSM = mkIntegrand cosmology (normImf cosmology iKind) sfrd delayNSM (const 1) z
      integrand_Novae_Element m = mEjNovae * nNovae * alphaNovae * yield_Novae yields * mkIntegrand cosmology (normImf cosmology iKind) sfrd delayNovae (const 1) z m
      integrand_SNe_Ia = normImf cosmology iKind
      integrand_HNe_Element m = (fractionHNe hneKind epsHNe0 (metalFractionAtZ cosmology z metalFraction) - epsMRSNe) * mkIntegrand cosmology (normImf cosmology iKind) sfrd (tauMS m) (yield_HNe yields) z m
      integrand_MRSNe_Element m = epsMRSNe * 1
      integrand_WR m = 1
      integrand_PISNe m = 1

      -- Numerical integration
      integrator = makeIntegrator prec

      e_AGB = integrator integrand_AGB mAGBd mAGBu
      e_AGB_Element = integrator integrand_AGB_Element mAGBd mAGBu
      o_Wind = epsW * integrator integrand_Wind (mDown cosmology z) mUp
      e_CCSN = integrator integrand_CCSN (mDown cosmology z) mUp
      e_CCSN_Element = integrator integrand_CCSN_Element (mDown cosmology z) mUp
      e_ECSN_Element = integrator integrand_ECSN_Element mECSNd mECSNu
      e_HNe_Element = integrator integrand_HNe_Element (mDown cosmology z) mUp
      e_MRSNe_Element = epsMRSNe * e_MRSNe_Element
      e_Novae_Element = integrator integrand_Novae_Element mNovaeD mNovaeU
      e_NSM = alphaNSM * integrator integrand_NSM mNSMd mNSMu
      e_NSM_Element = yield_nsm yields * e_NSM

      first_term =
        bRG
          * integrator (\m -> normImfSN cosmology mDLRG mDURG m * sfrd (tauMS m)) mDLRG mDURG
      second_term =
        bMS
          * integrator (\m -> normImfSN cosmology mDLMS mDUMS m * sfrd (tauMS m)) mDLMS mDUMS
      e_SNe_Ia =
        mCO
          * integrator integrand_SNe_Ia (maximum [mPL, mDynamicalRedshift cosmology z]) mPU
          * (first_term + second_term)
      e_SNe_Ia_Element = yield_ia yields * e_SNe_Ia
   in fmap
        (* yrGyr)
        MkEjectaOutflow
          { e_AGB,
            e_AGB_Element,
            e_CCSN,
            e_CCSN_Element,
            e_HNe_Element,
            e_MRSNe_Element,
            e_SNe_Ia,
            e_SNe_Ia_Element,
            e_Novae_Element,
            e_NSM,
            e_NSM_Element,
            o_Wind
          }

-- | Extract rates from the EjectaOutflow dataclass and convert them into an appropriate form
-- to be used further in the ODE solver
computeRates :: EjectaOutflow Double -> (Double, Double, Double, Double, Double, Double)
computeRates rates =
  ( e_AGB + e_SNe_Ia + e_NSM + e_CCSN,
    e_AGB_Element + e_CCSN_Element + e_SNe_Ia_Element + e_NSM_Element,
    eps_sn * e_AGB,
    o_Wind,
    eps_sn * e_AGB_Element,
    eps_sn * e_AGB + o_Wind
  )
  where
    eps_sn = 0.005
    MkEjectaOutflow {e_AGB, e_AGB_Element, e_CCSN, e_CCSN_Element, e_SNe_Ia, e_SNe_Ia_Element, e_NSM, e_NSM_Element, o_Wind} = rates

-- Separate IO function for IGM/ISM terms
igmTermsIO ::
  ReferenceStarFormationConfig ->
  ReferenceCosmology ->
  PowerSpectrum ->
  RemnantKind ->
  IMFKind ->
  SMFKind ->
  HMFKind ->
  WKind ->
  HNeKind ->
  Metallicity ->
  Mhalo ->
  SFRD ->
  Redshift ->
  Element ->
  IO (Double, Double, Double, Double, Double, Double)
igmTermsIO sfCfg cosmology pk rKind iKind sKind hKind wKind hneKind metalFraction mh_min sfrd z elem =
  do
    snIaYieldRetrieved <- retrieveYieldIa sfCfg elem
    ccsnYieldRetrieved <- retrieveYieldCCSN sfCfg elem (metalFraction z)
    agbYieldRetrieved <- retrieveYieldAGB sfCfg elem (metalFraction z)
    ecsnYieldRetrieved <- retrieveYieldECSN sfCfg elem
    hneYieldRetrieved <- retrieveYieldHNe sfCfg elem (metalFraction z)
    let yields =
          MkStarFormation
            { yield_ia = snIaYieldRetrieved,
              yield_ccsn = ccsnYieldRetrieved,
              yield_hne = hneYieldRetrieved,
              yield_ecsn = ecsnYieldRetrieved,
              yield_agb = agbYieldRetrieved,
              yield_sagb = ([0], [0]),
              yield_nsm = 0,
              yield_psn = ([0], [0]),
              yield_co = 0,
              yield_one = 0
            }
        result = computeRates $ interGalacticMediumTerms cosmology yields pk rKind iKind sKind hKind wKind hneKind metalFraction mh_min sfrd z
    return result

-- Separate IO function for SN terms
snTermsIO ::
  ReferenceStarFormationConfig ->
  ReferenceCosmology ->
  PowerSpectrum ->
  RemnantKind ->
  IMFKind ->
  SMFKind ->
  HMFKind ->
  WKind ->
  HNeKind ->
  Metallicity ->
  Mhalo ->
  SFRD ->
  Redshift ->
  Element ->
  IO (Double, Double, Double, Double)
snTermsIO sfCfg cosmology pk rKind iKind sKind hKind wKind hneKind metalFraction mh_min sfrd z elem =
  do
    snIaYieldRetrieved <- retrieveYieldIa sfCfg elem
    ccsnYieldRetrieved <- retrieveYieldCCSN sfCfg elem (metalFraction z)
    agbYieldRetrieved <- retrieveYieldAGB sfCfg elem (metalFraction z)
    ecsnYieldRetrieved <- retrieveYieldECSN sfCfg elem
    hneYieldRetrieved <- retrieveYieldHNe sfCfg elem (metalFraction z)
    let yields =
          MkStarFormation
            { yield_ia = snIaYieldRetrieved,
              yield_ccsn = ccsnYieldRetrieved,
              yield_hne = hneYieldRetrieved,
              yield_ecsn = ecsnYieldRetrieved,
              yield_agb = agbYieldRetrieved,
              yield_sagb = ([0], [0]),
              yield_nsm = 0,
              yield_psn = ([0], [0]),
              yield_co = 0,
              yield_one = 0
            }
        result =
          let MkEjectaOutflow {e_AGB, e_AGB_Element, e_CCSN, e_CCSN_Element, e_SNe_Ia, e_SNe_Ia_Element, e_NSM, e_NSM_Element, o_Wind} =
                interGalacticMediumTerms cosmology yields pk rKind iKind sKind hKind wKind hneKind metalFraction mh_min sfrd z
           in (e_CCSN, e_CCSN_Element, e_SNe_Ia, e_SNe_Ia_Element)

    return result

-- | Set up initial abundance for a given isotope
-- Big Bang Nucleosynthesis abundances are taken from the [Coc et al. 2014]
iniAbundance :: Element -> Double
iniAbundance elem =
  case element elem of
    "H" -> if isotope elem == 1 then (1 - fd) * fh else fd * fh
    "He" -> if isotope elem == 4 then (1 - fh3) * yp else fh3
    "Li" | isotope elem == 7 -> fli7
    _ -> 0
  where
    yp = 0.2464
    fh = 1 - yp
    fd = 2.64e-5
    fh3 = 1.05e-5 * (1 - fd) * fh
    fli7 = 5.18e-10 * (1 - fd) * fh

-- | Solve four coupled first-order differential equations that govern the evolution of:
--    * rho_IGM (1)
--    * rho_ISM (3)
--    * Xi_IGM  (4)
--    * Xi_ISM  (5)
-- with all equations being taken from the [Daigne et al. 2004]
{-# INLINE igmIsmEvolution #-}
igmIsmEvolution ::
  ReferenceStarFormationConfig ->
  ReferenceCosmology ->
  PowerSpectrum ->
  RemnantKind ->
  IMFKind ->
  SMFKind ->
  HMFKind ->
  WKind ->
  HNeKind ->
  Element ->
  Mhalo ->
  IO ([Double], [V.Vector Double])
igmIsmEvolution sfCfg cosmology@MkCosmology {h0, om0, ob0, gn} pk rKind iKind sKind hKind wKind hneKind elem mh_min =
  do
    let -- Constructing stellar yields datatype as a function of metallicity of ISM and given isotope
        zs = [20.0, 20.0 - 0.5 .. 0]

        mar =
          parMap rpar (\z -> 1e9 * baryonFormationRateDensity cosmology pk hKind wKind z) zs
        ts =
          parMap rpar (\z -> cosmicTime cosmology z) zs
        sfrd =
          parMap rpar (\z -> starFormationRateDensity cosmology pk sKind hKind wKind z mh_min) zs
        rhoCr z = 3 * h0 ** 2 / (8 * pi * gn) * (1 + z) ** 3

        -- Unpack outflow/inflow rates and interpolate over our redshift range
        (interpMAR, interpSFRD, interpZ, interpT) = (makeInterp zs mar, makeInterp ts sfrd, makeInterp ts zs, makeInterp zs ts)

        -- ICs are set assuming very small baryon fraction in the structures,
        -- with rho_ISM/rho_IGM ~ 0.01 (stellar mass is negligible at this redshift)
        -- following the prescription of [Daigne et al. 2006]
        -- Finally, we also adopt the BBN abundances for H (He),
        -- such that the ICs for Xi_ISM/Xi_IGM = 0.76 (0.24) * M_ISM/M_IGM.
        (nSteps, tInit, aInit, rhoTot, igmInit, ismInit, xiIgmInit, xiIsmInit, ejectaInit, outflowInit) =
          (10 :: Int, interpT (maximum zs), 0.01, rhoCr (maximum zs) * ob0, (1 - aInit) * rhoTot, aInit * rhoTot, iniAbundance elem, xiIgmInit, 0, 0)

        -- Convert Differential-Algebraic system into a system of ODEs
        odeSystem :: History -> Double -> V.Vector Double -> IO (V.Vector Double)
        odeSystem history t y = do
          let (times, metals) =
                unzip $
                  [(t, metal) | (t, v) <- history, V.length v > 2, let metal = (v V.! 4) / v V.! 1]

              interpMetal t =
                if length times < 1
                  then 0
                  else
                    (makeInterp times metals) t

              -- Unpack all rates at z(t)
              z = interpZ t

          (e_tot, e_tot_Element, o_SNe, o_Wind, o_SNe_Element, o_tot) <- igmTermsIO sfCfg cosmology pk rKind iKind sKind hKind wKind hneKind interpMetal mh_min interpSFRD z elem

          pure $
            V.fromList
              [ -interpMAR z + o_tot,
                (-1e9 * interpSFRD t + e_tot) + (interpMAR z - o_tot),
                1 / (y V.! 0) * (o_Wind * (y V.! 3 - y V.! 2) + (o_SNe_Element - o_SNe * y V.! 2)),
                1 / (y V.! 1) * ((e_tot_Element - e_tot * y V.! 3) + interpMAR z * (y V.! 2 - y V.! 3) - (o_SNe_Element - o_SNe * y V.! 3)),
                e_tot,
                o_tot
              ]

    -- Solve the system and unpack values
    zippedHistory <-
      rk4SolveHistIO odeSystem tInit ((interpT 0 - tInit) / fromIntegral nSteps) nSteps (V.fromList [igmInit, ismInit, xiIgmInit, xiIsmInit, ejectaInit, outflowInit])

    -- M_star = rho_tot - M_IGM - Mstar from the conservation equation
    -- In addition, we also normalise each mass by the total mass
    let (times, masses) = unzip zippedHistory
        result =
          ( \v ->
              V.zipWith ($) (V.fromList [(/ rhoTot), (/ rhoTot), (* 1), (* 1), (/ rhoTot), (/ rhoTot), (/ rhoTot)]) $ V.snoc v (rhoTot - v V.! 0 - v V.! 1)
          )
            <$> masses
    return (times, result)

{-
-- | Derive the metallicity of the IGM/ISM from outflow/inflow rates of metals,
-- currently using an approach presented in [Tan et al. 2018]
{-# INLINE igmIsmMetallicity #-}
igmIsmMetallicity ::
  ReferenceCosmology ->
  PowerSpectrum ->
  RemnantKind ->
  IMFKind ->
  SMFKind ->
  HMFKind ->
  WKind ->
  Mhalo ->
  ([Double], [Double], [Double])
igmIsmMetallicity cosmology pk rKind iKind sKind hKind wKind yield_ii yield_ia yield_nsm mh_min =
  let elem = Element "H" 1

      (times, densities) =
        bimap id (transpose . (fmap V.toList)) $
          igmIsmEvolution cosmology pk rKind iKind sKind hKind wKind yield_ii yield_ia yield_nsm elem mh_min

      result_igm = zipWith (/) (densities !! 5) (densities !! 0)
      result_ism = zipWith (/) (densities !! 4) (densities !! 1)
   in (times, result_igm, result_ism)

-}
