{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}

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

data IMF_kind
  = Salpeter
  | Kroupa
  | Chabrier
  | SN_Ia
  deriving (Eq, Show, Read)

data Remnant_Kind
  = Pereira
  | WoosleyWeaver
  deriving (Eq, Show, Read)

data HNe_Kind
  = Constant_HNe
  | Grimmett
  deriving (Eq, Show, Read)

data PNS_Kind
  = Arcones
  deriving (Eq, Show, Read)

data EjectaOutflow a
  = MkEjectaOutflow
  { e_AGB :: a,
    e_AGB_Element :: a,
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
initialMassFunction :: IMF_kind -> Mstar -> Double
initialMassFunction i_kind m =
  let (alpha0, alpha1, alpha2, m1, m2, k0, k1, k2) =
        (-0.3, -1.3, -2.3, 1, 0.08, 0.5, k0 * m1 ** (alpha0 - alpha1), k1 * m2 ** (alpha1 - alpha2))
      (a_Ch, b_Ch, center_Ch, sigma_Ch) =
        (0.85, 0.24, 0.079, 0.69)
   in case i_kind of
        Salpeter -> m ** (-2.35)
        Kroupa
          | m < m1 -> k0 * m ** alpha0
          | m >= m2 && m < m2 -> k1 * m ** alpha1
          | otherwise -> k2 * m ** alpha2
        Chabrier
          | m < 1 -> a_Ch * exp (-(log m - log center_Ch) ** 2 / (2 * sigma_Ch ** 2))
          | otherwise -> b_Ch * m ** (-1.3)
        SN_Ia -> m ** (-1.35)

-- | Next two functions normalise IMF between m_inf = 0.1 Msol and m_sup = 100 Msol
-- so that it can acts as a PDF in that range
imfNormalisation :: ReferenceCosmology -> IMF_kind -> Mstar -> Mstar -> Double
imfNormalisation MkCosmology {prec} i_kind m_down m_up =
  let m_arr = [m_down, m_down + 0.1 .. m_up]
      integrand m = m * initialMassFunction i_kind m
      integrator = makeIntegrator prec
   in integrator integrand m_down m_up

normalisedInitialMassFunction :: ReferenceCosmology -> IMF_kind -> Mstar -> Mstar -> Mstar -> Double
normalisedInitialMassFunction cosmology i_kind m_down m_up m =
  let norm = imfNormalisation cosmology i_kind m_down m_up
   in (initialMassFunction i_kind m) / norm

-- | Mass of a remnant produced by the supernova,
-- calculated according to various works in the units of [Msol].
-- Note that below a certain mass, stellar lifetimes are higher than the
-- t_0 (age of the universe), so there will be no remnant
massRemnant :: Mstar -> Remnant_Kind -> Double -> Double
massRemnant m r_kind metal_frac =
  case r_kind of
    -- Remnant masses taken from the [Woosley & Weaver 95]
    -- and [Hoek & Groenewegen 1996]
    WoosleyWeaver
      | m < 0.9 -> 0
      | m <= 8 ->
          let (mi, mr) = (unzip . M.toList) (remnantMediumMass metal_frac)
           in makeInterp mi mr $ m
      | m <= 40 ->
          let (mi, mr) = (unzip . M.toList) (remnantHighMass metal_frac)
           in makeInterp mi mr $ m
      | otherwise ->
          extrapolate
            m
            ( (\x y -> [x, y])
                <$> M.findWithDefault 0.0 35
                <*> M.findWithDefault 0.0 40
                $ remnantHighMass
                  metal_frac
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
massDynamical :: CosmicTime -> Mstar -> Mstar
massDynamical t m_up =
  let mass_range = [0.1, 0.1 + 0.5 .. m_up]
      interp_tau = makeInterp (tauMS <$> mass_range) mass_range
   in interp_tau t

-- | Progenitor mass - proto NS mass relation from []
massPNS :: PNS_Kind -> Mstar -> Mstar
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
fractionHNe :: HNe_Kind -> Double -> Double -> Double
fractionHNe hne_kind eps_hne0 metal_frac =
  case hne_kind of
    Constant_HNe -> eps_hne0
    Grimmett -> maximum (eps_hne0 * exp (-metal_frac / 0.001), 0.001)

-- Construct yield of SNe type II by interpolating over AGB, SAGB, ECSNe and CCSNe yields
-- If neutrino-driven yields are turned on, add yields to SN II yields following the
-- initial mass - NS mass relation
yieldII :: ReferenceStarFormationModel -> Double -> Double
yieldII yields m =
  y1 m + y2 m
  where
    (mSAGB, ySAGB) = yield_sagb yields
    (mECSN, yECSN) = yield_ecsn yields
    (mII, yII)
      | m < 8.0 = yield_agb yields
      | m > 8.8 && m < 9.0 && not (null yECSN) = (mECSN, yECSN)
      | m < 11.0 && null ySAGB = yield_ccsn yields
      | m < 11.0 = (mSAGB, ySAGB)
      | otherwise = yield_ccsn yields
    y1 = uncurry makeInterp $ yield_psn yields
    y2 = uncurry makeInterp $ (mII, yII)

-- | Some of the terms (IGM/ISM outflows for all mass and a specific element yield, SFRD),
-- to be used in the next function
interGalacticMediumTerms ::
  ReferenceCosmology ->
  ReferenceStarFormationModel ->
  PowerSpectrum ->
  Remnant_Kind ->
  IMF_kind ->
  SMF_kind ->
  HMF_kind ->
  W_kind ->
  HNe_Kind ->
  Metallicity ->
  Mhalo ->
  SFRD ->
  Redshift ->
  EjectaOutflow Double
interGalacticMediumTerms cosmology@MkCosmology {prec} yields pk r_kind i_kind s_kind h_kind w_kind hne_kind metal_frac mh_min sfrd z =
  let (m_CO, eps_w, eps_sn, yr_Gyr, kms_ergMsol, energy, mUp, m_pu, m_pl, m_du_rg, m_dl_rg, m_du_ms, m_dl_ms, b_rg, b_ms, m_NSM_d, m_NSM_u, alpha_NSM, delta_t_NSM, eps_hne0, eps_MRSNe, delta_t_Novae, alpha_Novae, eps_CO, m_ej_Novae, n_Novae, m_Novae_d, m_Novae_u) =
        ( 1.38, -- Mass of the CO white dwarf
          0.02, -- Fraction of the mass that contributes to the galactic winds
          0.005, -- Fraction of the mass that contributes to the SNe outflow to the IGM
          1e9, -- Conversion factor from [yr] to [Gyr]
          1.989 * 1e43, -- Conversion factor from [km^2/s^2] to [erg/Msol]
          2 * 1e51, -- Typical kinetic energy in [erg] released by a supernovae explosion
          100, -- Upper mass limit for an IMF, in [Msol]
          8.0, -- MS SN Ia progenitor upper mass, all variables below are in [Msol]
          3.0, -- Same, but lower mass
          1.5, -- RG+WD pair progenitor upper mass
          0.9, -- Same,  but lower mass
          2.6, -- MS+WD pair progenitor upper mass
          1.8, -- Same, but lower mass
          0.02, -- Fraction of primary progenitors that produce SN Ia for RG+WD pair
          0.04, -- Same but for MS+WD pair
          9.0, -- Minimum mass of a star that can leave NSM as a remnant in [Msol]
          30, -- Similarly, maximum mass of a star that can leave an NSM in [Msol]
          0.018, -- Fraction of NS that will give rise to NS-NS system and eventually coalesce
          1e7, -- Time delay between the formation of NSM and its collapse in [yr]
          0.5, -- Fraction of HNe among CCSN with M >= 20 Msol
          0.03, -- Fraction of HNe that are turned into MRSNe
          2e9, -- Time delay between the formation of a WD and a Novae explosion in [yr]
          0.01, -- Fraction of WDs that belong to Nova systems
          0.7, -- Fraction of CO WDs that contribute towards the Nova ejecta
          2e-5, -- Average ejected mass per Nova explosion
          1e4, -- Number of Nova explosions that fit in the lifetime of a WD
          0.8, -- Lower mass bound of a star that will end up as WD
          8.0 -- Same, but upper mass
        )

      zs, ts :: [Double]
      zs = [20, 20 - 0.5 .. 0]
      ts = parMap rpar (\z -> cosmicTime cosmology z) zs

      interpT :: Double -> Double
      interpT =
        makeInterp zs ts

      zTarget :: Double -> Double -> Double
      zTarget m = makeInterp zs ts
        where
          zs = [20, 20 - 0.5 .. 0]

      vescSq :: Double
      vescSq =
        escapeVelocitySq cosmology pk h_kind w_kind mh_min z

      mDyn = massDynamical (interpT z) mUp
      mDown z = maximum [8, mDyn]

      normImf :: Double -> Double
      normImf = normalisedInitialMassFunction cosmology i_kind 0.1 mUp

      normImfSN :: Double -> Double -> Double -> Double
      normImfSN = normalisedInitialMassFunction cosmology SN_Ia

      metalFraction :: Double -> Double -> Double
      metalFraction = (metal_frac .) . zTarget

      mkIntegrand :: (Double -> Double) -> Double -> (Double -> Double) -> Double -> Double
      mkIntegrand norm offset yield m =
        norm m
          * sfrd (zTarget z m - offset)
          * yield m

      -- Helper functions
      integrand0 :: (Double -> Double) -> Double -> Double
      integrand0 = mkIntegrand normImf 0

      massEjectedAGB, massEjectedAGB_Element :: Double -> Double
      massEjectedAGB m = m - massRemnant m r_kind (metalFraction z m)
      massEjectedAGB_Element m =
        massEjectedAGB m - yieldII yields m

      yield_Novae :: Double
      yield_Novae = eps_CO * (yield_co yields) + (1 - eps_CO) * (yield_one yields)

      eps_HNe :: Double -> Double
      eps_HNe m = fractionHNe hne_kind eps_hne0 (metalFraction z m)

      -- Integrands from all contributors towards ejecta and outflows
      integrand_AGB = integrand0 massEjectedAGB
      integrand_AGB_Element = integrand0 massEjectedAGB_Element
      integrand_Wind = mkIntegrand normImf 0 (\m -> 2 * energy / (kms_ergMsol * vescSq))
      integrand_CCSN_Element m = (1 - eps_HNe m - eps_MRSNe) * mkIntegrand normImf 0 (yieldII yields) m
      integrand_NSM = mkIntegrand normImf delta_t_NSM (const 1)
      integrand_Novae_Element m = m_ej_Novae * n_Novae * alpha_Novae * yield_Novae * mkIntegrand normImf delta_t_Novae (const 1) m
      integrand_SNe_Ia_1 = normImf
      integrand_SNe_Ia_2 md mu m = mkIntegrand (normImfSN md mu) 0 (const 1)
      integrand_HNe m = (eps_HNe m - eps_MRSNe) * 1
      integrand_WR m = 1
      integrand_PISNe m = 1

      integrator = makeIntegrator prec

      -- Integrate the integrands provided above with the chosen precision
      e_AGB = integrator integrand_AGB (mDown z) mUp
      e_AGB_Element = integrator integrand_AGB_Element (mDown z) mUp
      o_Wind = eps_w * integrator integrand_Wind (mDown z) mUp
      e_CCSN_Element = integrator integrand_CCSN_Element (mDown z) mUp
      e_HNe_Element = integrator integrand_HNe (mDown z) mUp
      e_MRSNe_Element = eps_MRSNe * e_HNe_Element
      e_Novae_Element = integrator integrand_Novae_Element m_Novae_d m_Novae_u
      e_NSM = alpha_NSM * integrator integrand_NSM m_NSM_d m_NSM_u
      e_NSM_Element = yield_nsm yields * e_NSM

      first_term =
        b_rg
          * integrator (normImfSN m_dl_rg m_du_rg) (maximum [m_dl_rg, mDyn]) m_du_rg
          / integrator (normImfSN m_dl_rg m_du_rg) m_dl_rg m_du_rg
      second_term =
        b_ms
          * integrator (normImfSN m_dl_ms m_du_ms) (maximum [m_dl_ms, mDyn]) m_du_ms
          / integrator (normImfSN m_dl_ms m_du_ms) m_dl_ms m_du_ms
      e_SNe_Ia =
        m_CO
          * integrator integrand_SNe_Ia_1 (maximum [m_pl, mDyn]) m_pu
          * (first_term + second_term)
      e_SNe_Ia_Element = yield_ia yields * e_SNe_Ia
   in fmap
        (* yr_Gyr)
        MkEjectaOutflow
          { e_AGB,
            e_AGB_Element,
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
  ( e_AGB + e_SNe_Ia + e_NSM,
    e_AGB_Element + e_CCSN_Element + e_SNe_Ia_Element + e_NSM_Element,
    eps_sn * e_AGB,
    o_Wind,
    eps_sn * e_AGB_Element,
    eps_sn * e_AGB + o_Wind
  )
  where
    eps_sn = 0.005
    MkEjectaOutflow {e_AGB, e_AGB_Element, e_CCSN_Element, e_SNe_Ia, e_SNe_Ia_Element, e_NSM, e_NSM_Element, o_Wind} = rates

-- Separate IO function for IGM/ISM terms
igmTermsIO ::
  ReferenceCosmology ->
  PowerSpectrum ->
  Remnant_Kind ->
  IMF_kind ->
  SMF_kind ->
  HMF_kind ->
  W_kind ->
  HNe_Kind ->
  Metallicity ->
  Mhalo ->
  SFRD ->
  Redshift ->
  Element ->
  IO (Double, Double, Double, Double, Double, Double)
igmTermsIO cosmology pk r_kind i_kind s_kind h_kind w_kind hne_kind metal_frac mh_min sfrd z elem =
  do
    let sf_cfg = MkStarFormationCfg {model_ia = "iwamoto99/WDD1", model_ccsn = "WW95"}
    sn_ia_yield <- retrieveYieldIa sf_cfg elem
    ccsn_yield <- retrieveYieldCCSN sf_cfg elem (metal_frac z)
    let yields =
          MkStarFormation
            { yield_ia = sn_ia_yield,
              yield_ccsn = ccsn_yield,
              yield_hne = ([1], [1]),
              yield_ecsn = ([1], [1]),
              yield_agb = ([1], [1]),
              yield_sagb = ([1], [1]),
              yield_nsm = 1,
              yield_psn = ([1], [1]),
              yield_co = 1,
              yield_one = 1
            }
        result = computeRates $ interGalacticMediumTerms cosmology yields pk r_kind i_kind s_kind h_kind w_kind hne_kind metal_frac mh_min sfrd z
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
  ReferenceCosmology ->
  PowerSpectrum ->
  Remnant_Kind ->
  IMF_kind ->
  SMF_kind ->
  HMF_kind ->
  W_kind ->
  HNe_Kind ->
  Element ->
  Mhalo ->
  IO ([Double], [V.Vector Double])
igmIsmEvolution cosmology@MkCosmology {h0, om0, ob0, gn} pk r_kind i_kind s_kind h_kind w_kind hne_kind elem mh_min =
  do
    let -- Constructing stellar yields datatype as a function of metallicity of ISM and given isotope

        zs = [20.0, 20.0 - 0.5 .. 0]

        mar =
          parMap rpar (\z -> 1e9 * baryonFormationRateDensity cosmology pk h_kind w_kind z) zs
        ts =
          parMap rpar (\z -> cosmicTime cosmology z) zs
        sfrd =
          parMap rpar (\z -> starFormationRateDensity cosmology pk s_kind h_kind w_kind z mh_min) zs
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

          (e_tot, e_tot_Element, o_SNe, o_Wind, o_SNe_Element, o_tot) <- igmTermsIO cosmology pk r_kind i_kind s_kind h_kind w_kind hne_kind interpMetal mh_min interpSFRD z elem

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
    zipped_result <-
      rk4SolveHistIO odeSystem tInit ((interpT 0 - tInit) / fromIntegral nSteps) nSteps (V.fromList [igmInit, ismInit, xiIgmInit, xiIsmInit, ejectaInit, outflowInit])

    -- M_star = rho_tot - M_IGM - Mstar from the conservation equation
    -- In addition, we also normalise each mass by the total mass
    let (times, masses) = unzip zipped_result
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
  Remnant_Kind ->
  IMF_kind ->
  SMF_kind ->
  HMF_kind ->
  W_kind ->
  Mhalo ->
  ([Double], [Double], [Double])
igmIsmMetallicity cosmology pk r_kind i_kind s_kind h_kind w_kind yield_ii yield_ia yield_nsm mh_min =
  let elem = Element "H" 1

      (times, densities) =
        bimap id (transpose . (fmap V.toList)) $
          igmIsmEvolution cosmology pk r_kind i_kind s_kind h_kind w_kind yield_ii yield_ia yield_nsm elem mh_min

      result_igm = zipWith (/) (densities !! 5) (densities !! 0)
      result_ism = zipWith (/) (densities !! 4) (densities !! 1)
   in (times, result_igm, result_ism)

-}
