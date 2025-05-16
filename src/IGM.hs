{-# LANGUAGE MultiWayIf #-}

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

import Control.Parallel.Strategies
import Cosmology
import Data.Bifunctor
import qualified Data.Map as M
import qualified Data.Vector as V
import HMF
import Helper
import Lookup
import Math.GaussianQuadratureIntegration
import SMF

data IMF_kind
  = Salpeter
  | Kroupa
  | Chabrier
  | SN_Ia
  deriving (Eq, Show)

data Remnant_Kind
  = Pereira
  | WoosleyWeaver
  deriving (Eq, Show)

initialMassFunction :: IMF_kind -> Mstar -> Double
initialMassFunction i_kind m =
  let (alpha0, alpha1, alpha2, m1, m2, k0, k1, k2) =
        (-0.3, -1.3, -2.3, 1, 0.08, 0.5, k0 * m1 ** (alpha0 - alpha1), k1 * m2 ** (alpha1 - alpha2))
      (a_Ch, b_Ch, center_Ch, sigma_Ch) =
        (0.85, 0.24, 0.079, 0.69)
   in case i_kind of
        -- Salpeter et al. 1955 IMF (single power-law)
        Salpeter -> m ** (-2.35)
        -- Kroupa et al. 2001 IMF (broken power-law)
        Kroupa
          | m < m1 -> k0 * m ** alpha0
          | m >= m2 && m < m2 -> k1 * m ** alpha1
          | otherwise -> k2 * m ** alpha2
        -- Chabrier et al. 2003 (log-normal) converted to [Mpc^-3]
        Chabrier
          | m < 1 -> a_Ch * exp (-(log m - log center_Ch) ** 2 / (2 * sigma_Ch ** 2))
          | otherwise -> b_Ch * m ** (-1.3)
        SN_Ia -> m ** (-0.35)

-- | Next two functions normalise IMF between m_inf = 0.1 Msol and m_sup = 100 Msol
-- so that it can acts as a PDF in that range
imfNormalisation :: ReferenceCosmology -> IMF_kind -> Mstar -> Mstar -> Double
imfNormalisation cosmology i_kind m_down m_up =
  let (_, _, _, _, _, _, _, prec) = unpackCosmology cosmology
      m_arr = [m_down, m_down + 0.1 .. m_up]
      integrand m = m * initialMassFunction i_kind m
      integrator = makeIntegrator (Precision prec)
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

-- | Some of the terms (IGM/ISM outflows for all mass and a specific element yield, SFRD),
-- to be used in the next function
interGalacticMediumTerms :: ReferenceCosmology -> PowerSpectrum -> Remnant_Kind -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Yield_II -> Yield_Ia -> Metallicity -> Mhalo -> SFRD -> Redshift -> (Double, Double, Double, Double, Double, Double)
interGalacticMediumTerms cosmology pk r_kind i_kind s_kind h_kind w_kind yield_ii yield_ia metal_frac mh_min sfrd z =
  let (_, _, _, _, _, _, _, prec) = unpackCosmology cosmology

      (m_CO, eps_w, eps_sn, kms_ergMsol, energy, m_up, m_pu, m_pl, m_du_rg, m_dl_rg, m_du_ms, m_dl_ms, b_rg, b_ms) =
        ( 1.38, -- Mass of the CO white dwarf
          0.02, -- Fraction of the mass that contributes to the galactic winds
          0.005, -- Fraction of the mass that contributes to the SNe outflow to the IGM
          1.989 * 1e43, -- Conversion factor from [km^2/s^2] to [erg/Msol]
          2 * 1e51, -- Typical kinetic energy in [erg] released by a supernovae explosion
          100, -- Upper mass limit for an IMF, in [Msol]
          8, -- MS SN Ia progenitor upper mass, all variables in below are in [Msol]
          3, -- Same, but lower mass
          1.5, -- RG+WD pair progenitor upper mass
          0.9, -- Same,  but lower mass
          2.6, -- MS+WD pair progenitor upper mas
          1.8, -- Same, but lower mass
          0.02, -- Fraction of primary progenitors that produce SN Ia for RG+WD pair
          0.05 -- Same but for MS+WD pair
        )

      z_arr = [20, 20 - 0.5 .. 0]

      vesc_sq =
        escapeVelocitySq cosmology pk h_kind w_kind mh_min z
      t_arr =
        parMap rpar (\z -> cosmicTime cosmology z) z_arr

      interp_t = makeInterp z_arr t_arr
      z_target z m = interp_t (cosmicTime cosmology z - tauMS m)

      mdyn = massDynamical (interp_t z) m_up
      m_down z = maximum [8, mdyn]
      norm_imf = normalisedInitialMassFunction cosmology i_kind 0.1 m_up
      norm_imf_sn md mu = normalisedInitialMassFunction cosmology SN_Ia md mu

      -- Ejecta by mass loss and SNe II
      integrand_loss m =
        norm_imf m
          * sfrd (z_target z m)
          * (1 - massRemnant m r_kind (metal_frac (z_target z m)))

      -- Ejecta per element from mass loss
      integrand_loss_Element m =
        norm_imf m
          * sfrd (z_target z m)
          * (1 - massRemnant m r_kind (metal_frac (z_target z m)) - yield_ii m (metal_frac (z_target z m)) / m)
          * metal_frac (z_target z m) -- ISM metal fraction at t - tauMS

      -- Ejecta per element from SNe II
      integrand_SNe_II_Element m =
        norm_imf m
          * sfrd (z_target z m)
          * yield_ii m (metal_frac (z_target z m))
          / m

      -- Ejecta from SNe Ia
      integrand_SNe_Ia_1 m =
        norm_imf m
          * (1 / m)
      integrand_SNe_Ia_2 md mu m =
        norm_imf_sn md mu m
          * sfrd (z_target z m)
          * (1 / m)

      -- IGM outflows from galactic winds
      integrand_Wind m =
        norm_imf m
          * sfrd (z_target z m)
          * (2 * energy / (kms_ergMsol * vesc_sq))

      integrator = makeIntegrator (Precision prec)

      -- Integrate the integrands provided above with the chosen precision
      e_loss =
        integrator integrand_loss (m_down z) m_up
      e_loss_Element =
        integrator integrand_loss_Element (m_down z) m_up
      o_Wind =
        eps_w * integrator integrand_Wind (m_down z) m_up
      e_SNe_II_Element = integrator integrand_SNe_II_Element (m_down z) m_up
      e_SNe_Ia =
        m_CO
          * integrator integrand_SNe_Ia_1 (maximum [m_pl, mdyn]) m_pu
          * ( b_rg * integrator (integrand_SNe_Ia_2 m_dl_rg m_du_rg) (maximum [m_dl_rg, mdyn]) m_du_rg
                + b_ms * integrator (integrand_SNe_Ia_2 m_dl_ms m_du_ms) (maximum [m_dl_ms, mdyn]) m_du_ms
            )
      e_SNe_Ia_Element =
        yield_ia / m_CO * e_SNe_Ia
   in (e_loss, e_loss_Element, e_SNe_II_Element, e_SNe_Ia, e_SNe_Ia_Element, o_Wind)

-- | Solve four copled first-order differential equations that govern the evolution of:
--    * M_IGM     (1)
--    * M_ISM     (3)
--    * Xi_IGM    (4)
--    * Xi_ISM    (5)
-- with all equations being taken from the [Daigne et al. 2004]
igmIsmEvolution :: ReferenceCosmology -> PowerSpectrum -> Remnant_Kind -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Yield_II -> Yield_Ia -> Element -> Mhalo -> ([Double], [V.Vector Double])
igmIsmEvolution cosmology pk r_kind i_kind s_kind h_kind w_kind yield_ii yield_ia elem mh_min =
  let (h0, om0, ob0, _, gn, _, _, _) = unpackCosmology cosmology
      z_arr = [20.0, 20.0 - 0.5 .. 0]

      terms_arr metal_frac sfrd z =
        interGalacticMediumTerms cosmology pk r_kind i_kind s_kind h_kind w_kind yield_ii yield_ia metal_frac mh_min sfrd z
      mar_arr =
        parMap rpar (\z -> 1e9 * baryonFormationRateDensity cosmology pk h_kind w_kind z) z_arr
      t_arr =
        parMap rpar (\z -> cosmicTime cosmology z) z_arr
      sfrd_arr =
        parMap rpar (\z -> 1e9 * starFormationRateDensity cosmology pk s_kind h_kind w_kind z) z_arr

      -- Unpack outflow/inflow rates and interpolate over our redshift range
      interp_mar = makeInterp z_arr mar_arr
      interp_sfrd = makeInterp z_arr sfrd_arr
      interp_z = makeInterp t_arr z_arr
      interp_t = makeInterp z_arr t_arr

      -- ICs are set assuming very small baryon fraction in the structures,
      -- with M_ISM/M_IGM ~ 0.01 (stellar mass is negligible at this redshift)
      -- following the prescription of [Daigne et al. 2006]
      -- Finally, we also adopt the BBN abundances for H (He),
      -- such that the ICs for Xi_ISM/Xi_IGM = 0.76 (0.24) * M_ISM/M_IGM.
      (n_steps, t_init, a_ini, mass_tot, igm_ini, ism_ini, xi_igm_ini, xi_ism_ini) =
        -- Nucleosynthesis abundances are taken from the [Coc et al. 2014]
        let yp = 0.2464 -- Primordial Helium abundance
            fh = 1 - yp -- Primordial hydrogen abudance
            fd = 2.64 * 1e-5 -- Deuterium fraction
            fh3 = 1.05 * 1e-5 * (1 - fd) * fh
            fli7 = 5.18 * 1e-10 * (1 - fd) * fh
            ini_abundance
              | element elem == "H" =
                  if
                    | isotope elem == 1 -> (1 - fd) * fh
                    | otherwise -> fd * fh
              | element elem == "He" =
                  if
                    | isotope elem == 4 -> (1 - fh3) * yp
                    | otherwise -> fh3
              | element elem == "Li" && isotope elem == 7 = fli7
              | otherwise = 0
         in (10 :: Int, interp_t (maximum z_arr), 0.01, 1e11, (1 - a_ini) * mass_tot, a_ini * mass_tot, ini_abundance * igm_ini, ini_abundance * ism_ini)

      -- Convert Differential-Algebraic system into a system of ODEs
      igm_ode :: History -> Double -> V.Vector Double -> V.Vector Double
      igm_ode history t y =
        let (times, metals) =
              unzip $
                [(t, metal) | (t, v) <- history, V.length v > 2, let metal = v V.! 2]
            interp_metal t =
              if length times < 1
                then 0
                else
                  (makeInterp times metals) t

            -- Unpack all rates at z(t)
            z = interp_z t
            (e_loss, e_loss_Element, e_SNe_II_Element, e_SNe_Ia, e_SNe_Ia_Element, o_Wind) =
              (terms_arr interp_metal interp_sfrd z)

            -- Ejecta
            e_tot = e_loss + e_SNe_Ia
            e_tot_Element = e_loss_Element + e_SNe_II_Element + e_SNe_Ia_Element

            -- Outflows
            eps_sn = 0.005
            o_SNe = eps_sn * e_tot
            o_SNe_Element = eps_sn * e_tot_Element
            o_tot = o_SNe + o_Wind
         in V.fromList
              [ -interp_mar z + o_tot,
                (-interp_sfrd z + e_tot) + (interp_mar z - o_tot),
                1 / (y V.! 0) * (o_Wind * (y V.! 3 - y V.! 2) + (o_SNe_Element - o_SNe * y V.! 2)),
                1 / (y V.! 1) * ((e_tot_Element - e_tot * y V.! 3) + interp_mar z * (y V.! 2 - y V.! 3) - (o_SNe_Element - o_SNe * y V.! 3))
              ]

      -- Solve the system and unpack values
      (times, masses) =
        unzip $
          rk4SolveHist igm_ode t_init ((interp_t 0 - t_init) / fromIntegral n_steps) n_steps (V.fromList [igm_ini, ism_ini, xi_igm_ini, xi_ism_ini])

      -- M_star = rho_tot - M_IGM - Mstar from the conservation equation
      -- In addition, we also normalise each mass by the total mass
      result = (\v -> V.map (/ mass_tot) $ V.snoc v (mass_tot - v V.! 0 - v V.! 1)) <$> masses
   in (times, result)

{-
-- | Derive the metallicity of the IGM from outflow/inflow rates of metals,
-- currently using an approach presented in [Tan et al. 2018]
igmMetallicity :: ReferenceCosmology -> PowerSpectrum -> Remnant_Kind -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Yield_Ia -> Yield_II -> Metallicity -> Mhalo -> ([Double], [Double])
igmMetallicity cosmology pk r_kind i_kind s_kind h_kind w_kind yield_ia yield_ii metal_frac mh_min =
  let (h0, om0, ob0, _, gn, _, _, _) = unpackCosmology cosmology
      z_arr = [20.0, 20.0 - 0.25 .. 0]
      rho_cr = 3 * h0 ** 2 / (8 * pi * gn)

      terms_arr =
        interGalacticMediumTerms cosmology pk r_kind i_kind s_kind h_kind w_kind yield metal_frac mh_min z_arr
      mar_arr =
        parMap rpar (\z -> 1e9 * baryonFormationRateDensity cosmology pk h_kind w_kind z) z_arr
      t_arr =
        parMap rpar (\z -> cosmicTime cosmology z) z_arr

      (_, osn_arr, _, ow_arr, _, _) = terms_arr
      o_arr = zipWith (+) osn_arr ow_arr

      numerator = cumulativeTrapezoid t_arr o_arr
      denominator = zipWith (\x y -> rho_cr * ob0 + x - y) (cumulativeTrapezoid t_arr mar_arr) numerator
      result = zipWith (\x y -> x / y) numerator denominator
   in (z_arr, result)

-}
