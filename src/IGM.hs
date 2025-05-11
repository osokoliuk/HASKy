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

-- | Next two functions normalise IMF between m_inf = 0.1 Msol and m_sup = 100 Msol
-- so that it can acts as a PDF in that range
imfNormalisation :: IMF_kind -> Mstar -> Double
imfNormalisation i_kind mup =
  let m_arr = [0.1, 0.1 + 0.1 .. mup]
      integrand m = m * initialMassFunction i_kind m
   in nIntegrate256 integrand (head m_arr) (last m_arr)

normalisedInitialMassFunction :: IMF_kind -> Mstar -> Mstar -> Double
normalisedInitialMassFunction i_kind mup m =
  let norm = imfNormalisation i_kind mup
   in (initialMassFunction i_kind m) / norm

-- | Mass of a remnant produced by the supernova,
-- calculated according to various works in the units of [Msol].
-- Note that below a certain mass, stellar lifetimes are higher than the
-- t_0 (age of the universe), so there will be no remnant
massRemnant :: Mstar -> Metallicity -> Double
massRemnant m metal_frac
  | m < 0.9 = 0
  | m >= 0.9 && m <= 8 =
      let (mi, mr) = (unzip . M.toList) (remnantMediumMass metal_frac)
       in makeInterp mi mr $ m
  | m <= 40 =
      let (mi, mr) = (unzip . M.toList) (remnantHighMass metal_frac)
       in makeInterp mi mr $ m
  | otherwise =
      extrapolate
        m
        ( (\x y -> [x, y])
            <$> M.findWithDefault 0.0 35
            <*> M.findWithDefault 0.0 40
            $ remnantHighMass
              metal_frac
        )
        [35, 40]

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

-- | Some of the terms (IGM/ISM outflows for all mass and a specific element yield, SFRD),
-- to be used in the next function
interGalacticMediumTerms :: ReferenceCosmology -> PowerSpectrum -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Yield -> Mhalo -> [Redshift] -> ([Double], [Double], [Double], [Double], [Double], [Double])
interGalacticMediumTerms cosmology pk i_kind s_kind h_kind w_kind yield mh_min z_arr =
  let (e_w, e_sn, kms_ergMsol, yr_Gyr, energy, m_up) =
        (0.02, 0.005, 1.989 * 1e43, 1e9, 2 * 1e51, 100)

      t_arr =
        parMap rpar (\z -> cosmicTime cosmology z) z_arr
      sfrd_arr =
        parMap rpar (\z -> yr_Gyr * starFormationRateDensity cosmology pk s_kind h_kind w_kind z) z_arr
      vesc_sq_arr =
        parMap rpar (\z -> escapeVelocitySq cosmology pk h_kind w_kind mh_min z) z_arr

      sfrd = makeInterp z_arr sfrd_arr
      vesc_sq = makeInterp z_arr vesc_sq_arr
      time = makeInterp t_arr z_arr

      massDynamical t =
        let mass_range = [0.1, 0.1 + 0.5 .. 100]
            interp_tau = makeInterp (tauMS <$> mass_range) mass_range
         in interp_tau t

      z_target z m = time (cosmicTime cosmology z - tauMS m)
      m_down z = maximum [8, (massDynamical (cosmicTime cosmology z))]
      norm_imf = normalisedInitialMassFunction i_kind m_up

      integrand_SNe z m =
        norm_imf m
          * sfrd (z_target z m)
          * (m - massRemnant m 0.1)
      integrand_SNe_Element z m =
        norm_imf m
          * sfrd (z_target z m)
          * yield m
          * (m - massRemnant m 0.1)
      integrand_Wind z m =
        norm_imf m
          * sfrd (z_target z m)
          * (2 * energy / (kms_ergMsol * vesc_sq z))
      integrand_ISM_Element = integrand_SNe_Element

      result_SNe =
        parMap rpar (\z -> e_sn * nIntegrate256 (integrand_SNe z) (m_down z) m_up) z_arr
      result_SNe_Element =
        parMap rpar (\z -> e_sn * nIntegrate256 (integrand_SNe_Element z) (m_down z) m_up) z_arr
      result_Wind =
        parMap rpar (\z -> e_w * nIntegrate256 (integrand_Wind z) (m_down z) m_up) z_arr
      result_ISM =
        parMap rpar (\z -> nIntegrate256 (integrand_SNe z) (massDynamical (cosmicTime cosmology z)) m_up) z_arr
      result_ISM_Element =
        parMap rpar (\z -> nIntegrate256 (integrand_ISM_Element z) (massDynamical (cosmicTime cosmology z)) m_up) z_arr
   in (sfrd_arr, result_SNe, result_SNe_Element, result_Wind, result_ISM, result_ISM_Element)

-- | Solve four copled first-order differential equations that govern the evolution of:
--    * M_IGM     (1)
--    * M_ISM     (3)
--    * Xi_IGM    (4)
--    * Xi_ISM    (5)
-- with all equations being taken from the [Daigne et al. 2004]
igmIsmEvolution :: ReferenceCosmology -> PowerSpectrum -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Yield -> Mhalo -> ([Double], [V.Vector Double])
igmIsmEvolution cosmology pk i_kind s_kind h_kind w_kind yield mh_min =
  let (h0, om0, ob0, _, gn, _, _) = unpackCosmology cosmology
      z_arr = [19.0, 19.0 - 0.1 .. 0]
      rho_mean z = 3 * h0 ** 2 * ob0 / (8 * pi * gn) * (1 + z) ** (-4)
      rho_tot = 1e11

      terms_arr =
        interGalacticMediumTerms cosmology pk i_kind s_kind h_kind w_kind yield mh_min z_arr
      mar_arr =
        parMap rpar (\z -> 1e9 * baryonFormationRateDensity cosmology pk h_kind w_kind z) z_arr
      t_arr =
        parMap rpar (\z -> cosmicTime cosmology z) z_arr

      (interp_sfrd, interp_osn, interp_osni, interp_ow, interp_e, interp_ei) =
        mapTuple6 (makeInterp z_arr) terms_arr
      interp_mar = makeInterp z_arr mar_arr
      interp_z = makeInterp t_arr z_arr
      interp_t = makeInterp z_arr t_arr
      interp_o = (+) <$> interp_osn <*> interp_ow

      -- ICs are set assuming very small baryon fraction in the structures,
      -- with M_ISM/M_IGM ~ 0.01 (stellar mass is negligible at this redshift)
      -- following the prescription of [Daigne et al. 2006]
      (n_steps, t_init, a_ini, igm_ini, ism_ini, xi_igm_ini, xi_ism_ini) =
        (20 :: Int, interp_t (maximum z_arr), 0.01, (1 - a_ini) * rho_tot, a_ini * rho_tot, 0, 0)

      -- Convert Differential-Algebraic system into ODE via Lagrangian multipliers
      -- We are using an approach laid out in [van der Houwen & de Swart 1997]
      igm_ode t y =
        let z = interp_z t
         in V.fromList
              [ -interp_mar z + interp_o z,
                (-interp_sfrd z + interp_e z) + (interp_mar z - interp_o z),
                1 / (y V.! 0) * (interp_ow z * (y V.! 3 - y V.! 2) + (interp_osni z - interp_osn z * y V.! 2)),
                1 / (y V.! 1) * ((interp_ei z - interp_e z * y V.! 3) + interp_mar z * (y V.! 2 - y V.! 3) - (interp_osni z - interp_osn z * y V.! 3))
              ]

      (times, masses) =
        unzip $
          rk4Solve igm_ode t_init ((interp_t 0 - t_init) / fromIntegral n_steps) n_steps (V.fromList [igm_ini, ism_ini, xi_igm_ini, xi_ism_ini])

      -- M_ISM = M_IGM - Mstar - rho_tot from the conservation equation
      result = (\v -> V.snoc v (rho_tot - v V.! 0 - v V.! 1)) <$> masses
   in (times, result)
