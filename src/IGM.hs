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

import Cosmology
import Data.List (unzip5)
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
  let alpha0, alpha1, alpha2, k0, k1, k2 :: Double
      (alpha0, alpha1, alpha2, m1, m2, k0, k1, k2) =
        (-0.3, -1.3, -2.3, 1, 0.08, 0.5, k0 * m1 ** (alpha0 - alpha1), k1 * m2 ** (alpha1 - alpha2))

      a_Ch, b_Ch, center_Ch, sigma_Ch :: Double
      (a_Ch, b_Ch, center_Ch, sigma_Ch) = (0.85, 0.24, 0.079, 0.69)
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
imfNormalisation :: IMF_kind -> Double
imfNormalisation i_kind =
  let m_arr = (10 **) <$> [0.1, 0.1 + 0.1 .. 100]
      integrand m = m * initialMassFunction i_kind m
   in nIntegrate512 integrand (head m_arr) (last m_arr)

normalisedInitialMassFunction :: IMF_kind -> Mstar -> Double
normalisedInitialMassFunction i_kind m =
  let norm = imfNormalisation i_kind
   in (initialMassFunction i_kind m) / norm

-- | Mass of a remnant produced by the supernova,
-- calculated according to [Iben & Tutukov 1984] in the units of [Msol]
massRemnant :: Mstar -> Metallicity -> Double
massRemnant m metal_frac
  | m >= 0.9 && m <= 8 = (mapLookup $ remnantMediumMass metal_frac) m
  | m <= 40 = (mapLookup $ remnantHighMass metal_frac) m
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

-- | Mass of a star that dies at the age t,
-- essentially an inverse of a function tauMS
interGalacticMediumTerms :: FilePath -> ReferenceCosmology -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Element -> Mhalo -> [Redshift] -> IO ([Double], [Double], [Double], [Double], [Double])
interGalacticMediumTerms filepath cosmology i_kind s_kind h_kind w_kind elem mh_min z_arr =
  let t_arr = (\z -> cosmicTime cosmology z) <$> z_arr

      e_w, e_sn :: Double
      (e_w, e_sn) = (0.02, 0.005)
   in do
        sfrd_arr <-
          mapM (\z -> starFormationRateDensity filepath cosmology s_kind h_kind w_kind z) z_arr

        vesc_sq_arr <-
          mapM (\z -> escapeVelocitySq filepath cosmology h_kind w_kind mh_min z) z_arr
        (mass_arr, yield_arr) <- yieldsHighMass 1 elem

        let interp_sfrd :: Double -> Double
            interp_sfrd = mapLookup (M.fromList (zip z_arr sfrd_arr))

            interp_vesc :: Double -> Double
            interp_vesc = mapLookup (M.fromList (zip z_arr vesc_sq_arr))

            interp_yields :: Double -> Double
            interp_yields = mapLookup (M.fromList (zip mass_arr yield_arr))

            interp_time :: Double -> Double
            interp_time = mapLookup (M.fromList (zip t_arr z_arr))

            massDynamical :: Double -> Double
            massDynamical t =
              let mass_range = [0.1, 0.1 + 0.5 .. 100]
                  interp_tau = mapLookup (M.fromList (zip (tauMS <$> mass_range) mass_range))
               in interp_tau t

            z_target :: Double -> Double -> Double
            z_target z m = interp_time (cosmicTime cosmology z - tauMS m)

            m_down :: Double -> Double
            m_down z = maximum [1e8, (massDynamical (cosmicTime cosmology z))]

            energy :: Double
            energy = 2 * 1e51

            kms_ergMsol :: Double
            kms_ergMsol = 1.989 * 1e43

            integrand_SNe :: Double -> Double -> Double
            integrand_SNe z m =
              normalisedInitialMassFunction i_kind m
                * interp_sfrd (z_target z m)
                * (m - massRemnant m 0.1)

            integrand_SNe_Element :: Double -> Double -> Double
            integrand_SNe_Element z m =
              normalisedInitialMassFunction i_kind m
                * interp_sfrd (z_target z m)
                * interp_yields m
                * (m - massRemnant m 0.1)

            integrand_Wind :: Double -> Double -> Double
            integrand_Wind z m =
              normalisedInitialMassFunction i_kind m
                * interp_sfrd (z_target z m)
                * (2 * energy / (kms_ergMsol * interp_vesc z))

            integrand_ISM_Element :: Double -> Double -> Double
            integrand_ISM_Element = integrand_SNe_Element

            result_SNe :: [Double]
            result_SNe =
              (\z -> e_sn * nIntegrate1024 (integrand_SNe z) (m_down z) 100) <$> z_arr

            result_SNe_Element :: [Double]
            result_SNe_Element =
              (\z -> e_sn * nIntegrate1024 (integrand_SNe_Element z) (m_down z) 100) <$> z_arr

            result_Wind :: [Double]
            result_Wind =
              (\z -> e_w * nIntegrate1024 (integrand_Wind z) (m_down z) 100) <$> z_arr

            result_ISM :: [Double]
            result_ISM =
              (\z -> nIntegrate1024 (integrand_SNe z) (massDynamical (cosmicTime cosmology z)) 100) <$> z_arr

            result_ISM_Element :: [Double]
            result_ISM_Element =
              (\z -> nIntegrate1024 (integrand_ISM_Element z) (massDynamical (cosmicTime cosmology z)) 100) <$> z_arr

        return $ (result_SNe, result_SNe_Element, result_Wind, result_ISM, result_ISM_Element)

-- | Coupled solver ...
igmIsmEvolution :: FilePath -> ReferenceCosmology -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Element -> Mhalo -> IO [(Double, V.Vector Double)]
igmIsmEvolution filepath cosmology i_kind s_kind h_kind w_kind elem mh_min =
  let z_arr = [20.0, 20.0 - 0.1 .. 0]
      m_tot = 6 * 1e22
   in do
        terms_arr <-
          interGalacticMediumTerms filepath cosmology i_kind s_kind h_kind w_kind elem mh_min z_arr

        sfrd_arr <-
          mapM (\z -> starFormationRateDensity filepath cosmology s_kind h_kind w_kind z) z_arr

        let (h0, om0, ob0, c, gn) = unpackCosmology cosmology

            (osn_arr, osni_arr, ow_arr, e_arr, ei_arr) = terms_arr

            baryon_mar :: [Double]
            baryon_mar = (\z -> ob0 / om0 * massAccretionRate cosmology m_tot z) <$> z_arr

            interp_osn = mapLookup $ M.fromList (zip z_arr osn_arr)
            interp_osni = mapLookup $ M.fromList (zip z_arr osni_arr)
            interp_ow = mapLookup $ M.fromList (zip z_arr ow_arr)
            interp_e = mapLookup $ M.fromList (zip z_arr e_arr)
            interp_ei = mapLookup $ M.fromList (zip z_arr ei_arr)
            interp_sfrd = mapLookup $ M.fromList (zip z_arr sfrd_arr)
            interp_o = (+) <$> interp_osn <*> interp_ow
            interp_mar = mapLookup $ M.fromList (zip z_arr baryon_mar)

            -- M_ISM = y!0
            -- M_IGM = y!1
            -- Xi_ISM = y!2
            -- Xi_IGM = y!3
            igm_ode z y =
              V.fromList
                [ -interp_mar z + interp_o z,
                  (-interp_sfrd z + interp_e z) + (interp_mar z - interp_o z),
                  1 / (y V.! 0) * (interp_ow z * (y V.! 3 - y V.! 2) + (interp_osni z - interp_osn z * y V.! 2)),
                  1 / (y V.! 1) * ((interp_ei z - interp_e z * y V.! 3) + interp_mar z * (y V.! 2 - y V.! 3) - (interp_osni z - interp_osn z * y V.! 3))
                ]
            result = rk4Solve igm_ode (minimum z_arr) 0.1 (length z_arr) (V.fromList [1, 1, 0.01, 0.01])

        return result

main_IGM :: IO ()
main_IGM = do
  x <- igmIsmEvolution "data/CAMB_Pk_z=0.txt" planck18 Kroupa DoublePower ST Smooth (Element "C" 12) 1e6
  print $ x
