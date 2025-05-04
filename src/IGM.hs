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

interGalacticMediumTerms :: ReferenceCosmology -> PowerSpectrum -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Yield -> Mhalo -> [Redshift] -> ([Double], [Double], [Double], [Double], [Double], [Double])
interGalacticMediumTerms cosmology pk i_kind s_kind h_kind w_kind yield mh_min z_arr =
  let t_arr = (\z -> cosmicTime cosmology z) <$> z_arr

      (e_w, e_sn) = (0.02, 0.005)
      kms_ergMsol = 1.989 * 1e43
      energy = 2 * 1e51

      sfrd_arr =
        (\z -> starFormationRateDensity cosmology pk s_kind h_kind w_kind z) <$> z_arr
      vesc_sq_arr =
        (\z -> escapeVelocitySq cosmology pk h_kind w_kind mh_min z) <$> z_arr

      sfrd = mapLookup (M.fromList (zip z_arr sfrd_arr))
      vesc_sq = mapLookup (M.fromList (zip z_arr vesc_sq_arr))

      massDynamical :: Double -> Double
      massDynamical t =
        let mass_range = [0.1, 0.1 + 0.5 .. 100]
            interp_tau = mapLookup (M.fromList (zip (tauMS <$> mass_range) mass_range))
         in interp_tau t

      interp_time :: Double -> Double
      interp_time = mapLookup (M.fromList (zip t_arr z_arr))

      z_target :: Double -> Double -> Double
      z_target z m = interp_time (cosmicTime cosmology z - tauMS m)

      m_down :: Double -> Double
      m_down z = maximum [1e8, (massDynamical (cosmicTime cosmology z))]

      kms_ergMsol :: Double

      integrand_SNe :: Double -> Double -> Double
      integrand_SNe z m =
        normalisedInitialMassFunction i_kind m
          * sfrd (z_target z m)
          * (m - massRemnant m 0.1)

      integrand_SNe_Element :: Double -> Double -> Double
      integrand_SNe_Element z m =
        normalisedInitialMassFunction i_kind m
          * sfrd (z_target z m)
          * yield m
          * (m - massRemnant m 0.1)

      integrand_Wind :: Double -> Double -> Double
      integrand_Wind z m =
        normalisedInitialMassFunction i_kind m
          * sfrd (z_target z m)
          * (2 * energy / (kms_ergMsol * vesc_sq z))

      integrand_ISM_Element :: Double -> Double -> Double
      integrand_ISM_Element = integrand_SNe_Element

      result_SNe :: [Double]
      result_SNe =
        (\z -> e_sn * nIntegrate256 (integrand_SNe z) (m_down z) 100) <$> z_arr

      result_SNe_Element :: [Double]
      result_SNe_Element =
        (\z -> e_sn * nIntegrate256 (integrand_SNe_Element z) (m_down z) 100) <$> z_arr

      result_Wind :: [Double]
      result_Wind =
        (\z -> e_w * nIntegrate256 (integrand_Wind z) (m_down z) 100) <$> z_arr

      result_ISM :: [Double]
      result_ISM =
        (\z -> nIntegrate256 (integrand_SNe z) (massDynamical (cosmicTime cosmology z)) 100) <$> z_arr

      result_ISM_Element :: [Double]
      result_ISM_Element =
        (\z -> nIntegrate256 (integrand_ISM_Element z) (massDynamical (cosmicTime cosmology z)) 100) <$> z_arr
   in (sfrd_arr, result_SNe, result_SNe_Element, result_Wind, result_ISM, result_ISM_Element)

-- | Solve four copled first-order differential equations that govern the evolution of:
--    * M_IGM   (1)
--    * M_ISM   (2)
--    * Xi_IGM  (3)
--    * Xi_ISM  (4)
-- with all equations being taken from the [Daigne et al. 2004]
igmIsmEvolution :: ReferenceCosmology -> PowerSpectrum -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Yield -> Mhalo -> [(Double, V.Vector Double)]
igmIsmEvolution cosmology pk i_kind s_kind h_kind w_kind yield mh_min =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology
      z_arr = [20.0, 20.0 - 0.1 .. 0]
      m_tot = 6 * 1e22

      terms_arr =
        interGalacticMediumTerms cosmology pk i_kind s_kind h_kind w_kind yield mh_min z_arr

      baryon_mar :: [Double]
      baryon_mar = (\z -> ob0 / om0 * massAccretionRate cosmology m_tot z) <$> z_arr

      (interp_sfrd, interp_osn, interp_osni, interp_ow, interp_e, interp_ei) =
        mapTuple7 (makeInterp z_arr) terms_arr

      interp_o = (+) <$> interp_osn <*> interp_ow
      interp_mar = mapLookup $ M.fromList (zip z_arr baryon_mar)

      igm_ode z y =
        V.fromList
          [ -interp_mar z + interp_o z,
            (-interp_sfrd z + interp_e z) + (interp_mar z - interp_o z),
            1 / (y V.! 0) * (interp_ow z * (y V.! 3 - y V.! 2) + (interp_osni z - interp_osn z * y V.! 2)),
            1 / (y V.! 1) * ((interp_ei z - interp_e z * y V.! 3) + interp_mar z * (y V.! 2 - y V.! 3) - (interp_osni z - interp_osn z * y V.! 3))
          ]
      result = rk4Solve igm_ode (minimum z_arr) 0.1 (length z_arr) (V.fromList [1, 1, 0.01, 0.01])
   in result

main_IGM :: IO ()
main_IGM =
  do
    (k_arr, pk_arr) <- powerSpectrum "data/CAMB_Pk_z=0.txt"

    (mass_arr, yield_arr) <- yieldsHighMass 1 $ Element "C" 12

    let interp_pk :: PowerSpectrum
        interp_pk = mapLookup (M.fromList (zip k_arr pk_arr))

        interp_yield :: Yield
        interp_yield = mapLookup (M.fromList (zip mass_arr yield_arr))

        x = igmIsmEvolution planck18 interp_pk Kroupa DoublePower ST Smooth interp_yield 1e6
    print x
