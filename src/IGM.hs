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
import qualified Data.Map as M
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
massDynamical :: CosmicTime -> Double
massDynamical t =
  let logt = log10 t
   in 138 * exp (1.4 * 1e-3 * sqrt (8.8753 * 1e6 + 3.2628 * 1e6 * logt))

-- | Mass of a star that dies at the age t,
-- essentially an inverse of a function tauMS
interGalacticMediumTerms :: FilePath -> ReferenceCosmology -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Mhalo -> Redshift -> IO (Double, Double)
interGalacticMediumTerms filepath cosmology i_kind s_kind h_kind w_kind mh_min z =
  let z_arr = [0, 0.25 .. 10]
      t_arr = (\z -> cosmicTime cosmology z) <$> z_arr

      e_w, e_sn :: Double
      (e_w, e_sn) = (0.02, 0.005)
   in do
        sfrd_arr <-
          mapM (\z -> starFormationRateDensity filepath cosmology s_kind h_kind w_kind z) z_arr

        vesc_sq <- escapeVelocitySq filepath cosmology h_kind w_kind mh_min z

        let interp_sfrd :: Double -> Double
            interp_sfrd = mapLookup (M.fromList (zip z_arr sfrd_arr))

            interp_time :: Double -> Double
            interp_time = mapLookup (M.fromList (zip t_arr z_arr))

            time_at_z :: Double
            time_at_z = cosmicTime cosmology z

            z_target :: Double -> Double
            z_target m = interp_time (time_at_z - tauMS m)

            m_down :: Double
            m_down = maximum [1e8, massDynamical time_at_z]

            energy :: Double
            energy = 2 * 1e51

            kms_ergMsol :: Double
            kms_ergMsol = 1.989 * 1e43

            integrand_SNe :: Double -> Double
            integrand_SNe m =
              normalisedInitialMassFunction i_kind m
                * interp_sfrd (z_target m)
                * (m - massRemnant m 0.1)

            integrand_Wind :: Double -> Double
            integrand_Wind m =
              normalisedInitialMassFunction i_kind m
                * interp_sfrd (z_target m)
                * (2 * energy / (kms_ergMsol * vesc_sq))

            result_SNe :: Double
            result_SNe = nIntegrate1024 integrand_SNe m_down 100

            result_Wind :: Double
            result_Wind = nIntegrate1024 integrand_Wind m_down 100

            result_ISM :: Double
            result_ISM = nIntegrate1024 integrand_SNe (massDynamical time_at_z) 100

        return $ (e_sn * result_SNe + e_w * result_Wind, result_ISM)

intergalacticMediumEvolution :: FilePath -> ReferenceCosmology -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Mhalo -> Redshift -> IO Double
intergalacticMediumEvolution filepath cosmology i_kind m0 =
  let z_arr = [20, 20 - 0.1 .. 0]
   in do
        terms_arr <-
          mapM (\z -> interGalacticMediumTerms filepath cosmology i_kind s_kind h_kind w_kind mh_min z) <$> z_arr

        let (o_arr, e_arr) = unzip terms_arr

            baryon_mar :: [Double]
            baryon_mar = (\z -> ob0 / om0 * massAccretionRate cosmology mh z) <$> z_arr

        return 1

main_IGM :: IO ()
main_IGM = do
  x <- interGalacticMediumTerms "../data/CAMB_Pk_z=0.txt" planck18 Kroupa DoublePower ST Smooth 1e6 0
  print $ x
