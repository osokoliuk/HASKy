module IGM where

import Cosmology
import qualified Data.Map as M
import HMF
import Helper
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

      a_Ch, center_Ch, sigma_Ch :: Double
      (a_Ch, center_Ch, sigma_Ch) = (0.086, 0.57, 0.22)
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
          | otherwise -> m ** (-1.3)

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
massRemnant :: Mstar -> Double
massRemnant m
  | m <= 6.8 = 0.11 * m + 0.4
  | otherwise = 1.5

-- | Lifetime of a main sequence star in relation to it's mass,
-- taken from the work of [Reid et al. 2002]
tauMS :: Mstar -> CosmicTime
tauMS m = 10 ** (1.015 - 3.491 * log10 m + 0.8157 * (log10 m) ** 2)

-- | Mass of a star that dies at the age t,
-- essentially an inverse of a function tauMS
massDynamical :: CosmicTime -> Double
massDynamical t =
  let logt = log10 t
   in 138 * exp (1.4 * 1e-3 * sqrt (8.8753 * 1e6 + 3.2628 * 1e6 * logt))

-- | Mass of a star that dies at the age t,
-- essentially an inverse of a function tauMS
supernovaEjecta :: FilePath -> ReferenceCosmology -> IMF_kind -> SMF_kind -> HMF_kind -> W_kind -> Redshift -> IO Double
supernovaEjecta filepath cosmology i_kind s_kind h_kind w_kind z =
  let z_arr = [0, 0.25 .. 10]
      t_arr = (\z -> cosmicTime cosmology z) <$> z_arr
   in do
        sfrd_arr <-
          mapM (\z -> starFormationRateDensity filepath cosmology s_kind h_kind w_kind z) z_arr

        let interp_sfrd :: Double -> Double
            interp_sfrd = mapLookup (M.fromList (zip z_arr sfrd_arr))

            interp_time :: Double -> Double
            interp_time = mapLookup (M.fromList (zip t_arr z_arr))

            time_at_z :: Double
            time_at_z = cosmicTime cosmology z

            z_target :: Double -> Double
            z_target m = interp_time (time_at_z - tauMS m)

            integrand m =
              normalisedInitialMassFunction i_kind m
                * interp_sfrd (z_target m)
                * (m - massRemnant m)

            result :: Double
            result = nIntegrate1024 integrand (massDynamical time_at_z) 100

        return $ result

main_IGM :: IO ()
main_IGM = do
  x <- supernovaEjecta "data/CAMB_Pk_z=0.txt" planck18 Kroupa DoublePower ST Smooth 0
  print $ x
