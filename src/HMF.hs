module HMF where

import Cosmology
import Data.Maybe
import qualified Data.Vector as V
import Numeric.Tools.Differentiation
import Numeric.Tools.Integration
import Numeric.Tools.Interpolation

-- Define an HMF datatype, consisting of two choices
data HMF_kind
  = Tinker
  | ST
  deriving (Eq, Show)

-- Similarly, define the options for the window function
data W_kind
  = TopHat
  | Smooth
  deriving (Eq, Show)

-- Some simple types added for convenience
type Wavenumber = Double

type Mhalo = Double

type Rhalo = Double

type Redshift = Double

-- Unpack the values in the cosmology_record into variables
unpackCosmology :: ReferenceCosmology -> (Double, Double, Double, Double, Double)
unpackCosmology cosmology =
  (,,,,) <$> h0' <*> om0' <*> ob0' <*> c' <*> gn' $ cosmology

-- | Define a radius for the uniform density sphere in terms of it's mass
rh :: ReferenceCosmology -> Mhalo -> Rhalo
rh cosmology mh =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology
   in 1 / c * (2 * mh * gn / (om0 * h0 ** 2)) ** (1 / 3)

-- | Produces a matter power spectrum at the linear level
-- To be imported from CAMB
powerSpectrum :: Wavenumber -> Redshift -> Double
powerSpectrum k z =
  let xPoints = [1, 2, 3, 4]
      yPoints = [2, 4, 6, 8]
      interp = linearInterp $ tabulate (uniformMesh (1e-4, 1e4) 4) (V.fromList xPoints)
   in interp `at` k

-- | Window function, used to derive the cosmic variance
-- You have a choice of two different ones, namely:
--    * Top-Hat
--    * Smooth-k
windowFunction :: ReferenceCosmology -> W_kind -> Wavenumber -> Mhalo -> Double
windowFunction cosmology w_kind k mh =
  let r = rh cosmology mh
      kr = k * r
      beta = 4.8
   in case w_kind of
        TopHat -> 3 * (sin (kr) - kr * cos (kr)) / (kr) ** 3
        Smooth -> (1 + kr ** beta) ** (-1)

-- | Cosmic variance squared, usually referred to as sigma^2(R,z)
-- Derived by integrating a matter power spectrum and a window function
cosmicVarianceSq :: ReferenceCosmology -> Mhalo -> Redshift -> W_kind -> Double
cosmicVarianceSq cosmology mh z w_kind =
  let integrand :: Wavenumber -> Double
      integrand k =
        (k ** 2 / (2 * pi ** 2))
          * (powerSpectrum k z)
          * (windowFunction cosmology w_kind k mh) ** 2

      params :: QuadParam
      params = QuadParam {quadPrecision = 1e-9, quadMaxIter = 20}

      result :: Maybe Double
      result = quadRes $ quadSimpson params (1e-4, 1e4) integrand
   in fromMaybe 0.0 result

-- | First-crossing distribution, crucial for the derivation of a
-- Halo Mass Function afterwards, again you have a choice of two:
--    * Sheth-Tormen
--    * Tinker
-- We are planning to add more options in the near future
firstCrossing :: ReferenceCosmology -> HMF_kind -> W_kind -> Mhalo -> Redshift -> Double
firstCrossing cosmology h_kind w_kind mh z =
  let sigma :: Double
      sigma = sqrt $ cosmicVarianceSq cosmology mh z w_kind

      nu :: Double
      nu = 1.686 / sigma

      a_T, a_ST, p, a, b, c :: Double
      (a_T, a_ST, p, a, b, c) = (0.186, 0.3222, 0.3, 1.47, 2.57, 1.19)
   in case h_kind of
        Tinker -> a_T * ((sigma / b) ** (-a) + 1) * exp (-c / sigma ** 2)
        ST -> a_ST * sqrt (2 * nu ** 2 / pi) * (1 + nu ** (-2 * p)) * exp (-nu ** 2 / 2)

-- | The Halo Mass Function (HMF) itself, uses most of the functions
-- defined within this module and a differentiation library
haloMassFunction :: ReferenceCosmology -> HMF_kind -> W_kind -> [Mhalo] -> Redshift -> [Double]
haloMassFunction cosmology h_kind w_kind mh_arr z =
  let (h0, om0, omb0, c, gn) = unpackCosmology cosmology

      diff_func :: Double -> Double
      diff_func mh = log $ sqrt (cosmicVarianceSq cosmology mh z w_kind)
      dsdm = (\mh -> diffRes $ diffRichardson diff_func 1.0 mh) <$> mh_arr

      dsdm :: [Double]
      dsdlogm = zipWith (*) mh_arr dsdm

      rho_mean :: Double
      rho_mean = 3 * h0 ** 2 * om0 / (8 * pi * gn)

      first_crossing :: [Double]
      first_crossing = (\mh -> firstCrossing cosmology h_kind w_kind mh z) <$> mh_arr

      fdsdlogm :: [Double]
      fdsdlogm = zipWith (*) dsdlogm first_crossing
   in zipWith (*) fdsdlogm ((* (-rho_mean)) <$> mh_arr)

main :: IO ()
main = print $ haloMassFunction (MkCosmology {h0' = 69, om0' = 0.305, ob0' = 0.05, c' = 3e5, gn' = 1e-11}) ST TopHat [1e10, 1e11] 0
