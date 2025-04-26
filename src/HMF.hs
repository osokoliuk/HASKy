module HMF where

import Cosmology
import Data.Maybe
import Numeric.Tools.Differentiation
import Numeric.Tools.Integration

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

-- Initialise background cosmology from the arguments, taken from Cosmology.hs
cosmology_record = initialiseCosmology []

-- Unpack the values in the cosmology_record into variables
h0, om0, c, gn :: Double
(h0, om0, c, gn) = (,,,) <$> h0' <*> om0' <*> c' <*> gn' $ cosmology_record

-- | Define a radius for the uniform density sphere in terms of it's mass
rh :: Mhalo -> Rhalo
rh mh = 1 / c * (2 * mh * gn / (om0 * h0 ** 2)) ** (1 / 3)

-- | Produces a matter power spectrum at the linear level
-- To be imported from CAMB
powerSpectrum :: Wavenumber -> Redshift -> Double
powerSpectrum k z = 1e0

-- | Window function, used to derive the cosmic variance
-- You have a choice of two different ones, namely:
--    * Top-Hat
--    * Smooth-k
windowFunction :: W_kind -> Wavenumber -> Mhalo -> Double
windowFunction w_kind k mh =
  let r = rh mh
      kr = k * r
      beta = 4.8
   in case w_kind of
        TopHat -> 3 * (sin (kr) - kr * cos (kr)) / (kr) ** 3
        Smooth -> (1 + kr ** beta) ** (-1)

-- | Cosmic variance squared, usually referred to as sigma^2(R,z)
-- Derived by integrating a matter power spectrum and a window function
cosmicVarianceSq :: Mhalo -> Redshift -> W_kind -> Double
cosmicVarianceSq mh z w_kind =
  let integrand :: Wavenumber -> Double
      integrand k =
        (k ** 2 / (2 * pi ** 2))
          * (powerSpectrum k z)
          * (windowFunction w_kind k mh) ** 2

      params :: QuadParam
      params = QuadParam {quadPrecision = 1e-9, quadMaxIter = 20}

      result :: Maybe Double
      result = quadRes $ quadSimpson params (1e-4, 1e4) sin
   in fromMaybe 0.0 result

-- | First-crossing distribution, crucial for the derivation of a
-- Halo Mass Function afterwards, again you have a choice of two:
--    * Sheth-Tormen
--    * Tinker
-- We are planning to add more options in the near future
firstCrossing :: HMF_kind -> W_kind -> Mhalo -> Redshift -> Double
firstCrossing h_kind w_kind mh z =
  let sigma :: Double
      sigma = sqrt $ cosmicVarianceSq mh z w_kind

      nu :: Double
      nu = 1.686 / sigma

      a_T, a_ST, p, a, b, c :: Double
      (a_T, a_ST, p, a, b, c) = (0.186, 0.3222, 0.3, 1.47, 2.57, 1.19)
   in case h_kind of
        Tinker -> a_T * ((sigma / b) ** (-a) + 1) * exp (-c / sigma ** 2)
        ST -> a_ST * sqrt (2 * nu ** 2 / pi) * (1 + nu ** (-2 * p)) * exp (-nu ** 2 / 2)

-- | The Halo Mass Function (HMF) itself, uses most of the functions
-- defined within this module and a differentiation library
haloMassFunction :: HMF_kind -> W_kind -> [Mhalo] -> Redshift -> [Double]
haloMassFunction h_kind w_kind mh_arr z =
  let diff_func :: Double -> Double
      diff_func mh = log $ sqrt (cosmicVarianceSq mh z w_kind)
      dsdm = (\mh -> diffRes $ diffRichardson diff_func 1.0 mh) <$> mh_arr

      dsdm :: [Double]
      dsdlogm = zipWith (*) mh_arr dsdm

      rho_mean :: Double
      rho_mean = 3 * h0 ** 2 * om0 / (8 * pi * gn)

      first_crossing :: [Double]
      first_crossing = (\mh -> firstCrossing h_kind w_kind mh z) <$> mh_arr

      fdsdlogm :: [Double]
      fdsdlogm = zipWith (*) dsdlogm first_crossing
   in zipWith (*) fdsdlogm ((* (-rho_mean)) <$> mh_arr)
