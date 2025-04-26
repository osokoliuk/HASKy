module HMF where

import Cosmology
import Data.Maybe
import Numeric.Tools.Integration

data HMF_kind
  = Tinker
  | ST
  deriving (Eq, Show)

data W_kind
  = TopHat
  | Smooth
  deriving (Eq, Show)

type Wavenumber = Double

type Mhalo = Double

type Rhalo = Double

type Redshift = Float

cosmology_record = initialiseCosmology []

h0 = h0' cosmology_record

om0 = om0' cosmology_record

c = c' cosmology_record

gn = gn' cosmology_record

rh :: Mhalo -> Rhalo
rh mh = 1 / c * (2 * mh * gn / (om0 * h0 ** 2)) ** (1 / 3)

powerSpectrum :: Wavenumber -> Redshift -> Double
powerSpectrum k z = 1e0

windowFunction :: W_kind -> Wavenumber -> Mhalo -> Double
windowFunction w_kind k mh =
  let r = rh mh
      kr = k * r
      beta = 4.8 -- Best-fit from the N-body simulation [Leo et al. 2018]
   in case w_kind of
        TopHat -> 3 * (sin (kr) - kr * cos (kr)) / (kr) ** 3
        Smooth -> (1 + kr ** beta) ** (-1)

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
