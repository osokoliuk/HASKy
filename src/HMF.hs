{-# LANGUAGE OverloadedStrings #-}

module HMF where

{-
Module      : HASKy.HMF
Description : Halo Mass Function
Copyright   : (c) Oleksii Sokoliuk, 20256
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

This module defines a bunch of routines that in the end
yield a Halo Mass Function for a given cosmology (i.e., values of
Hubble parameter H0, Omega_m0, Omega_b0)
-}

import Cosmology
import qualified Data.Map as M
import qualified Data.Map.Strict as M'
import Data.Maybe
import qualified Data.Text as T
import qualified Data.Vector as V
import Helper
import Math.GaussianQuadratureIntegration
import Numeric.Tools.Differentiation
import Numeric.Tools.Integration
import qualified Safe
import System.IO

-- Define an HMF datatype, consisting of five possible choices
data HMF_kind
  = Tinker
  | ST
  | Angulo
  | Jenkins
  | Warren
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

type PowerSpectrum = Redshift -> Wavenumber -> Double

-- | Define a radius for the uniform density sphere in terms of it's mass
rh :: ReferenceCosmology -> W_kind -> Mhalo -> Rhalo
rh cosmology w_kind mh =
  let (h0, om0, ob0, _, gn, _, _) = unpackCosmology cosmology
      rho_mean = 3 * h0 ** 2 * om0 / (8 * pi * gn)
      c_smooth = 3.3
   in case w_kind of
        TopHat -> (3 * mh / (4 * pi * rho_mean)) ** (1 / 3)
        Smooth -> (3 * mh / (4 * pi * rho_mean * c_smooth ** 3)) ** (1 / 3)

-- | Produces a matter power spectrum at the linear level
-- To be imported from CAMB
powerSpectrum :: FilePath -> IO ([Double], [Double])
powerSpectrum filepath =
  do
    content <- readFile filepath
    let linesOfFile = lines content
        parsedLines = map parseLine linesOfFile
        (xValues, yValues) = unzip parsedLines
    return (xValues, yValues)

-- | Window function, used to derive the cosmic variance
-- You have a choice of two different ones, namely:
--    * Top-Hat
--    * Smooth-k
windowFunction :: ReferenceCosmology -> W_kind -> Wavenumber -> Mhalo -> Double
windowFunction cosmology w_kind k mh =
  let r = rh cosmology w_kind mh
      kr = k * r
      beta = 4.8
   in case w_kind of
        TopHat -> 3 * (sin (kr) - kr * cos (kr)) / (kr) ** 3
        Smooth -> (1 + kr ** beta) ** (-1)

-- | Cosmic variance squared, usually referred to as sigma^2(R,z)
-- Derived by integrating a matter power spectrum and a window function
cosmicVarianceSq :: ReferenceCosmology -> PowerSpectrum -> Mhalo -> Redshift -> W_kind -> Double
cosmicVarianceSq cosmology pk mh z w_kind =
  let integrand :: Wavenumber -> Double
      integrand k =
        (k ** 2 / (2 * pi ** 2))
          * (pk z k)
          * (windowFunction cosmology w_kind k mh) ** 2

      result = nIntegrate256 integrand 1e-3 1e3
   in result

-- | First-crossing distribution, crucial for the derivation of a
-- Halo Mass Function afterwards, again you have a choice of two:
--    * Sheth-Tormen
--    * Tinker
-- We are planning to add more options in the near future
firstCrossing :: ReferenceCosmology -> PowerSpectrum -> HMF_kind -> W_kind -> Mhalo -> Redshift -> Double
firstCrossing cosmology pk h_kind w_kind mh z =
  let (h0, om0, ob0, _, _, _, _) = unpackCosmology cosmology
      sigma = sqrt $ cosmicVarianceSq cosmology pk mh z w_kind

      -- Critical linear overdensity threshold with
      -- redshift corrections from [Kitayama et al. 1996]
      delta_c :: Redshift -> Double
      delta_c z = 1.686 * (om0 * (1 + z) ** 3) ** 0.0055
      nu = delta_c z / sigma
      a_T', a_ST', p, a_T, b_T, c_T, a_Ang, b_Ang, c_Ang, a_Jen, b_Jen, a_War, b_War, c_War :: Double
      (a_T', a_ST', p, a_T, b_T, c_T, a_Ang, b_Ang, c_Ang, a_Jen, b_Jen, a_War, b_War, c_War) =
        (0.186, 0.3222, 0.3, 1.47, 2.57, 1.19, 0.201, 2.08, -1.172, 0.315, 0.61, 0.7234, 0.2538, 1.1982)
   in case h_kind of
        Tinker -> a_T' * ((sigma / b_T) ** (-a_T) + 1) * exp (-c_T / sigma ** 2)
        ST -> a_ST' * sqrt (2 * nu ** 2 / pi) * (1 + nu ** (-2 * p)) * exp (-nu ** 2 / 2)
        Angulo -> a_Ang * (b_Ang / sigma + 1) ** 1.7 * exp (c_Ang / sigma ** 2)
        Jenkins -> a_Jen * exp (-abs (log sigma ** (-1) + b_Jen) ** 3.8)
        Warren -> a_War * (sigma ** (-1.625) + b_War) * exp (-c_War / sigma)

-- | The Halo Mass Function (HMF) itself, uses most of the functions
-- defined within this module and a differentiation library
haloMassFunction :: ReferenceCosmology -> PowerSpectrum -> HMF_kind -> W_kind -> Mhalo -> Redshift -> Double
haloMassFunction cosmology pk h_kind w_kind mh z =
  let (h0, om0, ob0, _, gn, _, _) = unpackCosmology cosmology
      rho_mean = 3 * h0 ** 2 * om0 / (8 * pi * gn)

      sigma = \mh -> cosmicVarianceSq cosmology pk mh z w_kind
      first_crossing = \mh -> firstCrossing cosmology pk h_kind w_kind mh z

      diff_func mh = log . sqrt $ sigma mh
      dsdm = diffRes $ diffRichardson diff_func 100 mh
      dsdlogm = dsdm / mh
      fdsdlogm = dsdlogm * first_crossing mh
   in -rho_mean * fdsdlogm * mh

-- | Escape velocity squared of a star from a halo of mass M and radius R at the redshift z,
-- in the units of [km^2 s^-2], taken from the [Tan et al. 2018]
escapeVelocitySq :: ReferenceCosmology -> PowerSpectrum -> HMF_kind -> W_kind -> Mhalo -> Redshift -> Double
escapeVelocitySq cosmology pk h_kind w_kind mh_min z =
  let (h0, om0, ob0, _, gn, _, _) = unpackCosmology cosmology

      mh_arr = (10 **) <$> [5, 5 + 0.25 .. 17]

      hmf_arr =
        (\mh -> haloMassFunction cosmology pk h_kind w_kind mh z) <$> mh_arr
      dndmh = zipWith (/) hmf_arr mh_arr

      (_, n_CDM) =
        (unzip . M'.toList)
          (cumulativeTrapezoidMap $ M'.fromList (zip mh_arr dndmh))
      n_interp = makeInterp mh_arr n_CDM

      integrand_1 mh =
        mh * (2 * gn * mh / rh cosmology w_kind mh) * n_interp mh
      integrand_2 mh =
        mh * n_interp mh

      result =
        (nIntegrate256 integrand_1 mh_min 1e18)
          / (nIntegrate256 integrand_2 mh_min 1e18)
   in result
