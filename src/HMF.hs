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
import Data.Maybe
import qualified Data.Text as T
import qualified Data.Vector as V
import Helper
import Math.GaussianQuadratureIntegration
import Numeric.Tools.Differentiation
import Numeric.Tools.Integration
import qualified Safe
import System.IO

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

type PowerSpectrum = Wavenumber -> Double

-- | Define a radius for the uniform density sphere in terms of it's mass
rh :: ReferenceCosmology -> W_kind -> Mhalo -> Rhalo
rh cosmology w_kind mh =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology
      rho_mean :: Double
      rho_mean = 3 * h0 ** 2 * om0 / (8 * pi * gn)

      c_smooth :: Double
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
cosmicVarianceSq :: FilePath -> ReferenceCosmology -> PowerSpectrum -> Mhalo -> Redshift -> W_kind -> IO Double
cosmicVarianceSq filepath cosmology pk mh z w_kind =
  do
    let integrand :: Wavenumber -> Double
        integrand k =
          (k ** 2 / (2 * pi ** 2))
            * (pk k)
            * (windowFunction cosmology w_kind k mh) ** 2

        result :: Double
        result = nIntegrate1024 integrand 1e-4 1e4
    return $ result

-- | First-crossing distribution, crucial for the derivation of a
-- Halo Mass Function afterwards, again you have a choice of two:
--    * Sheth-Tormen
--    * Tinker
-- We are planning to add more options in the near future
firstCrossing :: FilePath -> ReferenceCosmology -> PowerSpectrum -> HMF_kind -> W_kind -> Mhalo -> Redshift -> IO Double
firstCrossing filepath cosmology pk h_kind w_kind mh z =
  do
    sigma <- cosmicVarianceSq filepath cosmology pk mh z w_kind

    let (h0, om0, ob0, c, gn) = unpackCosmology cosmology

        sigma_eval :: Double
        sigma_eval = sqrt sigma

        -- Critical linear overdensity threshold with
        -- redshift corrections from [Kitayama et al. 1996]
        delta_c :: Redshift -> Double
        delta_c z = 1.686 * (om0 * (1 + z) ** 3) ** 0.0055

        nu :: Double
        nu = delta_c z / sigma_eval -- CHANGE DELTA_c!!!!
        a_T', a_ST', p, a_T, b_T, c_T :: Double
        (a_T', a_ST', p, a_T, b_T, c_T) = (0.186, 0.3222, 0.3, 1.47, 2.57, 1.19)
    return $ case h_kind of
      Tinker -> a_T' * ((sigma_eval / b_T) ** (-a_T) + 1) * exp (-c_T / sigma_eval ** 2)
      ST -> a_ST' * sqrt (2 * nu ** 2 / pi) * (1 + nu ** (-2 * p)) * exp (-nu ** 2 / 2)

-- | The Halo Mass Function (HMF) itself, uses most of the functions
-- defined within this module and a differentiation library
haloMassFunction :: FilePath -> ReferenceCosmology -> HMF_kind -> W_kind -> [Mhalo] -> Redshift -> IO [Double]
haloMassFunction filepath cosmology h_kind w_kind mh_arr z =
  do
    -- Map arrays for power spectrum, variance and first crossing distributions
    (k_arr, pk_arr) <- powerSpectrum filepath

    let interp_pk :: PowerSpectrum
        interp_pk = mapLookup (M.fromList (zip k_arr pk_arr))

    sigma_arr <- mapM (\mh -> cosmicVarianceSq filepath cosmology interp_pk mh z w_kind) mh_arr
    first_crossing_arr <- mapM (\mh -> firstCrossing filepath cosmology interp_pk h_kind w_kind mh z) mh_arr

    let (h0, om0, ob0, c, gn) = unpackCosmology cosmology

        interp_sigma :: Double -> Double -- Interpolate cosmic variance
        interp_sigma = mapLookup (M.fromList (zip mh_arr sigma_arr))

        diff_func :: Double -> Double
        diff_func mh = log . sqrt $ interp_sigma mh
        dsdm = (\mh -> diffRes $ diffRichardson diff_func 1000 mh) <$> mh_arr

        dsdlogm :: [Double]
        dsdlogm = zipWith (/) dsdm mh_arr

        rho_mean :: Double
        rho_mean = 3 * h0 ** 2 * om0 / (8 * pi * gn)

        fdsdlogm :: [Double]
        fdsdlogm = zipWith (*) dsdlogm first_crossing_arr
    return $ zipWith (*) fdsdlogm ((* (-rho_mean)) <$> mh_arr)

main_HMF :: IO ()
main_HMF = do
  x <- haloMassFunction "../data/CAMB_Pk_z=0.txt" planck18 ST Smooth (map (\x -> 10 ** x) [9, 9 + 0.2 .. 15]) 0
  print $ x
