{-# LANGUAGE OverloadedStrings #-}

module HMF where

import Cosmology
import qualified Data.Map as M
import Data.Maybe
import qualified Data.Text as T
import qualified Data.Vector as V
import Numeric.Tools.Differentiation
import Numeric.Tools.Integration
import Numeric.Tools.Interpolation
import qualified Safe
import System.IO
import Text.Read (readMaybe)

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

-- | Helper function to linearly interpolate power spectrum
-- Taken from the https://cmears.id.au/articles/linear-interpolation.html
interpolate (a, av) (b, bv) x = av + (x - a) * (bv - av) / (b - a)

mapLookup :: M.Map Double Double -> Double -> Double
mapLookup m x =
  case (M.lookupLE x m, M.lookupGE x m) of
    (Just (a, av), Just (b, bv)) ->
      if a == b
        then av
        else interpolate (a, av) (b, bv) x
    (Nothing, Just (b, bv)) -> bv
    (Just (a, av), Nothing) -> av
    _ -> error "mapLookup"

-- | Helper function for reading the file into a table
parseLine :: String -> (Double, Double)
parseLine line =
  case mapM readMaybe (words line) of -- Try to read both values
    Just [x, y] -> (x, y) -- If both are parsed, return the tuple
    _ -> error ("Invalid line: " ++ line)

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
  let r = rh cosmology mh
      kr = k * r
      beta = 4.8
   in case w_kind of
        TopHat -> 3 * (sin (kr) - kr * cos (kr)) / (kr) ** 3
        Smooth -> (1 + kr ** beta) ** (-1)

-- | Cosmic variance squared, usually referred to as sigma^2(R,z)
-- Derived by integrating a matter power spectrum and a window function
cosmicVarianceSq :: FilePath -> ReferenceCosmology -> Mhalo -> Redshift -> W_kind -> IO Double
cosmicVarianceSq filepath cosmology mh z w_kind =
  do
    (k_arr, pk_arr) <- powerSpectrum filepath
    let interp_pk :: Double -> Double
        interp_pk k = mapLookup (M.fromList (zip k_arr pk_arr)) k

        integrand :: Wavenumber -> Double
        integrand k =
          (k ** 2 / (2 * pi ** 2))
            * (interp_pk k)
            * (windowFunction cosmology w_kind k mh) ** 2

        params :: QuadParam
        params = QuadParam {quadPrecision = 1e-9, quadMaxIter = 1000}

        result :: Maybe Double
        result = quadRes $ quadSimpson params (1e-4, 1e4) integrand
    return $ fromMaybe 0.0 result

-- | First-crossing distribution, crucial for the derivation of a
-- Halo Mass Function afterwards, again you have a choice of two:
--    * Sheth-Tormen
--    * Tinker
-- We are planning to add more options in the near future
firstCrossing :: FilePath -> ReferenceCosmology -> HMF_kind -> W_kind -> Mhalo -> Redshift -> IO Double
firstCrossing filepath cosmology h_kind w_kind mh z =
  do
    sigma <- cosmicVarianceSq filepath cosmology mh z w_kind

    let sigma_eval :: Double
        sigma_eval = sqrt sigma

        nu :: Double
        nu = 1.686 / sigma_eval

        a_T, a_ST, p, a, b, c :: Double
        (a_T, a_ST, p, a, b, c) = (0.186, 0.3222, 0.3, 1.47, 2.57, 1.19)
    return $ case h_kind of
      Tinker -> a_T * ((sigma_eval / b) ** (-a) + 1) * exp (-c / sigma_eval ** 2)
      ST -> a_ST * sqrt (2 * nu ** 2 / pi) * (1 + nu ** (-2 * p)) * exp (-nu ** 2 / 2)

-- | The Halo Mass Function (HMF) itself, uses most of the functions
-- defined within this module and a differentiation library
haloMassFunction :: FilePath -> ReferenceCosmology -> HMF_kind -> W_kind -> [Mhalo] -> Redshift -> IO [Double]
haloMassFunction filepath cosmology h_kind w_kind mh_arr z =
  do
    sigma_arr <- mapM (\mh -> cosmicVarianceSq filepath cosmology mh z w_kind) mh_arr
    first_crossing_arr <- mapM (\mh -> firstCrossing filepath cosmology h_kind w_kind mh z) mh_arr

    let (h0, om0, omb0, c, gn) = unpackCosmology cosmology

        interp_sigma :: Double -> Double
        interp_sigma mh = mapLookup (M.fromList (zip mh_arr sigma_arr)) mh

        diff_func :: Double -> Double
        diff_func mh = log . sqrt $ interp_sigma mh
        dsdm = (\mh -> diffRes $ diffRichardson diff_func 1.0 mh) <$> mh_arr

        dsdm :: [Double]
        dsdlogm = zipWith (*) mh_arr dsdm

        rho_mean :: Double
        rho_mean = 3 * h0 ** 2 * om0 / (8 * pi * gn)

        fdsdlogm :: [Double]
        fdsdlogm = zipWith (*) dsdlogm first_crossing_arr
    return $ zipWith (*) fdsdlogm ((* (-rho_mean)) <$> mh_arr)

main :: IO ()
main = do
  x <- haloMassFunction "camb.txt" (MkCosmology {h0' = 69, om0' = 0.305, ob0' = 0.05, c' = 3e5, gn' = 1e-11}) ST TopHat [1e10, 1e11] 0
  print x
