module Helper where

{-
Module      : HASKy.Helper
Description : Helper module
Copyright   : (c) Oleksii Sokoliuk, 20256
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

A module that defines lots of functions to be used by other modules.
Kind of useless by itself.
-}

import qualified Data.Map as M
import Text.Read (readMaybe)

-- | We define a log10 function (which is apparenly absent from the Prelude)
-- just for our convenience
log10 :: (Floating a) => a -> a
log10 x = log x / log 10

-- | Helper function to linearly interpolate power spectrum
-- Taken from the https://cmears.id.au/articles/linear-interpolation.html
extrapolate x [y1, y2] [x1, x2] = y1 + (x - x1) / (x2 - x1) * (y2 - y1)

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

-- | Runge-Kutta solver of 4th order
rungeKutta4 :: (Double -> Double -> Double) -> (Double, Double) -> Double -> (Double, Double)
rungeKutta4 f (t, y) t' = (t', y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0)
  where
    h = t' - t
    k1 = f t y
    k2 = f (t + 0.5 * h) (y + 0.5 * h * k1)
    k3 = f (t + 0.5 * h) (y + 0.5 * h * k2)
    k4 = f (t + 1.0 * h) (y + 1.0 * h * k3)
