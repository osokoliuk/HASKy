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
import qualified Data.Vector as V
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

vecAdd :: (Num a) => V.Vector a -> V.Vector a -> V.Vector a
vecAdd = V.zipWith (+)

vecMultiply :: (Num a) => a -> V.Vector a -> V.Vector a
vecMultiply x = V.map (x *)

rk4Step :: (Double -> V.Vector Double -> V.Vector Double) -> Double -> Double -> V.Vector Double -> V.Vector Double
rk4Step f t h y =
  vecAdd y $
    vecMultiply (1 / 6) $
      vecAdd k1 $
        vecAdd (vecMultiply 2 k2) $
          vecAdd (vecMultiply 2 k3) k4
  where
    k1 = vecMultiply h (f t y)
    k2 = vecMultiply h (f (t + h / 2) (vecAdd y (vecMultiply 0.5 k1)))
    k3 = vecMultiply h (f (t + h / 2) (vecAdd y (vecMultiply 0.5 k2)))
    k4 = vecMultiply h (f (t + h) (vecAdd y k3))

rk4Solve :: (Double -> V.Vector Double -> V.Vector Double) -> Double -> Double -> Int -> V.Vector Double -> [(Double, V.Vector Double)]
rk4Solve f t0 h n y0 = take (n + 1) $ iterate step (t0, y0)
  where
    step (t, y) =
      let y' = rk4Step f t h y
       in (t + h, y')

main :: IO ()
main = do
  let y0 = V.fromList [1.0, 0.0] -- Initial state [position, velocity]
      t0 = 0.0
      h = 0.1
      n = 100
      oscillator t y = V.fromList [y V.! 1, -(y V.! 0)]
      result = rk4Solve oscillator t0 h n y0
  mapM_ print result
