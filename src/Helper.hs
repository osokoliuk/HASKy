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

import Control.Monad (unless)
import Data.Bifunctor
import Data.Char (isDigit, isSpace, toLower, toUpper)
import Data.List (dropWhileEnd, transpose)
import qualified Data.Map as M
import qualified Data.Vector as V
import System.IO
import Text.Read (readMaybe)

data Element
  = Element
  { element :: String,
    isotope :: Int
  }
  deriving (Eq, Ord)

-- Create a read instance for the Element type so that it can ignore first uppercase letter
instance Read Element where
  readsPrec _ str =
    let (x : xs, isotope) = span (not . isDigit) str
     in case read isotope of
          i -> [(Element (toUpper x : xs) i, "")]
          _ -> []

instance Show Element where
  show (Element elem isotope) = (toLower <$> elem) ++ show isotope

data Table = Table
  { masses :: [Double],
    values :: [(Element, [Double])]
  }
  deriving (Show)

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

parseFileToTable :: FilePath -> IO Table
parseFileToTable path = do
  content <- readFile path
  let ls = lines content
  case ls of
    (_ : headerLine : _ : rows) -> do
      let headerWords = words headerLine
      case headerWords of
        ("#" : "M_init" : elems) -> do
          let tableRows = map words rows
              cols = transpose tableRows
          unless (length cols >= 1 + length elems) $
            error "Not enough columns for all elements"
          let massCol = map read (head cols)
              elemCols = map (map read) (tail cols)
              namedCols = zip (read <$> elems :: [Element]) elemCols
          return $ Table massCol namedCols
        _ -> error "Invalid header format"
    _ -> error "File too short"

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
  rows <- parseFileToTable "../data/test.txt"
  print rows
