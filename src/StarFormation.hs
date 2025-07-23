module StarFormation where

{-
Module      : HASKy.StarFormation
Description : StarFormation module
Copyright   : (c) Oleksii Sokoliuk, 2025
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

A module that defines reference star formation model. Pretty much used by
every other module in this library.
-}

import Helper
import Math.GaussianQuadratureIntegration
import System.Environment

-- | Define a star formation model datatype, which includes the following:
--    * yields_ia -> Element/Isotope yields model for SN type Ia, in [Msol]
--    * yields_ccsn -> Element/Isotope yields model for CCSN, in [Msol]
--    * yields_hne -> Element/Isotope yields model for HNe, in [Msol]
--    * yields_agb -> Element/Isotope yields model for AGB stars, in [MSol]
--    * yields_NSM  -> Element/Isotope yields model for NS-NS and NS-BH mergers, in [Msol]
data ReferenceStarFormationModel
  = MkStarFormation
  { model_ia :: [Char],
    model_ccsn :: [Char],
    model_hne :: [Char],
    model_agb :: [Char],
    model_sagb :: [Char],
    model_nsm :: [Char]
  }
  deriving (Eq, Show, Ord)

-- Functions to extract yields from Yield datatype
retrieveYieldIa :: Yield -> Element -> IO Double
retrieveYieldIa yield elem =
  let filepath = "data/Ia/" ++ show (model_ia yield) ++ "/" ++ (toLower <$> element elem) ++ ".dat"
   in do
        table <- parseFile_Ia filepath (toLower <$> show elem)
        return table

retrieveYieldCCSN :: Yield -> Element -> Double -> IO ([Double], [Double])
retrieveYieldCCSN yield elem metal_frac =
  let metal_str
        | metal_frac <= 0.001 = "z0001"
        | metal_frac <= 0.01 = "z001"
        | metal_frac <= 0.1 = "z01"
        | otherwise = "z1"
   in let filepath = "data/CCSNe/" ++ show (model_ccsn yield) ++ "/" ++ metal_str ++ "/" ++ (toLower <$> element elem) ++ ".dat"
       in do
            table <- parseFile_II filepath
            let yields_arr = lookup elem (values table)
            return $ (masses table, fromMaybe [] yields_arr)

retrieveYieldAGB :: Yield -> Element -> Double -> IO ([Double], [Double])
retrieveYieldAGB yield elem metal_frac =
  let metal_str
        | metal_frac <= 0.001 = "z0001"
        | metal_frac <= 0.01 = "z001"
        | metal_frac <= 0.1 = "z01"
        | otherwise = "z1"
   in let filepath = "data/AGB/" ++ show (model_agb yield) ++ "/" ++ metal_str ++ "/" ++ (toLower <$> element elem) ++ ".dat"
       in do
            table <- parseFile_II filepath
            let yields_arr = lookup elem (values table)
                yields_corrected = [if x >= 0 then x else 0 | x <- fromMaybe [] yields_arr]
            return $ (masses table, yields_corrected)

retrieveYieldSAGB :: Yield -> Element -> Double -> IO ([Double], [Double])
retrieveYieldSAGB yield elem metal_frac =
  let metal_str
        | metal_frac <= 0.001 = "z0001"
        | metal_frac <= 0.01 = "z001"
        | metal_frac <= 0.1 = "z01"
        | otherwise = "z1"
   in let filepath = "data/SAGB/" ++ show (model_agb yield) ++ "/" ++ metal_str ++ "/" ++ (toLower <$> element elem) ++ ".dat"
       in do
            table <- parseFile_II filepath
            let yields_arr = lookup elem (values table)
                yields_corrected = [if x >= 0 then x else 0 | x <- fromMaybe [] yields_arr]
            return $ (masses table, yields_corrected)

retrieveYieldNSM :: Yield -> Element -> IO Double
retrieveYieldNSM yield elem =
  pure 1.0

retrieveYieldHNe :: Yield -> Element -> IO Double
retrieveYieldHNe yield elem =
  pure 1.0
