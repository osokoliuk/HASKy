{-# LANGUAGE OverloadedStrings #-}

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

import Data.Char (toLower)
import Data.Maybe (fromMaybe)
import Helper
import Math.GaussianQuadratureIntegration
import System.Environment

-- | Define a star formation model datatype, which includes the following:
--    * yield_ia -> Element/Isotope yields model for SN type Ia, in [Msol]
--    * yield_ccsn -> Element/Isotope yields model for CCSN, in [Msol]
--    * yield_hne -> Element/Isotope yields model for HNe, in [Msol]
--    * yield_ecsn -> Element/Isotope yields model for ECSN, in [Msol]
--    * yield_agb -> Element/Isotope yields model for AGB stars, in [MSol]
--    * yield_NSM  -> Element/Isotope yields model for NS-NS and NS-BH mergers, in [Msol]
--    * yield_PSN -> Element/Isotope yields model for proto-NS neutrino-drived yields, in [Msol]
--    * yield_CO -> Element/Isotope yields model for CO WDs within a Nova system, in [Msol]
--    * yield_ONe -> Element/Isotope yields model for ONe WDs within a Nova system, in [Msol]
data ReferenceStarFormationModel
  = MkStarFormation
  { yield_ia :: Double,
    yield_ccsn :: ([Double], [Double]),
    yield_hne :: ([Double], [Double]),
    yield_ecsn :: ([Double], [Double]),
    yield_agb :: ([Double], [Double]),
    yield_sagb :: ([Double], [Double]),
    yield_nsm :: Double,
    yield_psn :: ([Double], [Double]),
    yield_co :: Double,
    yield_one :: Double
  }

data ReferenceStarFormationConfig
  = MkStarFormationCfg
  { model_ia :: [Char],
    model_ccsn :: [Char]
  }

-- Functions to extract yields from Yield datatype
retrieveYieldIa :: ReferenceStarFormationConfig -> Element -> IO Double
retrieveYieldIa cfg elem =
  let filepath = "data/" ++ model_ia cfg ++ "/" ++ (toLower <$> element elem) ++ ".dat"
   in do
        table <- parseFile_Ia filepath (toLower <$> show elem)
        return table

retrieveYieldCCSN :: ReferenceStarFormationConfig -> Element -> Double -> IO ([Double], [Double])
retrieveYieldCCSN cfg elem metal_frac =
  let metal_str
        | metal_frac <= 0.001 = "z0001"
        | metal_frac <= 0.01 = "z001"
        | metal_frac <= 0.1 = "z01"
        | otherwise = "z1"
   in let filepath = "data/" ++ model_ccsn cfg ++ "/" ++ metal_str ++ "/" ++ (toLower <$> element elem) ++ ".dat"
       in do
            table <- parseFile_II filepath
            let yields_arr = lookup elem (values table)
            return $ (masses table, fromMaybe [] yields_arr)

{-
retrieveYieldAGB :: ReferenceStarFormationModel -> Element -> Double -> IO ([Double], [Double])
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

retrieveYieldSAGB :: ReferenceStarFormationModel -> Element -> Double -> IO ([Double], [Double])
retrieveYieldSAGB yield elem metal_frac =
  let metal_str
        | metal_frac <= 0.001 = "z0001"
        | metal_frac <= 0.01 = "z001"
        | metal_frac <= 0.1 = "z01"
        | otherwise = "z1"
   in let filepath = "data/SAGB/" ++ show (model_sagb yield) ++ "/" ++ metal_str ++ "/" ++ (toLower <$> element elem) ++ ".dat"
       in do
            table <- parseFile_II filepath
            let yields_arr = lookup elem (values table)
                yields_corrected = [if x >= 0 then x else 0 | x <- fromMaybe [] yields_arr]
            return $ (masses table, yields_corrected)

retrieveYieldNSM :: ReferenceStarFormationModel -> Element -> IO Double
retrieveYieldNSM yield elem =
  pure 1.0

retrieveYieldHNe :: ReferenceStarFormationModel -> Element -> IO Double
retrieveYieldHNe yield elem =
  pure 1.0

-}
