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
import qualified Data.Map as M
import Data.Maybe (fromMaybe)
import Helper
import Math.GaussianQuadratureIntegration
import Numeric (showFFloat)
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
    yield_ecsn :: Double,
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
    model_ccsn :: [Char],
    model_agb :: [Char],
    model_ecsn :: [Char]
  }

-- Functions to extract yields from Yield datatype
retrieveYieldIa :: ReferenceStarFormationConfig -> Element -> IO Double
retrieveYieldIa cfg isotopeElem =
  let filepath = "data/SN_Ia/" <> model_ia cfg <> "/" <> (toLower <$> element isotopeElem) <> ".dat"
   in do
        table <- parseFile_Ia filepath (toLower <$> show isotopeElem)
        return table

retrieveYieldCCSN :: ReferenceStarFormationConfig -> Element -> Double -> IO ([Double], [Double])
retrieveYieldCCSN cfg isotopeElem metalFrac =
  let metalStr
        | metalFrac <= 0.001 = "z0001"
        | metalFrac <= 0.01 = "z001"
        | metalFrac <= 0.1 = "z01"
        | otherwise = "z1"
   in let filepath = "data/CCSN/" <> model_ccsn cfg <> "/" <> metalStr <> "/" <> (toLower <$> element isotopeElem) <> ".dat"
       in do
            table <- parseFile_II filepath
            let yields_arr = lookup isotopeElem (values table)
            return $ (masses table, fromMaybe [] yields_arr)

retrieveYieldAGB :: ReferenceStarFormationConfig -> Element -> Double -> IO ([Double], [Double])
retrieveYieldAGB cfg isotopeElem metalFrac =
  -- Create a map lookup for a metal fraction/model
  let modelMetallicity :: M.Map String [Double]
      modelMetallicity =
        M.fromList $
          [ ("Cristallo11", [0.01, 0.001, 0.0001, 0.02, 0.002, 0.003, 0.0003, 0.006, 0.008, 0.014]),
            ("Karakas10", [0.0001, 0.02, 0.004, 0.008]),
            ("Karakas16", [0.03, 0.007, 0.014, 0.0028]),
            ("Ventura13", [0.001, 0.002, 0.0003, 0.04, 0.004, 0.008, 0.0014])
          ]
      metalStr =
        let lookupList = fromMaybe [] $ M.lookup (model_agb cfg) modelMetallicity
            idx = fromMaybe 0 $ findClosestList metalFrac lookupList
         in showFFloat Nothing (lookupList !! idx) ""
   in let filepath = "data/AGB/" <> model_agb cfg <> "/" <> metalStr <> "/" <> (toLower <$> element isotopeElem) <> ".dat"
       in do
            parseFile_AGB filepath

retrieveYieldECSN :: ReferenceStarFormationConfig -> Element -> IO Double
retrieveYieldECSN cfg isotopeElem =
  let filepath = "data/ECSN/" <> model_ecsn cfg <> "/" <> "yields.dat"
   in do
        yieldTable <- parseFile_ECSN filepath
        let yield = fromMaybe 0.0 $ M.lookup (show isotopeElem) yieldTable
        return yield
