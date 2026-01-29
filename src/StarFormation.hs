{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

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
    model_ecsn :: [Char],
    model_hne :: [Char]
  }

-- | Define a datatype containing all of the tunable parameters used in the IGM.hs
--    * PhysicalConstants - some of the unit conversions
--    * MassLimits - limits of integrations of IMF for different sources of nucleosynthesis yields
--    * Efficiencies - fractions of things related to star formation, such as the fraction of CCSN that are HNe
--    * DelayTimes - time delays in the NSM and Novae systems
--    * NovaeParams - tunable parameters for Novae systems
--    * ECSNParams - tunable parameters for ECSNe
--    * IGMParams - all of the aforementioned parameters compiled in a single datatype
data PhysicalConstants = PhysicalConstants
  { zMax :: Double,
    mCO :: Double,
    yrGyr :: Double,
    kmsErgMsol :: Double,
    snEnergy :: Double
  }

data MassLimits = MassLimits
  { mUp :: Double,
    mPU :: Double,
    mPL :: Double,
    mDURG :: Double,
    mDLRG :: Double,
    mDUMS :: Double,
    mDLMS :: Double,
    mNSMd :: Double,
    mNSMu :: Double,
    mNovaeD :: Double,
    mNovaeU :: Double,
    mAGBd :: Double,
    mAGBu :: Double,
    mECSNd :: Double,
    mECSNu :: Double
  }

data Efficiencies = Efficiencies
  { epsW :: Double,
    epsSN :: Double,
    bRG :: Double,
    bMS :: Double,
    epsHNe0 :: Double,
    epsMRSNe :: Double,
    alphaNSM :: Double,
    alphaNovae :: Double,
    epsCO :: Double
  }

data DelayTimes = DelayTimes
  { delayNSM :: Double,
    delayNovae :: Double
  }

data NovaeParams = NovaeParams
  { mEjNovae :: Double,
    nNovae :: Double
  }

data ECSNParams = ECSNParams
  { mEjECSN :: Double
  }

data IGMParams = IGMParams
  { phys :: PhysicalConstants,
    masses :: MassLimits,
    effs :: Efficiencies,
    delays :: DelayTimes,
    novae :: NovaeParams,
    ecsn :: ECSNParams
  }

defaultIGMParams :: IGMParams
defaultIGMParams =
  IGMParams
    { phys =
        PhysicalConstants
          { zMax = 20,
            mCO = 1.38,
            yrGyr = 1e9,
            kmsErgMsol = 1.989e43,
            snEnergy = 2e51
          },
      masses =
        MassLimits
          { mUp = 100,
            mPU = 8.0,
            mPL = 3.0,
            mDURG = 1.5,
            mDLRG = 0.9,
            mDUMS = 2.6,
            mDLMS = 1.8,
            mNSMd = 9.0,
            mNSMu = 30.0,
            mNovaeD = 0.8,
            mNovaeU = 8.0,
            mAGBd = 1.3,
            mAGBu = 8.0,
            mECSNd = 8.8,
            mECSNu = 9.0
          },
      effs =
        Efficiencies
          { epsW = 0.02,
            epsSN = 0.005,
            bRG = 0.02,
            bMS = 0.04,
            epsHNe0 = 0.5,
            epsMRSNe = 0.03,
            alphaNSM = 0.018,
            alphaNovae = 0.01,
            epsCO = 0.7
          },
      delays =
        DelayTimes
          { delayNSM = 1e7,
            delayNovae = 2e9
          },
      novae =
        NovaeParams
          { mEjNovae = 2e-5,
            nNovae = 1e4
          },
      ecsn =
        ECSNParams
          { mEjECSN = 1.14e-2
          }
    }

-- Functions to extract yields from ReferenceStarFormationConfig datatype
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
            return $ (indices table, fromMaybe [] yields_arr)

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

retrieveYieldHNe :: ReferenceStarFormationConfig -> Element -> Double -> IO ([Double], [Double])
retrieveYieldHNe cfg isotopeElem metalFrac =
  -- Create a map lookup for a metal fraction/model
  let modelMetallicity :: M.Map String [Double]
      modelMetallicity =
        M.fromList $
          [ ("Kobayashi06", [0.0, 0.02, 0.004, 0.001])
          ]
      metalStr =
        let lookupList = fromMaybe [] $ M.lookup (model_hne cfg) modelMetallicity
            idx = fromMaybe 0 $ findClosestList metalFrac lookupList
         in (\x -> if x == '.' then 'p' else x) <$> showFFloat Nothing (lookupList !! idx) ""
   in let filepath = "data/HNe/" <> model_hne cfg <> "/" <> "yields_z" <> metalStr <> ".dat"
       in do
            yieldTable <- parseFile_HNe filepath
            let yields = fromMaybe (replicate 4 0.0) $ M.lookup (show isotopeElem) yieldTable
            return ([20, 25, 30, 40], yields)
