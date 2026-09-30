{-# LANGUAGE RecordWildCards #-}

module Main where

import Control.Concurrent.Async (mapConcurrently)
import Control.Lens
import Control.Monad (when)
import Control.Parallel.Strategies
import Cosmology
import DSP.Basic (logspace)
import Data.Char (toLower)
import Data.Colour
import Data.Colour.Names
import Data.Default.Class
import Data.List
import qualified Data.Vector as V
import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Cairo
import Graphics.Rendering.Chart.Easy
import Graphs
import HMF
import Helper
import IGM
import Lookup
import Pk
import SMF
import StarFormation
import System.Directory (
    doesFileExist,
    removeFile,
 )

args :: [Double]
args = []

main :: IO ()
main = do
    -- Fix this interpolation later, this code is not prod ready...
    let pk = powerSpectrumEisensteinHu planck18
        elem = ["H", "Fe", "Si"]
        sfCfg = MkStarFormationCfg{model_ia = "iwamoto99/WDD1", model_ccsn = "WW95", model_agb = "Cristallo11", model_ecsn = "Wanajo13", model_hne = "Kobayashi06"}
        IGMParams{..} = defaultIGMParams
        PhysicalConstants{..} = phys
        MassLimits{..} = masses
        Efficiencies{..} = effs
        DelayTimes{..} = delays
        NovaeParams{..} = novae
        ECSNParams{..} = ecsn

    -- parsed <- parseFileColumns "data/CCSN/WW95/z0001/al.dat"
    -- print parsed

    (times, redshifts, chosenIsotopes, groupedChosenIsotopes, abundances, abundancesIsotopes) <- igmIsmEvolution sfCfg planck18 pk Pereira Kroupa DoublePower Tinker Smooth Constant_HNe elem mHaloMin

    plotMassFractions times abundances "igmism.png"
    plotIsotopeAbundances times abundancesIsotopes elem groupedChosenIsotopes "abundances.png"
    plotMetallicity times abundances "metallicities.png"
