{-# LANGUAGE RecordWildCards #-}

module Main where

import Control.Concurrent.Async (mapConcurrently)
import Control.Lens
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
import System.Directory (doesFileExist)

args :: [Double]
args = []

main :: IO ()
main = do
  -- Fix this interpolation later, this code is not prod ready...
  let pk = powerSpectrumEisensteinHu planck18
      elem = [Element {element = "Fe", isotope = 56}, Element {element = "Si", isotope = 28}, Element {element = "C", isotope = 12}, Element {element = "O", isotope = 16}]
      sfCfg = MkStarFormationCfg {model_ia = "iwamoto99/WDD1", model_ccsn = "WW95", model_agb = "Cristallo11", model_ecsn = "Wanajo13", model_hne = "Kobayashi06"}
      IGMParams {..} = defaultIGMParams
      PhysicalConstants {..} = phys
      MassLimits {..} = masses
      Efficiencies {..} = effs
      DelayTimes {..} = delays
      NovaeParams {..} = novae
      ECSNParams {..} = ecsn

  (times, redshifts, abundances) <- igmIsmEvolution sfCfg planck18 pk Pereira Kroupa DoublePower Tinker Smooth Constant_HNe elem 1e6
  print (length $ head abundances)

  plotMassFractions times abundances "igmism.png"

  plotIsotopeAbundances times abundances elem "abundances.png"

  plotMetallicity times abundances "metallicities.png"
