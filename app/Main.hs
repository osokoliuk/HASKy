module Main where

import Cosmology
import Data.Char (toLower)
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
main =
  do
    -- Fix this interpolation later, this code is not prod ready...
    let 
      pk = powerSpectrumEisensteinHu planck18
      elem = Element {element = "C", isotope = 12}
 
    mass_time <- igmIsmEvolution planck18 pk Pereira Kroupa Behroozi ST Smooth Constant_HNe elem 1e6
    print mass_time
