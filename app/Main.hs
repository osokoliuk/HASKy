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
    x <- retrieveYieldIa sf_cfg elem
    print x
    y <- doesFileExist filepath
    print filepath
    print y
  where
    pk = powerSpectrumEisensteinHu planck18
    sf_model =
      MkStarFormation
        { yield_ia = retrieveYieldIa sf_cfg elem,
          yield_ccsn = ([1], [1]),
          yield_hne = ([1], [1]),
          yield_ecsn = ([1], [1]),
          yield_agb = ([1], [1]),
          yield_sagb = ([1], [1]),
          yield_nsm = 1,
          yield_psn = ([1], [1]),
          yield_co = 1,
          yield_one = 1
        }

    metal_frac x = 1.0
    sfrd x = 1.0
    elem = Element {element = "C", isotope = 12}

    -- x = interGalacticMediumTerms planck18 sf_model pk Pereira Kroupa Behroozi ST Smooth Constant_HNe metal_frac 1e6 sfrd 0
    -- x = igmIsmEvolution planck18 sf_model pk Pereira Kroupa Behroozi ST Smooth Constant_HNe (Element "Fe" 54) 1e6
    sf_cfg = MkStarFormationCfg {model_ia = "iwamoto99/WDD1"}
    filepath = "data/" ++ model_ia sf_cfg ++ "/" ++ (toLower <$> element elem) ++ ".dat"
