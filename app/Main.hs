module Main where

import Cosmology
import HMF
import Helper
import IGM
import Lookup
import Pk
import SMF

args :: [Double]
args = []

main :: IO ()
main =
  do
    (mass_arr, yield_arr) <- yieldsHighMass 1 $ Element "H" 1

    let interp_yield :: Yield
        interp_yield = makeInterp mass_arr yield_arr

        pk = powerSpectrumEisensteinHu planck18

        x = igmIsmEvolution planck18 (powerSpectrumEisensteinHu planck18) Kroupa DoublePower ST Smooth interp_yield 1e6

    print x
