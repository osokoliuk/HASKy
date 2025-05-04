module Main where

import Cosmology
import HMF
import Helper
import IGM
import Lookup
import SMF

args :: [Double]
args = []

main :: IO ()
main =
  do
    (k_arr, pk_arr) <- powerSpectrum "data/CAMB_Pk_z=0.txt"

    (mass_arr, yield_arr) <- yieldsHighMass 1 $ Element "C" 12

    let interp_pk :: PowerSpectrum
        interp_pk = makeInterp k_arr pk_arr

        interp_yield :: Yield
        interp_yield = makeInterp mass_arr yield_arr

        x = igmIsmEvolution planck18 interp_pk Kroupa DoublePower ST Smooth interp_yield 1e6
    print x
