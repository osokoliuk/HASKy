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
    (mass_arr, yield_arr) <- yieldsHighMass 1 $ Element "C" 12

    let interp_yield :: Yield
        interp_yield = makeInterp mass_arr yield_arr

        pk = powerSpectrumEisensteinHu planck18

        -- x = interGalacticMediumTerms planck18 (powerSpectrumEisensteinHu planck18) Kroupa DoublePower ST Smooth interp_yield 1e6 [20, 20 - 1 .. 0]
        x = (\mh -> 1.988 * 1e43 * escapeVelocitySq planck18 pk Tinker TopHat 1e11 mh) <$> [0, 0 + 0.5 .. 20]

    print $ x
