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
        coeff = 1.988 * 1e43
        -- x = interGalacticMediumTerms planck18 (powerSpectrumEisensteinHu planck18) Kroupa DoublePower ST Smooth interp_yield 1e6 [20, 20 - 1 .. 0]
        -- x = (\mh -> sqrt $ escapeVelocitySq planck18 pk ST Smooth mh 0) <$> ((10 **) <$> [6, 6 + 0.1 .. 16])
        x = igmIsmEvolution planck18 pk Kroupa Behroozi ST Smooth interp_yield 1e6

    print $ x
