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
    (mass_arr, yield_arr) <- yieldsHighMass 1 $ Element "O" 16

    let interp_yield :: Yield
        interp_yield = makeInterp mass_arr yield_arr

        pk = powerSpectrumEisensteinHu planck18
        coeff = 1.988 * 1e43
        z_arr = [20, 20 - 1 .. 0]
        -- x = map (\z -> 0.05 / 0.315 * massAccretionRate planck18 1e6 z) z_arr
        -- x = interGalacticMediumTerms planck18 (powerSpectrumEisensteinHu planck18) Kroupa DoublePower ST Smooth interp_yield 1e6 [20, 20 - 1 .. 0]
        -- x = (\mh -> sqrt $ escapeVelocitySq planck18 pk ST Smooth mh 0) <$> ((10 **) <$> [6, 6 + 0.1 .. 16])
        -- x = (\z -> baryonFormationRateDensity planck18 pk ST Smooth z) <$> z_arr
        -- x = (\m -> m - massRemnant m 0.1) <$> [0.1, 1.1 .. 100]
        x = igmIsmEvolution planck18 pk Kroupa Behroozi ST Smooth interp_yield 1e7
    -- x = (\z -> baryonAccretionRate planck18 pk ST Smooth 1e6 z) <$> z_arr

    print $ x
