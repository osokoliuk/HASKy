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
    let elem = Element "He" 4
    (mass_arr, yield_arr) <- yieldsHighMass 1 elem
    yield_ia <- yields_Ia elem

    let interp_yield_ii :: Yield_II
        interp_yield_ii m aa = (makeInterp mass_arr yield_arr) m

        pk = powerSpectrumEisensteinHu planck18
        coeff = 1.988 * 1e43
        z_arr = [20, 20 - 1 .. 0]
        -- x = map (\z -> 1e9 * baryonFormationRateDensity planck18 pk ST Smooth z) z_arr
        metal_frac :: Metallicity
        metal_frac x = 1
        sfrd :: SFRD
        sfrd x = 1e-3
        -- x = interGalacticMediumTerms planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield_ii yield_ia metal_frac 1e6 0
        -- x = (\mh -> sqrt $ escapeVelocitySq planck18 pk ST Smooth mh 0) <$> ((10 **) <$> [6, 6 + 0.1 .. 16])
        -- x = (\z -> baryonFormationRateDensity planck18 pk ST Smooth z) <$> z_arr
        -- x = (\m -> m - massRemnant m 0.1) <$> [0.1, 1.1 .. 100]
        x = igmIsmEvolution planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield_ii yield_ia elem 1e7

    -- x = igmMetallicity planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield 0 1e7
    -- x = (\z -> baryonAccretionRate planck18 pk ST Smooth 1e6 z) <$> z_arr

    print $ x
