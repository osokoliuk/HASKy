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
    let elem = Element "Fe" 56
        metal_fracs = [1e-3, 1e-2, 1e-1, 1e0]
        masses = [1, 5 .. 350]
        argument = [(x, elem) | x <- metal_fracs]
    yield_ii <- mapM (uncurry yieldsHighMass) argument

    yield_ia <- yields_Ia elem

    -- Fix this interpolation later, this code is not prod ready...
    let yield_nsm = 0
        interp_yield_ii m metal_frac =
          let (masses, yields)
                | metal_frac <= 0.001 = yield_ii !! 0
                | metal_frac <= 0.01 = yield_ii !! 1
                | metal_frac <= 0.1 = yield_ii !! 2
                | otherwise = yield_ii !! 3
              interp_mass = makeInterp masses yields
           in interp_mass m

        pk = powerSpectrumEisensteinHu planck18
        -- x = map (\z -> 1e9 * baryonFormationRateDensity planck18 pk ST Smooth z) z_arr
        -- x = interGalacticMediumTerms planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield_ii yield_ia metal_frac 1e6 sfrd 0
        -- x = (\mh -> sqrt $ escapeVelocitySq planck18 pk ST Smooth mh 0) <$> ((10 **) <$> [6, 6 + 0.1 .. 16])
        -- x = (\z -> baryonFormationRateDensity planck18 pk ST Smooth z) <$> z_arr
        -- x = (\m -> m - massRemnant m 0.1) <$> [0.1, 1.1 .. 100]

        x = igmIsmMetallicity planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield_ii yield_ia yield_nsm 1e11
     in print x
