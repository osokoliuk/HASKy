module Main where

import Cosmology
import HMF
import Helper
import IGM
import Lookup
import Pk
import SMF
import StarFormation

args :: [Double]
args = []

main :: IO ()
main =
  do
    -- Fix this interpolation later, this code is not prod ready...
    let pk = powerSpectrumEisensteinHu planck18
        sf_model =
          MkStarFormation
            { yield_ia = 1,
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
        x = interGalacticMediumTerms planck18 sf_model pk Pereira Kroupa Behroozi ST Smooth Constant_HNe metal_frac 1e6 sfrd 0
     in -- x = (\mh -> sqrt $ escapeVelocitySq planck18 pk ST Smooth mh 0) <$> ((10 **) <$> [6, 6 + 0.1 .. 16])
        -- x = (\z -> baryonFormationRateDensity planck18 pk ST Smooth z) <$> z_arr
        -- x = (\m -> m - massRemnant m 0.1) <$> [0.1, 1.1 .. 100]

        -- x = igmIsmMetallicity planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield_ii yield_ia yield_nsm 1e11
        print x
