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
    (mass_arr, yield_arr) <- yieldsHighMass 1 elem
    yield_ia <- yields_Ia elem

    let interp_yield_ii :: Yield_II
        interp_yield_ii m aa = (makeInterp mass_arr yield_arr) m

        pk = powerSpectrumEisensteinHu planck18
        coeff = 1.988 * 1e43
        z_arr = [20, 20 - 1 .. 0]
        sfrd_arr = [24530.09869163431, 35822.861675889355, 51850.3598353403, 74892.7980540064, 108355.51640473121, 157287.32839658318, 229153.6156463161, 334940.1723087513, 490698.47535302996, 719691.7858093071, 1055347.6083541322, 1545141.927004186, 2255765.7667246205, 3279455.2851683698, 4741912.106333696, 6811144.932782047, 9707011.147283336, 1.3710533499749003e7, 1.9168928298233483e7, 2.6495529196722653e7, 3.615617087670479e7, 4.8631717802108966e7, 6.435175425091211e7, 8.356729762959377e7, 1.061620409608078e8, 1.3138271944922253e8, 1.5752692355362806e8, 1.8173101178121063e8, 2.0012639867724532e8, 2.0870928288348246e8, 2.0495940936830974e8, 1.8930995103883353e8, 1.652400084569673e8, 1.375290218982621e8, 1.101254444916516e8, 8.516908548160467e7, 6.330773795331535e7, 4.4442662040499136e7, 2.80926955618159e7, 1.455996367536503e7, 5445438.810876669]
        -- x = map (\z -> 1e9 * baryonFormationRateDensity planck18 pk ST Smooth z) z_arr
        metal_frac :: Metallicity
        metal_frac x = 1
        sfrd :: SFRD
        t_arr = [0.0, 1.3658033847049145, 2.731606769409829, 4.097410154114743, 5.463213538819658, 6.829016923524573, 8.194820308229486, 9.560623692934401, 10.926427077639316, 12.29223046234423, 13.658033847049145]
        sfrd = makeInterp t_arr sfrd_arr
        -- x = interGalacticMediumTerms planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield_ii yield_ia metal_frac 1e6 sfrd 0
        -- x = (\mh -> sqrt $ escapeVelocitySq planck18 pk ST Smooth mh 0) <$> ((10 **) <$> [6, 6 + 0.1 .. 16])
        -- x = (\z -> baryonFormationRateDensity planck18 pk ST Smooth z) <$> z_arr
        -- x = (\m -> m - massRemnant m 0.1) <$> [0.1, 1.1 .. 100]
        x = igmIsmEvolution planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield_ii yield_ia elem 1e7

    -- x = igmMetallicity planck18 pk Pereira Kroupa Behroozi ST Smooth interp_yield 0 1e7
    -- x = (\z -> baryonAccretionRate planck18 pk ST Smooth 1e6 z) <$> z_arr

    print $ x
