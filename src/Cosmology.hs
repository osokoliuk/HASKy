module Cosmology where

{-
Module      : HASKy.Helper
Description : Helper module
Copyright   : (c) Oleksii Sokoliuk, 20256
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

A module that defines reference cosmology. Pretty much used by
every other module in this library.
-}

-- Some of the imports

import Math.GaussianQuadratureIntegration
import System.Environment

-- Define a cosmology datatypes, which includes the following:
--    * h0'   -> Hubble parameter in [km s^-1 Mpc^-1]
--    * om0'  -> mass fraction of baryons + CDM
--    * ob0'  -> analoguously, mass fraction of baryons only
--    * tcmb0' -> temperature of the CMB at the current time in [K]
--    * gn'   -> Newton's gravitational constant, in [km^2 Mpc Msun^-1 s^-2]
--    * as' -> 1e9 * amplitude of the primordial scalar power spectrum
--    * ns' -> spectral index of the primordial scalar power spectrum
data ReferenceCosmology
  = MkCosmology
  { h0' :: Double,
    om0' :: Double,
    ob0' :: Double,
    tcmb0' :: Double,
    gn' :: Double,
    as' :: Double,
    ns' :: Double
  }
  deriving (Eq, Show)

type Redshift = Double

type CosmicTime = Double

-- | Unpack the values in the cosmology_record into variables
unpackCosmology :: ReferenceCosmology -> (Double, Double, Double, Double, Double, Double, Double)
unpackCosmology cosmology =
  (,,,,,,) <$> h0' <*> om0' <*> ob0' <*> tcmb0' <*> gn' <*> as' <*> ns' $ cosmology

-- | Example Planck 2018 cosmology
planck18 :: ReferenceCosmology
planck18 =
  MkCosmology
    { h0' = 67.66,
      om0' = (0.02242 / (67.66 / 100) ** 2 + 0.11933 / (67.66 / 100) ** 2),
      ob0' = 0.02242 / (67.66 / 100) ** 2,
      tcmb0' = 2.7255,
      gn' = 4.301 * 10 ** (-9),
      as' = 2.105,
      ns' = 0.9665
    }

-- | Initialise cosmological model from the args, to be provided from the CLI
initialiseCosmology :: [String] -> ReferenceCosmology
initialiseCosmology args =
  case args of
    [h0, om0, ob0, tcmb0, gn, as, ns] ->
      let [h0, om0, ob0, tcmb0, gn, as, ns] = read <$> args
       in MkCosmology h0 om0 ob0 tcmb0 gn as ns
    _ -> MkCosmology 0 0 0 0 0 0 0

-- | Fiducial Lambda CDM Hubble parameter, in the units of [km s^-1 Mpc^-1]
hubbleParameter :: ReferenceCosmology -> Redshift -> Double
hubbleParameter cosmology z =
  let (h0, om0, ob0, _, _, _, _) = unpackCosmology cosmology
   in h0 * sqrt (om0 * (1 + z) ** 3 + 1 - om0)

-- | dt/dz, will be used as an integrand to derive cosmic time
dtdz :: ReferenceCosmology -> Redshift -> Double
dtdz cosmology z =
  let (h0, om0, ob0, _, _, _, _) = unpackCosmology cosmology
      km_Mpc = 3.24 * 1e-20
      s_Gyr = 3.15 * 1e16
   in (km_Mpc * s_Gyr * (1 + z) * hubbleParameter cosmology z) ** (-1)

-- | Define cosmic time t(z) in the units of [Gyr]
cosmicTime :: ReferenceCosmology -> Redshift -> CosmicTime
cosmicTime cosmology z =
  nIntegrate128 (dtdz cosmology) z 20
