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
--    * c'    -> speed of light, in [km s^-1]
--    * gn'   -> Newton's gravitational constant, in [km^2 Mpc Msun^-1 s^-2]
data ReferenceCosmology
  = MkCosmology
  { h0' :: Double,
    om0' :: Double,
    ob0' :: Double,
    c' :: Double,
    gn' :: Double
  }
  deriving (Eq, Show)

type Redshift = Double

type CosmicTime = Double

-- | Unpack the values in the cosmology_record into variables
unpackCosmology :: ReferenceCosmology -> (Double, Double, Double, Double, Double)
unpackCosmology cosmology =
  (,,,,) <$> h0' <*> om0' <*> ob0' <*> c' <*> gn' $ cosmology

-- | Example Planck 2018 cosmology
planck18 :: ReferenceCosmology
planck18 =
  MkCosmology
    { h0' = 67.66,
      om0' = (0.02242 / (67.66 / 100) ** 2 + 0.11933 / (67.66 / 100) ** 2),
      ob0' = 0.02242 / (67.66 / 100) ** 2,
      c' = 299792.45800000057,
      gn' = 4.301 * 10 ** (-9)
    }

-- | Initialise cosmological model from the args, to be provided from the CLI
initialiseCosmology :: [String] -> ReferenceCosmology
initialiseCosmology args =
  case args of
    [h0, om0, ob0, c, gn] ->
      let [h0, om0, ob0, c, gn] = read <$> args
       in MkCosmology h0 om0 ob0 c gn
    _ -> MkCosmology 0 0 0 0 0

-- | Fiducial Lambda CDM Hubble parameter, in the units of [km s^-1 Mpc^-1]
hubbleParameter :: ReferenceCosmology -> Redshift -> Double
hubbleParameter cosmology z =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology
   in h0 * sqrt (om0 * (1 + z) ** 3 + 1 - om0)

-- | Define cosmic time t(z) in the units of [Gyr]
cosmicTime :: ReferenceCosmology -> Redshift -> Redshift -> CosmicTime
cosmicTime cosmology zinit z =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology

      km_Mpc :: Double
      km_Mpc = 3.24 * 1e-20

      s_Gyr = 3.15 * 1e16

      integrand :: Double -> Double
      integrand z =
        (km_Mpc * s_Gyr * (1 + z) * hubbleParameter cosmology z) ** (-1)
   in nIntegrate512 integrand z zinit
