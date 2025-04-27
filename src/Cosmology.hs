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
