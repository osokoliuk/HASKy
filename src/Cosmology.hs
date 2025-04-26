module Cosmology where

import System.Environment

data ReferenceCosmology
  = MkCosmology
  { h0' :: Double,
    om0' :: Double,
    ob0' :: Double,
    c' :: Double,
    gn' :: Double
  }
  deriving (Eq, Show)

initialiseCosmology :: [String] -> ReferenceCosmology
initialiseCosmology args =
  case args of
    [h0, om0, ob0, c, gn] ->
      let [h0, om0, ob0, c, gn] = read <$> args
       in MkCosmology h0 om0 ob0 c gn
