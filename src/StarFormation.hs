module StarFormation where

{-
Module      : HASKy.StarFormation
Description : StarFormation module
Copyright   : (c) Oleksii Sokoliuk, 2025
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

A module that defines reference star formation model. Pretty much used by
every other module in this library.
-}

import Helper
import Math.GaussianQuadratureIntegration
import System.Environment

-- Define a cosmology datatypes, which includes the following:
--    * yields_II -> Element/Isotope yields model for SN type II, in [Msol]
--    * yields_Ia -> Element/Isotope yields model for SN type Ia, in [Msol]
--    * yields_NSM  -> Element/Isotope yields model for NS-NS and NS-BH mergers, in [Msol]
data ReferenceStarFormationModel
  = MkStarFormation
  { yields_II :: [Char],
    yields_Ia :: [Char],
    yields_NSM :: [Char]
  }
  deriving (Eq, Show)
