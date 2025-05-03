module Lookup where

{-
Module      : HASKy.Lookup
Description : Lookup tables
Copyright   : (c) Oleksii Sokoliuk, 20256
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

This module stores all of the lookup tables used in the code.
-}

import Data.Map (fromList)
import qualified Data.Map as M

type Metallicity = Double

data Element
  = H1
  | He4
  | H2
  | He3
  | Li7
  | Be7
  | Be9
  | B10
  | B11
  | C11
  | C12
  | C13
  | C14
  | N14
  | N15
  | O16
  | O17
  | O18
  | F19
  | Ne20
  | Ne21
  | Ne22
  | Na22
  | Mg24
  | Mg25
  | Mg26
  | Al26
  | Al27
  | Si28
  | Si29
  | Si30
  | P31
  | S32
  | S33
  | S34
  | S35
  | S36
  | Cl35
  | Cl36
  | Cl37
  | Ar36
  | Ar37
  | Ar38
  | Ar40
  | K39
  | K40
  | K41
  | Ca40
  | Ca41
  | Ca42
  | Ca43
  | Ca44
  | Ca45
  | Ca46
  | Ca48
  | Sc43
  | Sc45
  | Ti44
  | Ti45
  | Ti46
  | Ti47
  | Ti48
  | Ti49
  | Ti50
  | V47
  | V48
  | V49
  | V50
  | V51
  | Cr48
  | Cr49
  | Cr50
  | Cr51
  | Cr52
  | Cr53
  | Cr54
  | Mn51
  | Mn52
  | Mn53
  | Mn54
  | Mn55
  | Fe52
  | Fe53
  | Fe54
  | Fe55
  | Fe56
  | Fe57
  | Fe58
  | Fe59
  | Fe60
  | Co55
  | Co56
  | Co57
  | Co58
  | Co59
  | Co60
  | Co61
  | Ni56
  | Ni57
  | Ni58
  | Ni59
  | Ni60
  | Ni61
  | Ni62
  | Ni63
  | Ni64
  | Ni65
  | Cu59
  | Cu60
  | Cu61
  | Cu62
  | Cu63
  | Cu64
  | Cu65
  | Cu66
  | Zn60
  | Zn61
  | Zn62
  | Zn63
  | Zn64
  | Zn65
  | Zn66
  | Zn67
  | Zn68
  | Zn69
  | Ga64
  | Ga65
  | Ga66
  | Ga67
  | Ga68
  | Ga69
  | Ga70
  | Ge64
  | Ge65
  | Ge66
  | Ge68
  | Ge69
  | Ge70
  | Ge71
  deriving (Eq, Show, Ord, Enum)

-- | Stellar remnant mass for a white dwarf, taken from the [Hoek & Groenewegen 1996]
remnantMediumMass :: Metallicity -> M.Map Double Double
remnantMediumMass metal_frac
  | metal_frac <= 0.001 =
      fromList $
        [(0.9, 0.56), (1.0, 0.57), (1.3, 0.59), (1.5, 0.62), (1.7, 0.64), (2.0, 0.68), (2.5, 0.74), (3.0, 0.82), (4.0, 0.91), (5.0, 0.98), (7.0, 1.13), (8.0, 1.84)]
  | metal_frac <= 0.004 =
      fromList $
        [(0.9, 0.58), (1.0, 0.58), (1.3, 0.59), (1.5, 0.60), (1.7, 0.60), (2.0, 0.62), (2.5, 0.64), (3.0, 0.66), (4.0, 0.91), (5.0, 0.98), (7.0, 1.12), (8.0, 1.20)]
  | metal_frac <= 0.008 =
      fromList $
        [(0.9, 0.58), (1.0, 0.58), (1.3, 0.59), (1.5, 0.59), (1.7, 0.60), (2.0, 0.61), (2.5, 0.63), (3.0, 0.64), (4.0, 0.90), (5.0, 0.97), (7.0, 1.11), (8.0, 1.20)]
  | metal_frac <= 0.02 =
      fromList $
        [(0.9, 0.57), (1.0, 0.58), (1.3, 0.59), (1.5, 0.59), (1.7, 0.59), (2.0, 0.60), (2.5, 0.60), (3.0, 0.62), (4.0, 0.79), (5.0, 0.92), (7.0, 1.07), (8.0, 1.15)]
  | otherwise =
      fromList $
        [(0.9, 0.55), (1.0, 0.55), (1.3, 0.56), (1.5, 0.57), (1.7, 0.58), (2.0, 0.59), (2.5, 0.60), (3.0, 0.61), (4.0, 0.68), (5.0, 0.77), (7.0, 0.98), (8.0, 1.09)]

-- | Stellar remnant masses for either black hole or a neutron star, taken from the [Woosley & Weaver 1995]
remnantHighMass :: Metallicity -> M.Map Double Double
remnantHighMass metal_frac
  | metal_frac <= 0.001 =
      fromList $
        [(12, 1.28), (13, 1.44), (15, 1.63), (18, 1.61), (20, 1.97), (22, 2.01), (25, 1.87), (30, 2.08), (35, 3.03), (40, 4.09)]
  | metal_frac <= 0.01 =
      fromList $
        [(12, 1.40), (13, 1.44), (15, 1.55), (18, 1.58), (20, 1.98), (22, 2.04), (25, 1.87), (30, 2.21), (35, 2.42), (40, 4.42)]
  | metal_frac <= 0.1 =
      fromList $
        [(12, 1.38), (13, 1.31), (15, 1.49), (18, 1.69), (20, 1.97), (22, 2.12), (25, 1.99), (30, 2.01), (35, 3.39), (40, 4.45)]
  | otherwise =
      fromList $
        [(12, 1.32), (13, 1.46), (15, 1.43), (18, 1.76), (19, 1.98), (20, 2.06), (22, 2.02), (25, 2.07), (30, 1.94), (35, 3.86), (40, 5.45)]

-- | Stellar lifetimes for the higher mass end from the [Schaerer 2002]
tauHighMass :: M.Map Double Double
tauHighMass =
  fromList $
    [(9, 20.22), (15, 10.40), (25, 6.459), (40, 3.864), (60, 3.464), (80, 3.012), (120, 2.521)]
