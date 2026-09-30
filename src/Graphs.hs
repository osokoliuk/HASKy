{-# LANGUAGE OverloadedStrings #-}

module Graphs where

import qualified Data.Map as M
import qualified Data.Vector as V

{-
Module      : HASKy.Graphs
Description : Graphs module
Copyright   : (c) Oleksii Sokoliuk, 20256
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

A module that defines a bunch of very useful routines for plotting specific
outputs of HASKy solver
-}

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Cairo
import qualified Graphics.Rendering.Chart.Easy as C
import Helper
import IGM

-- Plot general information

plotMassFractions ::
    [Double] ->
    [V.Vector Double] ->
    [Char] ->
    IO ()
plotMassFractions times abundances@(a : as) filename = toFile C.def filename $
    do
        layout_title C..= "ISM/IGM mass fractions"
        C.setColors [C.opaque C.blue, C.opaque C.red, C.opaque C.green, C.opaque C.orange]
        C.plot (C.line "IGM" [zip times $ ((\x -> x V.! 0) <$> abundances)])
        C.plot (C.line "ISM" [zip times $ ((\x -> x V.! 1) <$> abundances)])
        C.plot (C.line "Stars" [zip times $ ((\x -> x V.! (V.length a - 3)) <$> abundances)])

plotIsotopeAbundances ::
    [Double] ->
    M.Map String (M.Map Element (IsotopeInformation [Double])) ->
    [String] ->
    [[Element]] ->
    [Char] ->
    IO ()
plotIsotopeAbundances times abundances elems groupedIsotopes filename =
    toFile C.def filename $ do
        layout_title C..= "ISM abundances"
        C.setColors (jetColors (length (concat groupedIsotopes)))

        mapM_
            ( \i ->
                mapM_
                    ( \j ->
                        case M.lookup (elems !! i) abundances >>= M.lookup (groupedIsotopes !! i !! j) of
                            Nothing -> pure ()
                            Just isotopeInfo ->
                                C.plot $
                                    C.line
                                        (show (groupedIsotopes !! i !! j))
                                        [zip times (map (log10) (ismFraction isotopeInfo))]
                    )
                    [0 .. length (groupedIsotopes !! i) - 1]
            )
            [0 .. length elems - 1]

plotMetallicity :: [Double] -> [V.Vector Double] -> [Char] -> IO ()
plotMetallicity times abundances@(a : as) filename = toFile C.def filename $
    do
        layout_title C..= "ISM/IGM metallicity"
        C.setColors [C.opaque C.blue, C.opaque C.red]
        C.plot (C.line "Z_IGM" [zip times $ ((\x -> x V.! (V.length a - 1)) <$> abundances)])
        C.plot (C.line "Z_ISM" [zip times $ ((\x -> x V.! (V.length a - 2)) <$> abundances)])

-- Plot ejection / outflow rates from different sources

-- Plot log-scale isotope ratio phase space
