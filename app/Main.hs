{-# LANGUAGE RecordWildCards #-}

module Main where

import Control.Concurrent.Async (mapConcurrently)
import Control.Lens
import Control.Parallel.Strategies
import Cosmology
import DSP.Basic (logspace)
import Data.Char (toLower)
import Data.Colour
import Data.Colour.Names
import Data.Default.Class
import Data.List
import qualified Data.Vector as V
import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Cairo
import Graphics.Rendering.Chart.Easy
import HMF
import Helper
import IGM
import Lookup
import Pk
import SMF
import StarFormation
import System.Directory (doesFileExist)

args :: [Double]
args = []

main :: IO ()
main = do
  -- Fix this interpolation later, this code is not prod ready...
  let pk = powerSpectrumEisensteinHu planck18
      elem = [Element {element = "Fe", isotope = 56}, Element {element = "Si", isotope = 28}, Element {element = "C", isotope = 12}, Element {element = "O", isotope = 16}]
      sfCfg = MkStarFormationCfg {model_ia = "iwamoto99/WDD1", model_ccsn = "WW95", model_agb = "Cristallo11", model_ecsn = "Wanajo13", model_hne = "Kobayashi06"}
      ts =
        parMap rpar (\z -> cosmicTime planck18 z) zs
      sfrd = makeInterp ts (parMap rpar (\z -> starFormationRateDensity planck18 pk DoublePower ST TopHat z 1e6) zs)
      hmf = (\mh -> haloMassFunction planck18 pk ST TopHat mh 10) <$> ((\x -> 10 ** x) <$> [6.0, 6.5 .. 16])
      pk_approx = (\k -> pk 10 k) <$> ((\x -> 10 ** x) <$> [-3, -2.75 .. 3])
      ccsn_integrand t m = (normImf planck18 Kroupa m) * sfrd (t - tauMS m)
      IGMParams {..} = defaultIGMParams
      PhysicalConstants {..} = phys
      MassLimits {..} = masses
      Efficiencies {..} = effs
      DelayTimes {..} = delays
      NovaeParams {..} = novae
      ECSNParams {..} = ecsn

      first_term t =
        bRG
          * makeIntegrator P128 (\m -> normImfSN planck18 mDLRG mDURG m * sfrd (t - tauMS m)) mDLRG mDURG
      second_term t =
        bMS
          * makeIntegrator P128 (\m -> normImfSN planck18 mDLMS mDUMS m * sfrd (t - tauMS m)) mDLMS mDUMS
      snia z =
        makeIntegrator P128 (\m -> normImf planck18 Kroupa m) (maximum [mPL, mDynamicalRedshift planck18 z]) mPU
          * (first_term (interpT planck18 z) + second_term (interpT planck18 z))

  (times, redshifts, abundances) <- igmIsmEvolution sfCfg planck18 pk Pereira Kroupa DoublePower Tinker Smooth Constant_HNe elem 1e6
  -- imf <- pure $ parMap rpar (\m -> normImf planck18 Kroupa m) $ logspace (-2) 2 50
  -- ccsn <- pure $ makeIntegrator P128 (ccsn_integrand 0) (mDown planck18 0)
  -- print $ parMap rpar (\z -> (makeIntegrator P512 (\m -> ccsn_integrand (interpT planck18 z) m) (mDown planck18 z) 100)) zs
  -- print $ parMap rpar (\t -> sfrd t) ts

  toFile def "igmism.png" $
    do
      layout_title .= "ISM/IGM mass fractions"
      setColors [opaque blue, opaque red, opaque green, opaque orange]
      plot (line "IGM" [zip times $ ((\x -> x V.! 0) <$> abundances)])
      plot (line "ISM" [zip times $ ((\x -> x V.! 1) <$> abundances)])
      plot (line "Stars" [zip times $ ((\x -> x V.! length abundances - 1) <$> abundances)])

  toFile def "abundances.png" $
    do
      layout_title .= "ISM abundances"
      setColors (jetColors . length $ elem)
      let nCols = length elem
      mapM_ (\i -> plot $ line (show (elem !! i)) [zip times ((\x -> log10 x) <$> (\v -> v V.! (5 + 2 * i)) <$> abundances)]) [0 .. length elem - 1]
