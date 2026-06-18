module Main where

import Control.Concurrent.Async (mapConcurrently)
import Control.Parallel.Strategies
import Cosmology
import DSP.Basic (logspace)
import Data.Char (toLower)
import Data.List
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
main =
  do
    -- Fix this interpolation later, this code is not prod ready...
    let pk = powerSpectrumEisensteinHu planck18
        elem = Element {element = "Fe", isotope = 56}
        sfCfg = MkStarFormationCfg {model_ia = "iwamoto99/WDD1", model_ccsn = "WW95", model_agb = "Cristallo11", model_ecsn = "Wanajo13", model_hne = "Kobayashi06"}

        defaultIGMParams :: IGMParams
        defaultIGMParams =
          IGMParams
            { phys =
                PhysicalConstants
                  { mCO = 1.38,
                    yrGyr = 1e9,
                    kmsErgMsol = 1.989e43,
                    snEnergy = 2e51
                  },
              masses =
                MassLimits
                  { mUp = 100,
                    mPU = 8.0,
                    mPL = 3.0,
                    mDURG = 1.5,
                    mDLRG = 0.9,
                    mDUMS = 2.6,
                    mDLMS = 1.8,
                    mNSMd = 9.0,
                    mNSMu = 30.0,
                    mNovaeD = 0.8,
                    mNovaeU = 8.0,
                    mAGBd = 1.3,
                    mAGBu = 8.0,
                    mECSNd = 8.8,
                    mECSNu = 9.0
                  },
              effs =
                Efficiencies
                  { epsW = 0.02,
                    epsSN = 0.005,
                    bRG = 0.02,
                    bMS = 0.04,
                    epsHNe0 = 0.5,
                    epsMRSNe = 0.03,
                    alphaNSM = 0.018,
                    alphaNovae = 0.01,
                    epsCO = 0.7
                  },
              delays =
                DelayTimes
                  { delayNSM = 1e7,
                    delayNovae = 2e9
                  },
              novae =
                NovaeParams
                  { mEjNovae = 2e-5,
                    nNovae = 1e4
                  },
              ecsn =
                ECSNParams
                  { mEjECSN = 1.14e-2
                  }
            }
        zs = [20.0, 20.0 - 0.5 .. 0]
        ts =
          parMap rpar (\z -> cosmicTime planck18 z) zs
        sfrd = (parMap rpar (\z -> starFormationRateDensity planck18 pk DoublePower ST Smooth z 1e6) zs)
        hmf = (\mh -> haloMassFunction planck18 pk ST Smooth mh 10) <$> ((\x -> 10 ** x) <$> [6.0, 6.5 .. 16])
        pk_approx = (\k -> pk 10 k) <$> ((\x -> 10 ** x) <$> [-3, -2.75 .. 3])
    -- ccsn_integrand z m = mkIntegrand planck18 (normImf planck18 Kroupa) sfrd (tauMS m) (const 1) z m

    -- mass_time <- igmIsmEvolution sfCfg planck18 pk Pereira Kroupa Behroozi ST Smooth Constant_HNe elem 1e6
    -- igm <- mapConcurrently (\z -> snTermsIO sfCfg planck18 pk Pereira Kroupa Behroozi ST Smooth Constant_HNe (\x -> 1e-3) 1e6 sfrd z elem) zs
    -- imf <- pure $ parMap rpar (\m -> normImf planck18 Kroupa m) $ logspace (-2) 2 50
    -- ccsn <- pure $ parMap rpar (\z -> makeIntegrator P128 (ccsn_integrand z) (mDown planck18 z) 100) [20.0, 20.0 - 0.5 .. 0]
    -- print $ parMap rpar (\z -> ccsn_integrand z 8) zs
    print $ pk_approx
