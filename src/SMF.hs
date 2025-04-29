module SMF where

{-
Module      : HASKy.SMF
Description : Stellar Mass Function
Copyright   : (c) Oleksii Sokoliuk, 20256
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

This module defines a bunch of routines that in the end
yield a Stellar Mass Function for a given cosmology (i.e., values of
Hubble parameter H0, Omega_m0, Omega_b0)
-}

-- Import HMF, Cosmology modules (to be used to derive SMF)

import Cosmology
import qualified Data.Map as M
import Data.Maybe (fromMaybe)
import HMF
import Helper
import Math.GaussianQuadratureIntegration
import Numeric.Tools.Differentiation

-- As usual, specify all kinds of star formation efficiencies that we consider
data SMF_kind
  = DoublePower
  | Behroozi
  | EMERGE
  deriving (Eq, Show)

-- | Star formation efficiency, i.e. Star mass / Halo mass
-- There are three possible choices for the SFE model:
--    * Simple double power-law
--    * Behroozi et al. 2013 model
--    * EMERGE semi-analytical model
epsStar :: ReferenceCosmology -> SMF_kind -> Mhalo -> Redshift -> Double
epsStar cosmology s_kind mh z =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology

      -- A set of best fit parameter for the double power-law
      eps_0, mh_0, gamma_lo, gamma_hi :: Double
      (eps_0, mh_0, gamma_lo, gamma_hi) = (0.21, 2.8 * 1e11, 0.49, -0.61)

      -- A set of best fit parameters for the Behrozoi et al. 2013 SFE model
      e0, e1, e2, e3, m0, m1, m2, a0, a1, d0, d1, d2, g0, g1, g2 :: Double
      (e0, e1, e2, e3, m0, m1, m2, a0, a1, d0, d1, d2, g0, g1, g2) =
        ( -1.777,
          -0.006,
          -0.000,
          -0.119,
          11.514,
          -1.793,
          -0.251,
          -1.412,
          0.731,
          3.508,
          2.608,
          -0.043,
          0.316,
          1.319,
          0.279
        )

      -- Some helper parameters for Behroozi et al. 2013 SFE model
      a, nu, epsilon, m_1, alpha, delta, gamma :: Double
      a = 1 / (1 + z)
      nu = exp (-4 * a ** 2)
      epsilon = (10 **) $ e0 + (e1 * (a - 1) + e2 * z) * nu + e3 * (a - 1)
      m_1 = (10 **) $ m0 + (m1 * (a - 1) + m2 * z) * nu
      alpha = a0 + (a1 * (a - 1)) * nu
      delta = d0 + (d1 * (a - 1) + d2 * z) * nu
      gamma = g0 + (g1 * (a - 1) + g2 * z) * nu
      fBehroozi x = -log10 (10 ** (alpha * x) + 1) + delta * (log10 (1 + exp (x))) ** gamma / (1 + exp (10 ** (-x)))
   in case s_kind of
        DoublePower -> eps_0 / ((mh / mh_0) ** gamma_lo + (mh / mh_0) ** gamma_hi)
        Behroozi ->
          let mstar = 10 ** (log10 (epsilon * m_1) + fBehroozi (log10 (mh / m_1)) - fBehroozi (0))
           in (mstar / mh) / (ob0 / om0)

-- | Mass accretion history for the CDM halo with
-- parameters alpha_MAR and beta_MAR taken to be an average for the LCDM model
-- within the 1e8 <= Mh <= 1e14 bounds
massAccretionHistory :: Mhalo -> Redshift -> Double
massAccretionHistory mh z =
  let alpha_MAR = 0.24
      neta_MAR = -0.75
   in mh * (1 + z) ** (alpha_MAR) * exp (neta_MAR * z)

-- | Mass accretion rate, i.e., rate at which halo of mass Mh gains mass,
-- adopted in the units of a solar mass from the [Fakhouri et al. 2013] work
massAccretionRate :: ReferenceCosmology -> Mhalo -> Redshift -> Double
massAccretionRate cosmology mh z =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology
   in 25.3 * (mh / 1e12) ** 1.1
        + (1 + 1.65 * z) * sqrt (om0 * (1 + z) ** 3 + 1 - om0)

-- | Star formation rate, derived simply as a normalised halo mass accretion rate
-- eps_star converts baryonic mass to a stellar mass, while factor Om0/Ob0 converts
-- halo mass to baryonic mass
starFormationRate :: SMF_kind -> ReferenceCosmology -> Mhalo -> Redshift -> Double
starFormationRate s_kind cosmology mh z =
  let (h0, om0, ob0, c, gn) = unpackCosmology cosmology
      ep = epsStar cosmology s_kind mh z
   in ep * ob0 / om0 * massAccretionRate cosmology mh z

-- | Stellar mass function, derived from the HMF and SFE via a simple chain rule
stellarMassFunction :: FilePath -> ReferenceCosmology -> SMF_kind -> HMF_kind -> W_kind -> [Mhalo] -> Redshift -> IO ([Double], [Double])
stellarMassFunction filepath cosmology s_kind h_kind w_kind mh_arr z =
  do
    hmf_arr <- haloMassFunction filepath cosmology h_kind w_kind mh_arr z

    let (h0, om0, ob0, c, gn) = unpackCosmology cosmology

        ms :: Mhalo -> Double -- Function that gives stellar mass
        ms mh = mh * (ob0 / om0) * epsStar cosmology s_kind mh z

        ms_arr :: [Double]
        ms_arr = ms <$> mh_arr

        ln_factors :: [Double] -- Turns dMh/dMstar into dlnMh/dlog10Mstar
        ln_factors = zipWith (\x y -> x * log 10 / y) ms_arr mh_arr

        dmhdms :: [Double] -- Part of the chain rule to turn HMF into SMF
        dmhdms =
          zipWith (/) ln_factors $
            (\mh -> diffRes $ diffRichardson ms 10 mh) <$> mh_arr

    return $ (ms_arr, zipWith (\x y -> x * y) hmf_arr dmhdms)

-- | Star formation rate density,
-- to be used in the IGM/ISM mass fraction differential equations
starFormationRateDensity :: FilePath -> ReferenceCosmology -> SMF_kind -> HMF_kind -> W_kind -> Redshift -> IO Double
starFormationRateDensity filepath cosmology s_kind h_kind w_kind z =
  let mh_arr = (10 **) <$> [6, 6 + 0.1 .. 18]
   in do
        hmf_arr <- haloMassFunction filepath cosmology h_kind w_kind mh_arr z

        let dndmh :: [Double] -- Convert dn/dlnMh to dn/dMh
            dndmh = zipWith (/) hmf_arr mh_arr

            interp_hmf :: Double -> Double
            interp_hmf = mapLookup (M.fromList (zip mh_arr dndmh))

            integrand :: Double -> Double -- Function giving epsStar * Ob0/Om0 * dn/dt
            integrand mh =
              (interp_hmf mh)
                * (starFormationRate s_kind cosmology mh z)

            result :: Double
            result = nIntegrate512 integrand (minimum mh_arr) (maximum mh_arr)

        return $ result

main_SMF = do
  x <- starFormationRateDensity "data/CAMB_Pk_z=0.txt" planck18 DoublePower ST Smooth 0
  print $ x
