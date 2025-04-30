module IGM where

import Cosmology
import HMF
import Helper
import SMF

data IMF_kind
  = Salpeter
  | MillerScalo
  | Kroupa
  | Chabrier
  deriving (Eq, Show)

initialMassFunction :: IMF_kind -> Mhalo -> Double
initialMassFunction i_kind m =
  let alpha0, alpha1, alpha2, k0, k1, k2 :: Double
      (alpha0, alpha1, alpha2, m1, m2, k0, k1, k2) =
        (-0.3, -1.3, -2.3, 1, 0.08, 0.5, k0 * m1 ** (alpha0 - alpha1), k1 * m2 ** (alpha1 - alpha2))

      a_Ch, center_Ch, sigma_Ch :: Double
      (a_Ch, center_Ch, sigma_Ch) = (0.086, 0.57, 0.22)
   in case i_kind of
        -- Salpeter et al. 1955 IMF (single power-law)
        Salpeter -> m ** (-2.35)
        -- Kroupa et al. 2001 IMF (broken power-law)
        Kroupa
          | m < m1 -> k0 * m ** alpha0
          | m >= m2 && m < m2 -> k1 * m ** alpha1
          | otherwise -> k2 * m ** alpha2
        -- Chabrier et al. 2003 (log-normal) converted to [Mpc^-3]
        Chabrier
          | m < 1 -> a_Ch * exp (-(log m - log center_Ch) ** 2 / (2 * sigma_Ch ** 2))
          | otherwise -> m ** (-1.3)
