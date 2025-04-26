module SMF where

import HMF

data SMF_kind
  = DoublePower
  | Behroozi
  | EMERGE
  deriving (Eq, Show)

log10 :: (Floating a) => a -> a
log10 x = log x / log 10

epsStar :: SMF_kind -> Mhalo -> Redshift -> Double
epsStar s_kind mh z =
  let eps_0, mh_0, gamma_lo, gamma_hi :: Double
      (eps_0, mh_0, gamma_lo, gamma_hi) = (0.05, 2.8 * 1e11, 0.49, -0.61)

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

      a, nu, epsilon, m_1, alpha, delta, gamma :: Double
      a = 1 / (1 + z)
      nu = exp (-4 * a ** 2)
      epsilon = e0 + (e1 * (a - 1) + e2 * z) * nu + e3 * (a - 1)
      m_1 = m0 + (m1 * (a - 1) + m2 * z) * nu
      alpha = a0 + a1 * (a - 1) * nu
      delta = d0 + (d1 * (a - 1) + d2 * z) * nu
      gamma = g0 + (g1 * (a - 1) + g2 * z) * nu
      fBehroozi x = -log10 (10 ** (alpha * x) + 1) + delta * (log10 (1 + exp (x))) ** gamma / (1 + exp (10 ** (-x)))
   in case s_kind of
        DoublePower -> eps_0 / ((mh_0) ** gamma_lo + (mh / mh_0) ** gamma_hi)
        Behroozi -> 10 ** (log10 (epsilon * m_1) + fBehroozi (log10 (mh / m_1)) - fBehroozi (0))
