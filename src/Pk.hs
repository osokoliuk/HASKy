module Pk where

import Control.Parallel.Strategies
import Cosmology
import HMF
import Helper

-- | Eisentein-Hu transfer function with the BAO wiggles,
-- implementation is essentially adapted from the SkyPy package
transferEinsensteinHu :: Wavenumber -> ReferenceCosmology -> Double
transferEinsensteinHu k cosmology =
  let (h0, om0, ob0, tcmb0, gn, as, ns, _) = unpackCosmology cosmology
      little_h = h0 / 100
      om0h2 = om0 * little_h ** 2
      ob0h2 = ob0 * little_h ** 2
      f_baryon = ob0 / om0

      k_eq = 7.46e-2 * om0h2 * (tcmb0 / 2.7) ** (-2)
      z_eq = 2.5e4 * om0h2 * (tcmb0 / 2.7) ** (-4)

      z_drag_b1 = 0.313 * om0h2 ** (-0.419) * (1 + 0.607 * om0h2 ** 0.674)
      z_drag_b2 = 0.238 * om0h2 ** 0.223
      z_drag =
        1291.0
          * om0h2 ** 0.251
          / (1 + 0.659 * om0h2 ** 0.828)
          * (1 + z_drag_b1 * ob0h2 ** z_drag_b2)

      r_drag = 31.5 * ob0h2 * (tcmb0 / 2.7) ** (-4) * (1000.0 / z_drag)
      r_eq = 31.5 * ob0h2 * (tcmb0 / 2.7) ** (-4) * (1000.0 / z_eq)

      sound_horizon =
        2
          / (3 * k_eq)
          * sqrt (6 / r_eq)
          * log
            ( (sqrt (1 + r_drag) + sqrt (r_drag + r_eq))
                / (1 + sqrt (r_eq))
            )
      k_silk = 1.6 * ob0h2 ** 0.52 * om0h2 ** 0.73 * (1 + (10.4 * om0h2) ** (-0.95))

      alpha_c_a1 = (46.9 * om0h2) ** 0.670 * (1 + (32.1 * om0h2) ** (-0.532))
      alpha_c_a2 = (12.0 * om0h2) ** 0.424 * (1 + (45.0 * om0h2) ** (-0.582))
      alpha_c = alpha_c_a1 ** (-f_baryon) * alpha_c_a2 ** (-f_baryon ** 3)

      beta_c_b1 = 0.944 / (1 + (458 * om0h2) ** (-0.708))
      beta_c_b2 = (0.395 * om0h2) ** (-0.0266)
      beta_c = 1 / (1 + beta_c_b1 * ((1 - f_baryon) ** beta_c_b2 - 1))

      y = (1.0 + z_eq) / (1 + z_drag)
      alpha_b_G =
        y
          * ( -6 * sqrt (1 + y)
                + (2 + 3 * y)
                  * log ((sqrt (1 + y) + 1) / (sqrt (1 + y) - 1))
            )
      alpha_b = 2.07 * k_eq * sound_horizon * (1 + r_drag) ** (-0.75) * alpha_b_G

      beta_node = 8.41 * om0h2 ** 0.435
      beta_b = 0.5 + f_baryon + (3 - 2 * f_baryon) * sqrt ((17.2 * om0h2) ** 2 + 1.0)

      q = k / (13.41 * k_eq)
      ks = k * sound_horizon

      t_c_ln_beta = log (exp 1 + 1.8 * beta_c * q)
      t_c_ln_nobeta = log (exp 1 + 1.8 * q)
      t_c_C_alpha = 14.2 / alpha_c + 386.0 / (1 + 69.9 * q ** 1.08)
      t_c_C_noalpha = 14.2 + 386.0 / (1 + 69.9 * q ** 1.08)

      t_c_f = 1 / (1 + (ks / 5.4) ** 4)

      f = (\a b -> a / (a + b * q ** 2))

      t_c =
        t_c_f * f t_c_ln_beta t_c_C_noalpha
          + (1 - t_c_f) * f t_c_ln_beta t_c_C_alpha

      s_tilde = sound_horizon * (1 + (beta_node / ks) ** 3) ** (-1 / 3)
      ks_tilde = k * s_tilde

      t_b_T0 = f t_c_ln_nobeta t_c_C_noalpha
      t_b_1 = t_b_T0 / (1 + (ks / 5.2) ** 2)
      t_b_2 = alpha_b / (1 + (beta_b / ks) ** 3) * exp (-(k / k_silk) ** 1.4)
      t_b = sinc (ks_tilde / pi) * (t_b_1 + t_b_2)

      transfer = f_baryon * t_b + (1 - f_baryon) * t_c
   in transfer

powerSpectrumEisensteinHu :: ReferenceCosmology -> Redshift -> Wavenumber -> Double
powerSpectrumEisensteinHu cosmology z k =
  let (h0, om0, ob0, tcmb0, gn, as, ns, _) = unpackCosmology cosmology
      transfer = transferEinsensteinHu k cosmology

      om = om0 * (1 + z) ** (-3)
      ode = 1 - om

      expr = 2.5 * om / (1 + z)
      dz = expr / ((om ** (4.0 / 7.0) - ode + (1 + 0.5 * om) * (1.0 + ode / 70.0)))
      power_spectrum =
        sqrt (dz)
          * 1e-9
          * as
          * (2 * (k / (h0 / 299792.458)) ** 2 / (5 * om0)) ** 2
          * (transfer) ** 2
          * (k / 0.02) ** (ns - 1)
          * (2 * pi ** 2 / k ** 3)
   in power_spectrum
