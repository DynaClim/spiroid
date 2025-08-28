use crate::universe::particles::{Planet, Star};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum Inertial {
    #[default]
    Disabled,
    FrequencyAveraged,
    FrequencyDependant,
}

impl Inertial {
    pub fn tidal_quality(self, star: &Star, _planet: &Planet) -> f64 {
        match self {
            // When tides are disabled the tidal quality factor would be 1 / 0.
            // We set the tidal_quality to infinity such that 1 / infinity == 0 == disabled.
            Inertial::Disabled => f64::INFINITY,
            Inertial::FrequencyAveraged => Self::tidal_quality_frequency_averaged(star),
            Inertial::FrequencyDependant => {
                todo!()
            }
        }
    }

    // Dynamical tide
    // Equivalent Q' factor, tidal quality factor
    // This is the inverse of Eq. 4 of Bolmont & Mathis 2016,
    // which give the tidal quality factor Q' as a function of <D>w
    // (Eq. 1 of Mathis 2015, or Bolmont & Mathis 2016)
    // But here, keep in mind that the spin was removed out of <D>w to be then multiplied here
    // Q' = 3 / ( 2 * <D>w ) | Bolmont & Mathis
    // Q' = 3 / ( 2 * <D>w' * spin^2 ) | here
    fn tidal_quality_frequency_averaged(star: &Star) -> f64 {
        // Epsilon to ensure a smooth transition equilibrium / dynamical tide
        let epsilon_step = 1.0E-06;
        // dynamical_tide_quality_factors
        3. / (2. * star.dynamical_tide_dissipation * star.spin.powi(2)) * 2.
            / (1. + tanh!((star.tidal_frequency + 2. * star.spin) / epsilon_step))
    }
}

// References:
// Mathis 2015, https://doi.org/10.1051/0004-6361/201526472
// Bolmont & Mathis 2016, https://doi.org/10.1007/s10569-016-9690-3
