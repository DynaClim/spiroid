use serde::{Deserialize, Serialize};
use crate::universe::particles::Star;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum Inertial {
    #[default]
    Disabled,
    FrequencyAveraged,
    FrequencyDependant
}

impl Inertial {
    pub fn tidal_quality(&self, star: &Star) -> f64 {
        match self {
            Inertial::Disabled => 0.0,
            Inertial::FrequencyAveraged => {
                self.tidal_quality_frequency_averaged(star)
            }
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
    pub fn tidal_quality_frequency_averaged(&self, star: &Star) -> f64 {
        // Epsilon to ensure a smooth transition equilibrium / dynamical tide
        let epsilon_step = 1.0E-06;

        let dynamical_tide_quality_factors =
            3. / (2. * star.dynamical_tide_dissipation * star.spin.powi(2))
                / 0.5
                / (1. + tanh!((star.tidal_frequency + 2. * star.spin) / epsilon_step));

        dynamical_tide_quality_factors
    }
}