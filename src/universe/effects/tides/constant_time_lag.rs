use serde::{Deserialize, Serialize};
use crate::universe::particles::{Star, Planet};
use crate::constants::GRAVITATIONAL;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ConstantTimeLag {
    pub(crate) sigma_bar_star: f64, // dimensionless sigma_bar_star
}

impl ConstantTimeLag {
    
    pub fn new(sigma_bar_star: f64) -> Self {
        Self { sigma_bar_star }
    }

    // This is a re-write of Eq. 3 and 19 from Benbakoura et al. 2019
    // without the factors that are in the function semi_major_axis_13_div_2_derivative in physics.rs
    // The a^-6 is here to compensate the a^6 in physics.rs
    pub fn tidal_torque(&self, star: &Star, planet: &Planet) -> f64 {
        let tidal_quality = star.tidal_quality(self.sigma_bar_star);
        // Smoothing parameter when tidal frequency is 0
        let depth = 1E-08;
        -(9. / 4.)
            * planet.mass.powi(2)
            * GRAVITATIONAL
            * planet.semi_major_axis.powi(-6)
            * tanh!(star.tidal_frequency / depth)
            * star.radius.powi(5)
            / tidal_quality
    }
}