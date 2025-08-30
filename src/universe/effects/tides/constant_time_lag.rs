use crate::constants::GRAVITATIONAL;
use crate::universe::particles::{Planet, Star};
use serde::{Deserialize, Serialize};

pub(crate) mod equilibrium;
pub use equilibrium::Equilibrium;

pub(crate) mod inertial;
pub use inertial::Inertial;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub struct ConstantTimeLag {
    pub equilibrium: Equilibrium,
    pub inertial: Inertial,
}

impl ConstantTimeLag {
    // This is a re-write of Eq. 3 and 19 from Benbakoura et al. 2019
    // without the factors that are in the function semi_major_axis_13_div_2_derivative in physics.rs
    // The a^-6 is here to compensate the a^6 in physics.rs
    pub fn tidal_torque(&self, star: &Star, planet: &Planet) -> f64 {
        let total_dissipation = 1. / self.equilibrium.tidal_quality(star, planet)
            + 1. / self.inertial.tidal_quality(star, planet);

        let tidal_quality = 1. / total_dissipation;
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

#[cfg(test)]
pub mod tests;

// References:
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
