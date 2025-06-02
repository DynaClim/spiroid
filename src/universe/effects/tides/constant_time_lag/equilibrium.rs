use serde::{Deserialize, Serialize};
use crate::universe::particles::Star;
use crate::constants::{GRAVITATIONAL, SOLAR_MASS, SOLAR_RADIUS};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum Equilibrium {
    #[default]
    Disabled,
    SigmaBarStar {
        sigma_bar_star: f64, // dimensionless sigma_bar_star
    },
    Evolution
}

impl Equilibrium {
    pub fn tidal_quality(&self, star: &Star) -> f64 {
        match self {
            Equilibrium::Disabled => 0.0,
            Equilibrium::SigmaBarStar { sigma_bar_star } => {
                self.tidal_quality_sigma_bar_star(star, *sigma_bar_star)
            }
            Equilibrium::Evolution => {
                todo!()
            }
        }
    }

    // Equilibrium tide
    // Normalization constant for equilibrium tide (Bolmont & Mathis,  2016,  Eq. 8)
    pub fn tidal_quality_sigma_bar_star(&self, star: &Star, sigma_bar_star: f64) -> f64 {
        // Epsilon to ensure that equilibrium_tide_quality_factors stays finite
        let epsilon_secure = 1.0E-10;

        let normalisation_constant = sqrt!(GRAVITATIONAL / (SOLAR_MASS * SOLAR_RADIUS.powi(7)));
        // Tidal quality factors for the equilibrium tide
        // This is Eq. 22 of Benbakoura et al. 2019
        let equilibrium_tide_quality_factors = GRAVITATIONAL
            / (abs!(star.tidal_frequency + epsilon_secure)
                * sigma_bar_star
                * normalisation_constant
                * star.radius.powi(5));
                
        equilibrium_tide_quality_factors
    }
}