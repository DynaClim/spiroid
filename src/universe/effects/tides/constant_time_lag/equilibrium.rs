use serde::{Deserialize, Serialize};
use crate::universe::particles::{Star, Planet};
use crate::constants::{GRAVITATIONAL, SOLAR_MASS, SOLAR_RADIUS, PI};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum Equilibrium {
    #[default]
    Disabled,
    SigmaBarStar {
        sigma_bar_star: f64, // dimensionless sigma_bar_star
    },
    Zahn {
        f_prime: f64,
        c_f: f64,
        gamma_f: f64,
    },
    Evolution
}

impl Equilibrium {
    pub fn tidal_quality(&self, star: &Star, planet: &Planet) -> f64 {
        match self {
            Equilibrium::Disabled => f64::INFINITY,
            Equilibrium::SigmaBarStar { sigma_bar_star } => {
                self.tidal_quality_sigma_bar_star(star, *sigma_bar_star)
            }
            Equilibrium::Zahn { f_prime, c_f, gamma_f } => {
                self.tidal_quality_zahn(star, planet, *f_prime, *c_f, *gamma_f)
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
        // This is Eq. 22 of Benbakoura et al. 2019, omitting the factor 3/2 (which is due to a typo in the paper)
        let equilibrium_tide_quality_factors = GRAVITATIONAL
            / (abs!(star.tidal_frequency + epsilon_secure)
                * sigma_bar_star
                * normalisation_constant
                * star.radius.powi(5));
                
        equilibrium_tide_quality_factors
    }

    pub fn tidal_quality_zahn(&self, star: &Star, planet: &Planet, f_prime: f64, c_f: f64, gamma_f: f64) -> f64 {
        let f2 = f_prime * min!(1_f64, ((2. * PI) / (2. * planet.mean_motion * c_f * star.convective_turnover_time)).powf(gamma_f));
        let k2 = 1. / 27.
            / star.convective_turnover_time
            * (1. - star.radiative_mass / star.mass)
            * (star.mass + planet.mass) / star.mass
            / planet.mean_motion
            * (star.radius / planet.semi_major_axis).powi(3)
            * (2. * f2);

        3. / 2. / k2
    }
}