use crate::constants::{GRAVITATIONAL, PI, SOLAR_MASS, SOLAR_RADIUS};
use crate::universe::particles::{Planet, Star};
use serde::{Deserialize, Serialize};

// Order unity parameters (Mustill & Villaver, 2012)
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Zahn {
    f_prime: f64,
    c_f: f64,
    gamma_f: f64,
}

impl Default for Zahn {
    fn default() -> Self {
        // Values adopted from Mustill & Villaver (2012)
        Self {
            f_prime: 9. / 5., // multiplication factor for the dissipation
            c_f: 1.,
            gamma_f: 2., // exponent for the frequency dependence of the viscosity
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum Equilibrium {
    #[default]
    Disabled,
    SigmaBarStar(f64),
    Zahn(Zahn),
    Evolution,
}

impl Equilibrium {
    pub fn tidal_quality(&self, star: &Star, planet: &Planet) -> f64 {
        match self {
            // When tides are disabled the tidal quality factor would be 1 / 0.
            // We set the tidal_quality to infinity such that 1 / infinity == 0 == disabled.
            Equilibrium::Disabled => f64::INFINITY,
            Equilibrium::SigmaBarStar(sigma_bar_star) => {
                Self::tidal_quality_sigma_bar_star(star, *sigma_bar_star)
            }
            Equilibrium::Zahn(zahn) => {
                Self::tidal_quality_zahn(star, planet, zahn.f_prime, zahn.c_f, zahn.gamma_f)
            }
            Equilibrium::Evolution => {
                todo!()
            }
        }
    }

    // Equilibrium tide
    // Normalization constant for equilibrium tide (Bolmont & Mathis,  2016,  Eq. 8)
    fn tidal_quality_sigma_bar_star(star: &Star, sigma_bar_star: f64) -> f64 {
        // Epsilon to ensure that equilibrium_tide_quality_factors stays finite
        let epsilon_secure = 1.0E-10;

        let normalisation_constant = sqrt!(GRAVITATIONAL / (SOLAR_MASS * SOLAR_RADIUS.powi(7)));
        // Tidal quality factors for the equilibrium tide
        // This is Eq. 22 of Benbakoura et al. 2019
        // omitting the factor 3/2 (which is due to a typo in the paper)

        GRAVITATIONAL
            / (abs!(star.tidal_frequency + epsilon_secure)
                * sigma_bar_star
                * normalisation_constant
                * star.radius.powi(5))
    }

    // Equilibrium tide with Zahn prescription
    // following the parametrisation (Mustill & Villaver, 2012, Eq. 1)
    fn tidal_quality_zahn(
        star: &Star,
        planet: &Planet,
        f_prime: f64,
        c_f: f64,
        gamma_f: f64,
    ) -> f64 {
        let f2 = f_prime
            * min!(
                1_f64,
                ((2. * PI) / (2. * planet.mean_motion * c_f * star.convective_turnover_time))
                    .powf(gamma_f)
            );
        let k2 = 1. / 27. / star.convective_turnover_time
            * (1. - star.radiative_mass / star.mass)
            * (star.mass + planet.mass)
            / star.mass
            / planet.mean_motion
            * (star.radius / planet.semi_major_axis).powi(3)
            * (2. * f2);

        3. / 2. / k2
    }
}

// References:
// Bolmont & Mathis 2016, https://doi.org/10.1007/s10569-016-9690-3
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Mustill & Villaver 2012, http://doi.org/10.1088/0004-637X/761/2/121
