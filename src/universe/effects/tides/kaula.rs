use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

mod love_number;
mod polynomials;

use love_number::{LoveNumber, ParticleComposition};
use polynomials::Polynomials;

use crate::universe::particles::ParticleT;
use crate::utils::{factorial, kronecker_delta};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
#[serde(default)]
#[serde(deny_unknown_fields)]
pub struct Kaula {
    pub(crate) particle_type: ParticleComposition,
    // Inclination and eccentricty polynomials
    #[serde(skip)]
    polynomials: Polynomials,
    #[serde(skip)]
    love_number: LoveNumber,
}

impl Kaula {
    pub fn solid_file(&self) -> Option<&PathBuf> {
        match self.particle_type {
            ParticleComposition::Solid { ref solid_file }
            | ParticleComposition::Atmosphere { ref solid_file, .. }
            | ParticleComposition::SolidAtmosphere { ref solid_file, .. }
            | ParticleComposition::SolidOcean { ref solid_file, .. }
            | ParticleComposition::SolidAtmosphereOcean { ref solid_file, .. } => Some(solid_file),
            _ => None,
        }
    }

    pub fn ocean_file(&self) -> Option<&PathBuf> {
        match self.particle_type {
            ParticleComposition::SolidOcean { ref ocean_file, .. }
            | ParticleComposition::SolidAtmosphereOcean { ref ocean_file, .. } => Some(ocean_file),
            _ => None,
        }
    }

    pub fn interpolate_dir(&self) -> Option<&PathBuf> {
        match self.particle_type {
            ParticleComposition::Interpolate {
                ref interpolate_dir,
                ..
            }
            | ParticleComposition::InterpolateAtmosphere {
                ref interpolate_dir,
                ..
            } => Some(interpolate_dir),
            _ => None,
        }
    }

    pub fn initialise_love_number_solid(&mut self, love_solid: &[Vec<f64>]) {
        self.love_number
            .imaginary_solid
            .init(&love_solid[0], &love_solid[1]);
        self.love_number
            .real_solid
            .init(&love_solid[0], &love_solid[2]);
    }

    pub fn initialise_love_number_ocean(&mut self, love_ocean: &[Vec<f64>]) {
        self.love_number
            .imaginary_oceanic
            .init(&love_ocean[0], &love_ocean[1]);
    }

    pub(crate) fn refresh(
        &mut self,
        time: f64,
        planet: &impl ParticleT,
        star: &impl ParticleT,
    ) -> Result<()> {
        self.polynomials
            .refresh_cache(planet.eccentricity(), planet.inclination());
        self.love_number
            .refresh_cache(time, planet, star, &self.particle_type)?;
        Ok(())
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the semi-major-axis derivative.
    // Boue & Efroimksy (2019) Eq 116 and Revol et al. (2023) Eq A.1.
    pub(crate) fn summation_of_longitudinal_modes_semi_major_axis(&self) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut sum_over_m = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut sum_over_p = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut sum_over_q = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    let q_factor = p_factor + (q as f64 - 7.);
                    let imk2 = self.love_number.imaginary(m, p, q);
                    sum_over_q += q_factor * q_val * imk2;
                }
                sum_over_p += m_val[p] * sum_over_q;
            }
            sum_over_m +=
                sum_over_p * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        sum_over_m
    }

    // Summation over longitudinal modes m for the computation of the spin derivative.
    // Boue & Efroimksy (2019) Eq 123 and Revol et al. (2023) Eq A.3
    pub(crate) fn summation_of_longitudinal_modes_spin(&self) -> f64 {
        self.sum_over_m_imaginary_mfactor(
            &self.polynomials.eccentricity_2pq_squared,
            &self.polynomials.inclination_2mp_squared,
        )
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the eccentricity derivative.
    // Boue & Efroimksy (2019) Eq 117 and Revol et al. (2023) Eq A.3
    pub(crate) fn summation_of_longitudinal_modes_eccentricity(&self, eccentricity: f64) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut sum_over_m = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut sum_over_p = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut sum_over_q = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    let q_factor = p_factor + (q as f64 - 7.);
                    let term = q_factor * sqrt!(1.0 - eccentricity.powi(2)) - p_factor;
                    let imk2 = self.love_number.imaginary(m, p, q);
                    sum_over_q += imk2 * q_val * term;
                }
                sum_over_p += sum_over_q * m_val[p];
            }
            sum_over_m +=
                sum_over_p * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        sum_over_m
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the inclination derivative.
    // by Boue & Efroimksy (2019) Eq 118 and Revol et al. (2023) Eq A.7
    pub(crate) fn summation_of_longitudinal_modes_inclination(
        &self,
        inclination: f64,
        term1: f64,
        term3: f64,
    ) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut sum_over_m = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut sum_over_p = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut sum_over_q = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    let imk2 = self.love_number.imaginary(m, p, q);
                    sum_over_q += imk2 * q_val;
                }
                sum_over_p += sum_over_q
                    * m_val[p]
                    * (term1 * (m as f64 * cos!(inclination) - p_factor)
                        - ((p_factor * cos!(inclination) - m as f64) * term3));
            }
            sum_over_m +=
                sum_over_p * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        sum_over_m
    }

    fn summation_of_longitudinal_modes_triple_common(
        &self,
        term1: f64,
        term2: f64,
        term3: f64,
    ) -> f64 {
        (self.sum_over_m_real(
            &self.polynomials.eccentricity_2pq_squared,
            &self.polynomials.inclination_2mp_squared_derivative,
        ) * term1
            * 0.5)
            + (self.sum_over_m_imaginary_mfactor(
                &self.polynomials.eccentricity_2pq_squared,
                &self.polynomials.inclination_2mp_squared,
            ) * term2)
            + (self.sum_over_m_imaginary_pfactor(
                &self.polynomials.eccentricity_2pq_squared,
                &self.polynomials.inclination_2mp_squared,
            ) * term3)
    }

    // Summation over longitudinal modes m for the computation of the longitude of ascending node derivative.
    // Boue & Efroimksy (2019) Eq 121 and Revol et al. (2023) Eq A.9
    pub(crate) fn summation_of_longitudinal_modes_longitude_ascending_node(
        &self,
        term1: f64,
        term2: f64,
        term3: f64,
    ) -> f64 {
        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the spin axis inclination derivative.
    // Boue & Efroimksy (2019) Eq 122 and Revol et al. (2023) Eq A.12
    pub(crate) fn summation_of_longitudinal_modes_spin_axis_inclination(
        &self,
        longitude_ascending_node: f64,
        inclination: f64,
    ) -> f64 {
        let term1 = -sin!(longitude_ascending_node);
        let term2 = cos!(longitude_ascending_node) * 1. / tan!(inclination);
        let term3 = -(cos!(longitude_ascending_node) / sin!(inclination));

        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the eccentricity dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    pub(crate) fn summation_of_longitudinal_modes_pericentre_eccentricity(
        &self,
        eccentricity: f64,
    ) -> f64 {
        if eccentricity == 0. {
            0.0
        } else {
            self.sum_over_m_real(
                &self.polynomials.eccentricity_2pq_squared_derivative,
                &self.polynomials.inclination_2mp_squared,
            ) * 0.5
        }
    }

    // Summation over longitudinal modes m for the computation of the inclination dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    pub(crate) fn summation_of_longitudinal_modes_pericentre_inclination(
        &self,
        inclination: f64,
        spin_inclination: f64,
    ) -> f64 {
        if (inclination == 0.) || (spin_inclination == 0.) {
            0.0
        } else {
            self.sum_over_m_real(
                &self.polynomials.eccentricity_2pq_squared,
                &self.polynomials.inclination_2mp_squared_derivative,
            ) * 0.5
        }
    }

    // Iteration over the provided 2D arrays (outer 3x15 and inner 3x3), summing the contents of:
    // (love_number(m, p, q) * inner[p][q]) * outer[m]p] * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0))
    fn sum_over_m_real(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3]) -> f64 {
        let mut m_sum = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut p_sum = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut q_sum = 0.0;
                for (q, q_val) in p_val.iter().enumerate() {
                    q_sum += q_val * self.love_number.real(m, p, q);
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    fn sum_over_m_imaginary_pfactor(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3]) -> f64 {
        let mut m_sum = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut p_sum = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut q_sum = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    q_sum += q_val * self.love_number.imaginary(m, p, q) * p_factor;
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }

    // Precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    fn sum_over_m_imaginary_mfactor(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3]) -> f64 {
        let mut m_sum = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut p_sum = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut q_sum = 0.0;
                for (q, q_val) in p_val.iter().enumerate() {
                    q_sum += q_val * self.love_number.imaginary(m, p, q);
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum
                * (factorial(2 - m) / factorial(2 + m))
                * (2. - kronecker_delta(m, 0))
                * (m as f64);
        }
        m_sum
    }
}

#[cfg(test)]
pub mod tests;
