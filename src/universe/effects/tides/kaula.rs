use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

mod love_number;
mod polynomials;

use love_number::{LoveNumber, ParticleComposition};
use polynomials::Polynomials;

use crate::universe::particles::{ParticleT, Planet};

// Upper and lower bound for the m, p, q summation.
// Calculated at each time step, based on inclination and eccentricity.
struct Mpq {
    m_min: u8,
    m_max: u8,
    p_min: u8,
    p_max: u8,
    q_min: u8,
    q_max: u8,
}

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

    // cache
    #[serde(skip)]
    sum_over_m_real_2pq_2mp_dt: f64,
    #[serde(skip)]
    sum_over_m_real_2pq_dt_2mp: f64,
    #[serde(skip)]
    sum_over_m_imaginary_mfactor: f64,
    #[serde(skip)]
    sum_over_m_imaginary_pfactor: f64,
    #[serde(skip)]
    sum_over_m_imaginary_qfactor: f64,
    #[serde(skip)]
    sum_over_m_imaginary_eccentricity: f64,
    #[serde(skip)]
    sum_over_m_imaginary_inclination: f64,
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

    fn bound_q_by_eccentricity(eccentricity: f64) -> (u8, u8) {
        match () {
            // Select the order of the summation q over the eccentricity function G_lpq
            () if (eccentricity > 0.30) => (0, 15), // q: -7 <= q <= 7
            () if (eccentricity > 0.25) => (1, 14), // q: -6 <= q <= 6
            () if (eccentricity > 0.20) => (2, 13), // q: -5 <= q <= 5
            () if (eccentricity > 0.15) => (3, 12), // q: -4 <= q <= 4
            () if (eccentricity > 1e-8) => (5, 10), // q: -2 <= q <= 2
            () if (eccentricity > 0.0) => (6, 9),   // q: -1 <= q <= 1
            () if (eccentricity == 0.0) => (7, 8),  // q:  0 <= q <= 0
            () => unreachable!("eccentricity cannot be negative."),
        }
    }

    // All the calculations using the polynomials and love number are performed here
    // and stored in the `sum_over_xxx` caches.
    // The caches are used when the derivitaves are calculated for each keplerian element.
    pub(crate) fn refresh(
        &mut self,
        time: f64,
        planet: &impl ParticleT,
        star: &impl ParticleT,
    ) -> Result<()> {
        self.polynomials
            .refresh_cache(planet.eccentricity(), planet.inclination());
        let (q_min, q_max) = Self::bound_q_by_eccentricity(planet.eccentricity());

        if planet.inclination() <= 1e-8 {
            // If inclination is close to zero, only compute
            // m = 0, p = 1 and m = 2, p = 0
            let mpq_01q = Mpq {
                m_min: 0,
                m_max: 1,
                p_min: 1,
                p_max: 2,
                q_min,
                q_max,
            };

            let mpq_20q = Mpq {
                m_min: 2,
                m_max: 3,
                p_min: 0,
                p_max: 1,
                q_min,
                q_max,
            };

            self.love_number
                .refresh_cache(time, planet, star, &self.particle_type, &mpq_01q)?;

            self.love_number
                .refresh_cache(time, planet, star, &self.particle_type, &mpq_20q)?;

            self.sum_over_m_imaginary_mfactor = self.sum_over_m_imaginary_mfactor(&mpq_01q)
                + self.sum_over_m_imaginary_mfactor(&mpq_20q);
            self.sum_over_m_imaginary_qfactor = self.sum_over_m_imaginary_qfactor(&mpq_01q)
                + self.sum_over_m_imaginary_qfactor(&mpq_20q);

            if planet.eccentricity() != 0.0 {
                self.sum_over_m_imaginary_eccentricity = self
                    .sum_over_m_imaginary_eccentricity(planet, &mpq_01q)
                    + self.sum_over_m_imaginary_eccentricity(planet, &mpq_20q);
                self.sum_over_m_real_2pq_dt_2mp = self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared_derivative,
                    &self.polynomials.inclination_2mp_squared,
                    &mpq_01q,
                ) + self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared_derivative,
                    &self.polynomials.inclination_2mp_squared,
                    &mpq_20q,
                );
            }

            if planet.inclination() != 0.0 {
                self.sum_over_m_imaginary_pfactor = self.sum_over_m_imaginary_pfactor(&mpq_01q)
                    + self.sum_over_m_imaginary_pfactor(&mpq_20q);
                self.sum_over_m_real_2pq_2mp_dt = self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared,
                    &self.polynomials.inclination_2mp_squared_derivative,
                    &mpq_01q,
                ) + self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared,
                    &self.polynomials.inclination_2mp_squared_derivative,
                    &mpq_20q,
                );
            }

            if sin!(planet.inclination()) != 0.0 {
                self.sum_over_m_imaginary_inclination = self
                    .sum_over_m_imaginary_inclination(planet, &mpq_01q)
                    + self.sum_over_m_imaginary_inclination(planet, &mpq_20q);
            }
        } else {
            let mpq = Mpq {
                m_min: 0,
                m_max: 3,
                p_min: 0,
                p_max: 3,
                q_min,
                q_max,
            };

            self.love_number
                .refresh_cache(time, planet, star, &self.particle_type, &mpq)?;

            self.sum_over_m_imaginary_mfactor = self.sum_over_m_imaginary_mfactor(&mpq);
            self.sum_over_m_imaginary_qfactor = self.sum_over_m_imaginary_qfactor(&mpq);

            if planet.eccentricity() != 0.0 {
                self.sum_over_m_imaginary_eccentricity =
                    self.sum_over_m_imaginary_eccentricity(planet, &mpq);
                self.sum_over_m_real_2pq_dt_2mp = self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared_derivative,
                    &self.polynomials.inclination_2mp_squared,
                    &mpq,
                );
            }

            if planet.inclination() != 0.0 {
                self.sum_over_m_imaginary_pfactor = self.sum_over_m_imaginary_pfactor(&mpq);
                self.sum_over_m_real_2pq_2mp_dt = self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared,
                    &self.polynomials.inclination_2mp_squared_derivative,
                    &mpq,
                );
            }

            if sin!(planet.inclination()) != 0.0 {
                self.sum_over_m_imaginary_inclination =
                    self.sum_over_m_imaginary_inclination(planet, &mpq);
            }
        }

        Ok(())
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the semi-major-axis derivative.
    // Boue & Efroimksy (2019) Eq 116 and Revol et al. (2023) Eq A.1.
    pub(crate) fn summation_of_longitudinal_modes_semi_major_axis(&self) -> f64 {
        self.sum_over_m_imaginary_qfactor
    }

    // Summation over longitudinal modes m for the computation of the spin derivative.
    // Boue & Efroimksy (2019) Eq 123 and Revol et al. (2023) Eq A.3
    pub(crate) fn summation_of_longitudinal_modes_spin(&self) -> f64 {
        self.sum_over_m_imaginary_mfactor
    }

    pub(crate) fn summation_of_longitudinal_modes_eccentricity(&self) -> f64 {
        self.sum_over_m_imaginary_eccentricity
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the eccentricity derivative.
    // Boue & Efroimksy (2019) Eq 117 and Revol et al. (2023) Eq A.3
    fn sum_over_m_imaginary_eccentricity(&self, planet: &impl ParticleT, mpq: &Mpq) -> f64 {
        let semi_minor_axis_ratio = sqrt!(1. - planet.eccentricity().powi(2));
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut m_sum = 0.0;
        for (m, m_val) in outer
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
        {
            let mut p_sum = 0.0;
            for (p, p_val) in inner
                .iter()
                .enumerate()
                .take(mpq.p_max.into())
                .skip(mpq.p_min.into())
            {
                let mut q_sum = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val
                    .iter()
                    .enumerate()
                    .take(mpq.q_max.into())
                    .skip(mpq.q_min.into())
                {
                    let q_factor = p_factor + (q as f64 - 7.);
                    let term = q_factor * semi_minor_axis_ratio - p_factor;
                    let imk2 = self.love_number.imaginary(m, p, q);
                    q_sum += imk2 * q_val * term;
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }

    pub(crate) fn summation_of_longitudinal_modes_inclination(&self) -> f64 {
        self.sum_over_m_imaginary_inclination
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the inclination derivative.
    // by Boue & Efroimksy (2019) Eq 118 and Revol et al. (2023) Eq A.7
    fn sum_over_m_imaginary_inclination(&self, planet: &impl ParticleT, mpq: &Mpq) -> f64 {
        let semi_minor_axis_ratio = sqrt!(1. - planet.eccentricity().powi(2));

        let term1 = (planet.reduced_mass()
            * planet.mean_motion().powi(2)
            * planet.semi_major_axis().powi(2))
            / (planet.moment_of_inertia() * planet.spin());
        let term3 = planet.mean_motion() / semi_minor_axis_ratio;
        let cos_inc = cos!(planet.inclination());

        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut m_sum = 0.0;
        for (m, m_val) in outer
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
        {
            let mut p_sum = 0.0;
            for (p, p_val) in inner
                .iter()
                .enumerate()
                .take(mpq.p_max.into())
                .skip(mpq.p_min.into())
            {
                let mut q_sum = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val
                    .iter()
                    .enumerate()
                    .take(mpq.q_max.into())
                    .skip(mpq.q_min.into())
                {
                    let imk2 = self.love_number.imaginary(m, p, q);
                    q_sum += imk2 * q_val;
                }
                p_sum += q_sum
                    * m_val[p]
                    * (term1 * (m as f64 * cos_inc - p_factor)
                        - ((p_factor * cos_inc - m as f64) * term3));
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }

    fn summation_of_longitudinal_modes_triple_common(
        &self,
        term1: f64,
        term2: f64,
        term3: f64,
    ) -> f64 {
        (self.sum_over_m_real_2pq_2mp_dt * term1 * 0.5)
            + (self.sum_over_m_imaginary_mfactor * term2)
            + (self.sum_over_m_imaginary_pfactor * term3)
    }

    // Summation over longitudinal modes m for the computation of the longitude of ascending node derivative.
    // Boue & Efroimksy (2019) Eq 121 and Revol et al. (2023) Eq A.9
    pub(crate) fn summation_of_longitudinal_modes_longitude_ascending_node(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term1 = (1. / (planet.moment_of_inertia * planet.spin * planet.tan_inc))
            - (planet.cos_lan / (planet.moment_of_inertia * planet.spin * planet.tan_spin_inc))
            + (1.
                / (planet.reduced_mass
                    * planet.mean_motion
                    * planet.semi_major_axis.powi(2)
                    * planet.semi_minor_axis_ratio
                    * planet.sin_inc));

        let term2 = -(planet.sin_lan * cotan!(planet.inclination))
            / (planet.moment_of_inertia * planet.spin * planet.tan_spin_inc);

        let term3 = planet.sin_lan
            / (planet.moment_of_inertia * planet.spin * planet.tan_spin_inc * planet.sin_inc);

        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the spin axis inclination derivative.
    // Boue & Efroimksy (2019) Eq 122 and Revol et al. (2023) Eq A.12
    pub(crate) fn summation_of_longitudinal_modes_spin_axis_inclination(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term1 = -planet.sin_lan;
        let term2 = planet.cos_lan / planet.tan_inc;
        let term3 = -(planet.cos_lan / planet.sin_inc);

        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the eccentricity dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    pub(crate) fn summation_of_longitudinal_modes_pericentre_eccentricity(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term2 = planet.semi_minor_axis_ratio
            / (planet.mean_motion
                * planet.semi_major_axis.powi(2)
                * planet.eccentricity
                * planet.reduced_mass);
        self.sum_over_m_real_2pq_dt_2mp * 0.5 * term2
    }

    // Summation over longitudinal modes m for the computation of the inclination dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    pub(crate) fn summation_of_longitudinal_modes_pericentre_inclination(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term1 = -((1. / (planet.moment_of_inertia * planet.spin * planet.sin_inc))
            + (1.
                / (planet.mean_motion
                    * planet.semi_major_axis.powi(2)
                    * planet.semi_minor_axis_ratio
                    * planet.tan_inc
                    * planet.reduced_mass)));
        self.sum_over_m_real_2pq_2mp_dt * 0.5 * term1
    }

    // Iteration over the provided 2D arrays (outer 3x15 and inner 3x3), summing the contents of:
    // (love_number(m, p, q) * inner[p][q]) * outer[m]p] * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0))
    fn sum_over_m_real(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3], mpq: &Mpq) -> f64 {
        let mut m_sum = 0.0;
        for (m, m_val) in outer
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
        {
            let mut p_sum = 0.0;
            for (p, p_val) in inner
                .iter()
                .enumerate()
                .take(mpq.p_max.into())
                .skip(mpq.p_min.into())
            {
                let mut q_sum = 0.0;
                for (q, q_val) in p_val
                    .iter()
                    .enumerate()
                    .take(mpq.q_max.into())
                    .skip(mpq.q_min.into())
                {
                    let rek2 = self.love_number.real(m, p, q);
                    q_sum += rek2 * q_val;
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
    fn sum_over_m_imaginary_pfactor(&self, mpq: &Mpq) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut m_sum = 0.0;
        for (m, m_val) in outer
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
        {
            let mut p_sum = 0.0;
            for (p, p_val) in inner
                .iter()
                .enumerate()
                .take(mpq.p_max.into())
                .skip(mpq.p_min.into())
            {
                let mut q_sum = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val
                    .iter()
                    .enumerate()
                    .take(mpq.q_max.into())
                    .skip(mpq.q_min.into())
                {
                    let imk2 = self.love_number.imaginary(m, p, q);
                    q_sum += imk2 * q_val * p_factor;
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }

    // Precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    fn sum_over_m_imaginary_mfactor(&self, mpq: &Mpq) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;
        // Skip over the case of m = 0.
        let m_min = max!(1, mpq.m_min);
        let mut m_sum = 0.0;
        for (m, m_val) in outer
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(m_min.into())
        {
            let mut p_sum = 0.0;
            for (p, p_val) in inner
                .iter()
                .enumerate()
                .take(mpq.p_max.into())
                .skip(mpq.p_min.into())
            {
                let mut q_sum = 0.0;
                for (q, q_val) in p_val
                    .iter()
                    .enumerate()
                    .take(mpq.q_max.into())
                    .skip(mpq.q_min.into())
                {
                    let imk2 = self.love_number.imaginary(m, p, q);
                    q_sum += imk2 * q_val;
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

    // Precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    fn sum_over_m_imaginary_qfactor(&self, mpq: &Mpq) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut m_sum = 0.0;
        for (m, m_val) in outer
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
        {
            let mut p_sum = 0.0;
            for (p, p_val) in inner
                .iter()
                .enumerate()
                .take(mpq.p_max.into())
                .skip(mpq.p_min.into())
            {
                let mut q_sum = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val
                    .iter()
                    .enumerate()
                    .take(mpq.q_max.into())
                    .skip(mpq.q_min.into())
                {
                    let q_factor = p_factor + (q as f64 - 7.);
                    let imk2 = self.love_number.imaginary(m, p, q);
                    q_sum += imk2 * q_val * q_factor;
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }
}

// Precomputed factorials as f64 for 0..=4
const FACTORIALS: &[f64] = &[1.0, 1.0, 2.0, 6.0, 24.0];

pub(crate) fn factorial(x: usize) -> f64 {
    FACTORIALS[x]
}

// x == y -> 1.0
// x != y -> 0.0
pub(crate) fn kronecker_delta(x: usize, y: usize) -> f64 {
    f64!(x == y)
}

#[cfg(test)]
pub mod tests;
