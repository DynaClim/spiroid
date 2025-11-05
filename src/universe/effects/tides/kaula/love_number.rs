use crate::constants::{
    GAS_CONSTANT, GRAVITATIONAL, PI, SECONDS_IN_DAY, SECONDS_IN_YEAR, SOLAR_LUMINOSITY,
};
use crate::universe::effects::tides::kaula::Mpq;
use crate::universe::particles::ParticleT;
use sci_file::Interpolator;

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use anyhow::Result;
use serde_big_array::BigArray;

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub enum ParticleComposition {
    #[default]
    None,
    Solid {
        solid_file: PathBuf,
    },
    Atmosphere {
        solid_file: PathBuf,
        thermal_tide_model: ThermalTideModel,
    },
    SolidAtmosphere {
        solid_file: PathBuf,
        thermal_tide_model: ThermalTideModel,
    },
    SolidOcean {
        solid_file: PathBuf,
        ocean_file: PathBuf,
    },
    SolidAtmosphereOcean {
        solid_file: PathBuf,
        ocean_file: PathBuf,
        thermal_tide_model: ThermalTideModel,
    },
    Interpolate {
        interpolate_dir: PathBuf,
    },
    InterpolateAtmosphere {
        interpolate_dir: PathBuf,
        thermal_tide_model: ThermalTideModel,
    },
}

// Real and imaginary love numbers calculated each timestep.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub(crate) struct LoveNumber {
    // Cache
    #[serde(with = "BigArray")]
    real: [f64; 57],
    #[serde(with = "BigArray")]
    imaginary: [f64; 57],

    pub(crate) imaginary_atmosphere: Interpolator<f64>,
    pub(crate) imaginary_oceanic: Interpolator<f64>,
    pub(crate) imaginary_solid: Interpolator<f64>,
    pub(crate) real_solid: Interpolator<f64>,
    love_interpolator: Vec<Interpolator<f64>>,
}

impl Default for LoveNumber {
    fn default() -> Self {
        Self {
            real: [0.0; 57],
            imaginary: [0.0; 57],
            imaginary_atmosphere: Interpolator::new(),
            imaginary_oceanic: Interpolator::new(),
            imaginary_solid: Interpolator::new(),
            real_solid: Interpolator::new(),
            love_interpolator: vec![Interpolator::new()],
        }
    }
}

impl LoveNumber {
    /// Fetches the real love number from the cache for index of tuple (m, p, q).
    pub(crate) fn real(&self, m: usize, p: usize, q: usize) -> f64 {
        // The cached data is stored in a 1D array, so the 3D coordinates are mapped to the 1D index.
        let index = Self::get_index(m, p, q);
        self.real[index]
    }

    /// Fetches the imaginary love number from the cache for index of tuple (m, p, q).
    pub(crate) fn imaginary(&self, m: usize, p: usize, q: usize) -> f64 {
        // The cached data is stored in a 1D array, so the 3D coordinates are mapped to the 1D index.
        let index = Self::get_index(m, p, q);
        self.imaginary[index]
    }

    #[allow(clippy::cast_possible_truncation)]
    #[allow(clippy::cast_sign_loss)]
    fn get_index(m: usize, p: usize, q: usize) -> usize {
        // Truncation not applicable since the values are in [0..3).
        let pq_fac = Self::pq_fac(p as u8, q as u8);
        // Sign loss not applicable since the values are in [0..18]
        // Convert pq_fac from -9..=9 to 0..18
        let pq_fac = (pq_fac + 9) as usize;

        // Convert from 2d index of q_fac, m to 1d index (19 is number of possible q_fac values)
        pq_fac + m * 19
    }

    /// Sets the real and imaginary love number into the cache for index of tuple (m, p, q).
    fn set_k2(&mut self, m: u8, p: u8, q: u8, real: f64, imaginary: f64) {
        let index = Self::get_index(m.into(), p.into(), q.into());
        self.real[index] = real;
        self.imaginary[index] = imaginary;
    }

    /// Clear the cache for the tidal frequencies in the range to be populated.
    fn clear_cache(&mut self, mpq: Mpq) {
        for q in mpq.q_min..mpq.q_max {
            for p in mpq.p_min..mpq.p_max {
                for m in mpq.m_min..mpq.m_max {
                    self.set_k2(m, p, q, 0.0, 0.0);
                }
            }
        }
    }

    fn pq_fac(p: u8, q: u8) -> i32 {
        let p = i32::from(p);
        // Convert q from 0..=14 to -7..=7
        let q = i32::from(q) - 7;

        2 - 2 * p + q
    }

    #[allow(clippy::cast_possible_truncation)]
    fn tidal_excitation_frequency_mode_sigma_2mpq(
        planet: &impl ParticleT,
        m: u8,
        p: u8,
        q: u8,
    ) -> f64 {
        let pq_fac = f64!(Self::pq_fac(p, q));
        let m = f64!(m);
        pq_fac * planet.mean_motion() - m * planet.spin()
    }

    /// Recomputes all the love number values.
    // Called at each time step to cache love numbers for that iteration, to prevent duplicate calculations.
    pub(crate) fn refresh_cache(
        &mut self,
        time: f64,
        planet: &impl ParticleT,
        star: &impl ParticleT,
        particle_type: &ParticleComposition,
        mpq: Mpq,
    ) -> Result<()> {
        let mut k2_re;
        let mut k2_im;
        self.clear_cache(mpq);

        // This internal cache acts as an allocation free hash table, where the index (key) is derived from (q_fac, m) and value is the k2
        // for the relevant tidal frequency.
        // If the array value is zero, the cache is assumed to not been filled. Tidal frequency of zero is not calculated,
        // so the (possible) false cache miss is insignificant.
        // There are 19 possible q_fac values (-9..=9) for each m value (0..=2)
        // This cache means only 57 calculations are done for the full 135 (q = 15 * p = 3 * m = 3) iterations of the loop.
        for q in mpq.q_min..mpq.q_max {
            for p in mpq.p_min..mpq.p_max {
                for m in mpq.m_min..mpq.m_max {
                    k2_re = self.real(m.into(), p.into(), q.into());
                    k2_im = self.imaginary(m.into(), p.into(), q.into());
                    if k2_re == 0.0 && k2_im == 0.0 {
                        // Cache miss, compute k2
                        let w_2lmpq =
                            Self::tidal_excitation_frequency_mode_sigma_2mpq(planet, m, p, q);
                        k2_re = self.real_part(w_2lmpq, particle_type)?;
                        k2_im = self.imaginary_part(time, w_2lmpq, planet, star, particle_type)?;
                        // Add to cache
                        self.set_k2(m, p, q, k2_re, k2_im);
                    }
                }
            }
        }
        Ok(())
    }

    // Select the correct love number calculation based on the composition of the planet.
    // TODO add additional solid love numbers for different particle types when they become available.
    fn real_part(&self, freq: f64, particle_type: &ParticleComposition) -> Result<f64> {
        match particle_type {
            ParticleComposition::None => {
                todo!();
            }
            _ => Ok(self.real_solid(freq)?),
        }
    }

    // Real part of love number is always negative.
    fn real_solid(&self, freq: f64) -> Result<f64> {
        let (_, real_k2) = self.real_solid.interpolate(abs!(freq))?;
        Ok(-real_k2)
    }

    // Select the correct love number calculation based on the composition of the planet.
    fn imaginary_part(
        &self,
        time: f64,
        freq: f64,
        planet: &impl ParticleT,
        star: &impl ParticleT,
        particle_type: &ParticleComposition,
    ) -> Result<f64> {
        match particle_type {
            ParticleComposition::None => {
                todo!();
            }
            ParticleComposition::Solid { .. } => self.imaginary_solid(freq),
            ParticleComposition::Atmosphere {
                thermal_tide_model, ..
            } => Ok(thermal_tide_model.imaginary_atmosphere(freq, planet, star)),
            ParticleComposition::SolidAtmosphere {
                thermal_tide_model, ..
            } => Ok(self.imaginary_solid(freq)?
                + thermal_tide_model.imaginary_atmosphere(freq, planet, star)),
            ParticleComposition::SolidOcean { .. } => {
                Ok(self.imaginary_solid(freq)? + self.imaginary_oceanic(freq)?)
            }
            ParticleComposition::SolidAtmosphereOcean {
                thermal_tide_model, ..
            } => Ok(self.imaginary_solid(freq)?
                + thermal_tide_model.imaginary_atmosphere(freq, planet, star)
                + self.imaginary_oceanic(freq)?),
            ParticleComposition::Interpolate { .. } => {
                self.love_number_interpolated_by_frequency(time, freq)
            }
            ParticleComposition::InterpolateAtmosphere {
                thermal_tide_model, ..
            } => Ok(self.love_number_interpolated_by_frequency(time, freq)?
                + thermal_tide_model.imaginary_atmosphere(freq, planet, star)),
        }
    }

    fn imaginary_solid(&self, freq: f64) -> Result<f64> {
        let (_, im_k2) = self.imaginary_solid.interpolate(abs!(freq))?;
        Ok(freq.signum() * im_k2)
    }

    fn imaginary_oceanic(&self, freq: f64) -> Result<f64> {
        let (_, im_k2) = self.imaginary_oceanic.interpolate(abs!(freq))?;
        Ok(im_k2)
    }

    // Love number data stored across multiple files 1.0, 1.1, 1.2, ..., 4.0
    // The number represents the giga-year
    // Convert the time to giga-years, then index into the vector to access the relevant data
    // e.g. time ~= 1.0 gigayears: (1 - 1) * 10 == 0, so vec[0] contains relevant data
    // e.g. time ~= 3.5 gigayears: (3.5 - 1) * 10 == 25, so vec[25] contains relevant data.
    fn love_number_interpolated_by_frequency(&self, time: f64, freq: f64) -> Result<f64> {
        // Find which section of the love number data files to use, based on the "giga-year" and convert it to an index
        #[allow(clippy::cast_possible_truncation)]
        #[allow(clippy::cast_sign_loss)]
        let index = (time / 1e9 / SECONDS_IN_YEAR) as usize * 10 - 10;
        let love_interpolator = &self.love_interpolator[index];
        if freq == 0.0 {
            Ok(0.0)
        } else {
            let (_, im_k2) = love_interpolator.interpolate(abs!(freq))?;
            Ok(freq.signum() * im_k2)
        }
    }
}

#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub enum ThermalTideModel {
    Analytic,
    Auclair {
        surface_temperature: f64,
        radiative_frequency: f64, // [omega] Radiative thermal frequency of the atmosphere (s^-1).
    },
    AuclairScaling {
        surface_pressure: f64, // Surface pressure (Pa)
    },
    Leconte {
        thermal_tide_amplitude: f64, // [q_0] Amplitude of the atmospheric quadrupole (Pa)
        radiative_frequency: f64, // [omega] Radiative thermal frequency of the atmosphere (s^-1).
    },
}

impl ThermalTideModel {
    fn imaginary_atmosphere(
        &self,
        freq: f64,
        planet: &impl ParticleT,
        star: &impl ParticleT,
    ) -> f64 {
        match self {
            ThermalTideModel::Analytic => Self::imaginary_atmosphere_analytic(freq, star, planet),
            ThermalTideModel::Auclair {
                surface_temperature,
                radiative_frequency,
            } => Self::imaginary_atmosphere_auclair(
                *surface_temperature,
                *radiative_frequency,
                freq,
                star,
                planet,
            ),
            ThermalTideModel::AuclairScaling { surface_pressure } => {
                Self::imaginary_atmosphere_auclair_scaling(*surface_pressure, freq, star, planet)
            }
            ThermalTideModel::Leconte {
                thermal_tide_amplitude,
                radiative_frequency,
            } => Self::imaginary_atmosphere_leconte(
                *thermal_tide_amplitude,
                *radiative_frequency,
                freq,
                star,
                planet,
            ),
        }
    }

    // Auclair-desrotour 2017a Eq. 173
    // Values of physical parameter table 1
    fn imaginary_atmosphere_analytic(
        freq: f64,
        star: &impl ParticleT,
        planet: &impl ParticleT,
    ) -> f64 {
        // Related to the first adiabatic exponent of the gas.
        let kappa = 0.286;
        let surface_temperature = 737.;
        // Radiative thermal frequency of the atmosphere (s^-1).
        let omega = 2_f64 * 3.77e-7;
        // Effective fraction of power absorbed by the atmosphere.
        let epsilon = 0.04;
        // Shape factor defined on the spatial distribution of tidal heat sources.
        let alpha = 0.2;

        let ra = 191.; // Specific gas constant
        // Imaginary part of the thermal Love number.
        // Auclair-Desrotour 2017b Eq. 5 + 6
        -(epsilon * alpha * star.luminosity() * planet.semi_major_axis() * kappa)
            / (5. * star.mass() * ra * surface_temperature * planet.radius())
            * freq
            / (freq.powi(2) + omega.powi(2))
    }

    // Thermal Love number
    // Auclair-Desrotour et al. (2017b) based on Equation 5 and 6
    fn imaginary_atmosphere_auclair(
        surface_temperature: f64,
        omega: f64,
        freq: f64,
        star: &impl ParticleT,
        planet: &impl ParticleT,
    ) -> f64 {
        // Related to the first adiabatic exponent of the gas.
        let kappa = 0.286;
        // Effective fraction of power absorbed by the atmosphere.
        let epsilon = 0.04;
        // Shape factor defined on the spatial distribution of tidal heat sources.
        let alpha = 0.14;

        // Efficiency of dynamical (viscous) coupling between atmospheric layers
        let beta = 1.0;
        // Mean molar mass of the atmosphere
        let m_a = 43.45e-3;

        // Specific gas constant.
        let r_a = GAS_CONSTANT / m_a;
        let factor = -(4.0 / 32.0)
            * kappa
            * beta
            * alpha
            * epsilon
            * star.luminosity()
            * planet.semi_major_axis()
            / (r_a * surface_temperature * star.mass() * planet.radius());
        // Rescaled Radiative frequency
        let w_0 = omega * (star.luminosity() / SOLAR_LUMINOSITY).powf(0.75);
        // Maxwell-like frequency dependence
        let q_a = freq / (freq.powi(2) + w_0.powi(2));

        factor * q_a
    }

    // Thermal tide scaling model
    // Auclair-Desrotour et al. (2019) Sec. 5.3
    fn imaginary_atmosphere_auclair_scaling(
        surface_pressure: f64,
        freq: f64,
        star: &impl ParticleT,
        planet: &impl ParticleT,
    ) -> f64 {
        // Avoid division by zero NaN.
        if freq == 0. {
            0.
        } else {
            // scaling model for the atmospheric love number, Eq. 49
            let a1 = 0.734;
            let a2 = -1.;
            let b1 = 0.171;
            let b2 = -0.031;
            let btrans = -0.02;
            let d1 = 0.01;
            let d2 = 0.023;
            let chi1 = -0.277;
            let chi2 = 0.29;
            // Convert from Pa to Bar
            let surface_pressure_bar = surface_pressure / 100_000.;
            // Scaled thermal time scale and amplitude (using scaling formulation with fixed sma a = a_venus) Eq. 44 and 45
            // Scaled pressure Eq. 44
            let q_0 = 10_f64.powf(0.48 * log10!(surface_pressure_bar) + 2.87);
            // Scaled time-scale Eq. 45
            let tau_0 = 10.0_f64.powf(0.3 * log10!(surface_pressure_bar) + 0.038);

            // Scaled frequency Eq. 46
            let sigma = freq * SECONDS_IN_DAY;
            let chi = log10!(abs!(tau_0 * sigma));
            // Activation funtions Eq. 48
            let f1 = 1.0 / (1.0 + ((chi - chi1) / d1).exp());
            let f2 = 1.0 / (1.0 + (-(chi - chi2) / d2).exp());

            // Parametrized function Eq. 24 and 47
            let f_par = (a1 * chi + b1) * f1 + (a2 * chi + b2) * f2 + btrans * (1.0 - f1 - f2);

            // Imaginary part of the spherical harmonic of surface pressure variations Eq. 46
            let imaginary_delta_pressure_2 = q_0 * 10.0_f64.powf(f_par) * freq.signum();

            // Imaginary tidal love number associated with conversion factor derived from Leconte et al. 2015
            Self::imaginary_tidal_love_number_leconte(imaginary_delta_pressure_2, star, planet)
        }
    }

    // Imaginary part of the thermal Love number defined as Leconte et al. (2015)
    fn imaginary_atmosphere_leconte(
        thermal_tide_amplitude: f64,
        radiative_frequency: f64,
        freq: f64,
        star: &impl ParticleT,
        planet: &impl ParticleT,
    ) -> f64 {
        let tmp = freq / radiative_frequency;
        // Maxwell-like frequency dependence
        let q_a = thermal_tide_amplitude * (tmp / (1.0 + tmp.powi(2)));

        Self::imaginary_tidal_love_number_leconte(q_a, star, planet)
    }

    // Imaginary part of the thermal Love number defined as Leconte et al. (2015)
    fn imaginary_tidal_love_number_leconte(
        frequency_dependence: f64,
        star: &impl ParticleT,
        planet: &impl ParticleT,
    ) -> f64 {
        -sqrt!(32.0 * PI / 15.0) * (planet.semi_major_axis().powi(3) * planet.radius())
            / (GRAVITATIONAL * star.mass() * planet.mass())
            * frequency_dependence
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Auclair-Desrotour 2017a https://doi.org/10.1051/0004-6361/201628252
// Auclair-Desrotour 2017b https://doi.org/10.1051/0004-6361/201628701
// Auclair-Desrotour et al. 2019 https://doi.org/10.1051/0004-6361/201834685
// Leconte et al. 2015 https://doi.org/10.48550/arXiv.1502.01952
