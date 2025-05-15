use crate::constants::{
    GAS_CONSTANT, GRAVITATIONAL, PI, SECONDS_IN_DAY, SECONDS_IN_YEAR, SOLAR_LUMINOSITY,
};
use crate::universe::particles::ParticleT;
use crate::utils::map_3d_to_1d;
use sci_file::Interpolator;

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use anyhow::Result;

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
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub(crate) struct LoveNumber {
    // Cache
    real: Vec<f64>,
    imaginary: Vec<f64>,

    pub(crate) imaginary_atmosphere: Interpolator<f64>,
    pub(crate) imaginary_oceanic: Interpolator<f64>,
    pub(crate) imaginary_solid: Interpolator<f64>,
    pub(crate) real_solid: Interpolator<f64>,
    love_interpolator: Vec<Interpolator<f64>>,
}

impl LoveNumber {
    /// Fetches the real love number from the cache for index of tuple (m, p, q).
    pub(crate) fn real(&self, m: usize, p: usize, q: usize) -> f64 {
        // The cached data is stored in a 1D array, so the 3D coordinates are mapped to the 1D index.
        // m_max and p_max are both 3
        self.real[map_3d_to_1d(m, 3, p, 3, q)]
    }

    /// Fetches the real love number from the cache for index of tuple (m, p, q).
    pub(crate) fn imaginary(&self, m: usize, p: usize, q: usize) -> f64 {
        // The cached data is stored in a 1D array, so the 3D coordinates are mapped to the 1D index.
        // m_max and p_max are both 3
        self.imaginary[map_3d_to_1d(m, 3, p, 3, q)]
    }

    /// Recomputes all the love number values.
    // Called at each time step to cache love numbers for that iteration, to prevent duplicate calculations.
    pub(crate) fn refresh_cache(
        &mut self,
        time: f64,
        planet: &impl ParticleT,
        star: &impl ParticleT,
        particle_type: &ParticleComposition,
    ) -> Result<()> {
        self.real.clear();
        self.imaginary.clear();
        for q in 0..=14 {
            for p in 0..3 {
                let q_fac = f64!((2 - 2 * p) + q - 7);
                for m in 0..3 {
                    let w_2lmpq = planet.mean_motion() * q_fac - (planet.spin() * f64!(m));
                    // Loop is organised so arrays are filled in order: index == 0..=135
                    self.real.push(self.real_part(w_2lmpq, particle_type)?);
                    self.imaginary.push(self.imaginary_part(
                        time,
                        w_2lmpq,
                        planet,
                        star,
                        particle_type,
                    )?);
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
    // TODO this match is called 135 times at every timestep. Since the ParticleComposition is runtime static, this might be inefficient.
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
            } => Ok(imaginary_atmosphere(thermal_tide_model, freq, planet, star)),
            ParticleComposition::SolidAtmosphere {
                thermal_tide_model, ..
            } => Ok(self.imaginary_solid(freq)?
                + imaginary_atmosphere(thermal_tide_model, freq, planet, star)),
            ParticleComposition::SolidOcean { .. } => {
                Ok(self.imaginary_solid(freq)? + self.imaginary_oceanic(freq)?)
            }
            ParticleComposition::SolidAtmosphereOcean {
                thermal_tide_model, ..
            } => Ok(self.imaginary_solid(freq)?
                + imaginary_atmosphere(thermal_tide_model, freq, planet, star)
                + self.imaginary_oceanic(freq)?),
            ParticleComposition::Interpolate { .. } => {
                self.love_number_interpolated_by_frequency(time, freq)
            }
            ParticleComposition::InterpolateAtmosphere {
                thermal_tide_model, ..
            } => Ok(self.love_number_interpolated_by_frequency(time, freq)?
                + imaginary_atmosphere(thermal_tide_model, freq, planet, star)),
        }
    }

    fn imaginary_solid(&self, freq: f64) -> Result<f64> {
        if let 0.0 = freq {
            Ok(0.0)
        } else {
            let (_, im_k2) = self.imaginary_solid.interpolate(abs!(freq))?;
            Ok(freq.signum() * im_k2)
        }
    }

    fn imaginary_oceanic(&self, freq: f64) -> Result<f64> {
        if let 0.0 = freq {
            Ok(0.0)
        } else {
            let (_, im_k2) = self.imaginary_oceanic.interpolate(abs!(freq))?;
            Ok(im_k2)
        }
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
        let index = (time / 1.0E9 / SECONDS_IN_YEAR) as usize * 10 - 10;
        let love_interpolator = &self.love_interpolator[index];
        if let 0.0 = freq {
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
    Auclair,
    AuclairScaling,
    Leconte,
}

pub(crate) fn imaginary_atmosphere(
    thermal_tide_model: &ThermalTideModel,
    freq: f64,
    planet: &impl ParticleT,
    star: &impl ParticleT,
) -> f64 {
    match thermal_tide_model {
        ThermalTideModel::Analytic => imaginary_atmosphere_analytic(freq, star, planet),
        ThermalTideModel::Auclair => imaginary_atmosphere_auclair(freq, star, planet),
        ThermalTideModel::AuclairScaling => {
            imaginary_atmosphere_auclair_scaling(freq, star, planet)
        }
        ThermalTideModel::Leconte => imaginary_atmosphere_leconte(freq, star, planet),
    }
}

// Auclair-desrotour 2017a Eq. 173
// Values of physical parameter table 1
fn imaginary_atmosphere_analytic(freq: f64, star: &impl ParticleT, planet: &impl ParticleT) -> f64 {
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
fn imaginary_atmosphere_auclair(freq: f64, star: &impl ParticleT, planet: &impl ParticleT) -> f64 {
    // Related to the first adiabatic exponent of the gas.
    let kappa = 0.286;
    let surface_temperature = 737.;
    // Radiative thermal frequency of the atmosphere (s^-1).
    let omega = 2_f64 * 3.77e-7;
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
        // Scaled thermal time scale and amplitude (using scaling formulation with fixed sma a = a_venus) Eq. 44 and 45
        // Surface pressure in bar
        let p_s = 0.27;
        // Scaled pressure Eq. 44
        let q_0 = 10.0_f64.powf(0.48 * log10!(p_s) + 2.87);
        // Scaled time-scale Eq. 45
        let tau_0 = 10.0_f64.powf(0.3 * log10!(p_s) + 0.038);

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
        imaginary_tidal_love_number_leconte(imaginary_delta_pressure_2, star, planet)
    }
}

// Imaginary part of the thermal Love number defined as Leconte et al. (2015)
fn imaginary_atmosphere_leconte(freq: f64, star: &impl ParticleT, planet: &impl ParticleT) -> f64 {
    // Amplitude of the atmospheric quadrupole (Pa)
    let q_0 = 201.0;
    // Radiative thermal frequency of the atmosphere (s^-1)
    let omega = 2_f64 * 3.77e-7;
    let tmp = freq / omega;
    // Maxwell-like frequency dependence
    let q_a = q_0 * (tmp / (1.0 + tmp.powi(2)));

    imaginary_tidal_love_number_leconte(q_a, star, planet)
}

// Imaginary part of the thermal Love number defined as Leconte et al. (2015)
fn imaginary_tidal_love_number_leconte(
    frequency_dependence: f64,
    star: &impl ParticleT,
    planet: &impl ParticleT,
) -> f64 {
    -sqrt!(32.0_f64 * PI / 15.0) * (planet.semi_major_axis().powi(3) * planet.radius())
        / (GRAVITATIONAL * star.mass() * planet.mass())
        * frequency_dependence
}

#[cfg(test)]
pub mod tests;

// References:
// Auclair-Desrotour 2017a https://doi.org/10.1051/0004-6361/201628252
// Auclair-Desrotour 2017b https://doi.org/10.1051/0004-6361/201628701
// Auclair-Desrotour et al. 2019 https://doi.org/10.1051/0004-6361/201834685
// Leconte et al. 2015 https://doi.org/10.48550/arXiv.1502.01952
