pub(crate) mod star_csv;
use crate::constants::{
    GRAVITATIONAL, PI, ROSSBY_SATURATION, ROSSBY_SUN, SECONDS_IN_YEAR, SOLAR_ANGULAR_VELOCITY,
    SOLAR_MASS, SOLAR_MASS_LOSS_RATE, SOLAR_RADIUS, TWO_PI,
};
use crate::universe::particles::{ParticleT, Planet};
use serde::{Deserialize, Serialize};
pub use star_csv::StarCsv;
use std::path::PathBuf;

use anyhow::Result;
use sci_file::Interpolator;

#[derive(Deserialize, Serialize, PartialEq, Clone, Default)]
enum Evolution {
    #[default]
    Disabled,
    Starevol {
        star_file_path: PathBuf,
        #[serde(skip)]
        interpolator: Interpolator<Vec<f64>>,
    },
    Mesa {
        star_file_path: PathBuf,
        #[serde(skip)]
        interpolator: Interpolator<Vec<f64>>,
    },
}

// Custom debug implementation to only print the stellar evolution file name instead of a data dump.
impl std::fmt::Debug for Evolution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Evolution::Disabled => write!(f, "Disabled"),
            Evolution::Starevol { star_file_path, .. } => {
                write!(f, "Starevol: \"{}\"", &star_file_path.display())
            }
            Evolution::Mesa { star_file_path, .. } => {
                write!(f, "Mesa: \"{}\"", &star_file_path.display())
            }
        }
    }
}

#[derive(Debug, Deserialize, Serialize, PartialEq, Default, Clone)]
#[serde(deny_unknown_fields)]
#[serde(default)]
pub struct Star {
    // Input parameters
    pub(crate) mass: f64,                            // (kg)
    pub(crate) spin: f64,                            // (rad.s-1)
    pub(crate) core_envelope_coupling_constant: f64, // (s)

    // Evolution model of the star (if enabled).
    evolution: Evolution,
    // Evolving parameters
    age: f64,                       // (s)
    pub(crate) radius: f64,         // (m)
    convective_radius: f64,         // (m)
    pub(crate) radiative_mass: f64, // (kg)
    pub(crate) convective_moment_of_inertia_derivative: f64,
    pub(crate) convective_moment_of_inertia: f64, // (kg.m2)
    pub(crate) radiative_moment_of_inertia: f64,  // (kg.m2)
    radiative_mass_derivative: f64,
    pub(crate) luminosity: f64, // solar units

    // Calculated internally
    pub(crate) dynamical_tide_dissipation: f64,
    pub(crate) convective_turnover_time: f64,
    convective_turnover_time_sun: f64,
    pub(crate) angular_momentum_redistribution: f64,
    pub(crate) mass_transfer_envelope_to_core_torque: f64, // structural evolution
    pub(crate) rossby: f64,
    mass_loss_rate: f64, // (kg.s-1)
    magnetic_field: f64,
    pub(crate) wind_torque: f64,
    pub(crate) alfven_radius: f64,

    tidal_quality: f64,
    pub(crate) tidal_frequency: f64,
    pub(crate) magnetic_torque: f64,
    pub(crate) tidal_torque_convective: f64,

    // Evolved parameters
    pub(crate) evolved_wind_torque: f64,
    pub(crate) evolved_change_semi_major_axis: f64,
    // Additional mass loss rate during the evolved phase of the star.
    evolved_mass_loss_rate: f64, // (kg.s-1)
    terminal_wind_speed: f64,    // (m.s-1)
    mass_accretion_efficiency: f64,
    wind_orbital_angular_momentum_loss: f64,

    // Integration parameters
    pub(crate) convective_zone_angular_momentum: f64, // (kg.m^2.s-1)
    pub(crate) radiative_zone_angular_momentum: f64,  // (kg.m^2.s-1)
}

impl ParticleT for Star {
    fn semi_major_axis(&self) -> f64 {
        todo!();
    }
    fn mean_motion(&self) -> f64 {
        todo!();
    }
    fn mass(&self) -> f64 {
        self.mass
    }
    fn radius(&self) -> f64 {
        self.radius
    }
    fn spin(&self) -> f64 {
        self.spin
    }
    fn eccentricity(&self) -> f64 {
        todo!()
    }
    fn inclination(&self) -> f64 {
        todo!()
    }
    fn luminosity(&self) -> f64 {
        self.luminosity
    }
}

impl Star {
    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        Self::default()
    }

    // Returns `true` if evolution is enabled for the star `Evolution`.
    pub fn evolution_enabled(&self) -> bool {
        !matches!(self.evolution, Evolution::Disabled)
    }

    // Initialise stellar values from the stellar evolution file if evolution is interpolated.
    pub fn initialise_evolution(&mut self, star_ages: &[f64], star_values: &[Vec<f64>]) {
        match self.evolution {
            Evolution::Disabled => {}
            Evolution::Starevol {
                ref mut interpolator,
                ..
            }
            | Evolution::Mesa {
                ref mut interpolator,
                ..
            } => {
                interpolator.init(star_ages, star_values);
            }
        }
    }

    // Provide a reference to the stellar evolution file if evolution is interpolated.
    pub fn evolution_file(&mut self) -> Option<&PathBuf> {
        match self.evolution {
            Evolution::Starevol {
                ref star_file_path, ..
            }
            | Evolution::Mesa {
                ref star_file_path, ..
            } => Some(star_file_path),
            Evolution::Disabled => None,
        }
    }

    // Evolves the star by interpolating the stellar evolution values by time (if evolution is enabled).
    fn stellar_evolution(&mut self, time: f64) -> Result<()> {
        match self.evolution {
            Evolution::Disabled => Ok(()),
            Evolution::Mesa {
                ref interpolator, ..
            }
            | Evolution::Starevol {
                ref interpolator, ..
            } => {
                let (age, values) = interpolator.interpolate(time)?;
                // Update the star properties with the interpolated values.
                self.age = age;
                self.radius = values[1];
                self.mass = values[2];
                self.convective_radius = values[3];
                self.radiative_mass = values[4];
                self.radiative_moment_of_inertia = values[5];
                self.convective_moment_of_inertia = values[6];
                self.luminosity = values[7];
                self.radiative_mass_derivative = values[8];
                self.convective_moment_of_inertia_derivative = values[9];

                if matches!(self.evolution, Evolution::Mesa { .. }) {
                    self.convective_turnover_time = values[10];
                    self.evolved_mass_loss_rate = values[11];
                }

                self.dynamical_tide_dissipation = self.dynamical_tide_dissipation();

                Ok(())
            }
        }
    }

    pub(crate) fn initialise(&mut self, time: f64) -> Result<()> {
        self.stellar_evolution(time)?;
        // 0.02 is the convection zone mass of the Sun divided by its total mass.
        // Christensen-Dalsgaard et al. 1991
        self.convective_turnover_time_sun = Self::convective_turnover_time(0.02);

        Ok(())
    }

    // Recompute independent star values and interpolate if required.
    // ***WARNING!***
    // Stateful function.
    // The order of these calculations is important.
    // Lower order calculations depend on previous values.
    // ***WARNING!***
    pub(crate) fn refresh(
        &mut self,
        time: f64,
        radiative_zone_angular_momentum: f64,
        convective_zone_angular_momentum: f64,
        disk_is_dissipated: bool,
    ) -> Result<()> {
        self.radiative_zone_angular_momentum = radiative_zone_angular_momentum;
        self.convective_zone_angular_momentum = convective_zone_angular_momentum;

        self.stellar_evolution(time)?;
        // Update the spin only after the disk has dissipated.
        if disk_is_dissipated {
            self.spin = self.spin(); // requires convective_zone_angular_momentum, convective_moment_of_inertia
        }

        self.angular_momentum_redistribution = self.angular_momentum_redistribution(); // requires convective_moment_of_inertia, radiative_moment_of_inertia, convective_zone_angular_momentum, radiative_zone_angular_momentum
        self.mass_transfer_envelope_to_core_torque = self.mass_transfer_envelope_to_core_torque(); // requires convective_radius, radiative_mass_derivative, spin

        if matches!(self.evolution, Evolution::Mesa { .. }) {
            self.core_envelope_coupling_constant = self.evolving_core_envelope_coupling_constant(); // requres mass, spin
        } else {
            // Only used by tides and magnetism
            let convective_zone_mass_ratio = (self.mass - self.radiative_mass) / self.mass;
            self.convective_turnover_time =
                Self::convective_turnover_time(convective_zone_mass_ratio);
        }

        self.rossby = self.rossby(); // requires convective_turnover_time, spin
        self.mass_loss_rate = self.mass_loss_rate(); // requires mass, rossby

        // Zero the torques. They will be calculated if associated effects are enabled.
        self.tidal_torque_convective = 0.0;
        self.magnetic_torque = 0.0;
        self.wind_torque = 0.0;
        self.alfven_radius = 0.0;

        Ok(())
    }

    // Recompute star values that depend on planet (tidal_frequency used by tidal and magnetic torque).
    pub(crate) fn refresh_tidal_frequency(&mut self, planet: &Planet) {
        self.tidal_frequency = self.tidal_frequency(planet);
    }

    // Update the tidal torque.
    pub(crate) fn update_tidal_torque(&mut self, tidal_torque_convective: f64) {
        self.tidal_torque_convective = tidal_torque_convective;
    }

    // Update the magnetic torque.
    pub(crate) fn update_magnetic_torque(&mut self, magnetic_torque: f64) {
        self.magnetic_torque = magnetic_torque;
    }

    // Update the wind torque.
    pub(crate) fn update_wind_torque(&mut self, enabled: bool) {
        if enabled {
            self.wind_torque = self.wind_torque(); // requires mass, radius
            // alfven_radius is recalculated with the updated wind_torque
            self.alfven_radius = self.alfven_radius_estimate(); // requires mass_loss_rate, wind_torque
            self.evolved_wind_torque = self.evolved_wind_torque(); // requires spin, evolved_mass_loss_rate, radius
        }
    }

    pub(crate) fn update_evolved_change_semi_major_axis(&mut self, enabled: bool, planet: &Planet) {
        if enabled {
            self.evolved_change_semi_major_axis = self.evolved_change_semi_major_axis(planet); // requires mass, radius, planet mass, semi_major_axis, mean_motion
        }
    }

    // TODO Only works for circular orbits, can be extended for eccentric orbits.
    // The change in semi-major axis is computed using the terminal wind speed and orbital velocity.
    // Depending on these two velocities, the mass accretion efficiency and wind orbital angular momentum loss are computed.
    // These quantities are then used to compute the change in semi-major axis.
    fn evolved_change_semi_major_axis(&mut self, planet: &Planet) -> f64 {
        self.terminal_wind_speed = self.terminal_wind_speed(); // requires mass, radius
        let orbital_velocity = planet.mean_motion * planet.semi_major_axis; // requires semi_major_axis, mean_motion
        let mass_ratio = planet.mass / self.mass; // requires mass, planet mass; Esseldeurs et al. 2025, below Eq. 18
        self.mass_accretion_efficiency =
            self.mass_accretion_efficiency(mass_ratio, orbital_velocity); // requires terminal_wind_speed
        self.wind_orbital_angular_momentum_loss =
            self.wind_orbital_angular_momentum_loss(mass_ratio, orbital_velocity); // requires terminal_wind_speed

        // Esseldeurs et al. 2025, Eq. 20
        2. * planet.semi_major_axis * self.evolved_mass_loss_rate / self.mass
            * (1.
                - self.mass_accretion_efficiency / mass_ratio
                - self.wind_orbital_angular_momentum_loss
                    * (1. - self.mass_accretion_efficiency)
                    * (1. + mass_ratio)
                    / mass_ratio
                - (1. - self.mass_accretion_efficiency) / 2. / (1. + mass_ratio))
    }

    // Computes the terminal wind speed.
    // This is the speed at which the stellar wind reaches its maximum velocity.
    // Esseldeurs et al. 2025, Eq. 22
    fn terminal_wind_speed(&self) -> f64 {
        let alpha_wind = 1. / 8.; // Esseldeurs et al. 2025, below Eq. 22
        sqrt!(2. * alpha_wind * GRAVITATIONAL * self.mass / self.radius)
    }

    // Computes the mass accretion efficiency.
    // This is the fraction of the stellar wind that is accreted by the planet.
    // Esseldeurs et al. 2025, Eq. 20, based on Saladino et al. 2019
    fn mass_accretion_efficiency(&self, mass_ratio: f64, orbital_velocity: f64) -> f64 {
        // Bondi-Hoyle-Lyttleton accretion, see Edgar 2004
        let mass_accretion_efficiency_bhl = mass_ratio.powi(2) / (1. + mass_ratio).powi(2)
            * orbital_velocity.powi(4)
            / (self.terminal_wind_speed
                * (self.terminal_wind_speed.powi(2) + orbital_velocity.powi(2)).powf(1.5));
        // Esseldeurs et al. 2025, below Eq. 20
        let mass_accretion_efficiency = (0.75
            + 1.0
                / (1.7
                    + 0.3 / mass_ratio
                    + ((0.5 + 0.2 / mass_ratio) * self.terminal_wind_speed / orbital_velocity)
                        .powi(5)))
            * mass_accretion_efficiency_bhl;
        min!(mass_accretion_efficiency, 0.3_f64, 1.4 * mass_ratio.powi(2))
    }

    // Computes the wind orbital angular momentum loss.
    // This is the fraction of the orbital angular momentum that is lost due to the stellar wind.
    // Esseldeurs et al. 2025, Eq. 21, based on Saladino et al. 2019
    fn wind_orbital_angular_momentum_loss(&self, mass_ratio: f64, orbital_velocity: f64) -> f64 {
        let wind_orbital_angular_momentum_loss_iso = mass_ratio.powi(2) / (1. + mass_ratio).powi(2);
        let wind_orbital_angular_momentum_loss = 1.0
            / (max!(mass_ratio.powi(-1), 0.6 * mass_ratio.powf(-1.7))
                + ((1.5 + 0.3 / mass_ratio) * self.terminal_wind_speed / orbital_velocity).powi(3))
            + wind_orbital_angular_momentum_loss_iso;
        min!(wind_orbital_angular_momentum_loss, 0.6)
    }

    fn dynamical_tide_dissipation(&self) -> f64 {
        // Computing the dynamical_tide_dissipation (cf. Ogilvie 2013,  Mathis 2015)
        // Critical angular velocity of the star (for instance, page 2 of Mathis 2015)
        let omega_crit = sqrt!(GRAVITATIONAL * self.mass / self.radius.powi(3));
        // Radius aspect ratio
        let alpha = self.convective_radius / self.radius;
        // Mass aspect ratio
        // Both mass and radius aspect ratio can be zero before the convective core appears on the PMS
        // But it's only a problem if beta is zero (gamma has a 1/beta), so the 1e-20 is there to prevent NaNs
        let beta = max!(1E-20_f64, self.radiative_mass / self.mass);
        // Gamma parameter from Mathis 2015, Eq.2
        let gamma = max!(
            1E-20_f64,
            alpha.powi(3) * (1. - beta) / (beta * (1. - alpha.powi(3)))
        );
        // Frequency-averaged tidal dissipation (Eq. B3 of Ogilvie 2013, or Eq. 1 of Mathis 2015)
        // but without the Spin^2, which was taken out here and multiplied back in fn tidal_quality
        let dynamical_tide_dissipation = omega_crit.powi(-2) * 100. * PI / 63.
            * (alpha.powi(5) / (1. - alpha.powi(5)))
            * (1. - gamma).powi(2)
            * (1. - alpha).powi(4)
            * (1. + 2. * alpha + 3. * alpha.powi(2) + 1.5 * alpha.powi(3)).powi(2)
            * (1. + ((1. - gamma) / gamma) * alpha.powi(3))
            * (1.
                + 1.5 * gamma
                + 2.5 / gamma * (1. + 0.5 * gamma - 1.5 * gamma.powi(2)) * alpha.powi(3)
                - 2.25 * (1. - gamma) * alpha.powi(5))
            .powi(-2);

        // 1e-20 is to prevent the dissipation to go to zero.
        // Maybe in case spin = 0.
        max!(1E-20_f64, dynamical_tide_dissipation)
    }

    // Angular momentum redistribution. See MacGregor & Brenner 1991,  Eq. 1
    fn angular_momentum_redistribution(&self) -> f64 {
        (self.convective_moment_of_inertia * self.radiative_zone_angular_momentum
            - self.radiative_moment_of_inertia * self.convective_zone_angular_momentum)
            / (self.convective_moment_of_inertia + self.radiative_moment_of_inertia)
    }

    // Adjust the mass in each layer (radiative and convective) based on the stellar evolution model.
    // Benbakoura et al. 2019, Eq 2.
    // TODO does this need another reference due to changes?
    fn mass_transfer_envelope_to_core_torque(&self) -> f64 {
        // Takes into account the structural evolution of the star and the torques applied on both radiative and convective zones.
        if self.radiative_mass_derivative >= 0.0 {
            // If the radiative mass is increasing, the rotation of the convective zone is transferred to the radiative zone.
            (2. / 3.) * self.convective_radius.powi(2) * self.spin * self.radiative_mass_derivative
        } else if self.radiative_zone_angular_momentum == 0.0
            || self.radiative_moment_of_inertia == 0.0
        {
            0.0
        } else {
            let minspin = min!(
                self.spin,
                self.radiative_zone_angular_momentum / self.radiative_moment_of_inertia
            );
            (2. / 3.) * self.convective_radius.powi(2) * minspin * self.radiative_mass_derivative
        }
    }

    // Computes the Rossby number.
    // Ardestani et al. 2017
    fn rossby(&self) -> f64 {
        (TWO_PI / self.spin) / self.convective_turnover_time
    }

    // Calculate spin from the ratio of angular momentum and moment inertia.
    fn spin(&self) -> f64 {
        self.convective_zone_angular_momentum / self.convective_moment_of_inertia
    }

    // Mass loss rate in the stellar wind.
    // Matt et al. 2015, Eq. 4
    fn mass_loss_rate(&self) -> f64 {
        // Mass loss rate due to stellar wind
        let mass_loss = SOLAR_MASS_LOSS_RATE
            * (max!(self.rossby, ROSSBY_SATURATION) / ROSSBY_SUN).powi(-2)
            * (self.mass / SOLAR_MASS).powi(4);

        mass_loss * SOLAR_MASS / SECONDS_IN_YEAR
    }

    // Stellar wind torque.
    // Matt et al. 2015, Eq. 3
    fn wind_torque(&self) -> f64 {
        // Torque applied on the envelope by the wind
        // Solar wind torque, in Joule (Matt et al. 2015)
        // There is a debate in the community about the value of solar_wind_torque_sun.
        // The best estimate so far is from Finley et al. (2018), giving 2.9e30 erg = 2.9e23 J
        // Most scaling laws were adjusted with this constant as 8e23 to recover the Sun.
        // A clean study should be made again before changing this.

        let gamma =
            8e23 * (self.radius / SOLAR_RADIUS).powf(3.1) * (self.mass / SOLAR_MASS).powf(0.5);
        // Wind braking torque in Joules, following (Matt et al. 2015)
        if self.rossby > ROSSBY_SATURATION {
            -gamma
                * (self.convective_turnover_time / self.convective_turnover_time_sun).powi(2)
                * (self.spin / SOLAR_ANGULAR_VELOCITY).powi(3)
        } else {
            -gamma * (ROSSBY_SUN / ROSSBY_SATURATION).powi(2) * (self.spin / SOLAR_ANGULAR_VELOCITY)
        }
    }

    // Stellar wind torque during the evolved phases of the star.
    // Dust wind torque
    // Madappatt et al 2016, Eq. 2
    fn evolved_wind_torque(&self) -> f64 {
        -2. / 3. * self.spin * self.evolved_mass_loss_rate * self.radius.powi(2)
    }

    // Gallet & Delorme 2019, Eq. 1.
    // Evolving core-envelope coupling constant that reduces
    // the coupling efficiency during the evolved phases.
    fn evolving_core_envelope_coupling_constant(&self) -> f64 {
        74.6e6
            * SECONDS_IN_YEAR
            * (self.mass / SOLAR_MASS).powf(-3.83)
            * (abs!(self.spin) / SOLAR_ANGULAR_VELOCITY).powf(-0.69)
    }

    // Alfven radius estimate from the stellar wind torque and mass loss rate.
    // Benbakoura et al. 2019, Eq. 7
    fn alfven_radius_estimate(&self) -> f64 {
        sqrt!(abs!(self.wind_torque) / (self.mass_loss_rate * abs!(self.spin)))
    }

    // Ardestani et al. 2017 Eq. A1
    fn convective_turnover_time(adjusted_convective_mass: f64) -> f64 {
        10_f64.powf(
            8.79 - 2. * abs!(log10!(adjusted_convective_mass)).powf(0.349)
                - 0.0194 * abs!(log10!(adjusted_convective_mass)).powi(2)
                - 1.62 * min!(log10!(adjusted_convective_mass) + 8.55, 0.),
        )
    }

    // Tidal frequency for a coplanar circular orbit
    // Efroimsky 2012, Eq. 103 (for l = m = 2, p = q = 0)
    fn tidal_frequency(&self, planet: &Planet) -> f64 {
        2. * (self.spin - planet.mean_motion)
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Ardestani et al. 2017, https://doi.org/10.1093/mnras/stx2039
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Christensen-Dalsgaard et al. 1991 https://doi.org/10.1086/170441
// Edgar 2004, https://doi.org/10.1016/j.newar.2004.06.001
// Efroimsky 2012, https://doi.org/10.1007/s10569-011-9397-4
// Esseldeurs et al. 2025, submitted TODO update when published
// Finley et al. 2018 https://doi.org/10.3847/1538-4357/aad7b6
// Gallet & Delorme 2019, https://doi.org/10.1051/0004-6361/201834898
// MacGregor & Brenner 1991, https://doi.org/10.1086/170269
// Madappatt et al, 2016, https://doi.org/10.1093/mnras/stw2025
// Mathis 2015, https://doi.org/10.1051/0004-6361/201526472
// Matt et al. 2015, https://doi.org/10.1088/2041-8205/799/2/L23
// Ogilvie 2013, https://doi.org/10.1093/mnras/sts362
// Saladino et al. 2019, https://doi.org/10.1051/0004-6361/201834598
