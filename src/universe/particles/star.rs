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
    Interpolated {
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
            Evolution::Interpolated { star_file_path, .. } => {
                write!(f, "Interpolated: \"{}\"", &star_file_path.display())
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
    pub(crate) footpoint_conductance: f64,           // (Ohm-1)

    // Evolution model of the star (if enabled).
    evolution: Evolution,
    // Evolving parameters
    age: f64,               // (s)
    pub(crate) radius: f64, // (m)
    convective_radius: f64, // (m)
    convective_mass: f64,   // (kg)
    pub(crate) convective_moment_of_inertia_derivative: f64,
    pub(crate) convective_moment_of_inertia: f64, // (kg.m2)
    pub(crate) radiative_moment_of_inertia: f64,  // (kg.m2)
    radiative_mass_derivative: f64,
    pub(crate) luminosity: f64, // solar units

    // Calculated internally
    dynamical_tide_dissipation: f64,
    convective_turnover_time: f64,
    convective_turnover_time_sun: f64,
    pub(crate) angular_momentum_redistribution: f64,
    pub(crate) mass_transfer_envelope_to_core_torque: f64, // structural evolution
    pub(crate) rossby: f64,
    mass_loss_rate: f64,
    magnetic_field: f64,
    pub(crate) wind_torque: f64,
    pub(crate) alfven_radius: f64,

    tidal_quality: f64,
    pub(crate) tidal_frequency: f64,
    pub(crate) magnetic_torque: f64,
    pub(crate) tidal_torque: f64,

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
        if let Evolution::Interpolated {
            ref mut interpolator,
            ..
        } = self.evolution
        {
            interpolator.init(star_ages, star_values);
        }
    }

    // Provide a reference to the stellar evolution file if evolution is interpolated.
    pub fn evolution_file(&mut self) -> Option<&PathBuf> {
        match self.evolution {
            Evolution::Interpolated {
                ref star_file_path, ..
            } => Some(star_file_path),
            Evolution::Disabled => None,
        }
    }

    // Evolves the star by interpolating the stellar evolution values by time (if evolution is enabled).
    fn stellar_evolution(&mut self, time: f64) -> Result<()> {
        if let Evolution::Interpolated {
            ref interpolator, ..
        } = self.evolution
        {
            let (age, values) = interpolator.interpolate(time)?;
            // Update the star properties with the interpolated values.
            self.age = age;
            self.radius = values[1];
            self.mass = values[2];
            self.convective_radius = values[3];
            self.convective_mass = values[4];
            self.radiative_moment_of_inertia = values[5];
            self.convective_moment_of_inertia = values[6];
            self.luminosity = values[7];
            self.radiative_mass_derivative = values[8];
            self.convective_moment_of_inertia_derivative = values[9];

            self.dynamical_tide_dissipation = self.dynamical_tide_dissipation();
        }

        Ok(())
    }

    pub(crate) fn initialise(&mut self, time: f64) -> Result<()> {
        if let Evolution::Interpolated { .. } = self.evolution {
            self.stellar_evolution(time)?;
        }
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
        let radiative_zone_mass_ratio = (self.mass - self.convective_mass) / self.mass;
        self.convective_turnover_time = Self::convective_turnover_time(radiative_zone_mass_ratio);
        self.mass_transfer_envelope_to_core_torque = self.mass_transfer_envelope_to_core_torque(); // requires convective_radius, radiative_mass_derivative, spin

        // TODO Might be possible to guard these calculations based on the tides (i.e. avoid unused calculations)
        self.rossby = self.rossby(); // requires convective_turnover_time, spin
        self.mass_loss_rate = self.mass_loss_rate(); // requires mass, rossby
        self.wind_torque = self.wind_torque(); // requires mass, radius
        self.alfven_radius = self.alfven_radius_estimate(); // requires mass_loss_rate, wind_torque

        Ok(())
    }

    // Recompute star values that depend on planet (tidal_frequency used by tidal and magnetic torque).
    pub(crate) fn refresh_tidal_frequency(&mut self, planet: &Planet) {
        self.tidal_frequency = self.tidal_frequency(planet);
    }

    // Update the tidal and magnetic torque.
    pub(crate) fn update_torques(&mut self, tidal_torque: f64, magnetic_torque: f64) {
        self.tidal_torque = tidal_torque;
        self.magnetic_torque = magnetic_torque;
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
        let beta = max!(1E-20_f64, self.convective_mass / self.mass);
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
    fn mass_transfer_envelope_to_core_torque(&self) -> f64 {
        // Takes into account the structural evolution of the star and the torques applied on both radiative and convective zones.
        (2. / 3.) * self.convective_radius.powi(2) * self.spin * self.radiative_mass_derivative
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

    // Alfven radius estimate from the stellar wind torque and mass loss rate.
    // Benbakoura et al. 2019 Eq. 7
    fn alfven_radius_estimate(&self) -> f64 {
        sqrt!(abs!(self.wind_torque) / (self.mass_loss_rate * self.spin))
    }

    // This is a re-write of Eq. 3 and 19 from Benbakoura et al. 2019
    // without the factors that are in the function semi_major_axis_13_div_2_derivative in physics.rs
    // The a^-6 is here to compensate the a^6 in physics.rs
    pub fn tidal_torque_ctl(&self, equilibrium_tide_dissipation: f64, planet: &Planet) -> f64 {
        let tidal_quality = self.tidal_quality(equilibrium_tide_dissipation);
        // Smoothing parameter when tidal frequency is 0
        let depth = 1E-08;
        -(9. / 4.)
            * planet.mass.powi(2)
            * GRAVITATIONAL
            * planet.semi_major_axis.powi(-6)
            * tanh!(self.tidal_frequency / depth)
            * self.radius.powi(5)
            / tidal_quality
    }

    // Calculates the equivalent tidal quality factor as in Mathis 2015 and Bolmont & Mathis 2016.
    pub fn tidal_quality(&self, equilibrium_tide_dissipation: f64) -> f64 {
        // Epsilon to ensure that equilibrium_tide_quality_factors stays finite
        let epsilon_secure = 1.0E-10;
        // Epsilon to ensure a smooth transition equilibrium / dynamical tide
        let epsilon_step = 1.0E-06;

        // Equilibrium tide
        // Normalization constant for equilibrium tide (Bolmont & Mathis,  2016,  Eq. 8)
        let normalisation_constant = sqrt!(GRAVITATIONAL / (SOLAR_MASS * SOLAR_RADIUS.powi(7)));
        // Tidal quality factors for the equilibrium tide
        // This is Eq. 22 of Benbakoura et al. 2019
        let equilibrium_tide_quality_factors = GRAVITATIONAL
            / (abs!(self.tidal_frequency + epsilon_secure)
                * equilibrium_tide_dissipation
                * normalisation_constant
                * self.radius.powi(5));

        // Dynamical tide
        // Equivalent Q' factor, tidal quality factor
        // This is the inverse of Eq. 4 of Bolmont & Mathis 2016,
        // which give the tidal quality factor Q' as a function of <D>w
        // (Eq. 1 of Mathis 2015, or Bolmont & Mathis 2016)
        // But here, keep in mind that the spin was removed out of <D>w to be then multiplied here
        // Q' = 3 / ( 2 * <D>w ) | Bolmont & Mathis
        // Q' = 3 / ( 2 * <D>w' * spin^2 ) | here
        let dynamical_tide_quality_factors =
            3. / (2. * self.dynamical_tide_dissipation * self.spin.powi(2));

        // Total dissipation
        // The multiplicative factor after dynamical_tide_quality_factors allows to pass continuously above 0 when inertia waves are raised
        // Benbakoura et al. 2019 Eq. 20
        let total_dissipiation = 1. / equilibrium_tide_quality_factors
            + (1. / dynamical_tide_quality_factors)
                * 0.5
                * (1. + tanh!((self.tidal_frequency + 2. * self.spin) / epsilon_step));

        1. / total_dissipiation
    }

    // Ardestani et al. 2017 Eq. A1
    fn convective_turnover_time(adjusted_convective_mass: f64) -> f64 {
        10_f64.powf(
            8.79 - 2. * abs!(log10!(adjusted_convective_mass)).powf(0.349)
                - 0.0194 * abs!(log10!(adjusted_convective_mass)).powi(2)
                - 1.62 * min!(log10!(adjusted_convective_mass) + 8.55, 0.),
        )
    }

    fn tidal_frequency(&self, planet: &Planet) -> f64 {
        // Tidal frequency for a coplanar circular orbit
        // Eq. 103 of Efroimsky 2012, for l=m=2, p=q=0
        2. * (self.spin - planet.mean_motion)
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Ardestani et al. 2017, https://doi.org/10.1093/mnras/stx2039
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Christensen-Dalsgaard et al. 1991 https://doi.org/10.1086/170441
// Efroimsky 2012, https://doi.org/10.1007/s10569-011-9397-4
// Finley et al. 2018 https://doi.org/10.3847/1538-4357/aad7b6
// MacGregor & Brenner 1991, https://doi.org/10.1086/170269
// Mathis 2015, https://doi.org/10.1051/0004-6361/201526472
// Matt et al. 2015, https://doi.org/10.1088/2041-8205/799/2/L23
// Ogilvie 2013, https://doi.org/10.1093/mnras/sts362
