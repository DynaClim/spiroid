use crate::constants::GRAVITATIONAL;
use crate::universe::particles::{ParticleT, magnetic_pressure};

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Default, Clone, PartialEq)]
#[serde(deny_unknown_fields)]
#[serde(default)]
pub struct Planet {
    // Input parameters
    pub(crate) mass: f64,            // Kg
    pub(crate) radius: f64,          // m
    pub(crate) semi_major_axis: f64, // m
    pub(crate) is_destroyed: bool,   // ()
    pub(crate) magnetic_field: f64,  // (T)

    // Input parameters, only if kaula tides enabled.
    pub(crate) spin: f64,                     // rad.s
    pub(crate) eccentricity: f64,             // ()
    pub(crate) inclination: f64,              // rad
    pub(crate) longitude_ascending_node: f64, // rad
    pub(crate) pericentre_omega: f64,         // rad
    pub(crate) spin_inclination: f64,         // rad
    pub(crate) radius_of_gyration_2: f64,     // ()

    // Calculated internally
    pub(crate) magnetic_pressure: f64,
    pub(crate) mean_motion: f64, // rad.s
    roche_limit: f64,            // m
    orbit_lower_limit: f64,      // m
    density_ratio: f64,          // ()

    // Calculated internally, only if kaula tides enabled.
    pub(crate) moment_of_inertia: f64, // kg.m^2
}

impl ParticleT for Planet {
    fn semi_major_axis(&self) -> f64 {
        self.semi_major_axis
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
        self.eccentricity
    }
    fn inclination(&self) -> f64 {
        self.inclination
    }
    fn mean_motion(&self) -> f64 {
        self.mean_motion
    }
    fn luminosity(&self) -> f64 {
        todo!()
    }
}

impl Planet {
    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Initialise the units and set initial tidal and magnetic values.
    pub(crate) fn initialise(&mut self) {
        // 1.0E-4 is a Gauss to Tesla unit conversion, to go back to S.I.
        self.magnetic_pressure = magnetic_pressure(self.magnetic_field * 1.0E-4);
        // Only relevant if kaula tides enabled.
        // Hut 1981, text just after Eq. 11
        self.moment_of_inertia = self.radius_of_gyration_2 * self.mass * self.radius.powi(2);
    }

    // Refresh internal planet values.
    // ***WARNING!***
    // Stateful function.
    // The order of these calculations is important.
    // Lower order calculations depend on previous values.
    // ***WARNING!***
    pub(crate) fn refresh(&mut self, semi_major_axis: f64, star: &impl ParticleT) {
        self.semi_major_axis = semi_major_axis;
        self.mean_motion = self.mean_motion(star.mass());
        self.density_ratio = self.density_ratio(star.mass(), star.radius());

        self.roche_limit = self.roche_limit(star.radius()); // requires density_ratio
        self.orbit_lower_limit = self.orbit_lower_limit(star.radius()); // requires density_ratio, roche_limit

        // Destroy the planet if it is too close to the star (or negative semi major axis).
        if self.crossed_orbital_lower_limit() {
            self.destroy();
        }
    }

    // Refresh additional elements, used by kaula tides.
    pub(crate) fn refresh_orbital_elements(
        &mut self,
        spin: f64,
        eccentricity: f64,
        inclination: f64,
        longitude_ascending_node: f64,
        pericentre_omega: f64,
        spin_inclination: f64,
    ) {
        // If the planet's spin is close enough to the synchronization state
        // then set the spin equal to the mean motion for numerical stability.
        if abs!(1. - abs!(spin / self.mean_motion)) < 1.0E-9 {
            self.spin = self.mean_motion;
        } else {
            self.spin = spin;
        }
        // Invert the exponent of e^2 to normalise the eccentricity.
        self.eccentricity = eccentricity;
        self.longitude_ascending_node = longitude_ascending_node;
        self.pericentre_omega = pericentre_omega;
        // If inclination < 1e-4 degrees, it is close to zero
        // so set it and the spin to zero to avoid the computation of the derivative and singularities.
        if abs!(inclination) < 1.7453e-6 {
            self.inclination = 0.0;
            self.spin_inclination = 0.0;
        } else {
            self.inclination = inclination;
            self.spin_inclination = spin_inclination;
        }
    }

    // Calculates the density_ratio (planet / star).
    // Bulk density ratio (M/R^3)
    fn density_ratio(&self, star_mass: f64, star_radius: f64) -> f64 {
        self.mass / star_mass * (star_radius / self.radius).powi(3)
    }

    // Planetary orbit's mean motion with Kepler's third law from the semi-major axis.
    // Kepler
    fn mean_motion(&self, star_mass: f64) -> f64 {
        // Note: semi_major_axis must be positive or below will be NaN
        sqrt!(GRAVITATIONAL * (star_mass + self.mass) / self.semi_major_axis.powi(3))
    }

    // Determines whether the planet is inside the alfven radius of the star.
    pub(crate) fn inside_alfven_radius(&self, alfven_radius: f64) -> bool {
        alfven_radius >= self.semi_major_axis
    }

    // Roche limit under which the planet may be tidally disrupted.
    fn roche_limit(&self, star_radius: f64) -> f64 {
        // This is from Benbakoura et al. 2019, Eq. 6 but the numerical factor
        // is not quite the same (2.42 here against 2.44) for the fluid case in Strugarek+14.
        // Antoine: This is a small change in the composition assumed for the planet.
        // We could want to do this more properly and generically.
        2.42_f64 * star_radius * self.density_ratio.powf(-1. / 3.)
    }

    // Limit of the orbit as defined in Zhang & Penev 2014.
    // Determines whether the lower limit of the orbit is the Roche limit or
    // the star radius according to the observations of Metzger et al. (2012).
    fn orbit_lower_limit(&self, star_radius: f64) -> f64 {
        if self.density_ratio > 5. {
            star_radius
        } else {
            self.roche_limit
        }
    }

    // Tests whether the planet has crossed the orbital lower limit (and will be destroyed).
    fn crossed_orbital_lower_limit(&self) -> bool {
        self.semi_major_axis.is_nan() || (self.semi_major_axis <= self.orbit_lower_limit)
    }

    // Destroy the planet. No further interaction is possible.
    fn destroy(&mut self) {
        self.is_destroyed = true;
        self.mass = 0.0;
        self.radius = 0.0;
        self.semi_major_axis = 0.0;
        self.mean_motion = 0.0;
        self.roche_limit = 0.0;
        self.mean_motion = 0.0;
        self.density_ratio = 0.0;
        self.magnetic_pressure = 0.0;
        self.roche_limit = 0.0;
        self.orbit_lower_limit = 0.0;
        self.density_ratio = 0.0;
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Hut 1981 (no DOI for that one), https://ui.adsabs.harvard.edu/abs/1981A%26A....99..126H/abstract
// Metzger et al. 2012, https://doi.org/10.1111/j.1365-2966.2012.21444.x
// Zhang & Penev 2014, https://doi.org/10.1088/0004-637X/787/2/131
