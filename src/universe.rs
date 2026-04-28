pub(crate) mod effects;
pub(crate) mod particles;
pub(crate) mod physics;

pub use effects::Kaula;
pub use particles::{Particle, ParticleType, Planet, Star, StarCsv};

use anyhow::Result;
use derive_more::{Add, Div, Mul, Sub};
use serde::{Deserialize, Serialize};
use simulation::DopriNumOps;

// Returns true if number is denormal
fn denormal_check(num: f64) -> bool {
    !(num == 0.0 || num.is_normal())
}
#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub enum UniverseKind {
    StarPlanet,
    BinaryStar, //TODO unimplemented
    PlanetMoon, //TODO unimplemented
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct Universe {
    #[serde(default)]
    time: f64,
    disk_lifetime: f64,
    #[serde(default)]
    pub disk_is_dissipated: bool,
    #[serde(default)]
    pub derivatives: UniverseIntegral,
    /// The central body of the simulation, e.g. the star.
    pub central_body: Particle,
    /// The oribiting body of the simulation, e.g. the planet.
    pub orbiting_body: Particle,
}

#[derive(
    Deserialize, Serialize, Debug, Clone, Copy, PartialEq, Default, Add, Mul, Div, Sub, DopriNumOps,
)]
#[mul(forward)]
#[div(forward)]
pub struct UniverseIntegral {
    central_body: StarIntegral,
    orbiting_body: PlanetIntegral,
}

// Manually implementing Mul<f64>
// Deriving Mul not possible for both Mul<Self> and Mul<T>
// https://github.com/JelteF/derive_more/issues/361
// https://github.com/JelteF/derive_more/pull/450
impl<T> std::ops::Mul<T> for UniverseIntegral
where
    T: Clone + Copy,
    f64: std::ops::Mul<T, Output = f64>,
{
    type Output = UniverseIntegral;
    fn mul(self, scalar: T) -> UniverseIntegral {
        UniverseIntegral {
            central_body: self.central_body * scalar,
            orbiting_body: self.orbiting_body * scalar,
        }
    }
}
impl<T> std::ops::Add<T> for UniverseIntegral
where
    T: Clone + Copy,
    f64: std::ops::Add<T, Output = f64>,
{
    type Output = UniverseIntegral;
    fn add(self, scalar: T) -> UniverseIntegral {
        UniverseIntegral {
            central_body: self.central_body + scalar,
            orbiting_body: self.orbiting_body + scalar,
        }
    }
}

impl UniverseIntegral {
    fn zero(&mut self) {
        self.central_body.zero();
        self.orbiting_body.zero();
    }
    fn denormal_check(&self) -> bool {
        self.central_body.denormal_check() || self.orbiting_body.denormal_check()
    }
}

#[derive(
    Deserialize, Serialize, Debug, Clone, Copy, PartialEq, Default, Add, Mul, Div, Sub, DopriNumOps,
)]
#[mul(forward)]
#[div(forward)]
struct StarIntegral {
    radiative_zone_angular_momentum: f64,
    convective_zone_angular_momentum: f64,
}

impl StarIntegral {
    // Creates initial quantities to be integrated for the particle.
    fn new(star: &Star) -> StarIntegral {
        StarIntegral {
            radiative_zone_angular_momentum: star.spin * star.radiative_moment_of_inertia,
            convective_zone_angular_momentum: star.spin * star.convective_moment_of_inertia,
        }
    }

    fn zero(&mut self) {
        self.radiative_zone_angular_momentum = 0.0;
        self.convective_zone_angular_momentum = 0.0;
    }
    fn denormal_check(&self) -> bool {
        denormal_check(self.radiative_zone_angular_momentum)
            || denormal_check(self.convective_zone_angular_momentum)
    }
}
impl<T> std::ops::Mul<T> for StarIntegral
where
    T: Clone + Copy,
    f64: std::ops::Mul<T, Output = f64>,
{
    type Output = StarIntegral;
    fn mul(self, scalar: T) -> StarIntegral {
        StarIntegral {
            radiative_zone_angular_momentum: self.radiative_zone_angular_momentum * scalar,
            convective_zone_angular_momentum: self.convective_zone_angular_momentum * scalar,
        }
    }
}
impl<T> std::ops::Add<T> for StarIntegral
where
    T: Clone + Copy,
    f64: std::ops::Add<T, Output = f64>,
{
    type Output = StarIntegral;
    fn add(self, scalar: T) -> StarIntegral {
        StarIntegral {
            radiative_zone_angular_momentum: self.radiative_zone_angular_momentum + scalar,
            convective_zone_angular_momentum: self.convective_zone_angular_momentum + scalar,
        }
    }
}

// Only if kaula tides are enabled on the planet:
#[derive(
    Deserialize, Serialize, Debug, Clone, Copy, PartialEq, Default, Add, Mul, Div, Sub, DopriNumOps,
)]
#[mul(forward)]
#[div(forward)]
struct PlanetIntegral {
    // semi-major axis^6.5
    semi_major_axis: f64,
    // spin
    spin: f64,
    // orbital eccentricity^2
    eccentricity: f64,
    // orbital inclination (with respect to the planet equatorial plane)
    inclination: f64,
    // longitude of ascending node
    longitude_ascending_node: f64,
    // argument of periapsis
    pericentre_omega: f64,
    // spin axis inclination (with respect to the total angular momentum)
    spin_inclination: f64,
}
impl PlanetIntegral {
    // Creates initial quantities to be integrated for the particle.
    fn new(planet: &Planet) -> PlanetIntegral {
        PlanetIntegral {
            // The actual integrated quantity is sma^6.5.
            // See comment for planet_semi_major_axis_13_div_2_derivative
            semi_major_axis: planet.semi_major_axis.powf(6.5),
            spin: planet.spin,
            // The actual integrated quantity is eccentricity^2 to avoid singularities when eccentricity goes to 0.
            eccentricity: planet.eccentricity.powi(2),
            inclination: planet.inclination,
            longitude_ascending_node: planet.longitude_ascending_node,
            pericentre_omega: planet.pericentre_omega,
            spin_inclination: planet.spin_inclination,
        }
    }
    fn zero(&mut self) {
        self.semi_major_axis = 0.0;
        self.spin = 0.0;
        self.eccentricity = 0.0;
        self.inclination = 0.0;
        self.longitude_ascending_node = 0.0;
        self.pericentre_omega = 0.0;
        self.spin_inclination = 0.0;
    }
    fn denormal_check(&self) -> bool {
        denormal_check(self.semi_major_axis)
            || denormal_check(self.spin)
            || denormal_check(self.eccentricity)
            || denormal_check(self.inclination)
            || denormal_check(self.longitude_ascending_node)
            || denormal_check(self.pericentre_omega)
            || denormal_check(self.spin_inclination)
    }
}
impl<T> std::ops::Mul<T> for PlanetIntegral
where
    T: Clone + Copy,
    f64: std::ops::Mul<T, Output = f64>,
{
    type Output = PlanetIntegral;
    fn mul(self, scalar: T) -> PlanetIntegral {
        PlanetIntegral {
            semi_major_axis: self.semi_major_axis * scalar,
            spin: self.spin * scalar,
            eccentricity: self.eccentricity * scalar,
            inclination: self.inclination * scalar,
            longitude_ascending_node: self.longitude_ascending_node * scalar,
            pericentre_omega: self.pericentre_omega * scalar,
            spin_inclination: self.spin_inclination * scalar,
        }
    }
}
impl<T> std::ops::Add<T> for PlanetIntegral
where
    T: Clone + Copy,
    f64: std::ops::Add<T, Output = f64>,
{
    type Output = PlanetIntegral;
    fn add(self, scalar: T) -> PlanetIntegral {
        PlanetIntegral {
            semi_major_axis: self.semi_major_axis + scalar,
            spin: self.spin + scalar,
            eccentricity: self.eccentricity + scalar,
            inclination: self.inclination + scalar,
            longitude_ascending_node: self.longitude_ascending_node + scalar,
            pericentre_omega: self.pericentre_omega + scalar,
            spin_inclination: self.spin_inclination + scalar,
        }
    }
}

impl Universe {
    /// Initialise the `Star` and `Planet`.
    /// # Errors
    ///
    /// Will return `Err` if evolution is enabled on the `Star`(s)
    /// and time is outside of the interpolation range.
    pub fn initialise(&mut self, time: f64) -> Result<()> {
        self.central_body.initialise(time)?;
        self.orbiting_body.initialise(time)?;

        Ok(())
    }

    // Determines the kind of universe, based on the particles.
    pub fn universe_kind(&self) -> UniverseKind {
        if self.central_body.is_star() && self.orbiting_body.is_planet() {
            UniverseKind::StarPlanet
        } else if self.central_body.is_star() && self.orbiting_body.is_star() {
            UniverseKind::BinaryStar
        } else if self.central_body.is_planet() && self.orbiting_body.is_planet() {
            UniverseKind::PlanetMoon
        } else {
            unreachable!()
        }
    }

    // Creates a vector of initial quantities to be integrated, depending on the simulation configuration.
    pub fn integration_quantities(&self) -> UniverseIntegral {
        match self.universe_kind() {
            UniverseKind::StarPlanet => {
                let ParticleType::Star(star) = &self.central_body.kind else {
                    unreachable!()
                };
                let ParticleType::Planet(planet) = &self.orbiting_body.kind else {
                    unreachable!()
                };

                UniverseIntegral {
                    central_body: StarIntegral::new(star),
                    orbiting_body: PlanetIntegral::new(planet),
                }
            }
            UniverseKind::BinaryStar => todo!(),
            UniverseKind::PlanetMoon => todo!(),
        }
    }

    // Provide values to bound the subsequent step size
    pub(crate) fn interpolation_step_size_hint(&self, time: f64) -> Option<f64> {
        let ParticleType::Star(star) = &self.central_body.kind else {
            unreachable!()
        };
        star.stellar_evolution_step_size_hint(time)
    }

    fn disk_is_dissipated(&mut self) {
        self.disk_is_dissipated = self.disk_is_dissipated || (self.time > self.disk_lifetime);
    }

    fn update_time(&mut self, time_in_seconds: f64) {
        self.time = time_in_seconds;
    }

    // Update routine for a Star Planet simulation.
    fn update_star_planet(&mut self, new_state: &UniverseIntegral) -> Result<()> {
        // ***WARNING!***
        // Stateful function
        // The order of these calculations is important.
        // Lower order calculations depend on previous values.
        // ***WARNING!***
        let &mut ParticleType::Star(ref mut star) = &mut self.central_body.kind else {
            unreachable!()
        };
        let &mut ParticleType::Planet(ref mut planet) = &mut self.orbiting_body.kind else {
            unreachable!()
        };
        // Update radiative zone and convective zone angular momentum
        // and recompute independent values.
        star.refresh(
            self.time,
            new_state.central_body.radiative_zone_angular_momentum,
            new_state.central_body.convective_zone_angular_momentum,
            self.disk_is_dissipated,
        )?;

        // Invert the exponent of sma^6.5 to normalise the semi major axis.
        // Recompute planet values, including those depending on star.
        planet.refresh(new_state.orbiting_body.semi_major_axis.powf(2. / 13.), star);

        // No torques during disk lifetime.
        if !self.disk_is_dissipated {
            return Ok(());
        }

        // Compute the stellar wind torque (if enabled).
        star.update_wind_torque(self.central_body.wind.wind_torque());

        // No planetary torques after the planet is destroyed.
        if planet.is_destroyed() {
            return Ok(());
        }

        // Recompute star values that depend on planet (tidal and magnetic torque).
        star.refresh_tidal_frequency(planet);

        // Compute the enabled effects dependent on a planet (magnetism, planet tides)
        star.update_tidal_torque(self.central_body.tides.tidal_torque(star, planet));
        star.update_magnetic_torque(self.central_body.magnetism.magnetic_torque(planet, star)); // Requires wind torque to be calculated first.
        star.update_evolved_change_semi_major_axis(self.central_body.wind.wind_torque(), planet);

        if self.orbiting_body.tides.kaula_enabled() {
            //(spin, eccentricity, inclination, longitude_ascending_node, pericentre_omega, spin_inclination)
            // Invert the exponent of e^2 to normalise the eccentricity.
            planet.refresh_orbital_elements(
                new_state.orbiting_body.spin,
                sqrt!(new_state.orbiting_body.eccentricity),
                new_state.orbiting_body.inclination,
                new_state.orbiting_body.longitude_ascending_node,
                new_state.orbiting_body.pericentre_omega,
                new_state.orbiting_body.spin_inclination,
            );
            // Recompute the kaula tidal effects.
            self.orbiting_body
                .tides
                .refresh_kaula(self.time, star, planet)?;
        } else if self.central_body.general_relativity.is_enabled() {
            // General relativity only evolves pericentre_omega
            // If general relativity and kaula tides on the planet are both enabled,
            // pericentre_omega is updated in `planet.refresh_orbital_elements` for kaula
            planet.pericentre_omega = new_state.orbiting_body.pericentre_omega;
        }

        Ok(())
    }

    // Update the planet and star values from the integrator prior to the derivation step.
    pub(crate) fn update(&mut self, time: f64, new_state: &UniverseIntegral) -> Result<()> {
        // Set the time.
        self.update_time(time);
        // Calculate the dissipation status of the disk.
        self.disk_is_dissipated();

        match self.universe_kind() {
            UniverseKind::StarPlanet => self.update_star_planet(new_state)?,
            UniverseKind::BinaryStar => todo!(),
            UniverseKind::PlanetMoon => todo!(),
        }

        Ok(())
    }

    pub(crate) fn clear_destroyed_particles(&mut self) {
        match self.universe_kind() {
            UniverseKind::StarPlanet => {
                let &mut ParticleType::Planet(ref mut planet) = &mut self.orbiting_body.kind else {
                    unreachable!()
                };
                planet.destroy();
            }
            UniverseKind::BinaryStar => todo!(),
            UniverseKind::PlanetMoon => todo!(),
        }
    }
}

#[cfg(test)]
pub mod tests;
