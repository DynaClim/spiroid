pub(crate) mod effects;
pub(crate) mod particles;
pub(crate) mod physics;

pub use effects::Kaula;
pub use particles::{Particle, ParticleType, Planet, Star, StarCsv};

use anyhow::Result;
use serde::{Deserialize, Serialize};

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
    disk_is_dissipated: bool,
    pub orbiting_body: Particle,
    pub central_body: Particle,
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

    // Creates a vector of initial quantities to be integrated for the particle.
    fn integration_quantities_per_particle(particle: &Particle) -> Vec<f64> {
        let mut vec = vec![];

        match &particle.kind {
            ParticleType::Star(star) => vec.append(&mut vec![
                star.spin * star.radiative_moment_of_inertia,
                star.spin * star.convective_moment_of_inertia,
            ]),
            ParticleType::Planet(planet) => {
                // The actual integrated quantity is sma^6.5.
                // See comment for planet_semi_major_axis_13_div_2_derivative
                vec.append(&mut vec![planet.semi_major_axis.powf(6.5)]);

                if particle.tides.kaula_enabled() {
                    vec.append(&mut vec![
                        planet.spin,
                        // The actual integrated quantity is eccentricity^2 to avoid singularities when eccentricity goes to 0.
                        planet.eccentricity.powi(2),
                        planet.inclination,
                        planet.longitude_ascending_node,
                        planet.pericentre_omega,
                        planet.spin_inclination,
                    ]);
                }
            }
        }

        vec
    }

    // Creates a vector of initial quantities to be integrated, depending on the simulation configuration.
    pub fn integration_quantities(&self) -> Vec<f64> {
        let mut vec = vec![];

        vec.append(&mut Self::integration_quantities_per_particle(
            &self.central_body,
        ));
        vec.append(&mut Self::integration_quantities_per_particle(
            &self.orbiting_body,
        ));

        vec
    }

    fn disk_is_dissipated(&mut self) {
        self.disk_is_dissipated = self.disk_is_dissipated || (self.time > self.disk_lifetime);
    }

    fn update_time(&mut self, time_in_seconds: f64) {
        self.time = time_in_seconds;
    }

    // Update routine for a Star Planet simulation.
    fn update_star_planet(&mut self, y: &[f64]) -> Result<()> {
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
        // Update radiative zone (y[0]) and convective zone (y[1]) angular momentum
        // and recompute independent values.
        star.refresh(self.time, y[0], y[1], self.disk_is_dissipated)?;

        // Invert the exponent of sma^6.5 to normalise the semi major axis.
        // Recompute planet values, including those depending on star.
        planet.refresh(y[2].powf(2. / 13.), star);

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
            planet.refresh_orbital_elements(y[3], sqrt!(y[4]), y[5], y[6], y[7], y[8]);
            // Recompute the kaula tidal effects.
            self.orbiting_body
                .tides
                .refresh_kaula(self.time, star, planet)?;
        }

        Ok(())
    }

    // Update the planet and star values from the integrator prior to the derivation step.
    pub(crate) fn update(&mut self, time: f64, y: &[f64]) -> Result<()> {
        // Set the time.
        self.update_time(time);
        // Calculate the dissipation status of the disk.
        self.disk_is_dissipated();

        match self.universe_kind() {
            UniverseKind::StarPlanet => self.update_star_planet(y)?,
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
