use crate::constants::MAGNETIC_PERMEABILITY_OF_VACUUM;
use crate::universe::effects::tides::TidalModel;
use crate::universe::effects::{MagneticModel, WindModel};
pub(crate) mod planet;
pub(crate) mod star;

pub use planet::Planet;
pub use star::{Star, StarCsv};

use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub enum ParticleType {
    Planet(Planet),
    Star(Star),
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct Particle {
    pub kind: ParticleType,
    #[serde(default)]
    pub tides: TidalModel,
    #[serde(default)]
    pub(crate) magnetism: MagneticModel,
    #[serde(default)]
    pub(crate) wind: WindModel,
}

impl Particle {
    pub(crate) fn is_star(&self) -> bool {
        matches!(self.kind, ParticleType::Star(_))
    }

    pub(crate) fn is_planet(&self) -> bool {
        matches!(self.kind, ParticleType::Planet(_))
    }

    pub(crate) fn initialise(&mut self, time: f64) -> Result<()> {
        match &mut self.kind {
            ParticleType::Star(star) => star.initialise(time)?,
            ParticleType::Planet(planet) => planet.initialise(),
        }

        Ok(())
    }
}

// Common properties of both Star and Planet.
// Enables making functions generic over impl ParticleT.
pub trait ParticleT {
    fn semi_major_axis(&self) -> f64;
    fn mass(&self) -> f64;
    fn radius(&self) -> f64;
    fn spin(&self) -> f64;
    fn spin_inclination(&self) -> f64;
    fn eccentricity(&self) -> f64;
    fn inclination(&self) -> f64;
    fn luminosity(&self) -> f64;
    fn mean_motion(&self) -> f64;
    fn moment_of_inertia(&self) -> f64;
    fn reduced_mass(&self) -> f64;
}

// https://en.wikipedia.org/wiki/Magnetic_pressure
pub(crate) fn magnetic_pressure(magnetic_field: f64) -> f64 {
    magnetic_field.powi(2) / (2. * MAGNETIC_PERMEABILITY_OF_VACUUM)
}
