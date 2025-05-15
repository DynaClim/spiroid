use crate::constants::MAGNETIC_PERMEABILITY_OF_VACUUM;
use crate::universe::effects::MagneticModel;
use crate::universe::effects::tides::TidalModel;
pub(crate) mod planet;
pub(crate) mod star;

pub use planet::Planet;
pub use star::{Star, StarCsv};

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
    pub magnetism: MagneticModel,
}

// Common properties of both Star and Planet.
// Enables making functions generic over impl ParticleT.
pub trait ParticleT {
    fn semi_major_axis(&self) -> f64;
    fn mass(&self) -> f64;
    fn radius(&self) -> f64;
    fn spin(&self) -> f64;
    fn eccentricity(&self) -> f64;
    fn inclination(&self) -> f64;
    fn luminosity(&self) -> f64;
    fn mean_motion(&self) -> f64;
}

// https://en.wikipedia.org/wiki/Magnetic_pressure
pub fn magnetic_pressure(magnetic_field: f64) -> f64 {
    magnetic_field.powi(2) / (2. * MAGNETIC_PERMEABILITY_OF_VACUUM)
}
