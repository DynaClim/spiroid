use super::*;
use crate::universe::effects::magnetism::{IsothermalWind, MagneticModel};
use crate::universe::effects::tides::ConstantTimeLag;
use crate::universe::effects::tides::constant_time_lag::Equilibrium;
use crate::universe::effects::tides::constant_time_lag::Inertial;
use crate::universe::particles::TidalModel;
use crate::universe::particles::planet::tests::{test_planet, test_planet_magnetic};
use crate::universe::tests::{DISK_IS_DISSIPATED, TEST_TIME};

use pretty_assertions::assert_eq;
use sci_file::read_csv_rows_from_file;

fn add_interpolate_to_test_star(star: &mut Star) {
    star.evolution = Evolution::Starevol {
        star_file_path: "examples/data/star/evolution/savgol_08.csv".into(),
        interpolator: Interpolator::new(),
    };
    // Load stellar evolution data from file.
    if let Some(star_file_path) = star.evolution_file() {
        let mut stellar_data = read_csv_rows_from_file::<StarCsv>(star_file_path).unwrap();
        // Configure the stellar evolution interpolator.
        let (star_ages, star_values) = StarCsv::initialise(&mut stellar_data);
        star.initialise_evolution(&star_ages, &star_values);
    }
}

pub fn test_star_evolving() -> Star {
    let mut star = test_star();
    add_interpolate_to_test_star(&mut star);
    star
}

pub fn test_star() -> Star {
    let mut star = Star::default();

    star.mass = 1.5909177014856084e30;
    star.core_envelope_coupling_constant = 369539496.0e6;
    // The tolerance is very high: any value in range 4.7e5 to 9.3e5 will pass the existing tests.
    star.radius = 544588072.4685764;
    star.radiative_mass = 1.5048623991131647e30;
    star.convective_radius = 374606632.43479675;

    star.convective_moment_of_inertia = 1.0420137774656348e46;
    star.radiative_moment_of_inertia = 3.605010133078022e46;
    star.radiative_mass_derivative = 23849190556.112328;

    star.convective_turnover_time_sun = Star::convective_turnover_time(0.02);
    star.spin = 1.2583862403723232e-6;
    star.angular_momentum_redistribution = 2.4499591272215565e37;

    //    star.convective_moment_of_inertia_derivative = 2.380403432967787e27;
    //    star.dynamical_tide_dissipation = 3037.223911926055;

    let radiative_zone_angular_momentum = 4.547421124942826e40;
    let convective_zone_angular_momentum = 1.3112557998411429e40;
    star.refresh(
        TEST_TIME,
        radiative_zone_angular_momentum,
        convective_zone_angular_momentum,
        DISK_IS_DISSIPATED,
    )
    .unwrap();

    star.update_wind_torque(true);

    star
}

// Tests below

#[test]
fn _angular_momentum_redistribution() {
    let expected = 2.4499591272215565e37;
    let star = test_star();
    let result = star.angular_momentum_redistribution();
    assert_eq!(expected, result);
}

#[test]
fn _wind_torque() {
    let expected = -9.356968580603306e22;
    let star = test_star();
    let result = star.wind_torque();
    assert_eq!(expected, result);
}

#[test]
fn _evolved_wind_torque() {
    let expected = -7.119008297630572e27;
    let mut star = test_star();
    star.evolved_mass_loss_rate = 2.8612812361645612e16;
    let result = star.evolved_wind_torque();
    assert_eq!(expected, result);
}

#[test]
fn _mass_transfer_envelope_to_core_torque() {
    let expected = 2.80767781316852e21;
    let star = test_star();
    let result = star.mass_transfer_envelope_to_core_torque();
    assert_eq!(expected, result);
}

#[test]
fn _mass_loss_rate() {
    let expected = 439989563.7273058;
    let star = test_star();
    let result = star.mass_loss_rate();
    assert_eq!(expected, result);
}

#[test]
fn _alfven_radius_estimate() {
    let expected = 12999882215.232397;
    let star = test_star();
    let result = star.alfven_radius_estimate();
    assert_eq!(expected, result);
}

#[test]
fn _dynamical_tide_dissipation() {
    let expected = 3034.8061299412066;
    let star = test_star();
    let result = star.dynamical_tide_dissipation();
    assert_eq!(expected, result);
}

#[test]
fn _rossby() {
    let expected = 1.2930654606068248;
    let star = test_star();
    let result = star.rossby();
    assert_eq!(expected, result);
}

#[test]
fn _convective_turnover_time() {
    let adjusted_convective_mass = 0.02;
    let expected = 2126270.90231897;
    let result = Star::convective_turnover_time(adjusted_convective_mass);
    assert_eq!(expected, result);
}

#[test]
fn _tidal_frequency() {
    let expected = -0.00024612822573271255;
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let result = Star::tidal_frequency(&star, &planet);
    assert_eq!(expected, result);
}

#[test]
fn _magnetic_torque_enabled() {
    let expected = 4.648379104022687e22;
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7.0e4;
    let mut magnetism = MagneticModel::Wind(wind);

    let result = magnetism.magnetic_torque(&planet, &star);
    assert_eq!(expected, result);
}

#[test]
// No magnetic torque if magnetism is disabled.
fn _magnetic_torque_disabled() {
    let expected = 0.0;
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let mut magnetism = MagneticModel::Disabled;
    star.magnetic_torque = magnetism.magnetic_torque(&planet, &star);

    let result = magnetism.magnetic_torque(&planet, &star);
    assert_eq!(expected, result);
}

#[test]
// No tidal torque if tides are disabled.
fn _tidal_torque_disabled() {
    let expected = 0.0;
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let tides = TidalModel::Disabled;
    let result = tides.tidal_torque(&star, &planet);
    assert_eq!(expected, result);
}

#[test]
fn _tidal_torque_enabled() {
    let expected = 6.325284391272144e23;
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let tides = TidalModel::ConstantTimeLag(ConstantTimeLag {
        equilibrium: Equilibrium::SigmaBarStar(1e-6),
        inertial: Inertial::FrequencyAveraged,
    });
    let result = tides.tidal_torque(&star, &planet);
    assert_eq!(expected, result);
}
