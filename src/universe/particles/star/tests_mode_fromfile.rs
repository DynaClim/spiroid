use super::*;
use crate::universe::effects::magnetism::{IsothermalWind, MagneticModel};
use crate::universe::effects::tides::ConstantTimeLag;
use crate::universe::effects::tides::constant_time_lag::Equilibrium;
use crate::universe::effects::tides::constant_time_lag::Inertial;
use crate::universe::particles::TidalModel;
use crate::universe::particles::planet::tests::{test_planet, test_planet_magnetic};
use crate::universe::tests::{DISK_IS_DISSIPATED, TEST_TIME_mode_fromfile};

use pretty_assertions::assert_eq;
use sci_file::read_csv_rows_from_file;

// Same structure and tests as `tests_mode_ardestani`, for the `FromFile` turnover-time mode and
// a 1 Msun star. The structural values of `test_star` are those of the 2026 STAREVOL track of
// Louis Amard (starevol_m1p00_2026.csv) at the solar age, 4.567 Gyr; the star rotates at the
// solar rate. NOTE: that track is not tracked by git yet.

fn add_interpolate_to_test_star(star: &mut Star) {
    star.evolution = Evolution::Starevol {
        star_file_path: "examples/data/star/evolution/starevol_m1p00_2026.csv".into(),
        convective_turnover_time_mode: TurnoverTimeMode::FromFile,
        interpolator: Interpolator1D::new(),
    };
    // Load stellar evolution data from file.
    if let Some(star_file_path) = star.evolution_file() {
        let mut stellar_data = read_csv_rows_from_file::<StarCsv>(star_file_path).unwrap();
        // Configure the stellar evolution interpolator.
        StarCsv::initialise(&mut stellar_data);
        let star_ages = StarCsv::ages(&stellar_data);
        star.initialise_evolution(&star_ages, &stellar_data)
            .unwrap();
    }
}

pub fn test_star_evolving() -> Star {
    let mut star = test_star();
    add_interpolate_to_test_star(&mut star);
    star
}

pub fn test_star() -> Star {
    let mut star = Star::default();

    star.mass = 1.9884158605722266e30;
    // Gallet & Bouvier 2013 slow rotator (30 Myr).
    star.core_envelope_coupling_constant = 946728e6;
    star.radius = 713711897.4193548;
    star.radiative_mass = 1.9638428890823805e30;
    star.convective_radius = 538312237.8387097;

    star.convective_moment_of_inertia = 5.710775131281403e45;
    star.radiative_moment_of_inertia = 6.201607952241131e46;
    star.radiative_mass_derivative = -255255910999.14395;

    star.convective_turnover_time_sun = CONVECTIVE_TURNOVER_TIME_SUN_STAREVOL_2026;
    star.wind_torque_prefactor = 8e23;
    // Solar rotation rate, 2 pi / 25.38 d.
    star.spin = 2.8653290845717256e-6;
    star.angular_momentum_redistribution = -3.260979477440655e24;

    // Both zones corotating at the solar rate: spin * moment of inertia.
    let radiative_zone_angular_momentum = 1.7769647636667813e41;
    let convective_zone_angular_momentum = 1.6363250079109518e40;
    star.refresh(
        TEST_TIME_mode_fromfile,
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
    let expected = -3.260979477440655e24;
    let star = test_star();
    let result = star.angular_momentum_redistribution();
    assert_eq!(expected, result);
}

#[test]
fn _wind_torque() {
    let expected = -7.368258251123122e23;
    let star = test_star();
    let result = star.wind_torque();
    assert_eq!(expected, result);
}

#[test]
fn _evolved_wind_torque() {
    let expected = -2.784131017369365e28;
    let mut star = test_star();
    star.evolved_mass_loss_rate = 2.8612812361645612e16;
    let result = star.evolved_wind_torque();
    assert_eq!(expected, result);
}

#[test]
fn _mass_transfer_envelope_to_core_torque() {
    let expected = -1.4129525029179846e23;
    let star = test_star();
    let result = star.mass_transfer_envelope_to_core_torque();
    assert_eq!(expected, result);
}

#[test]
fn _mass_loss_rate() {
    let expected = 1233085092.120401;
    let star = test_star();
    let result = star.mass_loss_rate();
    assert_eq!(expected, result);
}

#[test]
fn _alfven_radius_estimate() {
    let expected = 14441045684.524355;
    let star = test_star();
    let result = star.alfven_radius_estimate();
    assert_eq!(expected, result);
}

#[test]
fn _dynamical_tide_dissipation() {
    let expected = 1300.7613258982299;
    let star = test_star();
    let result = star.dynamical_tide_dissipation();
    assert_eq!(expected, result);
}

#[test]
fn _rossby() {
    let expected = 1.3411832983570864;
    let star = test_star();
    let result = star.rossby();
    assert_eq!(expected, result);
}

#[test]
fn _convective_turnover_time() {
    // The turnover time `refresh` stores for the evolving star in FromFile mode: the file column
    // (tauc_hp of the 2026 track) interpolated at TEST_TIME_mode_fromfile, not the Ardestani fit.
    // Calling test_star_evolving is necessary for reading convective_turnover_time from the file.
    let expected = 1754129.196;
    let mut star = test_star_evolving();
    star.refresh(
        TEST_TIME_mode_fromfile,
        star.radiative_zone_angular_momentum,
        star.convective_zone_angular_momentum,
        DISK_IS_DISSIPATED,
    )
    .unwrap();
    let result = star.convective_turnover_time;
    assert_eq!(expected, result);
    let fit = Star::convective_turnover_time((star.mass - star.radiative_mass) / star.mass);
    assert_ne!(fit, result);
}

#[test]
fn _tidal_frequency() {
    let expected = -0.00024291434004431378;
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let result = Star::tidal_frequency(&star, &planet);
    assert_eq!(expected, result);
}

#[test]
fn _magnetic_torque_enabled() {
    let expected = 8.362759311054455e22;
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;
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
    let expected = 9.330821869632422e24;
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
