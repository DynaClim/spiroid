use super::*;
use crate::universe::particles::planet::tests::test_planet_magnetic;
use crate::universe::particles::star::tests_mode_ardestani::test_star;

use pretty_assertions::assert_eq;

#[test]
fn _init_weber_davis() {
    let expected = IsothermalWind {
        // input
        footpoint_conductance: 7e4,
        // intermediate
        speed_of_sound: 154526.12347777415,
        critical_radius: 2223403684.8240876,
        critical_radius_div_alfven_radius: 0.15847855533221863,
        radial_magnetic_field: 2.226577824572609e-5,
        magnetic_pressure: 0.00019725857851166478,
        integration_constant: 0.4553635280344644,
        wind_velocity: 0.3094177219897195,
        surface_wind_velocity: 0.008056559716731736,
        wind_density: 1.5205555914799212e-16,
        alfvenic_mach: 0.16498899414711332,
        azimuthal_velocity: 2392.2989794266414,
        alfven_speed_at_alfven_radius: 407510.82539476344,
        interaction: MagneticInteraction::Dipolar,
        // output
        magnetic_torque: 5.068552639495255e22,
    };
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;
    wind.init_weber_davis(&planet, &star);
    assert_eq!(expected, wind);
}

#[test]
fn _radial_magnetic_field() {
    let expected = 2.226577824572609e-5;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let surface_magnetic_field =
        IsothermalWind::magnetic_field(star.mass, star.rossby, star.rossby_sun_code());
    let result = IsothermalWind::radial_magnetic_field(
        surface_magnetic_field,
        star.radius,
        planet.semi_major_axis,
    );
    assert_eq!(expected, result);
}

#[test]
fn _magnetic_pressure() {
    let expected = 0.00019725857851166478;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let surface_magnetic_field =
        IsothermalWind::magnetic_field(star.mass, star.rossby, star.rossby_sun_code());
    let magnetic_field = IsothermalWind::radial_magnetic_field(
        surface_magnetic_field,
        star.radius,
        planet.semi_major_axis,
    );
    let result = magnetic_pressure(magnetic_field);
    assert_eq!(expected, result);
}

#[test]
fn _density_profile() {
    let expected = 1.5205555914799212e-16;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let coronal_density =
        IsothermalWind::coronal_density(star.mass, star.rossby, star.rossby_sun_code());

    let result = wind.density_profile(star.radius, coronal_density, planet.semi_major_axis);
    assert_eq!(expected, result);
}

#[test]
fn _alfvenic_mach() {
    let expected = 0.16498899414711332;
    let star = test_star();
    let planet = test_planet_magnetic();
    let keplerian_velocity = sqrt!(GRAVITATIONAL * star.mass / planet.semi_major_axis);
    let star = test_star();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let result = wind.alfvenic_mach(keplerian_velocity);
    assert_eq!(expected, result);
}

#[test]
fn _integration_constant() {
    let expected = 0.4553635280344644;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let result = wind.integration_constant(&star);
    assert_eq!(expected, result);
}

#[test]
fn _weber_davis_velocity_profile() {
    let expected = 0.008056559716731736;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let result = wind.weber_davis_velocity_profile(star.radius, &star);
    assert_eq!(expected, result);
}

#[test]
fn _alfven_speed_at_alfven_radius() {
    let expected = 407510.82539476344;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let result = wind.alfven_speed_at_alfven_radius(&star);
    assert_eq!(expected, result);
}

#[test]
fn _magnetic_torque() {
    let expected = 5.068552639495255e22;
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;
    wind.init_weber_davis(&planet, &star);
    let result = wind.magnetic_torque(&planet, &star);
    assert_eq!(expected, result);
}

#[test]
fn _magnetic_field_magnetic() {
    let expected = 0.000271334249863209;
    let star = test_star();
    let result = IsothermalWind::magnetic_field(star.mass, star.rossby, star.rossby_sun_code());
    assert_eq!(expected, result);
}
