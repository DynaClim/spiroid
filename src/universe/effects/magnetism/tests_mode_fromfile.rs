use super::*;
use crate::universe::particles::planet::tests::test_planet_magnetic;
use crate::universe::particles::star::tests_mode_fromfile::test_star;

use pretty_assertions::assert_eq;

// Same tests as `tests_mode_ardestani` with the FromFile-mode test star (1 Msun at the solar
// age, solar reference CONVECTIVE_TURNOVER_TIME_SUN_STAREVOL_2026, see star::tests_mode_fromfile).

#[test]
fn _init_weber_davis() {
    let expected = IsothermalWind {
        // input
        footpoint_conductance: 7e4,
        // intermediate
        speed_of_sound: 156665.93675243173,
        critical_radius: 2703538045.6409526,
        critical_radius_div_alfven_radius: 0.18721206931283238,
        radial_magnetic_field: 2.6001826742856518e-5,
        magnetic_pressure: 0.0002690096507232475,
        integration_constant: 0.44669805455408906,
        wind_velocity: 0.2530473913021395,
        surface_wind_velocity: 0.013080451485249047,
        wind_density: 8.103557070876562e-16,
        alfvenic_mach: 0.3404900406229801,
        azimuthal_velocity: 5447.233627343446,
        alfven_speed_at_alfven_radius: 395413.6412459326,
        interaction: MagneticInteraction::Dipolar,
        // output
        magnetic_torque: 8.362759311054455e22,
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
    let expected = 2.6001826742856518e-5;
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
    let expected = 0.0002690096507232475;
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
    let expected = 8.103557070876562e-16;
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
    let expected = 0.3404900406229801;
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
    let expected = 0.44669805455408906;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let result = wind.integration_constant(&star);
    assert_eq!(expected, result);
}

#[test]
fn _weber_davis_velocity_profile() {
    let expected = 0.013080451485249047;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let result = wind.weber_davis_velocity_profile(star.radius, &star);
    assert_eq!(expected, result);
}

#[test]
fn _alfven_speed_at_alfven_radius() {
    let expected = 395413.6412459326;
    let star = test_star();
    let planet = test_planet_magnetic();
    let mut wind = IsothermalWind::default();
    wind.init_weber_davis(&planet, &star);
    let result = wind.alfven_speed_at_alfven_radius(&star);
    assert_eq!(expected, result);
}

#[test]
fn _magnetic_torque() {
    let expected = 8.362759311054455e22;
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
    let expected = 0.00018448496386013329;
    let star = test_star();
    let result = IsothermalWind::magnetic_field(star.mass, star.rossby, star.rossby_sun_code());
    assert_eq!(expected, result);
}
