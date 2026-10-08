use super::*;
use crate::Universe;
use crate::universe::effects::general_relativity::GeneralRelativityModel;
use crate::universe::effects::magnetism::{IsothermalWind, MagneticModel};
use crate::universe::effects::tides::ConstantTimeLag;
use crate::universe::effects::tides::TidalModel;
use crate::universe::effects::tides::constant_time_lag::Equilibrium;
use crate::universe::effects::tides::constant_time_lag::Inertial;
use crate::universe::effects::tides::kaula::tests::test_kaula;
use crate::universe::effects::wind::WindModel;
use crate::universe::particles::planet::tests::{
    test_planet, test_planet_kaula, test_planet_magnetic, test_planet_mercury,
};
use crate::universe::particles::star::tests_mode_fromfile::{test_star, test_star_evolving};
use crate::universe::tests::{DISK_IS_DISSIPATED, TEST_DISK_LIFETIME, TEST_TIME_mode_fromfile};
use crate::universe::{Particle, ParticleType};
use pretty_assertions::assert_eq;

// Same tests as `tests_mode_ardestani` for the FromFile turnover-time mode: the 1 Msun star of
// `star::tests_mode_fromfile` (with the 2026 STAREVOL track for the evolving tests), evaluated at
// TEST_TIME_mode_fromfile.

#[test]
fn _derivatives_magnetic() {
    let star = test_star_evolving();
    let planet = test_planet_magnetic();

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);

    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Wind(wind),
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME_mode_fromfile,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };
    universe.update(TEST_TIME_mode_fromfile, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -3.892224542501776e27;
    expected.central_body.convective_zone_angular_momentum = 3.891815235879581e27;
    expected.orbiting_body.semi_major_axis = -3.4278571902979133e43;

    assert_eq!(expected, result);
}

#[test]
fn _derivatives_tides() {
    let star = test_star_evolving();
    let planet = test_planet();

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::ConstantTimeLag(ConstantTimeLag {
                equilibrium: Equilibrium::SigmaBarStar(1e-6),
                inertial: Inertial::FrequencyAveraged,
            }),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME_mode_fromfile,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME_mode_fromfile, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -3.892224542501776e27;
    expected.central_body.convective_zone_angular_momentum = 4.0381625899274537e27;
    expected.orbiting_body.semi_major_axis = -4.1111956874199373e46;

    assert_eq!(expected, result);
}

#[test]
fn _derivatives_magnetic_tides() {
    let star = test_star_evolving();
    let planet = test_planet_magnetic();

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);

    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::ConstantTimeLag(ConstantTimeLag {
                equilibrium: Equilibrium::SigmaBarStar(1e-6),
                inertial: Inertial::FrequencyAveraged,
            }),
            magnetism: MagneticModel::Wind(wind),
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME_mode_fromfile,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME_mode_fromfile, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -3.892224542501776e27;
    expected.central_body.convective_zone_angular_momentum = 4.038284714120102e27;
    expected.orbiting_body.semi_major_axis = -4.114623544610235e46;

    assert_eq!(expected, result);
}

// Note: Tests are currently tailored for amd64 and may fail (with slight numerical discrepancies)
// on other platforms due to differences in cpu architecture.
#[test]
fn _derivatives_kaula() {
    let star = test_star();
    let planet = test_planet_kaula();

    let mut y = UniverseIntegral::default();
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);
    y.orbiting_body.spin = 8.062093352143078e-7;
    y.orbiting_body.eccentricity = 2.500000000179822e-5;
    y.orbiting_body.inclination = 0.34999207817863753;
    y.orbiting_body.longitude_ascending_node = 1.0465602799892118;
    y.orbiting_body.pericentre_omega = -0.11536773671287792;
    y.orbiting_body.spin_inclination = 0.31581363067032314;

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::KaulaTides(test_kaula()),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME_mode_fromfile,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME_mode_fromfile, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.orbiting_body.semi_major_axis = 4.757609277202951e49;
    expected.orbiting_body.spin = -2.6970753910788047e-9;
    expected.orbiting_body.eccentricity = 8.015556953720482e-16;
    expected.orbiting_body.inclination = 0.0012249201319543006;
    expected.orbiting_body.longitude_ascending_node = -0.003249149194701831;
    expected.orbiting_body.pericentre_omega = 9.142954948700943e-6;
    expected.orbiting_body.spin_inclination = 0.0006158499873966261;

    assert_eq!(expected, result);
}

#[test]
fn _star_radiative_zone_angular_momentum_derivative() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = star_radiative_zone_angular_momentum_derivative(&star);
    let expected = -1.4129525028835398e23;
    assert_eq!(expected, result);
}

#[test]
fn _star_convective_zone_angular_momentum_derivative() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = star_convective_zone_angular_momentum_derivative(&star, DISK_IS_DISSIPATED);
    let expected = -5.955305748239582e23;
    assert_eq!(expected, result);
}

#[test]
fn _planet_semi_major_axis_13_div_2_derivative() {
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let tides = TidalModel::ConstantTimeLag(ConstantTimeLag {
        equilibrium: Equilibrium::SigmaBarStar(1e-6),
        inertial: Inertial::FrequencyAveraged,
    });
    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;
    let mut magnetism = MagneticModel::Wind(wind);
    let tidal_torque_convective = tides.tidal_torque(&star, &planet);
    let magnetic_torque = magnetism.magnetic_torque(&planet, &star);
    let wind_torque = WindModel::Enabled.wind_torque();

    star.update_wind_torque(wind_torque);
    star.update_tidal_torque(tidal_torque_convective);
    star.update_magnetic_torque(magnetic_torque);

    let result = planet_semi_major_axis_13_div_2_derivative(&planet, &star);
    let expected = -2.6425057626752656e45;
    assert_eq!(expected, result);
}

#[test]
fn _kaula_planet_semi_major_axis_13_div_2_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);

    let tides = TidalModel::ConstantTimeLag(ConstantTimeLag {
        equilibrium: Equilibrium::SigmaBarStar(1e-6),
        inertial: Inertial::FrequencyAveraged,
    });
    let mut magnetism = MagneticModel::Wind(IsothermalWind::default());
    let tidal_torque_convective = tides.tidal_torque(&star, &planet);
    let magnetic_torque = magnetism.magnetic_torque(&planet, &star);
    let wind_torque = WindModel::Enabled.wind_torque();

    star.update_wind_torque(wind_torque);
    star.update_tidal_torque(tidal_torque_convective);
    star.update_magnetic_torque(magnetic_torque);

    let kaula = test_kaula();

    let result = kaula_planet_semi_major_axis_13_div_2_derivative(&planet, &star, &kaula);
    let expected = -3.2280952749067313e52;
    assert_eq!(expected, result);
}

#[test]
fn _planet_spin_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_derivative(&planet, &star, &kaula);
    let expected = 1.9529585509471158e-6;
    assert_eq!(expected, result);
}

#[test]
fn _planet_eccentricity_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_eccentricity_derivative(&planet, &star, &kaula);
    let expected = -3.3142375458188136e-13;
    assert_eq!(expected, result);
}

#[test]
fn _planet_inclination_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.3337470586838154;
    assert_eq!(expected, result);
}

#[test]
fn _planet_longitude_ascending_node_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_longitude_ascending_node_derivative(&planet, &star, &kaula);
    let expected = 1.0369274014892227;
    assert_eq!(expected, result);
}

#[test]
fn _planet_argument_pericentre_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_argument_pericentre_derivative(&planet, &star, &kaula);
    let expected = -65686541.65272909;
    assert_eq!(expected, result);
}

#[test]
fn _planet_spin_axis_inclination_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_axis_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.2563721522272757;
    assert_eq!(expected, result);
}

#[test]
fn _general_relativity_pericentre_precession_rate() {
    let planet = test_planet_mercury();
    let result = general_relativity_pericentre_precession_rate(&planet);
    let expected = 4.725582990934671e-14;
    assert_eq!(expected, result);
}
