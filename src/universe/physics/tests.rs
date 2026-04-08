use super::*;
use crate::Universe;
use crate::universe::effects::general_relativity::GRModel;
use crate::universe::effects::magnetism::{IsothermalWind, MagneticModel};
use crate::universe::effects::tides::ConstantTimeLag;
use crate::universe::effects::tides::TidalModel;
use crate::universe::effects::tides::constant_time_lag::Equilibrium;
use crate::universe::effects::tides::constant_time_lag::Inertial;
use crate::universe::effects::tides::kaula::tests::test_kaula;
use crate::universe::effects::wind::WindModel;
use crate::universe::particles::planet::tests::{
    test_planet, test_planet_kaula, test_planet_magnetic,
};
use crate::universe::particles::star::tests::{test_star, test_star_evolving};
use crate::universe::tests::{DISK_IS_DISSIPATED, TEST_DISK_LIFETIME, TEST_TIME};
use crate::universe::{Particle, ParticleType};
use pretty_assertions::assert_eq;

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
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Wind(wind),
            wind: WindModel::Enabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };
    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -6.348994811695528e22;
    expected.central_body.convective_zone_angular_momentum = 1.6351930535408648e22;
    expected.orbiting_body.semi_major_axis = -1.4634701453519956e43;

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
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::ConstantTimeLag(ConstantTimeLag {
                equilibrium: Equilibrium::SigmaBarStar(1e-6),
                inertial: Inertial::FrequencyAveraged,
            }),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Enabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -6.348994811695528e22;
    expected.central_body.convective_zone_angular_momentum = 6.020027165936562e23;
    expected.orbiting_body.semi_major_axis = -1.9848639097150575e44;

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
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::ConstantTimeLag(ConstantTimeLag {
                equilibrium: Equilibrium::SigmaBarStar(1e-6),
                inertial: Inertial::FrequencyAveraged,
            }),
            magnetism: MagneticModel::Wind(wind),
            wind: WindModel::Enabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -6.348994811695528e22;
    expected.central_body.convective_zone_angular_momentum = 6.486208599049978e23;
    expected.orbiting_body.semi_major_axis = -2.131210924250257e44;

    assert_eq!(expected, result);
}

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
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Enabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.orbiting_body.semi_major_axis = 3.0436830855775857e49;
    expected.orbiting_body.spin = -1.5432986377090881e-9;
    expected.orbiting_body.eccentricity = 5.129250165061513e-16;
    expected.orbiting_body.inclination = 0.0007011714730044864;
    expected.orbiting_body.longitude_ascending_node = -0.0018601514572935667;
    expected.orbiting_body.pericentre_omega = 5.876804234917257e-6;
    expected.orbiting_body.spin_inclination = 0.0003527171243323208;

    assert_eq!(expected, result);
}

#[test]
fn _star_radiative_zone_angular_momentum_derivative() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = star_radiative_zone_angular_momentum_derivative(&star);
    let expected = -6.348994811695822e22;
    assert_eq!(expected, result);
}

#[test]
fn _star_convective_zone_angular_momentum_derivative() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = star_convective_zone_angular_momentum_derivative(&star, DISK_IS_DISSIPATED);
    let expected = -3.0079737689074846e22;
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
    let expected = -2.1307551258578705e44;
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
    let expected = -2.3102674836008928e52;
    assert_eq!(expected, result);
}

#[test]
fn _planet_spin_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_derivative(&planet, &star, &kaula);
    let expected = 1.2501842317327892e-6;
    assert_eq!(expected, result);
}

#[test]
fn _planet_eccentricity_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_eccentricity_derivative(&planet, &star, &kaula);
    let expected = -2.371917363949444e-13;
    assert_eq!(expected, result);
}

#[test]
fn _planet_inclination_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.26702870058883815;
    assert_eq!(expected, result);
}

#[test]
fn _planet_longitude_ascending_node_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_longitude_ascending_node_derivative(&planet, &star, &kaula);
    let expected = 0.6637879161156786;
    assert_eq!(expected, result);
}

#[test]
fn _planet_argument_pericentre_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_argument_pericentre_derivative(&planet, &star, &kaula);
    let expected = -42049166.159453586;
    assert_eq!(expected, result);
}

#[test]
fn _planet_spin_axis_inclination_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_axis_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.1641163464603492;
    assert_eq!(expected, result);
}

// // Unit test for the GR precession rate function.
// // dω/dt = 3 n³ a² / (c² (1 - e²))
// #[test]
// fn _gr_pericentre_precession_rate() {
//     let planet = test_planet_kaula();
//     let result = gr_pericentre_precession_rate(&planet);
//     let expected = todo!(); // fill in: 3 * n^3 * a^2 / (c^2 * (1 - e^2))
//     assert_eq!(expected, result);
// }

// // Integration test: GR enabled, Kaula disabled.
// // Only pericentre_omega should be nonzero among orbital element derivatives.
// #[test]
// fn _derivatives_gr_only() {
//     let star = test_star();
//     let planet = test_planet_kaula();
//
//     let mut y = UniverseIntegral::default();
//     y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);
//     y.orbiting_body.pericentre_omega = planet.pericentre_omega;
//
//     let mut universe = Universe {
//         orbiting_body: Particle {
//             kind: ParticleType::Planet(planet),
//             tides: TidalModel::Disabled,
//             magnetism: MagneticModel::Disabled,
//             wind: WindModel::Disabled,
//         },
//         central_body: Particle {
//             kind: ParticleType::Star(star),
//             tides: TidalModel::Disabled,
//             magnetism: MagneticModel::Disabled,
//             wind: WindModel::Enabled,
//             general_relativity: GRModel::Enabled,
//         },
//         perturbing_body: None,
//         time: TEST_TIME,
//         disk_lifetime: TEST_DISK_LIFETIME,
//         disk_is_dissipated: DISK_IS_DISSIPATED,
//         derivatives: UniverseIntegral::default(),
//     };
//
//     universe.update(TEST_TIME, &y).unwrap();
//     let mut result = UniverseIntegral::default();
//     let _ = force(
//         &universe.central_body,
//         &universe.orbiting_body,
//         universe.disk_is_dissipated,
//         &mut result,
//     )
//     .unwrap();
//
//     let mut expected = UniverseIntegral::default();
//     expected.orbiting_body.pericentre_omega = todo!(); // fill in GR precession rate
//     // All other orbital element derivatives should be zero:
//     // spin, eccentricity, inclination, longitude_ascending_node, spin_inclination = 0
//     assert_eq!(expected, result);
// }

// // Integration test: both Kaula and GR enabled.
// // pericentre_omega derivative should be the sum of Kaula and GR contributions.
// #[test]
// fn _derivatives_kaula_and_gr() {
//     let star = test_star();
//     let planet = test_planet_kaula();
//
//     let mut y = UniverseIntegral::default();
//     y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);
//     y.orbiting_body.spin = 8.062093352143078e-7;
//     y.orbiting_body.eccentricity = 2.500000000179822e-5;
//     y.orbiting_body.inclination = 0.34999207817863753;
//     y.orbiting_body.longitude_ascending_node = 1.0465602799892118;
//     y.orbiting_body.pericentre_omega = -0.11536773671287792;
//     y.orbiting_body.spin_inclination = 0.31581363067032314;
//
//     let mut universe = Universe {
//         orbiting_body: Particle {
//             kind: ParticleType::Planet(planet),
//             tides: TidalModel::KaulaTides(test_kaula()),
//             magnetism: MagneticModel::Disabled,
//             wind: WindModel::Disabled,
//         },
//         central_body: Particle {
//             kind: ParticleType::Star(star),
//             tides: TidalModel::Disabled,
//             magnetism: MagneticModel::Disabled,
//             wind: WindModel::Enabled,
//             general_relativity: GRModel::Enabled,
//         },
//         perturbing_body: None,
//         time: TEST_TIME,
//         disk_lifetime: TEST_DISK_LIFETIME,
//         disk_is_dissipated: DISK_IS_DISSIPATED,
//         derivatives: UniverseIntegral::default(),
//     };
//
//     universe.update(TEST_TIME, &y).unwrap();
//     let mut result = UniverseIntegral::default();
//     let _ = force(
//         &universe.central_body,
//         &universe.orbiting_body,
//         universe.disk_is_dissipated,
//         &mut result,
//     )
//     .unwrap();
//
//     let mut expected = UniverseIntegral::default();
//     // Same as _derivatives_kaula except pericentre_omega += GR contribution:
//     expected.orbiting_body.semi_major_axis = 3.0436830855775857e49;
//     expected.orbiting_body.spin = -1.5432986377090881e-9;
//     expected.orbiting_body.eccentricity = 5.129250165061513e-16;
//     expected.orbiting_body.inclination = 0.0007011714730044864;
//     expected.orbiting_body.longitude_ascending_node = -0.0018601514572935667;
//     expected.orbiting_body.pericentre_omega = todo!(); // fill in: kaula value (5.876804234917257e-6) + GR value
//     expected.orbiting_body.spin_inclination = 0.0003527171243323208;
//     assert_eq!(expected, result);
// }
