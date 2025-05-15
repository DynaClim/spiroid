use super::*;
use crate::universe::effects::tides::kaula::love_number::tests::test_love_number;
use crate::universe::effects::tides::kaula::polynomials::tests::test_polynomials;
use crate::universe::particles::planet::tests::test_planet_kaula;

use pretty_assertions::assert_eq;

#[cfg(test)]
pub fn test_kaula() -> Kaula {
    Kaula {
        particle_type: ParticleComposition::Solid {
            solid_file: "dummy".into(),
        },
        polynomials: test_polynomials(),
        love_number: test_love_number(),
    }
}

#[test]
fn _summation_of_longitudinal_modes_semi_major_axis() {
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_semi_major_axis();
    let expected = 0.5630675119283106;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_spin() {
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_spin();
    let expected = 0.5641312760456983;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_eccentricity() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_eccentricity(planet.eccentricity);
    let expected = 7.73439335780939e-5;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_inclination(planet.inclination, 2., 3.);
    let expected = 0.03973365404734674;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_longitude_ascending_node() {
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_longitude_ascending_node(2., 3., 4.);
    let expected = 3.9613491613661527;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_spin_axis_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_spin_axis_inclination(
        planet.longitude_ascending_node,
        planet.inclination,
    );
    let expected = -0.05924449315462055;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_pericentre_eccentricity() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_pericentre_eccentricity(planet.eccentricity);
    let expected = -137114722399835.89;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_pericentre_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_pericentre_inclination(
        planet.inclination,
        planet.spin_inclination,
    );
    let expected = 0.008511407400842122;
    assert_eq!(expected, result);
}
