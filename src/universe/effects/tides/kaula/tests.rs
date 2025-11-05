use super::*;
use crate::universe::effects::tides::kaula::love_number::tests::test_love_number;
use crate::universe::effects::tides::kaula::polynomials::tests::test_polynomials;
use crate::universe::particles::planet::tests::test_planet_kaula;
use crate::universe::particles::star::tests::test_star;
use crate::universe::tests::TEST_TIME;

use pretty_assertions::assert_eq;

#[cfg(test)]
pub fn test_kaula() -> Kaula {
    let mut kaula = Kaula {
        particle_type: ParticleComposition::Solid {
            solid_file: "dummy".into(),
        },
        polynomials: test_polynomials(),
        love_number: test_love_number(),
        summation: Summation::default(),

        prev_mean_motion: f64::NAN,
        prev_spin: f64::NAN,
        prev_inclination: f64::NAN,
        prev_eccentricity: f64::NAN,
    };
    let mpq = test_mpq();
    let planet = test_planet_kaula();

    kaula.summation.real_2pq_2mp_dt = kaula.sum_over_m_real(
        &kaula.polynomials.eccentricity_2pq_squared,
        &kaula.polynomials.inclination_2mp_squared_derivative,
        mpq,
    );
    kaula.summation.real_2pq_dt_2mp = kaula.sum_over_m_real(
        &kaula.polynomials.eccentricity_2pq_squared_derivative,
        &kaula.polynomials.inclination_2mp_squared,
        mpq,
    );

    kaula.summation.imaginary_mfactor = kaula.sum_over_m_imaginary_mfactor(mpq);
    kaula.summation.imaginary_pfactor = kaula.sum_over_m_imaginary_pfactor(mpq);
    kaula.summation.imaginary_qfactor = kaula.sum_over_m_imaginary_qfactor(mpq);

    kaula.summation.imaginary_inclination = kaula.sum_over_m_imaginary_inclination(&planet, mpq);
    kaula.summation.imaginary_eccentricity = kaula.sum_over_m_imaginary_eccentricity(&planet, mpq);

    kaula
}

pub(crate) fn test_mpq() -> Mpq {
    Mpq {
        m_min: 0,
        m_max: 3,
        p_min: 0,
        p_max: 3,
        q_min: 0,
        q_max: 15,
    }
}

#[test]
fn _refresh_high_ecc() {
    let mut kaula = test_kaula();
    let star = test_star();
    let mut planet = test_planet_kaula();
    planet.eccentricity = 0.3;

    kaula.refresh(TEST_TIME, &planet, &star).unwrap();
    let expected = Summation {
        real_2pq_2mp_dt: -2.4457628278616497e-6,
        real_2pq_dt_2mp: -0.00910296147744826,
        imaginary_mfactor: -0.002293444338943969,
        imaginary_pfactor: -0.002442067516179405,
        imaginary_qfactor: -0.0043923613097381084,
        imaginary_eccentricity: -0.0017479781239668566,
        imaginary_inclination: 64.60961289373019,
    };
    let result = kaula.summation;
    assert_eq!(expected, result);
}

#[test]
fn _refresh_low_ecc() {
    let mut kaula = test_kaula();
    let star = test_star();
    let mut planet = test_planet_kaula();
    planet.inclination = 0.0;
    planet.refresh(planet.semi_major_axis, &star);
    kaula.refresh(TEST_TIME, &planet, &star).unwrap();
    let expected = Summation {
        real_2pq_2mp_dt: 0.017022814801684244,
        real_2pq_dt_2mp: -5.2982137061585994e-5,
        imaginary_mfactor: -0.0007413858810131647,
        imaginary_pfactor: 0.5629831296068434,
        imaginary_qfactor: -0.0007415623709616942,
        imaginary_eccentricity: -1.672203609571526e-7,
        imaginary_inclination: -7423.755336647935,
    };
    let result = kaula.summation;
    assert_eq!(expected, result);
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
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_eccentricity();
    let expected = 7.73439335780939e-5;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_inclination() {
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_inclination();
    let expected = -7423.755336647935;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_longitude_ascending_node() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_longitude_ascending_node(&planet);
    let expected = 5.076723278969285e-33;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_spin_axis_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_spin_axis_inclination(&planet);
    let expected = -0.05924449315462055;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_pericentre_eccentricity() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_pericentre_eccentricity(&planet);
    let expected = -3.2159666539736226e-25;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_pericentre_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_pericentre_inclination(&planet);
    let expected = -5.258897087031897e-34;
    assert_eq!(expected, result);
}
