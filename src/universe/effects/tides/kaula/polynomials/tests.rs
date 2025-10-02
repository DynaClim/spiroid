use super::*;
use pretty_assertions::assert_eq;

#[cfg(test)]
pub fn test_polynomials() -> Polynomials {
    Polynomials {
        inclination_2mp_squared: [
            [
                0.0019441123903527048,
                0.169592269793094,
                0.0019441123903527048,
            ],
            [
                0.24875653835984152,
                0.2334467410593072,
                0.0002431018222873341,
            ],
            [
                7.957335142227162,
                0.031105798245643277,
                7.599675858850527e-6,
            ],
        ],

        inclination_2mp_squared_derivative: [
            [
                0.021303678127986697,
                -0.3979485529163214,
                0.021303678127986697,
            ],
            [
                1.2749785238785618,
                1.1086309462370176,
                0.0037763044283222027,
            ],
            [
                -5.6277016901034065,
                0.34085885004778715,
                0.00017192998808024731,
            ],
        ],

        eccentricity_2pq_squared: [
            [
                8.951155653314441e-35,
                1.9290123456790128e-30,
                3.910782952816371e-26,
                6.7819213887956415e-22,
                6.781917154315517e-18,
                0.0,
                6.249960937764488e-6,
                0.9998750049218077,
                0.00030621636852872246,
                4.515434092363666e-8,
                4.84170951316148e-12,
                4.334291873777955e-16,
                3.452767780162516e-20,
                2.5351301409874433e-24,
                1.7519717957478673e-28,
            ],
            [
                1.2075026123849163e-30,
                2.3913118839263917e-26,
                4.684250683722252e-22,
                9.047011818090832e-18,
                1.7145172690212558e-13,
                3.1641855519440455e-9,
                5.62531642025811e-5,
                1.0000750037501562,
                5.62531642025811e-5,
                3.1641855519440455e-9,
                1.7145172690212558e-13,
                9.047011818090832e-18,
                4.684250683722252e-22,
                2.3913118839263917e-26,
                1.2075026123849163e-30,
            ],
            [
                1.7519717957478673e-28,
                2.5351301409874433e-24,
                3.452767780162516e-20,
                4.334291873777955e-16,
                4.84170951316148e-12,
                4.515434092363666e-8,
                0.00030621636852872246,
                0.9998750049218077,
                6.249960937764488e-6,
                0.0,
                6.781917154315517e-18,
                6.7819213887956415e-22,
                3.910782952816371e-26,
                1.9290123456790128e-30,
                8.951155653314441e-35,
            ],
        ],

        eccentricity_2pq_squared_derivative: [
            [
                250632357.9314715,
                4629629.62962963,
                29744241.200387478,
                547010.6336805556,
                1957820.7411699824,
                0.0,
                4120.284288194444,
                -869213.6284722221,
                -38894276.42578125,
                163066561.1979167,
                18004583056534.527,
                126315129241.94336,
                462890276611781.3,
                6084312338369.864,
                490552102809402.75,
            ],
            [
                3381007314677.764,
                57391485214.2334,
                5608456177.740459,
                166040649.4140625,
                104574627275.08397,
                8260539.0625,
                46132312.21758541,
                9.04224537037037e-5,
                46132312.21758541,
                8260539.0625,
                104574627275.08397,
                166040649.4140625,
                5608456177.740459,
                57391485214.2334,
                3381007314677.764,
            ],
            [
                490552102809402.75,
                6084312338369.864,
                462890276611781.3,
                126315129241.94336,
                18004583056534.527,
                163066561.1979167,
                -38894276.42578125,
                -869213.6284722221,
                4120.284288194444,
                0.0,
                1957820.7411699824,
                547010.6336805556,
                29744241.200387478,
                4629629.62962963,
                250632357.9314715,
            ],
        ],
    }
}

#[test]
fn _inclination_polynomials() {
    let mut polynomials = Polynomials::default();
    let expected = test_polynomials();
    let inclination: f64 = 0.35;
    polynomials.inclination_polynomials(cos!(inclination), sin!(inclination));
    assert_eq!(
        expected.inclination_2mp_squared,
        polynomials.inclination_2mp_squared
    );
}

#[test]
fn _inclination_polynomials_derivatives() {
    let mut polynomials = Polynomials::default();
    let expected = test_polynomials();
    let inclination: f64 = 0.35;
    polynomials.inclination_polynomials_derivatives(cos!(inclination), sin!(inclination));
    assert_eq!(
        expected.inclination_2mp_squared_derivative,
        polynomials.inclination_2mp_squared_derivative
    );
}

#[test]
fn _eccentricity_polynomials() {
    let mut polynomials = Polynomials::default();
    let expected = test_polynomials();
    let eccentricity = 0.005;
    polynomials.eccentricity_polynomials(eccentricity);
    assert_eq!(
        expected.eccentricity_2pq_squared,
        polynomials.eccentricity_2pq_squared
    );
}

#[test]
fn _eccentricity_polynomials_derivatives() {
    let mut polynomials = Polynomials::default();
    let expected = test_polynomials();
    let eccentricity = 5.;
    polynomials.eccentricity_polynomials_derivatives(eccentricity);
    assert_eq!(
        expected.eccentricity_2pq_squared_derivative,
        polynomials.eccentricity_2pq_squared_derivative
    );
}
