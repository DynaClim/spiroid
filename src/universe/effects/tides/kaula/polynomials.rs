use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub(crate) struct Polynomials {
    pub(crate) inclination_2mp_squared: [[f64; 3]; 3],
    pub(crate) inclination_2mp_squared_derivative: [[f64; 3]; 3],
    pub(crate) eccentricity_2pq_squared: [[f64; 15]; 3],
    pub(crate) eccentricity_2pq_squared_derivative: [[f64; 15]; 3],
}

impl Polynomials {
    pub(crate) fn refresh_cache(&mut self, eccentricity: f64, inclination: f64) {
        let cos_inc = cos!(inclination);
        let sin_inc = sin!(inclination);
        self.inclination_polynomials(cos_inc, sin_inc);
        self.inclination_polynomials_derivatives(cos_inc, sin_inc);
        self.eccentricity_polynomials(eccentricity);
        self.eccentricity_polynomials_derivatives(eccentricity);
    }

    // Inclination polynomials defined as Kaula (1964)
    fn inclination_polynomials(&mut self, cos_inc: f64, sin_inc: f64) {
        let sin_inc_2 = sin_inc.powi(2);

        let [m0p, m1p, m2p] = &mut self.inclination_2mp_squared;

        m0p[0] = -(3. / 8.) * sin_inc_2;
        m0p[1] = (3. / 4.) * sin_inc_2 - (1. / 2.);
        m0p[2] = m0p[0];

        for val in m0p.iter_mut() {
            *val = val.powi(2);
        }

        m1p[0] = (3. / 4.) * sin_inc * (1. + cos_inc);
        m1p[1] = -(3. / 2.) * sin_inc * cos_inc;
        m1p[2] = (3. / 4.) * sin_inc * (cos_inc - 1.);

        for val in m1p.iter_mut() {
            *val = val.powi(2);
        }

        m2p[0] = (3. / 4.) * (1. + cos_inc).powi(2);
        m2p[1] = (3. / 2.) * sin_inc_2;
        m2p[2] = (3. / 4.) * (1. - cos_inc).powi(2);

        for val in m2p.iter_mut() {
            *val = val.powi(2);
        }
    }

    // Partial derivative of inclination polynomials (defined as Kaula 1964) with respect to the inclination.
    fn inclination_polynomials_derivatives(&mut self, cos_inc: f64, sin_inc: f64) {
        let [m0p, m1p, m2p] = &mut self.inclination_2mp_squared_derivative;

        m0p[0] = (9. / 16.) * sin_inc.powi(3) * cos_inc;
        m0p[1] = ((3. / 2.) * sin_inc.powi(2) - 1.) * (3. / 2.) * cos_inc * sin_inc;
        m0p[2] = m0p[0];

        m1p[0] =
            (9. / 8.) * sin_inc * (1. + cos_inc) * (cos_inc * (1. + cos_inc) - sin_inc.powi(2));
        m1p[1] = (9. / 2.) * sin_inc * cos_inc * (cos_inc.powi(2) - sin_inc.powi(2));
        m1p[2] = (9. / 8.) * sin_inc * cos_inc * (cos_inc - 1.).powi(2)
            - sin_inc.powi(3) * (cos_inc - 1.);

        m2p[0] = -(9. / 4.) * sin_inc * (1. + cos_inc).powi(3);
        m2p[1] = 9. * sin_inc.powi(3) * cos_inc;
        m2p[2] = (9. / 4.) * sin_inc * (1. - cos_inc).powi(3);
    }

    // Partial derivative of eccentricity polynomials (see Kaula 1964, defined in the table of Cayley 1861) with respect to the eccentricity.
    #[allow(clippy::similar_names)]
    fn eccentricity_polynomials_derivatives(&mut self, eccentricity: f64) {
        let ecc_2 = eccentricity.powi(2);
        let ecc_3 = eccentricity.powi(3);
        let ecc_5 = eccentricity.powi(5);
        let ecc_7 = eccentricity.powi(7);
        let ecc_8 = eccentricity.powi(8);
        let ecc_9 = eccentricity.powi(9);
        let ecc_11 = eccentricity.powi(11);
        let ecc_13 = eccentricity.powi(13);

        let [p0q, p1q, p2q] = &mut self.eccentricity_2pq_squared_derivative;

        p0q[0] = ECC_DY_P0Q0 * ecc_13;
        p0q[1] = ECC_DY_P0Q1 * ecc_11;
        p0q[2] = ECC_DY_P0Q2_A * ecc_9 + ECC_DY_P0Q2_B * ecc_11 + ECC_DY_P0Q2_C * ecc_13;
        p0q[3] = ECC_DY_P0Q3_A * ecc_7 + ECC_DY_P0Q3_B * ecc_9 + ECC_DY_P0Q3_C * ecc_11;
        p0q[4] = ECC_DY_P0Q4_A * ecc_5
            + ECC_DY_P0Q4_B * ecc_7
            + ECC_DY_P0Q4_C * ecc_9
            + ECC_DY_P0Q4_D * ecc_11
            + ECC_DY_P0Q4_E * ecc_13;
        //p0q[5] = 0.;
        p0q[6] = ECC_DY_P0Q6_A * eccentricity
            + ECC_DY_P0Q6_B * ecc_3
            + ECC_DY_P0Q6_C * ecc_5
            + ECC_DY_P0Q6_D * ecc_7;
        p0q[7] = ECC_DY_P0Q7_A * eccentricity
            + ECC_DY_P0Q7_B * ecc_3
            + ECC_DY_P0Q7_C * ecc_5
            + ECC_DY_P0Q7_D * ecc_7;

        p0q[8] = ECC_DY_P0Q8_A * eccentricity
            + ECC_DY_P0Q8_B * ecc_3
            + ECC_DY_P0Q8_C * ecc_5
            + ECC_DY_P0Q8_D * ecc_7;
        p0q[9] = ECC_DY_P0Q9_A * ecc_3 + ECC_DY_P0Q9_B * ecc_5 + ECC_DY_P0Q9_C * ecc_7;
        p0q[10] = ECC_DY_P0Q10_A * ecc_5
            + ECC_DY_P0Q10_B * ecc_7
            + ECC_DY_P0Q10_C * ecc_9
            + ECC_DY_P0Q10_D * ecc_11
            + ECC_DY_P0Q10_E * ecc_13;
        p0q[11] = ECC_DY_P0Q11_A * ecc_7 + ECC_DY_P0Q11_B * ecc_8 + ECC_DY_P0Q11_C * ecc_9;
        p0q[12] = ECC_DY_P0Q12_A * ecc_9 + ECC_DY_P0Q12_B * ecc_11 + ECC_DY_P0Q12_C * ecc_13;
        p0q[13] = ECC_DY_P0Q13 * ecc_11;
        p0q[14] = ECC_DY_P0Q14 * ecc_13;

        p1q[0] = ECC_DY_P1Q0 * ecc_13;
        p1q[1] = ECC_DY_P1Q1 * ecc_11;
        p1q[2] = ECC_DY_P1Q2_A * ecc_9 + ECC_DY_P1Q2_B * ecc_11 + ECC_DY_P1Q2_C * ecc_13;
        p1q[3] = ECC_DY_P1Q3_A * ecc_7 + ECC_DY_P1Q3_B * ecc_9;
        p1q[4] = ECC_DY_P1Q4_A * ecc_5
            + ECC_DY_P1Q4_B * ecc_7
            + ECC_DY_P1Q4_C * ecc_9
            + ECC_DY_P1Q4_D * ecc_11
            + ECC_DY_P1Q4_E * ecc_13;
        p1q[5] = ECC_DY_P1Q5_A * ecc_3 + ECC_DY_P1Q5_B * ecc_5 + ECC_DY_P1Q5_C * ecc_7;
        p1q[6] = ECC_DY_P1Q6_A * eccentricity
            + ECC_DY_P1Q6_B * ecc_3
            + ECC_DY_P1Q6_C * ecc_5
            + ECC_DY_P1Q6_D * ecc_7;

        p1q[7] = ECC_DY_P1Q7 * eccentricity * (1. - ecc_2).powi(-4);
        // array reflected here
        let (left, right) = p1q.split_at_mut(8);
        right.copy_from_slice(&left[0..=6]);
        right.reverse();

        // Copy in reverse order
        p2q.clone_from(p0q);
        p2q.reverse();
    }

    // Eccentricity polynomials (see Kaula 1964, defined in the table of Cayley 1861).
    fn eccentricity_polynomials(&mut self, eccentricity: f64) {
        let ecc = eccentricity;
        let ecc_2 = eccentricity.powi(2);
        let ecc_3 = eccentricity.powi(3);
        let ecc_4 = eccentricity.powi(4);
        let ecc_5 = eccentricity.powi(5);
        let ecc_6 = eccentricity.powi(6);
        let ecc_7 = eccentricity.powi(7);

        let [p0q, p1q, p2q] = &mut self.eccentricity_2pq_squared;

        p0q[0] = ECC_P0Q0 * ecc_7;
        p0q[1] = ECC_P0Q1 * ecc_6;
        p0q[2] = ECC_P0Q2_A * ecc_5 + ECC_P0Q2_B * ecc_7;
        p0q[3] = ECC_P0Q3_A * ecc_4 + ECC_P0Q3_B * ecc_6;
        p0q[4] = ECC_P0Q4_A * ecc_3 + ECC_P0Q4_B * ecc_5 + ECC_P0Q4_C * ecc_7;
        // p0q[5] stays as 0.;
        p0q[6] = ECC_P0Q6_A * ecc + ECC_P0Q6_B * ecc_3 + ECC_P0Q6_C * ecc_5 + ECC_P0Q6_D * ecc_7;
        p0q[7] = 1. + ECC_P0Q7_A * ecc_2 + ECC_P0Q7_B * ecc_4 + ECC_P0Q7_C * ecc_6;
        p0q[8] = ECC_P0Q8_A * ecc + ECC_P0Q8_B * ecc_3 + ECC_P0Q8_C * ecc_5 + ECC_P0Q8_D * ecc_7;
        p0q[9] = ECC_P0Q9_A * ecc_2 + ECC_P0Q9_B * ecc_4 + ECC_P0Q9_C * ecc_6;
        p0q[10] = ECC_P0Q10_A * ecc_3 + ECC_P0Q10_B * ecc_5 + ECC_P0Q10_C * ecc_7;
        p0q[11] = ECC_P0Q11_A * ecc_4 + ECC_P0Q11_B * ecc_6;
        p0q[12] = ECC_P0Q12_A * ecc_5 + ECC_P0Q12_B * ecc_7;
        p0q[13] = ECC_P0Q13 * ecc_6;
        p0q[14] = ECC_P0Q14 * ecc_7;

        for val in p0q.iter_mut() {
            *val = val.powi(2);
        }

        p1q[0] = ECC_P1Q0 * ecc_7;
        p1q[1] = ECC_P1Q1 * ecc_6;
        p1q[2] = ECC_P1Q2_A * ecc_5 + ECC_P1Q2_B * ecc_7;
        p1q[3] = ECC_P1Q3_A * ecc_4 + ECC_P1Q3_B * ecc_6;
        p1q[4] = ECC_P1Q4_A * ecc_3 + ECC_P1Q4_B * ecc_5 + ECC_P1Q4_C * ecc_7;
        p1q[5] = ECC_P1Q5_A * ecc_2 + ECC_P1Q5_B * ecc_4 + ECC_P1Q5_C * ecc_6;
        p1q[6] = ECC_P1Q6_A * ecc + ECC_P1Q6_B * ecc_3 + ECC_P1Q6_C * ecc_5 + ECC_P1Q6_D * ecc_7;
        p1q[7] = (1. - ecc_2).powf(ECC_P1Q7);
        // array refelected here
        let (left, right) = p1q.split_at_mut(8);
        right.copy_from_slice(&left[0..=6]);
        right.reverse();

        for val in p1q.iter_mut() {
            *val = val.powi(2);
        }

        // Copy in reverse order
        p2q.clone_from(p0q);
        p2q.reverse();
    }
}

// Eccentricity factors for polynomials
// for p = 0, q = [-7, 7] as [0, 14]
const ECC_P0Q0: f64 = 15625. / 129_024.;

const ECC_P0Q1: f64 = 4. / 45.;

const ECC_P0Q2_A: f64 = 81. / 1280.;
const ECC_P0Q2_B: f64 = 81. / 2048.;

const ECC_P0Q3_A: f64 = 1. / 24.;
const ECC_P0Q3_B: f64 = 7. / 240.;

const ECC_P0Q4_A: f64 = 1. / 48.;
const ECC_P0Q4_B: f64 = 11. / 768.;
const ECC_P0Q4_C: f64 = 313. / 30_720.;

const ECC_P0Q6_A: f64 = -1. / 2.;
const ECC_P0Q6_B: f64 = 1. / 16.;
const ECC_P0Q6_C: f64 = -5. / 384.;
const ECC_P0Q6_D: f64 = -143. / 18_432.;

const ECC_P0Q7_A: f64 = -5. / 2.;
const ECC_P0Q7_B: f64 = 13. / 16.;
const ECC_P0Q7_C: f64 = -35. / 288.;

const ECC_P0Q8_A: f64 = 7. / 2.;
const ECC_P0Q8_B: f64 = -123. / 16.;
const ECC_P0Q8_C: f64 = 489. / 128.;
const ECC_P0Q8_D: f64 = -1763. / 2048.;

const ECC_P0Q9_A: f64 = 17. / 2.;
const ECC_P0Q9_B: f64 = -115. / 16.;
const ECC_P0Q9_C: f64 = 601. / 48.;

const ECC_P0Q10_A: f64 = 845. / 48.;
const ECC_P0Q10_B: f64 = -32525. / 768.;
const ECC_P0Q10_C: f64 = 208_225. / 6144.;

const ECC_P0Q11_A: f64 = 533. / 16.;
const ECC_P0Q11_B: f64 = -13827. / 160.;

const ECC_P0Q12_A: f64 = 228_347. / 3840.;
const ECC_P0Q12_B: f64 = -3_071_075. / 18432.;

const ECC_P0Q13: f64 = 73369. / 720.;
const ECC_P0Q14: f64 = 12_144_273. / 71680.;

// for p = 1, q = [-7, 7] as [0, 14]
const ECC_P1Q0: f64 = 432_091. / 30720.;

const ECC_P1Q1: f64 = 3167. / 320.;

const ECC_P1Q2_A: f64 = 1773. / 256.;
const ECC_P1Q2_B: f64 = 4987. / 6144.;

const ECC_P1Q3_A: f64 = 77. / 16.;
const ECC_P1Q3_B: f64 = 129. / 160.;

const ECC_P1Q4_A: f64 = 53. / 16.;
const ECC_P1Q4_B: f64 = 393. / 256.;
const ECC_P1Q4_C: f64 = 24753. / 10240.;

const ECC_P1Q5_A: f64 = 9. / 4.;
const ECC_P1Q5_B: f64 = 7. / 4.;
const ECC_P1Q5_C: f64 = 141. / 64.;

const ECC_P1Q6_A: f64 = 3. / 2.;
const ECC_P1Q6_B: f64 = 27. / 16.;
const ECC_P1Q6_C: f64 = 261. / 128.;
const ECC_P1Q6_D: f64 = 14309. / 6144.;

const ECC_P1Q7: f64 = -3. / 2.;

// Eccentricity derivatives for polynomials
// for p = 0, q = [-7, 7] as [0, 14]
const ECC_DY_P0Q0: f64 = 3_417_968_750. / 16_647_192_600.;

const ECC_DY_P0Q1: f64 = 64. / 675.;

const ECC_DY_P0Q2_A: f64 = 6561. / 163_840.;
const ECC_DY_P0Q2_B: f64 = 19683. / 327_680.;
const ECC_DY_P0Q2_C: f64 = 45927. / 2_097_152.;

const ECC_DY_P0Q3_A: f64 = 1. / 72.;
const ECC_DY_P0Q3_B: f64 = 7. / 288.;
const ECC_DY_P0Q3_C: f64 = 49. / 4800.;

const ECC_DY_P0Q4_A: f64 = 1. / 384.;
const ECC_DY_P0Q4_B: f64 = 11. / 2304.;
const ECC_DY_P0Q4_C: f64 = 619. / 98304.;
const ECC_DY_P0Q4_D: f64 = 3443. / 983_040.;
const ECC_DY_P0Q4_E: f64 = 685_783. / 471_859_200.;

const ECC_DY_P0Q6_A: f64 = 1. / 2.;
const ECC_DY_P0Q6_B: f64 = -1. / 4.;
const ECC_DY_P0Q6_C: f64 = 13. / 128.;
const ECC_DY_P0Q6_D: f64 = 113. / 2304.;

const ECC_DY_P0Q7_A: f64 = -10.;
const ECC_DY_P0Q7_B: f64 = 63. / 2.;
const ECC_DY_P0Q7_C: f64 = -155. / 6.;
const ECC_DY_P0Q7_D: f64 = -2921. / 288.;

const ECC_DY_P0Q8_A: f64 = 49. / 2.;
const ECC_DY_P0Q8_B: f64 = -861. / 4.;
const ECC_DY_P0Q8_C: f64 = 65925. / 128.;
const ECC_DY_P0Q8_D: f64 = -132_635. / 256.;

const ECC_DY_P0Q9_A: f64 = 289.;
const ECC_DY_P0Q9_B: f64 = -5865. / 8.;
const ECC_DY_P0Q9_C: f64 = 203_147. / 96.;

const ECC_DY_P0Q10_A: f64 = 714_025. / 384.;
const ECC_DY_P0Q10_B: f64 = -27_483_625. / 2304.;
const ECC_DY_P0Q10_C: f64 = 2_936_126_875. / 98_304.;
const ECC_DY_P0Q10_D: f64 = -6_772_518_125. / 196_608.;
const ECC_DY_P0Q10_E: f64 = 303_503_554_375. / 18_874_368.;

const ECC_DY_P0Q11_A: f64 = 284_089. / 32.;
const ECC_DY_P0Q11_B: f64 = -66_328_119. / 1280.;
const ECC_DY_P0Q11_C: f64 = 191_185_929. / 2560.;

const ECC_DY_P0Q12_A: f64 = 52_142_352_409. / 1_474_560.;
const ECC_DY_P0Q12_B: f64 = -140_254_152_605. / 589_824.;
const ECC_DY_P0Q12_C: f64 = 66_020_511_589_375. / 169_869_312.;

const ECC_DY_P0Q13: f64 = 5_383_010_161. / 43200.;
const ECC_DY_P0Q14: f64 = 147_483_366_698_529. / 367_001_600.;

// for p = 1, q = [-7, 7] as [0, 14]
const ECC_DY_P1Q0: f64 = 1_306_918_425_967. / 471_859_200.;

const ECC_DY_P1Q1: f64 = 30_089_667. / 25600.;

const ECC_DY_P1Q2_A: f64 = 15_717_645. / 32768.;
const ECC_DY_P1Q2_B: f64 = -8_841_951. / 65536.;
const ECC_DY_P1Q2_C: f64 = 174_091_183. / 18_874_368.;

const ECC_DY_P1Q3_A: f64 = 5929. / 32.;
const ECC_DY_P1Q3_B: f64 = 9933. / 128.;

const ECC_DY_P1Q4_A: f64 = 8427. / 128.;
const ECC_DY_P1Q4_B: f64 = 20829. / 256.;
const ECC_DY_P1Q4_C: f64 = 6_019_881. / 32768.;
const ECC_DY_P1Q4_D: f64 = 29_183_787. / 327_680.;
const ECC_DY_P1Q4_E: f64 = 4_288_977_063. / 52_428_800.;

const ECC_DY_P1Q5_A: f64 = 81. / 4.;
const ECC_DY_P1Q5_B: f64 = 189. / 4.;
const ECC_DY_P1Q5_C: f64 = 1661. / 16.;

const ECC_DY_P1Q6_A: f64 = 9. / 2.;
const ECC_DY_P1Q6_B: f64 = 81. / 4.;
const ECC_DY_P1Q6_C: f64 = 6885. / 128.;
const ECC_DY_P1Q6_D: f64 = 12_123_879. / 20608.;

const ECC_DY_P1Q7: f64 = 6.;

#[cfg(test)]
pub mod tests;

// References:
// Cayley 1861 https://ui.adsabs.harvard.edu/abs/1861MmRAS..29..191C
// Kaula 1964 10.1029/RG002i004p00661
