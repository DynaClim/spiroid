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
        self.inclination_polynomials(inclination);
        self.inclination_polynomials_derivatives(inclination);
        self.eccentricity_polynomials(eccentricity);
        self.eccentricity_polynomials_derivatives(eccentricity);
    }

    // Inclination polynomials defined as Kaula (1964)
    fn inclination_polynomials(&mut self, inclination: f64) {
        let cos_inc = cos!(inclination);
        let sin_inc = sin!(inclination);
        let sin_inc_2 = sin_inc.powi(2);

        let borrowed_slice = &mut self.inclination_2mp_squared[..];
        let [first, second, third] = borrowed_slice else {
            unreachable!()
        };

        first[0] = -(3. / 8.) * sin_inc_2;
        first[1] = (3. / 4.) * sin_inc_2 - (1. / 2.);
        first[2] = first[0];

        for val in first.iter_mut() {
            *val = val.powi(2);
        }

        second[0] = (3. / 4.) * sin_inc * (1. + cos_inc);
        second[1] = -(3. / 2.) * sin_inc * cos_inc;
        second[2] = (3. / 4.) * sin_inc * (cos_inc - 1.);

        for val in second.iter_mut() {
            *val = val.powi(2);
        }

        third[0] = (3. / 4.) * (1. + cos_inc).powi(2);
        third[1] = (3. / 2.) * sin_inc_2;
        third[2] = (3. / 4.) * (1. - cos_inc).powi(2);

        for val in third.iter_mut() {
            *val = val.powi(2);
        }
    }

    // Partial derivative of inclination polynomials (defined as Kaula 1964) with respect to the inclination.
    fn inclination_polynomials_derivatives(&mut self, inclination: f64) {
        let cos_inc = cos!(inclination);
        let sin_inc = sin!(inclination);

        let borrowed_slice = &mut self.inclination_2mp_squared_derivative[..];
        let [first, second, third] = borrowed_slice else {
            unreachable!()
        };

        first[0] = (9. / 16.) * sin_inc.powi(3) * cos_inc;
        first[1] = ((3. / 2.) * sin_inc.powi(2) - 1.) * (3. / 2.) * cos_inc * sin_inc;
        first[2] = first[0];

        second[0] =
            (9. / 8.) * sin_inc * (1. + cos_inc) * (cos_inc * (1. + cos_inc) - sin_inc.powi(2));
        second[1] = (9. / 2.) * sin_inc * cos_inc * (cos_inc.powi(2) - sin_inc.powi(2));
        second[2] = (9. / 8.) * sin_inc * cos_inc * (cos_inc - 1.).powi(2)
            - sin_inc.powi(3) * (cos_inc - 1.);

        third[0] = -(9. / 4.) * sin_inc * (1. + cos_inc).powi(3);
        third[1] = 9. * sin_inc.powi(3) * cos_inc;
        third[2] = (9. / 4.) * sin_inc * (1. - cos_inc).powi(3);
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

        let borrowed_slice = &mut self.eccentricity_2pq_squared_derivative[..];
        let [first, second, third] = borrowed_slice else {
            unreachable!()
        };

        first[0] = (3_417_968_750. / 16_647_192_600.) * ecc_13;
        first[1] = (64. / 675.) * ecc_11;
        first[2] = (6561. / 163_840.) * ecc_9
            + (19683. / 327_680.) * ecc_11
            + (45927. / 2_097_152.) * ecc_13;
        first[3] = (1. / 72.) * ecc_7 + (7. / 288.) * ecc_9 + (49. / 4800.) * ecc_11;
        first[4] = (1. / 384.) * ecc_5
            + (11. / 2304.) * ecc_7
            + (619. / 98304.) * ecc_9
            + (3443. / 983_040.) * ecc_11
            + (685_783. / 471_859_200.) * ecc_13;
        first[5] = 0.;
        first[6] =
            0.5 * eccentricity - (1. / 4.) * ecc_3 + (13. / 128.) * ecc_5 + (113. / 2304.) * ecc_7;
        first[7] = -10.0 * eccentricity + (63. / 2.) * ecc_3
            - (155. / 6.) * ecc_5
            - (2921. / 288.) * ecc_7;

        first[8] = (49. / 2.) * eccentricity - (861. / 4.) * ecc_3 + (65925. / 128.) * ecc_5
            - (132_635. / 256.) * ecc_7;
        first[9] = 289. * ecc_3 - (5865. / 8.) * ecc_5 + (203_147. / 96.) * ecc_7;
        first[10] = (714_025. / 384.) * ecc_5 - (27_483_625. / 2304.) * ecc_7
            + (2_936_126_875. / 98304.) * ecc_9
            - (6_772_518_125. / 196_608.) * ecc_11
            + (303_503_554_375. / 18_874_368.) * ecc_13;
        first[11] = (284_089. / 32.) * ecc_7 - (66_328_119. / 1280.) * ecc_8
            + (191_185_929. / 2560.) * ecc_9;
        first[12] = (52_142_352_409. / 1_474_560.) * ecc_9 - (140_254_152_605. / 589_824.) * ecc_11
            + (66_020_511_589_375. / 169_869_312.) * ecc_13;
        first[13] = (5_383_010_161. / 43200.) * ecc_11;
        first[14] = (147_483_366_698_529. / 367_001_600.) * ecc_13;

        second[0] = (1_306_918_425_967. / 471_859_200.) * ecc_13;
        second[1] = (30_089_667. / 25600.) * ecc_11;
        second[2] = (15_717_645. / 32768.) * ecc_9 - (8_841_951. / 65536.) * ecc_11
            + (174_091_183. / 18_874_368.) * ecc_13;
        second[3] = (5929. / 32.) * ecc_7 + (9933. / 128.) * ecc_9;
        second[4] = (8427. / 128.) * ecc_5
            + (20829. / 256.) * ecc_7
            + (6_019_881. / 32768.) * ecc_9
            + (29_183_787. / 327_680.) * ecc_11
            + (4_288_977_063. / 52_428_800.) * ecc_13;
        second[5] = (81. / 4.) * ecc_3 + (189. / 4.) * ecc_5 + (1661. / 16.) * ecc_7;
        second[6] = (9. / 2.) * eccentricity
            + (81. / 4.) * ecc_3
            + (6885. / 128.) * ecc_5
            + (12_123_879. / 20608.) * ecc_7;

        second[7] = 6. * eccentricity * (1. - ecc_2).powi(-4);
        // array reflected here
        let (left, right) = second.split_at_mut(8);
        right.copy_from_slice(&left[0..=6]);
        right.reverse();

        // Copy in reverse order
        third.copy_from_slice(&first[..]);
        third.reverse();
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

        let borrowed_slice = &mut self.eccentricity_2pq_squared[..];
        let [first, second, third] = borrowed_slice else {
            unreachable!()
        };

        first[0] = (15625. / 129_024.) * ecc_7;
        first[1] = (4. / 45.) * ecc_6;
        first[2] = (81. / 1280.) * ecc_5 + (81. / 2048.) * ecc_7;
        first[3] = (1. / 24.) * ecc_4 + (7. / 240.) * ecc_6;
        first[4] = (1. / 48.) * ecc_3 + (11. / 768.) * ecc_5 + (313. / 30720.) * ecc_7;
        first[5] = 0.;
        first[6] = -0.5 * ecc + (1. / 16.) * ecc_3 - (5. / 384.) * ecc_5 - (143. / 18432.) * ecc_7;
        first[7] = 1. - (5. / 2.) * ecc_2 + (13. / 16.) * ecc_4 - (35. / 288.) * ecc_6;
        first[8] = (7. / 2.) * ecc - (123. / 16.) * ecc_3 + (489. / 128.) * ecc_5
            - (1763. / 2048.) * ecc_7;
        first[9] = (17. / 2.) * ecc_2 - (115. / 16.) * ecc_4 + (601. / 48.) * ecc_6;
        first[10] = (845. / 48.) * ecc_3 - (32525. / 768.) * ecc_5 + (208_225. / 6144.) * ecc_7;
        first[11] = (533. / 16.) * ecc_4 - (13827. / 160.) * ecc_6;
        first[12] = (228_347. / 3840.) * ecc_5 - (3_071_075. / 18432.) * ecc_7;
        first[13] = (73369. / 720.) * ecc_6;
        first[14] = (12_144_273. / 71680.) * ecc_7;

        for val in first.iter_mut() {
            *val = val.powi(2);
        }

        second[0] = (432_091. / 30720.) * ecc_7;
        second[1] = (3167. / 320.) * ecc_6;
        second[2] = (1773. / 256.) * ecc_5 + (4987. / 6144.) * ecc_7;
        second[3] = (77. / 16.) * ecc_4 + (129. / 160.) * ecc_6;
        second[4] = (53. / 16.) * ecc_3 + (393. / 256.) * ecc_5 + (24753. / 10240.) * ecc_7;
        second[5] = (9. / 4.) * ecc_2 + (7. / 4.) * ecc_4 + (141. / 64.) * ecc_6;
        second[6] = (3. / 2.) * ecc
            + (27. / 16.) * ecc_3
            + (261. / 128.) * ecc_5
            + (14309. / 6144.) * ecc_7;
        second[7] = (1. - ecc_2).powf(-3. / 2.);
        // array refelected here
        let (left, right) = second.split_at_mut(8);
        right.copy_from_slice(&left[0..=6]);
        right.reverse();

        for val in second.iter_mut() {
            *val = val.powi(2);
        }

        // Copy in reverse order
        third.copy_from_slice(&first[..]);
        third.reverse();
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Cayley 1861 https://ui.adsabs.harvard.edu/abs/1861MmRAS..29..191C
// Kaula 1964 10.1029/RG002i004p00661
