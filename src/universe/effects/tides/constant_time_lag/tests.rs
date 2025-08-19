use super::*;
use crate::universe::particles::planet::tests::test_planet;
use crate::universe::particles::star::tests::test_star;

use pretty_assertions::assert_eq;

pub fn test_constant_time_lag() -> ConstantTimeLag {
    ConstantTimeLag {
        equilibrium: Equilibrium::SigmaBarStar(1e-6),
        inertial: Inertial::FrequencyAveraged,
    }
}

#[test]
// This function is only called if tides are enabled.
fn _tidal_torque() {
    let expected = 6.325284391272144e23;
    let mut star = test_star();
    let planet = test_planet();
    let tidal_model = test_constant_time_lag();
    star.refresh_tidal_frequency(&planet);
    let result = tidal_model.tidal_torque(&star, &planet);
    assert_eq!(expected, result);
}
