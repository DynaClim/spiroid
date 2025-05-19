use crate::constants::{SECONDS_IN_YEAR, SOLAR_LUMINOSITY, SOLAR_MASS, SOLAR_RADIUS};
use serde::{Deserialize, Serialize};

// Interpolation values deserialized from user provided CSV.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub struct StarCsv {
    age: f64,                          // (s)
    radius: f64,                       // (m)
    mass: f64,                         // (kg)
    convective_radius: f64,            // (m)
    convective_mass: f64,              // (kg)
    radiative_moment_of_inertia: f64,  // (kg.m2)
    convective_moment_of_inertia: f64, // (kg.m2)
    luminosity: f64,                   // (J.s-1)

    // Calculated internally, not included in the CSV.
    #[serde(default)]
    convective_moment_of_inertia_derivative: f64,
    #[serde(default)]
    radiative_mass_derivative: f64,
}

impl StarCsv {
    pub fn initialise(stars: &mut [Self]) -> (Vec<f64>, Vec<Vec<f64>>) {
        stars.iter_mut().for_each(Self::convert_units);
        Self::compute_derivatives(stars);

        // Split the values into a vector of ages and a nested vector of remaining values.
        // The ages are used as the index to interpolate remaining values, based on time.
        let ages = stars
            .iter()
            .map(|starcsv| starcsv.age)
            .collect::<Vec<f64>>();
        let rest = stars.iter().map(StarCsv::to_vec).collect::<Vec<Vec<f64>>>();

        (ages, rest)
    }

    // Initialise the input values with unit conversion.
    fn convert_units(&mut self) {
        self.age *= SECONDS_IN_YEAR;
        self.radius *= SOLAR_RADIUS;
        self.mass *= SOLAR_MASS;
        self.luminosity *= SOLAR_LUMINOSITY;
        self.convective_radius *= SOLAR_RADIUS;
        self.convective_mass *= SOLAR_MASS;
        self.radiative_moment_of_inertia *= self.mass * self.radius.powi(2);
        self.convective_moment_of_inertia *= self.mass * self.radius.powi(2);
    }

    // Calcultes the radiative_mass_derivative and convective_moment_of_inertia_derivative for each record.
    fn compute_derivatives(stars: &mut [Self]) {
        let stars_len = stars.len();

        // Derivative is zero for first and last timesteps.
        stars[0].radiative_mass_derivative = 0.;
        stars[0].convective_moment_of_inertia_derivative = 0.;
        stars[stars_len - 1].radiative_mass_derivative = 0.;
        stars[stars_len - 1].convective_moment_of_inertia_derivative = 0.;

        for i in 1..stars_len - 1 {
            // Unpack values of the star at three consecutive timesteps to compute the derivatives.
            let [prev, curr, next] = &mut stars[i - 1..=i + 1] else {
                unreachable!()
            };
            curr.radiative_mass_derivative =
                (next.convective_mass - prev.convective_mass) / (next.age - prev.age);
            curr.convective_moment_of_inertia_derivative = (next.convective_moment_of_inertia
                - prev.convective_moment_of_inertia)
                / (next.age - prev.age);
        }
    }

    fn to_vec(&self) -> Vec<f64> {
        vec![
            self.age,
            self.radius,
            self.mass,
            self.convective_radius,
            self.convective_mass,
            self.radiative_moment_of_inertia,
            self.convective_moment_of_inertia,
            self.luminosity,
            self.radiative_mass_derivative,
            self.convective_moment_of_inertia_derivative,
        ]
    }
}
