use crate::constants::{
    BOLTZMANN_CONST, GRAVITATIONAL, MAGNETIC_PERMEABILITY_OF_VACUUM, PI, PROTON_MASS,
    ROSSBY_SATURATION, ROSSBY_SUN, SECONDS_IN_DAY, SOLAR_CORONA_DENSITY, SOLAR_CORONA_TEMPERATURE,
    SOLAR_MASS, SOLAR_SURFACE_MAGNETIC_FIELD, TWO_PI,
};
use crate::universe::particles::{Planet, Star, magnetic_pressure};

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub enum MagneticModel {
    #[default]
    Disabled,
    Wind(IsothermalWind),
}

impl MagneticModel {
    // Compute the magnetic torque if magnetism is enabled and the planet is inside the alfven radius.
    pub(crate) fn magnetic_torque(&mut self, planet: &Planet, star: &Star) -> f64 {
        match self {
            MagneticModel::Disabled => 0.0,
            MagneticModel::Wind(wind) => {
                // Avoid performing magnetic torque calculation if:
                // - the planet is destroyed, or
                // - the planet has no magnetic field, or
                // - the planet is too far away from the star
                if planet.is_destroyed()
                    || planet.magnetic_pressure == 0.0
                    || !planet.inside_alfven_radius(star.alfven_radius)
                {
                    0.0
                } else {
                    // Computation of the useful wind parameters.
                    wind.init_weber_davis(planet, star);

                    wind.magnetic_torque
                }
            }
        }
    }
}

// Output flag indicating the nature of the magnetic interaction
#[derive(Debug, PartialEq, Default, Serialize, Deserialize, Clone)]
pub enum MagneticInteraction {
    #[default]
    None, // Only exists for Default trait.
    Unipolar,
    Dipolar,
}

#[derive(Debug, PartialEq, Default, Serialize, Deserialize, Clone)]
#[serde(default)]
// Computation of a 1D magnetized isothermal wind to quantify angular momentum exchange in Star-Planet Magnetic Interactions
pub struct IsothermalWind {
    // User specified
    pub(crate) footpoint_conductance: f64, // (Ohm-1)

    // Calculated internally
    speed_of_sound: f64,                    // (m.s-1)
    critical_radius: f64,                   // (m)
    critical_radius_div_alfven_radius: f64, // ()
    magnetic_torque: f64,                   // (erg)
    radial_magnetic_field: f64,             // (T)
    magnetic_pressure: f64,                 // (erg.cm-3)
    integration_constant: f64,              // ()
    wind_velocity: f64,                     // (alfven speed at alfven radius)
    surface_wind_velocity: f64,             // (alfven speed at alfven radius)
    wind_density: f64,                      // (g.cm-3)

    alfvenic_mach: f64,                 // ()
    azimuthal_velocity: f64,            // (cm.s-1)
    alfven_speed_at_alfven_radius: f64, // (cm.s-1)

    interaction: MagneticInteraction, // Nature of the magnetic interaction (unipolar or dipolar)
}

impl IsothermalWind {
    // ***WARNING!***
    // Stateful function.
    // The order of these calculations is important.
    // Lower order calculations depend on previous values.
    // ***WARNING!***
    // Calculates the characteristics of the stellar wind at a given distance from the star.
    // The characteristics are computed following the magnetized model of Weber & Davis (1967)
    fn init_weber_davis(&mut self, planet: &Planet, star: &Star) {
        let surface_magnetic_field = Self::magnetic_field(star.mass, star.rossby);
        self.radial_magnetic_field = Self::radial_magnetic_field(
            surface_magnetic_field,
            star.radius,
            planet.semi_major_axis,
        );

        let coronal_temperature = Self::coronal_temperature(star.mass, star.rossby);
        self.speed_of_sound = Self::speed_of_sound(coronal_temperature);

        self.critical_radius = GRAVITATIONAL * star.mass / (2. * self.speed_of_sound.powi(2));
        self.critical_radius_div_alfven_radius = self.critical_radius / star.alfven_radius;

        self.alfven_speed_at_alfven_radius = self.alfven_speed_at_alfven_radius(star);

        self.integration_constant = self.integration_constant(star); // Requires alfven_speed_at_alfven_radius
        self.magnetic_pressure = magnetic_pressure(self.radial_magnetic_field); // Requires radial_magnetic_field
        self.wind_velocity = self.weber_davis_velocity_profile(planet.semi_major_axis, star); // Requires alfven_speed_at_alfven_radius
        self.surface_wind_velocity = self.weber_davis_velocity_profile(star.radius, star); // Requires alfven_speed_at_alfven_radius

        let coronal_density = Self::coronal_density(star.mass, star.rossby);
        self.wind_density =
            self.density_profile(star.radius, coronal_density, planet.semi_major_axis);

        self.azimuthal_velocity = star.spin * planet.semi_major_axis;
        let keplerian_velocity = sqrt!(GRAVITATIONAL * star.mass / planet.semi_major_axis);
        self.alfvenic_mach = self.alfvenic_mach(keplerian_velocity);
        self.magnetic_torque = self.magnetic_torque(planet, star);
    }

    // Estimate the stellar surface magnetic field based on scaling laws.
    // Ahuir et al. 2020, Eq. 67
    fn magnetic_field(star_mass: f64, star_rossby: f64) -> f64 {
        SOLAR_SURFACE_MAGNETIC_FIELD
            * (max!(star_rossby, ROSSBY_SATURATION) / ROSSBY_SUN).powi(-1)
            * (star_mass / SOLAR_MASS).powf(-1.76)
    }

    // Estimate the stellar coronal density based on scaling laws.
    // Ahuir et al. 2020, Eq. 66
    fn coronal_density(star_mass: f64, star_rossby: f64) -> f64 {
        let max_rossby = max!(star_rossby, ROSSBY_SATURATION);

        SOLAR_CORONA_DENSITY
            * (ROSSBY_SUN / max_rossby).powf(1.07)
            * (star_mass / SOLAR_MASS).powf(1.97)
    }

    // Estimate the stellar coronal temperature based on scaling laws.
    // Ahuir et al. 2020, Eq. 65
    fn coronal_temperature(star_mass: f64, star_rossby: f64) -> f64 {
        let max_rossby = max!(star_rossby, ROSSBY_SATURATION);

        SOLAR_CORONA_TEMPERATURE
            * (ROSSBY_SUN / max_rossby).powf(0.11)
            * (star_mass / SOLAR_MASS).powf(0.12)
    }

    // Speed of sound at the base of the corona, based on coronal temperature.
    fn speed_of_sound(coronal_temperature: f64) -> f64 {
        sqrt!(2. * BOLTZMANN_CONST * coronal_temperature / PROTON_MASS)
    }

    // Magnetic field strength at a given distance from the star.
    // Assumes an open magnetic field structure, with amplitude decreases like 1/r^2 to keep Div(B)=0.
    fn radial_magnetic_field(
        surface_magnetic_field: f64,
        star_radius: f64,
        distance_from_stellar_center: f64,
    ) -> f64 {
        // Open-field configuration
        surface_magnetic_field * (star_radius / distance_from_stellar_center).powi(2)
    }

    // Compute density at a distance from the star.
    // The conservation of mass is used, i.e. mp*n(Rstar)*v(Rstar)*Rstar^2 = mp*n(r)*v(r)*r^2.
    fn density_profile(
        &self,
        star_radius: f64,
        star_coronal_density: f64,
        distance_from_stellar_center: f64,
    ) -> f64 {
        PROTON_MASS * star_coronal_density * self.surface_wind_velocity * star_radius.powi(2)
            / (self.wind_velocity * distance_from_stellar_center.powi(2))
    }

    // Compute Alfvénic Mach number at a given distance, which corresponds to a given keplerian velocity.
    // Ma = v/va with v = sqrt(vr^2+vphi^2) and va = B/sqrt(mu0 * rho).
    // The radial component of the wind speed, vr, has to be multiplied buy the alfven speed at the alfven radius for units purposes here.
    // The relative longitudinal velocity is v_keplerian - v_phi.
    fn alfvenic_mach(&self, keplerian_velocity: f64) -> f64 {
        // Open-field configuration
        sqrt!(
            (self.wind_velocity * self.alfven_speed_at_alfven_radius).powi(2)
                + (keplerian_velocity - self.azimuthal_velocity).powi(2)
        ) * sqrt!(MAGNETIC_PERMEABILITY_OF_VACUUM * self.wind_density)
            / self.radial_magnetic_field
    }

    // Integration constant when solving for the Weber-Davis wind solution.
    // It corresponds to the energy flux computed at the Alfven radius.
    fn integration_constant(&self, star: &Star) -> f64 {
        // Integration constant used in the implicit equation which describes the WD profile
        0.5 - (GRAVITATIONAL * star.mass
            / (star.alfven_radius * self.alfven_speed_at_alfven_radius.powi(2)))
            + (star.spin.powi(2) * star.alfven_radius.powi(2)
                / (2. * self.alfven_speed_at_alfven_radius.powi(2)))
    }

    // Get wind speed at given distance from star.
    fn velocity_profile(&self, distance_from_stellar_center: f64, star: &Star) -> f64 {
        // Intial guess of the speed value at distance from star, to start the Newton-Raphson solver.
        let mut velocity = if (distance_from_stellar_center < 1.)
            && (distance_from_stellar_center >= self.critical_radius_div_alfven_radius)
        {
            // Supersonic and subalfvenic regime
            1.
        } else if distance_from_stellar_center >= 1. {
            // Superalfvenic regime
            1.5
        } else {
            // Subsonic regime
            1e-7
        };

        // energy_flux_difference corresponds to total energy flux F (Eq. 24 in Weber & Davis 1967) for initial guess 'velocity', minus the
        // input energy flux (integration constant).
        // F should be equal to 0 when the wind solution is found.
        let mut energy_flux_difference = self.total_energy_flux_minus_constant(
            self.alfven_speed_at_alfven_radius,
            velocity,
            distance_from_stellar_center,
            self.integration_constant,
            star,
        );
        while abs!(energy_flux_difference) >= 1e-7 {
            // Computation of the derivative of F for the down-gradient Newton-Raphson method.
            let d_energy_flux_by_dv =
                self.derivative_d_energy_flux_by_dv(velocity, distance_from_stellar_center, star);
            // Estimate of velocity from the previous Newton-Raphson step.
            velocity = Self::velocity_midpoint(
                velocity,
                velocity - energy_flux_difference / d_energy_flux_by_dv,
            );

            energy_flux_difference = self.total_energy_flux_minus_constant(
                self.alfven_speed_at_alfven_radius,
                velocity,
                distance_from_stellar_center,
                self.integration_constant,
                star,
            );
        }

        velocity
    }

    // New value of alfven_speed when the Alfven speed is negative.
    fn velocity_midpoint(prev_velocity: f64, mut velocity: f64) -> f64 {
        while velocity <= 0. {
            velocity = f64::midpoint(prev_velocity, velocity);
        }
        velocity
    }

    // Total energy flux F (Eq. 24 in Weber & Davis 1967) minus the integration constant.
    // When the WD solution is met, this function should return 0.
    fn total_energy_flux_minus_constant(
        &self,
        alfven_speed: f64,
        velocity: f64,
        radius: f64,
        integration_constant: f64,
        star: &Star,
    ) -> f64 {
        // Function F to cancel, as a function of the Alfven radius (F(speed_of_sound / alfven_speed_at_alfven_radius, critical_radius / alfven_radius) = 0 by knowing alfven_speed_at_alfven_radius as a function of integration_constant from F(1, 1) = 0)
        0.5 * velocity.powi(2)
            - ln!(velocity) * (self.speed_of_sound / alfven_speed).powi(2)
            - 2. * ln!(radius) * (self.speed_of_sound / alfven_speed).powi(2)
            - GRAVITATIONAL * star.mass / (star.alfven_radius * radius * alfven_speed.powi(2))
            + (star.spin.powi(2) * star.alfven_radius.powi(2)
                / (2. * radius.powi(2) * alfven_speed.powi(2)))
                * (1.
                    + ((-1. + 2. * velocity * radius.powi(2)) * (1. - radius.powi(2)).powi(2)
                        / (-1. + velocity * radius.powi(2)).powi(2)))
            - integration_constant
    }

    // Derivate of total_energy_flux_minus_constant with respect to velocity, used for the down-gradient Newton-Raphson method.
    fn derivative_d_energy_flux_by_dv(&self, velocity: f64, radius: f64, star: &Star) -> f64 {
        velocity
            - (1. / velocity) * (self.speed_of_sound / self.alfven_speed_at_alfven_radius).powi(2)
            + (star.spin.powi(2) * star.alfven_radius.powi(2)
                / (radius * self.alfven_speed_at_alfven_radius.powi(2)))
                * (1. - radius.powi(2)).powi(2)
                * radius.powi(2)
                * (1. / (velocity * radius.powi(2) - 1.).powi(2)
                    - (2. * velocity * radius.powi(2) - 1.)
                        / (velocity * radius.powi(2) - 1.).powi(3))
    }

    // Test if a solution is looked for too close to the sonic point, that diverges using the method coded here.
    // A safeguard of 1.e-3 Rstar is put into place around the sonic point.
    fn near_sonic_point(
        &self,
        width_around_sonic_point: f64,
        normalised_radius: f64,
        star: &Star,
    ) -> bool {
        //1e-3 is emprically determined to transit from unipolar to dipolar smoothly.
        (self.critical_radius_div_alfven_radius + 1e-3 - width_around_sonic_point
            > star.radius / star.alfven_radius)
            && (abs!(
                (2. * (self.critical_radius_div_alfven_radius + 1e-3) - width_around_sonic_point)
                    - 2. * normalised_radius
            ) <= width_around_sonic_point)
    }

    // Computation of the derivative, by calculating the derivative of F * alfven_speed_at_alfven_radius^2
    fn estimate_d_energy_flux_by_dv(&self, alfven_speed: f64, star: &Star) -> f64 {
        (self.speed_of_sound.powi(2) / alfven_speed) - alfven_speed
            + (star.spin.powi(2)
                * star.alfven_radius.powi(2)
                * (1. - self.critical_radius_div_alfven_radius.powi(2)).powi(2)
                / (self.critical_radius_div_alfven_radius.powi(2)))
                * (1.
                    / (self.speed_of_sound * self.critical_radius_div_alfven_radius.powi(2)
                        - alfven_speed)
                    + alfven_speed
                        * (2.
                            * self.speed_of_sound
                            * self.critical_radius_div_alfven_radius.powi(2)
                            - alfven_speed)
                        / (self.speed_of_sound * self.critical_radius_div_alfven_radius.powi(2)
                            - alfven_speed)
                            .powi(3))
    }

    // Value of the alfven speed at the alfven radius.
    // Used for computing the Weber-Davis solution at any other radius.
    fn alfven_speed_at_alfven_radius(&self, star: &Star) -> f64 {
        // Uses a Newton-Raphson method based on the fact that the wind profile will pass through the sonic and the alfvenic points.
        // Initial guess for the integration constant.
        let mut integration_constant = {
            if (-2. * GRAVITATIONAL * star.mass / star.alfven_radius)
                + star.spin.powi(2) * star.alfven_radius.powi(2)
                < 0.
            {
                0.445
            } else {
                0.6
            }
        };

        // Initial guess of the Alfven speed ensuring a super-alfvenic wind,  F(1, 1) = 0.
        let mut alfven_speed = sqrt!(
            ((-2. * GRAVITATIONAL * star.mass / star.alfven_radius)
                + star.spin.powi(2) * star.alfven_radius.powi(2))
                / (2. * integration_constant - 1.)
        );

        // Function F to cancel, as a function of the Alfven radius (F(speed_of_sound / alfven_speed_at_alfven_radius, critical_radius / alfven_radius) = 0 by knowing alfven_speed_at_alfven_radius as a function of integration_constant from F(1, 1) = 0)
        let mut energy_flux_difference = self.total_energy_flux_minus_constant(
            alfven_speed,
            self.speed_of_sound / alfven_speed,
            self.critical_radius_div_alfven_radius,
            integration_constant,
            star,
        );

        // Newton-Raphson method to find the Alfven speed.
        while abs!(energy_flux_difference) >= 1e-7 {
            let estimate_d_energy_flux_by_dv =
                self.estimate_d_energy_flux_by_dv(alfven_speed, star);
            let d_energy_flux_by_dv = estimate_d_energy_flux_by_dv / alfven_speed.powi(2)
                - 2. * energy_flux_difference / alfven_speed;

            alfven_speed = Self::velocity_midpoint(
                alfven_speed,
                alfven_speed - energy_flux_difference / d_energy_flux_by_dv,
            );

            integration_constant = 0.5
                - (GRAVITATIONAL * star.mass / (star.alfven_radius * alfven_speed.powi(2)))
                + (star.spin.powi(2) * star.alfven_radius.powi(2) / (2. * alfven_speed.powi(2)));

            // Update of the function to cancel.
            energy_flux_difference = self.total_energy_flux_minus_constant(
                alfven_speed,
                self.speed_of_sound / alfven_speed, // Update of the ratio speed_of_sound / alfven_speed_at_alfven_radius.
                self.critical_radius_div_alfven_radius,
                integration_constant,
                star,
            );
        }

        alfven_speed
    }

    // Computes a Weber-Davis velocity profile by using a linear interpolation near the critical points.
    fn weber_davis_velocity_profile(&self, radius: f64, star: &Star) -> f64 {
        let normalised_radius = radius / star.alfven_radius;

        // Width of the non-confidence area, where the Newton-Raphson method is not performed.
        let width_around_sonic_point = if TWO_PI / (star.spin * SECONDS_IN_DAY) < 1. {
            1.
        } else {
            7e-2
        };

        let width_around_alfven_radius = 1e-2;
        let x_prev;
        let x_next;
        let v_prev;
        let v_next;

        // The solution is not computed near the critical points due to non-crossing contours of the function F.
        // Newton-Raphson method
        if (abs!(
            (2. * (self.critical_radius_div_alfven_radius + 1e-3) - width_around_sonic_point)
                - 2. * normalised_radius
        ) > width_around_sonic_point)
            && (abs!(normalised_radius - 1.) > width_around_alfven_radius)
        {
            let velocity = self.velocity_profile(normalised_radius, star);
            #[allow(clippy::float_cmp)]
            if velocity != 10. {
                return velocity;
            }
        }
        if abs!(normalised_radius - 1.) <= width_around_alfven_radius {
            // Near the alfvenic point.
            x_prev = 1. - width_around_alfven_radius - 1e-2;
            x_next = 1. - width_around_alfven_radius;
            v_prev = self.velocity_profile(x_prev, star);
            v_next = self.velocity_profile(x_next, star);
        } else if self.near_sonic_point(width_around_sonic_point, normalised_radius, star) {
            // Near the sonic point.
            x_prev = self.critical_radius_div_alfven_radius + 1e-3 - width_around_sonic_point;
            x_next = self.critical_radius_div_alfven_radius + 1e-3;
            v_prev = self.velocity_profile(x_prev, star);
            v_next = self.velocity_profile(x_next, star);
        } else {
            x_prev = star.radius / star.alfven_radius;
            x_next = self.critical_radius_div_alfven_radius + 1e-3;
            v_next = self.velocity_profile(x_next, star);
            let radius_over_sonic_point =
                star.radius / (self.critical_radius_div_alfven_radius * star.alfven_radius);
            let mut velocity: f64 = if radius_over_sonic_point >= 1. {
                10.
            } else {
                1e-7
            };

            loop {
                let energy_flux_difference = velocity.powi(2)
                    - 2. * ln!(velocity)
                    - 4. / radius_over_sonic_point
                    - 4. * ln!(radius_over_sonic_point)
                    + 3.;
                if abs!(energy_flux_difference) < 1e-7 {
                    break;
                }

                let previous_velocity = velocity;
                velocity = velocity - 0.5 * energy_flux_difference / (velocity - 1. / velocity);

                velocity = Self::velocity_midpoint(previous_velocity, velocity);
            }

            v_prev = self.speed_of_sound * velocity / self.alfven_speed_at_alfven_radius;
        }

        // Linear interpolation of velocity, in units of alfven radius speed.
        ((v_next - v_prev) / (x_next - x_prev)) * (normalised_radius - x_prev) + v_prev
    }

    /// Computes the magnetic torque, taking into account unipolar and dipolar interaction between the planet and the star, following Strugarek et al. (2017)
    fn magnetic_torque(&mut self, planet: &Planet, star: &Star) -> f64 {
        // tidal_frequency difference omega_convective - omega_orbital

        // Ratio planet magnetic pressure / wind magnetic pressure, to check if a magnetosphere could be sustained
        // If lambda <= 1., Unipolar interaction, without magnetosphere.
        // Otherwise, Dipolar interaction, with magnetosphere.
        let lambda = planet.magnetic_pressure / self.magnetic_pressure;

        self.interaction = if lambda <= 1. {
            // Unipolar interaction,  without magnetosphere
            MagneticInteraction::Unipolar
        } else {
            // Dipolar interaction,  with magnetosphere
            MagneticInteraction::Dipolar
        };
        // Unipolar torque
        let torque_unipolar = 8.
            * planet.radius.powi(2)
            * planet.semi_major_axis.powi(2)
            * star.tidal_frequency
            * self.radial_magnetic_field.powi(2)
            * self.footpoint_conductance;

        // Dipolar torque
        let torque_dipolar = 10.8
            * PI
            * (self.alfvenic_mach / sqrt!(1. + self.alfvenic_mach.powi(2)))
            * self.magnetic_pressure
            * self.alfvenic_mach.powf(-0.56)
            * planet.radius.powi(2)
            * planet.semi_major_axis
            * lambda.powf(0.28);

        // Magnetic torque, interpolated between unipolar and dipolar torque with a tanh
        // to avoid discontinuities close to the dipolar/unipolar transisition (lambda = 1).
        // The values of alpha, depth1 and depth have been empirically tested.
        let alpha = 0.4;
        let depth1 = -3. * ln!(alpha);
        let depth = 1e-8;
        let sign = -tanh!(star.tidal_frequency / depth);

        let torque = abs!(
            0.5 * (torque_dipolar - torque_unipolar) * tanh!(10. * (ln!(lambda) + depth1) / depth1)
                + 0.5 * (torque_unipolar + torque_dipolar)
        );

        sign * torque
    }
}

#[cfg(test)]
mod tests;

// References:
// Ahuir et al. (2020), doi: https://doi.org/10.1051/0004-6361/201936974
// Weber & Davis (1967), doi: https://doi.org/10.1086/149138
// Strugarek et al. (2017), doi: https://doi.org/10.3847/2041-8213/aa8d70
