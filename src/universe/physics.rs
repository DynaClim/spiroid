use crate::constants::GRAVITATIONAL;
use crate::universe::effects::tides::TidalModel;
use crate::universe::{Kaula, ParticleType, Planet, Star, Universe};
use anyhow::{Result, bail};

pub(crate) fn force(dy: &mut [f64], universe: &mut Universe) -> Result<()> {
    dy.fill(0.0);

    let ParticleType::Star(star) = &mut universe.central_body.kind else {
        todo!();
    };

    let ParticleType::Planet(planet) = &mut universe.orbiting_body.kind else {
        todo!();
    };

    // Star derivatives
    dy[0] = star_radiative_zone_angular_momentum_derivative(star);
    dy[1] = star_convective_zone_angular_momentum_derivative(star, universe.disk_is_dissipated);

    // If the planet is destroyed, only the star derivatives are computed.
    if planet.is_destroyed {
        return Ok(());
    }

    // Constant time lag semi major axis derivative.
    // Is 0 if tides are disabled on the star.
    dy[2] = planet_semi_major_axis_13_div_2_derivative(planet, star);

    // Immutable borrow of kaula properties if kaula planet tides enabled.
    if let TidalModel::KaulaTides { ref kaula } = universe.orbiting_body.tides {
        // Sum the semi major axis derivative to account for both CTL star tide (if enabled) and Kaula planet tide.
        dy[2] += kaula_planet_semi_major_axis_13_div_2_derivative(planet, star, kaula);

        dy[3] = planet_spin_derivative(planet, star, kaula);
        dy[4] = planet_eccentricity_derivative(planet, star, kaula);
        dy[5] = planet_inclination_derivative(planet, star, kaula);
        dy[6] = planet_longitude_ascending_node_derivative(planet, star, kaula);
        dy[7] = planet_argument_pericentre_derivative(planet, star, kaula);
        dy[8] = planet_spin_axis_inclination_derivative(planet, star, kaula);
    }

    for val in dy.iter() {
        if *val != 0.0 && !val.is_normal() {
            let msg = format!("{:?}, {:?}, dy: {:?}", &star, &planet, &dy);
            eprintln!("{}", &msg);
            bail!(
                "error in computation of derivatives: Houston, we have a NaN...infinity, and beyond! {msg}"
            );
        }
    }

    Ok(())
}

// Rate of change in the angular momentum in the convective zone.
// Ahuir et al. 2021, Eq. 2.
fn star_convective_zone_angular_momentum_derivative(star: &Star, disk_is_dissipated: bool) -> f64 {
    // If the disk has not dissipated, the spin of the star has not evolved.
    if !disk_is_dissipated {
        return star.spin * star.convective_moment_of_inertia_derivative;
    }
    star.angular_momentum_redistribution / star.core_envelope_coupling_constant
        - star.mass_transfer_envelope_to_core_torque
        + star.wind_torque
        + star.magnetic_torque
        + star.tidal_torque
}

// Rate of change in the angular momentum in the radiative zone.
// Benbakoura et al. 2019, Eq. 2
fn star_radiative_zone_angular_momentum_derivative(star: &Star) -> f64 {
    -star.angular_momentum_redistribution / star.core_envelope_coupling_constant
        + star.mass_transfer_envelope_to_core_torque
}

// This loosely comes from Eq. 1 from Ahuir et al. 2021 (for the sum of tidal and magnetic components)
// Also Benbakoura et al. 2019, Eq. 3.
// It is the derivative of semi major axis (a) to the power 13/2
// This is obtained by moving the 1/a^6 dependency of the tidal torque to the left of Eq. 3, alongside the a^(1/2)
// this means that what we call here the tidal torque is not exactly the tidal torque, but the tidal torque * a^6
// or tidal torque without the semi-major axis dependency
pub(crate) fn planet_semi_major_axis_13_div_2_derivative(planet: &Planet, star: &Star) -> f64 {
    -13. * sqrt!((star.mass + planet.mass) / GRAVITATIONAL)
        * (1. / (star.mass * planet.mass))
        * planet.semi_major_axis.powi(6)
        * (star.magnetic_torque + star.tidal_torque)
}

// Semi-major axis derivative.
// Boue & Efroimksy (2019) Eq. 116 and Revol et al. (2023) Eq A.1
pub(crate) fn kaula_planet_semi_major_axis_13_div_2_derivative(
    planet: &Planet,
    star: &Star,
    kaula: &Kaula,
) -> f64 {
    -13. * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
        * (star.mass / planet.mass)
        * planet.radius.powi(5)
        * kaula.summation_of_longitudinal_modes_semi_major_axis()
}

// Spin derivative.
// Boue & Efroimksy (2019) Eq. 123 and Revol et al. (2023) Eq A.3
fn planet_spin_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    let planet_tidal_torque = (GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5))
        / (planet.semi_major_axis.powi(6));

    (planet_tidal_torque / planet.moment_of_inertia) * kaula.summation_of_longitudinal_modes_spin()
}

// Eccentricity derivative.
// Boue & Efroimksy (2019) Eq. 117 and Revol et al. (2023) Eq A.3
fn planet_eccentricity_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if planet.eccentricity == 0. {
        0.
    } else {
        -2.0 * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
            * (planet.radius.powi(5) / planet.semi_major_axis.powf(6.5))
            * (star.mass / planet.mass)
            * sqrt!(1.0 - planet.eccentricity.powi(2))
            * kaula.summation_of_longitudinal_modes_eccentricity(planet.eccentricity)
    }
}

// Inclination derivative.
// Boue & Efroimksy (2019) Eq. 118 and Revol et al. (2023) Eq A.7
// The inclination refers to the angle between the orbital planet and the planet's equatorial plane (i.e. obliquity)
fn planet_inclination_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if sin!(planet.inclination) == 0.0 {
        0.0
    } else {
        let beta = (star.mass * planet.mass) / (star.mass + planet.mass);
        let term1 = (beta * planet.mean_motion.powi(2) * planet.semi_major_axis.powi(2))
            / (planet.moment_of_inertia * planet.spin);
        let term3 = planet.mean_motion / sqrt!(1. - planet.eccentricity.powi(2));

        (1. / sin!(planet.inclination))
            * (star.mass / planet.mass)
            * (planet.radius / planet.semi_major_axis).powi(5)
            * kaula.summation_of_longitudinal_modes_inclination(planet.inclination, term1, term3)
    }
}

// Longitude of ascending node derivative.
// Boue & Efroimksy (2019) Eq. 121 and Revol et al. (2023) Eq A.9
fn planet_longitude_ascending_node_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if (planet.inclination == 0.) || (planet.spin_inclination == 0.) {
        0.0
    } else {
        let beta = (star.mass * planet.mass) / (star.mass + planet.mass);
        let term1 = (1. / (planet.moment_of_inertia * planet.spin * tan!(planet.inclination)))
            - (cos!(planet.longitude_ascending_node)
                / (planet.moment_of_inertia * planet.spin * tan!(planet.spin_inclination)))
            + (1.
                / (beta
                    * planet.mean_motion
                    * planet.semi_major_axis.powi(2)
                    * sqrt!(1. - planet.eccentricity.powi(2))
                    * sin!(planet.inclination)));

        let term2 = -(sin!(planet.longitude_ascending_node) * cotan!(planet.inclination))
            / (planet.moment_of_inertia * planet.spin * tan!(planet.spin_inclination));
        let term3 = sin!(planet.longitude_ascending_node)
            / (planet.moment_of_inertia
                * planet.spin
                * tan!(planet.spin_inclination)
                * sin!(planet.inclination));

        ((GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5))
            / planet.semi_major_axis.powi(6))
            * kaula.summation_of_longitudinal_modes_longitude_ascending_node(term1, term2, term3)
    }
}

// Longitude of pericentre derivative.
// Boue & Efroimksy (2019) Eq. 120 and Revol et al. (2023) Eq A.11
fn planet_argument_pericentre_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if (planet.eccentricity == 0.)
        && ((planet.inclination == 0.) || (planet.spin_inclination == 0.))
    {
        // no eccentricity and no inclination.
        return 0.0;
    }

    let beta = (star.mass * planet.mass) / (star.mass + planet.mass);

    // inclination
    let summation_of_longitudinal_modes_pericentre_inclination = if planet.inclination == 0. {
        0.
    } else {
        let term1 = -((1. / (planet.moment_of_inertia * planet.spin * sin!(planet.inclination)))
            + (1.
                / (planet.mean_motion
                    * planet.semi_major_axis.powi(2)
                    * sqrt!(1. - planet.eccentricity.powi(2))
                    * tan!(planet.inclination)
                    * beta)));
        kaula.summation_of_longitudinal_modes_pericentre_inclination(
            planet.inclination,
            planet.spin_inclination,
        ) * term1
    };

    // eccentricity
    let summation_of_longitudinal_modes_pericentre_eccentricity = if planet.eccentricity == 0. {
        0.
    } else {
        let term2 = sqrt!(1. - planet.eccentricity.powi(2))
            / (planet.mean_motion * planet.semi_major_axis.powi(2) * planet.eccentricity * beta);
        kaula.summation_of_longitudinal_modes_pericentre_eccentricity(planet.eccentricity) * term2
    };

    ((GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5)) / planet.semi_major_axis.powi(6))
        * (summation_of_longitudinal_modes_pericentre_eccentricity
            + summation_of_longitudinal_modes_pericentre_inclination)
}

// Spin axis inclination derivative.
// Boue & Efroimksy (2019) Eq 122 and Revol et al. (2023) Eq A.13
// The spin axis inclination refers to the inclination of the planet's rotational vector
// with respect to the total angular momentum.
fn planet_spin_axis_inclination_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if (planet.inclination == 0.0) || (planet.spin_inclination == 0.0) {
        0.
    } else {
        (GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5))
            / (planet.semi_major_axis.powi(6) * planet.moment_of_inertia * planet.spin)
            * kaula.summation_of_longitudinal_modes_spin_axis_inclination(
                planet.longitude_ascending_node,
                planet.inclination,
            )
    }
}

#[cfg(test)]
mod tests;

// References:
// Ahuir et al. 2021, https://doi.org/10.1051/0004-6361/202040173
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Boué and Efroimsky 2019, https://doi.org/10.1007/s10569-019-9908-2
// Revol et al. 2023, https://doi.org/10.1051/0004-6361/202245790
