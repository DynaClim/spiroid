use anyhow::Result;
use rayon::prelude::*;
use riptide::Layer;
use sci_file::{DataStore, read_csv_rows_from_file, read_json_from_file};
use spiroid_lib::{ParticleType, SECONDS_IN_YEAR, Simulation, StarCsv, Universe, UniverseIntegral};
fn main() -> Result<()> {
    // Parse the command line arguments into one (or more) simulations.
    // Universe is the struct containig the state of the universe and corresponding derivation functions.
    // UniverseIntegral contains the quantities that are to be integrated.
    let simulations = Simulation::<Universe, UniverseIntegral>::new()?;

    // Read in necessary data tables and launch the simulations in parallel.
    simulations
        .into_par_iter()
        .map(|mut simulation| {
            if let ParticleType::Star(star) = &mut simulation.system.central_body.kind {
                // Load stellar evolution data from file if stellar evolution is enabled.
                if let Some(star_file) = star.evolution_file() {
                    // Maps every row of the csv file into a `StarCsv`.
                    let mut stellar_data = read_csv_rows_from_file::<StarCsv>(star_file)?;
                    // Configure the stellar evolution interpolator.
                    StarCsv::initialise(&mut stellar_data);
                    let star_ages = StarCsv::ages(&stellar_data);
                    star.initialise_evolution(&star_ages, &stellar_data)?;
                }
            }

            // Load love number data from file(s) if kaula tides are enabled.
            if let Some(kaula) = simulation.system.orbiting_body.tides.kaula_get_mut() {
                if kaula.interpolation_mode() {
                    // Solid planet
                    if let Some((solid_file, solid_k2)) = kaula.solid_file() {
                        *solid_k2 = read_json_from_file(solid_file)?;
                        if let DataStore::Interpolate2D(interpolator_2d) = solid_k2 {
                            // Transpose the interpolator from (x:years, y:tidal_frequency) to (x:tidal_frequency, y:years)
                            interpolator_2d.transpose();
                            // Convert the time dimension of the 2D interpolation from years to seconds.
                            interpolator_2d
                                .y_vals_mut()
                                .iter_mut()
                                .for_each(|year| *year *= SECONDS_IN_YEAR);
                        }
                        solid_k2.dimension_check()?;
                    }
                    // Liquid planet
                    if let Some((liquid_file, liquid_k2)) = kaula.liquid_file() {
                        *liquid_k2 = read_json_from_file(liquid_file)?;
                        if let DataStore::Interpolate2D(interpolator_2d) = liquid_k2 {
                            // Transpose the interpolator from (x:years, y:tidal_frequency) to (x:tidal_frequency, y:years)
                            interpolator_2d.transpose();
                            // Convert the time dimension of the 2D interpolation from years to seconds.
                            interpolator_2d
                                .y_vals_mut()
                                .iter_mut()
                                .for_each(|year| *year *= SECONDS_IN_YEAR);
                        }
                        liquid_k2.dimension_check()?;
                    }
                } else {
                    if let Some(data_file) = kaula.internal_structure_file() {
                        let layers = read_csv_rows_from_file::<Layer>(data_file)?;
                        kaula.initialise_internal_structure(&layers);
                    }
                }
                if let ParticleType::Star(star) = &simulation.system.central_body.kind
                    && let ParticleType::Planet(planet) = &simulation.system.orbiting_body.kind
                {
                    kaula.initialise_cache(simulation.initial_time, star, planet)?;
                }
            }

            // Initialise the universe (star, planet, etc).
            simulation.system.initialise(simulation.initial_time)?;

            // Initialise the values to integrate.
            let y = simulation.system.integration_quantities();
            // y[0] = Star radiative zone angular momentum
            // y[1] = Star convective zone angular momentum
            // y[2] = Planet semi-major axis^6.5

            // Only if kaula tides are enabled on the planet:
            // y[3] = Planet spin
            // y[4] = Planet orbital eccentricity^2
            // y[5] = Planet orbital inclination (with respect to the planet equatorial plane)
            // y[6] = Planet longitude of ascending node
            // y[7] = Planet argument of periapsis
            // y[8] = Planet spin axis inclination (with respect to the total angular momentum)

            simulation.launch(simulation.initial_time, simulation.final_time, &[y])?;
            Ok(())
        })
        .collect::<Result<()>>()?;
    Ok(())
}
