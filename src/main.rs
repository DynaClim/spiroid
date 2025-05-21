use anyhow::Result;
use rayon::prelude::*;
use sci_file::{read_csv_columns_from_file, read_csv_rows_from_dir, read_csv_rows_from_file};
use spiroid_lib::{ParticleType, Simulation, StarCsv, Universe};

fn main() -> Result<()> {
    let simulations = Simulation::<Universe>::new()?;
    simulations
        .into_par_iter()
        .map(|mut simulation| {
            let initial_time = simulation.initial_time;
            let final_time = simulation.final_time;

            if let ParticleType::Star(star) = &mut simulation.system.central_body.kind {
                // Load stellar evolution data from file if stellar evolution is enabled.
                if let Some(star_file) = star.evolution_file() {
                    // Maps every row of the csv file into a `StarCsv`.
                    let mut stellar_data = read_csv_rows_from_file::<StarCsv>(star_file)?;
                    // Configure the stellar evolution interpolator.
                    let (star_ages, star_values) = StarCsv::initialise(&mut stellar_data);
                    star.initialise_evolution(&star_ages, &star_values);
                }
            }

            // Load love number data from file(s) if kaula tides are enabled.
            if let Some(kaula) = simulation.system.orbiting_body.tides.kaula_get_mut() {
                if let Some(solid_file) = kaula.solid_file() {
                    // Maps each column of love number data into a vector.
                    let love_solid = read_csv_columns_from_file::<f64>(solid_file)?;
                    kaula.initialise_love_number_solid(&love_solid);
                }
                if let Some(ocean_file) = kaula.ocean_file() {
                    let love_ocean = read_csv_columns_from_file::<f64>(ocean_file)?;
                    kaula.initialise_love_number_ocean(&love_ocean);
                }
                if let Some(interpolate_dir) = kaula.interpolate_dir() {
                    let _interpolation_3d = read_csv_rows_from_dir::<f64>(interpolate_dir)?;
                    todo!();
                }
            }

            // Initialise the universe (star, planet, etc).
            simulation.system.initialise(initial_time)?;

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

            simulation.launch(initial_time, final_time, &y)?;
            Ok(())
        })
        .collect::<Result<()>>()?;
    Ok(())
}
