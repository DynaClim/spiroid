use super::*;
use pretty_assertions::assert_eq;
use sci_file::{
    OutputFile, deserialize_csv_column_vectors_from_path, deserialize_csv_rows_from_path,
    deserialize_json_from_path,
};
use simulation::InputConfig;
use simulation::{Integrator, System};
use std::path::{Path, PathBuf};

struct Test {
    pub data: Universe,
}

// Mock System implementation without writing any output.
impl System for Test {
    type Data = Universe;
    type Output = OutputFile;

    fn new(_: OutputFile, data: Universe) -> Self {
        Self { data }
    }

    fn derive(
        &mut self,
        time: f64,
        y: &[f64],
        dy: &mut [f64],
    ) -> Result<(), Box<(dyn std::error::Error + Send + Sync + 'static)>> {
        // Update the state of the universe based on the current integration values.
        self.data.update(time, y)?;
        // Compute the derivatives using the updated values.
        force(dy, &mut self.data)?;
        Ok(())
    }

    fn solout(
        &mut self,
        _time: f64,
        _y: &[f64],
    ) -> Result<(), Box<(dyn std::error::Error + Send + Sync + 'static)>> {
        Ok(())
    }
}

fn test_simulation(config: PathBuf) -> Universe {
    // Parse the config file.
    let mut config: InputConfig<Universe> = deserialize_json_from_path(&config).unwrap();

    // Load stellar evolution data from file.
    if let ParticleType::Star(star) = &mut config.universe.central_body.kind {
        // Load stellar evolution data from file if stellar evolution is enabled.
        if let Some(star_file) = star.evolution_file() {
            let mut stellar_data = deserialize_csv_rows_from_path::<StarCsv>(star_file).unwrap();
            // Configure the stellar evolution interpolator.
            let (star_ages, star_values) = StarCsv::initialise(&mut stellar_data);
            star.initialise_evolution(&star_ages, &star_values);
        };
    }

    // Load love number data from file(s) if kaula tides are enabled.
    if let Some(kaula) = config.universe.orbiting_body.tides.kaula_get_mut() {
        if let Some(solid_file) = kaula.solid_file() {
            let love_solid = deserialize_csv_column_vectors_from_path::<f64>(solid_file).unwrap();
            kaula.initialise_love_number_solid(&love_solid);
        }
        if let Some(ocean_file) = kaula.ocean_file() {
            let love_ocean = deserialize_csv_column_vectors_from_path::<f64>(ocean_file).unwrap();
            kaula.initialise_love_number_ocean(&love_ocean);
        }
    }

    // Initialise the universe (star, planet, etc).
    config.universe.initialise(config.initial_time).unwrap();

    // Initial values for the integrator.
    let y = config.universe.integration_quantities();
    let mut system = Test {
        data: config.universe,
    };
    config
        .integrator
        .initialise(config.initial_time, config.final_time, &y);

    // Run the full integration.
    let _ = config.integrator.integrate(&mut system).unwrap();

    system.data
}

fn compare_or_create(path: impl AsRef<Path> + std::fmt::Display, result: &Universe) {
    match deserialize_json_from_path::<Universe>(&path) {
        Ok(expected) => {
            // Saved file exists, compare the results.
            // We roundtrip our `Universe` through serde before comparison
            // to reset fields that are not serialized (serde skip_serializing)
            // (i.e. interpolation data read from file, internal buffers).
            let tmp = serde_json::to_string(&result).unwrap();
            let result: Universe = serde_json::from_str(&tmp).unwrap();
            assert_eq!(expected, result);
        }
        Err(err) => {
            match err {
                sci_file::Error::FileIo(_) => {
                    // Saved file does not exist save the results.
                    let mut writer = OutputFile::new(&path).unwrap();
                    writer.write_json(&result).unwrap();
                    panic!("comparison file `{path}` did not exist, so it was created");
                }
                _ => {
                    dbg!(&err);
                    panic!(
                        "the comparison file `{path}` is corrupt or has invalid structure. if it contains 'null' values, the value was probably NaN or inifinity"
                    );
                }
            }
        }
    }
}

#[test]
fn example_no_effects() {
    let result = test_simulation("examples/no_effects.conf".into());
    compare_or_create("examples/no_effects.expected", &result);
}

#[test]
fn example_tides() {
    let result = test_simulation("examples/tides.conf".into());
    compare_or_create("examples/tides.expected", &result);
}

#[test]
fn example_magnetic() {
    let result = test_simulation("examples/magnetic.conf".into());
    compare_or_create("examples/magnetic.expected", &result);
}

#[test]
fn example_magnetic_tides() {
    let result = test_simulation("examples/magnetic_tides.conf".into());
    compare_or_create("examples/magnetic_tides.expected", &result);
}

#[test]
fn example_kaula_solid() {
    //    let result = test_simulation("examples/kaula_solid.conf".into());
    //    compare_or_create("examples/kaula_solid.expected", &result);
}

#[test]
fn example_all_effects() {
    let result = test_simulation("examples/all_effects.conf".into());
    compare_or_create("examples/all_effects.expected", &result);
}
