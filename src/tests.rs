use super::*;
use postcard;
use pretty_assertions::assert_eq;
use sci_file::{
    OutputWriter, read_csv_columns_from_file, read_csv_rows_from_file, read_json_from_file,
};
use simulation::Integrator;
use simulation::simulation::InputConfig;
use std::path::Path;

fn test_simulation(config: impl AsRef<Path>) -> Universe {
    // Parse the config file.
    let mut config: InputConfig<Universe> = read_json_from_file(&config).unwrap();
    // Load stellar evolution data from file.
    if let ParticleType::Star(star) = &mut config.system.central_body.kind {
        // Load stellar evolution data from file if stellar evolution is enabled.
        if let Some(star_file) = star.evolution_file() {
            // Maps every row of the csv file into a `StarCsv`.
            let mut stellar_data = read_csv_rows_from_file::<StarCsv>(star_file).unwrap();
            // Configure the stellar evolution interpolator.
            let (star_ages, star_values) = StarCsv::initialise(&mut stellar_data);
            star.initialise_evolution(&star_ages, &star_values);
        }
    }

    // Load love number data from file(s) if kaula tides are enabled.
    if let Some(kaula) = config.system.orbiting_body.tides.kaula_get_mut() {
        if let Some(solid_file) = kaula.solid_file() {
            let love_solid = read_csv_columns_from_file::<f64>(solid_file).unwrap();
            kaula.initialise_love_number_solid(&love_solid);
        }
        if let Some(ocean_file) = kaula.ocean_file() {
            let love_ocean = read_csv_columns_from_file::<f64>(ocean_file).unwrap();
            kaula.initialise_love_number_ocean(&love_ocean);
        }
    }

    // Initialise the universe (star, planet, etc).
    config.system.initialise(config.initial_time).unwrap();

    // Initial values for the integrator.
    let y = config.system.integration_quantities();
    config
        .integrator
        .initialise(config.initial_time, config.final_time, &y);

    // Run the full integration.
    let _ = config.integrator.integrate(&mut config.system).unwrap();

    config.system
}

fn compare_or_create(path: impl AsRef<Path> + std::fmt::Display, result: &Universe) {
    match read_json_from_file::<Universe>(&path) {
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
                    let mut writer = OutputWriter::new(&path).unwrap();
                    writer.write(&result).unwrap();
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
fn serde_roundtrip() {
    let (config, _) = make_testcase_paths("all_effects");

    let config: InputConfig<Universe> = read_json_from_file(&config).unwrap();
    let universe = config.system;

    let tmp = serde_json::to_string(&universe).unwrap();
    let serde_universe: Universe = serde_json::from_str(&tmp).unwrap();

    assert_eq!(universe, serde_universe)
}

#[test]
fn postcard_roundtrip() {
    let (config, _) = make_testcase_paths("all_effects");

    let config: InputConfig<Universe> = read_json_from_file(&config).unwrap();
    let universe = config.system;

    let tmp = postcard::to_stdvec(&universe).unwrap();
    let postcard_universe: Universe = postcard::from_bytes(&tmp).unwrap();

    assert_eq!(universe, postcard_universe)
}

#[test]
fn postcard_vs_serde() {
    let (config, _) = make_testcase_paths("all_effects");

    let universe = test_simulation(config);

    let tmp = serde_json::to_string(&universe).unwrap();
    let serde_universe: Universe = serde_json::from_str(&tmp).unwrap();

    let tmp = postcard::to_stdvec(&universe).unwrap();
    let postcard_universe: Universe = postcard::from_bytes(&tmp).unwrap();

    assert_eq!(serde_universe, postcard_universe)
}

fn make_testcase_paths(testcase: &str) -> (String, String) {
    let config = format!("examples/{}.conf", testcase);
    let expected = format!("examples/{}_{}.json", testcase, "expected");
    (config, expected)
}

#[test]
fn example_no_effects_starevol() {
    let (config, expected) = make_testcase_paths("no_effects_starevol");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}

#[test]
fn example_star_equilibrium_tides() {
    let (config, expected) = make_testcase_paths("star_equilibrium_tides");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}

#[test]
fn example_magnetic() {
    let (config, expected) = make_testcase_paths("magnetic");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}

#[test]
fn example_magnetic_tides() {
    let (config, expected) = make_testcase_paths("magnetic_tides");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}

#[test]
fn example_planet_kaula_solid_tides() {
    let (config, expected) = make_testcase_paths("planet_kaula_solid_tides");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}

#[test]
fn example_all_effects() {
    let (config, expected) = make_testcase_paths("all_effects");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}

#[test]
fn example_mesa() {
    let (config, expected) = make_testcase_paths("mesa");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}

#[test]
fn example_zahn() {
    let (config, expected) = make_testcase_paths("star_zahn_tides_mesa");
    let result = test_simulation(&config);
    compare_or_create(&expected, &result);
}
