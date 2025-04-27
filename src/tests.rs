use super::*;
use pretty_assertions::assert_eq;
use sci_file::{deserialize_postcard_from_path,
create_buffered_file_reader,
    OutputFile, deserialize_csv_column_vectors_from_path, deserialize_csv_rows_from_path,
    deserialize_json_from_path,
};
use simulation::simulation::InputConfig;
use simulation::Integrator;
use std::path::{Path, PathBuf};
use std::io::BufRead;
use postcard;

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
    config
        .integrator
        .initialise(config.initial_time, config.final_time, &y);

    // Run the full integration.
    let _ = config.integrator.integrate(&mut config.universe).unwrap();

    config.universe
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

// Ensure that the same fields are being serialized/deserialized (serde) as encoded/decoded (bincode)
#[test]
fn bincode_vs_serde() {
    let path = "examples/all_effects.conf";
    let universe = test_simulation(path.into());

    let tmp = serde_json::to_string(&universe).unwrap();
    let serde_universe: Universe = serde_json::from_str(&tmp).unwrap();

    let tmp = bincode::encode_to_vec(&universe, bincode::config::standard()).unwrap();
    let (bincode_universe, _): (Universe, usize) =
        bincode::decode_from_slice(&tmp, bincode::config::standard()).unwrap();

    assert_eq!(bincode_universe, serde_universe)
}

#[test]
fn serde_roundtrip() {
    let path = "examples/all_effects.conf";
    let config: InputConfig<Universe> = deserialize_json_from_path(&path).unwrap();
    let universe = config.universe;

    let tmp = serde_json::to_string(&universe).unwrap();
    let serde_universe: Universe = serde_json::from_str(&tmp).unwrap();

    assert_eq!(universe, serde_universe)

}

#[test]
fn postcard_roundtrip() {
    let path = "examples/all_effects.conf";
    let config: InputConfig<Universe> = deserialize_json_from_path(&path).unwrap();
    let universe = config.universe;

    let tmp = postcard::to_stdvec(&universe).unwrap();
    let postcard_universe: Universe = postcard::from_bytes(&tmp).unwrap();

    assert_eq!(universe, postcard_universe)

}

#[test]
fn postcard_vs_serde() {
    let path = "examples/all_effects.conf";
    let universe = test_simulation(path.into());

    let tmp = serde_json::to_string(&universe).unwrap();
    let serde_universe: Universe = serde_json::from_str(&tmp).unwrap();

    let tmp = postcard::to_stdvec(&universe).unwrap();
    let postcard_universe: Universe = postcard::from_bytes(&tmp).unwrap();

    assert_eq!(serde_universe, postcard_universe)
}


#[test]
fn bincode_roundtrip() {
    let path = "examples/all_effects.conf";
    let config: InputConfig<Universe> = deserialize_json_from_path(&path).unwrap();
    let universe = config.universe;

    let tmp = bincode::encode_to_vec(&universe, bincode::config::standard()).unwrap();
    let (bincode_universe, _): (Universe, usize) =
        bincode::decode_from_slice(&tmp, bincode::config::standard()).unwrap();

    assert_eq!(universe, bincode_universe)
}

#[test]
fn bincode_multiread() {
    let path = "output/run_2/magnetic.bincode2";
    let mut reader = create_buffered_file_reader(&path).unwrap();
    let mut count = 0;
    loop {
        let peek = reader.fill_buf().unwrap();
        dbg!(peek.len());
        if peek.is_empty() {
            break;
        }
        let data: Universe = match bincode::decode_from_std_read(&mut reader, bincode::config::standard()) {
            Ok(data) => data,
            Err(err) => {
                dbg!(count);
                dbg!(&err);
                panic!("Error deserializing data: {}", err);

            },
        };
        // Do something with the deserialized data
        count += 1;
//        dbg!(&data, count);
    }
}

#[test]
fn postcard_multiread() {
    let path = "output/run_1/magnetic.postcard";
    let data: Vec<Universe> = deserialize_postcard_from_path(&path).unwrap();
    assert_eq!(data.len(), 25401);
}

#[test]
fn postcard_size() {
    let path = "output/run_1/magnetic.postcard";
    let data: Vec<Universe> = deserialize_postcard_from_path(&path).unwrap();

    let mut outfile = OutputFile::new(&"output/run_1/magnetic_loop.pcsz").unwrap();
    for d in &data {
        outfile.write(&d).unwrap();
    }
    assert!(data.len() != 0);
}

#[test]
fn jsonl_multiread() {
    let path = "output/run_0/magnetic.jsonl";
    let mut reader = create_buffered_file_reader(&path).unwrap();
    let data: Vec<Universe> = reader.lines().map(|line|
        serde_json::from_str(&line.unwrap()).unwrap()
    ).collect();
    assert_eq!(data.len(), 25401);
}

#[test]
fn example_no_effects() {
    let result = test_simulation("examples/no_effects.conf".into());
    compare_or_create("examples/no_effects_expected.json", &result);
}

#[test]
fn example_tides() {
    let result = test_simulation("examples/tides.conf".into());
    compare_or_create("examples/tides_expected.json", &result);
}

#[test]
fn example_magnetic() {
    let result = test_simulation("examples/magnetic.conf".into());
    compare_or_create("examples/magnetic_expected.json", &result);
}

#[test]
fn example_magnetic_tides() {
    let result = test_simulation("examples/magnetic_tides.conf".into());
    compare_or_create("examples/magnetic_tides_expected.json", &result);
}

#[test]
fn example_kaula_solid() {
    let result = test_simulation("examples/kaula_solid.conf".into());
    compare_or_create("examples/kaula_solid_expected.json", &result);
}

#[test]
fn example_all_effects() {
    let result = test_simulation("examples/all_effects.conf".into());
    compare_or_create("examples/all_effects_expected.json", &result);
}
