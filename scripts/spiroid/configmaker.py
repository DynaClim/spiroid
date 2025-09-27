"""
This script generates JSON config files for input into spiroid.

A simulation input config file is generated for all combinations
of the values of the star and planet.

"""

import sys

sys.dont_write_bytecode = True

import json
from pathlib import Path
from spiroid.config import generate_all_configs, generate_all_effect_combinations


def make_config_files(simulation_name, all_configs, output_path):
    """Generate a simulation input config for all combinations of the star and planet values."""
    for i, config in enumerate(all_configs):
        # Append the name of the stellar evolution model, if used
        evolution_model = config["universe"]["central_body"]["kind"]["Star"].get("evolution", None)
        if evolution_model != "Disabled":
            evolution_model = (
                f"-{list(evolution_model.keys())[0].lower()}" if evolution_model else ""
            )
        else:
            evolution_model = ""

        config_name = f"{output_path}/{simulation_name}{evolution_model}_{i}.json.conf"
        print(f"Making config: {config_name}")
        with open(config_name, "x") as f:
            f.write(json.dumps(config, indent=4))


def make_configs(simulator_setup, effect_setup, planet_setup, star_setup, integrator_setup):
    """Generates a simulation configuration file for each combination of planets and stars."""
    if len(sys.argv) != 2:
        print("usage: python3 setup.py path/to/output/folder")
        return

    output_path = sys.argv[-1]
    # Create the output directory path if it doesn't already exist.
    Path(output_path).mkdir(parents=True, exist_ok=True)
    # Initialise the simulation properties.
    (simulation, disk_lifetime) = simulator_setup()

    # Initialise the integrator.
    integrator = integrator_setup()

    # Generates all combinations of enabled effects.
    all_effects = generate_all_effect_combinations(effect_setup())

    # Create a simulation for all combinations of planet, star, and effect values.
    for effects, effect_label in all_effects:
        simulation_name = f'{simulation["name"]}_{effect_label}'
        planet_base = planet_setup(effects)
        star_base = star_setup(effects)
        all_configs = generate_all_configs(
            simulation["start_time"],
            simulation["final_time"],
            disk_lifetime,
            planet_base,
            star_base,
            effects,
            integrator,
        )
        make_config_files(simulation_name, all_configs, output_path)
