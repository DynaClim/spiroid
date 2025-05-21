"""
Initial conditions for input into spiroid.

A simulation input config file is generated for all combinations
of the values of the star, planet and effects.

"""

import sys

sys.dont_write_bytecode = True
from spiroid.configmaker import make_configs
from units import AU, SECONDS_IN_YEAR, SOLAR_MASS


def simulator_setup():
    ##############################################################
    ####################### SIMULATOR SETUP ######################
    ##############################################################

    # The prefix simulation name.
    # Produces test_1.json.conf, test_2.json.conf, test_3.json.conf etc.
    simulation = {
        "name": "test",
        # Simulation start time, seconds (from years)
        "start_time": SECONDS_IN_YEAR * 1.0e6,
        # Simulation end time, seconds (from years)
        "final_time": SECONDS_IN_YEAR * 1.0e9,
    }

    # seconds (from years)
    disk_lifetime = SECONDS_IN_YEAR * 2.482e6

    return (simulation, disk_lifetime)


def effect_setup():
    # Enables or disables certain effects for all simulations.
    # Must be [True], [False] or [True, False]
    effects = {
        "MAGNETIC_EFFECT_ENABLED": [True, False],
        "STAR_EVOLUTION_ENABLED": [True, False],
        # Constant Time Lag stellar tide
        "STAR_TIDES_ENABLED": [True, False],
        # Kaula planetary tides
        "PLANET_TIDES_ENABLED": [False],
    }

    return effects


def planet_setup(effects):
    ##############################################################
    ####################### PLANET SETUP #########################
    ##############################################################
    planet_base = {
        # kg
        "mass": [1.898e26],
        # m
        "radius": [3.255e7],
        # m (from AU)
        "semi_major_axis": [AU * x for x in [0.019]],
        # Gauss
        "magnetic_field": [10.0],
        "is_destroyed": [False],
    }

    if effects["PLANET_TIDES_ENABLED"]:
        # For Kaula
        planet_base.update(
            {
                # rad.s
                "spin": [8.093879511357418e-07],
                # No units
                "eccentricity": [0.005],
                # rad
                "inclination": [0.3490658503988659],
                # rad
                "longitude_ascending_node": [1.0],
                # rad
                "pericentre_omega": [0.0],
                # rad
                "spin_inclination": [0.34906584951436426],
                # No units
                "radius_of_gyration": [0.33070368308499226],
                "type_and_file": [
                    (
                        "Solid",
                        "examples/data/planet/tides/kaula/leconte2015_steinberger.csv",
                    )
                ],
            }
        )

    return planet_base


def star_setup(effects):
    ##############################################################
    ####################### STAR SETUP ###########################
    ##############################################################
    star_base = {
        # kg (from Msun)
        "mass": [SOLAR_MASS * x for x in [0.8]],
        # rad.s-1
        "spin": [5.194e-05],
        # seconds (from years)
        "core_envelope_coupling_constant": [SECONDS_IN_YEAR * x for x in [1.171e7]],
        # Ohm-1
        "footpoint_conductance": [5.8e4],
        "star_file_path": [None],  # Do not edit.
        "sigma_bar": [None],  # Do not edit.
    }

    if effects["STAR_EVOLUTION_ENABLED"]:
        star_base["star_file_path"] = ["examples/data/star/evolution/savgol_08.csv"]
    else:
        # Set the initial star values that would otherwise be provided by savgol data if evolution were enabled.
        # Must be non-zero (to avoid NaN).
        # No units
        star_base["radiative_moment_of_inertia"] = [1.0]
        star_base["convective_moment_of_inertia"] = [1.0]

    if effects["STAR_TIDES_ENABLED"]:
        star_base["sigma_bar"] = [1.0e-6]

    return star_base



def integrator_setup():
    ##############################################################
    #################### INTEGRATOR SETUP ########################
    ##############################################################

    # Output filtered each timestep that occurs greater than `absolute_tolerance`
    # or if the percentage difference between two timesteps exceeds `relative_tolerance`,
    # in which case the `relative_tolerance` is increased by the additive `scaling_factor`
    # Used to simulate logscale output.
    filter = {
        "Filtered": {
            "absolute_tolerance": SECONDS_IN_YEAR * 1e6,
            "relative_tolerance": 1.15,
            "incremental_scaling_factor": 0.002,
            "decremental_scaling_factor": 0.0
        }
    }

    odex = {
        "Odex": {
            "step_size_reduction_factor": 0.66666666666666666666,
            "step_size_selection_b": 2.0,
            "step_size_max": SECONDS_IN_YEAR * 5e5,
            "max_integration_steps": 100000000,
            # Uncomment the entire solution_output for Dense output
            "solution_output": filter,

        }
    }

    # Kaula tides on the planet modifies the integrator defaults to avoid timestep issue.
    odex_kaula = {
        "Odex": {
            "step_size_reduction_factor": 0.66666666666666666666,
            "step_size_selection_b": 2.0,
            "step_size_max": SECONDS_IN_YEAR * 1e3,
            "max_integration_steps": 500000000,
            "step_control_safety_a": 0.05,
            "step_control_safety_b": 0.2,
        }
    }

    dopri853 = {
        "Dopri853": {
            "step_size_controller": {
                "relative_tolerance": 1e-10,
                "absolute_tolerance": 1e-10,
                "step_size_factor_min": 0.3333333333333333,
                "step_size_factor_max": 6.0,
                "step_size_error_factor": 0.9,
                "step_size_max": SECONDS_IN_YEAR * 5e5,
                "alpha": 0.125,
                "beta": 0.0,
            },
            "step_size_underflow": None,
            "stiffness_test": "Disabled",
            "max_integration_steps": 100000000,
            # Uncomment the entire solution_output for Dense output
            "solution_output": filter,
        }
    }

    # Uncomment only the desired integrator.
    return odex
#    return odex_kaula
#    return dopri853


if __name__ == "__main__":
    make_configs(
        simulator_setup, effect_setup, planet_setup, star_setup, integrator_setup
    )
