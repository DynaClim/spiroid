"""
Initial conditions to reproduce the Mercury GR apsidal precession.

A simulation input config file is generated for the Sun-Mercury system
with only General Relativity 1PN precession enabled. All other effects
(tides, magnetism, wind) are disabled.

Expected result: pericentre_omega precesses at ~42.98 arcsec/century,
consistent with the GR prediction (Einstein 1915).

"""

##############################################################
########################!!! WARNING !!!#######################
##############################################################
"""
Do _NOT_ change the order of dictionary keys in this file.
The values are unpacked in an order dependent way in config.py
to create the initial conditions.
"""
##############################################################
########################!!! WARNING !!!#######################
##############################################################


import sys

sys.dont_write_bytecode = True
from spiroid.configmaker import make_configs
from units import (
    AU,
    SECONDS_IN_YEAR,
    SOLAR_MASS,
    SOLAR_RADIUS,
)


def simulator_setup():
    ##############################################################
    ####################### SIMULATOR SETUP ######################
    ##############################################################

    simulation = {
        # The prefix simulation name.
        "name": "mercury",
        # Description of the science case.
        "decription": "Sun-Mercury system: GR apsidal precession (~42.98 arcsec/century)",
        # Simulation start time, seconds (from years).
        # Start at 1 Gyr so the disk is long dissipated.
        "start_time": SECONDS_IN_YEAR * 1.0e9,
        # Simulation end time, seconds (from years).
        # Integrate for 1 Myr to accumulate measurable precession.
        "final_time": SECONDS_IN_YEAR * (1.0e9 + 1.0e6),
    }

    # seconds (from years) — short disk, dissipated well before start_time
    disk_lifetime = SECONDS_IN_YEAR * 1.0e6

    return (simulation, disk_lifetime)


def effect_setup():
    # Enables or disables certain effects for all simulations.
    # Must be [True], [False] or [True, False].
    effects = {
        "MAGNETIC_EFFECT_ENABLED": [False],
        "STAR_EVOLUTION_ENABLED": [True],
        # Constant Time Lag stellar tide
        "STAR_TIDES_ENABLED": [False],
        # Kaula planetary tides
        "PLANET_TIDES_ENABLED": [False],
        # Disable wind
        "WIND_ENABLED": [False],
        # General Relativity 1PN apsidal precession
        "GR_ENABLED": [True],
    }

    return effects


def planet_setup(effects):
    ##############################################################
    ####################### PLANET SETUP #########################
    ##############################################################
    # Source: IAU / NASA planetary fact sheets
    planet_base = {
        # kg
        "mass": [3.3011e23],
        # m
        "radius": [2.4397e6],
        # m (from AU)
        "semi_major_axis": [AU * 0.387098],
        "magnetic_field": [None],  # Do not edit.
    }

    if effects.get("GR_ENABLED"):
        planet_base.update(
            {
                # No units
                "eccentricity": [0.20563],
                # rad (initial value)
                "pericentre_omega": [0.0],
            }
        )

    return planet_base


def star_setup(_effects):
    ##############################################################
    ####################### STAR SETUP ###########################
    ##############################################################
    star_base = {
        "mass": [None],  # Do not edit.
        "radius": [None],  # Do not edit.
        # rad.s-1 (solar rotation rate; irrelevant here — no tides/magnetism)
        "spin": [2.865e-6],
        # seconds (from years)
        "core_envelope_coupling_constant": [SECONDS_IN_YEAR * 1.171e7],
        "footpoint_conductance": [None],  # Do not edit.
        "evolution": [None],  # Do not edit.
        "sigma_bar": [None],  # Do not edit.
    }

    # Star evolution disabled: set initial stellar values manually.
    # Must be non-zero (to avoid NaN).
    star_base["mass"] = [SOLAR_MASS]
    star_base["radius"] = [SOLAR_RADIUS]
    star_base["radiative_moment_of_inertia"] = [1.0]
    star_base["convective_moment_of_inertia"] = [1.0]

    return star_base


def integrator_setup():
    ##############################################################
    #################### INTEGRATOR SETUP ########################
    ##############################################################

    dopri853 = {
        "Dopri853": {
            "step_controller": {
                "relative_tolerance": 1e-10,
                "absolute_tolerance": 1e-10,
                "step_size_factor_min": 0.3333333333333333,
                "step_size_factor_max": 6.0,
                "step_size_error_factor": 0.9,
                # step_size_max: 1000 yr — fine enough to resolve GR precession
                "step_size_max": SECONDS_IN_YEAR * 1e3,
                "alpha": 0.125,
                "beta": 0.0,
            },
            "step_size_underflow": None,
            "stiffness_test": "Disabled",
            "max_integration_steps": 100000000,
        }
    }

    return dopri853


if __name__ == "__main__":
    make_configs(simulator_setup, effect_setup, planet_setup, star_setup, integrator_setup)
