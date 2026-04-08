AU = 149597870700.0
SECONDS_IN_YEAR = 31557600.0
SOLAR_MASS = 1.9884158605722266e30
SOLAR_RADIUS = 6.957e8
SOLAR_MASS_OVER_JUPITER_MASS = 1.047_348_644e3
JUPITER_MASS = SOLAR_MASS / SOLAR_MASS_OVER_JUPITER_MASS
GRAVITATIONAL = 6.674_28e-11
GRAVITATIONAL_EARTH_MASS = 3.986_004_415e14
EARTH_MASS = GRAVITATIONAL_EARTH_MASS / GRAVITATIONAL
EARTH_RADIUS = 6.378_136_6e6
JUPITER_RADIUS = 7.1492e7
MERCURY_MASS = 3.3011e23
MERCURY_RADIUS = 2.4397e6

# Plots will _not_ be made for these quantities.
# Uncomment a quantity to prevent it being used in plots.
# Comment to use the quantity.
IGNORED_KEYS = [
    "time",
    "disk_is_dissipated",
    "disk_lifetime",
    "planet_density_ratio",
    "planet_eccentricity",
    "planet_inclination",
    "planet_is_destroyed",
    "planet_atmosphere",
    "planet_general_relativity",
    "planet_tides",
    "planet_tides_kaulatides_particle_type_solid_solid_file",
    "planet_longitude_ascending_node",
    "planet_luminosity",
    "planet_magnetic_field",
    "planet_magnetic_pressure",
    "planet_magnetism",
    "planet_mass",
    "planet_mean_motion",
    "planet_moment_of_inertia",
    "planet_orbit_lower_limit",
    # "planet_pericentre_omega",
    "planet_radius",
    "planet_radius_of_gyration_2",
    "planet_reduced_mass",
    "planet_roche_limit",
    "planet_semi_major_axis",
    "planet_spin",
    "planet_spin_inclination",
    "planet_wind",
    "star_age",
    "star_alfven_radius",
    "star_angular_momentum_redistribution",
    "star_convective_mass",
    "star_convective_moment_of_inertia",
    "star_convective_moment_of_inertia_derivative",
    "star_convective_radius",
    "star_convective_turnover_time",
    "star_convective_turnover_time_sun",
    "star_convective_zone_angular_momentum",
    "star_core_envelope_coupling_constant",
    "star_dynamical_tide_dissipation",
    "star_evolution",
    "star_evolved_change_semi_major_axis",
    "star_evolved_mass_loss_rate",
    "star_evolved_wind_torque",
    "star_evolution_interpolated_star_file_path",
    "star_luminosity",
    "star_magnetic_field",
    "star_magnetic_torque",
    "star_magnetism_wind_alfvenic_mach",
    "star_magnetism_wind_alfven_speed_at_alfven_radius",
    "star_magnetism_wind_azimuthal_velocity",
    "star_magnetism_wind_critical_radius",
    "star_magnetism_wind_critical_radius_div_alfven_radius",
    "star_magnetism_wind_footpoint_conductance",
    "star_magnetism_wind_integration_constant",
    "star_magnetism_wind_interaction",
    "star_magnetism_wind_magnetic_pressure",
    "star_magnetism_wind_magnetic_torque",
    "star_magnetism_wind_radial_magnetic_field",
    "star_magnetism_wind_speed_of_sound",
    "star_magnetism_wind_surface_wind_velocity",
    "star_magnetism_wind_wind_density",
    "star_magnetism_wind_wind_velocity",
    "star_mass",
    "star_mass_accretion_efficiency",
    "star_mass_loss_rate",
    "star_mass_transfer_envelope_to_core_torque",
    "star_radiative_mass",
    "star_radiative_mass_derivative",
    "star_radiative_moment_of_inertia",
    "star_radiative_zone_angular_momentum",
    "star_radius",
    "star_rossby",
    "star_spin",
    "star_tidal_frequency",
    "star_tidal_quality",
    "star_tidal_torque_convective",
    "star_tidal_torque",
    "star_tides",
    "star_tides_constanttimelag",
    "star_wind",
    "star_wind_orbital_angular_momentum_loss",
    "star_wind_torque",
]


def filter_keys(keys):
    """Remove INGORED_KEYS from the provided list of keys."""
    for key in IGNORED_KEYS:
        if key in keys:
            keys.remove(key)

    return keys


# Quantities and multiplicative inverse unit conversion
UNITS = {
    "time": {"label": "years", "unit": SECONDS_IN_YEAR},
    "planet_semi_major_axis": {"label": "AU", "unit": AU},
    "star_mass": {"label": "Msun", "unit": SOLAR_MASS},
    "star_radiative_zone_angular_momentum": {
        "label": "kg.m^2.s-1",
    },
    "star_convective_zone_angular_momentum": {
        "label": "kg.m^2.s-1",
    },
    "planet_spin": {
        "label": "rad.s",
    },
    "planet_eccentricity": {
        "label": "",
    },
    "planet_inclination": {
        "label": "rad",
    },
    "planet_longitude_ascending_node": {
        "label": "rad",
    },
    "planet_pericentre_omega": {
        "label": "rad",
    },
    "planet_spin_inclination": {
        "label": "rad",
    },
}


def get_units_label(key):
    """Returns the unit label for a given quantity."""
    return UNITS.get(key, {}).get("label", "")


def convert_units(key, val):
    """Applies unit conversion to the given quantity."""
    unit = UNITS.get(key, {}).get("unit", None)
    if not unit:
        return val
    else:
        return val / unit


KEY_PREFIXES = {
    "orbiting_body_kind_": "",
    "orbiting_body": "Planet",
    "central_body_kind_": "",
    "central_body": "Star",
    "derivatives_0": "star_radiative_zone_dt",
    "derivatives_1": "star_convective_zone_dt",
    "derivatives_2": "planet_spin_dt",
    "derivatives_3": "planet_orbital_eccentricity_dt",
    "derivatives_4": "planet_orbital_inclination_dt",
    "derivatives_5": "planet_longitude_ascending_node_dt",
    "derivatives_6": "planet_argument_of_periapsis_dt",
    "derivatives_8": "planet_spin_axis_inclination_dt",
}


def sanitise_key(key):
    """Replaces verbose prefixes from the flattened json into more readable names."""
    for prefix in KEY_PREFIXES:
        if key.startswith(prefix):
            key = key.replace(prefix, KEY_PREFIXES.get(prefix)).lower()
            break
    return key


def partition_keys(keys):
    """Partition keys into obvious plot categories, based on the key prefix."""
    prefixes = {x.split("_")[0]: [] for x in keys}
    for key in keys:
        for prefix in prefixes.keys():
            if key.startswith(prefix):
                prefixes[prefix].append(key)
    return prefixes
