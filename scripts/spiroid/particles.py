from units import (
    AU,
    JUPITER_MASS,
    JUPITER_RADIUS,
    MERCURY_MASS,
    MERCURY_RADIUS,
)

# Source: IAU / NASA planetary fact sheets
mercury = {
    # kg
    "mass": [MERCURY_MASS],
    # m
    "radius": [MERCURY_RADIUS],
    # m (from AU)
    "semi_major_axis": [AU * 0.387098],
    "magnetic_field": [10.0],
    # No units
    "eccentricity": [0.20563],
    # rad (initial value)
    "pericentre_omega": [0.0],
}

jupiter = {
    # kg
    "mass": [JUPITER_MASS],
    # m
    "radius": [JUPITER_RADIUS],
    # m (from AU)
    "semi_major_axis": [AU * 5.2],
    # Gauss
    "magnetic_field": [10.0],
    # rad.s
    "spin": [1.76e-4],
    # No units
    "eccentricity": [0.0484],
    # rad
    "inclination": [0.022776547],
    # rad
    "longitude_ascending_node": [1.755],
    # rad
    "pericentre_omega": [0.257],
    # rad
    "spin_inclination": [0.02278],
    # No units
    "radius_of_gyration": [0.3307],
}
