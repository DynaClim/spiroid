# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "mesa-reader",
#     "numpy",
#     "pandas",
#     "astropy",
#     "tqdm",
# ]
# ///

"""
This script processes stellar evolution data from MESA output into spiroid compatible CSV files.
It computes numerical derivatives for physical quantities and filters out rows where changes are insignificant.
"""

import numpy as np
import mesa_reader as mr
import pandas as pd
from tqdm import tqdm
import os
import sys

from astropy.constants import R_sun, M_sun, L_sun

SOLAR_RADIUS_CGS = R_sun.cgs.value  # Solar radius in CGS units
SOLAR_MASS_CGS = M_sun.cgs.value  # Solar mass in CGS units
SOLAR_LUMINOSITY_CGS = L_sun.cgs.value  # Solar luminosity in CGS units


def usage():
    print("""usage: convert_MESA_to_csv.py star_mass mesa_directory
       mesa_directory (string): path to your mesa directory to load the LOGS from""")


def load_mesa(log: mr.MesaLogDir, profile_number: int):
    """Load the data from the mesa log files and the profile snapshot

    Args:
        log (MesaLogDir): directory to the mesa log files
        profile_number (int): profile number


    Units: star_age in years, others in CGS units

    Returns:
        Tuple:
            star_mass (float): Mass of the star
            star_radius (float): Radius of the star
            star_luminosity (float): Luminosity of the star
            star_age (float): Age of the star in years
            radius (NDArray[single]): radius of the datapoints from the stellar evolution simulation
            density (NDArray[single]): density of the datapoints from the stellar evolution simulation
            mass (NDArray[single]): mass inside the sphere of the datapoints from the stellar evolution simulation
            bv_frequency (NDArray[single]): Brunt-Väisälä frequency of the datapoints from the stellar evolution simulation
            radiative_envelope_interface_i (int): index of the interface of the convective to the radiative zone (convective core radiative envelope)
            convective_envelope_interface_i (int): index of the interface of the convective to the radiative zone (radiative core convective envelope)
            mixing_length (NDArray[single]): mixing length of the datapoints from the stellar evolution simulation
            convective_velocity (NDArray[single]): convective velocity of the datapoints from the stellar evolution
    """

    # Setting the mesa log director
    history = log.history_data
    profile = log.profile_data(profile_number=profile_number)

    # Retreiving the mass and radius of the star at the given profile snapshot
    star_mass = (
        history.data_at_model_number("star_mass", log.model_with_profile_number(profile_number))
        * SOLAR_MASS_CGS
    )
    star_radius = (
        history.data_at_model_number("radius", log.model_with_profile_number(profile_number))
        * SOLAR_RADIUS_CGS
    )
    star_luminosity = (
        history.data_at_model_number("luminosity", log.model_with_profile_number(profile_number))
        * SOLAR_LUMINOSITY_CGS
    )
    star_age = history.data_at_model_number(
        "star_age", log.model_with_profile_number(profile_number)
    )

    star_mass_loss_rate = -history.data_at_model_number(
        "star_mdot", log.model_with_profile_number(profile_number)
    )

    # Retreiving the internal structure profiles for the given profile snapshot
    # All in CGS units
    radius = np.flip(profile.rmid) * SOLAR_RADIUS_CGS
    density = np.flip(profile.rho)
    mass = np.flip(profile.mass) * SOLAR_MASS_CGS
    bv_frequency_2 = np.flip(profile.brunt_N2)
    bv_frequency = np.sqrt(np.where(bv_frequency_2 < 0, 0, bv_frequency_2))
    mixing_length = np.flip(profile.mlt_mixing_length)
    convective_velocity = np.flip(profile.conv_vel)

    convective_envelope_interface_i, radiative_envelope_interface_i = calculate_tri_layer(
        bv_frequency_2, radius, star_radius
    )

    return (
        star_mass,
        star_radius,
        star_luminosity,
        star_age,
        star_mass_loss_rate,
        radius,
        density,
        mass,
        bv_frequency,
        radiative_envelope_interface_i,
        convective_envelope_interface_i,
        mixing_length,
        convective_velocity,
    )


def calculate_tri_layer(bv_frequency_2, radius, star_radius):
    """Calculate the indices of the interfaces between convective and radiative zones

    Args:
        bv_frequency_2 (NDArray[single]): squared Brunt-Väisälä frequency
        radius (NDArray[single]): radius of the datapoints
        star_radius (float): radius of the star

    Returns:
        Tuple[int, int]: indices of the interfaces between convective and radiative zones
    """

    # calculate the interaction layers
    radiative_envelope_interface_i = 0
    convective_envelope_interface_i = 0

    sign = np.sign(np.array(bv_frequency_2))
    signroll = np.roll(sign, 1)
    # to catch if the core is convective or radiative
    signroll[0] = 1
    sign_compare = ((signroll - sign) != 0).astype(int)
    # check where the sign of the bv_frequency_2 changes
    sign_change = np.where(sign_compare == 1)[0]
    # fully convective
    if len(sign_change) == 0:
        # indices already zero-initialised
        pass
    # radiative core
    elif sign_change[0] == 0:
        # fully radiative
        if len(sign_change) == 1:
            radiative_envelope_interface_i = len(radius) - 1
            convective_envelope_interface_i = len(radius) - 1
        # convective envelope
        else:
            for j in range(1, len(sign_change) - 1, 2):
                # convective core should be larger than 1e-4 of the star
                if (
                    radius[sign_change[j]] > 1e-4 * star_radius
                    and (radius[sign_change[j + 1]] - radius[sign_change[j]])
                    / radius[sign_change[j]]
                    > 0.3
                ):
                    radiative_envelope_interface_i = sign_change[j]
                    break
            for i in range(j + 1, len(sign_change), 2):
                # take the first one that is not due to numerical instabilities
                if (radius[sign_change[i]] - radius[sign_change[i - 1]]) / radius[
                    sign_change[i]
                ] > 0.4:
                    convective_envelope_interface_i = sign_change[i]
                    break
    # convective core
    else:
        sign_change = list(sign_change)
        sign_change.append(len(radius) - 1)
        # When len(sign_change) is too low, the loop will not be executed.
        # So i must be initialised to zero as is reused in the next loop.
        i = 0
        for i in range(len(sign_change) - 1, 1, -1):
            # take the first one that is not due to numerical instabilities
            if (radius[sign_change[i]] - radius[sign_change[i - 1]]) / radius[
                sign_change[i - 1]
            ] > 0.3:
                if i % 2 == 1:
                    i = i - 1
                for j in range(i, 0, -2):
                    # take the first one inside the numerical instability
                    if (radius[sign_change[j]] - radius[sign_change[j - 1]]) / radius[
                        sign_change[j - 1]
                    ] > 0.3:
                        i = j
                        break
                    else:
                        i = j - 2
                convective_envelope_interface_i = sign_change[i]
                break
        if convective_envelope_interface_i == 0:
            convective_envelope_interface_i = sign_change[0]
        for i in range(i, 0, -2):
            # take the first one that is not due to numerical instabilities
            if (radius[sign_change[i]] - radius[sign_change[i - 1]]) / radius[
                sign_change[i - 1]
            ] > 0.3:
                for j in range(i, 0, -2):
                    # check if there is a numerical instability in the shell
                    if j > 1:
                        if (radius[sign_change[j - 1]] - radius[sign_change[j - 2]]) / radius[
                            sign_change[j - 1]
                        ] > 0.3:
                            if j == 2:
                                # If j = 2, the numerical error was inside the convective shell
                                radiative_envelope_interface_i = sign_change[i - 1]
                            else:
                                # If j > 2, the numerical error was outside the convective shell
                                radiative_envelope_interface_i = sign_change[j - 1]
                            break
                    else:
                        radiative_envelope_interface_i = sign_change[i - 1]
                        break

    return convective_envelope_interface_i, radiative_envelope_interface_i


def convert_values(log):
    """Convert the values from raw MESA output into a dataframe matching the CSV struture, in SI units."""
    star = {
        "age": np.zeros(len(log.profile_numbers)),
        "radius": np.zeros(len(log.profile_numbers)),
        "mass": np.zeros(len(log.profile_numbers)),
        "convective_radius": np.zeros(len(log.profile_numbers)),
        "radiative_mass": np.zeros(len(log.profile_numbers)),
        "radiative_moment_of_inertia": np.zeros(len(log.profile_numbers)),
        "convective_moment_of_inertia": np.zeros(len(log.profile_numbers)),
        "luminosity": np.zeros(len(log.profile_numbers)),
        "convective_turnover_time": np.zeros(len(log.profile_numbers)),
        "mass_loss_rate": np.zeros(len(log.profile_numbers)),
    }

    for profile in tqdm(log.profile_numbers):
        (
            star_mass,
            star_radius,
            star_luminosity,
            star_age,
            star_mass_loss_rate,
            radius,
            density,
            mass,
            bv_frequency,
            radiative_envelope_interface_i,
            convective_envelope_interface_i,
            mixing_length,
            convective_velocity,
        ) = load_mesa(log, profile)

        # calculate the convective turnover time throughout the convective envelope
        if (
            convective_envelope_interface_i < len(radius) - 1
            and convective_envelope_interface_i != 0
        ):
            # calculate the convective turnover time at the center of the convective envelope
            i_cent = np.where(
                radius >= 0.5 * (star_radius + radius[convective_envelope_interface_i])
            )[0][0]
            while (
                convective_velocity[i_cent] == 0.0 or bv_frequency[i_cent] > 0.0
            ) and i_cent > convective_envelope_interface_i:
                # go down until we find a convective point
                i_cent -= 1
            tc_out = mixing_length[i_cent] / convective_velocity[i_cent]
        else:
            tc_out = 1e99

        index = profile - 1
        star["age"][index] = star_age
        star["radius"][index] = star_radius / SOLAR_RADIUS_CGS
        star["mass"][index] = star_mass / SOLAR_MASS_CGS
        star["convective_radius"][index] = (
            radius[convective_envelope_interface_i] / SOLAR_RADIUS_CGS
        )
        # MESA radiative mass
        star["radiative_mass"][index] = (mass[convective_envelope_interface_i]) / SOLAR_MASS_CGS
        integrand = density * radius**4
        star["radiative_moment_of_inertia"][index] = max(
            8
            * np.pi
            / 3
            * np.trapezoid(
                integrand[:convective_envelope_interface_i],
                radius[:convective_envelope_interface_i],
            ),
            1e44 * 1e7,
        ) / (star_mass * star_radius**2)
        star["convective_moment_of_inertia"][index] = max(
            8
            * np.pi
            / 3
            * np.trapezoid(
                integrand[convective_envelope_interface_i:],
                radius[convective_envelope_interface_i:],
            ),
            1e44 * 1e7,
        ) / (star_mass * star_radius**2)
        star["luminosity"][index] = star_luminosity / SOLAR_LUMINOSITY_CGS
        star["convective_turnover_time"][index] = tc_out
        star["mass_loss_rate"][index] = star_mass_loss_rate

    return pd.DataFrame(star)


def filter_values(df):
    """Reduce the dataset size by keeping only rows with meaningful evolution."""

    # Adjust the period based on the length of the DataFrame
    periods = len(df) // 100

    # Compute numerical derivatives for all columns except the first, with respect to the index
    derivatives = df.iloc[:, 1:].diff(periods=periods)
    derivatives.columns = [f"d_{col}/d_index" for col in df.columns[1:]]

    # Exclude convective_turnover_time from filtering criteria
    columns = [col for col in derivatives.columns if col[2:-8] != "convective_turnover_time"]

    # Calculate the maximum relative derivative for each row
    derivative_max = np.zeros(len(derivatives))
    for i in range(len(derivatives)):
        derivative_max[i] = max(
            [
                abs(derivatives[col][i] / df[col[2:-8]][i])
                for col in columns
                if df[col[2:-8]][i] != 0
            ]
        )

    keep_indices = derivative_max > 0.001
    # Always keep the first 'periods' rows
    keep_indices[:periods] = True
    for i in range(periods, len(derivatives)):
        # Ensure at least one row is kept in each segment
        if np.sum(keep_indices[i - periods : i]) == 0:
            keep_indices[i] = True

    # Filter the DataFrame and recompute derivatives for the filtered data
    filtered_df = df.iloc[keep_indices].reset_index(drop=True)
    filtered_derivatives = filtered_df.iloc[:, 1:].diff(periods=periods)
    filtered_derivatives.columns = [f"d_{col}/d_index" for col in filtered_df.columns[1:]]
    filtered_derivative_max = np.zeros(len(filtered_derivatives))
    for i in range(len(filtered_derivatives)):
        filtered_derivative_max[i] = max(
            [
                abs(filtered_derivatives[col][i] / filtered_df[col[2:-8]][i])
                for col in columns
                if filtered_df[col[2:-8]][i] != 0
            ]
        )

    return df


def save_mesa_to_csv(df):
    """Save the MESA data as CSV, ready for spiroid."""
    output_csv = f"./examples/data/star/evolution/mesa_{10*df.mass[0]:02.0f}.csv"
    df.to_csv(output_csv, index=False)


def main():
    # user must provide star_mass and mesa_dir_path
    if len(sys.argv) != 2:
        usage()
        exit()
    elif len(sys.argv) == 2:
        # load the mesa log file
        mesa_dir = os.path.join(sys.argv[1], "LOGS")
        if not os.path.isdir(mesa_dir):
            print(f"invalid directory path: {sys.argv[1]}")
            usage()
            exit()
        log = mr.MesaLogDir(mesa_dir, memoize_profiles=False)

        # convert the units and data format
        df = convert_values(log)
        # remove data points that are insiginificant
        df = filter_values(df)
        # write out the CSV data file
        save_mesa_to_csv(df)


if __name__ == "__main__":
    main()
