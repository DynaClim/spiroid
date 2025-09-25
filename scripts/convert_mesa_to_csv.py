# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "mesa-reader",
#     "numpy",
#     "pandas",
#     "steda",
#     "tqdm",
# ]
# ///

"""This script processes stellar evolution data from MESA output into spiroid compatible CSV files.
It reads the data, computes numerical derivatives for physical quantities,
and filters out rows where changes are insignificant.
"""

import numpy as np
import mesa_reader as mr
import pandas as pd
from tqdm import tqdm
import os
import sys

from STEDA.MESA import loadMesa
from STEDA.constants import Rsun, Msun, Lsun


def usage():
    print("""usage: convert_MESA_to_csv.py star_mass mesa_directory
                    star_mass: float mass of star in Msun
                    mesa_directory: path to your mesa directory to load the LOG file from""")


def convert_values(mass, log):
    properties = {
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
            Ms,
            Rs,
            Ls,
            age,
            fase,
            r,
            rho,
            drho,
            rho_mean,
            mass,
            g,
            dg,
            dN2,
            N,
            dN,
            Kt,
            i_int_enve,
            i_int_core,
            lc,
            vc,
        ) = loadMesa(log, profile)

        h = log.history_data
        Mdot = -h.data_at_model_number(
            "star_mdot", log.model_with_profile_number(profile)
        )

        # Calculate the convective turnover time throughout the convective envelope
        if i_int_core < len(r) - 1 and i_int_core != 0:
            # # Calculate the convective turnover time at the interface between the radiative core and convective envelope
            # option 1
            # filter = [i >= i_int_core and i < i_int_core+30 for i in range(len(r))]
            # filter[-1] = False
            # tc_out = abs(np.percentile(lc[filter]/vc[filter],50))
            # option 2
            # filter = [i >= i_int_core and N[i] == 0. for i in range(len(r))]
            # filter[-1] = False
            # tc_out = abs(np.percentile(lc[filter]/vc[filter],50))

            # Calculate the convective turnover time at the center of the convective envelope
            i_cent = np.where(r >= 0.5 * (Rs + r[i_int_core]))[0][0]
            while (vc[i_cent] == 0.0 or N[i_cent] > 0.0) and i_cent > i_int_core:
                i_cent -= 1  # go down until we find a convective point
            tc_out = lc[i_cent] / vc[i_cent]
        else:
            tc_out = 1e99

        # adjusted_convective_mass = mass[i_int_core]/Ms
        # t_c_base_out = 10**(8.79 - 2. * abs(np.log10(adjusted_convective_mass))**(0.349) - 0.0194 * abs(np.log10(adjusted_convective_mass))**2 - 1.62 * min(np.log10(adjusted_convective_mass) + 8.55, 0.))

        properties["age"][profile - 1] = age  # yr
        properties["radius"][profile - 1] = Rs / Rsun
        properties["mass"][profile - 1] = Ms / Msun
        properties["convective_radius"][profile - 1] = r[i_int_core] / Rsun
        properties["radiative_mass"][profile - 1] = (
            (mass[i_int_core]) / Msun
        )  # MESA radiative mass
        integrand = rho * r**4
        properties["radiative_moment_of_inertia"][profile - 1] = max(
            8 * np.pi / 3 * np.trapz(integrand[:i_int_core], r[:i_int_core]), 1e44 * 1e7
        ) / (Ms * Rs**2)
        properties["convective_moment_of_inertia"][profile - 1] = max(
            8 * np.pi / 3 * np.trapz(integrand[i_int_core:], r[i_int_core:]), 1e44 * 1e7
        ) / (Ms * Rs**2)
        properties["luminosity"][profile - 1] = Ls / Lsun
        properties["convective_turnover_time"][profile - 1] = tc_out
        properties["mass_loss_rate"][profile - 1] = Mdot

    return pd.DataFrame(properties)


def save_mesa_to_csv(mass, df):
    output_csv = f"../examples/data/star/evolution/mesa_{10*mass:02.0f}.csv"
    df.to_csv(output_csv, index=False)


def filter_values(mass, df):
    # The goal is to reduce the dataset size by keeping only rows with meaningful evolution,
    # reducing computational cost for models that use this data.
    # The filtered data is then saved back to the original CSV file location.

    periods = len(df) // 100  # Adjust the period based on the length of the DataFrame

    # Compute numerical derivatives for all columns except the first, with respect to the index
    derivatives = df.iloc[:, 1:].diff(periods=periods)
    derivatives.columns = [f"d_{col}/d_index" for col in df.columns[1:]]

    # Exclude convective_turnover_time from filtering criteria
    columns = [
        col for col in derivatives.columns if col[2:-8] != "convective_turnover_time"
    ]

    # Calculate the maximum relative derivative for each row
    maxderiv = np.zeros(len(derivatives))
    for i in range(len(derivatives)):
        maxderiv[i] = max(
            [
                abs(derivatives[col][i] / df[col[2:-8]][i])
                for col in columns
                if df[col[2:-8]][i] != 0
            ]
        )

    keep_indices = maxderiv > 0.001
    keep_indices[:periods] = True  # Always keep the first 'periods' rows
    for i in range(periods, len(derivatives)):
        if (
            np.sum(keep_indices[i - periods : i]) == 0
        ):  # Ensure at least one row is kept in each segment
            keep_indices[i] = True

    # Filter the DataFrame and recompute derivatives for the filtered data
    filtered_df = df.iloc[keep_indices].reset_index(drop=True)
    filtered_derivatives = filtered_df.iloc[:, 1:].diff(periods=periods)
    filtered_derivatives.columns = [
        f"d_{col}/d_index" for col in filtered_df.columns[1:]
    ]
    filtered_maxderiv = np.zeros(len(filtered_derivatives))
    for i in range(len(filtered_derivatives)):
        filtered_maxderiv[i] = max(
            [
                abs(filtered_derivatives[col][i] / filtered_df[col[2:-8]][i])
                for col in columns
                if filtered_df[col[2:-8]][i] != 0
            ]
        )

    # Save the filtered DataFrame back to the original CSV file location
    output_csv = f"../examples/data/star/evolution/mesa_{10*mass:02.0f}.csv"
    filtered_df.to_csv(output_csv, index=False)


def main():
    if len(sys.argv) != 3:
        usage()
        exit()
    elif len(sys.argv) == 2:
        mesa_dir = sys.argv[1]
        mesa_file = os.path.join(mesa_dir, "LOGS")
        log = mr.MesaLogDir(mesa_file, memoize_profiles=False)
        mass = float(sys.argv[2])
        df = convert_values(mass, log)
        df = filter_values(mass, df)
        save_mesa_to_csv(mass, df)


if __name__ == "__main__":
    main()
