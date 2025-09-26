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

from astropy.constants import R_sun, M_sun, L_sun

Rsun = R_sun.cgs.value   # Solar radius in CGS units
Msun = M_sun.cgs.value   # Solar mass in CGS units
Lsun = L_sun.cgs.value   # Solar luminosity in CGS units

def usage():
    print("""usage: convert_MESA_to_csv.py star_mass mesa_directory
                    star_mass: float mass of star in Msun
                    mesa_directory: path to your mesa directory to load the LOG file from""")


# Load Mesa data
#---------------
def loadMesa(log: mr.MesaLogDir, profile: int):
    """ Load the data from the mesa log files and the profile snapshot

    Args:
        log (MesaLogDir): directory to the mesa log files
        profile (int): profile number

    Returns:
        Tuple:
            Ms (float): Mass of the star
            Rs (float): Radius of the star
            Ls (float): Luminosity of the star
            age (float): Age of the star
            r (NDArray[single]): radius of the datapoints from the stellar evolution simulation
            rho (NDArray[single]): density of the datapoints from the stellar evolution simulation
            mass (NDArray[single]): mass inside the sphere of the datapoints from the stellar evolution simulation
            N (NDArray[single]): Brunt-Väisälä frequency of the datapoints from the stellar evolution simulation
            i_int_enve (float): index of the interface of the convective to the radiative zone (convective core radiative envelope)
            i_int_core (float): index of the interface of the convective to the radiative zone (radiative core convective envelope)
            lc (NDArray[single]): mixing length of the datapoints from the stellar evolution simulation
            vc (NDArray[single]): convective velocity of the datapoints from the stellar evolution
    """

    # Setting the mesa log director
    h = log.history_data
    p = log.profile_data(profile_number=profile)

    # Retreiving the mass and radius of the star at the given profile snapshot
    Ms  = h.data_at_model_number("star_mass",log.model_with_profile_number(profile))*Msun   # cgs
    Rs  = h.data_at_model_number("radius",log.model_with_profile_number(profile))*Rsun      # cgs
    Ls  = h.data_at_model_number("luminosity",log.model_with_profile_number(profile))*Lsun  # cgs
    age = h.data_at_model_number("star_age",log.model_with_profile_number(profile))         # yr

    # Retreiving the internal structure profiles for the given profile snapshot
    r         = np.flip(p.rmid)*Rsun                              # cgs
    rho       = np.flip(p.rho)                                    # cgs
    mass      = np.flip(p.mass)*Msun                              # cgs
    N2        = np.flip(p.brunt_N2)                               # cgs
    N         = np.sqrt(np.where(N2 < 0, 0, N2))                  # cgs
    Kt        = np.flip(p.thermal_diffusivity)                    # cgs
    Kt        = np.where(N2 < 0, 0, Kt)                           # cgs
    lc        = np.flip(p.mlt_mixing_length)                      # cgs
    vc        = np.flip(p.conv_vel)                               # cgs

    i_int_core, i_int_enve = calcTriLayer(N2, r, Rs)

    Kt = np.where(r < r[i_int_core], Kt, 0.) # Only take into account the radiative zone

    return Ms, Rs, Ls, age, r, rho, mass, N, i_int_enve, i_int_core, lc, vc


def calcTriLayer(N2, r, Rs):
    """ Calculate the indices of the interfaces between convective and radiative zones

    Args:
        N2 (NDArray[single]): squared Brunt-Väisälä frequency
        r (NDArray[single]): radius of the datapoints
        Rs (float): radius of the star

    Returns:
        Tuple[int, int]: indices of the interfaces between convective and radiative zones
    """

    # Calculate the interaction layers
    i_int_core, i_int_enve = 0, 0
    a = np.array(N2)
    asign = np.sign(a)
    asignroll = np.roll(asign, 1)
    asignroll[0] = 1 # to catch if the core is convective or radiative
    signchange = ((asignroll - asign) != 0).astype(int)
    where = np.where(signchange==1)[0]
    if len(where) == 0: # fully convective
        i_int_enve = 0
        i_int_core = 0
    elif where[0] == 0: # radiative core
        if len(where) == 1: # fully radiative
            i_int_enve = len(r)-1
            i_int_core = len(r)-1
        else: # convective envelope
            i_int_enve = 0
            i_int_core = 0
            for j in range(1, len(where)-1, 2): # convective core should be larger than 1e-4 of the star
                if r[where[j]] > 1e-4*Rs and (r[where[j+1]] - r[where[j]]) / r[where[j]] > 0.3:
                    i_int_enve = where[j]
                    break
            for i in range(j+1, len(where), 2): # take the first one that is not due to numerical instabilities
                if (r[where[i]] - r[where[i-1]]) / r[where[i]] > 0.4:
                    i_int_core = where[i]
                    break
    else: # convective core
        i_int_enve = 0
        i_int_core = 0
        where = list(where)
        where.append(len(r)-1)
        i = 0
        for i in range(len(where)-1, 1, -1): # take the first one that is not due to numerical instabilities
            if (r[where[i]] - r[where[i-1]]) / r[where[i-1]] > 0.3:
                if (i%2 == 1): i = i-1
                for j in range(i, 0, -2): # Take the first one inside the numerical instability
                    if (r[where[j]] - r[where[j-1]]) / r[where[j-1]] > 0.3:
                        i = j
                        break
                    else:
                        i = j-2
                i_int_core = where[i]
                break
        if i_int_core == 0: i_int_core = where[0]
        for i in range(i, 0, -2): # take the first one that is not due to numerical instabilities
            if (r[where[i]] - r[where[i-1]]) / r[where[i-1]] > 0.3:
                for j in range(i, 0, -2): # check if there is a numerical instability in the shell
                    if j > 1:
                        if (r[where[j-1]] - r[where[j-2]]) / r[where[j-1]] > 0.3:
                            if j == 2: i_int_enve = where[i-1] # If j = 2, the numerical error was inside the convective shell
                            else:      i_int_enve = where[j-1] # If j > 2, the numerical error was outside the convective shell
                            break
                    else:
                        i_int_enve = where[i-1]
                        break

    return i_int_core, i_int_enve


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
            r,
            rho,
            mass,
            N,
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
            8 * np.pi / 3 * np.trapezoid(integrand[:i_int_core], r[:i_int_core]), 1e44 * 1e7
        ) / (Ms * Rs**2)
        properties["convective_moment_of_inertia"][profile - 1] = max(
            8 * np.pi / 3 * np.trapezoid(integrand[i_int_core:], r[i_int_core:]), 1e44 * 1e7
        ) / (Ms * Rs**2)
        properties["luminosity"][profile - 1] = Ls / Lsun
        properties["convective_turnover_time"][profile - 1] = tc_out
        properties["mass_loss_rate"][profile - 1] = Mdot

    return pd.DataFrame(properties)


def save_mesa_to_csv(mass, df):
    output_csv = f"./examples/data/star/evolution/mesa_{10*mass:02.0f}.csv"
    df.to_csv(output_csv, index=False)


def filter_values(df):
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

    return df


def main():
    if len(sys.argv) != 3:
        usage()
        exit()
    elif len(sys.argv) == 3:
        mass = float(sys.argv[1])
        mesa_dir = sys.argv[2]
        mesa_file = os.path.join(mesa_dir, "LOGS")
        log = mr.MesaLogDir(mesa_file, memoize_profiles=False)
        df = convert_values(mass, log)
        df = filter_values(df)
        save_mesa_to_csv(mass, df)


if __name__ == "__main__":
    main()
