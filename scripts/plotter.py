"""
This script generates plot figures for JSONL (JSON Lines) data containing only numeric values.

If the only argument is a directory, plots will be generated for all JSONL files found (recursively).
If the first (and subsequent) arguments are JSONL file containing data, the last argument must be an output destination.

The JSONL must include a time variable (used for the x axis).

Example:
`plotter.py data1.jsonl data2.jsonl`

If the format of data1 and data2 is:
{"time": 0, "mass:" 1, "speed": 0},
{"time": 1, "mass:" 1, "speed": 11},
{"time": 3, "mass:" 7, "speed": 2}

Two types of plots would be created (one for each variable): mass vs time and speed vs time.
A plot for each variable is produced in both regular and logscale units.

The output would be:
mytitle-mass.png
mytitle-mass-logscale.png
mytitle-speed.png
mytitle-speed-logscale.png

"""

import sys

sys.dont_write_bytecode = True

import glob
import json
try:
    import lz4tools
except ImportError:
    pass
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from pathlib import Path
from flatten_json import flatten

from units import (
    sanitise_key,
    partition_keys,
    convert_units,
    get_units_label,
    filter_keys,
)


def parse_values(file_path):
    print(f"reading data file: {file_path})")
    if file_path.endswith(".jsonl.lz4"):
        data = lz4tools.open(file_path).read().split(b'\n')
        return parse_jsonl(data)
    if file_path.endswith(".jsonl"):
        with open(file_path, "r") as data:
            return parse_jsonl(data)
    else:
        return parse_values_lines(file_path)


def parse_jsonl(data):
    """Reads a JSONL file and combines all entries into a single dictionary."""
    dict = {}

    # Replace booleans with 0/1 so they can be plotted.
    lookup = {False: 0, True: 1}

    for line in data:
        data = json.loads(line)
        data = flatten(data)
        # Update the combined dictionary with the current JSON object
        for key, value in data.items():
            key = sanitise_key(key)
            value = lookup.get(value, value)

            # Convert units for specific values.
            value = convert_units(key, value)


            if key in dict:
                dict[key].append(value)
            else:
                dict[key] = [value]
    return dict


def parse_values_lines(file_path):
    """Parses lines of key:value pairs to a dictionary of {key:[value0, value1, ..., valuen]}."""
    # Open the file containing the PLOTLINE data
    with open(file_path, "r") as file:
        # Initialize variables for storing the x and y values
        values = {}
        # Read each line of the file
        for line in file:
            if line == "\n":
                continue
            # Extract the field and value
            line = [x.strip() for x in line.split(":")]
            key = line[0].lower()
            val = line[1]
            if val == "true":
                val = 1.0
            elif val == "false":
                val = 0.0
            else:
                val = float(line[1])

            # Create the list or add the value to the list
            if values.get(key) is None:
                values[key] = [val]
            else:
                values[key].append(val)
    return values


def create_plot(title, x_label, y_label, subplots, logscale=False):
    """Creates a figure containing specified subplots."""
    plt.xlabel(f"{x_label} ({get_units_label(x_label)})")
    plt.ylabel(f"{y_label} ({get_units_label(y_label)})")

    # Colors for each subplot
    colors = list(matplotlib.colors.XKCD_COLORS.keys())

    # Plot all subplots, assigning each subplot a distinct color.
    for data in subplots:
        (name, x, y) = data
        plt.plot(x, y, "", label=name, color=colors.pop(), alpha=0.5)

    # Convert to logscale if required.
    if logscale:
        plt.yscale("log")
        plt.xscale("log")

    # Set the title of the plot.
    plt.title(title)

    # Make legend for the plot.
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))


def save_plot(title, output_path):
    """Saves the plot to png."""
    # Save the figure as a file.
    plt.savefig(
        f"{output_path}/{title.lower().replace(' ', '_').replace('\n', ':')}.png",
        dpi=500,
        bbox_inches="tight",
    )


def create_plots(x_label, y_label, subplots, output_path):
    """Creates two figures (normal and logscale)."""
    # Create normal and logscale figures.
    for logscale in [True, False]:
        # Create a title for the figure.
        title = f"{x_label} vs {y_label}"
        if logscale:
            title += "-logscale"
        create_plot(title, x_label, y_label, subplots, logscale)
        save_plot(title, output_path)
        plt.cla()


def create_merged_subplots(x_label, y_label, data_sources):
    """Creates figures for each quantity containing combined plots for each of the datafiles."""
    subplots = []
    # Append the name and the x and y values for each data source as a subplot.
    for name, dict in data_sources.items():
        if y_label in dict:
            subplots.append((name, dict[x_label], dict[y_label]))

    return subplots


def create_subplots(x_label, y_labels, data):
    """Creates subplots for each quantity."""
    subplots = []
    # Append the name and the x and y values for each data source as a subplot.
    for y_label in y_labels:
        if y_label != x_label:
            subplots.append((y_label, data[x_label], data[y_label]))

    return subplots


def main():
    if len(sys.argv) < 2:
        print(
            "usage: plotter.py input_data [input_data1] [input_data2] [...] output_path"
        )
        exit()
    elif len(sys.argv) == 2:
        output_path = None
        path = sys.argv[1]
        all_files = glob.glob(f"{path}/**/*.jsonl*", recursive=True)
    else:
        output_path = sys.argv.pop()
        # Create the output directory path if it doesn't already exist.
        Path(output_path).mkdir(parents=True, exist_ok=True)
        all_files = sys.argv[1:]

    # Parse the input data into a dictionary:
    # {
    #   filename1: {var1: [values], var2: [values], ...},
    #   filename2: {var1: [values], var2: [values], ...},
    #   ...
    # }
    all_data = {file: parse_values(file) for file in all_files}
    x_label = "time"
    if output_path:
        # Create plots for each quantitiy, containing data from all data files.
        all_keys = set()
        # Collect list of all quantities to plot.
        for _, d in all_data.items():
            all_keys.update(set(d.keys()))

        # Remove unwanted keys (keys specified as unworth for plotting).
        all_keys = filter_keys(all_keys)
        for y_label in all_keys:
            print(f"Making graph: {y_label}")
            subplots = create_merged_subplots(x_label, y_label, all_data)
            create_plots(x_label, y_label, subplots, output_path)
    else:
        # Create a plot for each quantity for each data file.
        for path, data in all_data.items():
            # Set the output directory to save the plots
            output_path = os.path.dirname(path)
            all_keys = {*data.keys()}
            print(f"processing file: {output_path}")
            # Remove unwanted keys (keys specified as unworth for plotting).
            all_keys = filter_keys(all_keys)
            # Create a merged plot for all grouped quantities for each data file.
            for y_label, key_set in partition_keys(all_keys).items():
                if len(key_set) == 1:
                    continue
                print(f"Making graph: {y_label}")
                subplots = create_subplots(x_label, key_set, data)
                create_plots(x_label, y_label, subplots, output_path)


if __name__ == "__main__":
    main()
