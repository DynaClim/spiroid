<img src="images/SNSF.png" align="left"> <img src="images/AGPLv3.png" align="right">

<img src="images/SPIROID.png"/>

> **S**tar **P**lanet **I**nteraction **R**otational and **O**rbital **I**ntegrator for **D**ynamics

**spiroid** is a fast numerical simulator written in Rust that models the evolution of a planet around a star. It is the successor to [ESPEM](https://doi.org/10.48550/arXiv.1811.06354), taking into account many physical processes:

- tide raised by the planet in the star following [Benbakoura et al. 2019](https://doi.org/10.1051/0004-6361/201833314)
- magnetic interactions between star and planet following [Ahuir et al. 2021](https://doi.org/10.1051/0004-6361/201936974)
- tide raised by the star in the planet following [Revol et al. 2023](https://doi.org/10.1051/0004-6361/202245790)


The physics has been validated by the following contributors:

- Stellar evolution and magnetic interaction by [Antoine Strugarek, CEA](https://tonione.github.io/antoinestrugarek.github.io/)
- Stellar and planetary formulations by [Émeline Bolmont, UniGE](https://www.unige.ch/sciences/astro/exoplanets/en/team/faculty-members/emeline-bolmont/) and [Leon Ka-Wang Kwok, UniGE](https://www.unige.ch/sciences/astro/exoplanets/en/team/scientific-collaborators/kwok-leon/)
- Tides in Kaula formalism by [Alexandre Revol, UniGE](https://www.unige.ch/sciences/astro/exoplanets/en/team/scientific-collaborators/revol-alexandre/)

With all effects enabled, `spiroid` can simulate between `2e6` and `2e7` years per second (i.e. a simulation of `1e7` years completes in 0.5 seconds).

## Installation
### Requirements

- Rust: [see rustup](https://www.rustup.rs/)

Clone the repository, build, and install:

```bash
cargo install --path .
```
The executable will be copied into `$HOME/.cargo/bin/`. You may need to add this directory to your `$PATH`.

> Note: This project uses Rust edition 2024. Building will require using the [nightly toolchain](https://rust-lang.github.io/rustup/concepts/channels.html) until ~Q2 2025. Hint: `rustup default nightly`

## Quickstart Example
```
# Create simulation case(s) into the "my_cases" directory.
python3 scripts/setup.py my_cases

# Launch all simulations from "my_cases" directory, putting results into "my_output" directory.
spiroid -b my_cases my_output

# Plot all simulations from "my_output" directory into individual plots.
python3 scripts/plotter.py my_output

# Plot all simulations from "my_output" directory into comparison plots.
python3 scripts/plotter.py my_output/*/*.jsonl my_output
```

## Usage

Initial conditions are specified into input configuration files (`.conf`) as JSON, which is read by `spiroid` to start a simulation. Output data from the simulation is in JSONL format. A simple python script is included to help in generating the JSON cases.

### Create a JSON case
Set the simulation, planet and star properties in `scripts/setup.py` to the desired value and enable or disable desired effects. Adding multiple values to a variable (e.g. `mass: [0.8, 1.0]`) will create a simulation for each value. Then create the `spiroid` config file(s):

`python3 scripts/setup.py output_directory`

This generates a JSON case (`.conf`) for each set of input conditions into the designated `output_directory`. The total number of simulations will be the combinatorial product of the variable values and effect states (enabled or disabled). 

### Start a simulation
Simulations can be launched individually, or in batch mode. The user must specify the location of the input file(s) (`.conf`) and the desired output location (`output_directory`).

#### Start a single simulation
`spiroid input_file.conf output_directory`

With input arguments:

- `input_file.conf`: Path to the desired JSON input configuration file (e.g. `./relative/path/to/config.conf` or `/home/$USER/simulations/example.conf`). 
- `output_directory`: Path to the desired output destination. Will be automatically created if it does not already exist.

#### Start multiple simulations (batch mode)
Batch mode launches simulations in parallel.

`spiroid -b input_directory output_directory`

With input argument:

- `input_directory`: Path to the directory containing one or more input configuration files (`.conf`).

### Output
Launching a simulation produces three output files:

- `input_file.conf`: A copy of the input configuration file used to launch this simulation (JSON).
- `simulation.log`: A log file containing status information for each simulation. Will be appended to if it already exists.
- `input_file.jsonl`: Output for the simulation (JSONL). 

For example:
```bash
spiroid example.conf simulations
```
will produce the following structure:
(A `-n` sub directory will be created for each simulation, with `n` as the smallest non pre-existing numerical suffix.)
```
simulations/
└── simulation.log
└── example
    ├── example.conf
    └── example.jsonl
```

### Analyse a simulation
The output data is available in `/output/simulation_name/simulation_name.jsonl`

#### Plotting
A simple plotting script is provided to view output quantities over time.

To create merged plots to compare output from one or more simulations:

`python3 scripts/plotter.py input_data [input_data1] [input_data2] [...] output_directory`

With input parameters:

- `input_data*`: Path(s) to spiroid output file(s). (e.g. `/path/to/simulations/simulation_name/simulation_name.jsonl`).
- `output_directory`: Path to the desired output destination. Will be automatically created if it does not already exist.

Example:

`python3 scripts/plotter.py simulation/{simulation_name,other_simulation_name}/*.jsonl plots`

Will create merged plots for each quantity inside the `plots` directory.

It is also possible to plot all simulations independently:

`python3 scripts/plotter.py input_directory`

With input parameters:
- `input_directory`: Path to the spiroid simulation output directory.

Example: 
`python3 scripts/plotter.py simulations`

Will create plots for all quantities inside the respective `simulations/simulation_name/` directories.

> Note: Be mindful when plotting large simulations, or a large number of simulations, as python is slow.
 
## Modes of operation
The following modes of operation can be toggled independently, in any combination:

- Evolution of the `Star`
- Magnetic interaction between the `Star` and the `Planet`
- Constant Time Lag tidal interaction of the `Planet` on the `Star`
- Kaula tidal interaction of the `Star` on the `Planet`

> Note: Operation with both stellar evolution and planetary tides (kaula) disabled is untested.

### Stellar evolution (`Star`)
Evolution of the star can be enabled by setting the `evolution` property of the `Star` to `Interpolated` and provide the appropriate data file.
(e.g. `examples/data/star/evolution/savgol_08.csv`).
The format of the file must be CSV (Comma Separated Values) with the following header and fields: 

- `age` (years)
- `radius` (solar radius)
- `mass` (solar mass)
- `convective_radius` (solar radius)
- `convective_mass` (solar mass)
- `radiative_moment_of_inertia`  (stellar mass * stellar radius^2)
- `convective_moment_of_inertia` (stellar mass * stellar radius^2)
- `luminosity` (solar luminosity)

### Magnetism (`Particle`)
Magnetic interaction (initiated by the star) can be toggled into the following states:

- Disabled
- Wind

#### Wind
Set the `magnetism` property of the `central_body` `Particle` to `Wind`.
Set the desired values of `magnetic_field` for `Planet` and `Star`.

### Tides

Tides for each particle can be toggled independently into the following states:

- Disabled
- Constant Time Lag
- Kaula Tides

> Note: Currently CTL is not implemented for the Planet, and Kaula is not implemented for the Star.

#### Constant Time Lag (star)
Set the `tides` property of the `central_body` particle to `ConstantTimeLag` and provide the `equilibrium_tide_dissipation` factor. 

#### Kaula tides (planet)
Set the `tides` property of the `orbiting_body` particle to `KaulaTides`, specify the `particle_type` (e.g. `Solid`) and provide the appropriate love number data file.
(e.g. `examples/data/planet/tides/kaula/leconte2015_steinberger.csv`).
The format of the file must be CSV (Comma Separated Values) with the following header and fields:

- `tidal_frequency`
- `imaginary_love_number`
- `real_love_number`

## Testing
Run the built in tests:
``` cargo test```

Create a new test case when functionality is added or modified. Use `cargo tarpaulin` to ensure the test coverage touches the additional functionality.

> Note: Tests are currently tailored for `amd64` and may fail (with slight numerical discrepancies) on other platforms due to differences in cpu architecture.

## Design and Structure

Spiroid uses the [simulation](https://github.com/DynaClim/simulation) crate for creating a numerical integration simulator. At each timestep the properties of the particles (e.g. planet and star) are updated prior to calculating the derivatives for the integration quantities. The general sequence of operation is:

1. Initial conditions are parsed from the input config, applying unit conversions. Any auxiliary data files are also read.
1. The simulation is launched using initial integration quantities.
1. At each timestep of the simulation, the integrator calls the `derive` function, to solve the ODEs.
1. Before the ODEs are solved, the state of the universe is updated based on the current values of the integrator. This calculates the current star properties, current planet properties and, if appropriate, magnetic and tidal interactions of the star, and the kaula values for the planetary tide.
1. The ODEs are solved and the output is returned to the integrator.
1. If the integrator accepts the integration step, it calls the `solout` function, which appends the current integration values to the output file.
1. The simulation terminates after the maximum duration is reached, maximum number of iterations are reached, or a numerical error occurs (attempt to interpolate out of range, computation of infinity or NaN, etc).

### File structure
```
src
├── lib.rs
├── main.rs
├── physics.rs
├── universe
│   ├── effects
│   │   ├── magnetism.rs
│   │   ├── tides
│   │   │   ├── constant_time_lag.rs
│   │   │   ├── kaula
│   │   │   │   ├── love_number.rs
│   │   │   │   └── polynomials.rs
│   │   │   └── kaula.rs
│   │   └── tides.rs
│   ├── effects.rs
│   ├── particles
│   │   ├── planet.rs
│   │   ├── star
│   │   │   └── star_csv.rs
│   │   └── star.rs
│   └── particles.rs
├── universe.rs
└── utils.rs
```

- lib.rs: Specifies the structure of the `Universe` for the config file and implements the `Integrator::System` trait required by the integrator.
- main.rs: Robust example of using the spiroid `lib` crate.
- physics.rs: `derive` function called at every iteration by the integrator (`system.derive()`) to compute the derivatives for the integration quantities.
- universe.rs: Public interface for the system containing particles.
- universe/effects.rs: Public interface for effects (`Tides` and `Magnetism`).
- universe/effects/magnetism.rs: Magnetic SPI (physics).
- universe/effects/tides.rs: Tidal SPI (physics).
- universe/particles.rs: Public interface for particles (`Star` and `Planet`) and generic `Particle` types.
- universe/particles/planet.rs: `Planet` model. 
- universe/particles/star.rs: `Star` model. 
- universe/particles/star/star_csv.rs: Structure of the stellar evolution from `STAREVOL` models. 
- utils.rs: Misc helper functions, precomputed factorial, etc.

# Known issues and limitations

- Due to issues with the integrator, certain simulations can become stuck in a seeminly endless loop when the timestep becomes too small. These have been observed with simulations enabling Kaula tides on the planet.
- Data is output at every accepted timestep. Files can become quite large for long running simulations.
- Simulations with both stellar evolution and planetary tides (kaula) disabled are untested.

# Notes

Development is funded by the Swiss National Science Foundation.
