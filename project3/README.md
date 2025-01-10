README.md
=========

# Molecular Dynamics Simulation Program

## Description
This program simulates the dynamics of atoms in a system using the Lennard-Jones (LJ) potential. It features modules for reading input, computing interatomic distances, calculating kinetic and potential energy, enforcing boundary conditions, and running the dynamics simulation.

The program reads atomic data from an input file, performs the simulation, and outputs a trajectory file in `.xyz` format for visualization.

## Directory Structure
- `/src`: Contains all source files.
  - `md_io.f90`: Handles reading input and writing the trajectory output file.
  - `md_math.f90`: Computes interatomic distances.
  - `md_potential.f90`: Calculates total LJ potential energy.
  - `md_kinetics.f90`: Computes total kinetic energy, acceleration matrix, and enforces boundary conditions.
  - `md_simulation.f90`: Combines all modules into a subroutine that runs the simulation and writes the trajectory output file.
  - `main.f90`: Main program file that initializes variables and calls the simulation routines.

## Compilation
The compilation instructions assume a Linux-based system with `gfortran` installed.

### Steps:
1. Navigate to the `/src` directory.
2. Run the following commands:

```bash
gfortran -c md_io.f90 md_math.f90 md_potential.f90 md_kinetics.f90 md_simulation.f90 main.f90
gfortran md_io.o md_math.o md_potential.o md_kinetics.o md_simulation.o main.o -o md_program
```

Alternatively, use the `compileall.sh` script located in the root directory to automate the process.

## Usage
1. Copy the compiled executable (`md_program`) to a directory containing an input file.
2. Run the program:
   ```bash
   ./md_program
   ```
3. The program will read the input file, perform the simulation, and write the output trajectory to a file (default: `trajectory.xyz`).

## Input File Format
Input files should follow this format:
- First line: Number of atoms (`nat`).
- Subsequent lines: Atomic data in the format `x y z mass`, where:
  - `x`, `y`, `z` are the coordinates of the atom in nm.
  - `mass` is the atomic mass in amu.

## Configurable Parameters
The following parameters are hardcoded in `main.f90` and can be changed by the user:
- `input` (string): Input file name (default: `inp.txt`).
- `output` (string): Output trajectory file name (default: `trajectory.xyz`).
- `atom` (string): Atom name (default: `Ar`).
- `epsilon`, `sigma` (real): LJ potential parameters (default values for Ar in J/mol and nm, respectively).
- `timestep` (real): Simulation time step in ps (recommended: `sqrt(m) * 0.0005`, where `m` is mass in amu).
- `maxstep` (integer): Maximum number of simulation steps.
- `Lx`, `Ly`, `Lz` (real): Box dimensions in nm (default: 8x8x8 nm³).

### Notes on Boundary Conditions
The boundary conditions enforce a hard wall from `-4` to `L` nm in all three directions. Ensure that the initial geometry fits within these bounds.

## Recommendations
- Recompile the program using `compileall.sh` after any changes to `main.f90` parameters.
- Avoid overwriting existing trajectory files by specifying a new output file name for each simulation.
- Use tools like ASE for visualization:
  ```bash
  ase gui trajectory.xyz
  ```

## Example Inputs
The `example_inputs` directory contains sample input files for 3, 20, and 50 Ar atom systems.

## Tests
The `tests` directory includes tests for individual modules. To run tests, ensure an `inp.txt` file is present in the directory and execute the corresponding test scripts.