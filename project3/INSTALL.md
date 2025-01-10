INSTALL.md
===========

# Installation Instructions

## Requirements
The program requires the following software:
- `gfortran` (GNU Fortran compiler)
- Linux-based operating system is recommended for these instructions.

## Compilation

1. Navigate to the `/src` directory.
2. Compile the program by running:
   ```bash
   gfortran -c md_io.f90 md_math.f90 md_potential.f90 md_kinetics.f90 md_simulation.f90 main.f90
   gfortran md_io.o md_math.o md_potential.o md_kinetics.o md_simulation.o main.o -o md_program
   ```

Alternatively, use the `compileall.sh` script in the root directory:
   ```bash
   ./compileall.sh
   ```

## Post-Compilation
- Copy the compiled executable (`md_program`) to a folder containing the input file.
- Example inputs are provided in the `example_inputs` directory.
- To recompile after changing `main.f90`, rerun the compilation steps or use the `compileall.sh` script.

The program is now ready to use!