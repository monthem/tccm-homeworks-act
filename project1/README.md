README.md
=========

# Electronic Energy Calculation Program

## Description
This program computes the Hartree-Fock (HF) energy and MP2 correlation energy for closed-shell systems using integrals and orbital energies stored in a TREXIO file. The total energy, including the MP2 correction, is also computed and displayed.

The program reads the necessary data from a TREXIO `.h5` file and uses:
- `hf_energy.c` to compute HF energy
- `mp2_energy.c` to compute MP2 energy
- Wrapped TREXIO functions in `TREXIO_functions.c`

## Directory Structure
- `/src`: Contains all source files and headers.
  - `main.c`: Main program
  - `TREXIO_functions.c`: TREXIO function wrappers
  - `hf_energy.c`: HF energy computation
  - `mp2_energy.c`: MP2 energy computation
  - Corresponding `.h` files for linking
- `/tests` : Contains testing files from program development
- `/data` : Contains example molecule.h5 TREXIO files the user can test the program on

## Compilation
The compilation instructions assume a Linux-based system. To compile, run:

```bash
cd /src
gcc -o electronic_energy main.c TREXIO_functions.c hf_energy.c mp2_energy.c -I/usr/local/include -L/usr/local/lib -ltrexio -lm
```

This will create an executable named `electronic_energy`.

## Usage
1. Move the compiled executable to the directory containing the TREXIO `.h5` files. For example:
   ```bash
   cp electronic_energy /path/to/trexio/files/
   cd /path/to/trexio/files/
   ```

2. Run the program:
   ```bash
   ./electronic_energy
   ```

3. The program will prompt for a file name. Enter the name of the TREXIO `.h5` file (e.g., `h2o.h5`) and press enter.

4. The program will compute and display:
   - HF energy
   - MP2 energy correction
   - Total energy (HF + MP2 correction)

## Example
```bash
$ ./electronic_energy
Enter TREXIO file name: h2o.h5
Hartree-Fock energy: -76.026798708 a.u.
MP2 energy correction: -0.203959974 a.u.
Total energy: -76.230758682 a.u.
```

## Requirements
- GCC compiler
- TREXIO library
- HDF5 library

Refer to the `INSTALL.md` for detailed installation instructions.