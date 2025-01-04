#ifndef HF_ENERGY_H
#define HF_ENERGY_H

#include <stdio.h>
#include <stdlib.h>
#include "TREXIO_functions.h"

// Macro for 4D indexing in a 1D array
#define INDEX4D(i, j, k, l, dim) ((i)(dim)(dim)(dim) + (j)(dim)(dim) + (k)(dim) + (l))

// Function prototype for computing Hartree-Fock energy
double compute_HF_energy(double* one_electron_integrals, double* two_electron_integrals, 
                         int num_occupied_orbitals, int num_molecular_orbitals, double nuclear_repulsion_energy);

#endif // HF_ENERGY_H