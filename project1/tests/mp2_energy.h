#ifndef MP2_ENERGY_H
#define MP2_ENERGY_H

#include <stdio.h>
#include <stdlib.h>
#include "TREXIO_functions.h"

// Macro for 4D indexing in a 1D array
#define INDEX4D(i, j, k, l, dim) ((i)(dim)(dim)(dim) + (j)(dim)(dim) + (k)(dim) + (l))

// Function prototype for computing MP2 energy
double compute_MP2_energy(int64_t n_integrals, int32_t* index, double* value, int mo_num, int n_occ, double* orbital_energies);

#endif // MP2_ENERGY_H