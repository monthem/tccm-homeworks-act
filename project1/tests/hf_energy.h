#ifndef HF_ENERGY_H
#define HF_ENERGY_H

#include <stdio.h>
#include <stdlib.h>
#include "TREXIO_functions.h"

// Macro for 4D indexing in a 1D array
#define INDEX4D(i, j, k, l, dim) ((i)(dim)(dim)(dim) + (j)(dim)(dim) + (k)(dim) + (l))

// Function prototype for computing Hartree-Fock energy
double compute_HF_energy(double E_NN, double* one_e_integrals, int64_t n_integrals,
                         int32_t* index, double* value, int mo_num, int n_occ);

#endif // HF_ENERGY_H