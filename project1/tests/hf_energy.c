#include <stdio.h>
#include <stdlib.h>
#include "TREXIO_functions.h"
#include "hf_energy.h"

// A helper macro to do 4D indexing in a 1D array
#define MO_TEI(i,j,k,l, arr, mo_num) \
    (arr[ (((size_t)(i) * (mo_num) + (j)) * (mo_num) + (k)) * (mo_num) + (l) ])

double compute_HF_energy(double E_NN, double* one_e_integrals, int64_t n_integrals, int32_t* index, double* value, int mo_num, int n_occ) {
    
    double hf_energy = E_NN;
    double sum_one_e_i = 0.0;

    // Sum the one-electron contributions
    for (int i = 0; i < n_occ; i++) {
        sum_one_e_i += 2 * one_e_integrals[i * mo_num + i];
    }

    size_t total_size = (size_t) mo_num * mo_num * mo_num * mo_num;
    double* mo_tei = (double*) calloc(total_size, sizeof(double));
    if (mo_tei == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for mo_tei.\n");
        exit(1);
    }

    // Retrieve all stored 2e integrals
    for (int n = 0; n < n_integrals; n++) {
        int i = index[4*n + 0];
        int j = index[4*n + 1];
        int k = index[4*n + 2];
        int l = index[4*n + 3];


        // filter only HF relevant integrals i,j,k,l element [0, n_occ)
        if (i >= n_occ || j >= n_occ || k >= n_occ || l >= n_occ) {
        continue; // Skips virtual orbitals
        }

        double integral = value[n];

        // Store all 8 equal variations of (ij|kl)
        MO_TEI(i,j,k,l, mo_tei, mo_num) = integral;
        MO_TEI(i,l,k,j, mo_tei, mo_num) = integral;
        MO_TEI(k,l,i,j, mo_tei, mo_num) = integral;
        MO_TEI(k,j,i,l, mo_tei, mo_num) = integral;

        MO_TEI(j,i,l,k, mo_tei, mo_num) = integral;
        MO_TEI(l,i,j,k, mo_tei, mo_num) = integral;
        MO_TEI(l,k,j,i, mo_tei, mo_num) = integral;
        MO_TEI(j,k,l,i, mo_tei, mo_num) = integral;
    }

    // Sum all 2e contributions
    double two_e_sum = 0.0;
    for (int i = 0; i < n_occ; i++) {
        for (int j = i; j < n_occ; j++) {
            double coulomb  = MO_TEI(i, j, i, j, mo_tei, mo_num);
            double exchange = MO_TEI(i, j, j, i, mo_tei, mo_num);
        if (i == j) {
            two_e_sum += (2.0 * coulomb - exchange);
        } else {
            two_e_sum += 2.0 * (2.0 * coulomb - exchange);
        }
        }   
    }

    // Final summation
    hf_energy += sum_one_e_i + two_e_sum;

    // Free allocated memory
    free(mo_tei);
    
    // Return the final HF energy
    return hf_energy;
}