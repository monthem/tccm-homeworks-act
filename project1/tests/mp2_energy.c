#include <stdio.h>
#include <stdlib.h>
#include "TREXIO_functions.h"
#include "mp2_energy.h"


// A helper macro to do 4D indexing in a 1D array
#define MO_TEI(i,j,k,l, arr, mo_num) \
    ((arr[ (((size_t)(i) * (mo_num) + (j)) * (mo_num) + (k)) * (mo_num) + (l) ]))


double compute_MP2_energy(int64_t n_integrals, int32_t* index, double* value, int mo_num, int n_occ, double* orbital_energies) {

    double MP2_energy = 0.0;

    size_t total_size = (size_t) mo_num * mo_num * mo_num * mo_num;
    double* mo_tei = (double*) calloc(total_size, sizeof(double));
    if (mo_tei == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for mo_tei.\n");
        exit(1);
    }

    for (int n = 0; n < n_integrals; n++) {
        int i = index[4*n + 0];
        int j = index[4*n + 1];
        int k = index[4*n + 2];
        int l = index[4*n + 3];

        // Filter only MP2 relevant integrals < a b | i j >
        if (i < n_occ || j < n_occ || k >= n_occ || l >= n_occ) {
        continue; 
        }

        double integral = value[n];

        // Store all 8 index variations that are equal 
        MO_TEI(i,j,k,l, mo_tei, mo_num) = integral;
        MO_TEI(i,l,k,j, mo_tei, mo_num) = integral;
        MO_TEI(k,l,i,j, mo_tei, mo_num) = integral;
        MO_TEI(k,j,i,l, mo_tei, mo_num) = integral;

        MO_TEI(j,i,l,k, mo_tei, mo_num) = integral;
        MO_TEI(l,i,j,k, mo_tei, mo_num) = integral;
        MO_TEI(l,k,j,i, mo_tei, mo_num) = integral;
        MO_TEI(j,k,l,i, mo_tei, mo_num) = integral;
    }

    for (int i = 0; i < n_occ; i++) { // occupied orbitals
        for (int j = 0; j < n_occ; j++) { // occupied orbitals
            for (int a = n_occ; a < mo_num; a++) { // virtual orbitals
                for (int b = n_occ; b < mo_num; b++) { // virtual orbitals
                    double coulomb = MO_TEI(a, b, i, j, mo_tei, mo_num);
                    double exchange = MO_TEI(b, a, i, j, mo_tei, mo_num);
                    double denominator = orbital_energies[i] + orbital_energies[j] - orbital_energies[a] - orbital_energies[b];
                    double numerator = coulomb * (2*coulomb - exchange);
                    MP2_energy += numerator / denominator;
                }
            }
        }
    }
    free(mo_tei);
    return MP2_energy;
}