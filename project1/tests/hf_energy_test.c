#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "hf_energy.h"
#include "TREXIO_functions.h"

void test_compute_HF_energy() {

    // tests the HF energy function on the h2o.h5 file

    const char* filename = "h2o.h5";

    trexio_t* trexio_file = open_trexio_file(filename);

    int mo_num = read_mo_num(trexio_file);
    int n_occ = get_number_of_occupied_orbitals(trexio_file);
    double E_NN = read_nuclear_repulsion(trexio_file);
    int64_t n_integrals;
    int32_t* index;
    double* value; 
    double* one_e_integrals = read_one_electron_integrals(trexio_file, mo_num);
    read_two_electron_integrals(trexio_file, &n_integrals, &index, &value);

    double hf_energy = compute_HF_energy(E_NN, one_e_integrals, n_integrals, index, value, mo_num, n_occ);
    // Validate the result
    double diff = fabs(hf_energy + 76.0267987);
    assert(diff < 1e-7);
    printf("HF energy test passed\n");
}

int main() {
    test_compute_HF_energy();
    return 0;
}