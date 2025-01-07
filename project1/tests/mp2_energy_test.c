#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mp2_energy.h"
#include "TREXIO_functions.h"

void test_compute_MP2_energy() {

    // tests the MP2 energy function on the h2o.h5 file

    const char* filename = "h2o.h5";

    trexio_t* trexio_file = open_trexio_file(filename);

    int mo_num = read_mo_num(trexio_file);
    int n_occ = get_number_of_occupied_orbitals(trexio_file);
    int64_t n_integrals;
    int32_t* index;
    double* value; 
    double* one_e_integrals = read_one_electron_integrals(trexio_file, mo_num);
    read_two_electron_integrals(trexio_file, &n_integrals, &index, &value);
    double* orbital_energies = read_orbital_energies(trexio_file, mo_num);
    double MP2_energy = compute_MP2_energy(n_integrals, index, value, mo_num, n_occ, orbital_energies);
    // Validate the result
    double diff = fabs(MP2_energy + 0.20395997);
    assert(diff < 1e-7);
    printf("MP2 energy test passed\n");
}

int main() {
    test_compute_MP2_energy();
    return 0;
}