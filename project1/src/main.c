#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hf_energy.h"
#include "mp2_energy.h"



int main(){
    char filename[256];
    printf("Enter a TREXIO file name or a path to a TREXIO file: ");
    scanf("%255s", filename);
    trexio_t* trexio_file = open_trexio_file(filename);

    int mo_num = read_mo_num(trexio_file);
    int n_occ = get_number_of_occupied_orbitals(trexio_file);
    double E_NN = read_nuclear_repulsion(trexio_file);
    int64_t n_integrals;
    int32_t* index;
    double* value; 
    double* one_e_integrals = read_one_electron_integrals(trexio_file, mo_num);
    read_two_electron_integrals(trexio_file, &n_integrals, &index, &value);
    double* orbital_energies = read_orbital_energies(trexio_file, mo_num);

    double hf_energy = compute_HF_energy(E_NN, one_e_integrals, n_integrals, index, value, mo_num, n_occ);
    double MP2_energy = compute_MP2_energy(n_integrals, index, value, mo_num, n_occ, orbital_energies);
    printf("Computed HF energy:  %.9f a.u.\n", hf_energy);
    printf("Computed MP2 energy: %.9f a.u.\n", MP2_energy);
    printf("HF energy with the MP2 correction: %.9f a.u.\n", hf_energy+MP2_energy);

    free(index);
    free(value);
    free(one_e_integrals);
    free(orbital_energies);
    close_trexio_file(trexio_file);
    return 0;
}