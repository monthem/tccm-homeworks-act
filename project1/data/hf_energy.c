#include <stdio.h>
#include <stdlib.h>
#include "TREXIO_functions.h"
#include "hf_energy.h"

// A helper macro to do 4D indexing in a 1D array
#define MO_TEI(i,j,k,l, arr, mo_num) \
    (arr[ (((size_t)(i) * (mo_num) + (j)) * (mo_num) + (k)) * (mo_num) + (l) ])

//start HF part 

// void process_hf_data(const char* filename){
//     trexio_t* trexio_file = open_trexio_file(filename);
    
//     int mo_num = read_mo_num(trexio_file); //to read the # of MO 
//     printf("Number of Molecular Orbitals: %d\n", mo_num);

//     int n_occ = get_number_of_occupied_orbitals(trexio_file);  //to read the # n_occ
//     printf("Number of Occupied Orbitals: %d\n", n_occ); 

//     int n_virt = mo_num - n_occ;
//     printf("Number of Virtual Orbitals: %d\n", n_virt);

//     double E_NN = read_nuclear_repulsion(trexio_file);
//     printf("Nuclear Repulsion Energy: %.6f a.u.\n", E_NN); //to read the nuclear repulsion energy 

//     double* one_e_integrals = read_one_electron_integrals(trexio_file, mo_num); //alocate memory for one-electron integrals 

//     // printf("One-electron Integrals (matrix):\n");
//     // for (int i = 0; i < mo_num; i++){
//     //     for(int j = 0; j < mo_num; j++){
//     //         printf("%10.6f ", one_e_integrals[i * mo_num + j]);

//     //     }
//     //     printf("\n");
//     // }

//     int64_t n_integrals;    // # non zero integrals
//     int32_t* index;         // index array for sparse integrals
//     double* value;          // Value array for sparse integrals
//     read_two_electron_integrals(trexio_file, &n_integrals, &index, &value);

//     free(one_e_integrals); 
//     free(index);
//     free(value);
//     close_trexio_file(trexio_file);
// }


double compute_HF_energy(double E_NN, double* one_e_integrals, int64_t n_integrals, int32_t* index, double* value, int mo_num, int n_occ) {
    
    double hf_energy = E_NN;
    printf("Starting now. HF = %.6f au, just ENN\n", hf_energy);
    printf("number of mos is %d\n",mo_num);
    double sum_one_e_i = 0.0;
    // Add one-electron contributions
    printf("One electronintegralsfor occupied orbitals:\n");
    for (int i = 0; i < n_occ; i++) {
        printf("h[%d][%d] = %.6f\n", i, i, one_e_integrals[i * mo_num + i]);
        sum_one_e_i += 2 * one_e_integrals[i * mo_num + i];
    }

    printf("i computed sum of one electron integrals to be %.6f\n", sum_one_e_i);


    size_t total_size = (size_t) mo_num * mo_num * mo_num * mo_num;
    double* mo_tei = (double*) calloc(total_size, sizeof(double));
    if (mo_tei == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for mo_tei.\n");
        exit(1);
    }

    // 4) Fill in all symmetrical permutations of each integral
    printf("Two-electron integrals now: \n");
    for (int n = 0; n < n_integrals; n++) {
        int i = index[4*n + 0];
        int j = index[4*n + 1];
        int k = index[4*n + 2];
        int l = index[4*n + 3];

        if (i >= n_occ || j >= n_occ || k >= n_occ || l >= n_occ) {
        continue; // Skip virtual orbitals
        }

        double integral = value[n];
        // printf("<%d %d | %d %d> = %.6f\n", i, j, k, l, integral);

        // Place 'val' into all 8 permutations
        MO_TEI(i,j,k,l, mo_tei, mo_num) = integral;
        MO_TEI(i,l,k,j, mo_tei, mo_num) = integral;
        MO_TEI(k,l,i,j, mo_tei, mo_num) = integral;
        MO_TEI(k,j,i,l, mo_tei, mo_num) = integral;

        MO_TEI(j,i,l,k, mo_tei, mo_num) = integral;
        MO_TEI(l,i,j,k, mo_tei, mo_num) = integral;
        MO_TEI(l,k,j,i, mo_tei, mo_num) = integral;
        MO_TEI(j,k,l,i, mo_tei, mo_num) = integral;
    }

    printf("MO_TEI array values:\n");
    for (int i = 0; i < mo_num; i++) {
        for (int j = 0; j < mo_num; j++) {
            for (int k = 0; k < mo_num; k++) {
                for (int l = 0; l < mo_num; l++) {
                    double value = MO_TEI(i, j, k, l, mo_tei, mo_num);
                    if (value != 0.0) {
                        printf("MO_TEI[%d][%d][%d][%d] = %.6f\n", i, j, k, l, value);
                    }
                }
            }
        }
    }


    // 5) Compute the two-electron part of the HF energy
    //    E_2e = 0.5 * sum_{i=1..n_occ} sum_{j=1..n_occ} [2 * <ij|ij> - <ij|ji>]
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

        printf("Occupied pair (%d, %d): Coulomb = %.6f, Exchange = %.6f, Contribution = %.6f\n", 
                i, j, coulomb, exchange, (i == j ? (2.0 * coulomb - exchange) : 2.0 * (2.0 * coulomb - exchange)));

        }

    
    }

    // Final summation
    hf_energy += sum_one_e_i + two_e_sum;

    // 6) Free allocated memory
    free(mo_tei);

    // 7) Return the final HF energy
    return hf_energy;
}





// int main(){
    
//     const char* filename = "h2o.h5";
//     trexio_t* trexio_file = open_trexio_file(filename);

//     double E_NN = read_nuclear_repulsion(trexio_file);
//     int mo_num = read_mo_num(trexio_file);
//     int n_occ = get_number_of_occupied_orbitals(trexio_file);

//     double* one_e_integrals = read_one_electron_integrals(trexio_file, mo_num); //ask chat why does this function take number of mo orbitals as input

//     int64_t n_integrals;
//     int32_t* index;
//     double* value; 
//     read_two_electron_integrals(trexio_file, &n_integrals, &index, &value);

//     double hf_energy = compute_HF_energy(E_NN, one_e_integrals, n_integrals, index, value, mo_num, n_occ);
//     printf("Computed HF Energy: %.6f a.u.\n", hf_energy);

//     free(one_e_integrals);
//     free(index);
//     free(value);
//     close_trexio_file(trexio_file);
//     return 0;

// }