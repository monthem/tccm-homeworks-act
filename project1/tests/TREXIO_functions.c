#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>
#include "TREXIO_functions.h"

// All the functions bellow were provided by the professors
// Here, they are just wrapped for our personal ease of use

//Wrapped function for openning a TREXIO file
trexio_t* open_trexio_file(const char* filename){
    trexio_exit_code rc;
    //open the TREXIO file in a read mode, not write 
    trexio_t* trexio_file = trexio_open(filename, 'r',  TREXIO_AUTO, &rc);

    //To check if the file opened correctly
    if (rc != TREXIO_SUCCESS){
        printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1); // if i see something 1 then it fails. To exit the program if openning fails
    }
    return trexio_file;
    
}

// Wrapped function that closes a trexio file 
void close_trexio_file(trexio_t* trexio_file){
    trexio_exit_code rc = trexio_close(trexio_file);
    // next line to check if the file closed correctly 
    if (rc != TREXIO_SUCCESS){
        printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }
}

//Reads the nuclear repulsion energy from the TREXIO file 
double read_nuclear_repulsion(trexio_t* trexio_file){
    double energy; //to store the nuclear repulsion energy
    trexio_exit_code rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
    //line read if the energy was read correctly 
    if (rc != TREXIO_SUCCESS){
        printf("TREXIO Error reading nuclear repulsion energy: %s\n", trexio_string_of_error(rc));
        exit(1); 
    }
    return energy; //nuclear repulsion energy
}

// Reads the # of occupied orbital =  # of up-spin electrons)
int get_number_of_occupied_orbitals(trexio_t* trexio_file) {
    int32_t n_up; // Variable to store the # of up-spin electrons
    trexio_exit_code rc = trexio_read_electron_up_num(trexio_file, &n_up);
    
    // to check # of electrons was read correctly
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error reading number of occupied orbitals: %s\n", trexio_string_of_error(rc));
        exit(1); 
    }
    
    return n_up; // Return the # of occupied orbitals
}

// Reads the # of MO
int read_mo_num(trexio_t* trexio_file) {
    int32_t mo_num; // Variable to store the # of MO
    trexio_exit_code rc = trexio_read_mo_num(trexio_file, &mo_num);
    
    // to check if the # of orbitals was read successfully
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error reading number of molecular orbitals: %s\n", trexio_string_of_error(rc));
        exit(1); 
    }
    
    return mo_num; // Return the # MO
}

// Reads the one-electron integrals (core Hamiltonian) into a dynamically allocated array
double* read_one_electron_integrals(trexio_t* trexio_file, int mo_num) {
    // Allocate memory for storing the integrals
    double* one_e_integrals = malloc(mo_num * mo_num * sizeof(double));
    if (one_e_integrals == NULL) {
        fprintf(stderr, "Error allocating memory for one-electron integrals\n");
        exit(1);
    }

    // Read the one-electron integrals
    trexio_exit_code rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, one_e_integrals);
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error reading one-electron integrals: %s\n", trexio_string_of_error(rc));
        free(one_e_integrals); // Free memory before exiting
        exit(1);
    }

    return one_e_integrals; // Return the array containing the integrals
}

// Reads two-electron integrals in sparse format
void read_two_electron_integrals(trexio_t* trexio_file, int64_t* n_integrals, int32_t** index, double** value) {
    // Get the #s of non-zero two-electron integrals
    trexio_exit_code rc = trexio_read_mo_2e_int_eri_size(trexio_file, n_integrals);
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error reading two-electron integral size: %s\n", trexio_string_of_error(rc));
        exit(1);
    }

    // Allocate memory for indices and values
    *index = malloc(4 * (*n_integrals) * sizeof(int32_t)); // Indices need 4 integers per integral
    *value = malloc((*n_integrals) * sizeof(double));      // Values are doubles

    if (*index == NULL || *value == NULL) {
        fprintf(stderr, "Error allocating memory for two-electron integrals\n");
        free(*index);
        free(*value);
        exit(1);
    }

    // Read the two-electron integrals
    rc = trexio_read_mo_2e_int_eri(trexio_file, 0, n_integrals, *index, *value);
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error reading two-electron integrals: %s\n", trexio_string_of_error(rc));
        free(*index);
        free(*value);
        exit(1);
    }
}

// Reads orbital energies from the TREXIO file
double* read_orbital_energies(trexio_t* trexio_file, int mo_num) {
    // Allocate memory for orbital energies
    double* orbital_energies = malloc(mo_num * sizeof(double));
    if (orbital_energies == NULL) {
        fprintf(stderr, "Error allocating memory for orbital energies\n");
        exit(1);
    }

    // Read the orbital energies
    trexio_exit_code rc = trexio_read_mo_energy(trexio_file, orbital_energies);
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error reading orbital energies: %s\n", trexio_string_of_error(rc));
        free(orbital_energies); // Free memory before exiting
        exit(1);
    }

    return orbital_energies; // Return the array containing orbital energies
}