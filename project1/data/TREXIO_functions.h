#ifndef TREXIO_FUNCTIONS_H   // Include guard: Start
#define TREXIO_FUNCTIONS_H

#include <trexio.h>         // Include TREXIO library for function prototypes

// Function to open a TREXIO file
trexio_t* open_trexio_file(const char* filename);

// Function to close a TREXIO file
void close_trexio_file(trexio_t* trexio_file);

// Function to read nuclear repulsion energy
double read_nuclear_repulsion(trexio_t* trexio_file);

// Function to get the number of occupied orbitals
int get_number_of_occupied_orbitals(trexio_t* trexio_file);

// Function to get the number of molecular orbitals
int read_mo_num(trexio_t* trexio_file);

// Function to read one-electron integrals
double* read_one_electron_integrals(trexio_t* trexio_file, int mo_num);

// Function to read two-electron integrals
void read_two_electron_integrals(trexio_t* trexio_file, int64_t* n_integrals, int32_t** index, double** value);

#endif // TREXIO_FUNCTIONS_H  // Include guard: End
