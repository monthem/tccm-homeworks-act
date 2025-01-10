INSTALL.md
===========

# Installation Instructions

## Requirements
The program requires the following software installed on your Linux-based system:
- GCC compiler
- TREXIO library (version 2.5.0 or later)
- HDF5 library

## Installing the HDF5 Library
For Ubuntu:
```bash
sudo apt install libhdf5-dev
```

For Arch Linux:
```bash
sudo pacman -S hdf5
```

For macOS:
```bash
brew install hdf5
```

## Installing TREXIO Library

1. Download the TREXIO source code:
   ```bash
   wget https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz
   ```

2. Extract the files:
   ```bash
   tar -zxvf trexio-2.5.0.tar.gz
   cd trexio-2.5.0
   ```

3. Build and install:
   ```bash
   ./configure
   make
   sudo make install
   ```

   The library will be installed in `/usr/local/lib`, and header files will be in `/usr/local/include`.

4. Verify the installation:
   ```bash
   ls /usr/local/lib | grep trexio
   ```

## Compiling the Program
1. Navigate to the `/src` directory where the source files are located.

2. Compile the program:
   ```bash
   gcc -o electronic_energy main.c TREXIO_functions.c hf_energy.c mp2_energy.c -I/usr/local/include -L/usr/local/lib -ltrexio -lm
   ```

3. Copy the compiled executable to the directory containing the TREXIO `.h5` files:
   ```bash
   cp electronic_energy /path/to/trexio/files/
   ```

The program is now ready to use!
