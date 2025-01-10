gfortran -c md_io.f90 md_math.f90 md_potential.f90 md_kinetics.f90 md_simulation.f90 main.f90
gfortran md_io.o md_math.o md_potential.o md_kinetics.o md_simulation.o main.o -o main
