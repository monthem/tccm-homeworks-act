program small_molecular_dynamics
        
        use md_simulation
        implicit none

        ! Variables to initialize the system from the input file
        integer :: natoms
        double precision, allocatable :: coords(:,:), masses(:)
        ! Variables needed to run the dynamics
        integer :: maxstep
        double precision :: epsilon, sigma, timestep, Lx, Ly, Lz
        ! Variables needed for file I/O
        character(len=256) :: input, output
        character(len=2) :: atom
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! USER SPECIFIC SECTION, EDIT THE VALUES ACCORDINGLY (WITH CARE!)
        input = "inp.txt" ! Replace with wanted input file if necessary
        output = "trajectory.xyz" ! Replace with wanted output file if necessary
        ! System specific values - DEFAULT SYSTEM: Ar atoms
        atom = "Ar"
        epsilon = 0.06610d0 ! If changing atom type, replace with appropriate value in J/mol
        sigma = 0.33450d0 ! If changing atom type, replace with appropraite value in nm
        timestep = 0.003d0 ! Recommended value: sqrt(mass in amu) * 0.0005 ps
        maxstep = 10000 ! Change freely. Trajectory is updated every 10 steps
        ! Hard wall cell dimensions. Actual dimensions are + 4 nm in each direction
        Lx = 4.00d0
        Ly = 4.00d0
        Lz = 4.00d0
        ! END OF USER SPECIFIC SECTION. HAVE FUN!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call read_input_file(input, natoms, coords, masses)
        call run_dynamics(output, atom, natoms, masses, coords, epsilon, sigma, Lx, Ly, Lz, timestep, maxstep)
        
end program small_molecular_dynamics
