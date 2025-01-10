program test_md_simulation
        use md_simulation
        implicit none

        integer :: i, natoms, maxstep
        character(len=2) :: atom
        double precision, allocatable :: coords(:,:), masses(:)
        double precision :: epsilon, sigma, timestep, Lx, Ly, Lz
        call read_input_file("inp.txt", natoms, coords, masses)        
        epsilon = 0.06610d0 ! J / mol
        sigma = 0.33450d0 ! nm
        Lx = 10.00d0
        Ly = 10.00d0
        Lz = 10.00d0
        atom = "Ar"
        timestep = 0.003d0
        maxstep = 100000
        call run_dynamics("test_dynamics2.xyz", atom, natoms, masses, coords, epsilon, sigma, Lx, Ly, Lz, timestep, maxstep)
        
end program test_md_simulation
