program test_md_potential
        use md_io
        use md_math
        use md_potential
        implicit none

        integer :: i, natoms
        double precision, allocatable :: coords(:,:), masses(:)
        double precision, allocatable :: distances(:,:)
        double precision :: epsilon, sigma, v_total

        call read_input_file("inp.txt", natoms, coords, masses)

        allocate(distances(natoms, natoms))
        
        epsilon = 0.0661d0
        sigma = 0.3345d0
        call compute_interatomic_distances(natoms, coords, distances)
        v_total = compute_lj_potential(natoms, distances, epsilon, sigma)
        print *, "Total Lennard-Jones potential energy: ", v_total

        ! Memory clean up
        deallocate(coords, masses, distances)
end program test_md_potential
