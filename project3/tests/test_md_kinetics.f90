program test_md_kinetics
        use md_io
        use md_math
        use md_potential
        use md_kinetics
        implicit none

        integer :: i, natoms
        double precision, allocatable :: coords(:,:), masses(:)
        double precision, allocatable :: distances(:,:)
        double precision :: epsilon, sigma, v_total
        double precision, allocatable :: velocities(:,:)
        double precision :: T_total
        double precision :: E_total
        double precision, allocatable :: acceleration_matrix(:,:)

        allocate(distances(natoms, natoms))
        call read_input_file("inp.txt", natoms, coords, masses)        
        epsilon = 0.0661
        sigma = 0.3345
        call compute_interatomic_distances(natoms, coords, distances)
        v_total = compute_lj_potential(natoms, distances, epsilon, sigma)
        print *, "Total Lennard-Jones potential energy: ", v_total
        
        ! kinetics now        
        allocate(velocities(natoms, 3))
        velocities = 0.0
        
        T_total = compute_total_ke(natoms, masses, velocities)
        
        E_total = v_total + T_total
        print *, "Total system energy = ", E_total
        print *, "Kinetic energy = ", T_total, " and LJ potential energy = ", v_total
                
        allocate(acceleration_matrix(natoms, 3))
        call compute_acceleration_vectors(natoms, epsilon, sigma, masses, coords, distances, acceleration_matrix)
        print *, "Acceleration matrix:"
        print *, acceleration_matrix
        
        ! Memory clean up
        deallocate(coords, masses, distances, velocities, acceleration_matrix)
end program test_md_kinetics
