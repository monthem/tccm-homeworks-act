module md_simulation
        use md_io
        use md_math
        use md_potential
        use md_kinetics

        implicit none
        contains
        
        subroutine run_dynamics(traj_filename, atom, natoms, masses, coords, epsilon, sigma, Lx, Ly, Lz, timestep, maxstep)
                ! This subroutine implements the Verlet algorithm
                ! This algorithm lets atomic positions evolve over time 
                ! Coordinates are updated every timestep current_time
                ! But are printed every 100 steps

                ! Arguments
                character(len=*), intent(in) :: traj_filename
                character(2), intent(in) :: atom
                integer, intent(in) :: natoms, maxstep
                double precision, intent(in) :: epsilon, sigma, timestep
                double precision, intent(in) :: Lx, Ly, Lz
                double precision, intent(in) :: masses(natoms)
                double precision, intent(inout) :: coords(natoms, 3)

                ! Local variables
                integer :: n
                integer :: i
                double precision :: E, T, V, current_time
                double precision, allocatable :: distances(:,:)
                double precision, allocatable :: velocities(:,:)
                double precision, allocatable :: acceleration_matrix(:,:)
                allocate(distances(natoms, natoms))
                allocate(velocities(natoms, 3))
                allocate(acceleration_matrix(natoms, 3))
                
                current_time = 0.0
                velocities = 0.0d0 ! initialize velocities to 0
                call compute_interatomic_distances(natoms, coords, distances) ! compute initial interatomic distances
                V = compute_lj_potential(natoms, distances, epsilon, sigma) ! calculate initial LJ potential energy
                T = compute_total_ke(natoms, masses, velocities)
                E = T + V
                call compute_acceleration_vectors(natoms, epsilon, sigma, masses, coords, distances, acceleration_matrix) ! initial
                coords(:,:) = coords(:,:) + 4.00d0
                call write_xyz_trajectory(traj_filename, current_time, atom, natoms, coords, E, T, V)
                coords(:,:) = coords(:,:) - 4.00d0

                do n = 1, maxstep
                        current_time = current_time + timestep

                        coords = coords + timestep * velocities + 0.50d0 * timestep**2 * acceleration_matrix ! update coords
         
                        call compute_interatomic_distances(natoms, coords, distances) ! update distances
                        call enforce_boundary_conditions(natoms, coords, velocities, Lx, Ly, Lz)
                        velocities = velocities + 0.50d0 * timestep * acceleration_matrix ! first update to v
                        call compute_acceleration_vectors(natoms, epsilon, sigma, masses, coords, distances, acceleration_matrix) !
                        velocities = velocities + 0.50d0 * timestep * acceleration_matrix
                        
                        ! update energetic terms
                        V = compute_lj_potential(natoms, distances, epsilon, sigma)
                        T = compute_total_ke(natoms, masses,  velocities)
                        E = T + V

                        ! update trajectory file every 10 steps
                        if (mod(n, 10) == 0) then
                                coords(:,:) = coords(:,:) + 4.00d0
                                call write_xyz_trajectory(traj_filename, current_time, atom, natoms, coords, E, T, V)
                                coords(:,:) = coords(:,:) - 4.00d0
                        end if
                end do
        
        deallocate(distances, velocities, acceleration_matrix)

        end subroutine run_dynamics
end module md_simulation

                
                  

                
                

