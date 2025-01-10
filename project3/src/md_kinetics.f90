module md_kinetics
        implicit none
        contains

        function compute_total_ke(natoms, masses, velocities) result(T)
                ! Calculates total kinetic energy of the system
                ! T = sum over i (0.5 * m_i * v_i^2)
                ! Arguments
                integer, intent(in) :: natoms
                double precision, intent(in) :: masses(natoms)
                double precision, intent(in) :: velocities(natoms, 3)

                ! Local variables
                integer :: i

                ! Output
                double precision :: T

                T = 0.5d0 * sum(masses * (velocities(:, 1)**2 + velocities(:,2)**2 + velocities(:,3)**2))
        end function compute_total_ke

        subroutine compute_acceleration_vectors(natoms, epsilon, sigma, masses, coords, distances, acceleration_matrix)
                ! Calculates acceleration vectors for each atom

                ! Arguments
                integer, intent(in) :: natoms
                double precision, intent(in) :: sigma, epsilon
                double precision, intent(in) :: masses(natoms), coords(natoms, 3)
                double precision, intent(in) :: distances(natoms, natoms)
                double precision, intent(out) :: acceleration_matrix(natoms, 3)

                ! Local variables
                integer :: i, j
                double precision :: force_x, force_y, force_z, r_inv
                double precision :: mass_inv, U, r, norm_comp_x, norm_comp_y, norm_comp_z
                double precision, parameter :: r_min = 1.0d-5
                
                ! Initializion of the acceleration matrix
                acceleration_matrix = 0.0d0

                do i = 1, natoms
                        do j = i + 1, natoms
                                r = distances(i, j)
                                if (r > r_min) then
                                        r_inv = 1 / r                                        
                                        U = (24.0d0 * epsilon * r_inv ) * ((sigma * r_inv)**6 - 2.0d0 * (sigma * r_inv)**12)
                                        force_x = U * (coords(i, 1) - coords(j, 1)) * r_inv
                                        force_y = U * (coords(i, 2) - coords(j, 2)) * r_inv
                                        force_z = U * (coords(i, 3) - coords(j, 3)) * r_inv
                                
                                        ! Update acceleration for atom i
                                        acceleration_matrix(i, 1) = acceleration_matrix(i, 1) + force_x
                                        acceleration_matrix(i, 2) = acceleration_matrix(i, 2) + force_y
                                        acceleration_matrix(i, 3) = acceleration_matrix(i, 3) + force_z

                               
                                        ! Update acceleration for atom j (Newton's third law)
                                        acceleration_matrix(j, 1) = acceleration_matrix(j, 1) - force_x
                                        acceleration_matrix(j, 2) = acceleration_matrix(j, 2) - force_y
                                        acceleration_matrix(j, 3) = acceleration_matrix(j, 3) - force_z
                                end if
                        end do 
                end do
                ! Conversion of accumulated forces to accelerations 
                do i = 1, natoms
                        acceleration_matrix(i, :) = acceleration_matrix(i, :) / masses(i)
                end do
        end subroutine compute_acceleration_vectors

        subroutine enforce_boundary_conditions(natoms, coords, velocities, Lx, Ly, Lz)
                ! Enforces hard wall boundary conditions with elastic collisions 
                
                ! Arguments
                integer, intent(in) :: natoms
                double precision, intent(inout) :: coords(natoms, 3), velocities(natoms, 3)
                double precision, intent(in) :: Lx, Ly, Lz
                ! Local variables
                integer :: i, d
                double precision, parameter :: lower_bound = -4.0d0 ! IN NANOMETERS
                double precision, parameter :: reflection_factor = 2.0d0
                double precision :: box(3)

                ! Define the box dimensions
                box(1) = Lx
                box(2) = Ly
                box(3) = Lz 

                do i = 1, natoms
                        do d = 1, 3
                                if (coords(i, d) < lower_bound) then
                                        coordS(i, d) = -coords(i, d)
                                        velocities(i, d) = -velocities(i, d)
                                else if (coords(i, d) > box(d)) then
                                        coords(i, d) = reflection_factor * box(d) - coords(i, d)
                                        velocities(i, d) = -velocities(i, d)
                                end if
                        end do
                end do
        end subroutine enforce_boundary_conditions

end module md_kinetics


                                
                                 





                         
                        
                        
