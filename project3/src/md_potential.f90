module md_potential
        implicit none
        contains

        function compute_lj_potential(natoms, distances, epsilon, sigma) result(V)
                ! Calculates the total Lennard-Jones potential energy for the system
                ! V = sum over all r of 4 * eps * [(sigma/r)^12 - (sigma/r)^6]

                ! Arguments
                integer, intent(in) :: natoms 
                double precision, intent(in) :: distances(natoms, natoms)
                double precision, intent(in) :: epsilon, sigma

                ! Output
                double precision :: V

                ! Local variables
                integer :: i, j
                double precision :: r, v_lj
                double precision, parameter :: r_min = 1.0d-5

                V = 0.0d0

                ! Loop over all unique pairs of atoms
                do i = 1, natoms
                        do j = i + 1, natoms
                                r = distances(i, j)
                                if (r > r_min) then ! don't divide by close to 0 distances
                                       v_lj =  4.0d0 * epsilon * ((sigma/r)**12 - (sigma/r)**6)
                                       ! print *, "atoms ", i, " and", j, "and v_lj ij = ", v_lj
                                       V = V + v_lj ! In units of J / mol 
                                end if
                        end do
                end do
        end function compute_lj_potential
end module md_potential
