module md_math
        implicit none
        contains

        subroutine compute_interatomic_distances(natoms, coords, distances)
                ! Computes pairwise distances between all atoms

                ! Arguments
                integer, intent(in) :: natoms 
                double precision, intent(in) :: coords(natoms, 3)
                double precision, intent(out) :: distances(natoms, natoms) ! distances_ij = r(atom_i, atom_j)

                ! Local variables
                integer :: i, j
                double precision :: dx, dy, dz

                ! Initialize all elements of distances to be equal zero
                distances = 0.0d0

                ! Compute the distances 
                do i = 1, natoms
                        do j = i + 1, natoms ! To prevent calculating r(1,2) when it will calculate (2,1) first
                                dx = coords(i, 1) - coords(j, 1)
                                dy = coords(i, 2) - coords(j, 2)
                                dz = coords(i, 3) - coords(j, 3)
                                distances(i, j) = sqrt(dx**2 + dy**2 + dz**2) ! CONVERSION to nm 
                                distances(j, i) = distances(i, j) ! distance between i and j is equal to between j and i
                        end do
                end do
        end subroutine compute_interatomic_distances
end module md_math

