program test_md_math
        use md_io
        use md_math
        implicit none

        integer :: i, natoms
        double precision, allocatable :: coords(:,:), masses(:)
        double precision, allocatable :: distances(:,:)

        call read_input_file("inp.txt", natoms, coords, masses)

        allocate(distances(natoms, natoms))

        call compute_interatomic_distances(natoms, coords, distances)

        print *, "Distance matrix (from inp.txt):"
        do i = 1, natoms
                print "(3(F8.3,1X))", distances(i, :)
        end do

        ! Memory clean up
        deallocate(coords, masses, distances)
end program test_md_math
