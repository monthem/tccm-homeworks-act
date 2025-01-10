program test_md_io
        use md_io 
        implicit none

        integer :: i, natoms
        double precision, allocatable :: coords(:,:), masses(:)

        ! Call the read_input_file subroutine 
        call read_input_file("inp.txt", natoms, coords, masses)

        ! Print them back to terminal for verification
        print *, "Number of atoms: ", natoms
        do i = 1, natoms
                print '(F8.3, F8.3, F8.3, F8.3)', coords(i,1), coords(i,2), coords(i,3), masses(i)
        end do
end program test_md_io
