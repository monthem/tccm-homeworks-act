module md_io
        implicit none
        contains

        subroutine read_input_file(filename, natoms, coords, masses)
                ! This subroutine reads an input file and extracts the number of atoms, their coordinates and masses
                
                ! Arguments
                character(len=*), intent(in) :: filename ! Input file name
                integer, intent(out) :: natoms
                double precision, allocatable, intent(out) :: coords(:,:), masses(:) ! Arrays for coordinates and masses
                
                ! Local variables
                integer :: i, io_status
                logical :: file_exists
                double precision :: mass
                double precision, parameter :: amu_2_kgpermol = 0.00166054d0
                double precision :: coordinates(3)

                ! Check file existence
                inquire(file=filename, exist=file_exists)
                if (.not. file_exists) then
                        print *, "Error: Input file not found: ", filename 
                        stop
                end if

                ! Open the file for reading
                open(unit=10, file=filename, status='old', action='read', iostat = io_status)
                if (io_status /= 0) then
                        print *, "Error: Unable to open the file: ", filename
                        stop
                end if

                ! Read the number of atoms
                read(10, *) natoms ! Reads the first line into natoms variable

                ! Memory allocation according to the number of atoms
                allocate(coords(natoms, 3), masses(natoms)) ! a nat x 3 array and a 1D array with natoms elements
                ! Read the coordinates and masses
                do i = 1, natoms
                        read(10, *, iostat=io_status) coordinates, mass
                        
                        if (io_status /=0) then
                                print *, "Error reading atom data at line ", i + 1
                                stop
                        end if

                        coords(i, :) = coordinates ! Stores the i-th atom coordinates in the i-th row of coords array
                        !masses(i) = mass
                        masses(i) = mass * amu_2_kgpermol ! Stores the i-th atom mass at the i-th place in the masses array
                end do

                ! Close the file
                close(10)
        end subroutine read_input_file

        subroutine write_xyz_trajectory(filename, timestep, atom, natoms, coords, E_t, E_k, E_p)

                character(len=*), intent(in) :: filename
                character(2), intent(in) :: atom
                integer, intent(in) :: natoms 
                double precision, intent(in) :: timestep, E_t, E_k, E_p
                double precision, intent(in) :: coords(natoms, 3)
                integer :: i
                
                open(unit=20, file=filename, status='unknown', action='write', position='append')
                write(20, *) natoms ! first line 
                write(20, '(A, A, F10.2, A, F10.4, A, F10.4, A, F10.4)') "Lattice=""80.0 0.0 0.0 0.0 80.0 0.0 0.0 0.0 80.0""", " Current time: ", timestep, ", E = ", E_t, ", T = ", E_k, ", V = ", E_p
                do i = 1, natoms
                        write(20, '(A2, 3X, F10.3, 3X, F10.3, 3X, F10.3)') trim(adjustl(atom)), coords(i, 1) * 10.0d0, coords(i, 2) * 10.0d0, coords(i, 3) * 10.0d0
                end do

                close(20)
        end subroutine write_xyz_trajectory

end module md_io
 
