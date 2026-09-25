program mock_main
    use thermal_lines
    use read_field_config
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! ---------------------------------------------- MOCK MAIN ----------------------------------------------
    !
    ! Mock main file to check compilation and operators on a test configuration
    !
    ! Happy testing!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(real64), allocatable :: gauge_field(:,:,:,:)
    integer(int32) :: mu, site_idx

    ! Read parameters from input file
    call initialise_parameters("parameter_file.txt")

    ! Setup lattice pointers
    call setup_lattice()

    ! Allocate gauge field so that it only needs to be loaded once
    allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 4))
    call read_gauge_field("/Users/alexi/Work/phd/torelons/cnfg/run1_52x26x26x26nc2rADJnf2b2.300000m1.000000n134780", gauge_field)

    
    do mu=1,4
        do site_idx=1,LATTICE_VOLUME
            call print_matrix(gauge_field(:, :, site_idx, mu))
        enddo
    enddo

    contains

    ! Print matrix
    subroutine print_matrix(matrix)
        implicit none
        complex(real64), dimension(:,:), intent(in) :: matrix

        integer :: i, j
        integer, dimension(2) :: N

        N = shape(matrix)
        do i = 1, N(1)
            do j = 1, N(2)
                write(6, '(f0.6, a, f0.6, a)', advance='no') real(matrix(i,j)), "+",&
                aimag(matrix(i,j)), "i    "
            enddo
            print *, ""
        enddo
        print *, ""
    end subroutine

end program mock_main
