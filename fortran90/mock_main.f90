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

    complex(real64), allocatable :: gauge_field(:,:,:,:), blok_timeslice(:,:,:,:,:), A11(:,:)
    integer(int32) :: mu, site_idx, blevel, site_in, site_out


    ! Read parameters from input file
    call initialise_parameters("parameter_file.txt")

    ! Setup lattice pointers
    call setup_lattice()

    ! Allocate arrays
    allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 4))
    allocate(A11(NCOL, NCOL))

    A11 = get_I(2)
    call print_matrix(A11)

    ! Load gauge field
    call read_gauge_field("/Users/alexi/Work/phd/torelons/cnfg/run1_52x26x26x26nc2rADJnf2b2.300000m1.000000n134780", gauge_field)

    ! Get the blocked first time slice
    blok_timeslice = get_blocked_gauge_field(gauge_field(:, :, 1:SLICE_VOLUME, 1:3))

    ! get square pulses in `squares_up` and `squares_down`
    call get_square_pulses(blok_timeslice, 1)
    call get_plaquettes(blok_timeslice, 1)

    call loop_1(1, site_out, A11, 2)
    call print_matrix(A11)
    
    ! do blevel=1,MAX_BLOCKING_LEVEL
    !     do mu=1,2
    !         do site_idx=1,SLICE_VOLUME
    !             call print_matrix(squares_down(:, :, site_idx, mu, blevel))
    !         enddo
    !     enddo
    ! enddo

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
