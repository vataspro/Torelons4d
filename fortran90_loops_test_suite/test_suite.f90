program operator_test_suite
    use iso_fortran_env, only : real32, real64
    use parameters
    use read_field_config
    use lattice
    use thermal_lines
    use operator_dispatch
    implicit none

    character(len=512) :: gauge_filename
    character(len=512) :: output_filename
    integer :: argument_count
    integer :: x, y, z, mu, site
    integer :: site_x, site_y, site_z, time_slice
    integer :: first_loop, last_loop, tested_loops
    integer :: blocking_level, loop_index, loop_number
    integer :: in_site, out_site
    complex(real32), allocatable :: gauge_field_raw(:,:,:,:,:,:,:)
    complex(real64), allocatable :: gauge_field_slice(:,:,:,:)
    complex(real64), allocatable :: gauge_field_blocked(:,:,:,:,:)
    complex(real64), allocatable :: A11(:,:)

    call initialise_parameters('parameter_file.txt')
    call setup_lattice()

    argument_count = command_argument_count()
    if (argument_count >= 1) then
        call get_command_argument(1, gauge_filename)
    else
        gauge_filename = '../fortran77_scalar_looptest/run1_52x26x26x26nc2rADJnf2b2.300000m1.000000n1'
    endif

    if (argument_count >= 2) then
        call get_command_argument(2, output_filename)
    else
        output_filename = 'operator_outputs.dat'
    endif

    time_slice = 1
    site_x = 1
    site_y = 1
    site_z = 1
    first_loop = minval(LOOP_NUMBERS)
    last_loop = maxval(LOOP_NUMBERS)

    if (argument_count >= 3) call read_integer_argument(3, time_slice)
    if (argument_count >= 4) call read_integer_argument(4, site_x)
    if (argument_count >= 5) call read_integer_argument(5, site_y)
    if (argument_count >= 6) call read_integer_argument(6, site_z)
    if (argument_count >= 7) call read_integer_argument(7, first_loop)
    if (argument_count >= 8) call read_integer_argument(8, last_loop)

    call validate_inputs(time_slice, site_x, site_y, site_z, first_loop, last_loop)

    allocate(gauge_field_raw(NCOL, NCOL, LX1, LX2, LX3, LX4, 4))
    allocate(gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3))
    allocate(gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL))
    allocate(A11(NCOL, NCOL))

    call read_gauge_field(trim(gauge_filename), gauge_field_raw)

    do mu = 1, 3
        do z = 1, LX3
            do y = 1, LX2
                do x = 1, LX1
                    site = site_index(x, y, z)
                    gauge_field_slice(:,:,site,mu) = cmplx( &
                        real(gauge_field_raw(:,:,x,y,z,time_slice,mu), kind=real64), &
                        real(aimag(gauge_field_raw(:,:,x,y,z,time_slice,mu)), kind=real64), &
                        kind=real64 &
                    )
                enddo
            enddo
        enddo
    enddo

    gauge_field_blocked = get_blocked_gauge_field(gauge_field_slice)

    call get_square_pulses(gauge_field_blocked, 1)
    call get_plaquettes(gauge_field_blocked, 1)

    open(unit=77, file=trim(output_filename), status='replace', action='write')
    write(77,'(A)') 'Generated thermal_lines operator loop-by-loop A11 matrices'
    write(77,900) time_slice, site_x, site_y, site_z, first_loop, last_loop
900 format('Selection: N4=',I0,'  X=',I0,'  Y=',I0,'  Z=',I0,'  FIRST_LOOP=',I0,'  LAST_LOOP=',I0)

    in_site = site_index(site_x, site_y, site_z)
    tested_loops = 0

    do blocking_level = 1, MAX_BLOCKING_LEVEL
        do loop_index = 1, NUMBER_OF_LOOPS
            loop_number = LOOP_NUMBERS(loop_index)
            if (loop_number < first_loop .or. loop_number > last_loop) cycle

            call identity_matrix(A11)
            call call_loop(loop_number, in_site, out_site, A11, blocking_level, gauge_field_blocked)
            call dump_A11(77, time_slice, blocking_level, loop_number, in_site, out_site, A11)
            tested_loops = tested_loops + 1
        enddo
    enddo

    close(77)

    if (tested_loops == 0) error stop 'No generated loop subroutines matched the requested loop range'

contains

    subroutine read_integer_argument(argument_index, value)
        integer, intent(in) :: argument_index
        integer, intent(out) :: value
        character(len=512) :: argument_text
        integer :: io_status

        call get_command_argument(argument_index, argument_text)
        read(argument_text, *, iostat=io_status) value
        if (io_status /= 0) error stop 'Invalid integer command-line argument'
    end subroutine read_integer_argument

    subroutine validate_inputs(time_slice, site_x, site_y, site_z, first_loop, last_loop)
        integer, intent(in) :: time_slice, site_x, site_y, site_z, first_loop, last_loop

        if (time_slice < 1 .or. time_slice > LX4) error stop 'Requested time slice is outside 1:LX4'
        if (site_x < 1 .or. site_x > LX1) error stop 'Requested x coordinate is outside 1:LX1'
        if (site_y < 1 .or. site_y > LX2) error stop 'Requested y coordinate is outside 1:LX2'
        if (site_z < 1 .or. site_z > LX3) error stop 'Requested z coordinate is outside 1:LX3'
        if (first_loop > last_loop) error stop 'Requested first loop is larger than last loop'
    end subroutine validate_inputs

    subroutine identity_matrix(matrix)
        complex(real64), intent(out) :: matrix(:,:)
        integer :: colour_index

        matrix = cmplx(0.0_real64, 0.0_real64, kind=real64)
        do colour_index = 1, size(matrix, 1)
            matrix(colour_index, colour_index) = cmplx(1.0_real64, 0.0_real64, kind=real64)
        enddo
    end subroutine identity_matrix

    subroutine dump_A11(unit_number, time_slice, blocking_level, loop_number, in_site, out_site, A11)
        integer, intent(in) :: unit_number, time_slice, blocking_level, loop_number, in_site, out_site
        complex(real64), intent(in) :: A11(:,:)

        write(unit_number,'(A)') '----------------------------------------'
        write(unit_number,'(A,I0)') 'Loop ', loop_number
        write(unit_number,901) time_slice, blocking_level, in_site, out_site
        write(unit_number,'(A)') 'A11 ='
        write(unit_number,902) real(A11(1,1)), aimag(A11(1,1)), real(A11(1,2)), aimag(A11(1,2))
        write(unit_number,902) real(A11(2,1)), aimag(A11(2,1)), real(A11(2,2)), aimag(A11(2,2))
        write(unit_number,*)
901     format('N4=',I0,'  IBL=',I0,'  MN=',I0,'  OUT_SITE=',I0)
902     format('(',ES16.8,',',ES16.8,')  (',ES16.8,',',ES16.8,')')
    end subroutine dump_A11

end program operator_test_suite
