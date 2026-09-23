program test_suite
    use parameters
    use lattice
    use read_field_config
    implicit none

    logical :: error_code

    ! Read parameters from input file
    call initialise_parameters("parameter_file.txt")
    error_code = .true.

    ! Test parameters
    call test_dependant_parameters(error_code)
    call check_success(error_code)

    ! Test loading of gauge field
    call test_read_gauge_field(error_code)
    call check_success(error_code)

    ! Print success statement
    if (error_code) write(*, "(a)") "All tests passed."

    contains

    ! Quit program if a test has failed
    subroutine check_success(ierr)
        implicit none

        logical, intent(in) :: ierr

        if (.not.ierr) then
            error stop "Test failed. Program halting."
        end if
    end subroutine check_success

    ! Check parameters are correct for testing
    subroutine check_parameters()
        implicit none

        ! Check loaded parameters
        if (NCOL /= 2) then
            error stop "NCOL must be equal to 2"
        end if
        if (LX1 /= 26) then
            error stop "LX1 must be equal to 26"
        end if
        if (LX2 /= 26) then
            error stop "LX2 must be equal to 26"
        end if
        if (LX3 /= 26) then
            error stop "LX3 must be equal to 26"
        end if
        if (LX4 /= 52) then
            error stop "LX4 must be equal to 52"
        end if

        ! Output success
        write(*, "(a)") "Parameter check passed."
    end subroutine check_parameters

    ! Test dependant parameters
    subroutine test_dependant_parameters(ierr)
        implicit none
        logical, intent(out) :: ierr

        ! Check dependant parameters
        if (SLICE_VOLUME /= LX1 * LX2 * LX3) then
            ierr = .false.
            write(*, "(a)") "SLICE_VOLUME test failed"
            return
        end if
        if (LATTICE_VOLUME /= SLICE_VOLUME * LX4) then
            ierr = .false.
            write(*, "(a)") "LATTICE_VOLUME test failed"
            return
        end if
        if (NCOL2 /= NCOL * NCOL) then
            ierr = .false.
            write(*, "(a)") "NCOL2 test failed"
            return
        end if
        if (MAX_DELTA_T /= LX4 / 2) then
            ierr = .false.
            write(*, "(a)") "MAX_DELTA_T test failed"
            return
        end if

        ! Print success
        write(*, "(a)") "test_dependant_parameters passed"
    end subroutine test_dependant_parameters

    ! Test loading of field configurations
    subroutine test_read_gauge_field(ierr)
        implicit none
        logical, intent(out) :: ierr

        character(len=16) :: file_config_id
        character(len=256) :: directory
        complex(real64) :: gauge_field_conf(NCOL, NCOL, LATTICE_VOLUME, 4)
        complex(real64) :: gauge_field_correct(5), gauge_field_check(5)

        ! Set directory path
        write(file_config_id, "(i0)") CONFIG_START
        directory=trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)

        ! Load gauge field into memory
        call read_gauge_field(directory, gauge_field_conf)

        ! Set correct links
        gauge_field_correct = &
        [(0.545266926,0.530812383), &
        (-0.632124960,-0.146082729), &
        (0.632124960,-0.146082729), &
        (0.545266926,-0.530812383), &
        (-0.778977633,8.882141858E-02)]

        ! Get links that should match the set links from the loaded gauge field
        gauge_field_check(1) = gauge_field_conf(1,1,1,1)
        gauge_field_check(2) = gauge_field_conf(2,1,1,1)
        gauge_field_check(3) = gauge_field_conf(1,2,1,1)
        gauge_field_check(4) = gauge_field_conf(2,2,1,1)
        gauge_field_check(5) = gauge_field_conf(1,1,2,1)

        ! Check first 5 elements of loaded gauge field against correct values
        ierr = all(abs(gauge_field_check - gauge_field_correct) < epsilon(1.0))
        if (.not.ierr) then
            write(*, "(a)") "read_gauge_field test failed"
            print *, gauge_field_check
            return
        else
            write(*, "(a)") "read_gauge_field test passed"
        end if
    end subroutine test_read_gauge_field
end program test_suite
