program test_suite
    use parameters
    use lattice
    use read_field_config
    implicit none

    ! Read parameters from input file
    call initialise_parameters("parameter_file.txt")

    ! Check parameters are correct
    call check_parameters()

    ! Test parameters
    call check_success(test_dependant_parameters)

    ! Setup lattice pointers
    call setup_lattice()

    ! Test move function
    call check_success(test_move)

    ! Test loading of gauge field
    call check_success(test_read_gauge_field)

    ! Test unitarisation procedure
    call check_success(test_unitarise_SVD)

    ! Test smearing and blocking procedures
    call check_success(test_blocking_smearing)

    ! Print success statement
    write(*, "(a)") "----------ALL TESTS PASSED----------"

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
                write(6, '(f0.9, a, f0.9, a)', advance='no') real(matrix(i,j)), "+",&
                aimag(matrix(i,j)), "i    "
            enddo
            print *, ""
        enddo
        print *, ""
    end subroutine

    ! Quit program if a test has failed
    subroutine check_success(test_function)
        implicit none

        interface
            subroutine test_function(ierr)
                use parameters
                use lattice
                use read_field_config
                implicit none
                logical, intent(out) :: ierr
            end subroutine
        end interface

        logical :: error_code

        call test_function(error_code)
        if (.not.error_code) then
            stop "Test failed. Program halting."
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
        write(*, "(a)") "Parameter check PASSED."
        print *, ""
    end subroutine check_parameters

    ! Test dependant parameters
    subroutine test_dependant_parameters(ierr)
        implicit none
        logical, intent(out) :: ierr

        ierr = .true.

        ! Check dependant parameters
        if (SLICE_VOLUME /= LX1 * LX2 * LX3) then
            ierr = .false.
            write(*, "(a)") "SLICE_VOLUME test FAILED"
            return
        end if
        if (LATTICE_VOLUME /= SLICE_VOLUME * LX4) then
            ierr = .false.
            write(*, "(a)") "LATTICE_VOLUME test FAILED"
            return
        end if
        if (NCOL2 /= NCOL * NCOL) then
            ierr = .false.
            write(*, "(a)") "NCOL2 test FAILED"
            return
        end if
        if (MAX_DELTA_T /= LX4 / 2) then
            ierr = .false.
            write(*, "(a)") "MAX_DELTA_T test FAILED"
            return
        end if

        ! Print success
        write(*, "(a)") "Dependant parameter check PASSED"
        print *, ""
    end subroutine test_dependant_parameters

    ! Test move function
    subroutine test_move(ierr)
        implicit none
        logical, intent(out) :: ierr

        integer :: IUP_old(LATTICE_VOLUME, 4), IDN_old(LATTICE_VOLUME, 4), &
        IUPB_old(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1), IDNB_old(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1), &
        IUP_new(LATTICE_VOLUME, 4), IDN_new(LATTICE_VOLUME, 4), &
        IUPB_new(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1), IDNB_new(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1), &
        site, mu, blocking_level
        logical :: IUP_equal, IDN_equal, IUPB_equal, IDNB_equal

        ! Load old IUP, IDN from NEIGHBS.DAT file
        open(11, file='NEIGHBS.DAT', form='unformatted', status='old', access='stream')
        read(11) IUP_old
        read(11) IDN_old
        close(11)

        ! Load old IUPB, IDNB from BLOCKED_NEIGHBS.DAT file
        open(11, file='BLOCKED_NEIGHBS.DAT', form='unformatted', status='old', access='stream')
        read(11) IUPB_old
        read(11) IDNB_old
        close(11)

        ! Construct IUP and IDN using move function
        do mu = 1, 4
            do site = 1, LATTICE_VOLUME
                IUP_new(site, mu) = move(site, mu)
                IDN_new(site, mu) = move(site, -mu)
            enddo
        enddo

        ! Construct IUPB and IDNB using move function
        do blocking_level = 1, MAX_BLOCKING_LEVEL+1
            do mu = 1, 3
                do site = 1, SLICE_VOLUME
                    IUPB_new(site, mu, blocking_level) = move(site, mu, blocking_level)
                    IDNB_new(site, mu, blocking_level) = move(site, -mu, blocking_level)
                enddo
            enddo
        enddo

        ! Check if all moving arrays match
        IUP_equal = all(IUP_old == IUP_new)
        IDN_equal = all(IDN_old == IDN_new)
        IUPB_equal = all(IUPB_old == IUPB_new)
        IDNB_equal = all(IDNB_old == IDNB_new)

        ! Check if test has succeeded
        ierr = IUP_equal.and.IDN_equal.and.IUPB_equal.and.IDNB_equal
        if (ierr) then
            write(*, '(a)') "move test PASSED"
        else
            write(*, '(a)') "move test FAILED"
        endif
        print *, ""
    end subroutine

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
        ierr = all(abs(gauge_field_check - gauge_field_correct) <= epsilon(1.0))
        if (ierr) then
            write(*, "(a)") "read_gauge_field test PASSED"
        else
            write(*, "(a)") "read_gauge_field test FAILED"
            print *, gauge_field_check
        end if
        print *, ""
    end subroutine test_read_gauge_field

    ! Test unitarisation subroutine
    subroutine test_unitarise_SVD(ierr)
        implicit none
        logical, intent(out) :: ierr

        integer, parameter :: num_tests = 5
        character(len=16) :: file_config_id
        character(len=256) :: directory
        complex(real64) :: gauge_field(NCOL, NCOL, LATTICE_VOLUME, 4), &
        matrix_sum(NCOL, NCOL), trial_U(NCOL, NCOL), trial_identity(NCOL, NCOL)
        integer :: site, random_sites(num_tests), mu, random_directions(num_tests), i, j, k
        real :: random_numbers(2*num_tests)
        logical :: diagonal_good, off_diagonal_good, diagonal_mask(NCOL, NCOL), off_diagonal_mask(NCOL, NCOL)

        ! Set directory path
        write(file_config_id, "(i0)") CONFIG_START
        directory=trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)

        ! Load gauge field into memory
        call read_gauge_field(directory, gauge_field)

        ! Get 5 random sites and 5 random directions
        call random_number(random_numbers)
        random_numbers(1:num_tests) = random_numbers(1:num_tests) * LATTICE_VOLUME
        random_sites = ceiling(random_numbers(1:5))
        random_numbers(num_tests+1:2*num_tests) = random_numbers(num_tests+1:2*num_tests) * 4
        random_directions = ceiling(random_numbers(num_tests+1:2*num_tests))

        ! Get mask for diagonal and off-diagonal elements
        diagonal_mask = .false.
        do i = 1, NCOL
            diagonal_mask(i,i) = .true.
        enddo
        off_diagonal_mask = .not.diagonal_mask

        ! For each random site and direction, add 2 consectutive links from that site and in that direction. Apply the unitarisation procedure to the sum and test the result for unitarity. If any are not unitary, return ierr = .false. and halt execution.
        ierr = .true.
        do k = 1, num_tests
            site = random_sites(k)
            mu = random_directions(k)

            ! Sum consective links in direction mu from site
            matrix_sum = gauge_field(:,:,site,mu) + gauge_field(:,:,move(site,mu),mu)

            ! Apply unitaristation procedure
            trial_U = unitarise_SVD(matrix_sum)

            ! Form trial identity matrix
            trial_identity = matmul(herm(trial_U), trial_U)

            ! Test if trial identity matrix is sufficiently close to the identity
            diagonal_good = all((abs(trial_identity - cmplx(1.0, 0.0)) < TOL_SVD).or.off_diagonal_mask)
            off_diagonal_good = all((abs(trial_identity) < TOL_SVD).or.diagonal_mask)
            ierr = diagonal_good.and.off_diagonal_good
            if (.not.ierr) then
                write(*, "(a)") "unitarise_SVD test FAILED"
                write(*, "(a)") "Trial unitary matrix:"
                call print_matrix(trial_U)
                write(*, "(a)") "Trial identity matrix:"
                call print_matrix(trial_identity)
                return
            endif
        enddo
        write(*, "(a)") "unitarise_SVD test PASSED"
	print *, ""
    end subroutine

    ! Test smearing and blocking of configurations
    subroutine test_blocking_smearing(ierr)
        implicit none
        logical, intent(out) :: ierr

        complex(real32) :: smear_check(NCOL, NCOL, SLICE_VOLUME, 3), &
        blok_check(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        complex(real64) :: gauge_field(NCOL, NCOL, LATTICE_VOLUME, 4), &
        gauge_field_smeared(NCOL, NCOL, SLICE_VOLUME, 3), &
        gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        character(len=16) :: file_config_id
        character(len=256) :: directory
        logical :: blok_equal, smear_equal

        !! Load data from old code
        ! Load smear_check from SMEARED_SAVE.DAT file
        open(11, file='SMEARED_SAVE.DAT', form='unformatted', status='old', access='stream')
        read(11) smear_check
        close(11)

        ! Load blok_check from BLOCKED_SAVE.DAT file
        open(11, file='BLOCKED_SAVE.DAT', form='unformatted', status='old', access='stream')
        read(11) blok_check
        close(11)


        !! Load field configuration from file given in parameter file
        ! Set directory path
        write(file_config_id, "(i0)") CONFIG_START
        directory=trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)

        ! Load gauge field into memory
        call read_gauge_field(directory, gauge_field)


        !! Smear and block configuration using new functions
        ! Block configuration
        gauge_field_blocked = get_blocked_gauge_field(gauge_field(:, :, 1:SLICE_VOLUME, 1:3))

        ! Smear configuration at 4th blocking level
        gauge_field_smeared = get_smeared_gauge_field(gauge_field_blocked(:,:,:,:,4), 4)


        !! Check if smeared and blocked configurations are equal to those produced during old calculation
        smear_equal = all(abs(gauge_field_smeared - smear_check) <= max(1e-4, epsilon(1.0), TOL_SVD))
        blok_equal = all(abs(gauge_field_blocked - blok_check) <= max(1e-4, epsilon(1.0), TOL_SVD))
        ierr = smear_equal.and.blok_equal

        ! Output results
        if (ierr) then
            write(*, "(a)") "smearing and blocking tests PASSED"
        else
            write(*, "(a)") "smearing and blocking tests FAILED"
            if (smear_equal) then
                write(*, "(a)") "smearing test PASSED"
            else
                write(*, "(a)") "smearing test FAILED"
            endif
            if (blok_equal) then
                write(*, "(a)") "blocking test PASSED"
            else
                write(*, "(a)") "blocking test FAILED"
            endif
        end if
        print *, ""
    end subroutine
end program test_suite
