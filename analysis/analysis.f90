module analysis
    use iso_fortran_env, only : int8, int32, int64, real32, real64
    use ieee_arithmetic, only: ieee_is_finite
    use linalg

    implicit none

    integer :: NUM_OPERATORS, MAX_DELTA_T, NUM_BINS, NCONFIG, NUM_STATES, LT

    contains

    ! Provide help with options
    subroutine options_help()
        write(*, '(a)') "[Help]    Usage: ./analysis <options>"
        write(*, '(a)') "[Help]"
        write(*, '(a)') "[Help]    Required options for all states:"
        write(*, '(a, 24x, a)') "[Help]     * '-J <spin>'", &
        "Spin of flux tube state for masses to be calculated"

        write(*, '(a, 7x, a)') "[Help]     * '-q <longitudinal momentum>'", &
        "Longitudinal momentum of flux tube state for masses to be calculated"

        write(*, '(a)') "[Help]"

        write(*, '(a)') "[Help]    Required options for q=0, J=0,2 states:"

        write(*, '(a, 11x, a)') "[Help]     * '-P <transverse parity>'", &
        "Transverse parity of flux tube state for masses to be calculated"

        write(*, '(a, 12x, a)') "[Help]     * '-R <intrinsic parity>'", &
        "Intrinsic parity of flux tube state for masses to be calculated"

        write(*, '(a)') "[Help]"

        write(*, '(a)') "[Help]    Required options for q=0, J=1 states:"

        write(*, '(a, 12x, a)') "[Help]     * '-R <intrinsic parity>'", &
        "Intrinsic parity of flux tube state for masses to be calculated"

        write(*, '(a)') "[Help]"

        write(*, '(a)') "[Help]    Required options for q/=0, J=0,2 states:"

        write(*, '(a, 11x, a)') "[Help]     * '-P <transverse parity>'", &
        "Transverse parity of flux tube state for masses to be calculated"

        write(*, '(a)') "[Help]"

        write(*, '(a)') "[Help]    No further options for q/=0, J=1 states"

        write(*, '(a)') "[Help]"

        write(*, '(a)') "[Help]    Other passable options:"

        write(*, '(a, 9x, a)') "[Help]     * '-n <number_of_states>'", &
        "Total of number of states calculated in GEVP. Defaults to n = 3 if not passed."

        write(*, '(a, 11x, a)') "[Help]     * '-runs <number of runs>'", &
        " Total number of runs to used in calculation of states. Defaults to runs = 10 if not passed."

        write(*, '(a, 4x, a)') "[Help]     * '--max-dt' <maximum delta t>", &
        "Report effective masses upto and including maximum delta t. Defaults to max_dt = 7 if not passed."

        write(*, '(a, 23x, a, a)') "[Help]     * '--lines'", &
        "Show correlation functions and effective masses for individual line operators. ", &
        "Defaults to FALSE if not passed."

        write(*, '(a, 14x, a, a)') "[Help]     * '--lines-filtered'", &
        "Show correlation functions and effective masses for individual line operators, ", &
        "after filtering bad and identical operators. Defaults to FALSE if not passed."

        write(*, '(a, 25x, a)') "[Help]     * '--fit'", &
        "Use numerical cosh fit in calculation of effective masses. Defaults to TRUE if not passed."

        write(*, '(a, 23x, a)') "[Help]     * '--exact'", &
        "Use exact cosh formula in calculation of effective masses. Defaults to FALSE if not passed."

        stop
    end subroutine

    ! Calculate errors on quantity from jacknife bins
    pure function jacknife_error(data, average) result(error)
        implicit none
        real(real64), intent(in) :: data(:)
        real(real64), intent(in) :: average
        real(real64) :: error

        integer :: bin, nbins

        nbins = size(data)

        error = 0.0d0
        do bin = 1, nbins
            error = error + (data(bin) - average)**2
        enddo
        error = sqrt(real(nbins-1, kind=real64) * error / real(nbins, kind=real64))
    end function

    ! Read vevs and correlation matrix into memory, and remove duplicate operators and zero operators
    subroutine read_clean_data(vevs_filename, corr_matrix_filename, vevs, corr_matrix, nruns, read_nconfig, lines_passed)
        implicit none
        character(len=*), intent(in) :: vevs_filename, corr_matrix_filename
        complex(real64), allocatable, intent(out) :: vevs(:,:), corr_matrix(:,:,:,:)
        integer, intent(in) :: nruns
        logical, intent(in), optional :: read_nconfig, lines_passed

        character(len=32) :: filepath, file_run_id
        complex(real64), allocatable :: vevs_from_file(:,:), corr_matrix_from_file(:,:,:,:), &
        diagonal_correlators(:,:), corr_matrix_from_file_avg(:,:,:), &
        corr_matrix_from_file_jackknife(:,:,:,:), vevs_from_file_avg(:)
        integer(int32) :: dim1, dim2, dim3, dim4, i, j, nconfig_check, dt, run
        integer, allocatable :: distinct_operators(:)
        logical, allocatable :: operator_comparison(:)
        logical :: read_nconfig_default

        ! Set read_nconfig
        if (present(read_nconfig)) then
            read_nconfig_default = read_nconfig
        else
            read_nconfig_default = .true.
        endif

        !! Read vevs from file from run 1
        ! Write run1 to filepath
        filepath = 'run1/'

        ! Open file
        open(11, file=trim(filepath) // trim(vevs_filename), form='unformatted', access='stream', status='old')

        ! Read number of configurations
        if (read_nconfig_default) then
            read(11) NCONFIG
        endif

        ! Read dimensions of array
        read(11) dim1 ! NUM_OPERATORS
        read(11) dim2 ! NUM_BINS (per run)

        ! Allocate vevs array to take data from nruns runs
        allocate(vevs_from_file(dim1, nruns * dim2))

        ! Read vevs data from file
        read(11) vevs_from_file(:, 1:dim2)

        ! Close file
        close(11)

        !! Read correlation matrix from file
        ! Open file
        open(11, file=trim(filepath) // trim(corr_matrix_filename), &
        form='unformatted', access='stream', status='old')

        ! Check number of configurations is the same as for the vevs
        if (read_nconfig_default) then
            read(11) nconfig_check
            if (nconfig_check /= NCONFIG) then
                write(*, '(a)') &
                "[Error]    Number of configurations not consistent between vevs and correlation matrix"
                stop
            endif
        endif

        ! Read dimensions of array
        read(11) dim1 ! NUM_OPERATORS
        read(11) dim2 ! NUM_OPERATORS
        read(11) dim3 ! MAX_DELTA_T + 1
        read(11) dim4 ! NUM_BINS (per run)

        ! Allocate correlation matrix array to take data from nruns runs
        allocate(corr_matrix_from_file(dim1, dim2, dim3, nruns * dim4))

        ! Read vevs data from file
        read(11) corr_matrix_from_file(:, :, :, 1:dim4)

        ! Close file
        close(11)


        !! Find values of NUM_OPERATORS, MAX_DELTA_T and NUM_BINS
        ! Assume operators contained in vevs(:,j), corr_matrix(:,:,k,l)
        NUM_OPERATORS = size(vevs_from_file, dim=1)

        ! Check there are the same number of operators in corr_matrix as in vevs
        if (size(corr_matrix_from_file, dim=1) /= NUM_OPERATORS) &
        error stop "vevs and corr_matrix incompatible: different number of operators"
        if (size(corr_matrix_from_file, dim=1) /= size(corr_matrix_from_file, dim=2)) &
        error stop "corr_matrix not square"

        ! Assume time seperations contained in corr_matrix(i,j,:,l)
        MAX_DELTA_T = size(corr_matrix_from_file, dim=3) - 1
        LT = 2 * MAX_DELTA_T

        ! Assume bins contained in vevs(i,:), corr_matrix(i,j,k,:)
        NUM_BINS = size(vevs_from_file, dim=2)

        ! Check there are the same number of bins in corr_matrix as in vevs
        if (size(corr_matrix_from_file, dim=4) /= NUM_BINS) &
        error stop "vevs and corr_matrix incompatible: different number of bins"

        ! Output dimensions found
        write(*, '(a, i0)') "[Info]    Number of operators from file: ", NUM_OPERATORS
        write(*, '(a, i0)') "[Info]    Number time seperations from file: ", MAX_DELTA_T + 1
        write(*, '(a, i0)') "[Info]    Number of bins from file: ", NUM_BINS


        !! Read data from other requested runs
        do run = 2, nruns
            ! Set file path
            write(file_run_id, '(i0)') run
            filepath = 'run' // trim(file_run_id) // '/'

            ! Read vev dimensions from file
            open(11, &
                file=trim(filepath) // trim(vevs_filename), &
                form='unformatted', &
                access='stream', &
                status='old'&
                )
            if (read_nconfig_default) then
                read(11) nconfig_check
                if (nconfig_check /= NCONFIG) then
                    write(*, '(a, i0)') &
                    "[Error]    Number of configurations not consistent between run 1 and run ", run
                    stop
                endif
            endif
            read(11) dim1
            read(11) dim2

            ! Check dimensions are correct
            if (dim1 /= NUM_OPERATORS) then
                write(*, '(a, i0)') &
                "[Error]    Number of operators in vevs file not consistent between run 1 and run ", run
                stop
            endif
            if (nruns * dim2 /= NUM_BINS) then
                write(*, '(a, i0)') &
                "[Error]    Number of bins in vevs file not consistent between run 1 and run ", run
                stop
            endif

            ! Read vevs from file
            read(11) vevs_from_file(:, 1 + (run-1)*dim2 : run*dim2)
            close(11)

            ! Read correlation matrix dimensions from file
            open(11, &
                file=trim(filepath) // trim(corr_matrix_filename), &
                form='unformatted', &
                access='stream', &
                status='old'&
                )
            if (read_nconfig_default) then
                read(11) nconfig_check
                if (nconfig_check /= NCONFIG) then
                    write(*, '(a, i0)') &
                    "[Error]    Number of configurations not consistent between run 1 and run ", run
                    stop
                endif
            endif
            read(11) dim1
            read(11) dim2
            read(11) dim3
            read(11) dim4

            ! Check dimensions are correct
            if (dim1 /= NUM_OPERATORS .or. dim2 /= NUM_OPERATORS) then
                write(*, '(a, i0)') &
                "[Error]    Number of operators in corr_matrix file not consistent between run 1 and run ", run
                stop
            endif
            if (dim3 /= MAX_DELTA_T+1) then
                write(*, '(a, i0)') &
                "[Error]    Number of time slices in corr_matrix file not consistent between run 1 and run ", &
                run
                stop
            endif
            if (nruns * dim4 /= NUM_BINS) then
                write(*, '(a, i0)') &
                "[Error]    Number of bins in corr_matrix file not consistent between run 1 and run ", run
                stop
            endif

            ! Read correlation matrix from file
            read(11) corr_matrix_from_file(:, :, :, 1 + (run-1)*dim4 : run*dim4)
            close(11)
        enddo


        ! Output diagonal correlators if requested
        print *, ""
        if (lines_passed) then
            allocate(corr_matrix_from_file_avg(NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1), &
            corr_matrix_from_file_jackknife(NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1, NUM_BINS), &
            vevs_from_file_avg(NUM_OPERATORS))
            corr_matrix_from_file_avg = sum(corr_matrix_from_file, dim=4)
            corr_matrix_from_file_jackknife = corr_matrix_from_file
            vevs_from_file_avg = sum(vevs_from_file, dim=2)
            call jackknife(vevs_from_file, corr_matrix_from_file_jackknife, vevs_from_file_avg, corr_matrix_from_file_avg)
            write(*, '(a)') " *"
            write(*, '(a)') " *******************************************************"
            write(*, '(a)') " *"
            do i = 1, NUM_OPERATORS
                write(*, '(a, i0)') "Operator ", i
                call output_effective_masses(corr_matrix_from_file_jackknife(i, i, :, :)%re, &
                corr_matrix_from_file_avg(i, i, :)%re, .true.)
                print *, ""
            enddo
        endif
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        print *, ""


        !! Find duplicate operators by comparing diagonal elements of correlation matrix at zero seperation
        ! Compute average diagonal correlation matrix elements at 0, 1 seperation
        allocate(diagonal_correlators(NUM_OPERATORS, 2))
        do i = 1, NUM_OPERATORS
            do dt = 1, 2
                diagonal_correlators(i, dt) = sum(corr_matrix_from_file(i, i, dt, :))
            enddo
            diagonal_correlators(i, 2) = diagonal_correlators(i, 2) / diagonal_correlators(i, 1)
        enddo

        ! Find non-duplicate operators by comparing diagonal correlators
        distinct_operators = [1]
        do i = 2, NUM_OPERATORS
            ! Exclude operators that are exactly 0
            if ((abs(diagonal_correlators(i, 1)) <= epsilon(1.0d0))) then
                write(*, '(a, i0, a, f12.10)') "[Info]    Operator ", i, " diagonal correlator at dt=0 is ", &
                diagonal_correlators(i, 1)%re
                cycle
            ! Exclude operators that have a correlator that vanishes at 1 seperation
            elseif ((diagonal_correlators(i, 2)%re <= epsilon(1.0d0))) then
                write(*, '(a, i0, a, f12.10)') "[Info]    Operator ", i, " diagonal correlator at dt=1 is ", &
                diagonal_correlators(i, 2)%re
                cycle
            ! Exclude identical operators
            else
                allocate(operator_comparison(i-1))
                do j = 1, i-1
                    operator_comparison(j) &
                    = (1.0 - abs(sum(corr_matrix_from_file(i,j,1,:))) &
                    / sqrt(diagonal_correlators(i, 1)%re * diagonal_correlators(j, 1)%re) &
                    <= epsilon(1.0d0))
                enddo
                if (.not.any(operator_comparison)) distinct_operators = [distinct_operators, i]
                deallocate(operator_comparison)
            endif
        enddo

        !! Fill correlation matrix with distinct operators only
        ! Reset number of operators to number of distinct operators
        NUM_OPERATORS = size(distinct_operators)
        write(*, '(a, i0)') "[Info]    Number of operators after filter: ", NUM_OPERATORS

        ! Allocate vevs and correlation matrix outputs
        allocate(vevs(NUM_OPERATORS, NUM_BINS), &
        corr_matrix(NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1, NUM_BINS))

        ! Fill vevs with distinct operators
        do i = 1, NUM_OPERATORS
            vevs(i, :) = vevs_from_file(distinct_operators(i), :)
        enddo

        ! Fill corr_matrix with distinct operators
        do j = 1, NUM_OPERATORS
            do i = 1, NUM_OPERATORS
                corr_matrix(i, j, :, :) &
                = corr_matrix_from_file(distinct_operators(i), distinct_operators(j), :, :)
            enddo
        enddo
    end subroutine

    ! Arrange bins into jacknife bins
    subroutine jackknife(vevs, corr_matrix, vevs_avg, corr_matrix_avg)
        implicit none
        complex(real64), intent(inout) :: vevs(:,:), corr_matrix(:,:,:,:)
        complex(real64), intent(out), allocatable :: vevs_avg(:), corr_matrix_avg(:,:,:)
        integer :: i, j, dt

        ! Allocate average containers
        allocate(vevs_avg(NUM_OPERATORS), corr_matrix_avg(NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1))

        ! Divide vevs and correlation matrix by number of measurements
        vevs = vevs / (NCONFIG * MAX_DELTA_T * 2)
        corr_matrix = corr_matrix / (NCONFIG * MAX_DELTA_T * 2)

        ! Sum vevs and correlation matrix over bins to get average values
        vevs_avg = sum(vevs, dim=2)
        corr_matrix_avg = sum(corr_matrix, dim=4)

        ! Get jacknife bins
        if (NUM_BINS > 1) then
            ! Vev bins
            do concurrent (i = 1:NUM_OPERATORS)
                vevs(i, :) = (vevs_avg(i) - vevs(i,:)) &
                * real(NUM_BINS, kind=real64) / real(NUM_BINS-1, kind=real64)
            enddo

            ! Correlation matrix bins
            do concurrent (i = 1:NUM_OPERATORS, j = 1:NUM_OPERATORS, dt = 1:MAX_DELTA_T+1)
                corr_matrix(i, j, dt, :) = (corr_matrix_avg(i, j, dt) - corr_matrix(i, j, dt, :)) &
                * real(NUM_BINS, kind=real64) / real(NUM_BINS-1, kind=real64)
            enddo
        endif
    end subroutine

    ! Remove vevs from correlation matrix
    subroutine subtract_vevs(vevs, corr_matrix, vevs_avg, corr_matrix_avg)
        implicit none
        complex(real64), intent(in) :: vevs(:,:), vevs_avg(:)
        complex(real64), intent(inout) :: corr_matrix(:,:,:,:), corr_matrix_avg(:,:,:)
        integer :: i, j, bin

        ! Subtract vevs from average correlation matrices to form connected correlation matrcies
        do concurrent (i = 1:NUM_OPERATORS, j = 1:NUM_OPERATORS)
            corr_matrix_avg(i, j, :) = corr_matrix_avg(i, j, :) - vevs_avg(i) * conjg(vevs_avg(j))
        enddo

        ! Subtract vevs from jacknife bins of correaltion matrices
        do concurrent (i = 1:NUM_OPERATORS, j = 1:NUM_OPERATORS, bin = 1:NUM_BINS)
            corr_matrix(i, j, :, bin) = corr_matrix(i, j, :, bin) - vevs(i, bin) * conjg(vevs(j, bin))
        enddo
    end subroutine

    ! Run GEVP on correlation matrix. Return mass upper bounds and eigenvectors used to weight operators.
    subroutine GEVP(corr_matrix, corr_matrix_avg, mass_bound, mass_bound_avg, eigenvectors, eigenvectors_avg)
        implicit none
        complex(real64), intent(in) :: corr_matrix(:,:,:,:), corr_matrix_avg(:,:,:)
        real(real64), intent(out), allocatable :: mass_bound(:,:), mass_bound_avg(:)
        complex(real64), intent(out), allocatable :: eigenvectors(:,:,:), eigenvectors_avg(:,:)

        real(real64), allocatable :: eigenvalues(:)
        integer :: bin

        ! Check dimensions
        if (size(corr_matrix, dim=1) /= NUM_OPERATORS .or. size(corr_matrix, dim=2) /= NUM_OPERATORS) &
        error stop "[Error][GEVP]    corr_matrix_avg must have NUM_OPERATORS as size of first two dimensions"

        ! Allocate mass_bound and eigen arrays
        allocate(mass_bound(NUM_STATES, NUM_BINS), mass_bound_avg(NUM_STATES), &
        eigenvalues(NUM_STATES), &
        eigenvectors(NUM_OPERATORS, NUM_STATES, NUM_BINS), eigenvectors_avg(NUM_OPERATORS, NUM_STATES))

        ! Run GEVP on correlation matrix average values
        call generalised_eigenpairs(corr_matrix_avg(:,:,2), corr_matrix_avg(:,:,1), &
        eigenvalues, eigenvectors_avg)

        ! Calculate upper bounds on masses of states from eigenvalues
        mass_bound_avg = -log(eigenvalues)

        ! Run GEVP and get mass upper bounds on all jackknife bins
        do bin = 1, NUM_BINS
            call generalised_eigenpairs(corr_matrix(:,:,2,bin), corr_matrix(:,:,1,bin), &
            eigenvalues, eigenvectors(:,:,bin))
            mass_bound(:,bin) = -log(eigenvalues)
        enddo
    end subroutine

    ! Construct normalised correlation functions of operators from eigenvectors
    subroutine construct_correlators(&
        corr_matrix, &
        corr_matrix_avg, &
        eigenvectors, &
        eigenvectors_avg, &
        cfn, &
        cfn_avg, &
        cfn_error&
        )
        implicit none
        complex(real64), intent(in) :: corr_matrix(:,:,:,:), corr_matrix_avg(:,:,:), &
        eigenvectors(:,:,:), eigenvectors_avg(:,:)
        real(real64), allocatable, intent(out) :: cfn(:,:,:), cfn_avg(:,:), cfn_error(:,:)

        complex(real64) :: scalar_product
        integer :: state, dt, bin

        ! Check dimensions
        if (any(shape(corr_matrix) /= [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1, NUM_BINS])) then
            write(*, '(a, a)') "[Error][construct_correlators]    corr_matrix must have shape", &
            " [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1, NUM_BINS]"
            error stop
        endif

        if (any(shape(corr_matrix_avg) /= [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1])) then
            write(*, '(a, a)') "[Error][construct_correlators]    corr_matrix_avg must have shape", &
            " [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1]"
            error stop
        endif

        if (any(shape(eigenvectors) /= [NUM_OPERATORS, NUM_STATES, NUM_BINS])) &
        error stop &
        "[Error][construct_correlators]    eigenvectors must have shape [NUM_OPERATORS, NUM_STATES, NUM_BINS]"

        if (any(shape(eigenvectors_avg) /= [NUM_OPERATORS, NUM_STATES])) &
        error stop &
        "[Error][construct_correlators]    eigenvectors_avg must have shape [NUM_OPERATORS, NUM_STATES]"

        ! Allocate correlation function and average correlation funcion arrays
        allocate(cfn(MAX_DELTA_T+1, NUM_BINS, NUM_STATES), cfn_avg(MAX_DELTA_T+1, NUM_STATES), &
        cfn_error(MAX_DELTA_T+1, NUM_STATES))

        ! Construct correlation functions of states by computing the scalar product v^dagger*C(t)*v
        do state = 1, NUM_STATES
            do dt = 1, MAX_DELTA_T+1
                ! Compute scalar product
                scalar_product = dot_product(eigenvectors_avg(:, state), &
                matmul(corr_matrix_avg(:, :, dt), eigenvectors_avg(:, state)))

                ! Scalar product is real, so the real part is the correlation function entry
                cfn_avg(dt, state) = scalar_product%re
            enddo
        enddo

        do state = 1, NUM_STATES
            do bin = 1, NUM_BINS
                do dt = 1, MAX_DELTA_T+1
                    ! Compute scalar product
                    scalar_product = dot_product(eigenvectors(:, state, bin), &
                    matmul(corr_matrix(:, :, dt, bin), eigenvectors(:, state, bin)))

                    ! Scalar product is real, so the real part is the correlation function entry
                    cfn(dt, bin, state) = scalar_product%re
                enddo
            enddo
        enddo

        ! Normalise correlation functions so that cfn(0) = 1 for all states
        do state = 1, NUM_STATES
            cfn_avg(:, state) = cfn_avg(:, state) / cfn_avg(1, state)
            do bin = 1, NUM_BINS
                cfn(:, bin, state) = cfn(:, bin, state) / cfn(1, bin, state)
            enddo
        enddo

        ! Calculate errors on correlation functions
        do concurrent (dt = 1:MAX_DELTA_T+1, state = 1:NUM_STATES)
            cfn_error(dt, state) = jacknife_error(cfn(dt, :, state), cfn_avg(dt, state))
        enddo
    end subroutine

    ! Calculate effective mass exactly from 3 neighbouring correlation function data points
    pure function get_effective_mass_exact(cfn) result(m_eff)
        implicit none
        real(real64), intent(in) :: cfn(-1:1)
        real(real64) :: m_eff, x

        ! Calculate argument of arcosh function
        x = (cfn(1) + cfn(-1))/(2*cfn(0))

        ! Determine if arcosh(x) exists. If it does not, assign m_eff = 0.0
        if (x < 1) then
            m_eff = 0.0
        else
            m_eff = acosh(x)
        endif
    end function

    ! Calculate effective mass numerically using ratio of 2 neighbouring correlation function data points
    pure function get_effective_mass_fit(cfn, t, rel_tol, abs_tol, max_iter) result(m_eff)
        implicit none
        integer, intent(in) :: t
        real(real64), intent(in) :: cfn(t-1:t)
        real(real64), intent(in), optional :: rel_tol, abs_tol
        integer, intent(in), optional :: max_iter
        real(real64) :: m_eff

        real(real64) :: dt, ratio, rtol, atol, x, x_trial, width, low, high, F, dF
        integer :: max_iterations, iter

        ! Check t is valid
        if (t < 1 .or. t > MAX_DELTA_T) then
            error stop "[Error]    t invalid - must be in range [1, N_t/2]"
        endif

        ! Calculate dt = L_t/2 - t
        dt = real(MAX_DELTA_T) - real(t)

        ! Compute correlation function ratio
        ratio = cfn(t-1) / cfn(t)

        ! Check correlation function is falling
        if (ratio < 1) then
            m_eff = 0.0
            return
        endif

        ! Set default parameters
        if (present(rel_tol)) then
            rtol = rel_tol
        else
            rtol = 1e-13
        endif
        if (present(abs_tol)) then
            atol = abs_tol
        else
            atol = 1e-14
        endif
        if (present(max_iter)) then
            max_iterations = max_iter
        else
            max_iterations = 100
        endif

        ! Set initial bracket boundaries
        low = log(ratio)
        high = acosh(ratio)

        ! Set initial estimate as midpoint of bracket
        x = (low + high) / 2

        ! Conduct safeguarded Newton-Raphson iteration
        do iter = 1, max_iterations
            ! Compute F(x)
            F = x + log(1.0 + exp(-2.0 * (dt + 1) * x)) - log(1.0 + exp(-2.0 * dt * x)) - log(ratio)

            ! Update bracket
            if (F < 0) then
                low = x
            elseif (F > 0) then
                high = x
            else
                m_eff = x
                return
            endif

            ! Compute bracket width
            width = high - low

            ! Check for convergence
            if (width <= 2.0 * (atol + rtol * (low + high) / 2.0)) then
                m_eff = (low + high) / 2
                return
            endif

            ! Compute F'(x)
            dF = (dt + 1) * tanh((dt + 1) * x) - dt * tanh(dt * x)

            ! Produce trial update to x using Newton's method
            x_trial = x - F / dF

            ! Accept trial update if it is within the bracket. If not, use bisection to update x.
            if (x_trial > low .and. x_trial < high) then
                x = x_trial
            else
                x = (low + high) / 2
            endif
        enddo

        ! Report error if process did not converge
        error stop "[Error]    Numerical solver for effective mass did not converge"
    end function

    ! Calculate effective masses and errors from normalised correlation functions
    subroutine get_effective_masses(cfn, cfn_avg, effective_mass_avg, effective_mass_error, fit_mode)
        implicit none
        real(real64), intent(in) :: cfn(:,:), cfn_avg(:)
        real(real64), intent(out) :: effective_mass_avg(:), effective_mass_error(:)
        logical, intent(in), optional :: fit_mode

        real(real64), allocatable :: effective_mass_bins(:,:)
        integer :: dt, bin

        ! Check dimensions
        if (any(shape(cfn) /= [MAX_DELTA_T+1, NUM_BINS])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    cfn must have shape", &
            " [MAX_DELTA_T+1, NUM_BINS]"
            error stop
        endif

        if (any(shape(cfn_avg) /= [MAX_DELTA_T+1])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    cfn_avg must have shape", &
            " [MAX_DELTA_T+1]"
            error stop
        endif

        if (any(shape(effective_mass_avg) /= [MAX_DELTA_T])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    effective_mass_avg must have shape", &
            " [MAX_DELTA_T]"
            error stop
        endif

        if (any(shape(effective_mass_error) /= [MAX_DELTA_T])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    effective_mass_error must have shape", &
            " [MAX_DELTA_T]"
            error stop
        endif

        ! Allocate effective mass array
        allocate(effective_mass_bins(MAX_DELTA_T, NUM_BINS))

        ! Calculate effective masses for all time steps and states
        ! Use effective mass mode requested - default to cosh fit
        if (present(fit_mode) .and. .not.fit_mode) then
            do concurrent (dt = 1:MAX_DELTA_T-1)
                effective_mass_avg(dt) = get_effective_mass_exact(cfn_avg(dt:dt+2))
            enddo

            do concurrent (dt = 1:MAX_DELTA_T-1, bin = 1:NUM_BINS)
                effective_mass_bins(dt, bin) = get_effective_mass_exact(cfn(dt:dt+2, bin))
            enddo

            effective_mass_avg(MAX_DELTA_T) = get_effective_mass_exact(&
            [cfn_avg(MAX_DELTA_T-1), cfn_avg(MAX_DELTA_T), cfn_avg(MAX_DELTA_T-1)])

            do concurrent (bin = 1:NUM_BINS)
                effective_mass_bins(MAX_DELTA_T, bin) = get_effective_mass_exact(&
                [cfn(MAX_DELTA_T-1, bin), cfn(MAX_DELTA_T, bin), cfn(MAX_DELTA_T-1, bin)])
            enddo
        else
            do concurrent (dt = 1:MAX_DELTA_T)
                effective_mass_avg(dt) = get_effective_mass_fit(cfn_avg(dt:dt+1), dt)
            enddo

            do concurrent (dt = 1:MAX_DELTA_T, bin = 1:NUM_BINS)
                effective_mass_bins(dt, bin) = get_effective_mass_fit(cfn(dt:dt+1, bin), dt)
            enddo
        endif

        ! Calculate errors on effective masses
        do concurrent (dt = 1:MAX_DELTA_T)
            effective_mass_error(dt) = jacknife_error(effective_mass_bins(dt, :), effective_mass_avg(dt))
        enddo
    end subroutine

    ! Get effective masses for all states
    subroutine get_effective_masses_all_states(cfn, cfn_avg, effective_mass_avg, effective_mass_error, fit_mode)
        implicit none
        real(real64), intent(in) :: cfn(:,:,:), cfn_avg(:,:)
        real(real64), allocatable, intent(out) :: effective_mass_avg(:,:), effective_mass_error(:,:)
        logical, intent(in), optional :: fit_mode

        integer :: i

        ! Check dimensions
        if (any(shape(cfn) /= [MAX_DELTA_T+1, NUM_BINS, NUM_STATES])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    cfn must have shape", &
            " [MAX_DELTA_T+1, NUM_BINS, NUM_STATES]"
            error stop
        endif

        if (any(shape(cfn_avg) /= [MAX_DELTA_T+1, NUM_STATES])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    cfn_avg must have shape", &
            " [MAX_DELTA_T+1, NUM_STATES]"
            error stop
        endif

        ! Allocate effecive mass arrays
        allocate(effective_mass_avg(MAX_DELTA_T, NUM_STATES), effective_mass_error(MAX_DELTA_T, NUM_STATES))

        ! Calculate effective masses
        do i = 1, NUM_STATES
            if (present(fit_mode)) then
                call get_effective_masses(cfn(:,:,i), cfn_avg(:,i), &
                effective_mass_avg(:,i), effective_mass_error(:,i), fit_mode)
            else
                call get_effective_masses(cfn(:,:,i), cfn_avg(:,i), &
                effective_mass_avg(:,i), effective_mass_error(:,i))
            endif
        enddo
    end subroutine

    ! Output effective masses from a given correlation function
    subroutine output_effective_masses(cfn, cfn_avg, fit_mode)
        implicit none
        real(real64), intent(in) :: cfn(:,:), cfn_avg(:)
        logical, intent(in), optional :: fit_mode

        real(real64), allocatable :: effective_mass_avg(:), effective_mass_errors(:), &
        cfn_norm(:,:), cfn_avg_norm(:), cfn_error(:)
        integer :: i, b, delta_t

        ! Allocate temporary arrays
        allocate(cfn_norm(size(cfn, dim=1), size(cfn, dim=2)), cfn_avg_norm(size(cfn_avg)), &
        cfn_error(size(cfn_avg)), effective_mass_avg(size(cfn_avg)-1), effective_mass_errors(size(cfn_avg)-1))

        ! Calculate effective masses and errors
        call get_effective_masses(cfn, cfn_avg, effective_mass_avg, effective_mass_errors, fit_mode)

        ! Normalise cfn bins
        cfn_avg_norm = cfn_avg / cfn_avg(1)
        do b = 1, size(cfn, dim=2)
            cfn_norm(:, b) = cfn(:, b) / cfn(1, b)
        enddo

        ! Calculate errors on cfn function
        do delta_t = 1, size(cfn_avg)
            cfn_error(delta_t) = jacknife_error(cfn_norm(delta_t, :), cfn_avg_norm(delta_t))
        enddo

        ! Output effective masses and errors
        write(6, '(a, i0, a, f10.8, a, f10.8)') "delta t = ", 0, "  cfn = ",&
        cfn_avg_norm(1), "  error = ", cfn_error(1)
        do i = 1, 7
            write(6, '(a, i0, 4(a, f10.8))') "delta t = ", i, "  cfn = ",&
            cfn_avg_norm(i+1), "  error = ", cfn_error(i+1), "  m_eff = ",&
            effective_mass_avg(i), "  error = ", effective_mass_errors(i)
        enddo
        do i = 1, 7
            write(6, '(a, i0, a, f10.8, a, f10.8, a)') &
            "(", i, ",", effective_mass_avg(i), ") +- (0,", effective_mass_errors(i), ")"
        enddo
    end subroutine

    ! Output correlation functions and effective masses with errors
    subroutine output_results(&
        cfn_avg, &
        cfn_error, &
        effective_mass_avg, &
        effective_mass_error, &
        mass_bound, &
        mass_bound_avg, &
        max_dt&
        )
        real(real64), intent(in) :: cfn_avg(:,:), cfn_error(:,:), &
        effective_mass_avg(:,:), effective_mass_error(:,:)
        real(real64), intent(in) :: mass_bound(:,:), mass_bound_avg(:)
        integer, intent(in), optional :: max_dt

        real(real64) :: mass_bound_error
        integer :: i, state, dt, maximum_seperation

        ! Set default max_dt to 7
        if (present(max_dt)) then
            maximum_seperation = max_dt
        else
            maximum_seperation = 7
        endif

        ! Check dimensions
        if (any(shape(mass_bound) /= [NUM_STATES, NUM_BINS])) then
            write(*, '(a)') "[Error][output_results]    mass_bound must have shape [NUM_STATES, NUM_BINS]"
            error stop
        endif

        if (size(mass_bound_avg) /= NUM_STATES) then
            write(*, '(a)') "[Error][output_results]    mass_bound_avg must have size NUM_STATES"
            error stop
        endif

        if (any(shape(cfn_avg) /= [MAX_DELTA_T+1, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    cfn_avg must have shape", &
            " [MAX_DELTA_T+1, k]"
            error stop
        endif

        if (any(shape(cfn_error) /= [MAX_DELTA_T+1, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    cfn_error must have shape", &
            " [MAX_DELTA_T+1, k]"
            error stop
        endif

        if (any(shape(effective_mass_avg) /= [MAX_DELTA_T, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    effective_mass_avg must have shape", &
            " [MAX_DELTA_T, k]"
            error stop
        endif

        if (any(shape(effective_mass_error) /= [MAX_DELTA_T, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    effective_mass_error must have shape", &
            " [MAX_DELTA_T, k]"
            error stop
        endif

        ! Output correlation functions and effective masses fo each state in turn. Iterate through arrays in reverse order, as the lightest state is the last in each array
        state = 0
        do i = NUM_STATES, 1, -1
            ! Output state identifier
            print *, ""
            write(*, '(a)') "--------------------------------------------------"
            if (state == 0) then
                write(*, '(a)') "Ground state"
            else
                write(*, '(a, i0)') "Excited state ", state
            endif

            ! Output mass upper bound and error
            mass_bound_error = jacknife_error(mass_bound(i, :), mass_bound_avg(i))
            write(*, '(2(a, f6.4))') "Mass upper bound:    ", mass_bound_avg(i), " +- ", mass_bound_error
            write(*, '(a)') "--------------------------------------------------"

            ! Output table header row
            write(*, '(a)') "dt    cfn         error       m_eff       error   "
            write(*, '(a)') "--------------------------------------------------"

            ! Output correlation function at 0 seperation with errors
            write(*, '(i2, 2(a, f8.6), a)') 0, "    ", cfn_avg(1, i), "    ", cfn_error(1, i), &
            "    ########    ########"

            ! Output correlation functions and effective masses with errors upto supplied maximum seperation
            do dt = 1, maximum_seperation
                write(*, '(i2, 4(a, f8.6))') dt, "    ", cfn_avg(dt+1, i), "    ", cfn_error(dt+1, i), &
                "    ", effective_mass_avg(dt, i), "    ", effective_mass_error(dt, i)
            enddo

            ! Increment state counter
            state = state + 1
        enddo
    end subroutine


    ! ---------------------------------- Correlated cosh fitting functions ---------------------------------- !

    ! Cosh function correlation functions should be fitted to
    pure function cosh_fit(t, m)
        implicit none

        integer, intent(in) :: t
        real(real64), intent(in) :: m
        real(real64) :: cosh_fit

        cosh_fit = exp(-m * real(t, kind=real64)) + exp(-m * (real(LT - t, kind=real64)))
    end function

    ! Derivate of cosh function wrt m
    pure function dcosh_fit(t, m)
        implicit none

        integer, intent(in) :: t
        real(real64), intent(in) :: m
        real(real64) :: dcosh_fit

        dcosh_fit = -(real(t, kind=real64) * exp(-m * real(t, kind=real64)) &
        + real(LT - t, kind=real64) * exp(-m * (real(LT - t, kind=real64))))
    end function

    ! Second derivate of cosh function wrt m
    pure function d2cosh_fit(t, m)
        implicit none

        integer, intent(in) :: t
        real(real64), intent(in) :: m
        real(real64) :: d2cosh_fit

        d2cosh_fit = (real(t, kind=real64)**2) * exp(-m * real(t, kind=real64)) &
        + (real(LT - t, kind=real64)**2) * exp(-m * (real(LT - t, kind=real64)))
    end function

    ! Construct covariance matrix from normalised correlation functions and Jackknife bins
    subroutine get_covariance_matrix(jackknife_bins, data, covariance_matrix)
        implicit none
        real(real64), intent(in) :: jackknife_bins(:,:), data(:)
        real(real64), intent(out) :: covariance_matrix(:,:)

        integer :: t1, t2, ndata

        ! Get number of data points
        ndata = size(data)

        ! Check sizes of arrays
        if (any(shape(jackknife_bins) /= [ndata, NUM_BINS])) then
            write(*, '(a, a)') "[Error][get_covariance_matrix]    ", &
            "jackknife_bins must have shape [ndata, NUM_BINS]"
            stop
        endif

        if (any(shape(covariance_matrix) /= [ndata, ndata])) then
            write(*, '(a, a)') "[Error][get_covariance_matrix]    ", &
            "covariance_matrix must have shape [ndata, ndata]"
        endif

        ! Construct covariance matrix from jackknife samples
        covariance_matrix = 0.0d0
        do concurrent (t1 = 1:ndata, t2 = 1:ndata)
            covariance_matrix(t1, t2) &
            = sum((jackknife_bins(t1, :) - data(t1)) * (jackknife_bins(t2, :) - data(t2)))
        enddo
        covariance_matrix = covariance_matrix * real(NUM_BINS - 1, kind=real64) / real(NUM_BINS, kind=real64)
    end subroutine

    ! Diagonalise a window of the covariance matrix
    subroutine diagonalise_covariance_matrix(covariance_matrix, t_min, t_max, eigenvalues, eigenvectors, base)
        implicit none
        real(real64), intent(in) :: covariance_matrix(:,:)
        integer, intent(in) :: t_min, t_max
        integer, intent(in), optional :: base
        real(real64), intent(out) :: eigenvectors(:,:), eigenvalues(:)

        integer :: window_size, ndata, t_min_, t_max_

        ! Get size of window from [t_min, t_max]
        window_size = t_max - t_min + 1
        if (window_size <= 0) then
            write(*, '(a, a)') "[Error][diagonalise_covariance_matrix]    ", &
            "window cannot have size <= 0"
        endif

        ! Modify t_min and t_max if covariance_matrix is not one-based
        if (present(base)) then
            t_min_ = t_min - base + 1
            t_max_ = t_max - base + 1
        else
            t_min_ = t_min
            t_max_ = t_max
        endif

        ! Check sizes of other arrays
        ndata = size(covariance_matrix, dim=1)
        if (any(shape(covariance_matrix) /= [ndata, ndata])) then
            write(*, '(a, a)') "[Error][diagonalise_covariance_matrix]    ", &
            "covariance_matrix must have shape [ndata, ndata]"
        endif

        if (any(shape(eigenvectors) /= [window_size, window_size])) then
            write(*, '(a, a)') "[Error][diagonalise_covariance_matrix]    ", &
            "eigenvectors must have shape [window size, window_size]"
        endif

        if (size(eigenvalues) /= window_size) then
            write(*, '(a, a)') "[Error][diagonalise_covariance_matrix]    ", &
            "eigenvalues must have shape [window size]"
        endif

        ! Diagnoalise covariance matrix
        call eigenpairs(covariance_matrix(t_min_:t_max_, t_min_:t_max_), eigenvalues, eigenvectors)

        ! Check eigenvalues are well-conditions
        if (any(.not. ieee_is_finite(eigenvalues))) then
            write(*, '(a, 4x, a)') "[Error][diagnonalise_covariance_matrix]", &
            "Non-finite covariance eigenvalue"
            error stop
        endif
        if (any(eigenvalues <= 0.0_real64)) then
            write(*, '(a, 4x, a)') "[Error][diagnoalise_covariance_matrix]", &
            "Covariance not positive definite"
            error stop
        endif
        if (minval(eigenvalues) < epsilon(1.0_real64) * maxval(eigenvalues)) then
            write(*, '(a, 4x, a)') "[Warning][diagonalise_covariance_matrix]", &
            "Condition number of covariance matrix of order machine precision"
        endif
    end subroutine

    ! Vector representing (derivative of) cosh function in given t-window in eigenbasis given
    function transformed_cosh_vector(t_min, t_max, m, eigenvectors, cosh_derivative)
        implicit none

        interface
            pure function cosh_derivative(time, mass)
                use iso_fortran_env, only: real64
                implicit none
                integer, intent(in) :: time
                real(real64), intent(in) :: mass
                real(real64) :: cosh_derivative
            end function
        end interface

        integer, intent(in) :: t_min, t_max
        real(real64), intent(in) :: m, eigenvectors(t_min:t_max,t_min:t_max)
        real(real64) :: transformed_cosh_vector(t_min:t_max)

        real(real64) :: cosh_derivative_vector(t_min:t_max)
        integer :: t

        ! Fill vector with (derivative of) cosh function evaluated at relevant t-values
        do concurrent (t = t_min:t_max)
            cosh_derivative_vector(t) = cosh_derivative(t, m)
        enddo

        ! Transform vector into eigenbasis: L^{1/2} u(m) = U^T g(m)
        transformed_cosh_vector = matmul(transpose(eigenvectors), cosh_derivative_vector)
    end function

    ! Vector representing correlation function in eigenbasis given
    function transformed_cfn(t_min, t_max, eigenvectors, cfn)
        implicit none

        integer, intent(in) :: t_min, t_max
        real(real64), intent(in) :: cfn(0:MAX_DELTA_T), eigenvectors(t_min:t_max,t_min:t_max)
        real(real64) :: transformed_cfn(t_min:t_max)

        ! Transform t-window of correlation function into eigenbasis: L^{1/2} y = U^T C
        transformed_cfn = matmul(transpose(eigenvectors), cfn(t_min:t_max))
    end function

    ! Compute the function F(m) = |u(m)| (y \cdot \hat{u}(m)) and its derivative. F(m) is the function whose roots correspond to minima of the correlated chi-squared (at the optimum cosh amplitude) for a given t-window.
    subroutine get_F(t_min, t_max, m, eigenvalues, eigenvectors, transformed_data, F, dF, F_scale, dF_scale)
        implicit none

        integer, intent(in) :: t_min, t_max
        real(real64), intent(in) :: m, eigenvalues(t_min:t_max), eigenvectors(t_min:t_max,t_min:t_max), &
        transformed_data(t_min:t_max)
        real(real64), intent(out) :: F
        real(real64), intent(out), optional :: dF, F_scale, dF_scale

        real(real64) :: transformed_model(t_min:t_max), dtransformed_model(t_min:t_max), &
        d2transformed_model(t_min:t_max), U, dU, d2U, Y, dY, d2Y, F_fraction_term, dF_fraction_term
        integer :: t

        ! Get model vector and derivatives in eigenbasis
        transformed_model = transformed_cosh_vector(t_min, t_max, m, eigenvectors, cosh_fit)
        dtransformed_model = transformed_cosh_vector(t_min, t_max, m, eigenvectors, dcosh_fit)

        ! Second derivatives only need to be calculated if F'(m) is required
        if (present(dF) .or. present(dF_scale)) then
            d2transformed_model = transformed_cosh_vector(t_min, t_max, m, eigenvectors, d2cosh_fit)
        endif

        ! Get variables U, U', U'', Y, Y', Y'' in which we express F(m) and F'(m)
        U = 0.0
        dU = 0.0
        d2U = 0.0
        Y = 0.0
        dY = 0.0
        d2Y = 0.0
        do t = t_min, t_max
            ! U = u(m) * u(m)
            U = U + transformed_model(t) ** 2 / eigenvalues(t)

            ! Y = y * u(m)
            Y = Y + transformed_data(t) * transformed_model(t) / eigenvalues(t)

            ! U' = 2 (u'(m) * u(m))
            dU = dU + 2 * dtransformed_model(t) * transformed_model(t) / eigenvalues(t)

            ! Y' = y * u'(m)
            dY = dY + transformed_data(t) * dtransformed_model(t) / eigenvalues(t)

            ! Second derivatives only need to be calculated if F'(m) is required
            if (present(dF) .or. present(dF_scale)) then
                ! U'' = 2 (u''(m) * u(m) + u'(m) * u'(m))
                d2U = d2U + 2 * (d2transformed_model(t) * transformed_model(t) + dtransformed_model(t) ** 2) &
                / eigenvalues(t)

                ! Y'' = y * u''(m)
                d2Y = d2Y + transformed_data(t) * d2transformed_model(t) / eigenvalues(t)
            endif
        enddo

        ! Calculate F(m)
        F_fraction_term = (dU * Y) / (2 * U)
        F = dY - F_fraction_term

        ! Calculate F'(m) if required
        if (present(dF)) then
            dF_fraction_term = (U * d2U * Y + U * dU * dY - dU**2 * Y) / (2 * U**2)
            dF = d2Y - dF_fraction_term
        endif

        ! Calculate scale of F(m)
        if (present(F_scale)) then
            F_scale = abs(dY) + abs(F_fraction_term)
        endif

        ! Calculate scale of F'(m)
        if (present(dF_scale)) then
            if (.not. present(dF)) then
                dF_fraction_term = (U * d2U * Y + U * dU * dY - dU**2 * Y) / (2 * U**2)
            endif
            dF_scale = abs(d2Y) + abs(dF_fraction_term)
        endif
    end subroutine

    ! Compute the correlated chi-squared of a given model state mass at the optimal cosh amplitude
    function get_chi_squared(t_min, t_max, m, eigenvalues, eigenvectors, transformed_data) &
        result(chi_squared)
        implicit none

        integer, intent(in) :: t_min, t_max
        real(real64), intent(in) :: m, eigenvalues(t_min:t_max), eigenvectors(t_min:t_max,t_min:t_max), &
        transformed_data(t_min:t_max)
        real(real64) :: chi_squared

        real(real64) :: transformed_model(t_min:t_max), U, Y, A
        integer :: t

        ! Get model vector in eigenbasis
        transformed_model = transformed_cosh_vector(t_min, t_max, m, eigenvectors, cosh_fit)

        ! Get variables U, Y, A in which we express chi-squared
        U = 0.0
        Y = 0.0
        do t = t_min, t_max
            ! U = u(m) * u(m)
            U = U + transformed_model(t) ** 2 / eigenvalues(t)

            ! Y = y * u(m)
            Y = Y + transformed_data(t) * transformed_model(t) / eigenvalues(t)
        enddo

        ! Calculate chi-squared in a numerically stable way
        A = Y / U
        chi_squared = sum( &
            (transformed_data - A * transformed_model)**2 / eigenvalues &
        )
    end function

    ! Find the mass value that minimises the correlated chi-squared over a given mass range by using a Newton-bisection algorithm to find a zero of F(m) within the given range
    subroutine minimise_chi_squared(t_min, t_max, eigenvalues, eigenvectors, transformed_data, mass_bracket, &
        mass, chi_squared, rel_tol, abs_tol, max_iter, info)
        implicit none

        integer, intent(in) :: t_min, t_max
        real(real64), intent(in) :: eigenvalues(t_min:t_max), eigenvectors(t_min:t_max, t_min:t_max), &
        transformed_data(t_min:t_max), mass_bracket(2)
        real(real64), intent(out) :: mass, chi_squared
        real(real64), intent(in), optional :: rel_tol, abs_tol
        integer, intent(in), optional :: max_iter
        integer, intent(out), optional :: info

        real(real64) :: rtol, atol, low, high, x, x_newton, x_bisection, width, F, dF, F_low, F_high
        integer :: max_iterations, iter
        logical :: converged

        ! Set default parameters
        if (present(rel_tol)) then
            rtol = rel_tol
        else
            rtol = 1e-13
        endif
        if (present(abs_tol)) then
            atol = abs_tol
        else
            atol = 1e-14
        endif
        if (present(max_iter)) then
            max_iterations = max_iter
        else
            max_iterations = 100
        endif

        ! Set initial bracket boundaries
        low = mass_bracket(1)
        high = mass_bracket(2)

        ! Set initial estimate as midpoint of bracket
        x = (low + high) / 2

        ! Determine if F is rising or falling in the given interval
        call get_F(t_min, t_max, low, eigenvalues, eigenvectors, transformed_data, F_low)
        call get_F(t_min, t_max, high, eigenvalues, eigenvectors, transformed_data, F_high)
        if (F_low * F_high > 0.0_real64) then
            write(*, '(a, 4x, a)') "[Error][minimise_chi_squared]", &
            "No root exists in mass bracket given"
            error stop
        endif

        ! Conduct safeguarded Newton-Raphson iteration
        converged = .false.
        do iter = 1, max_iterations
            ! Compute F(m) and F'(m)
            call get_F(t_min, t_max, x, eigenvalues, eigenvectors, transformed_data, F, dF)

            ! Update bracket
            if (sgn(F) == sgn(F_low)) then
                low = x
            elseif (sgn(F) == sgn(F_high)) then
                high = x
            else
                mass = x
                converged = .true.
                exit
            endif

            ! Compute bracket width
            width = high - low

            ! Check for convergence
            if (width <= 2.0 * (atol + rtol * (low + high) / 2.0)) then
                mass = (low + high) / 2
                converged = .true.
                exit
            endif

            ! Update x using either Newton's method or bisection
            ! Safeguard x by checking size of derivative and update to x via Newton's method. If Newton's update will give a NaN value, use bisection.
            x_bisection = (low + high) / 2.0_real64
            if (abs(dF) > epsilon(1.0_real64)) then
                ! Produce trial update to x using Newton's method
                x_newton = x - F / dF

                ! Accept trial update if it is within the bracket. If not, use bisection to update x.
                if (ieee_is_finite(x_newton)) then
                    if (x_newton > low + 0.05 * width .and. x_newton < high - 0.05 * width) then
                        x = x_newton
                    else
                        x = x_bisection
                    endif
                else
                    x = x_bisection
                endif
            else
                x = x_bisection
            endif
        enddo

        ! Report error if process did not converge, and output x-value found at last itertaion
        if (.not.converged) then
            write(*, '(a, 4x, a, f6.4, a, f6.4)') "[Warning][minimise_chi_squared]", &
            "Minimum of correlated chi squared could not be found inside the bracket: [", &
            mass_bracket(1), ", ", mass_bracket(2)
            mass = x
        endif

        ! Calculate correlated chi-squared and report results
        chi_squared = get_chi_squared(t_min, t_max, mass, eigenvalues, eigenvectors, transformed_data)

        ! Output info status if requested
        if (present(info)) then
            if (converged) then
                info = 0
            else
                info = 1
            endif
        endif
    end subroutine

    ! Fit a given correlation function over a given t-window to a cosh function by finding the global minimum of the correlated chi-squared function over a mass range between m=0 and a supplied mass upper bound, which should be as a result of the variational calculation
    subroutine fit_cfn(t_min, t_max, cfn, covariance_matrix, mass_upper_bound, &
        mass, chi_squared, course_scan_size)
        implicit none

        integer, intent(in) :: t_min, t_max
        real(real64), intent(in) :: cfn(0:MAX_DELTA_T), covariance_matrix(0:MAX_DELTA_T, 0:MAX_DELTA_T), &
        mass_upper_bound
        real(real64), intent(out) :: mass, chi_squared
        integer, intent(in), optional :: course_scan_size

        real(real64) :: eigenvalues(t_min:t_max), eigenvectors(t_min:t_max, t_min:t_max), &
        transformed_data(t_min:t_max)
        real(real64), allocatable :: mass_scan(:), F_scan(:,:), mass_brackets(:,:), &
        mass_values(:), chi_squared_values(:)
        integer, allocatable :: bracket_indices(:), bracket_steps(:), error_codes(:)
        integer :: course_scan_size_, t_min_, t_max_, num_minima, i, k, mass_loc, step

        ! Set default course scan size
        if (present(course_scan_size)) then
            course_scan_size_ = course_scan_size
        else
            course_scan_size_ = 200
        endif

        ! Validate input parameters
        ! Set minimum course scan size
        if (course_scan_size_ < 20) then
            course_scan_size_ = 20
        endif
        ! Check t_min >= 0
        t_min_ = t_min
        if (t_min_ < 0) then
            t_min_ = 0
        endif
        ! Check t_max <= MAX_DELTA_T
        t_max_ = t_max
        if (t_max_ > MAX_DELTA_T) then
            t_max_ = MAX_DELTA_T
        endif
        ! Check t_min < t_max, and that there are at least 3 data points in the fitting window. If not, throw and error.
        if (t_min_ + 2 > t_max_) then
            write(*, '(a, 4x, a)') "[Error][fit_cfn]", &
            "Must be at least 3 data points in fitting window"
            error stop
        endif
        ! Check mass upper bound is finite and strictly positive
        if (ieee_is_finite(mass_upper_bound)) then
            if (mass_upper_bound <= 0.0_real64) then
                write(*, '(a, 4x, a)') "[Error][fit_cfn]", &
                "mass_upper_bound must be positive"
                error stop
            endif
        else
            write(*, '(a, 4x, a)') "[Error][fit_cfn]", &
            "mass_upper_bound must be finite"
            error stop
        endif

        ! Allocate arrays for course scan
        allocate(mass_scan(course_scan_size_), F_scan(2, course_scan_size_))

        ! Get eigenvalues and eigenvectors of covariance matrix
        call diagonalise_covariance_matrix(covariance_matrix, t_min_, t_max_, eigenvalues, eigenvectors, base=0)

        ! Transform correlation function window into eigenbasis of the covariance matrix
        transformed_data = transformed_cfn(t_min_, t_max_, eigenvectors, cfn)

        ! Generate an equally-spaced course grid of mass values in the range [0, m_max]
        mass_scan = linspace(0.0_real64, mass_upper_bound, course_scan_size_)

        ! Evaluate F over all masses in course grid
        do i = 1, course_scan_size_
            call get_F(t_min_, t_max_, mass_scan(i), eigenvalues, eigenvectors, transformed_data, &
            F_scan(1, i), F_scale=F_scan(2, i))
        enddo


        !! Find mass brackets by finding pairs of mass values over which F(m) changes sign
        ! Find indices of mass values i where there is a minimum of chi-squared in the range [m(i), m(i+1)]
        bracket_indices = [0]
        bracket_steps = [0]
        do i = 1, course_scan_size_ - 1
            ! Check if a course grid point has almost exactly found a zero of F
            if (abs(F_scan(1, i)) <  1.0e-10_real64 * F_scan(2, i) .and. i /= 1) then
                ! If signs of F(m_i-1) and F(m_i+1) differ, add this as a bracket
                if ((F_scan(1, i-1) < 0.0_real64 .and. F_scan(1, i+1) > 0.0_real64) .or. &
                    (F_scan(1, i-1) > 0.0_real64 .and. F_scan(1, i+1) < 0.0_real64)) then
                    bracket_indices = [bracket_indices, i-1]
                    bracket_steps = [bracket_steps, 2]
                endif
            endif
            ! Check if interval [m(i), m(i+1)] contians a root of F
            if ((F_scan(1, i) < 0.0_real64 .and. F_scan(1, i+1) > 0.0_real64) .or. &
                (F_scan(1, i) > 0.0_real64 .and. F_scan(1, i+1) < 0.0_real64)) then
                bracket_indices = [bracket_indices, i]
                bracket_steps = [bracket_steps, 1]
            endif
        enddo

        ! Check the algorithm has found at least one minimum
        if (size(bracket_indices) == 1) then
            num_minima = 0
        else
            ! Get number of minima
            bracket_indices = bracket_indices(2:)
            bracket_steps = bracket_steps(2:)
            num_minima = size(bracket_indices)

            ! Evaluate mass brackets that contain chi-squared minima
            allocate(mass_brackets(2, num_minima))
            do k = 1, num_minima
                i = bracket_indices(k)
                step = bracket_steps(k)
                mass_brackets(1, k) = mass_scan(i)
                mass_brackets(2, k) = mass_scan(i+step)
            enddo
        endif


        !! Run a Newton-bisection algorithm on each of these brackets to identify a root of F(m), and hence a minimim of chi-squared
        ! Save mass and chi-squared values of end-points of search space
        allocate(mass_values(0:num_minima+1), chi_squared_values(0:num_minima+1), error_codes(0:num_minima+1))
        ! Starting mass
        mass_values(0) = 0.0_real64
        chi_squared_values(0) &
        = get_chi_squared(t_min_, t_max_, 0.0_real64, eigenvalues, eigenvectors, transformed_data)
        error_codes(0) = 0
        ! Stopping mass
        mass_values(num_minima+1) = mass_upper_bound
        chi_squared_values(num_minima+1) &
        = get_chi_squared(t_min_, t_max_, mass_upper_bound, eigenvalues, eigenvectors, transformed_data)
        error_codes(num_minima+1) = 0

        ! Save mass and chi-squared values of minima of chi-squared by using a Newton-bisection algorithm to solve for roots of F(m) in the intervals determined above. Discard points that did not converge by setting chi-squared value to be very large.
        do k = 1, num_minima
            call minimise_chi_squared(t_min_, t_max_, eigenvalues, eigenvectors, transformed_data, &
            mass_brackets(:,k), mass_values(k), chi_squared_values(k), info=error_codes(k))
            if (error_codes(k) /= 0) then
                chi_squared_values(k) = huge(1.0_real64)
            endif
        enddo


        !! Output the best-fit value of the mass (for this t-window) by selecting the minimum with the minimal chi-squared. Output a warning if this is found to be one of the end-points of the search space.
        ! Find location of minimal chi-squared in list of values
        mass_loc = minloc(chi_squared_values, 1) - 1

        ! Output warning if this is an end-point of the search space
        if (mass_loc == 0 .or. mass_loc == num_minima+1) then
            write(*, '(a, 4x, a)') "[Warning][fit_cfn]", &
            "Minimum chi-squared value found at boundary of input mass range"
        endif

        ! Output mass and chi-squared values at this minimum
        mass = mass_values(mass_loc)
        chi_squared = chi_squared_values(mass_loc)

        ! Output a warning if mass value with minimum chi-squared did not converge
        if (error_codes(mass_loc) /= 0) then
            write(*, '(a, 4x, a)') "[Warning][fit_cfn]", &
            "Mass value found by fit_cfn did not converge"
        endif
    end subroutine
end module


program analyse_lines
    use analysis
    implicit none

    ! ------------------------------------------ Declare variables ------------------------------------------ !

    ! Command line input variables
    character(len=31) :: input_state, arg
    integer :: nargs
    logical :: read_nconfig

    ! File access variables
    character(len=7) :: file_spin, file_q, file_p, file_r
    character(len=255) :: vevs_file_name, corr_matrix_file_name
    logical :: file_exists

    ! Options variables
    integer :: max_dt, SPIN, PP, PR, MOMENTUM, RUNS
    logical :: max_dt_passed, lines_passed, lines_filtered_passed, mode_passed, use_fit_mode, &
    num_states_passed, spin_passed, PP_passed, PR_passed, momentum_passed, runs_passed

    ! Analysis variables
    complex(real64), allocatable :: vevs(:,:), corr_matrix(:,:,:,:), vevs_avg(:), corr_matrix_avg(:,:,:), &
    eigenvectors(:,:,:), eigenvectors_avg(:,:)
    real(real64), allocatable :: mass_upper_bound(:,:), mass_upper_bound_avg(:), &
    cfn(:,:,:), cfn_avg(:,:), cfn_error(:,:), m_eff(:,:), m_eff_error(:,:)

    ! Loop variables
    integer :: i


    ! ------------------------------------------- Set up analysis ------------------------------------------- !

    ! Get number of states to use and state to analyse from command line
    read_nconfig = .true.
    spin_passed = .false.
    PP_passed = .false.
    PR_passed = .false.
    momentum_passed = .false.
    max_dt_passed = .false.
    lines_passed = .false.
    lines_filtered_passed = .false.
    mode_passed = .false.
    use_fit_mode = .true.
    num_states_passed = .false.
    runs_passed = .false.
    nargs = command_argument_count()
    i = 1
    do while (i <= nargs)
        call get_command_argument(i, arg)

        select case (trim(arg))
        case("-h")
            call options_help()
            stop
        case("--help")
            call options_help()
            stop
        case("-J")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) SPIN
            spin_passed = .true.
        case("-q")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) MOMENTUM
            momentum_passed = .true.
        case("-P")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) PP
            PP_passed = .true.
        case("-R")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) PR
            PR_passed = .true.
        case("-n")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) NUM_STATES
            num_states_passed = .true.
        case("--runs")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) RUNS
            runs_passed = .true.
        case("--configs")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) NCONFIG
            read_nconfig = .false.
        case("--max-dt")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) max_dt
            max_dt_passed = .true.
        case("--lines")
            lines_passed = .true.
        case("--lines-filtered")
            lines_filtered_passed = .true.
        case("--fit")
            if (mode_passed) then
                write(*, '(a)') "[Warning]    Effective mass mode ambiguous. Defaulting to numerical fitting."
                mode_passed = .false.
            else
                use_fit_mode = .true.
                mode_passed = .true.
            endif
        case("--exact")
            if (mode_passed) then
                write(*, '(a)') "[Warning]    Effective mass mode ambiguous. Defaulting to numerical fitting."
                mode_passed = .false.
            else
                use_fit_mode = .false.
                mode_passed = .true.
            endif
        case default
            write(*, '(a, a)') "[Error]    Unknown option:    ", trim(arg)
        end select

        i = i + 1
    enddo

    ! Ensure input state exactly determines and valid state
    if (spin_passed .and. momentum_passed) then
        write(file_spin, '(i0)') SPIN
        write(file_q, '(i0)') MOMENTUM

        if ((SPIN == 0 .or. SPIN == 2) .and. MOMENTUM == 0) then
            if (PP_passed .and. PR_passed) then
                select case(PP)
                case(1)
                    file_p = 'P'
                case(-1)
                    file_p = 'M'
                case default
                    write(*, '(a)') "[Error]    Transverse parity must be P=+-1"
                    stop
                end select

                select case(PR)
                case(1)
                    file_r = 'P'
                case(-1)
                    file_r = 'M'
                case default
                    write(*, '(a)') "[Error]    Intrinsic parity must be R=+-1"
                    stop
                end select

                input_state = "J" // trim(file_spin) // trim(file_p) // trim(file_r) // "q" // trim(file_q)
            else
                write(*, '(a)') &
                "[Error]    Both transverse and intrinsic parity must be defined for J=0,2, q=0 states"
                stop
            endif

        elseif (SPIN == 1 .and. MOMENTUM == 0) then
            if (PR_passed) then
                select case(PR)
                case(1)
                    file_r = 'P'
                case(-1)
                    file_r = 'M'
                case default
                    write(*, '(a)') "[Error]    Intrinsic parity must be R=+-1"
                    stop
                end select

                input_state = "J" // trim(file_spin) // trim(file_r) // "q" // trim(file_q)
            else
                write(*, '(a)') &
                "[Error]    Intrinsic parity must be defined for J=1, q=0 states"
                stop
            endif

            if (PP_passed) then
                if (PP /= 1) then
                    write(*, '(a)') "[Warning]    Only P=+ state exists for J=1, q=0"
                endif
            endif

        elseif ((SPIN == 0 .or. SPIN == 2) .and. (MOMENTUM == 1 .or. MOMENTUM == 2)) then
            if (PP_passed) then
                select case(PP)
                case(1)
                    file_p = 'P'
                case(-1)
                    file_p = 'M'
                case default
                    write(*, '(a)') "[Error]    Transverse parity must be P=+-1"
                    stop
                end select

                input_state = "J" // trim(file_spin) // trim(file_p) // "q" // trim(file_q)
            else
                write(*, '(a)') &
                "[Error]    Transverse parity must be defined for J=0,2, q/=0 states"
                stop
            endif

            if (PR_passed) then
                if (PR /= 1) then
                    write(*, '(a)') "[Warning]    Only R=+ state exists for J=0,2, q/=0"
                endif
            endif

        elseif (SPIN == 1 .and. (MOMENTUM == 1 .or. MOMENTUM == 2)) then
            input_state = "J" // trim(file_spin) // "q" // trim(file_q)

            if (PP_passed) then
                if (PP /= 1) then
                    write(*, '(a)') "[Warning]    Only P=+ state exists for J=1, q/=0"
                endif
            endif

            if (PR_passed) then
                if (PR /= 1) then
                    write(*, '(a)') "[Warning]    Only R=+ state exists for J=1, q/=0"
                endif
            endif
        else
            write(*, '(a, a)') "[Error]    Spin can only take values J=0,1,2. ", &
            "Momentum can only take values q=0,1,2."
            stop
        endif
    else
        write(*, '(a)') "[Error]    Spin and momentum value must be passed to define a state"
        stop
    endif

    ! Set default effective mass mode to use numerical fit
    if (.not.mode_passed) use_fit_mode = .true.

    ! Set default number of states (This is changed later to the full spectrum if there are insufficient states that pass through the filters imposed)
    if (.not. num_states_passed) NUM_STATES = 3

    ! Set default number of runs (default to 10)
    if (.not. runs_passed) RUNS = 10

    ! Check vev and correlation matrix files exist
    vevs_file_name = "vevs_" // trim(input_state) // ".dat"
    inquire(file='run1/' // trim(vevs_file_name), exist=file_exists)
    if (.not. file_exists) then
        write(*, '(3(a))') "[Error]    File ", trim(vevs_file_name), " could not be found"
        error stop
    endif
    corr_matrix_file_name = "corr_matrix_" // trim(input_state) // ".dat"
    inquire(file='run1/' // trim(corr_matrix_file_name), exist=file_exists)
    if (.not. file_exists) then
        write(*, '(3(a))') "[Error]    File ", trim(corr_matrix_file_name), " could not be found"
        error stop
    endif


    ! -------------------------------------- Prepare data for analysis -------------------------------------- !

    ! Read data from files and remove duplicate operators
    call read_clean_data(vevs_file_name, corr_matrix_file_name, vevs, corr_matrix, &
    RUNS, read_nconfig, lines_passed)

    ! Set default max_dt
    if (.not. max_dt_passed) max_dt = 7
    if (max_dt < 0 .or. max_dt > MAX_DELTA_T) max_dt = 7

    ! Update N_CONFIG to NCONFIG * RUNS
    NCONFIG = NCONFIG * RUNS

    ! Arrange data into jackkife samples
    call jackknife(vevs, corr_matrix, vevs_avg, corr_matrix_avg)

    ! Form connected correlation matrix
    call subtract_vevs(vevs, corr_matrix, vevs_avg, corr_matrix_avg)

    ! Output diagonal correlators after filtering of operators if requested
    print *, ""
    if (lines_filtered_passed) then
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        print *, ""
        do i = 1, NUM_OPERATORS
            write(*, '(a, i0)') "Operator ", i
            call output_effective_masses(corr_matrix(i, i, :, :)%re, corr_matrix_avg(i, i, :)%re, use_fit_mode)
            print *, ""
        enddo
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
    endif


    ! ----------------------------- Run analysis and output results to console ----------------------------- !

    ! Run GEVP on correlation matrix to find approximate wavefunctionals for states
    call GEVP(corr_matrix, corr_matrix_avg, mass_upper_bound, mass_upper_bound_avg, &
    eigenvectors, eigenvectors_avg)

    ! Construct correlation functions for states
    call construct_correlators(corr_matrix, corr_matrix_avg, eigenvectors, eigenvectors_avg, &
    cfn, cfn_avg, cfn_error)

    ! Calculate effective masses for states
    call get_effective_masses_all_states(cfn, cfn_avg, m_eff, m_eff_error, use_fit_mode)

    ! Output correlation functions and effective masses with errors to console
    call output_results(cfn_avg, cfn_error, m_eff, m_eff_error, &
    mass_upper_bound, mass_upper_bound_avg, max_dt=max_dt)
end program