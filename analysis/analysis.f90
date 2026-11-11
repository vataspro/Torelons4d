module analysis
    use iso_fortran_env, only : int8, int32, int64, real32, real64
    use linalg

    implicit none

    integer :: NUM_OPERATORS, MAX_DELTA_T, NUM_BINS, NCONFIG, NUM_STATES

    contains

    ! Provide help with options
    subroutine options_help()
        write(*, '(a)') "[Help]    Usage: ./analysis <options>"
        write(*, '(a)') "[Help]"
        write(*, '(a)') "[Help]    Required options for q=0 states:"

        write(*, '(a, 11x, a)') "[Help]     * '-par <parallel parity>'", &
        "Parallel parity of flux tube state for masses to be calculated"

        write(*, '(a, 6x, a)') "[Help]     * '-perp <perpendicular parity>'", &
        "Perpendicular parity of flux tube state for masses to be calculated"

        write(*, '(a)') "[Help]"
        write(*, '(a)') "[Help]    Required options for q/=0 states:"

        write(*, '(a, 22x, a)') "[Help]     * '-P <parity>'", &
        "Parity of flux tube state for masses to be calculated"

        write(*, '(a, 8x, a)') "[Help]     * '-q <longitudinal momentum>'", &
        "Longitudinal momentum of flux tube state for masses to be calculated"

        write(*, '(a)') "[Help]"
        write(*, '(a)') "[Help]    Other passable options:"

        write(*, '(a, 9x, a)') "[Help]     * '-n <number_of_states>'", &
        "Total of number_of_states calculated in GEVP. Defaults to n = 3 if not passed."

        write(*, '(a, 4x, a)') "[Help]     * '--max-dt' <maximum delta t>", &
        "Report effective masses upto and including maximum delta t. Defaults to max_dt = 7 if not passed."

        write(*, '(a, 23x, a, a)') "[Help]     * '--lines'", &
        "Show correlation functions and effective masses for individual line operators. ", &
        "Defaults to FALSE if not passed."

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
    subroutine read_clean_data(vevs_filename, corr_matrix_filename, vevs, corr_matrix, read_nconfig)
        implicit none
        character(len=*), intent(in) :: vevs_filename, corr_matrix_filename
        complex(real64), allocatable, intent(out) :: vevs(:,:), corr_matrix(:,:,:,:)
        logical, intent(in), optional :: read_nconfig

        complex(real64), allocatable :: vevs_from_file(:,:), corr_matrix_from_file(:,:,:,:), &
        diagonal_correlators(:)
        integer(int32) :: dim1, dim2, dim3, dim4, i, j, nconfig_check
        integer, allocatable :: distinct_operators(:)
        logical, allocatable :: operator_comparison(:)
        logical :: read_nconfig_default

        ! Set read_nconfig
        if (present(read_nconfig)) then
            read_nconfig_default = read_nconfig
        else
            read_nconfig_default = .true.
        endif

        !! Read vevs from file
        ! Open file
        open(11, file=trim(vevs_filename), form='unformatted', access='stream', status='old')

        ! Read number of configurations
        if (read_nconfig_default) then
            read(11) NCONFIG
        endif

        ! Read dimensions of array
        read(11) dim1
        read(11) dim2

        ! Allocate vevs array
        allocate(vevs_from_file(dim1, dim2))

        ! Read vevs data from file
        read(11) vevs_from_file

        ! Close file
        close(11)


        !! Read correlation matrix from file
        ! Open file
        open(11, file=trim(corr_matrix_filename), form='unformatted', access='stream', status='old')

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
        read(11) dim1
        read(11) dim2
        read(11) dim3
        read(11) dim4

        ! Allocate vevs array
        allocate(corr_matrix_from_file(dim1, dim2, dim3, dim4))

        ! Read vevs data from file
        read(11) corr_matrix_from_file

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

        ! Assume bins contained in vevs(i,:), corr_matrix(i,j,k,:)
        NUM_BINS = size(vevs_from_file, dim=2)

        ! Check there are the same number of bins in corr_matrix as in vevs
        if (size(corr_matrix_from_file, dim=4) /= NUM_BINS) &
        error stop "vevs and corr_matrix incompatible: different number of bins"

        ! Output dimensions found
        write(*, '(a, i0)') "[Info]    Number of operators from file: ", NUM_OPERATORS
        write(*, '(a, i0)') "[Info]    Number time seperations from file: ", MAX_DELTA_T + 1
        write(*, '(a, i0)') "[Info]    Number of bins from file: ", NUM_BINS


        !! Find duplicate operators by comparing diagonal elements of correlation matrix at zero seperation
        ! Compute average diagonal correlation matrix elements at 0 seperation
        allocate(diagonal_correlators(NUM_OPERATORS))
        do i = 1, NUM_OPERATORS
            diagonal_correlators(i) = sum(corr_matrix_from_file(i, i, 1, :))
        enddo

        ! Find non-duplicate operators by comparing diagonal correlators
        distinct_operators = [1]
        do i = 2, NUM_OPERATORS
            ! Check if operator is exactly 0
            if (abs(diagonal_correlators(i)) <= epsilon(1.0d0)) then
                write(*, '(a, i0, a, f12.10)') "[Info]    Operator ", i, " diagonal correlator is ", &
                diagonal_correlators(i)%re
                cycle
            else
                allocate(operator_comparison(i-1))
                do j = 1, i-1
                    operator_comparison(j) &
                    = (1.0 - abs(sum(corr_matrix_from_file(i,j,1,:))) &
                    / sqrt(diagonal_correlators(i)%re * diagonal_correlators(j)%re) &
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
    subroutine GEVP(corr_matrix_avg, mass_bound, eigenvectors)
        implicit none
        complex(real64), intent(in) :: corr_matrix_avg(:,:,:)
        real(real64), intent(out), allocatable :: mass_bound(:)
        complex(real64), intent(out), allocatable :: eigenvectors(:,:)

        real(real64), allocatable :: eigenvalues(:)

        ! Check dimensions
        if (size(corr_matrix_avg, dim=1) /= NUM_OPERATORS .or. size(corr_matrix_avg, dim=2) /= NUM_OPERATORS) &
        error stop "[Error][GEVP]    corr_matrix_avg must have NUM_OPERATORS as size of first two dimensions"

        ! Allocate mass_bound and eigen arrays
        allocate(mass_bound(NUM_STATES), eigenvalues(NUM_STATES), eigenvectors(NUM_OPERATORS, NUM_STATES))

        ! Run GEVP on correlation matrix
        call generalised_eigenpairs(corr_matrix_avg(:,:,2), corr_matrix_avg(:,:,1), eigenvalues, eigenvectors)

        ! Calculate upper bounds on massed of states from eigenvalues
        mass_bound = -log(eigenvalues)
    end subroutine

    ! Construct normalised correlation functions of operators from eigenvectors
    subroutine construct_correlators(corr_matrix, corr_matrix_avg, eigenvectors, cfn, cfn_avg, cfn_error)
        implicit none
        complex(real64), intent(in) :: corr_matrix(:,:,:,:), corr_matrix_avg(:,:,:), eigenvectors(:,:)
        real(real64), allocatable, intent(out) :: cfn(:,:,:), cfn_avg(:,:), cfn_error(:,:)

        complex(real64) :: scalar_product
        integer :: state, dt, bin

        ! Check dimensions
        if (all(shape(corr_matrix) /= [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1, NUM_BINS])) then
            write(*, '(a, a)') "[Error][construct_correlators]    corr_matrix must have shape", &
            " [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1, NUM_BINS]"
            error stop
        endif

        if (all(shape(corr_matrix_avg) /= [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1])) then
            write(*, '(a, a)') "[Error][construct_correlators]    corr_matrix_avg must have shape", &
            " [NUM_OPERATORS, NUM_OPERATORS, MAX_DELTA_T+1]"
            error stop
        endif

        if (all(shape(eigenvectors) /= [NUM_OPERATORS, NUM_STATES])) &
        error stop &
        "[Error][construct_correlators]    eigenvectors must have shape [NUM_OPERATORS, NUM_STATES]"

        ! Allocate correlation function and average correlation funcion arrays
        allocate(cfn(MAX_DELTA_T+1, NUM_BINS, NUM_STATES), cfn_avg(MAX_DELTA_T+1, NUM_STATES), &
        cfn_error(MAX_DELTA_T+1, NUM_STATES))

        ! Construct correlation functions of states by computing the scalar product v^dagger*C(t)*v
        do state = 1, NUM_STATES
            do dt = 1, MAX_DELTA_T+1
                ! Compute scalar product
                scalar_product = dot_product(eigenvectors(:, state), &
                matmul(corr_matrix_avg(:, :, dt), eigenvectors(:, state)))

                ! Scalar product is real, so the real part is the correlation function entry
                cfn_avg(dt, state) = scalar_product%re
            enddo
        enddo

        do state = 1, NUM_STATES
            do bin = 1, NUM_BINS
                do dt = 1, MAX_DELTA_T+1
                    ! Compute scalar product
                    scalar_product = dot_product(eigenvectors(:, state), &
                    matmul(corr_matrix(:, :, dt, bin), eigenvectors(:, state)))

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
        real(8), intent(in) :: cfn(-1:1)
        real(8) :: m_eff, x

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
        real(8), intent(in) :: cfn(t-1:t)
        real(8), intent(in), optional :: rel_tol, abs_tol
        integer, intent(in), optional :: max_iter
        real(8) :: m_eff

        real(8) :: dt, ratio, rtol, atol, x, x_trial, width, low, high, F, dF
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
        if (all(shape(cfn) /= [MAX_DELTA_T+1, NUM_BINS])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    cfn must have shape", &
            " [MAX_DELTA_T+1, NUM_BINS]"
            error stop
        endif

        if (all(shape(cfn_avg) /= [MAX_DELTA_T+1])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    cfn_avg must have shape", &
            " [MAX_DELTA_T+1]"
            error stop
        endif

        if (all(shape(effective_mass_avg) /= [MAX_DELTA_T])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    effective_mass_avg must have shape", &
            " [MAX_DELTA_T]"
            error stop
        endif

        if (all(shape(effective_mass_error) /= [MAX_DELTA_T])) then
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
        if (all(shape(cfn) /= [MAX_DELTA_T+1, NUM_BINS, NUM_STATES])) then
            write(*, '(a, a)') "[Error][get_effective_masses]    cfn must have shape", &
            " [MAX_DELTA_T+1, NUM_BINS, NUM_STATES]"
            error stop
        endif

        if (all(shape(cfn_avg) /= [MAX_DELTA_T+1, NUM_STATES])) then
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
        real(8), intent(in) :: cfn(:,:), cfn_avg(:)
        logical, intent(in), optional :: fit_mode

        real(8), allocatable :: effective_mass_avg(:), effective_mass_errors(:), &
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
    subroutine output_results(cfn_avg, cfn_error, effective_mass_avg, effective_mass_error, mass_bound, max_dt)
        real(real64), intent(in) :: cfn_avg(:,:), cfn_error(:,:), &
        effective_mass_avg(:,:), effective_mass_error(:,:)
        real(real64), intent(in), optional :: mass_bound(:)
        integer, intent(in), optional :: max_dt

        integer :: i, state, dt, maximum_seperation

        ! Set default max_dt to 7
        if (present(max_dt)) then
            maximum_seperation = max_dt
        else
            maximum_seperation = 7
        endif

        ! Check dimensions
        if (size(mass_bound) /= NUM_STATES) then
            write(*, '(a)') "[Error][output_results]    mass_bound must have size NUM_STATES"
            error stop
        endif

        if (all(shape(cfn_avg) /= [MAX_DELTA_T+1, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    cfn_avg must have shape", &
            " [MAX_DELTA_T+1, k]"
            error stop
        endif

        if (all(shape(cfn_error) /= [MAX_DELTA_T+1, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    cfn_error must have shape", &
            " [MAX_DELTA_T+1, k]"
            error stop
        endif

        if (all(shape(effective_mass_avg) /= [MAX_DELTA_T, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    effective_mass_avg must have shape", &
            " [MAX_DELTA_T, k]"
            error stop
        endif

        if (all(shape(effective_mass_error) /= [MAX_DELTA_T, NUM_STATES])) then
            write(*, '(a, a)') "[Error][output_results]    effective_mass_error must have shape", &
            " [MAX_DELTA_T, k]"
            error stop
        endif

        ! Output correlation functions and effective masses fo each state in turn. Iterate through arrays in reverse order, as the lightest state is the last in each array
        write(*, '(a)') " *******************************************************"
        state = 0
        do i = NUM_STATES, 1, -1
            ! Output state identifier
            write(*, '(a)') " *"
            if (state == 0) then
                write(*, '(a)') "Ground state"
            else
                write(*, '(a, i0)') "Excited state ", state
            endif

            ! Output mass upper bound if supplied
            if (present(mass_bound)) then
                write(*, '(a, f6.4)') "Mass upper bound:    ", mass_bound(i)
            endif
            write(*, '(a)') " *"

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

            write(*, '(a)') " *"
            write(*, '(a)') " *******************************************************"
        enddo
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
    character(len=255) :: vevs_file_name, corr_matrix_file_name
    logical :: file_exists

    ! Options variables
    integer :: max_dt
    logical :: max_dt_passed, lines_passed, mode_passed, use_fit_mode, num_states_passed

    ! Analysis variables
    complex(real64), allocatable :: vevs(:,:), corr_matrix(:,:,:,:), vevs_avg(:), corr_matrix_avg(:,:,:), &
    eigenvectors(:,:)
    real(real64), allocatable :: mass_upper_bound(:), cfn(:,:,:), cfn_avg(:,:), cfn_error(:,:), m_eff(:,:), m_eff_error(:,:)

    ! Loop variables
    integer :: i


    ! ------------------------------------------- Set up analysis ------------------------------------------- !

    ! Get number of states to use and state to analyse from command line
    read_nconfig = .true.
    max_dt_passed = .false.
    lines_passed = .false.
    mode_passed = .false.
    use_fit_mode = .true.
    num_states_passed = .false.
    nargs = command_argument_count()
    i = 1
    do while (i <= nargs)
        call get_command_argument(i, arg)

        select case (trim(arg))
        case("-n")
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) NUM_STATES
            num_states_passed = .true.
        case("--state")
            i = i + 1
            call get_command_argument(i, arg)
            input_state = trim(arg)
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

    ! Set default effective mass mode to use numerical fit
    if (.not.mode_passed) use_fit_mode = .true.

    ! Set default number of states (This is changed later to the full spectrum if there are insufficient states that pass through the filters imposed)
    if (.not. num_states_passed) NUM_STATES = 3

    ! Check vev and correlation matrix files exist
    vevs_file_name = "vevs_" // trim(input_state) // ".dat"
    inquire(file=trim(vevs_file_name), exist=file_exists)
    if (.not. file_exists) then
        write(*, '(3(a))') "File ", trim(vevs_file_name), " could not be found"
        error stop
    endif
    corr_matrix_file_name = "corr_matrix_" // trim(input_state) // ".dat"
    inquire(file=trim(corr_matrix_file_name), exist=file_exists)
    if (.not. file_exists) then
        write(*, '(3(a))') "File ", trim(corr_matrix_file_name), " could not be found"
        error stop
    endif


    ! -------------------------------------- Prepare data for analysis -------------------------------------- !

    ! Read data from files and remove duplicate operators
    call read_clean_data(vevs_file_name, corr_matrix_file_name, vevs, corr_matrix, read_nconfig)

    ! Set default max_dt
    if (.not. max_dt_passed) max_dt = 7
    if (max_dt < 0 .or. max_dt > MAX_DELTA_T) max_dt = 7

    ! Arrange data into jackkife samples
    call jackknife(vevs, corr_matrix, vevs_avg, corr_matrix_avg)

    ! Form connected correlation matrix
    call subtract_vevs(vevs, corr_matrix, vevs_avg, corr_matrix_avg)


    ! ----------------------------- Run analysis and output results to console ----------------------------- !

    if (lines_passed) then
        do i = 1, NUM_OPERATORS
            write(*, '(a, i0)') "Operator ", i
            call output_effective_masses(corr_matrix(i, i, :, :)%re, corr_matrix_avg(i, i, :)%re, use_fit_mode)
            print *, ""
        enddo
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        print *, ""
    endif

    ! Run GEVP on correlation matrix to find approximate wavefunctionals for states
    call GEVP(corr_matrix_avg, mass_upper_bound, eigenvectors)

    ! Construct correlation functions for states
    call construct_correlators(corr_matrix, corr_matrix_avg, eigenvectors, cfn, cfn_avg, cfn_error)

    ! Calculate effective masses for states
    call get_effective_masses_all_states(cfn, cfn_avg, m_eff, m_eff_error, use_fit_mode)

    ! Output correlation functions and effective masses with errors to console
    call output_results(cfn_avg, cfn_error, m_eff, m_eff_error, mass_bound=mass_upper_bound, max_dt=max_dt)
end program