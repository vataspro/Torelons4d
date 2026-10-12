module analysis
    use iso_fortran_env, only : int8, int32, int64, real32, real64

    implicit none

    integer :: NUM_OPERATORS, MAX_DELTA_T, NUM_BINS, NCONFIG

    interface generalised_eigenpairs
        module procedure generalised_eigenpairs_real_symmetric
        module procedure generalised_eigenpairs_complex_Hermitian
    end interface

    contains

    ! Solve real symmetric generalised eigenvalue problem using LAPACK routine DSYGVX. Eigenpairs are returned in ascending order (so eigenvalues(k) refers to the largest eigenvalue)
    subroutine generalised_eigenpairs_real_symmetric(A, B, eigenvalues, eigenvectors)
        implicit none

        real(real64), intent(in) :: A(:,:), B(:,:)
        real(real64), intent(out) :: eigenvalues(:), eigenvectors(:,:)

        real(real64), allocatable :: Acopy(:,:), Bcopy(:,:), work(:)
        integer, allocatable :: iwork(:), ifail(:)
        integer :: lwork, info, n, m, il, iu, k
        real(real64) :: abstol

        ! Get size of matrices
        n = size(A, 1)

        ! Get requested number of eigenpairs
        k = size(eigenvalues)

        ! Check dimensions
        if (size(A, 2) /= n .or. size(B, 1) /= n .or. size(B, 2) /= n) &
        error stop "[Error][generalised_eigenpairs]    A and B must be square and of equal size"

        if (k < 1 .or. k > n) &
        error stop "[Error][generalised_eigenpairs]    k must be between 1 and n"

        if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= k) &
        error stop "[Error][generalised_eigenpairs]    eigenvectors must have shape (n,k)"

        ! Allocate arrays
        allocate(Acopy(n,n), Bcopy(n,n))
        allocate(iwork(5*n), ifail(n))

        ! Copy matrices, since DSYGVX overwrites them
        Acopy = A
        Bcopy = B

        ! Select largest k eigenvalues
        il = n - k + 1
        iu = n

        ! Default absolute tolerance
        abstol = 0.0d0

        ! Workspace query
        allocate(work(1))
        lwork = -1

        call DSYGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, iwork, ifail, info)

        if (info /= 0) &
        error stop "[Error][generalised_eigenpairs]    DSYGVX workspace query failed"

        ! Allocate optimal workspace
        lwork = max(8*n, int(work(1)))
        deallocate(work)
        allocate(work(lwork))

        ! Solve generalized eigenvalue problem
        call DSYGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, iwork, ifail, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][generalised_eigenpairs]    DSYGVX: illegal argument ', -info
            else
                write(*,'(a,i0)') '[Error][generalised_eigenpairs]    DSYGVX failed, INFO = ', info
            endif
            error stop
        endif

        if (m /= k) &
        error stop "[Error][generalised_eigenpairs]    DSYGVX returned an unexpected number of eigenpairs"
    end subroutine

    subroutine generalised_eigenpairs_complex_Hermitian(A, B, eigenvalues, eigenvectors)
        implicit none

        complex(real64), intent(in) :: A(:,:), B(:,:)
        real(real64), intent(out) :: eigenvalues(:)
        complex(real64), intent(out) :: eigenvectors(:,:)

        complex(real64), allocatable :: Acopy(:,:), Bcopy(:,:), work(:)
        real(real64), allocatable :: rwork(:)
        integer, allocatable :: iwork(:), ifail(:)

        integer :: lwork, info, n, m, il, iu, k
        real(real64) :: abstol

        ! Get size of matrices
        n = size(A, 1)

        ! Get requested number of eigenpairs
        k = size(eigenvalues)

        ! Check dimensions
        if (size(A, 2) /= n .or. size(B, 1) /= n .or. size(B, 2) /= n) &
        error stop "[Error][generalised_eigenpairs]    A and B must be square and of equal size"

        if (k < 1 .or. k > n) &
        error stop "[Error][generalised_eigenpairs]    k must be between 1 and n"

        if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= k) &
        error stop "[Error][generalised_eigenpairs]    eigenvectors must have shape (n,k)"

        ! Allocate arrays
        allocate(Acopy(n,n), Bcopy(n,n))
        allocate(rwork(7*n), iwork(5*n), ifail(n))

        ! Copy matrices, since ZHEGVX overwrites them
        Acopy = A
        Bcopy = B

        ! Select largest k eigenvalues
        il = n - k + 1
        iu = n

        ! Default absolute tolerance
        abstol = 0.0d0

        ! Workspace query
        allocate(work(1))
        lwork = -1

        call ZHEGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, rwork, iwork, ifail, info)

        if (info /= 0) &
        error stop "ZHEGVX workspace query failed"

        ! Allocate optimal workspace
        lwork = max(2*n, int(real(work(1))))
        deallocate(work)
        allocate(work(lwork))

        ! Solve generalized eigenvalue problem
        call ZHEGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, rwork, iwork, ifail, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][generalised_eigenpairs]    ZHEGVX: illegal argument ', -info
            else if (info <= n) then
                write(*,'(a,i0)') &
                '[Error][generalised_eigenpairs]    ZHEGVX: eigenvectors failed to converge, INFO = ', info
            else
                write(*,'(a,i0)') &
                '[Error][generalised_eigenpairs]    ZHEGVX: B is not positive definite, INFO = ', info
            endif
            error stop
        endif

        if (m /= k) &
        error stop "[Error][generalised_eigenpairs]    ZHEGVX returned an unexpected number of eigenpairs"
    end subroutine

    ! Read vevs and correlation matrix into memory, and remove duplicate operators and zero operators
    subroutine read_clean_data(vevs_filename, corr_matrix_filename, vevs, corr_matrix)
        implicit none
        character(len=255), intent(in) :: vevs_filename, corr_matrix_filename
        complex(real64), allocatable, intent(out) :: vevs(:,:), corr_matrix(:,:,:,:)

        complex(real64), allocatable :: vevs_from_file(:,:), corr_matrix_from_file(:,:,:,:), &
        diagonal_correlators(:)
        integer(int32) :: dim1, dim2, dim3, dim4, i, j, nconfig_check
        integer, allocatable :: distinct_operators(:)
        logical, allocatable :: operator_comparison(:)

        !! Read vevs from file
        ! Open file
        open(11, file=trim(vevs_filename), form='unformatted', access='stream', status='old')

        ! Read number of configurations
        read(11) NCONFIG

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
        read(11) nconfig_check
        if (nconfig_check /= NCONFIG) then
            write(*, '(a)') &
            "[Error]    Number of configurations not consistent between vevs and correlation matrix"
            stop
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
        NUM_OPERATORS = size(vevs, dim=1)

        ! Check there are the same number of operators in corr_matrix as in vevs
        if (size(corr_matrix_from_file, dim=1) /= NUM_OPERATORS) &
        error stop "vevs and corr_matrix incompatible: different number of operators"
        if (size(corr_matrix_from_file, dim=1) /= size(corr_matrix, dim=2)) &
        error stop "corr_matrix not square"

        ! Assume time seperations contained in corr_matrix(i,j,:,l)
        MAX_DELTA_T = size(corr_matrix_from_file, dim=3) - 1

        ! Assume bins contained in vevs(i,:), corr_matrix(i,j,k,:)
        NUM_BINS = size(vevs_from_file, dim=2)

        ! Check there are the same number of bins in corr_matrix as in vevs
        if (size(corr_matrix_from_file, dim=4) /= NUM_BINS) &
        error stop "vevs and corr_matrix incompatible: different number of bins"

        ! Output dimensions found
        write(*, '(a, i0)') "[Info] Number of operators found: ", NUM_OPERATORS
        write(*, '(a, i0)') "[Info] Number time seperations found: ", MAX_DELTA_T + 1
        write(*, '(a, i0)') "[Info] Number of bins found: ", NUM_BINS


        !! Find duplicate operators by comparing diagonal elements of correlation matrix at zero seperation
        ! Compute average diagonal correlation matrix elements at 0 seperation
        allocate(diagonal_correlators(NUM_OPERATORS))
        do i = 1, NUM_OPERATORS
            diagonal_correlators(i) = sum(corr_matrix_from_file(i, i, 0, :))
        enddo

        ! Find non-duplicate operators by comparing diagonal correlators
        distinct_operators = [1]
        do i = 2, NUM_OPERATORS
            allocate(operator_comparison(i-1))
            do j = 1, i-1
                operator_comparison(j) &
                = (abs(diagonal_correlators(j) - diagonal_correlators(i)) <= epsilon(1.0d0))
            enddo
            if (.not.any(operator_comparison)) distinct_operators = [distinct_operators, i]
        enddo

        !! Fill correlation matrix with distinct operators only
        ! Reset number of operators to number of distinct operators
        NUM_OPERATORS = size(distinct_operators)

        ! Allocate vevs and correlation matrix outputs
        allocate(vevs(NUM_OPERATORS, NUM_BINS), &
        corr_matrix(NUM_OPERATORS, NUM_OPERATORS, 0:MAX_DELTA_T, NUM_BINS))

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
    subroutine jacknife(vevs, corr_matrix, vevs_avg, corr_matrix_avg)
        implicit none
        complex(real64), intent(inout) :: vevs(:,:), corr_matrix(:,:,:,:)
        complex(real64), intent(out) :: vevs_avg(:), corr_matrix_avg(:,:,:)
        integer :: i, j, dt

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
            corr_matrix_avg(i, j, :) = corr_matrix_avg(i, j, :) - vevs_avg(i) * vevs_avg(j)
        enddo

        ! Subtract vevs from jacknife bins of correaltion matrices
        do concurrent (i = 1:NUM_OPERATORS, j = 1:NUM_OPERATORS, bin = 1:NUM_BINS)
            corr_matrix(i, j, :, bin) = corr_matrix(i, j, :, bin) - vevs(i, bin) * vevs(j, bin)
        enddo
    end subroutine

    ! Run GEVP on correlation matrix

    ! Calculate effective masses
end module