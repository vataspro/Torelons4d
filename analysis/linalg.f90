module linalg
    use iso_fortran_env, only : int8, int32, int64, real32, real64
    implicit none

    interface eigenpairs
        module procedure eigenpairs_real_symmetric
    end interface

    interface generalised_eigenpairs
        module procedure generalised_eigenpairs_real_symmetric
        module procedure generalised_eigenpairs_complex_Hermitian
    end interface

    contains

    ! Return a list of evenly spaced real numbers in the interval [start, stop]
    function linspace(start, stop, num)
        implicit none

        real(real64), intent(in) :: start, stop
        integer, intent(in) :: num

        real(real64) :: linspace(num)

        real(real64) :: step
        integer :: i

        step = (stop - start) / real(num - 1, kind=real64)
        do i = 1, num
            linspace(i) = start + real(i - 1, kind=real64) * step
        enddo
    end function

    ! Solve real symmetric eigenvalue problem using LAPACK routine DSYEV.
    ! All eigenpairs are returned in ascending order of eigenvalue.
    ! The columns of eigenvectors contain the corresponding eigenvectors.
    subroutine eigenpairs_real_symmetric(A, eigenvalues, eigenvectors)
        implicit none

        real(real64), intent(in)  :: A(:,:)
        real(real64), intent(out) :: eigenvalues(:), eigenvectors(:,:)

        real(real64), allocatable :: Acopy(:,:), work(:)
        integer :: lwork, info, n

        ! Get size of matrix
        n = size(A, 1)

        ! Check dimensions
        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][eigenpairs_real_symmetric]    A must be a non-empty square matrix"

        if (size(eigenvalues) /= n) &
        error stop "[Error][eigenpairs_real_symmetric]    eigenvalues must have shape (n)"

        if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n) &
        error stop "[Error][eigenpairs_real_symmetric]    eigenvectors must have shape (n,n)"

        ! Copy matrix, since DSYEV overwrites it
        allocate(Acopy(n,n))
        Acopy = A

        ! Workspace query
        allocate(work(1))
        lwork = -1

        call DSYEV('V', 'U', n, Acopy, n, eigenvalues, work, lwork, info)

        if (info /= 0) &
        error stop "[Error][eigenpairs_real_symmetric]    DSYEV workspace query failed"

        ! Allocate optimal workspace
        lwork = max(1, 3*n-1, int(work(1)))
        deallocate(work)
        allocate(work(lwork))

        ! Restore matrix before solving
        Acopy = A

        ! Solve symmetric eigenvalue problem
        call DSYEV('V', 'U', n, Acopy, n, eigenvalues, work, lwork, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') &
                    '[Error][eigenpairs_real_symmetric]    DSYEV: illegal argument ', -info
            else
                write(*,'(a,i0)') &
                    '[Error][eigenpairs_real_symmetric]    DSYEV failed to converge, INFO = ', info
            endif
            error stop
        endif

        ! DSYEV returns eigenvectors as columns of Acopy
        eigenvectors = Acopy

    end subroutine eigenpairs_real_symmetric


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

        complex(8), intent(in) :: A(:,:), B(:,:)
        real(real64), intent(out) :: eigenvalues(:)
        complex(8), intent(out) :: eigenvectors(:,:)

        complex(8), allocatable :: Acopy(:,:), Bcopy(:,:), work(:)
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
end module