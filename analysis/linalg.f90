module linalg
    implicit none

    interface generalised_eigenpairs
        module procedure generalised_eigenpairs_real_symmetric
        module procedure generalised_eigenpairs_complex_Hermitian
    end interface

    contains

    ! Solve real symmetric generalised eigenvalue problem using LAPACK routine DSYGVX. Eigenpairs are returned in ascending order (so eigenvalues(k) refers to the largest eigenvalue)
    subroutine generalised_eigenpairs_real_symmetric(A, B, eigenvalues, eigenvectors)
        implicit none

        real(8), intent(in) :: A(:,:), B(:,:)
        real(8), intent(out) :: eigenvalues(:), eigenvectors(:,:)

        real(8), allocatable :: Acopy(:,:), Bcopy(:,:), work(:)
        integer, allocatable :: iwork(:), ifail(:)
        integer :: lwork, info, n, m, il, iu, k
        real(8) :: abstol

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
        real(8), intent(out) :: eigenvalues(:)
        complex(8), intent(out) :: eigenvectors(:,:)

        complex(8), allocatable :: Acopy(:,:), Bcopy(:,:), work(:)
        real(8), allocatable :: rwork(:)
        integer, allocatable :: iwork(:), ifail(:)

        integer :: lwork, info, n, m, il, iu, k
        real(8) :: abstol

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