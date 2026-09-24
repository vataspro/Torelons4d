module lattice
    use parameters
    implicit none

    ! Containers for lattice pointers by 1 lattice step and multiple lattice steps in the case of blocked links
    integer, allocatable :: lattice_pointers_up(:,:), lattice_pointers_down(:,:), &
    lattice_pointers_blocked_up(:,:,:), lattice_pointers_blocked_down(:,:,:)

    interface herm
        module procedure herm_real32
        module procedure herm_real64
    end interface

    contains

    ! Return the site index, given a set of coordinates x, y, z, t
    function site_index(x, y, z, t)
        implicit none
        integer, intent(in) :: x, y, z
        integer, intent(in), optional :: t
        integer :: site_index

        if (present(t)) then
            site_index = (t-1) * SLICE_VOLUME + (z-1) * LX1 * LX2 + (y-1) * LX1 + x
        else
            site_index = (z-1) * LX1 * LX2 + (y-1) * LX1 + x
        endif
    end function

    ! Set up lattice pointers and blocked lattice pointers
    subroutine setup_lattice()
        implicit none

        integer :: n, neighbour, LX12, x, y, z, t, blocking_level, direction

        ! Allocate lattice points with pattern (site, direction)
        allocate(lattice_pointers_up(LATTICE_VOLUME, 4), lattice_pointers_down(LATTICE_VOLUME, 4))

        ! Allocate blocked lattice pointers with pattern (site, direction, blocking level)
        allocate(lattice_pointers_blocked_up(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1), &
        lattice_pointers_blocked_down(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1))

        ! Get lattice size of x-y slice
        LX12 = LX1 * LX2

        ! Calculate indices of lattice points in all directions
        n = 0
        do t = 1, LX4
            do z = 1, LX3
                do y = 1, LX2
                    do x = 1, LX1
                        n = n + 1

                        ! Moving in x-direction
                        ! Move up
                        neighbour = n + 1
                        if ((x+1) > LX1) neighbour = neighbour - LX1
                        lattice_pointers_up(n, 1) = neighbour
                        ! Move down
                        neighbour = n - 1
                        if ((x-1) < 1) neighbour = neighbour + LX1
                        lattice_pointers_down(n, 1) = neighbour

                        ! Moving in y-direction
                        ! Move up
                        neighbour = n + LX1
                        if ((y+1) > LX2) neighbour = neighbour - LX12
                        lattice_pointers_up(n, 2) = neighbour
                        ! Move down
                        neighbour = n - LX1
                        if ((y-1) < 1) neighbour = neighbour + LX12
                        lattice_pointers_down(n, 2) = neighbour

                        ! Moving in z-direction
                        ! Move up
                        neighbour = n + LX12
                        if ((z+1) > LX3) neighbour = neighbour - SLICE_VOLUME
                        lattice_pointers_up(n, 3) = neighbour
                        ! Move down
                        neighbour = n - LX12
                        if ((z-1) < 1) neighbour = neighbour + SLICE_VOLUME
                        lattice_pointers_down(n, 3) = neighbour

                        ! Moving in t-direction
                        ! Move up
                        neighbour = n + SLICE_VOLUME
                        if ((t+1 > LX4)) neighbour = neighbour - LATTICE_VOLUME
                        lattice_pointers_up(n, 4) = neighbour
                        ! Move down
                        neighbour = n - SLICE_VOLUME
                        if ((t-1) < 1) neighbour = neighbour + LATTICE_VOLUME
                        lattice_pointers_down(n, 4) = neighbour
                    enddo
                enddo
            enddo
        enddo

        !! Calculate indices of lattice points in all directions on blocked lattice
        ! Set blocking level 1 as equal to the original lattice
        lattice_pointers_blocked_up(:,:,1) = lattice_pointers_up(1:SLICE_VOLUME, 1:3)
        lattice_pointers_blocked_down(:,:,1) = lattice_pointers_down(1:SLICE_VOLUME, 1:3)

        ! For each subsequent blocking level, move twice the distance as the previous blocking level in all directions from all sites
        do blocking_level = 2, MAX_BLOCKING_LEVEL+1
            do direction = 1, 3
                do n = 1, SLICE_VOLUME
                    ! Move up by 2 (blocked) units
                    lattice_pointers_blocked_up(n, direction, blocking_level) &
                    = lattice_pointers_blocked_up(&
                    lattice_pointers_blocked_up(n, direction, blocking_level-1), &
                    direction, blocking_level-1)

                    ! Move down by 2 (blocked) units
                    lattice_pointers_blocked_down(n, direction, blocking_level) &
                    = lattice_pointers_blocked_down(&
                    lattice_pointers_blocked_down(n, direction, blocking_level-1), &
                    direction, blocking_level-1)
                enddo
            enddo
        enddo
    end subroutine setup_lattice

    ! Move around lattice by either blocked or unblocked amounts
    integer function move(site, direction, blocking_level)
        implicit none
        integer, intent(in) :: site, direction
        integer, intent(in), optional :: blocking_level

        ! Check parameters are compatible and fall within the range of the lattice, then return requested lattice site
        if (present(blocking_level)) then
            ! Site must belong to a time slice and direction must point in a spatial direction
            if (site > SLICE_VOLUME) then
                write(*, '(a, a, i0)') "site must be on a spatial slice in move function with blocking_level", &
                "argument present, therefore must be less than or equal to ", SLICE_VOLUME
            endif
            if (abs(direction) > 3 .or. direction == 0) then
                write(*, '(a, a)') "direction must be spatial in move function with blocking_level argument", &
                "present, therefore must be less than or equal to 3"
                stop
            endif
            if (blocking_level > MAX_BLOCKING_LEVEL+1) then
                write(*, '(a, i0)') "blocking_level must be less than ", MAX_BLOCKING_LEVEL+1
                stop
            endif

            ! Return blocked neighbour
            if (direction > 0) then
                move = lattice_pointers_blocked_up(site, direction, blocking_level)
            else
                move = lattice_pointers_blocked_down(site, abs(direction), blocking_level)
            endif
        else
            ! Site must belong to the lattice and direction must point in 4d
            if (site > LATTICE_VOLUME) then
                write(*, '(a, a, i0)') "site must be on the lattice in move function without ", &
                "blocking_level argument present, therefore must be less than or equal to ", LATTICE_VOLUME
            endif
            if (abs(direction) > 4 .or. direction  == 0) then
                write(*, '(a, a)') "direction must be spatial or temporal in move function without ", &
                "blocking_level argument present, therefore must be less than or equal to 4"
                stop
            endif

            ! Return unblocked neighbour
            if (direction > 0) then
                move = lattice_pointers_up(site, direction)
            else
                move = lattice_pointers_down(site, abs(direction))
            endif
        endif
    end function move

    ! Compute hermitian conjugate or matrix
    function herm_real32(A)
        implicit none
        complex(real32), intent(in) :: A(:,:)
        complex(real32), allocatable :: herm_real32(:,:)
        integer :: n

        n = size(A, 1)
        if (size(A, 2) /= n) then
            stop "Input matrix to herm function must be square"
        endif

        herm_real32 = conjg(transpose(A))
    end function

    function herm_real64(A)
        implicit none
        complex(real64), intent(in) :: A(:,:)
        complex(real64), allocatable :: herm_real64(:,:)
        integer :: n

        n = size(A, 1)
        if (size(A, 2) /= n) then
            stop "Input matrix to herm function must be square"
        endif

        herm_real64 = conjg(transpose(A))
    end function

    ! Frobenius norm of off-diagonal elements
    function off(A)
        implicit none

        complex(real64), dimension(NCOL,NCOL), intent(in) :: A
        real(real64) :: off

        integer :: i, j

        off = 0.0
        do j = 1, NCOL
            do i = 1, j-1
                off = off + abs(A(i,j))**2
            enddo
            do i = j+1, NCOL
                off = off + abs(A(i,j))**2
            enddo
        enddo
        off = sqrt(off)
    end function

    ! Find unitarisation of matrix A
    function unitarise_SVD(A) result(U)
        implicit none

        complex(real64), dimension(NCOL,NCOL), intent(in) :: A
        complex(real64), dimension(NCOL,NCOL) :: U

        complex(real64), dimension(NCOL,NCOL) :: A_dagger_A, V
        real(real64), dimension(NCOL) :: Sigma_inv
        complex(real64), dimension(2,NCOL) :: new_rows
        complex(real64), dimension(NCOL,2) :: new_cols
        complex(real64), dimension(2,2) :: pq_block
        real(real64) :: delta, c, t, tau
        complex(real64) :: s
        integer, dimension(2) :: index_pair
        logical, dimension(NCOL,NCOL) :: tri_mask
        integer :: i, j, p, q, iter

        !! Calculate A dagger A
        A_dagger_A = matmul(herm(A), A)

        !! Carry out Jacobi iteration on A_dagger_A
        ! Set tolerance parameter as delta = tol_SVD * ||A_dagger_A||
        delta = 0.0
        do j = 1, NCOL
            do i = 1, NCOL
                delta = delta + abs(A_dagger_A(i,j))**2
            enddo
        enddo
        delta = sqrt(delta)
        delta = TOL_SVD*delta

        ! Initialise matrix of eigenvectors, V
        V = cmplx(0.0,0.0)
        do i = 1, NCOL
            V(i,i) = cmplx(1.0,0.0)
        enddo

        ! Define mask to be true only if i < j
        do j = 1, NCOL
            do i = 1, NCOL
                tri_mask(i,j) = (i < j)
            enddo
        enddo

        ! Carry out Jacobi algorithm until convergence is reached
        iter = 0
        do while (off(A_dagger_A) > delta)
            !! Choose (p,q) so that |a_pq| = max(i/=j) |a_ij|
            index_pair = maxloc(abs(A_dagger_A), mask=tri_mask)
            p = index_pair(1)
            q = index_pair(2)

            !! Find c, s which make up Jacobi rotation matrix
            if (abs(A_dagger_A(p,q)) < max(delta, epsilon(1.0d0)*maxval(abs(A_dagger_A)))) then
                ! If target is already diagonal, apply no rotation
                cycle
            else
                ! Solve for t = tan(theta)
                tau = (real(A_dagger_A(q,q), kind=real64) - real(A_dagger_A(p,p), kind=real64)) &
                / (2*abs(A_dagger_A(p,q)))
                if (tau >= 0) then
                    t = 1.0/(tau + sqrt(1 + tau**2))
                else
                    t = 1.0/(tau - sqrt(1 + tau**2))
                endif

                ! Find c, s from t
                c = 1.0/sqrt(1 + t**2)
                s = t*c*A_dagger_A(p,q)/abs(A_dagger_A(p,q))
            endif

            !! Update A_dagger_A using Jacobi rotation matrix
            ! Find new p'th and q'th rows
            new_rows(1,:) = c*A_dagger_A(p,:) - s*A_dagger_A(q,:)
            new_rows(2,:) = conjg(s)*A_dagger_A(p,:) + c*A_dagger_A(q,:)

            ! Find new p'th and q'th columns
            new_cols(:,1) = c*A_dagger_A(:,p) - conjg(s)*A_dagger_A(:,q)
            new_cols(:,2) = s*A_dagger_A(:,p) + c*A_dagger_A(:,q)

            ! Find new (p,q) block
            pq_block(1,1) = c*new_cols(p,1) - s*new_cols(q,1)
            pq_block(2,1) = c*new_cols(q,1) + conjg(s)*new_cols(p,1)
            pq_block(1,2) = c*new_cols(p,2) - s*new_cols(q,2)
            pq_block(2,2) = c*new_cols(q,2) + conjg(s)*new_cols(p,2)

            ! Update A_dagger_A
            A_dagger_A(p,:) = new_rows(1,:)
            A_dagger_A(q,:) = new_rows(2,:)
            A_dagger_A(:,p) = new_cols(:,1)
            A_dagger_A(:,q) = new_cols(:,2)
            A_dagger_A(p,p) = pq_block(1,1)
            A_dagger_A(q,p) = pq_block(2,1)
            A_dagger_A(p,q) = pq_block(1,2)
            A_dagger_A(q,q) = pq_block(2,2)

            !! Update eigenvectors
            new_cols(:,1) = c*V(:,p) - conjg(s)*V(:,q)
            new_cols(:,2) = s*V(:,p) + c*V(:,q)
            V(:,p) = new_cols(:,1)
            V(:,q) = new_cols(:,2)

            iter = iter + 1
        enddo

        !! Find unitary matrix U
        ! Find Sigma inverse
        do i = 1, NCOL
            Sigma_inv(i) = 1.0/sqrt(real(A_dagger_A(i,i), kind=real64))
        enddo

        ! Find U = A * V * Sigma_inv * V_dagger
        U = herm(V) ! U = V_dagger
        do i = 1, NCOL
            U(i,:) = Sigma_inv(i)*U(i,:) ! U = Sigma_inv * V_dagger
        enddo
        U = matmul(V, U) ! U = V * Sigma_inv * V_dagger
        U = matmul(A, U)
    end function

    ! Determinant of NxN matrix
    function det(matrix)
        implicit none
        complex(real64), intent(in) :: matrix(:,:)
        complex(real64) :: det
        real(real64) :: pivot_tol

        !! Use the Bareiss algorithm to compute matrix determinant
        complex(real64), allocatable :: M(:,:)
        complex(real64) :: denominator, pivot
        integer :: i, j, k, n

        ! Find size of input matrix
        n = size(matrix, 1)
        if (n /= size(matrix, 2)) stop "Matrix input to det is not square"

        ! Allocate M
        allocate(M(n,n))

        ! Setup copy of matrix
        M = matrix
        denominator = cmplx(1.0,0.0)

        ! Perform itertaion
        do k = 1, n-1
            ! Check if matrix is singular
            pivot = M(k,k)
            pivot_tol = epsilon(1.0d0) * maxval(abs(matrix))
            if (abs(pivot) < pivot_tol) then
                det = cmplx(0.0,0.0)
                exit
            endif

            ! Conduct Bareiss algorithm step
            do i = k+1, n
                do j = k+1, n
                    M(i,j) = (M(i,j)*M(k,k) - M(i,k)*M(k,j))/denominator
                enddo
            enddo
            denominator = pivot
        enddo

        ! Entry n,n now contains the determinant of matrix
        det = M(n,n)
    end function

    ! Normalise a link by projecting it back into SU(N)
    function normalise_link(matrix) result(U)
        implicit none
        complex(real64), intent(in) :: matrix(NCOL, NCOL)
        complex(real64) :: U(NCOL, NCOL), determinant

        ! Unitarise smeared link to sit in U(N)
        U = unitarise_SVD(matrix)

        ! Normalise so that smeared links had determinant = 1, and thus sits in SU(N)
        determinant = det(U)
        determinant = determinant/abs(determinant) ! ensure det(U) is a phase as expected
        ! Scale U to force det(U) = 1 whilst maintaining unitarity
        U = U / (determinant**(1.0/real(NCOL)))
    end function

    ! Get diagonal links (diagonal_links) and lattice pointers to sites diagonal links end on (lattice_pointers_diagonal)
    subroutine get_diagonal_links(gauge_field_slice, direction, blocking_level, diagonal_links, &
        lattice_pointers_diagonal)
        implicit none
        complex(real64), intent(in) :: gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: direction, blocking_level
        complex(real64), intent(out) :: diagonal_links(NCOL, NCOL, SLICE_VOLUME, 4)
        integer, intent(out) :: lattice_pointers_diagonal(SLICE_VOLUME, 4)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! diagonal_links(n, nu_ku) (renamed from UDD) stores the diagonal link in the plane (nu, ku) orthogonal to the direction mu, where:
        ! nu_ku = 1    -->    (nu, ku) = (+, +)
        ! nu_ku = 2    -->    (nu, ku) = (+, -)
        ! nu_ku = 3    -->    (nu, ku) = (-, -)
        ! nu_ku = 4    -->    (nu, ku) = (-, +)
        ! lattice_pointers_diagonal(n, nu_ku) (renamed from IDD) stores the lattice site index that the diagonal link diagonal_links(n, nu_ku) ends on
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        integer :: mu, nu, ku, site, path_1_site, path_2_site, path_site
        complex(real64) :: path_1(NCOL, NCOL), path_2(NCOL, NCOL), path_sum(NCOL, NCOL)

        ! Get orthogonal directions to the input direction, mu
        ! ____mu|____nu|____ku|
        !      1|     2|     3|
        !      2|     3|     1|
        !      3|     1|     2|
        mu = direction
        nu = mod(mu, 3) + 1
        ku = 6 - mu - nu

        do site = 1, SLICE_VOLUME
            !! Compute diagonal links in (+, +) direction

            ! Compute path_1, which starts in the +nu direction
            path_1_site = move(site, nu, blocking_level)
            path_1 = matmul(gauge_field_slice(:,:,site,nu), gauge_field_slice(:,:,path_1_site,ku))

            ! Compute path_2, which starts in the +ku direction
            path_2_site = move(site, ku, blocking_level)
            path_2 = matmul(gauge_field_slice(:,:,site,ku), gauge_field_slice(:,:,path_2_site,nu))

            ! Add sum to diagonal_links
            path_sum = path_1 + path_2
            diagonal_links(:, :, site, 1) = path_sum

            ! Find site these paths point to and save in lattice_pointers_diagonal
            path_site = move(path_1_site, ku, blocking_level)
            lattice_pointers_diagonal(site, 1) = path_site

            !! Compute diagonal links in (-, -) direction from path_site using matrices already calculated. This is the Hermitian conjugate of the (+, +) direction, but labelled from path_site rather than site.

            ! Add Hermitian conjugate of path sum to diagonal links element 3 from path_site
            diagonal_links(:, :, path_site, 3) = herm(path_sum)

            ! Add lattice pointer to path_site, 3 index
            lattice_pointers_diagonal(path_site, 3) = site

            !! Compute diagonal links in (+, -) direction

            ! Compute path 1, which starts in the +nu direction (path_1_site remains unchanged)
            path_site = move(path_1_site, -ku, blocking_level)
            path_1 = matmul(gauge_field_slice(:,:,site,nu), herm(gauge_field_slice(:,:,path_site,ku)))

            ! Compute path 2, which starts in the -ku direction
            path_2_site = move(site, -ku, blocking_level)
            path_2 = matmul(herm(gauge_field_slice(:,:,path_2_site,ku)), gauge_field_slice(:,:,path_2_site,nu))

            ! Add sum to diagonal links
            path_sum = path_1 + path_2
            diagonal_links(:, :, site, 2) = path_sum

            ! Save path_site in lattice_pointers_diagonal
            lattice_pointers_diagonal(site, 2) = path_site

            !! Compute diagonal links in (-, +) direction from path_site using matrices already calculated. This is analagous to how the (-, -) direction is calculated from the results of the (+, +) direction

            ! Add Hermitian conjugate of path sum to diagonal links element 4 from path_site
            diagonal_links(:, :, path_site, 4) = herm(path_sum)

            ! Add lattice pointer to path_site, 4 index
            lattice_pointers_diagonal(path_site, 4) = site
        enddo
    end subroutine

    ! Return smeared input file configuration at given blocking level
    function get_smeared_gauge_field(gauge_field_slice, blocking_level) result(smear)
        implicit none
        complex(real64), intent(in) :: gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: blocking_level
        complex(real64) :: smear(NCOL, NCOL, SLICE_VOLUME, 3)

        complex(real64) :: diagonal_links_mu(NCOL, NCOL, SLICE_VOLUME, 4), staple(NCOL, NCOL), determinant
        integer :: diagonal_pointers_mu(SLICE_VOLUME, 4), mu, nu, nu_ku, site, temp_site1, temp_site2, site_plus_mu

        ! Smear over each direction individually
        do mu = 1, 3 ! mu is the direction of the link to be smeared
            ! Get diagonal links to all links in plane perpendicular to mu
            call get_diagonal_links(gauge_field_slice, mu, blocking_level, diagonal_links_mu, diagonal_pointers_mu)

            ! Smear each link pointing in direction mu individually
            do site = 1, SLICE_VOLUME
                ! Get site the given link points to
                site_plus_mu = move(site, mu, blocking_level)
                !! smeared link = original link + STAPLE_WEIGHT * staples + DIAGONAL_STAPLE_WEIGHT * diagonal staples
                ! Initialise smeared link as original link
                smear(:, :, site, mu) = gauge_field_slice(:, :, site, mu)

                ! Get staples over all directions perpendicular to mu
                do nu = 1, 3
                    ! Ignore direction parallel to mu
                    if (nu == mu) cycle

                    ! Get staple in +nu direction
                    ! To do this, we need (site + nu) -> store in temp_site1
                    temp_site1 = move(site, nu, blocking_level)
                    staple = matmul(matmul(&
                        gauge_field_slice(:,:,site,nu), &
                        gauge_field_slice(:,:,temp_site1,mu)), &
                        herm(gauge_field_slice(:,:,site_plus_mu,nu)))

                    ! Add staple to smeared link with appropriate weight
                    smear(:, :, site, mu) = smear(:, :, site, mu) + STAPLE_WEIGHT * staple

                    ! Get staple in -nu direction
                    ! To do this, we need (site - nu) -> store in temp_site1, and (site - nu + mu) -> store in temp_site2
                    temp_site1 = move(site, -nu, blocking_level)
                    temp_site2 = move(site_plus_mu, -nu, blocking_level)
                    staple = matmul(matmul(&
                        herm(gauge_field_slice(:,:,temp_site1,nu)), &
                        gauge_field_slice(:,:,temp_site1,mu)), &
                        gauge_field_slice(:,:,temp_site2,nu))

                    ! Add staple to smeared link with appropriate weight
                    smear(:, :, site, mu) = smear(:, :, site, mu) + STAPLE_WEIGHT * staple
                enddo

                ! Get diagonal staples over all diagonal directions perpendicular to mu
                do nu_ku = 1, 4
                    ! Get staple starting in diagonal direction nu_ku from site
                    ! To do this, we need (site + nu_ku) -> store in temp_site1
                    temp_site1 = diagonal_pointers_mu(site, nu_ku)

                    ! Get diagonal staple
                    staple = matmul(matmul(&
                        diagonal_links_mu(:,:,site,nu_ku), &
                        gauge_field_slice(:,:,temp_site1,mu)), &
                        herm(diagonal_links_mu(:,:,site_plus_mu,nu_ku)))

                    ! Add staple to smeared link with appropriate weight
                    smear(:, :, site, mu) = smear(:, :, site, mu) + DIAGONAL_STAPLE_WEIGHT * staple
                enddo

                ! Unitarise smeared link to sit in SU(N)
                smear(:, :, site, mu) = normalise_link(smear(:, :, site, mu))
            enddo
        enddo
    end function get_smeared_gauge_field

    function get_blocked_gauge_field(gauge_field_slice) result(blok)
        implicit none
        complex(real64), intent(in) :: gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3)
        complex(real64) :: blok(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        integer :: current_blocking_level, next_blocking_level, mu, site
        complex(real64) :: smeared_gauge_field(NCOL, NCOL, SLICE_VOLUME, 3)

        ! Initialise first blocking level as original lattice
        blok(:, :, :, :, 1) = gauge_field_slice

        ! Construct each blocked lattice from the previous blocking level
        do current_blocking_level = 1, MAX_BLOCKING_LEVEL-1
            ! Get next blocking level
            next_blocking_level = current_blocking_level + 1

            ! Smear configuration at current blocking level
            smeared_gauge_field = get_smeared_gauge_field(blok(:, :, :, :, current_blocking_level), &
            current_blocking_level)

            ! Form blocked configuration at next blocking level by blocking this smeared configuration
            do mu = 1, 3
                do site = 1, SLICE_VOLUME
                    blok(:, :, site, mu, next_blocking_level) &
                    = matmul(blok(:, :, site, mu, current_blocking_level), &
                    blok(:, :, move(site,mu,current_blocking_level), mu, current_blocking_level))
                enddo
            enddo
        enddo
    end function get_blocked_gauge_field
end module lattice