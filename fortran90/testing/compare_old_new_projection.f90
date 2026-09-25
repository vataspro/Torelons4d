program compare_old_new_projection
  use iso_fortran_env, only: real32, real64
  implicit none

  integer, parameter :: NCOL = 2
  complex(real32) :: a_old(NCOL,NCOL), u_old(NCOL,NCOL)
  complex(real64) :: a_new(NCOL,NCOL), u_new(NCOL,NCOL), u_old_as64(NCOL,NCOL)
  real(real64) :: diff

  call fixed_random_matrix(a_old)
  a_new = cmplx(real(a_old, kind=real64), real(aimag(a_old), kind=real64), kind=real64)

  call old_projection(a_old, u_old)
  u_new = new_projection(a_new)

  u_old_as64 = cmplx(real(u_old, kind=real64), real(aimag(u_old), kind=real64), kind=real64)
  diff = maxval(abs(u_old_as64 - u_new))

  print *, "Input matrix A:"
  call print_matrix32(a_old)

  print *, "Old projection: RENORMBS + DODET, real32:"
  call print_matrix32(u_old)

  print *, "New projection: SVD/polar + determinant phase, real64:"
  call print_matrix64(u_new)

  print *, "Difference matrix: U_old - U_new:"
  call print_matrix64(u_old_as64 - u_new)

  print *, "max |U_old - U_new| = ", diff
  print *, ""

  print *, "Old projection checks:"
  print *, "max |U^dagger U - I| = ", unitarity_error64(u_old_as64)
  print *, "|det(U) - 1|       = ", abs(det64(u_old_as64) - cmplx(1.0_real64, 0.0_real64, kind=real64))
  print *, ""

  print *, "New projection checks:"
  print *, "max |U^dagger U - I| = ", unitarity_error64(u_new)
  print *, "|det(U) - 1|       = ", abs(det64(u_new) - cmplx(1.0_real64, 0.0_real64, kind=real64))

contains

  subroutine fixed_random_matrix(a)
    complex(real32), intent(out) :: a(NCOL,NCOL)
    integer :: i, j, n
    integer, allocatable :: seed(:)
    real(real32) :: r1, r2

    call random_seed(size=n)
    allocate(seed(n))
    seed = 12345
    call random_seed(put=seed)

    do j = 1, NCOL
      do i = 1, NCOL
        call random_number(r1)
        call random_number(r2)
        a(i,j) = cmplx(2.0_real32*r1 - 1.0_real32, &
                       2.0_real32*r2 - 1.0_real32, kind=real32)
      end do
    end do
  end subroutine fixed_random_matrix

  subroutine old_projection(a, u)
    complex(real32), intent(in) :: a(NCOL,NCOL)
    complex(real32), intent(out) :: u(NCOL,NCOL)
    complex(real32) :: determinant, cdet, cnorm, dncol
    integer :: i, jmat

    call renormbs_old(a, u)

    jmat = NCOL
    call detnant_old(jmat, determinant, u)

    dncol = cmplx(1.0_real32 / real(NCOL, kind=real32), 0.0_real32, kind=real32)
    cdet = determinant**dncol
    cnorm = 1.0_real32 / cdet

    do i = 1, NCOL
      u(:,i) = u(:,i) * cnorm
    end do
  end subroutine old_projection

  subroutine renormbs_old(uu1, uren11)
    complex(real32), intent(in) :: uu1(NCOL,NCOL)
    complex(real32), intent(out) :: uren11(NCOL,NCOL)
    complex(real32) :: csum, adum(NCOL,NCOL)
    complex(real32) :: a11(NCOL,NCOL), b11(NCOL,NCOL), c11(NCOL,NCOL)
    real(real32) :: sum_value, anorm
    integer :: n1, n2, n3, i, j

    adum = uu1

    do n2 = 1, NCOL
      do n3 = 1, n2 - 1
        csum = cmplx(0.0_real32, 0.0_real32, kind=real32)
        do n1 = 1, NCOL
          csum = csum + adum(n1,n2) * conjg(adum(n1,n3))
        end do
        do n1 = 1, NCOL
          adum(n1,n2) = adum(n1,n2) - csum * adum(n1,n3)
        end do
      end do

      sum_value = 0.0_real32
      do n1 = 1, NCOL
        sum_value = sum_value + real(adum(n1,n2) * conjg(adum(n1,n2)), kind=real32)
      end do
      anorm = 1.0_real32 / sqrt(sum_value)
      do n1 = 1, NCOL
        adum(n1,n2) = adum(n1,n2) * anorm
      end do
    end do

    c11 = adum

    do j = 1, NCOL
      do i = 1, NCOL
        a11(i,j) = conjg(uu1(j,i))
      end do
    end do

    call vmx_old(a11, c11, b11)
    call rungrb_old(b11, c11)

    uren11 = c11
  end subroutine renormbs_old

  subroutine rungrb_old(b11, c11)
    complex(real32), intent(inout) :: b11(NCOL,NCOL), c11(NCOL,NCOL)

    call subgrb_old(1, 2, b11, c11)
  end subroutine rungrb_old

  subroutine subgrb_old(ldu, ldl, b11, c11)
    integer, intent(in) :: ldu, ldl
    complex(real32), intent(inout) :: b11(NCOL,NCOL), c11(NCOL,NCOL)
    complex(real32) :: f11, f12, a11(NCOL,NCOL), s11(NCOL,NCOL), t11(NCOL,NCOL)
    real(real32) :: umag
    integer :: ij1

    a11 = cmplx(0.0_real32, 0.0_real32, kind=real32)

    f11 = (b11(ldu,ldu) + conjg(b11(ldl,ldl))) * 0.5_real32
    f12 = (b11(ldl,ldu) - conjg(b11(ldu,ldl))) * 0.5_real32
    umag = sqrt(real(f11*conjg(f11) + f12*conjg(f12), kind=real32))
    umag = 1.0_real32 / umag
    f11 = f11 * umag
    f12 = f12 * umag

    a11(ldu,ldu) = conjg(f11)
    a11(ldu,ldl) = conjg(f12)
    a11(ldl,ldu) = -f12
    a11(ldl,ldl) = f11

    do ij1 = 1, NCOL
      s11(ij1,ldu) = c11(ij1,ldu)*(a11(ldu,ldu)-1.0_real32) + c11(ij1,ldl)*a11(ldl,ldu)
      s11(ij1,ldl) = c11(ij1,ldu)*a11(ldu,ldl) + c11(ij1,ldl)*(a11(ldl,ldl)-1.0_real32)
      t11(ij1,ldu) = b11(ij1,ldu)*(a11(ldu,ldu)-1.0_real32) + b11(ij1,ldl)*a11(ldl,ldu)
      t11(ij1,ldl) = b11(ij1,ldu)*a11(ldu,ldl) + b11(ij1,ldl)*(a11(ldl,ldl)-1.0_real32)
    end do

    do ij1 = 1, NCOL
      c11(ij1,ldu) = c11(ij1,ldu) + s11(ij1,ldu)
      c11(ij1,ldl) = c11(ij1,ldl) + s11(ij1,ldl)
      b11(ij1,ldu) = b11(ij1,ldu) + t11(ij1,ldu)
      b11(ij1,ldl) = b11(ij1,ldl) + t11(ij1,ldl)
    end do
  end subroutine subgrb_old

  subroutine vmx_old(a, b, c)
    complex(real32), intent(in) :: a(NCOL,NCOL), b(NCOL,NCOL)
    complex(real32), intent(out) :: c(NCOL,NCOL)
    complex(real32) :: csum
    integer :: i, j, k

    do j = 1, NCOL
      do i = 1, NCOL
        csum = cmplx(0.0_real32, 0.0_real32, kind=real32)
        do k = 1, NCOL
          csum = csum + a(i,k) * b(k,j)
        end do
        c(i,j) = csum
      end do
    end do
  end subroutine vmx_old

  subroutine detnant_old(numop, determinant, aa)
    integer, intent(in) :: numop
    complex(real32), intent(in) :: aa(numop,numop)
    complex(real32), intent(out) :: determinant
    complex(real32) :: b(numop,numop)

    b = aa
    call detmat_old(b, numop, determinant)
  end subroutine detnant_old

  subroutine detmat_old(a, np, determinant)
    integer, intent(in) :: np
    complex(real32), intent(inout) :: a(np,np)
    complex(real32), intent(out) :: determinant
    integer :: indx(np), j
    real(real32) :: d
    complex(real32) :: dd

    call ludcmp_old(a, np, np, indx, d)
    dd = cmplx(d, 0.0_real32, kind=real32)
    do j = 1, np
      dd = dd * a(j,j)
    end do
    determinant = dd
  end subroutine detmat_old

  subroutine ludcmp_old(a, n, np, indx, d)
    integer, intent(in) :: n, np
    complex(real32), intent(inout) :: a(np,np)
    integer, intent(out) :: indx(n)
    real(real32), intent(out) :: d
    real(real32), parameter :: tiny = 1.0e-20_real32
    real(real32) :: vv(NCOL), aamax, dum
    complex(real32) :: sum_value, cdum
    integer :: i, imax, j, k

    d = 1.0_real32
    do i = 1, n
      aamax = 0.0_real32
      do j = 1, n
        if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))
      end do
      if (aamax == 0.0_real32) return
      vv(i) = 1.0_real32 / aamax
    end do

    do j = 1, n
      do i = 1, j - 1
        sum_value = a(i,j)
        do k = 1, i - 1
          sum_value = sum_value - a(i,k) * a(k,j)
        end do
        a(i,j) = sum_value
      end do

      aamax = 0.0_real32
      imax = j
      do i = j, n
        sum_value = a(i,j)
        do k = 1, j - 1
          sum_value = sum_value - a(i,k) * a(k,j)
        end do
        a(i,j) = sum_value
        dum = vv(i) * abs(sum_value)
        if (dum >= aamax) then
          imax = i
          aamax = dum
        end if
      end do

      if (j /= imax) then
        do k = 1, n
          cdum = a(imax,k)
          a(imax,k) = a(j,k)
          a(j,k) = cdum
        end do
        d = -d
        vv(imax) = vv(j)
      end if

      indx(j) = imax
      if (abs(a(j,j)) == 0.0_real32) a(j,j) = cmplx(tiny, 0.0_real32, kind=real32)
      if (j /= n) then
        cdum = 1.0_real32 / a(j,j)
        do i = j + 1, n
          a(i,j) = a(i,j) * cdum
        end do
      end if
    end do
  end subroutine ludcmp_old

  function new_projection(a) result(u)
    complex(real64), intent(in) :: a(NCOL,NCOL)
    complex(real64) :: u(NCOL,NCOL)
    complex(real64) :: determinant

    u = unitarise_svd_new(a)
    determinant = det64(u)
    determinant = determinant / abs(determinant)
    u = u / (determinant**(1.0_real64 / real(NCOL, kind=real64)))
  end function new_projection

  function unitarise_svd_new(a) result(u)
    complex(real64), intent(in) :: a(NCOL,NCOL)
    complex(real64) :: u(NCOL,NCOL)
    complex(real64) :: a_dagger_a(NCOL,NCOL), v(NCOL,NCOL)
    real(real64) :: sigma_inv(NCOL)
    complex(real64) :: new_rows(2,NCOL), new_cols(NCOL,2), pq_block(2,2)
    real(real64) :: delta, c, t, tau
    complex(real64) :: s
    integer :: index_pair(2)
    logical :: tri_mask(NCOL,NCOL)
    integer :: i, j, p, q
    real(real64), parameter :: tol_svd = 1.0e-8_real64

    a_dagger_a = matmul(herm64(a), a)

    delta = 0.0_real64
    do j = 1, NCOL
      do i = 1, NCOL
        delta = delta + abs(a_dagger_a(i,j))**2
      end do
    end do
    delta = sqrt(delta)
    delta = tol_svd * delta

    v = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do i = 1, NCOL
      v(i,i) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    end do

    do j = 1, NCOL
      do i = 1, NCOL
        tri_mask(i,j) = (i < j)
      end do
    end do

    do while (off64(a_dagger_a) > delta)
      index_pair = maxloc(abs(a_dagger_a), mask=tri_mask)
      p = index_pair(1)
      q = index_pair(2)

      if (abs(a_dagger_a(p,q)) < max(delta, epsilon(1.0_real64)*maxval(abs(a_dagger_a)))) then
        cycle
      else
        tau = (real(a_dagger_a(q,q), kind=real64) - real(a_dagger_a(p,p), kind=real64)) / &
              (2.0_real64 * abs(a_dagger_a(p,q)))
        if (tau >= 0.0_real64) then
          t = 1.0_real64 / (tau + sqrt(1.0_real64 + tau**2))
        else
          t = 1.0_real64 / (tau - sqrt(1.0_real64 + tau**2))
        end if

        c = 1.0_real64 / sqrt(1.0_real64 + t**2)
        s = t * c * a_dagger_a(p,q) / abs(a_dagger_a(p,q))
      end if

      new_rows(1,:) = c*a_dagger_a(p,:) - s*a_dagger_a(q,:)
      new_rows(2,:) = conjg(s)*a_dagger_a(p,:) + c*a_dagger_a(q,:)

      new_cols(:,1) = c*a_dagger_a(:,p) - conjg(s)*a_dagger_a(:,q)
      new_cols(:,2) = s*a_dagger_a(:,p) + c*a_dagger_a(:,q)

      pq_block(1,1) = c*new_cols(p,1) - s*new_cols(q,1)
      pq_block(2,1) = c*new_cols(q,1) + conjg(s)*new_cols(p,1)
      pq_block(1,2) = c*new_cols(p,2) - s*new_cols(q,2)
      pq_block(2,2) = c*new_cols(q,2) + conjg(s)*new_cols(p,2)

      a_dagger_a(p,:) = new_rows(1,:)
      a_dagger_a(q,:) = new_rows(2,:)
      a_dagger_a(:,p) = new_cols(:,1)
      a_dagger_a(:,q) = new_cols(:,2)
      a_dagger_a(p,p) = pq_block(1,1)
      a_dagger_a(q,p) = pq_block(2,1)
      a_dagger_a(p,q) = pq_block(1,2)
      a_dagger_a(q,q) = pq_block(2,2)

      new_cols(:,1) = c*v(:,p) - conjg(s)*v(:,q)
      new_cols(:,2) = s*v(:,p) + c*v(:,q)
      v(:,p) = new_cols(:,1)
      v(:,q) = new_cols(:,2)
    end do

    do i = 1, NCOL
      sigma_inv(i) = 1.0_real64 / sqrt(real(a_dagger_a(i,i), kind=real64))
    end do

    u = herm64(v)
    do i = 1, NCOL
      u(i,:) = sigma_inv(i) * u(i,:)
    end do
    u = matmul(v, u)
    u = matmul(a, u)
  end function unitarise_svd_new

  function herm64(a) result(h)
    complex(real64), intent(in) :: a(:,:)
    complex(real64) :: h(size(a,1),size(a,2))

    h = conjg(transpose(a))
  end function herm64

  function off64(a) result(off)
    complex(real64), intent(in) :: a(NCOL,NCOL)
    real(real64) :: off
    integer :: i, j

    off = 0.0_real64
    do j = 1, NCOL
      do i = 1, j - 1
        off = off + abs(a(i,j))**2
      end do
      do i = j + 1, NCOL
        off = off + abs(a(i,j))**2
      end do
    end do
    off = sqrt(off)
  end function off64

  function det64(matrix) result(determinant)
    complex(real64), intent(in) :: matrix(:,:)
    complex(real64) :: determinant
    complex(real64), allocatable :: m(:,:)
    complex(real64) :: denominator, pivot
    real(real64) :: pivot_tol
    integer :: i, j, k, n

    n = size(matrix, 1)
    allocate(m(n,n))
    m = matrix
    denominator = cmplx(1.0_real64, 0.0_real64, kind=real64)

    do k = 1, n - 1
      pivot = m(k,k)
      pivot_tol = epsilon(1.0_real64) * maxval(abs(matrix))
      if (abs(pivot) < pivot_tol) then
        determinant = cmplx(0.0_real64, 0.0_real64, kind=real64)
        return
      end if

      do i = k + 1, n
        do j = k + 1, n
          m(i,j) = (m(i,j)*m(k,k) - m(i,k)*m(k,j)) / denominator
        end do
      end do
      denominator = pivot
    end do

    determinant = m(n,n)
  end function det64

  function unitarity_error64(u) result(error)
    complex(real64), intent(in) :: u(NCOL,NCOL)
    real(real64) :: error
    complex(real64) :: identity(NCOL,NCOL)

    identity = cmplx(0.0_real64, 0.0_real64, kind=real64)
    identity(1,1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    identity(2,2) = cmplx(1.0_real64, 0.0_real64, kind=real64)

    error = maxval(abs(matmul(herm64(u), u) - identity))
  end function unitarity_error64

  subroutine print_matrix32(m)
    complex(real32), intent(in) :: m(NCOL,NCOL)
    integer :: i

    do i = 1, NCOL
      write(*,'(2("(",ES14.6,",",ES14.6,")",2X))') m(i,1), m(i,2)
    end do
  end subroutine print_matrix32

  subroutine print_matrix64(m)
    complex(real64), intent(in) :: m(NCOL,NCOL)
    integer :: i

    do i = 1, NCOL
      write(*,'(2("(",ES14.6,",",ES14.6,")",2X))') m(i,1), m(i,2)
    end do
  end subroutine print_matrix64

end program compare_old_new_projection
