module analysis
    use iso_fortran_env, only : int8, int32, int64, real32, real64

    implicit none

    integer :: NUM_OPERATORS, MAX_DELTA_T, NUM_BINS, NCONFIG

    contains

    ! Read vevs and correlation matrix into memory, and remove duplicate operators and zero operators
    subroutine read_clean_data(vevs_filename, corr_matrix_filename, vevs, corr_matrix)
        implicit none
        character(len=255), intent(in) :: vevs_filename, corr_matrix_filename
        complex(real64), allocatable, intent(out) :: vevs(:,:), corr_matrix(:,:,:,:)

        complex(real64), allocatable :: vevs_from_file(:,:), corr_matrix_from_file(:,:,:,:), &
        diagonal_correlators(:)
        integer(int32) :: dim1, dim2, dim3, dim4, vev_shape(2), corr_matrix_shape(4), i, j, nconfig_check
        integer, allocatable :: distinct_operators(:)
        logical, allocatable :: operator_comparison(:)
        integer :: indices_i(4), indices_j(4)

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

    ! Run GEVP on correlation matrix

    ! Calculate effective masses
end module