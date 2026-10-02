module analysis
    use iso_fortran_env, only : int8, int32, int64, real32, real64

    implicit none

    integer :: NUM_OPERATORS, MAX_DELTA_T, NUM_BINS

    contains

    ! Read vevs and correlation matrix into memory, and remove duplicate operators and zero operators
    subroutine read_clean_data(vevs_filename, corr_matrix_filename, vevs, corr_matrix)
        implicit none
        character(len=255), intent(in) :: vevs_filename, corr_matrix_filename
        complex(real64), allocatable, intent(out) :: vevs(:,:), corr_matrix(:,:,:,:)

        complex(real64), allocatable :: vevs_from_file(:,:), corr_matrix_from_file(:,:,:,:)
        integer(int32) :: dim1, dim2, dim3, dim4, vev_shape(2), corr_matrix_shape(4), i, j
        integer, allocatable :: duplicate_operators(:)
        integer :: indices_i(4), indices_j(4)

        !! Read vevs from file
        ! Open file
        open(11, file=trim(vevs_filename), form='unformatted', access='stream', status='old')

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
        if (size(corr_matrix, dim=1) /= NUM_OPERATORS) &
        error stop "vevs and corr_matrix incompatible: different number of operators"
        if (size(corr_matrix, dim=1) /= size(corr_matrix, dim=2)) &
        error stop "corr_matrix not square"

        ! Assume time seperations contained in corr_matrix(i,j,:,l)
        MAX_DELTA_T = size(corr_matrix, dim=3)

        ! Assume bins contained in vevs(i,:), corr_matrix(i,j,k,:)
        NUM_BINS = size(vevs, dim=2)

        ! Check there are the same number of bins in corr_matrix as in vevs
        if (size(corr_matrix, dim=4) /= NUM_BINS) &
        error stop "vevs and corr_matrix incompatible: different number of bins"

        ! Output dimensions found
        write(*, '(a, i0)') "[Info] Number of operators found: ", NUM_OPERATORS
        write(*, '(a, i0)') "[Info] Number time seperations found: ", MAX_DELTA_T + 1
        write(*, '(a, i0)') "[Info] Number of bins found: ", NUM_BINS


        !! Find duplicate operators by comparing diagonal elements of correlation matrix at zero seperation
        do j = 2, NUM_OPERATORS
            do i = 1, j-1
                
            enddo
        enddo
    end subroutine

    ! Arrange bins into jacknife bins

    ! Remove vevs from correlation matrix

    ! Run GEVP on correlation matrix

    ! Calculate effective masses
end module