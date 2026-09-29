module states_class
    use parameters
    use lattice

    implicit none

    private
    public :: torelon_state
    public :: torelon_momentum_state

    type torelon_state
        private

        ! Descriptor variables
        integer :: SPIN, PP, PR, NUM_OPERATORS
        integer, allocatable :: LINE_INDICES(:) ! Gives the indices of the line that project onto this state

        ! Measurement variables
        complex(real64), allocatable :: vevs(:,:), corr_matrix(:,:,:,:)

        contains

        ! Allocate all parameters
        procedure :: init=>init_torelon_state

        ! Output number of operators
        procedure :: get_num_operators=>get_num_operators_torelon
        
        ! Update vacuum expectation values of operators contributing to this state
        procedure :: update_vevs=>update_torelon_vevs

        ! Update correlation matrix of this state
        procedure :: update_corr_matrix=>update_torelon_corr_matrix

        ! Write results to disk
        procedure :: output_results=>output_torelon_results
    end type

    type torelon_momentum_state
        private

        ! Descriptor variables
        integer :: SPIN, PARITY, MOMENTUM, NUM_OPERATORS
        integer, allocatable :: LINE_INDICES(:) ! Gives the indices of the line that project onto this state

        ! Measurment variables
        complex(real64), allocatable :: vevs(:,:), corr_matrix(:,:,:,:)
        contains

        ! Allocate all parameters
        procedure :: init=>init_torelon_momentum_state

        ! Output number of operators
        procedure :: get_num_operators=>get_num_operators_torelon_momentum

        ! Update vacuum expectation values of operators contributing to this state
        procedure :: update_vevs=>update_torelon_momentum_vevs

        ! Update correlation matrix of this state
        procedure :: update_corr_matrix=>update_torelon_momentum_corr_matrix

        ! Write results to disk
        procedure :: output_results=>output_torelon_momentum_results
    end type

    interface time_slice_average
        module procedure time_slice_average_single
        module procedure time_slice_average_double
        module procedure time_slice_average_double_complex
    end interface

    contains

    ! Compute average of quantity phi defined on each time slice
    pure function time_slice_average_single(phi, t) result(time_slice_avg)
        implicit none

        integer, intent(in) :: t
        real(real64), dimension(LX4), intent(in) :: phi
        real(real64) :: time_slice_avg

        integer :: s, t_plus_s

        time_slice_avg = 0.0
        t_plus_s = t
        do s = 1, LX4
            t_plus_s = mod(t_plus_s, LX4) + 1
            time_slice_avg = time_slice_avg + phi(t_plus_s)*phi(s)
        enddo
    end function

    ! Compute average of quantity phi defined on each time slice
    pure function time_slice_average_double(phi_1, phi_2, t) result(time_slice_avg)
        implicit none

        integer, intent(in) :: t
        real(real64), dimension(LX4), intent(in) :: phi_1, phi_2
        real(real64) :: time_slice_avg

        integer :: s, t_plus_s

        time_slice_avg = 0.0
        t_plus_s = t
        do s = 1, LX4
            t_plus_s = mod(t_plus_s, LX4) + 1
            time_slice_avg = time_slice_avg + phi_1(t_plus_s)*phi_2(s)
        enddo
    end function

    ! Compute average of complex quantity phi defined on each time slice
    pure function time_slice_average_double_complex(phi_1, phi_2, t) result(time_slice_avg)
        implicit none

        integer, intent(in) :: t
        complex(real64), dimension(LX4), intent(in) :: phi_1, phi_2
        complex(real64) :: time_slice_avg

        integer :: s, t_plus_s

        time_slice_avg = cmplx(0.0,0.0)
        t_plus_s = t
        do s = 1, LX4
            t_plus_s = mod(t_plus_s, LX4) + 1
            time_slice_avg = time_slice_avg &
            + 0.5d0 * (phi_1(t))*conjg(phi_2(t_plus_s) + phi_1(t_plus_s)*conjg(phi_2(t)))
        enddo
    end function

    ! Allocate variables for torelon state
    subroutine init_torelon_state(torelon, spin, transverse_parity, parallel_parity, num_operators)
        implicit none
        class(torelon_state), intent(inout) :: torelon
        integer, intent(in) :: spin, transverse_parity, parallel_parity, num_operators

        ! Check variables take appropriate values
        if (spin < 0 .or. spin > 2) then
            error stop "Invalid spin passed to init_torelon_state. Must be 0, 1 or 2."
        endif
        if (abs(transverse_parity) /= 1) then
            error stop "Invalid parity pp passed to init_torelon_state. Must be +-1."
        endif
        if (abs(parallel_parity) /= 1) then
            error stop "Invalid parity pr passed to init_torelon_state. Must be +-1."
        endif
        if (num_operators < 0) then
            error stop "Invalid number of operators passed to init_torelon_state. Must be greater than 0."
        endif

        ! Set state variables based on input
        torelon%SPIN = spin
        torelon%PP = transverse_parity
        torelon%PR = parallel_parity
        torelon%NUM_OPERATORS = num_operators

        ! Check variables are compatible
        if (spin == 1 .and. transverse_parity == -1) then
            write(*, '(a)') "[Warning]    J=1, PP=- state is undefined. Setting transverse parity to +"
            torelon%PP = 1
        endif

        ! Set line indices array
        select case(spin)
        case(0)
            select case(transverse_parity)
            case(1)
                select case(parallel_parity)
                case(1)
                    torelon%LINE_INDICES = [1, 2, 5, 8, 11, 14, 17, 20, 23, 26, 32, 38, 44, 50, 53, 56, 62, &
                    68, 74, 80, 86, 98, 104, 110, 116, 122, 128, 134, 140, 146, 156, 166, 176, 186, 196, 206, &
                    216, 226]
                case(-1)
                    torelon%LINE_INDICES = [87, 147, 157, 167, 177, 187, 197, 207, 217, 227]
                end select
            case(-1)
                select case(parallel_parity)
                case(1)
                    torelon%LINE_INDICES = [63, 81, 88, 123, 129, 135, 148, 158, 168, 178, 188, 198, 208, 218, &
                    228]
                case(-1)
                    torelon%LINE_INDICES = [27, 33, 39, 45, 57, 69, 75, 89, 99, 105, 111, 117, 141, 149, 159, &
                    169, 179, 189, 199, 209, 219, 229]
                end select
            end select
        case(1)
            select case(parallel_parity)
            case(1)
                torelon%LINE_INDICES = [3, 6, 18, 21, 28, 34, 40, 46, 51, 58, 64, 65, 70, 76, 90, 91, 100, &
                106, 112, 118, 124, 125, 130, 131, 142, 150, 160, 170, 180, 190, 200, 210, 220, 230]
            case(-1)
                torelon%LINE_INDICES = [9, 12, 15, 24, 29, 35, 41, 47, 54, 59, 71, 77, 82, 83, 92, 93, 101, &
                107, 113, 136, 137, 143, 151, 161, 171, 181, 191, 201, 211, 221, 231, 119]
            end select
        case(2)
            select case(transverse_parity)
            case(1)
                select case(parallel_parity)
                case(1)
                    torelon%LINE_INDICES = [4, 7, 10, 13, 16, 19, 22, 25, 52, 55, 66, 72, 78, 84, 94, 126, &
                    132, 138, 152, 162, 172, 182, 192, 202, 212, 222, 232]
                case(-1)
                    torelon%LINE_INDICES = [30, 36, 42, 48, 60, 95, 102, 108, 114, 120, 144, 153, 163, 173, &
                    183, 193, 203, 213, 223, 233]
                end select
            case(-1)
                select case(parallel_parity)
                case(1)
                    torelon%LINE_INDICES = [31, 37, 43, 49, 61, 67, 85, 96, 103, 109, 115, 121, 127, 133, 139, &
                    145, 154, 164, 174, 184, 194, 204, 214, 224, 234]
                case(-1)
                    torelon%LINE_INDICES = [97, 155, 165, 175, 185, 195, 205, 215, 225, 235, 73, 79]
                end select
            end select
        end select

        ! Throw warning if number of operators passed not equal to number of operators inferred from LINE_INDICES
        if (num_operators /= size(torelon%LINE_INDICES) * MAX_BLOCKING_LEVEL) then
            write(*, '(a, a, i0)') &
            "[Warning]    Number of operators passed is not equal to that of line indices. ", &
            "Setting number of operators to ", size(torelon%LINE_INDICES) * MAX_BLOCKING_LEVEL
            torelon%NUM_OPERATORS = size(torelon%LINE_INDICES) * MAX_BLOCKING_LEVEL
        endif

        ! Allocate measurement variables
        allocate(torelon%vevs(torelon%NUM_OPERATORS, NUM_BINS), &
        torelon%corr_matrix(torelon%NUM_OPERATORS, torelon%NUM_OPERATORS, 0:MAX_DELTA_T, NUM_BINS))

        ! Set measurement variables to zero
        torelon%vevs = 0.0d0
        torelon%corr_matrix = 0.0d0
    end subroutine

    ! Allocate variables for torelon state
    subroutine init_torelon_momentum_state(torelon, spin, parity, momentum, num_operators)
        implicit none
        class(torelon_momentum_state), intent(inout) :: torelon
        integer, intent(in) :: spin, parity, momentum, num_operators

        ! Check variables take appropriate values
        if (spin < 0 .or. spin > 2) then
            error stop "Invalid spin passed to init_torelon_state. Must be 0, 1 or 2."
        endif
        if (abs(parity) /= 1) then
            error stop "Invalid parity passed to init_torelon_momentum_state, Must be +-1."
        endif
        if (abs(momentum) /= 0) then
            error stop "Invalid momentum passed to init_torelon_momentum_state. Must not be 0."
        endif
        if (num_operators < 0) then
            error stop &
            "Invalid number of operators passed to init_torelon_momentum_state. Must be greater than 0."
        endif

        ! Set state variables based on input
        torelon%SPIN = spin
        torelon%PARITY = parity
        torelon%MOMENTUM = momentum
        torelon%NUM_OPERATORS = num_operators

        ! Check variables are compatible
        if (spin == 1 .and. parity == -1) then
            write(*, '(a)') "[Warning]    J=1, P=- non-zero momentum state is undefined. Setting P to +"
            torelon%PARITY = 1
        endif

        ! Set line indices array
        select case(spin)
        case(0)
            select case(parity)
            case(1)
                torelon%LINE_INDICES = [2, 5, 8, 11, 14, 17, 20, 23, 26, 32, 38, 44, 50, 53, 56, 62, 68, 74, &
                80, 86, 87, 98, 227, 104, 110, 116, 122, 128, 134, 140, 146, 147, 156, 157, 166, 167, 176, &
                177, 186, 187, 196, 197, 206, 207, 216, 217, 226]
            case(-1)
                torelon%LINE_INDICES = [27, 33, 39, 45, 57, 63, 69, 75, 81, 88, 89, 105, 111, 117, 123, 129, &
                135, 141, 148, 149, 158, 159, 168, 169, 178, 179, 188, 189, 198, 199, 208, 209, 218, 219, 228, &
                229, 99]
            end select
        case(1)
            torelon%LINE_INDICES = [3, 6, 9, 12, 15, 18, 21, 24, 28, 29, 34, 35, 40, 41, 46, 47, 51, 54, 58, &
            59, 64, 65, 70, 71, 76, 77, 82, 83, 90, 91, 92, 93, 100, 101, 106, 107, 112, 113, 118, 119, 124, &
            125, 130, 131, 136, 137, 142, 143, 150, 151, 160, 161, 170, 171, 180, 181, 190, 191, 200, 201, &
            210, 211, 220, 221, 230, 231]
        case(2)
            select case(parity)
            case(1)
                torelon%LINE_INDICES = [4, 7, 10, 13, 16, 19, 22, 25, 30, 36, 42, 48, 52, 55, 60, 66, 72, 78, &
                84, 94, 95, 102, 108, 114, 120, 126, 132, 138, 144, 152, 153, 162, 163, 172, 173, 182, 183, &
                192, 193, 202, 203, 212, 213, 222, 223, 232, 233]
            case(-1)
                torelon%LINE_INDICES = [31, 37, 43, 49, 61, 67, 73, 79, 85, 96, 97, 103, 109, 115, 121, 127, &
                133, 139, 145, 154, 155, 164, 165, 174, 175, 184, 185, 194, 195, 204, 205, 214, 215, 224, 225, &
                234, 235]
            end select
        end select

        ! Throw warning if number of operators passed not equal to number of operators inferred from LINE_INDICES
        if (num_operators /= size(torelon%LINE_INDICES) * MAX_BLOCKING_LEVEL) then
            write(*, '(a, a, i0)') &
            "[Warning]    Number of operators passed is not equal to that of line indices. ", &
            "Setting number of operators to ", size(torelon%LINE_INDICES) * MAX_BLOCKING_LEVEL
            torelon%NUM_OPERATORS = size(torelon%LINE_INDICES) * MAX_BLOCKING_LEVEL
        endif

        ! Allocate measurement variables
        allocate(torelon%vevs(torelon%NUM_OPERATORS, NUM_BINS), &
        torelon%corr_matrix(torelon%NUM_OPERATORS, torelon%NUM_OPERATORS, 0:MAX_DELTA_T, NUM_BINS))

        ! Set measurement variables to zero
        torelon%vevs = 0.0d0
        torelon%corr_matrix = 0.0d0
    end subroutine

    ! Update vevs of torelon operators that contribute to the state from lines calculated in THERML1
    subroutine update_torelon_vevs(torelon, lines, bin_index)
        implicit none
        class(torelon_state), intent(inout) :: torelon
        complex(real64), intent(in) :: lines(LX4, MAX_BLOCKING_LEVEL, 235)
        integer, intent(in) :: bin_index

        integer :: id, line_index, blocking_level, i

        do concurrent (i = 1:size(torelon%LINE_INDICES), blocking_level = 1:MAX_BLOCKING_LEVEL)
            ! Find index of operator being accessed
            line_index = torelon%LINE_INDICES(i)
            id = (line_index - 1) * MAX_BLOCKING_LEVEL + blocking_level
            
            ! Update vevs
            torelon%vevs(id, bin_index) = torelon%vevs(id, bin_index) &
            + sum(lines(:, blocking_level, line_index))
        enddo
    end subroutine

    ! Output number of operators used in torelon state
    function get_num_operators_torelon(torelon) result(num_operators)
        implicit none
        class(torelon_state), intent(inout) :: torelon
        integer :: num_operators

        num_operators = torelon%NUM_OPERATORS
    end function

    ! Output number of operators used in torelon state with momentum
    function get_num_operators_torelon_momentum(torelon) result(num_operators)
        implicit none
        class(torelon_momentum_state), intent(inout) :: torelon
        integer :: num_operators

        num_operators = torelon%NUM_OPERATORS
    end function

    ! Update vevs of torelon operators with momentum that contribute to the state from lines calculated in THERML1
    subroutine update_torelon_momentum_vevs(torelon, momentum_lines, bin_index)
        implicit none
        class(torelon_momentum_state), intent(inout) :: torelon
        complex(real64), intent(in) :: momentum_lines(LX4, MAX_BLOCKING_LEVEL, 2, 2:235)
        integer, intent(in) :: bin_index

        integer :: id, line_index, blocking_level, i

        do concurrent (i = 1:size(torelon%LINE_INDICES), blocking_level = 1:MAX_BLOCKING_LEVEL)
            ! Find index of operator being accessed
            line_index = torelon%LINE_INDICES(i)
            id = (line_index - 1) * MAX_BLOCKING_LEVEL + blocking_level
            
            ! Update vevs
            torelon%vevs(id, bin_index) = torelon%vevs(id, bin_index) &
            + sum(momentum_lines(:, blocking_level, torelon%MOMENTUM, line_index))
        enddo
    end subroutine

    ! Update correlation matrix of torelon operators that contribute to the state from lines calculated in THERML1
    subroutine update_torelon_corr_matrix(torelon, lines, bin_index)
        implicit none
        class(torelon_state), intent(inout) :: torelon
        complex(real64), intent(in) :: lines(LX4, MAX_BLOCKING_LEVEL, 235)
        integer, intent(in) :: bin_index

        integer :: i, j, id1, id2, line_index1, line_index2, b1, b2, delta_t

        do concurrent &
            (delta_t = 0:MAX_DELTA_T, i = 1:size(torelon%LINE_INDICES), j = 1:size(torelon%LINE_INDICES), &
            b1 = 1:MAX_BLOCKING_LEVEL, b2=1:MAX_BLOCKING_LEVEL)
            ! Find id's of loops being accessed
            line_index1 = torelon%LINE_INDICES(i)
            line_index2 = torelon%LINE_INDICES(j)
            id1 = (line_index1 - 1) * MAX_BLOCKING_LEVEL + b1
            id2 = (line_index2 - 1) * MAX_BLOCKING_LEVEL + b2

            ! Update correlation matrices
            torelon%corr_matrix(id1, id2, delta_t, bin_index) &
            = torelon%corr_matrix(id1, id2, delta_t, bin_index) &
            + time_slice_average(lines(:, b1, line_index1), lines(:, b2, line_index2), delta_t)
        enddo
    end subroutine

    ! Update correlation matrix of torelon operators with momentum that contribute to the state from lines calculated in THERML1
    subroutine update_torelon_momentum_corr_matrix(torelon, momentum_lines, bin_index)
        implicit none
        class(torelon_momentum_state), intent(inout) :: torelon
        complex(real64), intent(in) :: momentum_lines(LX4, MAX_BLOCKING_LEVEL, 2, 2:235)
        integer, intent(in) :: bin_index

        integer :: i, j, id1, id2, line_index1, line_index2, b1, b2, delta_t

        do concurrent &
            (delta_t = 0:MAX_DELTA_T, i = 1:size(torelon%LINE_INDICES), j = 1:size(torelon%LINE_INDICES), &
            b1 = 1:MAX_BLOCKING_LEVEL, b2 = 1:MAX_BLOCKING_LEVEL)
            ! Find id's of loops being accessed
            line_index1 = torelon%LINE_INDICES(i)
            line_index2 = torelon%LINE_INDICES(j)
            id1 = (line_index1 - 1) * MAX_BLOCKING_LEVEL + b1
            id2 = (line_index2 - 1) * MAX_BLOCKING_LEVEL + b2

            ! Update correlation matrices
            torelon%corr_matrix(id1, id2, delta_t, bin_index) &
            = torelon%corr_matrix(id1, id2, delta_t, bin_index) &
            + time_slice_average(momentum_lines(:, b1, torelon%MOMENTUM, line_index1), &
            momentum_lines(:, b2, torelon%MOMENTUM, line_index2), delta_t)
        enddo
    end subroutine

    ! Write vevs and correlation matrix of torelon state to disk
    subroutine output_torelon_results(torelon)
        implicit none
        class(torelon_state), intent(inout) :: torelon

        character(len=15) :: file_J, file_PP, file_PR
        character(len=255) :: file_name

        ! Write spin and parities into strings
        write(file_J, '(i0)') torelon%SPIN
        if (torelon%SPIN /= 1) then ! J=1 states have undefined transverse parity
            select case(torelon%PP)
            case(1)
                file_PP = 'P'
            case(-1)
                file_PP = 'M'
            end select
        endif
        select case(torelon%PR)
        case(1)
            file_PR = 'P'
        case(-1)
            file_PR = 'M'
        end select

        ! Construct vevs file name
        if (torelon%SPIN == 1) then
            file_name = 'vevs_J' // trim(file_J) // trim(file_PR) // 'q0.dat'
        else
            file_name = 'vevs_J' // trim(file_J) // trim(file_PP) // trim(file_PR) // 'q0.dat'
        endif

        ! Write dimension of array and vevs to file
        open(11, file=trim(file_name), form='unformatted', access='stream', status='replace')
        write(11) size(torelon%vevs, dim=1, kind=int32)
        write(11) size(torelon%vevs, dim=2, kind=int32)
        write(11) torelon%vevs
        close(11)

        ! Construct correlation matrix file name
        if (torelon%SPIN == 1) then
            file_name = 'corr_matrix_J' // trim(file_J) // trim(file_PR) // 'q0.dat'
        else
            file_name = 'corr_matrix_J' // trim(file_J) // trim(file_PP) // trim(file_PR) // 'q0.dat'
        endif

        ! Write dimension of array and correlation matrix to file
        open(11, file=trim(file_name), form='unformatted', access='stream', status='replace')
        write(11) size(torelon%corr_matrix, dim=1, kind=int32)
        write(11) size(torelon%corr_matrix, dim=2, kind=int32)
        write(11) size(torelon%corr_matrix, dim=3, kind=int32)
        write(11) size(torelon%corr_matrix, dim=4, kind=int32)
        write(11) torelon%corr_matrix
        close(11)
    end subroutine

    ! Write vevs and correlation matrix of torelon state with momentum to disk
    subroutine output_torelon_momentum_results(torelon)
        implicit none
        class(torelon_momentum_state), intent(inout) :: torelon

        character(len=15) :: file_J, file_P, file_q
        character(len=255) :: file_name

        ! Write spin and parities into strings
        write(file_J, '(i0)') torelon%SPIN
        if (torelon%SPIN /= 1) then ! J=1 states have undefined parity
            select case(torelon%PARITY)
            case(1)
                file_P = 'P'
            case(-1)
                file_P = 'M'
            end select
        endif
        write(file_q, '(i0)') torelon%MOMENTUM

        ! Construct vevs file name
        if (torelon%SPIN == 1) then
            file_name = 'vevs_J' // trim(file_J) // 'q' // trim(file_q) // '.dat'
        else
            file_name = 'vevs_J' // trim(file_J) // trim(file_P) // 'q' // trim(file_q) // '.dat'
        endif

        ! Write dimension of array and vevs to file
        open(11, file=trim(file_name), form='unformatted', access='stream', status='replace')
        write(11) size(torelon%vevs, dim=1, kind=int32)
        write(11) size(torelon%vevs, dim=2, kind=int32)
        write(11) torelon%vevs
        close(11)

        ! Construct correlation matrix file name
        if (torelon%SPIN == 1) then
            file_name = 'corr_matrix_J' // trim(file_J) // 'q' // trim(file_q) // '.dat'
        else
            file_name = 'corr_matrix_J' // trim(file_J) // trim(file_P) // 'q' // trim(file_q) // '.dat'
        endif

        ! Write dimension of array and correlation matrix to file
        open(11, file=trim(file_name), form='unformatted', access='stream', status='replace')
        write(11) size(torelon%corr_matrix, dim=1, kind=int32)
        write(11) size(torelon%corr_matrix, dim=2, kind=int32)
        write(11) size(torelon%corr_matrix, dim=3, kind=int32)
        write(11) size(torelon%corr_matrix, dim=4, kind=int32)
        write(11) torelon%corr_matrix
        close(11)
    end subroutine
end module