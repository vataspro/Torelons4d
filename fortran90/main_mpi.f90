program main_mpi
    use mpi_f08
    use torelon_parameters
    use lattice
    use read_field_config
    use here_be_dragons
    use states_class

    implicit none

    ! ---------------------------------------- Initialise variables ---------------------------------------- !

    ! Parallelisation variables
    integer :: my_rank, num_ranks, mpierr, my_LX4, my_t_start, my_t_stop, local_LX4, num_lines, num_mom_lines
    integer, allocatable :: rank_lookup_table(:)
    type(mpi_status) :: send_status(2)
    type(mpi_status), allocatable :: recv_status(:)
    type(mpi_request) :: send_request(2)
    type(mpi_request), allocatable :: recv_request(:)

    ! File access variables
    logical :: file_exists
    character(len=16) :: file_config_id
    character(len=256) :: file_path

    ! Gauge field
    complex(real64), allocatable :: gauge_field(:,:,:,:), gauge_field_blocked(:,:,:,:,:)

    ! Measurement variables
    complex(real64), allocatable, asynchronous :: lines_slice(:,:), lines(:,:,:), &
    momentum_lines_slice(:,:,:), momentum_lines(:,:,:,:)
    integer :: bin_index

    ! State variables
    type(torelon_state) :: states(10)
    type(torelon_momentum_state) :: momentum_states(10)

    ! Loop variables
    integer :: config, t, blocking_level, state, rank


    ! ----------------------------------------- Set up calculation ----------------------------------------- !

    ! Initialise MPI
    call mpi_init(mpierr)

    ! Get this rank number and number of ranks
    call mpi_comm_rank(mpi_comm_world, my_rank, mpierr)
    call mpi_comm_size(mpi_comm_world, num_ranks, mpierr)

    ! Check there is more than one rank
    if (num_ranks < 2) then
        write(*, '(a)') "[Error]    This version requires more than one MPI rank. Aborting..."
        call mpi_abort(mpi_comm_world, 1, mpierr)
    endif

    ! Load all parameters from file
    call initialise_parameters("parameter_file_local.txt")

    ! On rank 0, check configurations to be read are present and calculate what ranks will be sending rank 0 information for a given time slice. On other ranks, assign time slices.
    if (my_rank == 0) then
        ! Check all configurations can be found
        do config = CONFIG_START, CONFIG_STOP, CONFIG_STEP
            write(file_config_id, "(i0)") config
            file_path = trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)
            inquire(file=trim(file_path), exist=file_exists)
            if (.not.file_exists) then
                write(*, '(a, a)') "[Error][File access]    Cannot access file ", file_path
                call mpi_abort(mpi_comm_world, 1, mpierr)
            endif
        enddo

        ! Allocate lookup table to map between time slices and ranks
        allocate(rank_lookup_table(LX4))
        rank_lookup_table = 0

        ! Calculate temporal extent of every rank and add into lookup table
        t = 0
        do rank = 1, num_ranks-1
            ! Calculate number of time slices calculated by this rank
            local_LX4 = LX4 / (num_ranks-1)
            if (rank <= mod(LX4, num_ranks-1)) then
                local_LX4 = local_LX4 + 1
            endif

            ! Fill this many values of the lookup table
            rank_lookup_table(t+1 : t+local_LX4) = rank

            ! Advance time slice
            t = t + local_LX4
        enddo
    else
        ! Assign time slices to ranks
        ! Rank 0 is responsible for constructing the correlation matrices, and therefore does not get assigned any time slices.
        ! Other ranks get floor(LX4 / (num_ranks-1)) by default. Any extra time slices are given to the first ranks until none are left.
        my_LX4 = LX4 / (num_ranks-1)
        my_t_start = (my_rank - 1) * (LX4 / (num_ranks - 1)) + 1
        if (my_rank <= mod(LX4, num_ranks-1)) then
            my_LX4 = my_LX4 + 1
            my_t_start = my_t_start + my_rank - 1
        else
            my_t_start = my_t_start + mod(LX4, num_ranks-1)
        endif
        my_t_stop = my_t_start + my_LX4 - 1
    endif

    ! Setup lattice movers
    call setup_lattice()

    if (my_rank == 0) then
        ! Allocate arrays
        allocate(lines(MAX_BLOCKING_LEVEL, 235, LX4), &
        momentum_lines(MAX_BLOCKING_LEVEL, 2, 2:235, LX4), &
        recv_status(2*LX4), recv_request(2*LX4))

        ! Initialise states
        call states(1)%init(0, 1, 1, 38*MAX_BLOCKING_LEVEL) ! J=0, PP=+, PR=+
        call states(2)%init(0, 1, -1, 10*MAX_BLOCKING_LEVEL) ! J=0, PP=+, PR=-
        call states(3)%init(0, -1, 1, 15*MAX_BLOCKING_LEVEL) ! J=0, PP=-, PR=+
        call states(4)%init(0, -1, -1, 22*MAX_BLOCKING_LEVEL) ! J=0, PP=-, PR=-
        call states(5)%init(1, 1, 1, 34*MAX_BLOCKING_LEVEL) ! J=1, PR=+
        call states(6)%init(1, 1, -1, 32*MAX_BLOCKING_LEVEL) ! J=1, PR=-
        call states(7)%init(2, 1, 1, 27*MAX_BLOCKING_LEVEL) ! J=2, PP=+, PR=+
        call states(8)%init(2, 1, -1, 20*MAX_BLOCKING_LEVEL) ! J=2, PP=+, PR=-
        call states(9)%init(2, -1, 1, 25*MAX_BLOCKING_LEVEL) ! J=2, PP=-, PR=+
        call states(10)%init(2, -1, -1, 12*MAX_BLOCKING_LEVEL) ! J=2, PP=-, PR=-
        call momentum_states(1)%init(0, 1, 1, 47*MAX_BLOCKING_LEVEL) ! J=0, P=+, q=1
        call momentum_states(2)%init(0, -1, 1, 37*MAX_BLOCKING_LEVEL) ! J=0, P=-, q=1
        call momentum_states(3)%init(1, 1, 1, 66*MAX_BLOCKING_LEVEL) ! J=1, q=1
        call momentum_states(4)%init(2, 1, 1, 47*MAX_BLOCKING_LEVEL) ! J=2, P=+, q=1
        call momentum_states(5)%init(2, -1, 1, 37*MAX_BLOCKING_LEVEL) ! J=2, P=-, q=1
        call momentum_states(6)%init(0, 1, 2, 47*MAX_BLOCKING_LEVEL) ! J=0, P=+, q=2
        call momentum_states(7)%init(0, -1, 2, 37*MAX_BLOCKING_LEVEL) ! J=0, P=-, q=2
        call momentum_states(8)%init(1, 1, 2, 66*MAX_BLOCKING_LEVEL) ! J=1, q=2
        call momentum_states(9)%init(2, 1, 2, 47*MAX_BLOCKING_LEVEL) ! J=2, P=+, q=2
        call momentum_states(10)%init(2, -1, 2, 37*MAX_BLOCKING_LEVEL) ! J=2, P=-, q=2
    else
        ! Allocate arrays
        allocate(gauge_field(NCOL, NCOL, SLICE_VOLUME, 3), &
        gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL), &
        lines_slice(MAX_BLOCKING_LEVEL, 235), &
        momentum_lines_slice(MAX_BLOCKING_LEVEL, 2, 2:235))
    endif

    ! Calculate number of array elements to send and recieve for lines and momentum lines
    num_lines = MAX_BLOCKING_LEVEL * 235
    num_mom_lines = MAX_BLOCKING_LEVEL * 2 * 234

    ! Output calculation info
    if (my_rank == 0) then
        call output_calc_info()
    endif

    call mpi_barrier(mpi_comm_world, mpierr)


    ! ------------------------------- Calculate vevs and correlation matrices ------------------------------- !

    ! Conduct main iteration over configurations
    bin_index = 0
    do config = CONFIG_START, CONFIG_STOP, CONFIG_STEP
        ! Set bin into which measurements will be placed
        bin_index = mod(bin_index, NUM_BINS) + 1

        if (my_rank == 0) then
            ! Recieve lines and momentum_lines for all time slices
            do t = 1, LX4
                ! Recieve lines (with tag 0)
                call mpi_irecv(lines(1, 1, t), num_lines, mpi_double_complex, &
                rank_lookup_table(t), 0, mpi_comm_world, recv_request(2*t-1), mpierr)

                ! Recieve momentum_lines (with tag 1, operator indexing starts at 2)
                call mpi_irecv(momentum_lines(1, 1, 2, t), num_mom_lines, mpi_double_complex, &
                rank_lookup_table(t), 1, mpi_comm_world, recv_request(2*t), mpierr)
            enddo
            call mpi_waitall(2*LX4, recv_request, recv_status, mpierr)
        else
            ! Set directory of file
            write(file_config_id, "(i0)") config
            file_path = trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)

            ! Measure over time slices assigned to this rank
            send_request = mpi_request_null
            do t = my_t_start, my_t_stop
                ! Load configuration from file
                call read_gauge_field_slice(file_path, gauge_field, t)

                ! Block configuration
                gauge_field_blocked = get_blocked_gauge_field(gauge_field)

                ! Wait for previous communication to finish
                call mpi_waitall(2, send_request, send_status, mpierr)

                ! Measure lines
                do blocking_level = 1, MAX_BLOCKING_LEVEL
                    call THERML1(gauge_field_blocked, blocking_level, lines_slice, momentum_lines_slice)
                enddo

                ! Send lines to rank 0 (with tag 0)
                call mpi_isend(lines_slice, num_lines, mpi_double_complex, &
                0, 0, mpi_comm_world, send_request(1), mpierr)

                ! Send momentum lines to rank 0 (with tag 1)
                call mpi_isend(momentum_lines_slice, num_mom_lines, mpi_double_complex, &
                0, 1, mpi_comm_world, send_request(2), mpierr)
            enddo
            call mpi_waitall(2, send_request, send_status, mpierr)
        endif

        ! Wait for all measurements on this configuration to have been recieved
        call mpi_barrier(mpi_comm_world, mpierr)

        ! Ranks 1,...,num_ranks-1 now start measuring the next configuration whilst rank 0 updates the vevs and correlation matrices. Lines on rank 0 are safe as the next measurements have not yet been recieved into lines.
        if (my_rank == 0) then
            ! Update vevs of all states
            do state = 1, 10
                call states(state)%update_vevs(lines, bin_index)
                call momentum_states(state)%update_vevs(momentum_lines, bin_index)
            enddo

            ! Update correlation matrices of all states
            do state = 1, 10
                call states(state)%update_corr_matrix(lines, bin_index)
                call momentum_states(state)%update_corr_matrix(momentum_lines, bin_index)
            enddo
        endif
    enddo


    ! ------------------------------------ Output results of calculation ------------------------------------ !

    ! Output results
    if (my_rank == 0) then
        do state = 1, 10
            call states(state)%output_results()
            call momentum_states(state)%output_results()
        enddo
    endif

    ! Finalise MPI
    call mpi_finalize(mpierr)

    contains

    ! Output calculation info
    subroutine output_calc_info()
        implicit none

        write(*, '(a)') "                                                "
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') "                                                "
        write(*, '(a)') " ******  ******  ***    ******  **      ******  *    *"
        write(*, '(a)') "   **    *    *  *  *   *       **      *    *  **   *"
        write(*, '(a)') "   **    *    *  * *    ******  **      *    *  * *  *"
        write(*, '(a)') "   **    *    *  *  *   *       **      *    *  *  * *"
        write(*, '(a)') "   **    ******  *   *  ******  ******  ******  *   **"
        write(*, '(a)') "                                                "
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        write(*, '(4(a, i0))') "[Info][Lattice Size]     V = ", LX1, "x", LX2, "x", LX3, "x", LX4
        write(*, '(a, i0)') "[Info][Colours]           Number of Colours = ", NCOL
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        write(*, '(a, i0)') "[Info][Measurements]      Number of measurements = ", NCONFIG
        write(*, '(a, a, i0)') &
        "[Info][Measurements]      Starting configuration: ", trim(FILENAME), CONFIG_START
        write(*, '(a, a, i0)') &
        "[Info][Measurements]      Final configuration:    ", trim(FILENAME), CONFIG_STOP
        write(*, '(a, i0)') "[Info][Measurements]      Measurements per bin = ", CONFIG_PER_BIN
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        write(*, '(a)') "[Info][Number of Operators q=0]      J   P   R   #"
        write(*, '(a)') "[Info][Number of Operators q=0]      ----------------"
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      0   +   +   ", states(1)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      0   +   -   ", states(2)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      0   -   +   ", states(3)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      0   -   -   ", states(4)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      1   #   +   ", states(5)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      1   #   -   ", states(6)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      2   +   +   ", states(7)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      2   +   -   ", states(8)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      2   -   +   ", states(9)%get_num_operators()
        write(*, '(a, i0)') "[Info][Number of Operators q=0]      2   -   -   ", states(10)%get_num_operators()
        write(*, '(a)') "[Info][Number of Operators q=0]      ----------------"
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        write(*, '(a)') "[Info][Number of Operators q=1]      J   P   #"
        write(*, '(a)') "[Info][Number of Operators q=1]      ------------"
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=1]      0   +   ", momentum_states(1)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=1]      0   -   ", momentum_states(2)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=1]      1   #   ", momentum_states(3)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=1]      2   +   ", momentum_states(4)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=1]      2   -   ", momentum_states(5)%get_num_operators()
        write(*, '(a)') "[Info][Number of Operators q=1]      ------------"
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
        write(*, '(a)') " *"
        write(*, '(a)') "[Info][Number of Operators q=2]      J   P   #"
        write(*, '(a)') "[Info][Number of Operators q=2]      ------------"
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=2]      0   +   ", momentum_states(6)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=2]      0   -   ", momentum_states(7)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=2]      1   #   ", momentum_states(8)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=2]      2   +   ", momentum_states(9)%get_num_operators()
        write(*, '(a, i0)') &
        "[Info][Number of Operators q=2]      2   -   ", momentum_states(10)%get_num_operators()
        write(*, '(a)') "[Info][Number of Operators q=2]      ------------"
        write(*, '(a)') " *"
        write(*, '(a)') " *******************************************************"
    end subroutine
end program