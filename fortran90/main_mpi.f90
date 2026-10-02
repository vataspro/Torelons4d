program main_mpi
    use mpi_f08
    use torelon_parameters
    use lattice

    implicit none

    ! ---------------------------------------- Initialise variables ---------------------------------------- !

    ! Parallelisation variables
    integer :: my_rank, num_ranks, mpierr, my_LX4, my_t_start, my_t_stop, local_LX4
    integer, allocatable :: rank_lookup_table(:)
    type(mpi_status) :: send_status
    type(mpi_status), allocatable :: recv_status(:)
    type(mpi_request) :: send_request
    type(mpi_request), allocatable :: recv_request(:)
    logical :: test_success

    ! File access variables
    logical :: file_exists
    character(len=16) :: file_config_id
    character(len=256) :: file_path

    ! Gauge field
    complex(real64), allocatable :: gauge_field(:,:,:,:), gauge_field_blocked(:,:,:,:,:)

    ! Measurement variables
    integer, allocatable :: mock_lines(:)
    integer :: mock_vev

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
                stop
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
        ! write(*, '(4(a, i0))') "Rank ", my_rank, ":    t_start:    ", my_time_slice_start, "    t_stop:    ", &
        ! my_time_slice_start + my_LX4 - 1, "    my_LX4:    ", my_LX4
    endif

    ! call mpi_barrier(mpi_comm_world, mpierr)

    ! if (my_rank == 0) then
    !     do t = 1, LX4
    !         write(*, '(2(a, i0))') "t = ", t, "    rank = ", rank_lookup_table(t)
    !     enddo
    ! endif

    ! Setup lattice movers
    call setup_lattice()

    ! ! Allocate arrays
    ! if (my_rank == 0) then
    !     allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 4))
    ! else
    !     allocate(gauge_field(NCOL, NCOL, SLICE_VOLUME, 3), &
    !     gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL))
    ! endif

    call mpi_barrier(mpi_comm_world, mpierr)


    ! ------------------------------- Calculate vevs and correlation matrices ------------------------------- !

    ! Have each rank 1,...,num_ranks-1 send its time slices to rank 0
    if (my_rank == 0) then
        allocate(mock_lines(LX4))
        allocate(recv_status(LX4), recv_request(LX4))
    else
        allocate(mock_lines(my_LX4))
    endif
    mock_lines = 0
    mock_vev = 0

    if (my_rank == 0) then
        ! Recieve messages for all time slices
        do t = 1, LX4
            call mpi_irecv(mock_lines(t), 1, mpi_int, rank_lookup_table(t), &
            0, mpi_comm_world, recv_request(t), mpierr)
        enddo
        call mpi_waitall(LX4, recv_request, recv_status)
    else
        ! Send messages to rank 0
        do t = my_t_start, my_t_stop
            call mpi_isend(t, 1, mpi_int, 0, 0, mpi_comm_world, send_request, mpierr)
            call mpi_wait(send_request, send_status)
        enddo
    endif

    if (my_rank == 0) then
        test_success = .true.
        do t = 1, LX4
            if (mock_lines(t) /= t) then
                test_success = .false.
                exit
            endif
        enddo

        if (test_success) then
            print *, "Success!"
        else
            print *, "Test failed"
        endif
    endif

    ! Finalise MPI
    call mpi_finalize(mpierr)
end program