program main
    use torelon_parameters
    use lattice
    use read_field_config
    use here_be_dragons
    use states_class

    implicit none

    ! ---------------------------------------- Initialise variables ---------------------------------------- !

    ! File access variables
    logical :: file_exists
    character(len=16) :: file_config_id
    character(len=256) :: file_path

    ! Gauge field
    complex(real64), allocatable :: gauge_field(:,:,:,:), gauge_field_slice_blocked(:,:,:,:,:)

    ! Measurement variables
    complex(real64), allocatable :: lines(:,:,:), momentum_lines(:,:,:,:)
    integer :: bin_index

    ! State variables
    type(torelon_state) :: states(10)
    type(torelon_momentum_state) :: momentum_states(10)

    ! Timing variables
    real(real32) :: start, finish, &
    avg_runtime_loading, avg_runtime_blocking, avg_runtime_measurement, avg_runtime_correlation

    ! Loop variables
    integer :: config, t, blocking_level, state


    ! ----------------------------------------- Set up calculation ----------------------------------------- !

    ! Load all parameters from file
    call initialise_parameters("parameters.txt")

    ! Check first configuration can be found
    write(file_config_id, "(i0)") CONFIG_START
    file_path = trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)
    inquire(file=trim(file_path), exist=file_exists)
    if (.not.file_exists) then
        write(*, '(a, a)') "[Error][File access] Cannot access file ", file_path
        stop
    endif

    ! Setup lattice movers
    call setup_lattice()

    ! Allocate variables
    allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 4), &
    gauge_field_slice_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL))
    allocate(lines(LX4, MAX_BLOCKING_LEVEL, 235), &
    momentum_lines(LX4, MAX_BLOCKING_LEVEL, 2, 2:235))

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

    ! Output calculation info
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
    write(*, '(4(a, i0))') "[Info][Lattice Size]                 V = ", LX1, "x", LX2, "x", LX3, "x", LX4
    write(*, '(a, i0)') "[Info][Colours]                      Number of Colours = ", NCOL
    write(*, '(a)') " *"
    write(*, '(a)') " *******************************************************"
    write(*, '(a)') " *"
    write(*, '(a, i0)') "[Info][Measurements]                 Number of measurements = ", NCONFIG
    write(*, '(a, a, i0)') &
    "[Info][Measurements]                 Starting configuration: ", FILENAME, CONFIG_START
    write(*, '(a, i0)') "[Info][Measurements]                 Measurements per bin = ", CONFIG_PER_BIN
    write(*, '(a)') " *"
    write(*, '(a)') " *******************************************************"
    write(*, '(a)') " *"
    write(*, '(a)') "[Info][Number of Operators q=0]      J   P   R      #"
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
    write(*, '(a)') "[Info][Number of Operators q=1]      J   P      #"
    write(*, '(a)') "[Info][Number of Operators q=1]      ------------"
    write(*, '(a, i0)') "[Info][Number of Operators q=1]      0   +   ", momentum_states(1)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=1]      0   -   ", momentum_states(2)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=1]      1   #   ", momentum_states(3)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=1]      2   +   ", momentum_states(4)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=1]      2   -   ", momentum_states(5)%get_num_operators()
    write(*, '(a)') "[Info][Number of Operators q=1]      ------------"
    write(*, '(a)') " *"
    write(*, '(a)') " *******************************************************"
    write(*, '(a)') " *"
    write(*, '(a)') "[Info][Number of Operators q=2]      J   P      #"
    write(*, '(a)') "[Info][Number of Operators q=2]      ------------"
    write(*, '(a, i0)') "[Info][Number of Operators q=2]      0   +   ", momentum_states(6)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=2]      0   -   ", momentum_states(7)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=2]      1   #   ", momentum_states(8)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=2]      2   +   ", momentum_states(9)%get_num_operators()
    write(*, '(a, i0)') "[Info][Number of Operators q=2]      2   -   ", momentum_states(10)%get_num_operators()
    write(*, '(a)') "[Info][Number of Operators q=2]      ------------"
    write(*, '(a)') " *"
    write(*, '(a)') " *******************************************************"


    ! ------------------------------- Calculate vevs and correlation matrices ------------------------------- !

    ! Initilaise runtime variables to zero
    avg_runtime_loading = 0.0
    avg_runtime_blocking = 0.0
    avg_runtime_measurement = 0.0
    avg_runtime_correlation = 0.0

    ! Conduct main iteration over configurations
    bin_index = 0
    do config = CONFIG_START, CONFIG_STOP, CONFIG_STEP
        ! Set bin into which measurements will be placed
        bin_index = mod(bin_index, NUM_BINS) + 1

        ! Set directory of file
        write(file_config_id, "(i0)") config
        file_path = trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)

        ! Load gauge field
        call cpu_time(start)
        call read_gauge_field(file_path, gauge_field)
        call cpu_time(finish)
        avg_runtime_loading = avg_runtime_loading + (finish - start)

        ! Conduct measurements over every time slice
        do t = 1, LX4
            ! Block gauge field
            call cpu_time(start)
            gauge_field_slice_blocked = &
            get_blocked_gauge_field(gauge_field(:, :, (t-1)*SLICE_VOLUME+1 : t*SLICE_VOLUME, 1:3))
            call cpu_time(finish)
            avg_runtime_blocking = avg_runtime_blocking + (finish - start)

            ! Measure thermal lines over all blocking levels
            do blocking_level = 1, MAX_BLOCKING_LEVEL
                call cpu_time(start)
                call THERML1(gauge_field_slice_blocked, t, blocking_level, lines, momentum_lines)
                call cpu_time(finish)
                avg_runtime_measurement = avg_runtime_measurement + (finish - start)
            enddo
        enddo

        ! Update vevs of all states
        do state = 1, 10
            call cpu_time(start)
            call states(state)%update_vevs(lines, bin_index)
            call momentum_states(state)%update_vevs(lines, bin_index)
            call cpu_time(finish)
            avg_runtime_correlation = avg_runtime_correlation + (finish - start)
        enddo

        ! Update correlation matrices of all states
        do state = 1, 10
            call cpu_time(start)
            call states(state)%update_corr_matrix(lines, bin_index)
            call states(state)%update_corr_matrix(lines, bin_index)
            call cpu_time(finish)
            avg_runtime_correlation = avg_runtime_correlation + (finish - start)
        enddo
    enddo


    ! ------------------------------------ Output results of calculation ------------------------------------ !

    ! Get average runtimes
    avg_runtime_loading = avg_runtime_loading / NCONFIG
    avg_runtime_blocking = avg_runtime_blocking / NCONFIG
    avg_runtime_measurement = avg_runtime_measurement / NCONFIG
    avg_runtime_correlation = avg_runtime_correlation / NCONFIG

    ! Output runtimes
    write(*, '(a)') " *"
    write(*, '(a, f0.2)') &
    "[Info][Runtimes]                     Average runtime for loading a configuration    ", avg_runtime_loading
    write(*, '(a, f0.2)') &
    "[Info][Runtimes]                     Average runtime for blocking a configuration    ", &
    avg_runtime_blocking
    write(*, '(a, f0.2)') &
    "[Info][Runtimes]                     Average runtime for measuring lines on a configuration    ", &
    avg_runtime_loading
    write(*, '(a, f0.2)') &
    "[Info][Runtimes]                     Average runtime for updating vevs and correlators    ", avg_runtime_correlation
    write(*, '(a)') " *"
    write(*, '(a)') " *******************************************************"
    write(*, '(a)') " *"

    ! Output results
    do state = 1, 10
        call states(state)%output_results()
        call momentum_states(state)%output_results()
    enddo
end program