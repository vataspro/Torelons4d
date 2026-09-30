program blocking_example
    use torelon_parameters
    use lattice
    use read_field_config

    implicit none

    complex(real64), allocatable :: gauge_field(:,:,:,:), &
    gauge_field_slice(:,:,:,:), gauge_field_slice_blocked(:,:,:,:,:)
    character(len=16) :: file_config_id
    character(len=256) :: directory
    real :: start, finish

    ! Load parameters
    call initialise_parameters("parameter_file_local.txt")

    ! Setup lattice
    call setup_lattice()

    ! Allocate arrays
    allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 4), &
    gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3), &
    gauge_field_slice_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL))

    ! Set directory path
    write(file_config_id, '(i0)') CONFIG_START
    directory=trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)

    ! Load gauge field into memory
    call read_gauge_field(directory, gauge_field)

    ! Take time slice from gauge field
    gauge_field_slice = gauge_field(:, :, 1:SLICE_VOLUME, 1:3)

    ! Block gauge field
    call cpu_time(start)
    gauge_field_slice_blocked = get_blocked_gauge_field(gauge_field_slice)
    call cpu_time(finish)

    ! Output time taken
    write(*, '(a, f7.4, a)') "Blocking time:    ", finish-start, " seconds"
end program