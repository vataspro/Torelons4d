program blocking_example
    use parameters
    use lattice
    use read_field_config

    implicit none

    complex(real64), allocatable :: gauge_field_read(:,:,:,:,:,:,:), gauge_field(:,:,:,:,:), &
    gauge_field_slice(:,:,:,:), gauge_field_slice_smeared(:,:,:,:)
    character(len=16) :: file_config_id
    character(len=256) :: directory

    ! Load parameters
    call initialise_parameters("parameter_file.txt")

    ! Allocate arrays
    allocate(gauge_field_read(NCOL, NCOL, LX1, LX2, LX3, LX4, 4), &
    gauge_field(NCOL, NCOL, SLICE_VOLUME, LX4, 4), &
    gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3), &
    gauge_field_slice_smeared(NCOL, NCOL, SLICE_VOLUME, 3))

    ! Set directory path
    write(file_config_id, '(i0)') CONFIG_START
    directory=trim(FILEPATH) // trim(FILENAME) // trim(file_config_id)

    ! Load gauge field into memory
    call read_gauge_field(directory, gauge_field_read)

    ! Reshape gauge field
    gauge_field = reshape(gauge_field_read, [NCOL, NCOL, SLICE_VOLUME, LX4, 4])

    ! Take time slice from gauge field
    gauge_field_slice = gauge_field(:, :, :, 1, 1:3)

    ! Smear gauge field
    gauge_field_slice_smeared = smear(gauge_field_slice, 1)
end program