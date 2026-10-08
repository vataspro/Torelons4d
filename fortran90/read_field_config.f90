module read_field_config
    use torelon_parameters
    use lattice, only: site_index, move, herm
    implicit none

    contains

    ! Read 32-bit integer stored in big endian format
    function read_be_int32(unit) result(value)
        implicit none
        integer, intent(in) :: unit
        integer(int32) :: value
        integer(int8) :: bytes(4)
        integer :: i

        read(unit) bytes
        value = 0_int32
        do i = 1, size(bytes)
        value = ior(shiftl(value, 8), &
                    iand(int(bytes(i), int32), int(z'ff', int32)))
        end do
    end function read_be_int32

    ! Read 64-bit integer stored in big-endian format
    function read_be_real64(unit) result(value)
        implicit none
        integer, intent(in) :: unit
        real(real64) :: value
        integer(int8) :: bytes(8)
        integer(int64) :: bits, octet
        integer :: i

        read(unit) bytes
        bits = 0_int64
        do i = 1, size(bytes)
        octet = iand(int(bytes(i), int64), int(z'ff', int64))
        bits = ior(shiftl(bits, 8), octet)
        end do
        value = transfer(bits, value)
    end function read_be_real64

    function be_bytes_to_real64(bytes) result(value)
        implicit none
        integer(int8), intent(in) :: bytes(8)
        real(real64) :: value
        integer(int64) :: bits, octet
        integer :: i

        bits = 0_int64
        do i = 1, size(bytes)
            octet = iand(int(bytes(i), int64), int(z'ff', int64))
            bits = ior(shiftl(bits, 8), octet)
        end do
        value = transfer(bits, value)
    end function be_bytes_to_real64

    ! Convert quaternion to matrix for reading field configuration
    function quaternion_to_matrix(quaternion)
        implicit none
        real(real64), dimension(4), intent(in) :: quaternion
        complex(real32), dimension(2, 2) :: quaternion_to_matrix

        quaternion_to_matrix(1,1) = real(quaternion(1), kind=real32) + &
        & real(quaternion(4), kind=real32)*(0, 1)
        quaternion_to_matrix(1,2) = -real(quaternion(3), kind=real32) + &
        & real(quaternion(2), kind=real32)*(0, 1)
        quaternion_to_matrix(2,1) = real(quaternion(3), kind=real32) + &
        & real(quaternion(2), kind=real32)*(0, 1)
        quaternion_to_matrix(2,2) = real(quaternion(1), kind=real32) - &
        & real(quaternion(4), kind=real32)*(0, 1)

        return
    end function quaternion_to_matrix

    ! Convert a compressed Sp(2N) gauge link to a full 2Nx2N link
    ! Compressed links are stored as Nx2N 
    function decompress_spN(compressed_link) result(full_link)
        implicit none
        complex(real64), intent(in) :: compressed_link(:,:)
        complex(real64) :: full_link(size(compressed_link,1)*2, &
                                                  size(compressed_link,2))
        integer(int32) :: i, j, newi, newj

        do i=1, NCOL/2
            do j=1, NCOL
                full_link(i, j) = compressed_link(i, j)

                newi = i + NCOL/2
                ! Reading A
                if (j <= NCOL/2) then
                    newj = NCOL/2 + j
                    full_link(newi,newj) = conjg(compressed_link(i,j))
                ! Reading B
                else
                    newj = j - NCOL/2
                    full_link(newi,newj) = -conjg(compressed_link(i,j))
                end if
            end do
        end do

        return
     end function decompress_spN


    ! Read an SpN gauge configuration from a configuration file
    subroutine read_spN_gauge_field(gauge_field_filename, gauge_field)
        implicit none
        character(len=*), intent(in) :: gauge_field_filename
        complex(real64), intent(out) :: gauge_field(:,:,:,:)

        integer :: ios
        character(len=256) :: iomsg

        integer(int32) :: nc_read, nx_read, ny_read, nz_read, nt_read
        integer :: t, x, y, z, dir, dir_target, iun, iq, jq, site
        integer :: mu, nu, site_plus_mu, site_plus_nu, k
        real(real64) :: plaquette_read, plaquette_calc, plaquette_sum, re, im
        complex(real64) :: plaquette_loop(NCOL, NCOL)

        complex(real64), allocatable :: compressed_link(:,:)
        complex(real64), allocatable :: full_link(:,:)

        allocate(compressed_link(NCOL/2, NCOL))
        allocate(full_link(NCOL, NCOL))


        ! Check gauge field variable is the correct size
        if (size(gauge_field, 1) /= NCOL) then
            error stop "gauge_field must have size NCOL in 1st dimension"
        end if
        if (size(gauge_field, 2) /= NCOL) then
            error stop "gauge_field must have size NCOL in 2nd dimension"
        end if
        if (size(gauge_field, 3) /= LATTICE_VOLUME) then
            error stop "gauge_field must have size LATTICE_VOLUME in 3rd dimension"
        end if
        if (size(gauge_field, 4) /= 4) then
            error stop "gauge_field must have size 4 in 4th dimension"
        end if

        ! Open gauge field file
        open(newunit=iun, file=gauge_field_filename, access="stream", form="unformatted", &
        status="old", action="read", iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            error stop "Could not open gauge file: " // trim(iomsg)
        end if

        ! Read file metadata
        nc_read      = read_be_int32(iun)
        nt_read      = read_be_int32(iun)
        nx_read      = read_be_int32(iun)
        ny_read      = read_be_int32(iun)
        nz_read      = read_be_int32(iun)
        plaquette_read = read_be_real64(iun)

        ! Output gauge field metadata
        write(6, "(a, f8.6)") "[I/O][Plaq]    Plaquette value: ", plaquette_read
        WRITE(6, "(a, i2.1)") "[I/O][Ncol]    Number of Colors:", nc_read
        WRITE(6, "(a, 4i3.2)") "[I/O][Dim]    T x X x Y x Z=", nt_read, nx_read, ny_read, nz_read

        ! Read gauge field configuration
        do t = 1, LX4
            do x = 1, LX1
                do y = 1, LX2
                    do z = 1, LX3
                        ! Get index of site referenced by these coordinates
                        site = site_index(x, y, z, t)

                        ! Import all links from this site
                        do dir = 1, 4
                            do iq = 1, NCOL/2
                                do jq = 1, NCOL
                                    re = read_be_real64(iun)
                                    im = read_be_real64(iun)
                                    compressed_link(iq, jq) = cmplx(re, im)
                                end do
                            end do

                            full_link = decompress_spN(compressed_link)

                            if (dir == 1) then
                                dir_target = 4
                            else
                                dir_target = dir-1
                            end if

                            gauge_field(:, :, site, dir_target) = full_link

                        end do
                    end do
                end do
            end do
        end do

        ! Recompute the average plaquette from the loaded links as a read check.
        plaquette_sum = 0.0_real64
        do mu = 1, 4
            do nu = mu + 1, 4
                do site = 1, LATTICE_VOLUME
                    site_plus_mu = move(site, mu)
                    site_plus_nu = move(site, nu)
                    plaquette_loop = matmul(matmul(matmul( &
                        gauge_field(:, :, site, mu), &
                        gauge_field(:, :, site_plus_mu, nu)), &
                        herm(gauge_field(:, :, site_plus_nu, mu))), &
                        herm(gauge_field(:, :, site, nu)))
                    do k = 1, NCOL
                        plaquette_sum = plaquette_sum + real(plaquette_loop(k, k), kind=real64)
                    end do
                end do
            end do
        end do
        plaquette_calc = plaquette_sum / &
            (6.0_real64 * real(NCOL, real64) * real(LATTICE_VOLUME, real64))
        write(6, '(a, f16.12)') '[I/O][Plaq]    Recomputed plaquette: ', plaquette_calc
        write(6, '(a, es12.4)') '[I/O][Plaq]    Difference from header: ', plaquette_calc - plaquette_read

        close(iun)
    end subroutine read_spN_gauge_field


    ! Read gauge configuration from configuration file - renamed from READ_GF
    subroutine read_gauge_field(gauge_field_filename, gauge_field)
        implicit none
        character(len=*), intent(in) :: gauge_field_filename
        complex(real64), intent(out) :: gauge_field(:,:,:,:)

        integer :: ios
        character(len=256) :: iomsg

        integer(int32) :: nc_read, nx_read, ny_read, nz_read, nt_read
        integer :: t, x, y, z, dir, dir_target, iun, iq, site
        real(real64) :: plaquette_read

        real(real64) :: quaternion(4)

        ! Check gauge field variable is the correct size
        if (size(gauge_field, 1) /= NCOL) then
            error stop "gauge_field must have size NCOL in 1st dimension"
        end if
        if (size(gauge_field, 2) /= NCOL) then
            error stop "gauge_field must have size NCOL in 2nd dimension"
        end if
        if (size(gauge_field, 3) /= LATTICE_VOLUME) then
            error stop "gauge_field must have size LATTICE_VOLUME in 3rd dimension"
        end if
        if (size(gauge_field, 4) /= 4) then
            error stop "gauge_field must have size 4 in 4th dimension"
        end if

        ! Open gauge field file
        open(newunit=iun, file=gauge_field_filename, access="stream", form="unformatted", &
        status="old", action="read", iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            error stop "Could not open gauge file: " // trim(iomsg)
        end if

        ! Read file metadata
        nc_read      = read_be_int32(iun)
        nt_read      = read_be_int32(iun)
        nx_read      = read_be_int32(iun)
        ny_read      = read_be_int32(iun)
        nz_read      = read_be_int32(iun)
        plaquette_read = read_be_real64(iun)

        ! Output gauge field metadata
        write(6, "(a, f8.6)") "[I/O][Plaq]    Plaquette value: ", plaquette_read
        WRITE(6, "(a, i2.1)") "[I/O][Ncol]    Number of Colors:", nc_read
        WRITE(6, "(a, 4i3.2)") "[I/O][Dim]    T x X x Y x Z=", nt_read, nx_read, ny_read, nz_read

        ! Protect against bad parameter file / loaded configuration mismatch
        if (nc_read /= NCOL .or. nt_read /= LX4 .or. nx_read /= LX1 .or. &
            ny_read /= LX2 .or. nz_read /= LX3) then
            close(iun)
            error stop "Gauge file dimensions do not match configured lattice"
        end if


        ! Read gauge field configuration
        do t = 1, LX4
            do x = 1, LX1
                do y = 1, LX2
                    do z = 1, LX3
                        ! Get index of site referenced by these coordinates
                        site = site_index(x, y, z, t)

                        ! Import all links from this site
                        do dir = 1, 4
                            do iq = 1, size(quaternion)
                                quaternion(iq) = read_be_real64(iun)
                            end do

                            if (dir == 1) then
                                dir_target = 4
                            else
                                dir_target = dir-1
                            end if

                            gauge_field(1, 1, site, dir_target) = &
                            cmplx(real(quaternion(1), kind=real64),real(quaternion(4), kind=real64), &
                            kind=real64)
                            gauge_field(1, 2, site, dir_target) = &
                            cmplx(-real(quaternion(3), kind=real64),real(quaternion(2), kind=real64), &
                            kind=real64)
                            gauge_field(2, 1, site, dir_target) = &
                            cmplx(real(quaternion(3), kind=real64),real(quaternion(2), kind=real64), &
                            kind=real64)
                            gauge_field(2, 2, site, dir_target) = &
                            cmplx(real(quaternion(1), kind=real64),-real(quaternion(4), kind=real64), &
                            kind=real64)

                        end do
                    end do
                end do
            end do
        end do

        close(iun)
    end subroutine read_gauge_field

    ! Read one spatial time slice of an Sp(2N) gauge configuration.
    ! As in read_gauge_field_slice, file directions 2:4 map to in-memory
    ! directions 1:3. Each compressed link contains NCOL/2 x NCOL complex
    ! values, with each complex value stored as two big-endian real64 values.
    subroutine read_spN_gauge_field_slice(gauge_field_filename, gauge_field_slice, time_slice)
        implicit none
        character(len=*), intent(in) :: gauge_field_filename
        complex(real64), intent(out) :: gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: time_slice

        integer :: ios, iun
        integer :: x, y, z, dir, iq, jq, site, dir_target, byte_start
        integer(int32) :: nc_read, nx_read, ny_read, nz_read, nt_read
        integer(int64) :: file_site_number, position, site_record_bytes, link_bytes
        integer(int8), allocatable :: spatial_link_bytes(:)
        integer(int8) :: real_bytes(8)
        character(len=256) :: iomsg
        real(real64) :: plaquette_read, re, im
        complex(real64), allocatable :: compressed_link(:,:), full_link(:,:)

        if (time_slice < 1 .or. time_slice > LX4) then
            error stop "time_slice must be between 1 and LX4"
        end if

        if (NCOL < 2 .or. mod(NCOL, 2) /= 0) then
            error stop "Sp(2N) requires an even NCOL of at least 2"
        end if

        allocate(compressed_link(NCOL/2, NCOL), full_link(NCOL, NCOL))
        link_bytes = int(NCOL, int64) * int(NCOL, int64) * 8_int64
        site_record_bytes = 4_int64 * link_bytes
        allocate(spatial_link_bytes(3 * link_bytes))

        open(newunit=iun, file=gauge_field_filename, access="stream", form="unformatted", &
            status="old", action="read", iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            error stop "Could not open gauge file: " // trim(iomsg)
        end if

        nc_read = read_be_int32(iun)
        nt_read = read_be_int32(iun)
        nx_read = read_be_int32(iun)
        ny_read = read_be_int32(iun)
        nz_read = read_be_int32(iun)
        plaquette_read = read_be_real64(iun)

        if (nc_read /= NCOL .or. nt_read /= LX4 .or. nx_read /= LX1 .or. &
            ny_read /= LX2 .or. nz_read /= LX3) then
            close(iun)
            error stop "Gauge file dimensions do not match configured lattice"
        end if

        do x = 1, LX1
            do y = 1, LX2
                do z = 1, LX3
                    site = site_index(x, y, z)
                    ! File order is t,x,y,z (z fastest); skip the temporal link.
                    file_site_number = int((time_slice - 1) * SLICE_VOLUME + &
                        (x - 1) * LX2 * LX3 + (y - 1) * LX3 + z - 1, int64)
                    ! Header is 28 bytes; spatial links start after link 1.
                    position = 29_int64 + file_site_number * site_record_bytes + link_bytes
                    read(iun, pos=position, iostat=ios, iomsg=iomsg) spatial_link_bytes
                    if (ios /= 0) then
                        close(iun)
                        error stop "Could not seek to Sp(2N) gauge field slice: " // trim(iomsg)
                    end if

                    do dir = 1, 3
                        dir_target = dir
                        do iq = 1, NCOL/2
                            do jq = 1, NCOL
                                byte_start = (dir - 1) * link_bytes + &
                                    ((iq - 1) * NCOL + jq - 1) * 16 + 1
                                real_bytes = spatial_link_bytes(byte_start:byte_start + 7)
                                re = be_bytes_to_real64(real_bytes)
                                real_bytes = spatial_link_bytes(byte_start + 8:byte_start + 15)
                                im = be_bytes_to_real64(real_bytes)
                                compressed_link(iq, jq) = cmplx(re, im, kind=real64)
                            end do
                        end do

                        full_link = decompress_spN(compressed_link)
                        gauge_field_slice(:, :, site, dir_target) = full_link
                    end do
                end do
            end do
        end do

        close(iun)
    end subroutine read_spN_gauge_field_slice

    ! Read single time slice of the gauge configuration
    subroutine read_gauge_field_slice(gauge_field_filename, gauge_field_slice, time_slice)
        implicit none
        character(len=*), intent(in) :: gauge_field_filename
        complex(real64), intent(out) :: gauge_field_slice(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: time_slice

        integer :: ios, iun
        integer :: x, y, z, dir, iq, site, dir_target, byte_start
        integer(int32) :: nc_read, nx_read, ny_read, nz_read, nt_read
        integer(int64) :: file_site_number, position
        integer(int8) :: spatial_link_bytes(96)
        integer(int8) :: real_bytes(8)
        character(len=256) :: iomsg
        real(real64) :: plaquette_read, quaternion(4)

        if (time_slice < 1 .or. time_slice > LX4) then
            error stop "time_slice must be between 1 and LX4"
        end if

        open(newunit=iun, file=gauge_field_filename, access="stream", form="unformatted", &
        status="old", action="read", iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            error stop "Could not open gauge file: " // trim(iomsg)
        end if

        ! The header is 5 big-endian 32-bit integers followed by one real64.
        nc_read = read_be_int32(iun)
        nt_read = read_be_int32(iun)
        nx_read = read_be_int32(iun)
        ny_read = read_be_int32(iun)
        nz_read = read_be_int32(iun)
        plaquette_read = read_be_real64(iun)

        if (nc_read /= NCOL .or. nt_read /= LX4 .or. nx_read /= LX1 .or. &
            ny_read /= LX2 .or. nz_read /= LX3) then
            close(iun)
            error stop "Gauge file dimensions do not match configured lattice"
        end if

        ! Each site occupies four directions x four real64 quaternion values.
        ! Skip file direction 1 (temporal); read only spatial directions 2:4.
        do x = 1, LX1
            do y = 1, LX2
                do z = 1, LX3
                    site = site_index(x, y, z)
                    ! The configuration file ordering used by read_gauge_field is
                    ! t, x, y, z (with z varying fastest), unlike site_index's
                    ! in-memory ordering where x varies fastest.
                    file_site_number = int((time_slice - 1) * SLICE_VOLUME + &
                        (x - 1) * LX2 * LX3 + (y - 1) * LX3 + z - 1, int64)
                    ! POS is one-based: header=28 bytes, then site records of 128 bytes;
                    ! spatial links begin 32 bytes into each site's record.
                    position = 29_int64 + file_site_number * 128_int64 + 32_int64
                    read(iun, pos=position, iostat=ios, iomsg=iomsg) spatial_link_bytes
                    if (ios /= 0) then
                        close(iun)
                        error stop "Could not seek to gauge field slice: " // trim(iomsg)
                    end if

                    do dir = 2, 4
                        dir_target = dir - 1
                        do iq = 1, 4
                            byte_start = (dir - 2) * 32 + (iq - 1) * 8 + 1
                            real_bytes = spatial_link_bytes(byte_start:byte_start + 7)
                            quaternion(iq) = be_bytes_to_real64(real_bytes)
                        end do
                        gauge_field_slice(1, 1, site, dir_target) = &
                            cmplx(quaternion(1), quaternion(4), kind=real64)
                        gauge_field_slice(1, 2, site, dir_target) = &
                            cmplx(-quaternion(3), quaternion(2), kind=real64)
                        gauge_field_slice(2, 1, site, dir_target) = &
                            cmplx(quaternion(3), quaternion(2), kind=real64)
                        gauge_field_slice(2, 2, site, dir_target) = &
                            cmplx(quaternion(1), -quaternion(4), kind=real64)
                    end do
                end do
            end do
        end do

        close(iun)
    end subroutine
end module read_field_config
