module read_field_config
    use parameters
    use lattice
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
        if (size(gauge_field, dim=1) /= NCOL) then
            error stop "gauge_field must have size NCOL in 1st dimension"
        end if
        if (size(gauge_field, dim=2) /= NCOL) then
            error stop "gauge_field must have size NCOL in 2nd dimension"
        end if
        if (size(gauge_field, dim=3) /= LATTICE_VOLUME) then
            error stop "gauge_field must have size LATTICE_VOLUME in 3rd dimension"
        end if
        if (size(gauge_field, dim=4) /= 4) then
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
end module read_field_config
