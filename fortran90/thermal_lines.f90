module thermal_lines
    use parameters
    use lattice
    implicit none

    complex(real64), allocatable :: squares_up(:,:,:,:,:), squares_down(:,:,:,:,:)
    integer :: direction_index(3)

    contains

    ! Calculate square pulse from given site in given positive direction
    function square_pulse_up(gauge_field, blocking_level, site, flux_direction, pulse_direction)
        implicit none
        complex(real64), intent(in) :: gauge_field(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: blocking_level, site, flux_direction, pulse_direction
        complex(real64) :: square_pulse_up(NCOL, NCOL)

        integer :: site_plus_flux, site_plus_pulse

        ! Check flux and pulse directions are compatibile
        if (abs(flux_direction) == abs(pulse_direction)) error stop &
        "Flux direction cannot equal pulse direction in square_pulse_up"

        ! Check flux points in positive direction
        if (flux_direction <= 0) error stop &
        "Flux direction must be passed as positive in square_pulse_up"

        ! Check deformation points in positive direction
        if (pulse_direction <= 0) error stop &
        "Pulse direction must be passed as positive in square_pulse_up"

        ! Check both directions point in spatial directions
        if ((flux_direction > 3).or.(pulse_direction > 3)) error stop &
        "Flux and pulse directions must both be spatial directions"

        ! Get next site along flux tube
        site_plus_flux = move(site, flux_direction, blocking_level)

        ! Get site in direction of square pulse deformation
        site_plus_pulse = move(site, pulse_direction, blocking_level)

        ! Get square pulse
        square_pulse_up = matmul(matmul(&
            gauge_field(:,:,site,pulse_direction), &
            gauge_field(:,:,site_plus_pulse,flux_direction)), &
            herm(gauge_field(:,:,site_plus_flux,pulse_direction)))
    end function

    ! Calculate square pulse from given site in given negative direction
    function square_pulse_down(gauge_field, blocking_level, site, flux_direction, pulse_direction)
        implicit none
        complex(real64), intent(in) :: gauge_field(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: blocking_level, site, flux_direction, pulse_direction
        complex(real64) :: square_pulse_down(NCOL, NCOL)

        integer :: site_minus_pulse, corner_site

        ! Check flux and pulse directions are compatibile
        if (abs(flux_direction) == abs(pulse_direction)) error stop &
        "Flux direction cannot equal pulse direction in square_pulse_down"

        ! Check flux points in positive direction
        if (flux_direction <= 0) error stop &
        "Flux direction must be passed as positive in square_pulse_down"

        ! Check deformation points in positive direction
        if (pulse_direction <= 0) error stop &
        "Pulse direction must be passed as positive in square_pulse_down"

        ! Check both directions point in spatial directions
        if ((flux_direction > 3).or.(pulse_direction > 3)) error stop &
        "Flux and pulse directions must both be spatial directions"

        ! Get site in direction of square pulse deformation
        site_minus_pulse = move(site, -pulse_direction, blocking_level)

        ! Get site on opposite corner of square pulse from site
        corner_site = move(site_minus_pulse, flux_direction, blocking_level)

        ! Get square pulse
        square_pulse_down = matmul(matmul(&
            herm(gauge_field(:,:,site_minus_pulse,pulse_direction)), &
            gauge_field(:,:,site_minus_pulse,flux_direction)), &
            gauge_field(:,:,corner_site,pulse_direction))
    end function

    ! Calculate all square pulses in orthogonal directions over all blocking levels on the lattice
    subroutine get_square_pulses(gauge_field_blocked, flux_direction)
        implicit none
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(in) :: flux_direction

        integer :: site, pulse_direction, blocking_level, dir_count

        ! Check if square pulse containers are allocated. If not, allocate them.
        if (.not.allocated(squares_up)) &
        allocate(squares_up(NCOL, NCOL, SLICE_VOLUME, 2, MAX_BLOCKING_LEVEL))
        if (.not.allocated(squares_down)) &
        allocate(squares_down(NCOL, NCOL, SLICE_VOLUME, 2, MAX_BLOCKING_LEVEL))

        ! Iterate over lattice and calculate all square pulses in directions orthogonal to the flux direction
        dir_count = 1
        direction_index = 0
        do pulse_direction = 1, 3
            ! Skip direction parallel to flux direction
            if (pulse_direction == flux_direction) cycle

            ! Set direction index. This is a small array to map between pulse_direction and the position of the relevant entry in the 4th dimension of the squares_up/down arrays
            direction_index(pulse_direction) = dir_count

            ! Iterate over all blocking levels
            do blocking_level = 1, MAX_BLOCKING_LEVEL
                ! Iterate over all sites
                do site = 1, SLICE_VOLUME
                    ! Calculate square pulses
                    squares_up(:, :, site, dir_count, blocking_level) &
                    = square_pulse_up(gauge_field_blocked(:, :, :, :, blocking_level), &
                                      blocking_level, site, flux_direction, pulse_direction)
                    squares_down(:, :, site, dir_count, blocking_level) &
                    = square_pulse_down(gauge_field_blocked(:, :, :, :, blocking_level), &
                                        blocking_level, site, flux_direction, pulse_direction)
                enddo
            enddo

            ! Iterate dir_count
            dir_count = dir_count + 1
        enddo
    end subroutine

    ! Return square pulse from given site, in given direction
    function square_pulse(site, pulse_direction, blocking_level)
        implicit none
        integer, intent(in) :: site, pulse_direction, blocking_level
        complex(real64) :: square_pulse(NCOL, NCOL)

        ! Check pulse direction is not parallel to flux direction
        if (direction_index(abs(pulse_direction)) == 0) error stop &
        "Pulse direction must not be parallel to flux direction in square_pulse"

        ! Check pulse direction is a valid direction
        if (pulse_direction == 0 .or. abs(pulse_direction) > 3) error stop &
        "Pulse direction must point in a spatial direction in square_pulse"

        ! Return square pulse in given direction by looking up value in squares_up or squares_down
        if (pulse_direction > 0) then
            square_pulse = squares_up(:, :, site, direction_index(pulse_direction), blocking_level)
        else
            square_pulse = squares_down(:, :, site, direction_index(abs(pulse_direction)), blocking_level)
        endif
    end function

    !##############################################
    !           SQUARE PULSES
    !        __
    ! UP: __|  |__    DOWN: __    __
    !                         |__|
    !
    !##############################################
    ! Loop 1: Up square pulse Y
    subroutine loop_1(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(A11, square_pulse(in_site, 2, blocking_level))
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    ! Loop 2: Up square pulse Z
    subroutine loop_2(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(A11, square_pulse(in_site, 3, blocking_level))
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    ! Loop 3: Down square pulse Y
    subroutine loop_3(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(A11, square_pulse(in_site, -2, blocking_level))
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    ! Loop 4: Down square pulse Z
    subroutine loop_4(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(A11, square_pulse(in_site, -3, blocking_level))
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    !##############################################
    !          UP-UP SQUARE PULSES
    !    __  __            ____
    ! __|  ||  |__  ->  __|    |__
    !
    !##############################################
    ! Loop 5: up-up square pulse Y
    subroutine loop_5(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: middle_site

        middle_site = move(in_site, 1, blocking_level)

        A11 = matmul(A11, matmul( &
            square_pulse(in_site, 2, blocking_level), &
            square_pulse(middle_site, 2, blocking_level) &
        ))

        out_site = move(middle_site, 1, blocking_level)
    end subroutine loop_5

    ! Loop 6: up-up square pulse Z
    subroutine loop_6(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: middle_site

        middle_site = move(in_site, 1, blocking_level)

        A11 = matmul(A11, matmul( &
            square_pulse(in_site, 3, blocking_level), &
            square_pulse(middle_site, 3, blocking_level) &
        ))

        out_site = move(middle_site, 1, blocking_level)
    end subroutine loop_6

    !##############################################
    !        DOWN-DOWN SQUARE PULSES
    ! __        __        __      __
    !   |__||__|     ->     |____|
    !
    !##############################################
    ! Loop 7: down-down square pulse Y
    subroutine loop_7(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: middle_site

        middle_site = move(in_site, 1, blocking_level)

        A11 = matmul(A11, matmul( &
            square_pulse(in_site, -2, blocking_level), &
            square_pulse(middle_site, -2, blocking_level) &
        ))

        out_site = move(middle_site, 1, blocking_level)
    end subroutine loop_7

    ! Loop 8: down-down square pulse Z
    subroutine loop_8(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: middle_site

        middle_site = move(in_site, 1, blocking_level)

        A11 = matmul(A11, matmul( &
            square_pulse(in_site, -3, blocking_level), &
            square_pulse(middle_site, -3, blocking_level) &
        ))

        out_site = move(middle_site, 1, blocking_level)
    end subroutine loop_8

    !##############################################
    !          UP-DOWN SQUARE PULSES
    !    __
    ! __|  |   __
    !      |__|
    !
    !##############################################
    ! Loop 9: up-down square pulse Y
    subroutine loop_9(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: i

        A11 = matmul(square_pulse(in_site, 2, blocking_level), &
        square_pulse(move(in_site, 1, blocking_level), -2, blocking_level))

        out_site = in_site
        do  i = 1, 2
            out_site = move(out_site, 1, blocking_level)
        enddo
    end subroutine

    ! Loop 10: up-down square pulse Z
    subroutine loop_10(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: middle_site

        middle_site = move(in_site, 1, blocking_level)

        A11 = matmul(A11, matmul( &
            square_pulse(in_site, 3, blocking_level), &
            square_pulse(middle_site, -3, blocking_level) &
        ))

        out_site = move(middle_site, 1, blocking_level)
    end subroutine loop_10

    !##############################################
    !          DOWN-UP SQUARE PULSES
    !       __
    ! __   |  |__
    !   |__|
    !
    !##############################################
    ! Loop 11: down-up square pulse Y
    subroutine loop_11(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: middle_site

        middle_site = move(in_site, 1, blocking_level)

        A11 = matmul(A11, matmul( &
            square_pulse(in_site, -2, blocking_level), &
            square_pulse(middle_site, 2, blocking_level) &
        ))

        out_site = move(middle_site, 1, blocking_level)
    end subroutine loop_11

    ! Loop 12: down-up square pulse Z
    subroutine loop_12(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer :: middle_site

        middle_site = move(in_site, 1, blocking_level)

        A11 = matmul(A11, matmul( &
            square_pulse(in_site, -3, blocking_level), &
            square_pulse(middle_site, 3, blocking_level) &
        ))

        out_site = move(middle_site, 1, blocking_level)
    end subroutine loop_12



end module