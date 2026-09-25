module thermal_lines
    use parameters
    use lattice
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! --------------------------------------- THERML1 UPDATE FUNCTIONS ---------------------------------------
    !
    ! This module contains functions to calculate all square pulses and plaquettes on the lattice, used for 
    ! creating Torelon loops, with torelon flux pointing in a given direction.
    ! 
    ! To make use of these functions:
    ! - Nothing needs to be called in the parent code before the first call of THERML1
    !
    ! - Before the momentum sum (loop over sites) is initiated, the get_square_pulses and get_plaquettes should
    !   be called. These take as input:
    !    * The fully-blocked gauge field on a single time slice
    !    * The direction of the Torelon flux (!to match with the old code, this should be set to 1!)
    !   These functions store all square pulses and plaquettes in a memory-efficient way in the global
    !   variables squares_up, squares_down and plaquettes. After these function calls, the functions
    !   square_pulse and plaquette become available.
    !
    ! - The function square_pulse returns the square pulse in the requested direction, from the requested site.
    !   It takes as input:
    !    * The index of the site at which flux enters the square pulse
    !    * The index of the square pulse
    !    * The blocking level of the square pulse
    !   This maps to variables defined in the old code in the following way:
    !    * SQUY1 -> square_pulse(<site>, +2, <blocking_level>)
    !    * SQUZ1 -> square_pulse(<site>, +3, <blocking_level>)
    !    * SQDY1 -> square_pulse(<site>, -2, <blocking_level>)
    !    * SQDZ1 -> square_pulse(<site>, -3, <blocking_level>)
    !    * SQUY2 -> square_pulse(move(<site>, +1, <blocking_level>), +2, <blocking_level>)
    !    * SQUZ2 -> square_pulse(move(<site>, +1, <blocking_level>), +3, <blocking_level>)
    !    * SQDY2 -> square_pulse(move(<site>, +1, <blocking_level>), -2, <blocking_level>)
    !    * SQDZ2 -> square_pulse(move(<site>, +1, <blocking_level>), -3, <blocking_level>)
    !
    ! - The function plaquette returns the plaquette in the plane orthogonal to the flux tube with the
    !   requested index (detailed below), from the requested site. It takes as input:
    !    * The index of the site at which flux enters the plaquette
    !    * The index of the plaquette
    !    * The blocking level of the plaquette
    !   This function uses the convention of anti-clockwise flux around the plaquette, consistent with a
    !   'right-hand rule' in the direction of the Torelon flux. Dentoing mu and nu to be the directions in the
    !   plane perpendicular to the flux direction, they are defined cyclically:
    !   flux direction|____mu|____nu|
    !                1|     2|     3|
    !                2|     3|     1|
    !                3|     1|     2|
    !   The rotations of plaquettes are then indexed in the following way:
    !   ____(mu, nu)|____index|
    !         (+, +)|        1|
    !         (-, +)|        2|
    !         (-, -)|        3|
    !         (+, -)|        4|
    !   Inputting the negative of these indices gives the equivilant (mu, nu) plaquette with opposite
    !   circulation.
    !   This maps to variables defined in the old code in the following way:
    !    * PLQ1  -> plaquette(<site>, +1, <blocking_level>)
    !    * PLQ2  -> plaquette(<site>, +2, <blocking_level>) 
    !    * PLQ3  -> plaquette(<site>, +3, <blocking_level>) 
    !    * PLQ4  -> plaquette(<site>, +4, <blocking_level>) 
    !    * PLQ5  -> plaquette(<site>, -1, <blocking_level>)
    !    * PLQ6  -> plaquette(<site>, -2, <blocking_level>) 
    !    * PLQ7  -> plaquette(<site>, -3, <blocking_level>) 
    !    * PLQ8  -> plaquette(<site>, -4, <blocking_level>)
    !    * DPLQ1 -> plaquette(move(<site>, +1, <blocking_level>), +1, <blocking_level>)
    !    * DPLQ2 -> plaquette(move(<site>, +1, <blocking_level>), +2, <blocking_level>) 
    !    * DPLQ3 -> plaquette(move(<site>, +1, <blocking_level>), +3, <blocking_level>) 
    !    * DPLQ4 -> plaquette(move(<site>, +1, <blocking_level>), +4, <blocking_level>) 
    !    * DPLQ5 -> plaquette(move(<site>, +1, <blocking_level>), -1, <blocking_level>)
    !    * DPLQ6 -> plaquette(move(<site>, +1, <blocking_level>), -2, <blocking_level>) 
    !    * DPLQ7 -> plaquette(move(<site>, +1, <blocking_level>), -3, <blocking_level>) 
    !    * DPLQ8 -> plaquette(move(<site>, +1, <blocking_level>), -4, <blocking_level>)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(real64), allocatable :: squares_up(:,:,:,:,:), squares_down(:,:,:,:,:), plaquettes(:,:,:,:)
    integer :: squares_direction_index(3), plaquette_direction_index(3)

    contains

    ! Calculate square pulse from given site in given positive direction
    function calculate_square_pulse_up(gauge_field, blocking_level, site, flux_direction, pulse_direction)
        implicit none
        complex(real64), intent(in) :: gauge_field(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: blocking_level, site, flux_direction, pulse_direction
        complex(real64) :: calculate_square_pulse_up(NCOL, NCOL)

        integer :: site_plus_flux, site_plus_pulse

        ! Check flux and pulse directions are compatibile
        if (abs(flux_direction) == abs(pulse_direction)) error stop &
        "Flux direction cannot equal pulse direction in calculate_square_pulse_up"

        ! Check flux points in positive direction
        if (flux_direction <= 0) error stop &
        "Flux direction must be passed as positive in calculate_square_pulse_up"

        ! Check deformation points in positive direction
        if (pulse_direction <= 0) error stop &
        "Pulse direction must be passed as positive in calculate_square_pulse_up"

        ! Check both directions point in spatial directions
        if ((flux_direction > 3).or.(pulse_direction > 3)) error stop &
        "Flux and pulse directions must both be spatial directions"

        ! Get next site along flux tube
        site_plus_flux = move(site, flux_direction, blocking_level)

        ! Get site in direction of square pulse deformation
        site_plus_pulse = move(site, pulse_direction, blocking_level)

        ! Get square pulse
        calculate_square_pulse_up = matmul(matmul(&
            gauge_field(:,:,site,pulse_direction), &
            gauge_field(:,:,site_plus_pulse,flux_direction)), &
            herm(gauge_field(:,:,site_plus_flux,pulse_direction)))
    end function

    ! Calculate square pulse from given site in given negative direction
    function calculate_square_pulse_down(gauge_field, blocking_level, site, flux_direction, pulse_direction)
        implicit none
        complex(real64), intent(in) :: gauge_field(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: blocking_level, site, flux_direction, pulse_direction
        complex(real64) :: calculate_square_pulse_down(NCOL, NCOL)

        integer :: site_minus_pulse, corner_site

        ! Check flux and pulse directions are compatibile
        if (abs(flux_direction) == abs(pulse_direction)) error stop &
        "Flux direction cannot equal pulse direction in calculate_square_pulse_down"

        ! Check flux points in positive direction
        if (flux_direction <= 0) error stop &
        "Flux direction must be passed as positive in calculate_square_pulse_down"

        ! Check deformation points in positive direction
        if (pulse_direction <= 0) error stop &
        "Pulse direction must be passed as positive in calculate_square_pulse_down"

        ! Check both directions point in spatial directions
        if ((flux_direction > 3).or.(pulse_direction > 3)) error stop &
        "Flux and pulse directions must both be spatial directions"

        ! Get site in direction of square pulse deformation
        site_minus_pulse = move(site, -pulse_direction, blocking_level)

        ! Get site on opposite corner of square pulse from site
        corner_site = move(site_minus_pulse, flux_direction, blocking_level)

        ! Get square pulse
        calculate_square_pulse_down = matmul(matmul(&
            herm(gauge_field(:,:,site_minus_pulse,pulse_direction)), &
            gauge_field(:,:,site_minus_pulse,flux_direction)), &
            gauge_field(:,:,corner_site,pulse_direction))
    end function

    ! Calculate plaquette from given site in (mu, nu) = (+, +) direction
    function calculate_plaquette(gauge_field, blocking_level, site, flux_direction)
        implicit none
        complex(real64), intent(in) :: gauge_field(NCOL, NCOL, SLICE_VOLUME, 3)
        integer, intent(in) :: blocking_level, site, flux_direction
        complex(real64) :: calculate_plaquette(NCOL, NCOL)

        integer :: mu, nu, site_plus_mu, site_plus_nu

        ! Check flux direction is compatible with directions stored in plaquette_direction_index
        if (flux_direction /= plaquette_direction_index(1)) error stop &
        "flux_direction incompatible with plaquette_direction_index in calculate_plaquette"

        ! Set mu and nu from plaquette_direction_index
        mu = plaquette_direction_index(2)
        nu = plaquette_direction_index(3)

        ! Get sites in directions +mu and +nu from site
        site_plus_mu = move(site, mu, blocking_level)
        site_plus_nu = move(site, nu, blocking_level)

        ! Calculate plaquette in positive (mu, nu) direction
        calculate_plaquette = matmul(matmul(matmul(&
                              gauge_field(:, :, site, mu), &
                              gauge_field(:, :, site_plus_mu, nu)), &
                              herm(gauge_field(:, :, site_plus_nu, mu))), &
                              herm(gauge_field(:, :, site, nu)))
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
        squares_direction_index = 0
        do pulse_direction = 1, 3
            ! Skip direction parallel to flux direction
            if (pulse_direction == flux_direction) cycle

            ! Set direction index. This is a small array to map between pulse_direction and the position of the relevant entry in the 4th dimension of the squares_up/down arrays
            squares_direction_index(pulse_direction) = dir_count

            ! Iterate over all blocking levels
            do blocking_level = 1, MAX_BLOCKING_LEVEL
                ! Iterate over all sites
                do site = 1, SLICE_VOLUME
                    ! Calculate square pulses
                    squares_up(:, :, site, dir_count, blocking_level) &
                    = calculate_square_pulse_up(gauge_field_blocked(:, :, :, :, blocking_level), &
                                                blocking_level, site, flux_direction, pulse_direction)
                    squares_down(:, :, site, dir_count, blocking_level) &
                    = calculate_square_pulse_down(gauge_field_blocked(:, :, :, :, blocking_level), &
                                                  blocking_level, site, flux_direction, pulse_direction)
                enddo
            enddo

            ! Iterate dir_count
            dir_count = dir_count + 1
        enddo
    end subroutine

    ! Calculate all plaquettes in plane orthogonal to the flux direction over all blocking levels on the lattice
    subroutine get_plaquettes(gauge_field_blocked, flux_direction)
        implicit none
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(in) :: flux_direction

        integer :: site, blocking_level

        ! Check if plaquette container is allocated. If not, allocate it.
        if (.not.allocated(plaquettes)) allocate(plaquettes(NCOL, NCOL, SLICE_VOLUME, MAX_BLOCKING_LEVEL))

        ! Get orthogonal directions to the flux direction in order of cyclic permutations
        ! flux_direction|____mu|____nu|
        !              1|     2|     3|
        !              2|     3|     1|
        !              3|     1|     2|
        plaquette_direction_index(1) = flux_direction
        plaquette_direction_index(2) = mod(flux_direction, 3) + 1 ! mu
        plaquette_direction_index(3) = 6 - plaquette_direction_index(1) - plaquette_direction_index(2) ! nu

        ! Iterate over lattice and calculate all plaquettes in plane orthogonal to the flux direction
        do blocking_level = 1, MAX_BLOCKING_LEVEL
            do site = 1, SLICE_VOLUME
                ! Calculate plaquette
                plaquettes(:, :, site, blocking_level) &
                = calculate_plaquette(gauge_field_blocked(:, :, :, :, blocking_level), &
                                      blocking_level, site, flux_direction)
            enddo
        enddo
    end subroutine

    ! Return square pulse from given site, in given direction
    function square_pulse(site, pulse_direction, blocking_level)
        implicit none
        integer, intent(in) :: site, pulse_direction, blocking_level
        complex(real64) :: square_pulse(NCOL, NCOL)

        ! Check pulse direction is not parallel to flux direction
        if (squares_direction_index(abs(pulse_direction)) == 0) error stop &
        "Pulse direction must not be parallel to flux direction in square_pulse"

        ! Check pulse direction is a valid direction
        if (pulse_direction == 0 .or. abs(pulse_direction) > 3) error stop &
        "Pulse direction must point in a spatial direction in square_pulse"

        ! Return square pulse in given direction by looking up value in squares_up or squares_down
        if (pulse_direction > 0) then
            square_pulse = &
            squares_up(:, :, site, squares_direction_index(pulse_direction), blocking_level)
        else
            square_pulse = &
            squares_down(:, :, site, squares_direction_index(abs(pulse_direction)), blocking_level)
        endif
    end function

    ! Return plaquette from given site, with given index
    function plaquette(site, plaquette_index, blocking_level)
        implicit none
        integer, intent(in) :: site, plaquette_index, blocking_level
        complex(real64) :: plaquette(NCOL, NCOL)

        integer :: mu, nu

        ! Set mu and nu from plaquette_direction_index
        mu = plaquette_direction_index(2)
        nu = plaquette_direction_index(3)

        ! Get correct plaquette based on plaquette_index
        select case(plaquette_index)
        case(1)
            ! (mu, nu) = (+, +)
            plaquette = plaquettes(:, :, site, blocking_level)
        case(2)
            ! (mu, nu) = (-, +), equivalent to plaquette at (site - mu)
            plaquette = plaquettes(:, :, &
                                   move(site, -mu, blocking_level), &
                                   blocking_level)
        case(3)
            ! (mu, nu) = (-, -), equivalent to plaquette at (site - mu - nu)
            plaquette = plaquettes(:, :, &
                                   move(move(site, -mu, blocking_level), -nu, blocking_level), &
                                   blocking_level)
        case(4)
            ! (mu, nu) = (+, -), equivalent to plaquette at (site - nu)
            plaquette = plaquettes(:, :, &
                                   move(site, -nu, blocking_level), &
                                   blocking_level)
        case(-1)
            ! PLQ5 = herm(PLQ2)
            plaquette = herm(plaquettes(:, :, &
                                        move(site, -mu, blocking_level), &
                                        blocking_level))
        case(-2)
            ! PLQ6 = herm(PLQ1)
            plaquette = herm(plaquettes(:, :, site, blocking_level))
        case(-3)
            ! PLQ7 = herm(PLQ4)
            plaquette = herm(plaquettes(:, :, &
                                        move(site, -nu, blocking_level), &
                                        blocking_level))
        case(-4)
            ! PLQ8 = herm(PLQ3)
            plaquette = herm(plaquettes(:, :, &
                                        move(move(site, -mu, blocking_level), -nu, blocking_level), &
                                        blocking_level))
        case default
            ! plaquette_index is not a valid index
            error stop "Plaquette index must be +-1, +-2, +-3 or +-4"
        end select
    end function

    ! Loop 1
    !##############################################
    !           SQUARE PULSES
    !        __
    ! UP: __|  |__    DOWN: __    __
    !                         |__|
    !
    !##############################################
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
