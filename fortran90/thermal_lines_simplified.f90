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

    ! Check if there is enough space to fit a given operator in the flux direction
    logical function check_required_segments(pulse_start_offsets, operator_blocking_level)
        implicit none
        integer, intent(in) :: pulse_start_offsets(:)
        integer, intent(in) :: operator_blocking_level

        integer :: blocked_link_length
        integer :: available_segments
        integer :: required_segments

        ! Physical size of one blocked link at this operator blocking level
        blocked_link_length = 2**(operator_blocking_level - 1)
        ! Number of complete blocked links fitting in the flux direction
        available_segments = LX1 / blocked_link_length
        ! Number of blocked flux segments spanned by this operator
        required_segments = maxval(pulse_start_offsets) + 1

        ! Allow the operator only if enough blocked flux segments are available
        check_required_segments = available_segments >= required_segments
    end function check_required_segments

    ! Create the square-pulse operators given pulse directions and flux offsets
    function create_operator(in_site, out_site, number_of_square_pulses, pulse_directions, pulse_start_offsets, &
                             blocking_level, operator_flag, gauge_field_blocked) result(operator_matrix)
        implicit none
        integer, intent(in) :: in_site, number_of_square_pulses, blocking_level
        integer, intent(out) :: out_site
        integer, intent(in) :: pulse_directions(number_of_square_pulses)
        integer, intent(in) :: pulse_start_offsets(number_of_square_pulses)
        character(len=*), intent(in), optional :: operator_flag
        complex(real64), intent(in), optional :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        complex(real64) :: operator_matrix(NCOL, NCOL)

        complex(real64) :: bridge(NCOL, NCOL)
        integer :: pulse_index, pulse_site, colour_index, step, operator_blocking_level, bridge_site, bridge_next_site
        integer :: number_of_flux_steps
        integer :: pulse_spacing

        ! By default, build the operator at the requested blocking level
        operator_blocking_level = blocking_level

        ! Wave-like operators use the previous blocking level, clipped to level 1
        if (present(operator_flag)) then
            if (trim(operator_flag) == "W") then
                operator_blocking_level = blocking_level - 1
                if (operator_blocking_level == 0) operator_blocking_level = 1
            else
                error stop "Unknown operator flag in create_operator"
            endif
        endif

        ! Initialise the output operator to the identity matrix
        operator_matrix = cmplx(0.0_real64, 0.0_real64, kind=real64)
        do colour_index = 1, NCOL
            operator_matrix(colour_index, colour_index) = cmplx(1.0_real64, 0.0_real64, kind=real64)
        enddo

        ! If the operator cannot be built, leave the output as identity at the input site
        out_site = in_site

        if (.not. check_required_segments(pulse_start_offsets, operator_blocking_level)) return
        ! We can have the following cases:
        ! 2 square pulses:
        !       - adjacent square pulses
        !       - two links-spaced square pulses
        ! 4 square adjacent pulses
        select case(number_of_square_pulses)
        case(2)
            ! Decide whether the two pulses are adjacent or separated by a bridge
            pulse_spacing = pulse_start_offsets(2) - pulse_start_offsets(1)

            select case(pulse_spacing)
            case(1)
                ! Multiply adjacent square pulses in the requested order
                do pulse_index = 1, number_of_square_pulses
                    pulse_site = in_site

                    ! Move to the flux-offset where this pulse starts
                    do step = 1, pulse_start_offsets(pulse_index)
                        pulse_site = move(pulse_site, 1, operator_blocking_level)
                    enddo

                    operator_matrix = matmul( &
                        operator_matrix, &
                        square_pulse(pulse_site, pulse_directions(pulse_index), operator_blocking_level) &
                    )
                enddo
            case(3)
                ! Multiply first square pulse, bridge link, then second square pulse
                if (.not. present(gauge_field_blocked)) error stop &
                "gauge_field_blocked must be passed for square-pulse bridge square-pulse operators"

                pulse_site = in_site
                ! Move to the first square-pulse starting site
                do step = 1, pulse_start_offsets(1)
                    pulse_site = move(pulse_site, 1, operator_blocking_level)
                enddo

                operator_matrix = matmul( &
                    operator_matrix, &
                    square_pulse(pulse_site, pulse_directions(1), operator_blocking_level) &
                )

                bridge_site = in_site
                ! Bridge starts immediately after the first square pulse
                do step = 1, pulse_start_offsets(1) + 1
                    bridge_site = move(bridge_site, 1, operator_blocking_level)
                enddo

                ! At level 1, reproduce the old-code two-link bridge; otherwise use one blocked bridge
                if (blocking_level == 1) then
                    bridge_next_site = move(bridge_site, 1, operator_blocking_level)
                    bridge = matmul( &
                        gauge_field_blocked(:, :, bridge_site, 1, 1), &
                        gauge_field_blocked(:, :, bridge_next_site, 1, 1) &
                    )
                else
                    bridge = gauge_field_blocked(:, :, bridge_site, 1, blocking_level)
                endif

                operator_matrix = matmul(operator_matrix, bridge)

                pulse_site = in_site
                ! Move to the second square-pulse starting site
                do step = 1, pulse_start_offsets(2)
                    pulse_site = move(pulse_site, 1, operator_blocking_level)
                enddo

                operator_matrix = matmul( &
                    operator_matrix, &
                    square_pulse(pulse_site, pulse_directions(2), operator_blocking_level) &
                )
            case default
                error stop "Two square pulses must be adjacent or separated by one bridge in create_operator"
            end select
        case(4)
            ! Multiply four adjacent square pulses in the requested order
            do pulse_index = 1, number_of_square_pulses
                if (pulse_index < number_of_square_pulses) then
                    if (pulse_start_offsets(pulse_index + 1) - pulse_start_offsets(pulse_index) /= 1) &
                        error stop "Four square-pulse operators must have adjacent start offsets in create_operator"
                endif

                pulse_site = in_site

                ! Move to the flux-offset where this pulse starts
                do step = 1, pulse_start_offsets(pulse_index)
                    pulse_site = move(pulse_site, 1, operator_blocking_level)
                enddo

                operator_matrix = matmul( &
                    operator_matrix, &
                    square_pulse(pulse_site, pulse_directions(pulse_index), operator_blocking_level) &
                )
            enddo
        case default
            error stop "create_operator only supports 2 square pulses, square-bridge-square, or 4 square pulses"
        end select

        ! Final site is shifted by the total number of flux segments spanned
        number_of_flux_steps = maxval(pulse_start_offsets) + 1

        out_site = in_site
        do step = 1, number_of_flux_steps
            out_site = move(out_site, 1, operator_blocking_level)
        enddo
    end function


    !##############################################
    !           SQUARE PULSES
    !        __
    ! UP: __|  |__    DOWN: __    __
    !                         |__|
    !
    !##############################################
    ! Loop 1: SQUY1
    subroutine loop_1(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = square_pulse(in_site, 2, blocking_level)
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    ! Loop 2: SQUZ1
    subroutine loop_2(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = square_pulse(in_site, 3, blocking_level)
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    ! Loop 3: SQDY1
    subroutine loop_3(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = square_pulse(in_site, -2, blocking_level)
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    ! Loop 4: SQDZ1
    subroutine loop_4(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = square_pulse(in_site, -3, blocking_level)
        out_site = move(in_site, 1, blocking_level)
    end subroutine

    !##############################################
    !          UP-UP SQUARE PULSES
    !    __  __            ____
    ! __|  ||  |__  ->  __|    |__
    !
    !##############################################
    ! Loop 5: SQUY1 * SQUY2
    subroutine loop_5(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_5

    ! Loop 6: SQUZ1 * SQUZ2
    subroutine loop_6(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_6

    !##############################################
    !        DOWN-DOWN SQUARE PULSES
    ! __        __        __      __
    !   |__||__|     ->     |____|
    !
    !##############################################
    ! Loop 7: SQDY1 * SQDY2
    subroutine loop_7(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_7

    ! Loop 8: SQDZ1 * SQDZ2
    subroutine loop_8(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_8

    !##############################################
    !          UP-DOWN WAVE PULSES
    !    __
    ! __|  |   __
    !      |__|
    !
    !##############################################
    ! Loop 9: SQUY1 * SQDY2
    subroutine loop_9(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_9

    ! Loop 13: WSQUY1 * WSQDY2
    subroutine loop_13(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_13    

    ! Loop 10: SQUZ1 * SQDZ2
    subroutine loop_10(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_10

    ! Loop 14: WSQUZ1 * WSQDZ2
    subroutine loop_14(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_14

    !##############################################
    !          DOWN-UP WAVE PULSES
    !       __
    ! __   |  |__
    !   |__|
    !
    !##############################################
    ! Loop 11: SQDY1 * SQUY2
    subroutine loop_11(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_11

    ! Loop 15: WSQDY1 * WSQUY2
    subroutine loop_15(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_15

    ! Loop 12: SQDZ1 * SQUZ2
    subroutine loop_12(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_12

    ! Loop 16: WSQDZ1 * WSQUZ2
    subroutine loop_16(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_16

    !##############################################
    !          UP WAVE-UP WAVE PULSES 
    !    __    __     
    ! __|  |  |  |   __
    !      |__|  |__|
    !
    !##############################################
    ! Loop 17: WVUY1 * WVUY2
    subroutine loop_17(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -2, 2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_17

    ! Loop 18: WVUZ1 * WVUZ2
    subroutine loop_18(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -3, 3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_18

    !##############################################
    !          DOWN-WAVE DOWN-WAVE PULSES
    !       __    __ 
    ! __   |  |  |  |__
    !   |__|  |__|
    !
    !##############################################
    ! Loop 19: WVDY1 * WVDY2
    subroutine loop_19(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 2, -2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_19

    ! Loop 20: WVDZ1 * WVDZ2
    subroutine loop_20(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 3, -3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_20

    !##############################################
    !          UP WAVE-DOWN WAVE PULSES
    !    __        __              __       __
    ! __|  |      |  |__   ->   __|  |     |  |__
    !      |__||__|                  |__ __|
    !
    !##############################################
    ! Loop 21: WVUY1 * WVDY2
    subroutine loop_21(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -2, -2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_21

    ! Loop 22: WVUZ1 * WVDZ2
    subroutine loop_22(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -3, -3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_22

    !##############################################
    !          DOWN WAVE-UP WAVE SQUARE PULSES
    !       __  __                    __ __
    ! __   |  ||  |   __   ->   __   |     |   __
    !   |__|      |__|            |__|     |__|
    !
    !##############################################
    ! Loop 23: WVDY1 * WVUY2
    subroutine loop_23(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 2, 2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_23

    ! Loop 24: WVDZ1 * WVUZ2
    subroutine loop_24(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 3, 3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_24

    !##############################################
    !          UP-BRIDGE-UP PULSES
    !    __       __ 
    ! __|  |__ __|  |__
    !
    !##############################################
    ! Loop 25: WSQUY1 * D11 * WSQUY4
    subroutine loop_25(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_25

    ! Loop 26: WSQUZ1 * D11 * WSQUZ4
    subroutine loop_26(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_26

    !##############################################
    !          DOWN-BRIDGE-DOWN PULSES
    !   __    __ __    __ 
    !     |__|     |__|
    !
    !##############################################
    ! Loop 27: WSQDY1 * D11 * WSQDY4
    subroutine loop_27(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_27

    ! Loop 28: WSQDZ1 * D11 * WSQDZ4
    subroutine loop_28(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_28

    !##############################################
    !          UP-BRIDGE-DOWN PULSES
    !     __         
    !  __|  |__ __    __
    !             |__|
    !
    !##############################################
    ! Loop 29: WSQUY1 * D11 * WSQDY4
    subroutine loop_29(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_29

    ! Loop 30: WSQUZ1 * D11 * WSQDZ4
    subroutine loop_30(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_30

    !##############################################
    !          DOWN-BRIDGE-UP PULSES
    !              __        
    !  __    __ __|  |__
    !    |__|
    !
    !##############################################
    ! Loop 31: WSQDY1 * D11 * WSQUY4
    subroutine loop_31(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_31

    ! Loop 32: WSQDZ1 * D11 * WSQUZ4
    subroutine loop_32(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_32

    !##############################################
    !          UP Y PULSE + UP Z PULSE
    !        ___      
    !      _|_  |
    !   __/ |/  |__
    !
    !##############################################
    ! Loop 33: SQUY1 * SQUZ2
    subroutine loop_33(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_33

    !##############################################
    !          UP Z PULSE + DOWN Y PULSE
    !     __      
    !  __|  |__   __
    !         /__/
    !
    !##############################################
    ! Loop 34: SQUZ1 * SQDY2
    subroutine loop_34(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_34

    !##############################################
    !          DOWN Y PULSE + DOWN Z PULSE   
    !  __       __
    !   /__/|__|
    !
    !##############################################
    ! Loop 35: SQDY1 * SQDZ2
    subroutine loop_35(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_35

    !##############################################
    !          DOWN Z PULSE + UP Y PULSE   
    !           ___
    !  __    __/  /__
    !    |__|
    !
    !##############################################
    ! Loop 36: SQDZ1 * SQUY2
    subroutine loop_36(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_36

    !##############################################
    !          UP Y PULSE + DOWN Z PULSE
    !      __
    !   __|  |__    __
    !           /__/
    !
    !##############################################
    ! Loop 37: SQUY1 * SQDZ2
    subroutine loop_37(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_37

    !##############################################
    !          UP Z PULSE + UP Y PULSE
    !    __  __
    !  _|  |/ /__
    !        
    !##############################################
    ! Loop 38: SQUZ1 * SQUY2
    subroutine loop_38(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_38

    !##############################################
    !          DOWN Y PULSE + UP Z PULSE
    !          __
    !  __   __|  |__
    !   /__/ 
    !
    !##############################################
    ! Loop 39: SQDY1 * SQUZ2
    subroutine loop_39(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_39

    !##############################################
    !          DOWN Z PULSE + DOWN Y PULSE    
    !   __   ___   __
    !    |__/__|  / 
    !      /_____/    
    !
    !##############################################
    ! Loop 40: SQDZ1 * SQDY2
    subroutine loop_40(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level &
        )
    end subroutine loop_40

    !##############################################
    !          MIXED WAVE PULSES    
    !
    !##############################################
    ! Loop 41: WVUY1 * WVUZ2
    subroutine loop_41(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -2, 3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_41

    ! Loop 42: WVUZ1 * WVDY2
    subroutine loop_42(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -3, -2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_42

    ! Loop 43: WVDY1 * WVDZ2
    subroutine loop_43(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 2, -3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_43

    ! Loop 44: WVDZ1 * WVUY2
    subroutine loop_44(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 3, 2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_44

    ! Loop 45: WVUZ1 * WVUY2
    subroutine loop_45(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -3, 2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_45

    ! Loop 46: WVDY1 * WVUZ2
    subroutine loop_46(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 2, 3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_46

    ! Loop 47: WVDZ1 * WVDY2
    subroutine loop_47(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 3, -2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_47

    ! Loop 48: WVUY1 * WVDZ2
    subroutine loop_48(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -2, -3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_48

    !##############################################
    !          MIXED BRIDGE SQUARE PULSES
    !     __                __
    !  __|  |__ __ __ __ __/ /__
    !
    !##############################################
    ! Loop 49: WSQUY1 * D11 * WSQUZ4
    subroutine loop_49(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_49

    ! Loop 50: WSQUZ1 * D11 * WSQDY4
    subroutine loop_50(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_50

    ! Loop 51: WSQDY1 * D11 * WSQDZ4
    subroutine loop_51(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_51

    ! Loop 52: WSQDZ1 * D11 * WSQUY4
    subroutine loop_52(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_52

    ! Loop 53: WSQUY1 * D11 * WSQDZ4
    subroutine loop_53(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_53

    ! Loop 54: WSQUZ1 * D11 * WSQUY4
    subroutine loop_54(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_54

    ! Loop 55: WSQDY1 * D11 * WSQUZ4
    subroutine loop_55(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_55

    ! Loop 56: WSQDZ1 * D11 * WSQDY4
    subroutine loop_56(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2]
        pulse_start_offsets = [0, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W", &
            gauge_field_blocked=gauge_field_blocked &
        )
    end subroutine loop_56

    !##############################################
    !                 T-T4 OPERATORS
    !
    !##############################################
    ! Loop 57: WSQUY1 * WSQUZ2 * WSQUY3 * WSQUZ4
    subroutine loop_57(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_57

    ! Loop 58: WSQUZ1 * WSQDY2 * WSQUZ3 * WSQDY4
    subroutine loop_58(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_58

    ! Loop 59: WSQDY1 * WSQDZ2 * WSQDY3 * WSQDZ4
    subroutine loop_59(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_59

    ! Loop 60: WSQDZ1 * WSQUY2 * WSQDZ3 * WSQUY4
    subroutine loop_60(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_60

    ! Loop 61: WSQDZ1 * WSQDY2 * WSQDZ3 * WSQDY4
    subroutine loop_61(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_61

    ! Loop 62: WSQUY1 * WSQDZ2 * WSQUY3 * WSQDZ4
    subroutine loop_62(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_62

    ! Loop 63: WSQUZ1 * WSQUY2 * WSQUZ3 * WSQUY4
    subroutine loop_63(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_63

    ! Loop 64: WSQDY1 * WSQUZ2 * WSQDY3 * WSQUZ4
    subroutine loop_64(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_64

    !##############################################
    !                 T-T5 OPERATORS
    !
    !##############################################
    ! Loop 65: WSQUY1 * WSQUY2 * WSQUY3 * WSQUY4
    subroutine loop_65(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 2, 2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_65

    ! Loop 66: WSQUZ1 * WSQUZ2 * WSQUZ3 * WSQUZ4
    subroutine loop_66(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 3, 3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_66

    ! Loop 67: WSQDY1 * WSQDY2 * WSQDY3 * WSQDY4
    subroutine loop_67(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -2, -2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_67

    ! Loop 68: WSQDZ1 * WSQDZ2 * WSQDZ3 * WSQDZ4
    subroutine loop_68(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -3, -3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_68

    !##############################################
    !                 T-T6 OPERATORS
    !
    !##############################################
    ! Loop 69: DUY1 * DDY2
    subroutine loop_69(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 2, -2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_69

    ! Loop 70: DUZ1 * DDZ2
    subroutine loop_70(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 3, -3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_70

    ! Loop 71: DDY1 * DUY2
    subroutine loop_71(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -2, 2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_71

    ! Loop 72: DDZ1 * DUZ2
    subroutine loop_72(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -3, 3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_72

    !##############################################
    !                 T-T7 OPERATORS
    !
    !##############################################
    ! Loop 73: DUY1 * DUZ2
    subroutine loop_73(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 2, 3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_73

    ! Loop 74: DUZ1 * DDY2
    subroutine loop_74(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 3, -2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_74

    ! Loop 75: DDY1 * DDZ2
    subroutine loop_75(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -2, -3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_75

    ! Loop 76: DDZ1 * DUY2
    subroutine loop_76(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -3, 2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_76

    ! Loop 77: DDZ1 * DDY2
    subroutine loop_77(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -3, -2, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_77

    ! Loop 78: DUY1 * DDZ2
    subroutine loop_78(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 2, -3, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_78

    ! Loop 79: DUZ1 * DUY2
    subroutine loop_79(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 3, 2, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_79

    ! Loop 80: DDY1 * DUZ2
    subroutine loop_80(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -2, 3, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_80

    !##############################################
    !                 T-T8 OPERATORS
    !
    !##############################################
    ! Loop 81: WSQUY1 * WSQUZ2 * WSQUZ3 * WSQUY4
    subroutine loop_81(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_81

    ! Loop 82: WSQUZ1 * WSQDY2 * WSQDY3 * WSQUZ4
    subroutine loop_82(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_82

    ! Loop 83: WSQDY1 * WSQDZ2 * WSQDZ3 * WSQDY4
    subroutine loop_83(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_83

    ! Loop 84: WSQDZ1 * WSQUY2 * WSQUY3 * WSQDZ4
    subroutine loop_84(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_84

    ! Loop 85: WSQDY1 * WSQUZ2 * WSQUZ3 * WSQDY4
    subroutine loop_85(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_85

    ! Loop 86: WSQUZ1 * WSQUY2 * WSQUY3 * WSQUZ4
    subroutine loop_86(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_86

    ! Loop 87: WSQUY1 * WSQDZ2 * WSQDZ3 * WSQUY4
    subroutine loop_87(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_87

    ! Loop 88: WSQDZ1 * WSQDY2 * WSQDY3 * WSQDZ4
    subroutine loop_88(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_88

    !##############################################
    !                 T-T9 OPERATORS
    !
    !##############################################
    ! Loop 89: WSQUY1 * WSQUZ2 * WSQUZ3 * WSQDY4
    subroutine loop_89(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_89

    ! Loop 90: WSQUZ1 * WSQDY2 * WSQDY3 * WSQDZ4
    subroutine loop_90(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_90

    ! Loop 91: WSQDY1 * WSQDZ2 * WSQDZ3 * WSQUY4
    subroutine loop_91(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_91

    ! Loop 92: WSQDZ1 * WSQUY2 * WSQUY3 * WSQUZ4
    subroutine loop_92(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_92

    ! Loop 93: WSQUY1 * WSQDZ2 * WSQDZ3 * WSQDY4
    subroutine loop_93(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_93

    ! Loop 94: WSQUZ1 * WSQUY2 * WSQUY3 * WSQDZ4
    subroutine loop_94(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_94

    ! Loop 95: WSQDY1 * WSQUZ2 * WSQUZ3 * WSQUY4
    subroutine loop_95(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_95

    ! Loop 96: WSQDZ1 * WSQDY2 * WSQDY3 * WSQUZ4
    subroutine loop_96(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_96

    !##############################################
    !                 T-T10 OPERATORS
    !
    !##############################################
    ! Loop 97: WSQUY1 * WSQUZ2 * WSQDZ3 * WSQUY4
    subroutine loop_97(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_97

    ! Loop 98: WSQUZ1 * WSQDY2 * WSQUY3 * WSQUZ4
    subroutine loop_98(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_98

    ! Loop 99: WSQDY1 * WSQDZ2 * WSQUZ3 * WSQDY4
    subroutine loop_99(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_99

    ! Loop 100: WSQDZ1 * WSQUY2 * WSQDY3 * WSQDZ4
    subroutine loop_100(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_100

    ! Loop 101: WSQDY1 * WSQUZ2 * WSQDZ3 * WSQDY4
    subroutine loop_101(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_101

    ! Loop 102: WSQUZ1 * WSQUY2 * WSQDY3 * WSQUZ4
    subroutine loop_102(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_102

    ! Loop 103: WSQUY1 * WSQDZ2 * WSQUZ3 * WSQUY4
    subroutine loop_103(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_103

    ! Loop 104: WSQDZ1 * WSQDY2 * WSQUY3 * WSQDZ4
    subroutine loop_104(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_104

    !##############################################
    !                 T-T11 OPERATORS
    !
    !##############################################
    ! Loop 105: WSQUY1 * WSQUZ2 * WSQDZ3 * WSQDY4
    subroutine loop_105(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_105

    ! Loop 106: WSQUZ1 * WSQDY2 * WSQUY3 * WSQDZ4
    subroutine loop_106(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_106

    ! Loop 107: WSQDY1 * WSQDZ2 * WSQUZ3 * WSQUY4
    subroutine loop_107(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_107

    ! Loop 108: WSQDZ1 * WSQUY2 * WSQDY3 * WSQUZ4
    subroutine loop_108(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_108

    ! Loop 109: WSQDY1 * WSQUZ2 * WSQDZ3 * WSQUY4
    subroutine loop_109(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_109

    ! Loop 110: WSQUZ1 * WSQUY2 * WSQDY3 * WSQDZ4
    subroutine loop_110(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_110

    ! Loop 111: WSQUY1 * WSQDZ2 * WSQUZ3 * WSQDY4
    subroutine loop_111(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_111

    ! Loop 112: WSQDZ1 * WSQDY2 * WSQUY3 * WSQUZ4
    subroutine loop_112(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_112

    !##############################################
    !                 T-T12 OPERATORS
    !
    !##############################################
    ! Loop 113: WSQUY1 * WSQUZ2 * WSQUY3 * WSQDZ4
    subroutine loop_113(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_113

    ! Loop 114: WSQUZ1 * WSQDY2 * WSQUZ3 * WSQUY4
    subroutine loop_114(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_114

    ! Loop 115: WSQDY1 * WSQDZ2 * WSQDY3 * WSQUZ4
    subroutine loop_115(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_115

    ! Loop 116: WSQDZ1 * WSQUY2 * WSQDZ3 * WSQDY4
    subroutine loop_116(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_116

    ! Loop 117: WSQUZ1 * WSQDY2 * WSQDZ3 * WSQDY4
    subroutine loop_117(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_117

    ! Loop 118: WSQDY1 * WSQDZ2 * WSQUY3 * WSQDZ4
    subroutine loop_118(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_118

    ! Loop 119: WSQDZ1 * WSQUY2 * WSQUZ3 * WSQUY4
    subroutine loop_119(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_119

    ! Loop 120: WSQUY1 * WSQUZ2 * WSQDY3 * WSQUZ4
    subroutine loop_120(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_120

    ! Loop 121: WSQDY1 * WSQUZ2 * WSQDY3 * WSQDZ4
    subroutine loop_121(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_121

    ! Loop 122: WSQUZ1 * WSQUY2 * WSQUZ3 * WSQDY4
    subroutine loop_122(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_122

    ! Loop 123: WSQUY1 * WSQDZ2 * WSQUY3 * WSQUZ4
    subroutine loop_123(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_123

    ! Loop 124: WSQDZ1 * WSQDY2 * WSQDZ3 * WSQUY4
    subroutine loop_124(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_124

    ! Loop 125: WSQUZ1 * WSQUY2 * WSQDZ3 * WSQUY4
    subroutine loop_125(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_125

    ! Loop 126: WSQUY1 * WSQDZ2 * WSQDY3 * WSQDZ4
    subroutine loop_126(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_126

    ! Loop 127: WSQDZ1 * WSQDY2 * WSQUZ3 * WSQDY4
    subroutine loop_127(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_127

    ! Loop 128: WSQDY1 * WSQUZ2 * WSQUY3 * WSQUZ4
    subroutine loop_128(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_128

    !##############################################
    !                 T-T13 OPERATORS
    !
    !##############################################
    ! Loop 129: WSQUY1 * WSQUZ2 * WSQDY3 * WSQDZ4
    subroutine loop_129(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3, -2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_129

    ! Loop 130: WSQUZ1 * WSQDY2 * WSQDZ3 * WSQUY4
    subroutine loop_130(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2, -3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_130

    ! Loop 131: WSQDY1 * WSQDZ2 * WSQUY3 * WSQUZ4
    subroutine loop_131(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3, 2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_131

    ! Loop 132: WSQDZ1 * WSQUY2 * WSQUZ3 * WSQDY4
    subroutine loop_132(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2, 3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_132

    ! Loop 133: WSQUZ1 * WSQUY2 * WSQDZ3 * WSQDY4
    subroutine loop_133(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2, -3, -2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_133

    ! Loop 134: WSQDY1 * WSQUZ2 * WSQUY3 * WSQDZ4
    subroutine loop_134(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3, 2, -3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_134

    ! Loop 135: WSQDZ1 * WSQDY2 * WSQUZ3 * WSQUY4
    subroutine loop_135(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2, 3, 2]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_135

    ! Loop 136: WSQUY1 * WSQDZ2 * WSQDY3 * WSQUZ4
    subroutine loop_136(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 4
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3, -2, 3]
        pulse_start_offsets = [0, 1, 2, 3]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_136

    !##############################################
    !                 T-T14 OPERATORS
    !
    !##############################################
    ! Loop 137: WSQUY1 * WSQUZ2
    subroutine loop_137(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, 3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_137

    ! Loop 138: WSQUZ1 * WSQDY2
    subroutine loop_138(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, -2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_138

    ! Loop 139: WSQDY1 * WSQDZ2
    subroutine loop_139(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, -3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_139

    ! Loop 140: WSQDZ1 * WSQUY2
    subroutine loop_140(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, 2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_140

    ! Loop 141: WSQUY1 * WSQDZ2
    subroutine loop_141(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [2, -3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_141

    ! Loop 142: WSQUZ1 * WSQUY2
    subroutine loop_142(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [3, 2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_142

    ! Loop 143: WSQDY1 * WSQUZ2
    subroutine loop_143(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-2, 3]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_143

    ! Loop 144: WSQDZ1 * WSQDY2
    subroutine loop_144(in_site, out_site, A11, blocking_level)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)

        pulse_directions = [-3, -2]
        pulse_start_offsets = [0, 1]

        A11 = create_operator( in_site, out_site, &
            number_of_square_pulses, pulse_directions, pulse_start_offsets, &
            blocking_level, "W" &
        )
    end subroutine loop_144

    ! Loop 145: Simple Polyakov line
    subroutine loop_145(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        integer, parameter :: number_of_square_pulses = 2
        integer :: pulse_directions(number_of_square_pulses)
        integer :: pulse_start_offsets(number_of_square_pulses)


        A11 = gauge_field_blocked(:, :, 1, in_site, blocking_level)
    end subroutine loop_145

end module
