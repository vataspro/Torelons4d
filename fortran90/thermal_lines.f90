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

    complex(real64), allocatable :: squares_up(:,:,:,:,:), squares_down(:,:,:,:,:), plaquettes(:,:,:,:,:)
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

        integer :: site, blocking_level, mu, nu, site_minus_mu, site_minus_nu, diagonal_site
        complex(real64) :: similarity_matrix(NCOL, NCOL)

        ! Check if plaquette container is allocated. If not, allocate it.
        if (.not.allocated(plaquettes)) allocate(plaquettes(NCOL, NCOL, SLICE_VOLUME, 4, MAX_BLOCKING_LEVEL))

        ! Get orthogonal directions to the flux direction in order of cyclic permutations
        ! flux_direction|____mu|____nu|
        !              1|     2|     3|
        !              2|     3|     1|
        !              3|     1|     2|
        plaquette_direction_index(1) = flux_direction
        plaquette_direction_index(2) = mod(flux_direction, 3) + 1 ! mu
        mu = plaquette_direction_index(2)
        plaquette_direction_index(3) = 6 - plaquette_direction_index(1) - plaquette_direction_index(2) ! nu
        nu = plaquette_direction_index(3)

        ! Iterate over lattice and calculate all plaquettes in plane orthogonal to the flux direction
        do blocking_level = 1, MAX_BLOCKING_LEVEL
            ! Calculate plaquette in (mu, nu) = (+, +) direction
            do site = 1, SLICE_VOLUME
                plaquettes(:, :, site, 1, blocking_level) &
                = calculate_plaquette(gauge_field_blocked(:, :, :, :, blocking_level), &
                                      blocking_level, site, flux_direction)
            enddo

            ! Calculate plaquettes in other directions by applying similarity transformation to (+,+) plaquettes
            do site = 1, SLICE_VOLUME
                ! Get relevant sites
                site_minus_mu = move(site, -mu, blocking_level)
                site_minus_nu = move(site, -nu, blocking_level)
                diagonal_site = move(site_minus_mu, -nu, blocking_level)

                ! (mu, nu) = (-, +)
                similarity_matrix = gauge_field_blocked(:, :, site_minus_mu, mu, blocking_level)
                plaquettes(:, :, site, 2, blocking_level) = matmul(matmul(&
                    herm(similarity_matrix), &
                    plaquettes(:, :, site_minus_mu, 1, blocking_level)), &
                    similarity_matrix)

                ! (mu, nu) = (-, -)
                similarity_matrix = matmul(&
                    gauge_field_blocked(:, :, diagonal_site, nu, blocking_level), &
                    gauge_field_blocked(:, :, site_minus_mu, mu, blocking_level))
                plaquettes(:, :, site, 3, blocking_level) = matmul(matmul(&
                    herm(similarity_matrix), &
                    plaquettes(:, :, diagonal_site, 1, blocking_level)), &
                    similarity_matrix)

                ! (mu, nu) = (+, -)
                similarity_matrix = gauge_field_blocked(:, :, site_minus_nu, nu, blocking_level)
                plaquettes(:, :, site, 4, blocking_level) = matmul(matmul(&
                    herm(similarity_matrix), &
                    plaquettes(:, :, site_minus_nu, 1, blocking_level)), &
                    similarity_matrix)
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
            plaquette = plaquettes(:, :, site, 1, blocking_level)
        case(2)
            ! (mu, nu) = (-, +), similar to plaquette at (site - mu)
            plaquette = plaquettes(:, :, site, 2, blocking_level)
        case(3)
            ! (mu, nu) = (-, -), similar to plaquette at (site - mu - nu)
            plaquette = plaquettes(:, :, site, 3, blocking_level)
        case(4)
            ! (mu, nu) = (+, -), similar to plaquette at (site - nu)
            plaquette = plaquettes(:, :, site, 4, blocking_level)
        case(-1)
            ! PLQ5 = herm(PLQ2)
            plaquette = herm(plaquettes(:, :, site, 2, blocking_level))
        case(-2)
            ! PLQ6 = herm(PLQ1)
            plaquette = herm(plaquettes(:, :, site, 1, blocking_level))
        case(-3)
            ! PLQ7 = herm(PLQ4)
            plaquette = herm(plaquettes(:, :, site, 4, blocking_level))
        case(-4)
            ! PLQ8 = herm(PLQ3)
            plaquette = herm(plaquettes(:, :, site, 3, blocking_level))
        case default
            ! plaquette_index is not a valid index
            error stop "Plaquette index must be +-1, +-2, +-3 or +-4"
        end select
    end function

    ! Construct deformation composed of square pulses and gauge links
    ! Inputs:
    ! - in_site: the site at which the deformation starts
    ! - gauge_field_blocked: the blocked gauge field configuration on the given time slice
    ! - blocking_level: the blocking level of the polyakov loop to be attached to the deformation
    ! - square_pulse_positions: the relative positions of sqaure pulses along the flux tube in units of the blocking level (if narrow_operator = false) or at one lower blocking level (if narrow operator = true)
    ! - square_pulse_directions: the directions of the square pulses corresponding to those at square_pulse_positions
    ! - narrow_operator: a boolean flag indicating if the deformation is to be calculated at one lower blocking level
    ! - check_length: a boolean flag indicating whether to check that the operator is too large for the lattice. If true, the function returns the identity matrix and out_site=in_site.
    ! Outputs:
    ! - advance_site: a boolean flag indicating whether out_site should be returned as in_site or the end of the loop. If true, out_site is returned as the end point of the operator. If false, out_site is returned as in_site.
    ! - out_site: the site are which the deformation ends
    ! - A11: the matrix representing the deformation
    subroutine loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
        square_pulse_positions, square_pulse_directions, narrow_operator, check_length, advance_site, &
        deformation)
        integer, intent(in) :: in_site, blocking_level, square_pulse_positions(:), square_pulse_directions(:)
        integer, intent(out) :: out_site
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        logical, intent(in) :: narrow_operator, check_length, advance_site
        complex(real64), intent(out) :: deformation(NCOL, NCOL)

        integer :: deformation_blocking_level, site, array_index, position, i, n_pulses, length

        ! Check square_pulse_positions and square_pulse_directions are compatible
        n_pulses = size(square_pulse_directions)
        if (size(square_pulse_positions) /= n_pulses .or. square_pulse_positions(1) /= 0) &
        error stop "square_pulse_positions and square_pulse_directions incompatible"

        ! Check that elements of square_pulse_positions always increase, so that this corresponds to a continuous operator
        do array_index = 1, n_pulses - 1
            if (square_pulse_positions(array_index) >= square_pulse_positions(array_index+1)) &
            error stop "square_pulse_positions does not correspond to a continuous operator"
        enddo

        ! Check narrow operator can fit into Polyakov loop at input blocking level
        if (narrow_operator .and. mod(square_pulse_positions(n_pulses), 2) /= 1) &
        error stop "Narrow operator does not fit into Polyakov loop"

        ! Assign blocking level of the deformation
        if (narrow_operator) then
            if (blocking_level == 1) then
                deformation_blocking_level = 1
            else
                deformation_blocking_level = blocking_level - 1
            endif
        else
            deformation_blocking_level = blocking_level
        endif

        ! Initialise deformation to identity matrix
        deformation = cmplx(0.0,0.0,kind=real64)
        do i = 1, NCOL
            deformation(i, i) = cmplx(1.0,0.0,kind=real64)
        enddo

        ! Check that the operator fits inside the lattice by finding the position in unblocked lattice units of the end of the operator
        if (check_length) then
            length = (square_pulse_positions(n_pulses) + 1) * 2**(deformation_blocking_level-1)
            if (length > LX1) then ! Operator does not fit along the flux tube
                ! Output identity matrix and return input site
                out_site = in_site
                return
            endif
        endif

        ! Build operator
        site = in_site
        array_index = 1
        if (narrow_operator) then
            ! Advance along the flux tube, adding square pulses at blocking level b-1, and filling in with gauge links at blocking level b. If b=1, fill in with two gauge links.
            position = 0
            do while (position <= square_pulse_positions(n_pulses))
                if (array_index > n_pulses) error stop "array_index exceeded number of pulses"
                if (square_pulse_positions(array_index) == position) then
                    ! In this case, add a square pulse at this position at the appropriate blocking level
                    deformation = matmul(deformation, &
                    square_pulse(site, square_pulse_directions(array_index), deformation_blocking_level))

                    ! Advance by half a step
                    site = move(site, 1, deformation_blocking_level)
                    array_index = array_index + 1
                    position = position + 1
                else
                    ! In this case, fill the gap between square pulses in with gauge links
                    if (blocking_level == 1) then
                        do i = 1, 2
                            deformation = matmul(deformation, &
                            gauge_field_blocked(:, :, site, 1, 1))

                            ! Advance by a full step (after iteration)
                            site = move(site, 1, 1)
                        enddo
                    else
                        deformation = matmul(deformation, &
                        gauge_field_blocked(:, :, site, 1, blocking_level))

                        ! Advance by a full step
                        site = move(site, 1, blocking_level)
                    endif
                    position = position + 2
                endif
            enddo
        else
            ! Advance along the flux tube, adding square pulses and filling in with gauge links, until we get to the last square pulse
            do position = 0, square_pulse_positions(n_pulses)
                if (array_index > n_pulses) then
                    error stop "array_index exceeded number of pulses"
                elseif (square_pulse_positions(array_index) == position) then
                    ! In this case, add a square pulse at this position
                    deformation = matmul(deformation, &
                    square_pulse(site, square_pulse_directions(array_index), blocking_level))
                    array_index = array_index + 1
                else
                    ! In this case, fill the gap between pulses in with a blocked gauge link
                    deformation = matmul(deformation, &
                    gauge_field_blocked(:, :, site, 1, blocking_level))
                endif
                ! Move along flux tube to next position
                site = move(site, 1, blocking_level)
            enddo
        endif

        ! Set output position
        if (advance_site) then
            out_site = site
        else
            out_site = in_site
        endif
    end subroutine

    ! Construct deformation composed of plaquettes and gauge links
    ! Inputs:
    ! - in_site: the site at which the deformation starts
    ! - gauge_field_blocked: the blocked gauge field configuration on the given time slice
    ! - blocking_level: the blocking level of the polyakov loop to be attached to the deformation
    ! - plaquette_positions: the relative positions of the plaquettes along the flux tube in units of the blocking level
    ! - plaquette_orientations: the orientations of the plaquettes corresponding to those at plaquette_positions
    ! Outputs:
    ! - check_length: a boolean flag indicating whether to check that the operator is too large for the lattice. If true, the function returns the identity matrix and out_site=in_site.
    ! - advance_site: a boolean flag indicating whether out_site should be returned as in_site or the end of the loop. If true, out_site is returned as the end point of the operator. If false, out_site is returned as in_site.
    ! Outputs:
    ! - out_site: the site are which the deformation ends
    ! - A11: the matrix representing the deformation
    subroutine loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
        plaquette_positions, plaquette_orientations, check_length, advance_site, deformation)
        integer, intent(in) :: in_site, blocking_level, plaquette_positions(:), plaquette_orientations(:)
        integer, intent(out) :: out_site
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        logical, intent(in) :: check_length, advance_site
        complex(real64), intent(out) :: deformation(NCOL, NCOL)

        integer :: site, array_index, position, i, n_plaquettes, length

        ! Check plaquette_positions and plaquette_orientations are compatible
        n_plaquettes = size(plaquette_orientations)
        if (size(plaquette_positions) /= n_plaquettes .or. plaquette_positions(1) /= 0) &
        error stop "plaquette_positions and plaquette_orientations incompatible"

        ! Check that elements of plaquette_positions always increase or stay the same, so that this corresponds to a continuous operator
        do array_index = 1, n_plaquettes - 1
            if (plaquette_positions(array_index) > plaquette_positions(array_index+1)) &
            error stop "plaquette_positions does not correspond to a continuous operator"
        enddo

        ! Initialise deformation to identity matrix
        deformation = cmplx(0.0,0.0,kind=real64)
        do i = 1, NCOL
            deformation(i, i) = cmplx(1.0,0.0,kind=real64)
        enddo

        ! Check that the operator fits inside the lattice by finding the position in unblocked lattice units of the end of the operator
        if (check_length) then
            length = plaquette_positions(n_plaquettes) * 2**(blocking_level-1)
            if (length > LX1) then ! Operator does not fit along the flux tube
                ! Output identity matrix and return input site
                out_site = in_site
                return
            endif
        endif

        ! Build operator
        ! Multiply all plaquettes at the starting site into deformation in the given order
        array_index = 1
        do
            if (array_index > n_plaquettes) then
                error stop "plaquette_index exceeded number of plaquettes"
            elseif (plaquette_positions(array_index) /= 0) then
                exit
            else
                deformation = matmul(deformation, &
                plaquette(in_site, plaquette_orientations(array_index), blocking_level))
                array_index = array_index + 1
            endif
        enddo
        ! Multiply other plaquettes at further sites along the flux tube into deformation in the given order, separated by gauge links
        site = in_site
        do position = 1, plaquette_positions(n_plaquettes)
            ! Advance to this new position by adding a gauge link
            deformation = matmul(deformation, &
            gauge_field_blocked(:, :, site, 1, blocking_level))
            site = move(site, 1, blocking_level)

            ! Multiply all plaquettes at the current site into deformation in the given order
            do
                if (array_index > n_plaquettes) exit
                if (plaquette_positions(array_index) /= position) exit

                deformation = matmul(deformation, &
                plaquette(site, plaquette_orientations(array_index), blocking_level))

                array_index = array_index + 1
            enddo
        enddo

        ! Set output position
        if (advance_site) then
            out_site = site
        else
            out_site = in_site
        endif
    end subroutine

    ! Construct deformation composed of square pulses, plaquettes and gauge links
    ! Inputs:
    ! - in_site: the site at which the deformation starts
    ! - gauge_field_blocked: the blocked gauge field configuration on the given time slice
    ! - blocking_level: the blocking level of the polyakov loop to be attached to the deformation
    ! - square_pulse_positions: the relative positions of sqaure pulses along the flux tube in units of the blocking level
    ! - square_pulse_directions: the directions of the square pulses corresponding to those at square_pulse_positions
    ! - plaquette_positions: the relative positions of the plaquettes along the flux tube in units of the blocking level
    ! - plaquette_orientations: the orientations of the plaquettes corresponding to those at plaquette_positions
    ! Outputs:
    ! - check_length: a boolean flag indicating whether to check that the operator is too large for the lattice. If true, the function returns the identity matrix and out_site=in_site.
    ! - advance_site: a boolean flag indicating whether out_site should be returned as in_site or the end of the loop. If true, out_site is returned as the end point of the operator. If false, out_site is returned as in_site.
    ! Outputs:
    ! - out_site: the site are which the deformation ends
    ! - A11: the matrix representing the deformation
    subroutine loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
        square_pulse_positions, square_pulse_directions, plaquette_positions, plaquette_orientations, &
        check_length, advance_site, deformation)
        integer, intent(in) :: in_site, blocking_level, square_pulse_positions(:), square_pulse_directions(:), &
        plaquette_positions(:), plaquette_orientations(:)
        integer, intent(out) :: out_site
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        logical, intent(in) :: check_length, advance_site
        complex(real64), intent(out) :: deformation(NCOL, NCOL)

        integer :: site, pulse_index, plaquette_index, position, i, n_pulses, n_plaquettes, length

        ! Check square_pulse_positions and square_pulse_directions are compatible
        n_pulses = size(square_pulse_directions)
        if (size(square_pulse_positions) /= n_pulses) &
        error stop "square_pulse_positions and square_pulse_directions incompatible"

        ! Check plaquette_positions and plaquette_orientations are compatible
        n_plaquettes = size(plaquette_orientations)
        if (size(plaquette_positions) /= n_plaquettes) &
        error stop "plaquette_positions and plaquette_orientations incompatible"

        ! Check that elements of square_pulse_positions always increase, so that this corresponds to a continuous operator
        do pulse_index = 1, n_pulses - 1
            if (square_pulse_positions(pulse_index) >= square_pulse_positions(pulse_index+1)) &
            error stop "square_pulse_positions does not correspond to a continuous operator"
        enddo

        ! Check that elements of plaquette_positions always increase or stay the same, so that this corresponds to a continuous operator
        do plaquette_index = 1, n_plaquettes - 1
            if (plaquette_positions(plaquette_index) > plaquette_positions(plaquette_index+1)) &
            error stop "plaquette_positions does not correspond to a continuous operator"
        enddo

        ! Check that operator starts at in_site
        if (square_pulse_positions(1) /= 0 .and. plaquette_positions(1) /= 0) &
        error stop "Operator does not begin at input site"

        ! Initialise deformation to identity matrix
        deformation = cmplx(0.0,0.0,kind=real64)
        do i = 1, NCOL
            deformation(i, i) = cmplx(1.0,0.0,kind=real64)
        enddo

        ! Check that the operator fits inside the lattice by finding the position in unblocked lattice units of the end of the operator
        if (check_length) then
            if (square_pulse_positions(n_pulses) >= plaquette_positions(n_plaquettes)) then
                length = (square_pulse_positions(n_pulses) + 1) * 2**(blocking_level-1)
            else
                length = plaquette_positions(n_plaquettes) * 2**(blocking_level-1)
            endif
            if (length > LX1) then ! Operator does not fit along the flux tube
                ! Output identity matrix and return input site
                out_site = in_site
                return
            endif
        endif

        ! Build operator
        ! If the last deformation is a plaquette, we need to ensure we do not add an extra gauge link to the end of the operator
        if (square_pulse_positions(n_pulses) >= plaquette_positions(n_plaquettes)) then
            site = in_site
            pulse_index = 1
            plaquette_index = 1
            do position = 0, square_pulse_positions(n_pulses)
                ! Multiply all plaquettes at current site into deformation in the given order
                do
                    if (plaquette_index > n_plaquettes) exit
                    if (plaquette_positions(plaquette_index) /= position) exit

                    deformation = matmul(deformation, &
                    plaquette(site, plaquette_orientations(plaquette_index), blocking_level))

                    plaquette_index = plaquette_index + 1
                enddo

                ! Multiply by square pulse or gauge link and advance to next site
                if (pulse_index <= n_pulses) then
                    if (square_pulse_positions(pulse_index) == position) then
                        deformation = matmul(deformation, &
                        square_pulse(site, square_pulse_directions(pulse_index), blocking_level))
                        pulse_index = pulse_index + 1
                    else
                        deformation = matmul(deformation, &
                        gauge_field_blocked(:, :, site, 1, blocking_level))
                    endif
                else ! All square pulses have been used
                    deformation = matmul(deformation, &
                    gauge_field_blocked(:, :, site, 1, blocking_level))
                endif
                site = move(site, 1, blocking_level)
            enddo
        else
            ! Multiply all plaquettes at the starting site into deformation in the given order
            plaquette_index = 1
            do
                if (plaquette_index > n_plaquettes) exit
                if (plaquette_positions(plaquette_index) /= 0) exit

                deformation = matmul(deformation, &
                plaquette(in_site, plaquette_orientations(plaquette_index), blocking_level))

                plaquette_index = plaquette_index + 1
            enddo

            ! Multiple other plaquettes at further sites along the flux tube into deformation in the given order, separated square pulses or by gauge links
            site = in_site
            pulse_index = 1
            do position = 1, plaquette_positions(n_plaquettes)
                ! Advance to this new position by adding a square pulse or a gauge link
                if (pulse_index <= n_pulses) then
                    if (square_pulse_positions(pulse_index) == position-1) then
                        deformation = matmul(deformation, &
                        square_pulse(site, square_pulse_directions(pulse_index), blocking_level))
                        pulse_index = pulse_index + 1
                    else
                        deformation = matmul(deformation, &
                        gauge_field_blocked(:, :, site, 1, blocking_level))
                    endif
                else ! All square pulses have been used
                    deformation = matmul(deformation, &
                    gauge_field_blocked(:, :, site, 1, blocking_level))
                endif
                site = move(site, 1, blocking_level)

                ! Multiply all plaquettes at the current site into deformation in the given order
                do
                    if (plaquette_index > n_plaquettes) exit
                    if (plaquette_positions(plaquette_index) /= position) exit

                    deformation = matmul(deformation, &
                    plaquette(site, plaquette_orientations(plaquette_index), blocking_level))
                    
                    plaquette_index = plaquette_index + 1
                enddo
            enddo
        endif

        ! Set output position
        if (advance_site) then
            out_site = site
        else
            out_site = in_site
        endif
    end subroutine

    ! Return the nxn identity matrix 
    function get_I(n) result(eye_matrix)
        implicit none
        integer(int32), intent(in) :: n
        integer(int32) ::  ic, jc
        complex(real64) :: eye_matrix(NCOL, NCOL)

        do ic=1,NCOL
            do jc=1,NCOL
                if (ic == jc) then
                    eye_matrix(ic, jc) = 1
                else
                    eye_matrix(ic, jc) = 0
                end if
            end do
        end do     

        return
    end function get_I

    ! Loop 1
    subroutine loop_1(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0], [2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_1

    ! Loop 2
    subroutine loop_2(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0], [3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_2

    ! Loop 3
    subroutine loop_3(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0], [-2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_3

    ! Loop 4
    subroutine loop_4(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0], [-3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_4

    ! Loop 5
    subroutine loop_5(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [2, 2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_5

    ! Loop 6
    subroutine loop_6(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [3, 3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_6

    ! Loop 7
    subroutine loop_7(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-2, -2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_7

    ! Loop 8
    subroutine loop_8(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-3, -3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_8

    ! Loop 9
    subroutine loop_9(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [2, -2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_9

    ! Loop 10
    subroutine loop_10(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [3, -3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_10

    ! Loop 11
    subroutine loop_11(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-2, 2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_11

    ! Loop 12
    subroutine loop_12(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-3, 3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_12

    ! Loop 13
    subroutine loop_13(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_13

    ! Loop 14
    subroutine loop_14(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_14

    ! Loop 15
    subroutine loop_15(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_15

    ! Loop 16
    subroutine loop_16(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_16

    ! Loop 17
    subroutine loop_17(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -2, 2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_17

    ! Loop 18
    subroutine loop_18(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -3, 3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_18

    ! Loop 19
    subroutine loop_19(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 2, -2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_19

    ! Loop 20
    subroutine loop_20(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 3, -3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_20

    ! Loop 21
    subroutine loop_21(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -2, -2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_21

    ! Loop 22
    subroutine loop_22(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -3, -3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_22

    ! Loop 23
    subroutine loop_23(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 2, 2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_23

    ! Loop 24
    subroutine loop_24(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 3, 3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_24

    ! Loop 25
    subroutine loop_25(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_25

    ! Loop 26
    subroutine loop_26(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_26

    ! Loop 27
    subroutine loop_27(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_27

    ! Loop 28
    subroutine loop_28(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_28

    ! Loop 29
    subroutine loop_29(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_29

    ! Loop 30
    subroutine loop_30(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_30

    ! Loop 31
    subroutine loop_31(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_31

    ! Loop 32
    subroutine loop_32(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_32

    ! Loop 33
    subroutine loop_33(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [2, 3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_33

    ! Loop 34
    subroutine loop_34(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [3, -2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_34

    ! Loop 35
    subroutine loop_35(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-2, -3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_35

    ! Loop 36
    subroutine loop_36(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-3, 2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_36

    ! Loop 37
    subroutine loop_37(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [2, -3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_37

    ! Loop 38
    subroutine loop_38(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [3, 2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_38

    ! Loop 39
    subroutine loop_39(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-2, 3], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_39

    ! Loop 40
    subroutine loop_40(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-3, -2], &
                                  narrow_operator = .false., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_40

    ! Loop 41
    subroutine loop_41(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -2, 3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_41

    ! Loop 42
    subroutine loop_42(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -3, -2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_42

    ! Loop 43
    subroutine loop_43(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 2, -3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_43

    ! Loop 44
    subroutine loop_44(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 3, 2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_44

    ! Loop 45
    subroutine loop_45(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -3, 2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_45

    ! Loop 46
    subroutine loop_46(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 2, 3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_46

    ! Loop 47
    subroutine loop_47(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 3, -2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_47

    ! Loop 48
    subroutine loop_48(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -2, -3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_48

    ! Loop 49
    subroutine loop_49(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_49

    ! Loop 50
    subroutine loop_50(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_50

    ! Loop 51
    subroutine loop_51(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_51

    ! Loop 52
    subroutine loop_52(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_52

    ! Loop 53
    subroutine loop_53(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_53

    ! Loop 54
    subroutine loop_54(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_54

    ! Loop 55
    subroutine loop_55(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_55

    ! Loop 56
    subroutine loop_56(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 3], [-3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_56

    ! Loop 57
    subroutine loop_57(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_57

    ! Loop 58
    subroutine loop_58(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_58

    ! Loop 59
    subroutine loop_59(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_59

    ! Loop 60
    subroutine loop_60(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_60

    ! Loop 61
    subroutine loop_61(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_61

    ! Loop 62
    subroutine loop_62(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_62

    ! Loop 63
    subroutine loop_63(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_63

    ! Loop 64
    subroutine loop_64(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_64

    ! Loop 65
    subroutine loop_65(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 2, 2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_65

    ! Loop 66
    subroutine loop_66(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 3, 3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_66

    ! Loop 67
    subroutine loop_67(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -2, -2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_67

    ! Loop 68
    subroutine loop_68(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -3, -3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_68

    ! Loop 69
    subroutine loop_69(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 2, -2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_69

    ! Loop 70
    subroutine loop_70(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 3, -3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_70

    ! Loop 71
    subroutine loop_71(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -2, 2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_71

    ! Loop 72
    subroutine loop_72(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -3, 3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_72

    ! Loop 73
    subroutine loop_73(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 2, 3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_73

    ! Loop 74
    subroutine loop_74(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 3, -2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_74

    ! Loop 75
    subroutine loop_75(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -2, -3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_75

    ! Loop 76
    subroutine loop_76(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -3, 2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_76

    ! Loop 77
    subroutine loop_77(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -3, -2, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_77

    ! Loop 78
    subroutine loop_78(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 2, -3, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_78

    ! Loop 79
    subroutine loop_79(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 3, 2, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_79

    ! Loop 80
    subroutine loop_80(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -2, 3, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_80

    ! Loop 81
    subroutine loop_81(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_81

    ! Loop 82
    subroutine loop_82(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_82

    ! Loop 83
    subroutine loop_83(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_83

    ! Loop 84
    subroutine loop_84(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_84

    ! Loop 85
    subroutine loop_85(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_85

    ! Loop 86
    subroutine loop_86(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_86

    ! Loop 87
    subroutine loop_87(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_87

    ! Loop 88
    subroutine loop_88(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_88

    ! Loop 89
    subroutine loop_89(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_89

    ! Loop 90
    subroutine loop_90(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_90

    ! Loop 91
    subroutine loop_91(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_91

    ! Loop 92
    subroutine loop_92(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_92

    ! Loop 93
    subroutine loop_93(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_93

    ! Loop 94
    subroutine loop_94(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_94

    ! Loop 95
    subroutine loop_95(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_95

    ! Loop 96
    subroutine loop_96(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_96

    ! Loop 97
    subroutine loop_97(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_97

    ! Loop 98
    subroutine loop_98(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_98

    ! Loop 99
    subroutine loop_99(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_99

    ! Loop 100
    subroutine loop_100(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_100

    ! Loop 101
    subroutine loop_101(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_101

    ! Loop 102
    subroutine loop_102(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_102

    ! Loop 103
    subroutine loop_103(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_103

    ! Loop 104
    subroutine loop_104(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_104

    ! Loop 105
    subroutine loop_105(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_105

    ! Loop 106
    subroutine loop_106(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_106

    ! Loop 107
    subroutine loop_107(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_107

    ! Loop 108
    subroutine loop_108(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_108

    ! Loop 109
    subroutine loop_109(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_109

    ! Loop 110
    subroutine loop_110(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_110

    ! Loop 111
    subroutine loop_111(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_111

    ! Loop 112
    subroutine loop_112(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_112

    ! Loop 113
    subroutine loop_113(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_113

    ! Loop 114
    subroutine loop_114(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_114

    ! Loop 115
    subroutine loop_115(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_115

    ! Loop 116
    subroutine loop_116(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_116

    ! Loop 117
    subroutine loop_117(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_117

    ! Loop 118
    subroutine loop_118(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_118

    ! Loop 119
    subroutine loop_119(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_119

    ! Loop 120
    subroutine loop_120(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_120

    ! Loop 121
    subroutine loop_121(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_121

    ! Loop 122
    subroutine loop_122(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_122

    ! Loop 123
    subroutine loop_123(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_123

    ! Loop 124
    subroutine loop_124(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_124

    ! Loop 125
    subroutine loop_125(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_125

    ! Loop 126
    subroutine loop_126(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_126

    ! Loop 127
    subroutine loop_127(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_127

    ! Loop 128
    subroutine loop_128(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_128

    ! Loop 129
    subroutine loop_129(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, 3, -2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_129

    ! Loop 130
    subroutine loop_130(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, -2, -3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_130

    ! Loop 131
    subroutine loop_131(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, -3, 2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_131

    ! Loop 132
    subroutine loop_132(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, 2, 3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_132

    ! Loop 133
    subroutine loop_133(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [3, 2, -3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_133

    ! Loop 134
    subroutine loop_134(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-2, 3, 2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_134

    ! Loop 135
    subroutine loop_135(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [-3, -2, 3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_135

    ! Loop 136
    subroutine loop_136(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1, 2, 3], [2, -3, -2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .false., &
                                  deformation = A11)
    end subroutine loop_136

    ! Loop 137
    subroutine loop_137(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_137

    ! Loop 138
    subroutine loop_138(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_138

    ! Loop 139
    subroutine loop_139(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_139

    ! Loop 140
    subroutine loop_140(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_140

    ! Loop 141
    subroutine loop_141(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [2, -3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_141

    ! Loop 142
    subroutine loop_142(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [3, 2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_142

    ! Loop 143
    subroutine loop_143(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-2, 3], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_143

    ! Loop 144
    subroutine loop_144(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_square_pulse(in_site, out_site, gauge_field_blocked, blocking_level, &
                                  [0, 1], [-3, -2], &
                                  narrow_operator = .true., &
                                  check_length = .true., &
                                  advance_site = .true., &
                                  deformation = A11)
    end subroutine loop_144

    ! Loop 197
    subroutine loop_197(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-3, -2], &
                                             [1], [4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_197

    ! Loop 198
    subroutine loop_198(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-3, 2], &
                                             [1], [-2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_198

    ! Loop 199
    subroutine loop_199(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [2, 3], &
                                             [1], [-1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_199

    ! Loop 200
    subroutine loop_200(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [3, -2], &
                                             [1], [-4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_200

    ! Loop 201
    subroutine loop_201(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-2, -3], &
                                             [1], [-3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_201

    ! Loop 202
    subroutine loop_202(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-2, -3], &
                                             [1], [-1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_202

    ! Loop 203
    subroutine loop_203(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [3, -2], &
                                             [1], [-2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_203

    ! Loop 204
    subroutine loop_204(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [2, 3], &
                                             [1], [-3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_204

    ! Loop 205
    subroutine loop_205(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-3, 2], &
                                             [1], [-4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_205

    ! Loop 206
    subroutine loop_206(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-3, -2], &
                                             [1], [2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_206

    ! Loop 207
    subroutine loop_207(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-2, 3], &
                                             [1], [1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_207

    ! Loop 208
    subroutine loop_208(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [3, 2], &
                                             [1], [4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_208

    ! Loop 209
    subroutine loop_209(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [2, -3], &
                                             [1], [3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_209

    ! Loop 210
    subroutine loop_210(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1], [1, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_210

    ! Loop 211
    subroutine loop_211(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1], [2, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_211

    ! Loop 212
    subroutine loop_212(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1], [3, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_212

    ! Loop 213
    subroutine loop_213(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1], [4, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_213

    ! Loop 214
    subroutine loop_214(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0], [-1, -2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_214

    ! Loop 215
    subroutine loop_215(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0], [-4, -1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_215

    ! Loop 216
    subroutine loop_216(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0], [-3, -4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_216

    ! Loop 217
    subroutine loop_217(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0], [-2, -3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_217

    ! Loop 218
    subroutine loop_218(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1], [-1, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_218

    ! Loop 219
    subroutine loop_219(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1], [-2, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_219

    ! Loop 220
    subroutine loop_220(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1], [-3, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_220

    ! Loop 221
    subroutine loop_221(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1], [-4, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_221

    ! Loop 222
    subroutine loop_222(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0], [1, 2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_222

    ! Loop 223
    subroutine loop_223(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0], [4, 1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_223

    ! Loop 224
    subroutine loop_224(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0], [3, 4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_224

    ! Loop 225
    subroutine loop_225(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0], [2, 3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_225

    ! Loop 226
    subroutine loop_226(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1, 1], [1, 2, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_226

    ! Loop 227
    subroutine loop_227(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1, 1], [2, 3, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_227

    ! Loop 228
    subroutine loop_228(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1, 1], [3, 4, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_228

    ! Loop 229
    subroutine loop_229(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1, 1], [4, 1, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_229

    ! Loop 230
    subroutine loop_230(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0, 0], [-4, -1, -2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_230

    ! Loop 231
    subroutine loop_231(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0, 0], [-3, -4, -1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_231

    ! Loop 232
    subroutine loop_232(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0, 0], [-2, -3, -4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_232

    ! Loop 233
    subroutine loop_233(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0, 0], [-1, -2, -3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_233

    ! Loop 234
    subroutine loop_234(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1, 1], [-1, -2, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_234

    ! Loop 235
    subroutine loop_235(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1, 1], [-2, -3, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_235

    ! Loop 236
    subroutine loop_236(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1, 1], [-3, -4, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_236

    ! Loop 237
    subroutine loop_237(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1, 1], [-4, -1, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_237

    ! Loop 238
    subroutine loop_238(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0, 0], [4, 1, 2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_238

    ! Loop 239
    subroutine loop_239(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0, 0], [3, 4, 1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_239

    ! Loop 240
    subroutine loop_240(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0, 0], [2, 3, 4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_240

    ! Loop 241
    subroutine loop_241(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0, 0], [1, 2, 3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_241

    ! Loop 242
    subroutine loop_242(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0], [2, -3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_242

    ! Loop 243
    subroutine loop_243(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0], [3, -2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_243

    ! Loop 244
    subroutine loop_244(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0], [4, -1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_244

    ! Loop 245
    subroutine loop_245(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0], [1, -4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_245

    ! Loop 246
    subroutine loop_246(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1], [4, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_246

    ! Loop 247
    subroutine loop_247(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1], [1, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_247

    ! Loop 248
    subroutine loop_248(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1], [2, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_248

    ! Loop 249
    subroutine loop_249(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1], [3, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_249

    ! Loop 250
    subroutine loop_250(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0], [-2, 3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_250

    ! Loop 251
    subroutine loop_251(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0], [-3, 2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_251

    ! Loop 252
    subroutine loop_252(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0], [-4, 1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_252

    ! Loop 253
    subroutine loop_253(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0], [-1, 4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_253

    ! Loop 254
    subroutine loop_254(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1], [-4, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_254

    ! Loop 255
    subroutine loop_255(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1], [-1, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_255

    ! Loop 256
    subroutine loop_256(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1], [-2, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_256

    ! Loop 257
    subroutine loop_257(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1], [-3, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_257

    ! Loop 258
    subroutine loop_258(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0], [2, 4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_258

    ! Loop 259
    subroutine loop_259(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0], [3, 1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_259

    ! Loop 260
    subroutine loop_260(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0], [4, 2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_260

    ! Loop 261
    subroutine loop_261(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0], [1, 3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_261

    ! Loop 262
    subroutine loop_262(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1], [-3, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_262

    ! Loop 263
    subroutine loop_263(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1], [-2, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_263

    ! Loop 264
    subroutine loop_264(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1], [-1, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_264

    ! Loop 265
    subroutine loop_265(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1], [-4, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_265

    ! Loop 266
    subroutine loop_266(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0], [-2, -4], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_266

    ! Loop 267
    subroutine loop_267(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0], [-3, -1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_267

    ! Loop 268
    subroutine loop_268(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0], [-4, -2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_268

    ! Loop 269
    subroutine loop_269(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0], [-1, -3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
    end subroutine loop_269

    ! Loop 270
    subroutine loop_270(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [1, 1], [3, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_270

    ! Loop 271
    subroutine loop_271(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [1, 1], [2, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_271

    ! Loop 272
    subroutine loop_272(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [1, 1], [1, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_272

    ! Loop 273
    subroutine loop_273(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [1, 1], [4, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_273

    ! Loop 274
    subroutine loop_274(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0, 1], [2, 3, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_274

    ! Loop 275
    subroutine loop_275(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0, 1], [3, 4, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_275

    ! Loop 276
    subroutine loop_276(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0, 1], [4, 1, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_276

    ! Loop 277
    subroutine loop_277(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0, 1], [1, 2, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_277

    ! Loop 278
    subroutine loop_278(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1, 1], [4, -4, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_278

    ! Loop 279
    subroutine loop_279(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1, 1], [1, -3, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_279

    ! Loop 280
    subroutine loop_280(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1, 1], [2, -2, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_280

    ! Loop 281
    subroutine loop_281(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1, 1], [3, -1, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_281

    ! Loop 282
    subroutine loop_282(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 0, 1], [-2, -3, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_282

    ! Loop 283
    subroutine loop_283(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 0, 1], [-3, -4, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_283

    ! Loop 284
    subroutine loop_284(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 0, 1], [-4, -1, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_284

    ! Loop 285
    subroutine loop_285(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 0, 1], [-1, -2, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_285

    ! Loop 286
    subroutine loop_286(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1, 1], [-4, 4, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_286

    ! Loop 287
    subroutine loop_287(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1, 1], [-1, 3, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_287

    ! Loop 288
    subroutine loop_288(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1, 1], [-2, 2, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_288

    ! Loop 289
    subroutine loop_289(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1, 1], [-3, 1, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_289

    ! Loop 290
    subroutine loop_290(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [2, 3, 4, 1, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_290

    ! Loop 291
    subroutine loop_291(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [3, 4, 1, 2, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_291

    ! Loop 292
    subroutine loop_292(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [4, 1, 2, 3, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_292

    ! Loop 293
    subroutine loop_293(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [1, 2, 3, 4, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_293

    ! Loop 294
    subroutine loop_294(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [-1, -2, -3, -4, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_294

    ! Loop 295
    subroutine loop_295(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [-4, -1, -2, -3, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_295

    ! Loop 296
    subroutine loop_296(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [-3, -4, -1, -2, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_296

    ! Loop 297
    subroutine loop_297(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [-2, -3, -4, -1, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_297

    ! Loop 298
    subroutine loop_298(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [-2, -3, -4, -1, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_298

    ! Loop 299
    subroutine loop_299(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [-3, -4, -1, -2, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_299

    ! Loop 300
    subroutine loop_300(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [-4, -1, -2, -3, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_300

    ! Loop 301
    subroutine loop_301(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 0, 1, 1], [-1, -2, -3, -4, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_301

    ! Loop 302
    subroutine loop_302(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [1, 2, 3, 4, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_302

    ! Loop 303
    subroutine loop_303(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [4, 1, 2, 3, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_303

    ! Loop 304
    subroutine loop_304(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [3, 4, 1, 2, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_304

    ! Loop 305
    subroutine loop_305(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1, 1, 1], [2, 3, 4, 1, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_305

    ! Loop 306
    subroutine loop_306(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1], [2, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_306

    ! Loop 307
    subroutine loop_307(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1], [3, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_307

    ! Loop 308
    subroutine loop_308(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1], [4, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_308

    ! Loop 309
    subroutine loop_309(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1], [1, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_309

    ! Loop 310
    subroutine loop_310(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1], [4, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_310

    ! Loop 311
    subroutine loop_311(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1], [1, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_311

    ! Loop 312
    subroutine loop_312(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1], [2, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_312

    ! Loop 313
    subroutine loop_313(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1], [3, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_313

    ! Loop 314
    subroutine loop_314(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1], [-2, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_314

    ! Loop 315
    subroutine loop_315(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1], [-3, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_315

    ! Loop 316
    subroutine loop_316(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1], [-4, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_316

    ! Loop 317
    subroutine loop_317(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1], [-1, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_317

    ! Loop 318
    subroutine loop_318(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1], [-4, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_318

    ! Loop 319
    subroutine loop_319(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1], [-1, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_319

    ! Loop 320
    subroutine loop_320(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1], [-2, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_320

    ! Loop 321
    subroutine loop_321(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1], [-3, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_321

    ! Loop 322
    subroutine loop_322(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [2, -3, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_322

    ! Loop 323
    subroutine loop_323(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [3, -2, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_323

    ! Loop 324
    subroutine loop_324(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [4, -1, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_324

    ! Loop 325
    subroutine loop_325(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [1, -4, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_325

    ! Loop 326
    subroutine loop_326(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [1, 4, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_326

    ! Loop 327
    subroutine loop_327(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [2, 1, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_327

    ! Loop 328
    subroutine loop_328(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [3, 2, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_328

    ! Loop 329
    subroutine loop_329(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [4, 3, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_329

    ! Loop 330
    subroutine loop_330(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [-2, 3, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_330

    ! Loop 331
    subroutine loop_331(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [-3, 2, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_331

    ! Loop 332
    subroutine loop_332(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [-4, 1, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_332

    ! Loop 333
    subroutine loop_333(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 0, 1], [-1, 4, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_333

    ! Loop 334
    subroutine loop_334(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [-1, -4, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_334

    ! Loop 335
    subroutine loop_335(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [-2, -1, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_335

    ! Loop 336
    subroutine loop_336(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [-3, -2, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_336

    ! Loop 337
    subroutine loop_337(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1, 1], [-4, -3, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_337

    ! Loop 145: normal Polyakov link
    subroutine loop_145(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        A11 = gauge_field_blocked(:, :, in_site, 1, blocking_level)
        out_site = in_site
    end subroutine loop_145

    ! Loop 146: plaquette in (+,+) direction x 1 polyakov link
    subroutine loop_146(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(plaquette(in_site, 1, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))


        out_site = in_site

    end subroutine loop_146


    ! Loop 147: plaquette in (-,+) direction x 1 polyakov link
    subroutine loop_147(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        out_site = in_site

        A11 = matmul(plaquette(in_site, 2, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))

    end subroutine loop_147


    ! Loop 148: plaquette in (+,-) direction x 1 polyakov link
    subroutine loop_148(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        out_site = in_site

        A11 = matmul(plaquette(out_site, 3, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))

    end subroutine loop_148


    ! Loop 149: plaquette in (-,-) direction x 1 polyakov link
    subroutine loop_149(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        out_site = in_site

        A11 = matmul(plaquette(out_site, 4, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))

    end subroutine loop_149


    ! Loop 150: hermitian plaquette in (+,+) direction x 1 polyakov link
    subroutine loop_150(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(plaquette(in_site, -1, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))

        out_site = in_site

    end subroutine loop_150

    ! Loop 151: hermitian plaquette in (-,+) direction x 1 polyakov link
    subroutine loop_151(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(plaquette(in_site, -2, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))

        out_site = in_site

    end subroutine loop_151

    ! Loop 152: hermitian plaquette in (+,-) direction x 1 polyakov link
    subroutine loop_152(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(plaquette(in_site, -3, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))

        out_site = in_site

    end subroutine loop_152

    ! Loop 153: hermitian plaquette in (-,-) direction x 1 polyakov link
    subroutine loop_153(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)

        A11 = matmul(plaquette(in_site, -4, blocking_level), &
                    gauge_field_blocked(:, :, in_site, 1, blocking_level))

        out_site = in_site

    end subroutine loop_153


    !##############################################
    !        DOUBLE  PLAQUETTES
    !
    !             |\   |\
    !            x|_\__| \
    !             \ |  \ |
    !              \|   \|
    !
    !##############################################
    ! Loop 154
    subroutine loop_154(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [1, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_154

    ! Loop 155
    subroutine loop_155(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [2, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_155

    ! Loop 156
    subroutine loop_156(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [3, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_156

    ! Loop 157
    subroutine loop_157(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [4, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_157

    ! Loop 158
    subroutine loop_158(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-1, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_158

    ! Loop 159
    subroutine loop_159(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-2, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_159

    ! Loop 160
    subroutine loop_160(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-3, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_160

    ! Loop 161
    subroutine loop_161(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-4, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_161

    ! Loop 162
    subroutine loop_162(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [1, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_162

    ! Loop 163
    subroutine loop_163(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [2, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_163

    ! Loop 164
    subroutine loop_164(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [3, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_164

    ! Loop 165
    subroutine loop_165(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [4, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_165

    ! Loop 166
    subroutine loop_166(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-1, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_166

    ! Loop 167
    subroutine loop_167(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-2, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_167

    ! Loop 168
    subroutine loop_168(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-3, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_168

    ! Loop 169
    subroutine loop_169(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-4, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_169

    ! Loop 170
    subroutine loop_170(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1], [1, -2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_170

    ! Loop 171
    subroutine loop_171(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1], [2, -1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_171

    ! Loop 172
    subroutine loop_172(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1], [3, -4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_172

    ! Loop 173
    subroutine loop_173(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1], [4, -3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_173

    ! Loop 174
    subroutine loop_174(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [3], &
                                             [0, 1], [-1, 2], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_174

    ! Loop 175
    subroutine loop_175(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [2], &
                                             [0, 1], [-2, 1], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_175

    ! Loop 176
    subroutine loop_176(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-3], &
                                             [0, 1], [-3, 4], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_176

    ! Loop 177
    subroutine loop_177(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0], [-2], &
                                             [0, 1], [-4, 3], &
                                             check_length = .true., &
                                             advance_site = .true., &
                                             deformation = A11)
    end subroutine loop_177

    ! Loop 178
    subroutine loop_178(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [1, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_178

    ! Loop 179
    subroutine loop_179(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [2, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_179

    ! Loop 180
    subroutine loop_180(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [3, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_180

    ! Loop 181
    subroutine loop_181(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [4, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_181

    ! Loop 182
    subroutine loop_182(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-1, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_182

    ! Loop 183
    subroutine loop_183(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-2, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_183

    ! Loop 184
    subroutine loop_184(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-3, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_184

    ! Loop 185
    subroutine loop_185(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-4, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_185

    ! Loop 186
    subroutine loop_186(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [1, 3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_186

    ! Loop 187
    subroutine loop_187(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [2, 4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_187

    ! Loop 188
    subroutine loop_188(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [3, 1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_188

    ! Loop 189
    subroutine loop_189(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [4, 2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_189

    ! Loop 190
    subroutine loop_190(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-1, -3], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_190

    ! Loop 191
    subroutine loop_191(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-2, -4], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_191

    ! Loop 192
    subroutine loop_192(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-3, -1], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_192

    ! Loop 193
    subroutine loop_193(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)

        call loop_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                             [0, 1], [-4, -2], &
                             check_length = .true., &
                             advance_site = .true., &
                             deformation = A11)
    end subroutine loop_193

    ! Loop 194
    subroutine loop_194(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [2, -3], &
                                             [1], [1], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_194

    ! Loop 195
    subroutine loop_195(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [3, 2], &
                                             [1], [2], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_195

    ! Loop 196
    subroutine loop_196(in_site, out_site, A11, blocking_level, gauge_field_blocked)
        implicit none
        integer, intent(in) :: in_site, blocking_level
        integer, intent(out) :: out_site
        complex(real64), intent(inout) :: A11(NCOL, NCOL)
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL)
        call loop_square_pulse_plaquette(in_site, out_site, gauge_field_blocked, blocking_level, &
                                             [0, 1], [-2, 3], &
                                             [1], [3], &
                                             check_length = .true., &
                                             advance_site = .false., &
                                             deformation = A11)
        out_site = move(in_site, 1, blocking_level)
    end subroutine loop_196

end module
