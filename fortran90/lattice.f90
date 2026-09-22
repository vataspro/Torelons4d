module lattice
    use parameters
    implicit none

    ! Containers for lattice pointers by 1 lattice step and multiple lattice steps in the case of blocked links
    integer, allocatable :: lattice_pointers_up(:,:), lattice_pointers_down(:,:), &
    lattice_pointers_blocked_up(:,:,:), lattice_pointers_blocked_down(:,:,:)

    contains

    ! Set up lattice pointers and blocked lattice pointers
    subroutine setup_lattice()
        implicit none

        integer :: n, neighbour, LX12, x, y, z, t, blocking_level, direction

        ! Allocate lattice points with pattern (site, direction)
        allocate(lattice_pointers_up(LATTICE_VOLUME, 4), lattice_pointers_down(LATTICE_VOLUME, 4))

        ! Allocate blocked lattice pointers with pattern (site, direction, blocking level)
        allocate(lattice_pointers_blocked_up(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1), &
        lattice_pointers_blocked_down(SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL+1))

        ! Get lattice size of x-y slice
        LX12 = LX1 * LX2

        ! Calculate indices of lattice points in all directions
        n = 0
        do t = 1, LX4
            do z = 1, LX3
                do y = 1, LX2
                    do x = 1, LX1
                        n = n + 1

                        ! Moving in x-direction
                        ! Move up
                        neighbour = n + 1
                        if ((x+1) > LX1) neighbour = neighbour - LX1
                        lattice_pointers_up(n, 1) = neighbour
                        ! Move down
                        neighbour = n - 1
                        if ((x-1) < 1) neighbour = neighbour + LX1
                        lattice_pointers_down(n, 1) = neighbour

                        ! Moving in y-direction
                        ! Move up
                        neighbour = n + LX1
                        if ((y+1) > LX2) neighbour = neighbour - LX12
                        lattice_pointers_up(n, 2) = neighbour
                        ! Move down
                        neighbour = n - LX1
                        if ((y-1) < 1) neighbour = neighbour + LX12
                        lattice_pointers_down(n, 2) = neighbour

                        ! Moving in z-direction
                        ! Move up
                        neighbour = n + LX12
                        if ((z+1) > LX3) neighbour = neighbour - SLICE_VOLUME
                        lattice_pointers_up(n, 3) = neighbour
                        ! Move down
                        neighbour = n - LX12
                        if ((z-1) < 1) neighbour = neighbour + SLICE_VOLUME
                        lattice_pointers_down(n, 3) = neighbour

                        ! Moving in t-direction
                        ! Move up
                        neighbour = n + SLICE_VOLUME
                        if ((t+1 > LX4)) neighbour = neighbour - LATTICE_VOLUME
                        lattice_pointers_up(n, 4) = neighbour
                        ! Move down
                        neighbour = n - SLICE_VOLUME
                        if ((t-1) < 1) neighbour = neighbour + LATTICE_VOLUME
                        lattice_pointers_down(n, 4) = neighbour
                    enddo
                enddo
            enddo
        enddo

        !! Calculate indices of lattice points in all directions on blocked lattice
        ! Set blocking level 1 as equal to the original lattice
        lattice_pointers_blocked_up(:,:,1) = lattice_pointers_up(1:SLICE_VOLUME, 1:3)
        lattice_pointers_blocked_down(:,:,1) = lattice_pointers_down(1:SLICE_VOLUME, 1:3)

        ! For each subsequent blocking level, move twice the distance as the previous blocking level in all directions from all sites
        do blocking_level = 2, MAX_BLOCKING_LEVEL+1
            do direction = 1, 3
                do n = 1, SLICE_VOLUME
                    ! Move up by 2 (blocked) units
                    lattice_pointers_blocked_up(n, direction, blocking_level) &
                    = lattice_pointers_blocked_up(&
                    lattice_pointers_blocked_up(n, direction, blocking_level-1), &
                    direction, blocking_level-1)

                    ! Move down by 2 (blocked) units
                    lattice_pointers_blocked_down(n, direction, blocking_level) &
                    = lattice_pointers_blocked_down(&
                    lattice_pointers_blocked_down(n, direction, blocking_level-1), &
                    direction, blocking_level-1)
                enddo
            enddo
        enddo
    end subroutine setup_lattice

    ! Move around lattice by either blocked or unblocked amounts
    integer function move(site, direction, blocking_level)
        implicit none
        integer, intent(in) :: site, direction
        integer, intent(in), optional :: blocking_level

        ! Check parameters are compatible and fall within the range of the lattice, then return requested lattice site
        if (present(blocking_level)) then
            ! Site must belong to a time slice and direction must point in a spatial direction
            if (site > SLICE_VOLUME) then
                write(*, '(a, i0)') "site must be on a spatial slice in move function with blocking_level argument present, therefore must be less than or equal to ", SLICE_VOLUME
            endif
            if (abs(direction) > 3 .or. direction == 0) then
                write(*, '(a)') "direction must be spatial in move function with blocking_level argument present, therefore must be less than or equal to 3"
                stop
            endif
            if (blocking_level > MAX_BLOCKING_LEVEL+1) then
                write(*, '(a, i0)') "blocking_level must be less than ", MAX_BLOCKING_LEVEL+1
                stop
            endif

            ! Return blocked neighbour
            if (direction > 0) then
                move = lattice_pointers_blocked_up(site, direction, blocking_level)
            else
                move = lattice_pointers_blocked_down(site, abs(direction), blocking_level)
            endif
        else
            ! Site must belong to the lattice and direction must point in 4d
            if (site > LATTICE_VOLUME) then
                write(*, '(a, i0)') "site must be on the lattice in move function without blocking_level argument present, therefore must be less than or equal to ", LATTICE_VOLUME
            endif
            if (abs(direction) > 4 .or. direction  == 0) then
                write(*, '(a)') "direction must be spatial or temporal in move function without blocking_level argument present, therefore must be less than or equal to 4"
                stop
            endif

            ! Return unblocked neighbour
            if (direction > 0) then
                move = lattice_pointers_up(site, direction)
            else
                move = lattice_pointers_down(site, abs(direction))
            endif
        endif
    end function move
end module lattice