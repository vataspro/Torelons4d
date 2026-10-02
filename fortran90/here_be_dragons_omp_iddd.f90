module here_be_dragons
    use torelon_parameters
    use lattice
    use thermal_lines
    implicit none

    contains

    subroutine THERML1(gauge_field_blocked, time_slice, blocking_level, lines, momentum_lines)
        implicit none

        ! Inputs
        complex(real64), intent(in) :: gauge_field_blocked(NCOL, NCOL, SLICE_VOLUME, 3, MAX_BLOCKING_LEVEL) ! UB11
        integer, intent(in) :: time_slice, blocking_level ! N4, IBLL

        ! Outputs
        complex(real64), intent(inout) :: lines(LX4, MAX_BLOCKING_LEVEL, 235) ! ALINE
        complex(real64), intent(inout) :: momentum_lines(LX4, MAX_BLOCKING_LEVEL, 2, 2:235) ! ALINEMOM

        ! Dummy variables
        complex(real64) :: gauge_field(NCOL, NCOL, SLICE_VOLUME, 3) ! UC11
        integer(int32) :: lcnt(MAX_BLOCKING_LEVEL), lb(MAX_BLOCKING_LEVEL), ix(3), ls(3)
        complex(real64) :: A11(NCOL, NCOL), B11(NCOL, NCOL), C11(NCOL, NCOL), &
        D11(NCOL, NCOL), E11(NCOL, NCOL), UINT11(NCOL, NCOL), REM11(NCOL, NCOL)

        ! Loop variables
        integer :: idd, ib, idg, iloop, i, ic, nc, ig, mu, ik

        ! Site indices
        integer :: site, site_ku, site_ju, site_iu, mn, ml, m1, m2, m3, m4, m5

        ! Direction indices
        integer :: iu, ju, ku

        ! Blocking level variables
        integer :: id, lrest, ids, idsm1, irem, li, ico, idsw

        ! Operator variables
        integer :: iddd, ieee
        ! Rotation/reflection index (ieee) of each operator iddd. It is a fixed function of iddd,
        ! so the iddd iterations are independent. (0: csumn, which has no index)
        integer, parameter :: ieee_of(337) = [ &
            1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, &
            1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, &
            1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, &
            5, 6, 7, 8, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, &
            1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, &
            5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, &
            9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, &
            5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, &
            8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, &
            4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, &
            8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, &
            12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, &
            16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, &
            4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, &
            8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, &
            12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, &
            16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 &
        ]

        ! Momentum variables
        real(real64) :: phase, adiv1, adiv2, adiv3, adivn

        ! Square pulses
        complex(real64) :: squy1(NCOL, NCOL), sqdy1(NCOL, NCOL), squz1(NCOL, NCOL), sqdz1(NCOL, NCOL), &
                            squy2(NCOL, NCOL), sqdy2(NCOL, NCOL), squz2(NCOL, NCOL), sqdz2(NCOL, NCOL), &
                            squdy1(NCOL, NCOL), sqduy1(NCOL, NCOL), squdz1(NCOL, NCOL), sqduz1(NCOL, NCOL)

        ! Wave-like square pulses
        complex(real64) ::  wsquy1(NCOL, NCOL), wsquy2(NCOL, NCOL), wsquy3(NCOL, NCOL), wsquy4(NCOL, NCOL), &
                            wsquz1(NCOL, NCOL), wsquz2(NCOL, NCOL), wsquz3(NCOL, NCOL), wsquz4(NCOL, NCOL), &
                            wsqdy1(NCOL, NCOL), wsqdy2(NCOL, NCOL), wsqdy3(NCOL, NCOL), wsqdy4(NCOL, NCOL), &
                            wsqdz1(NCOL, NCOL), wsqdz2(NCOL, NCOL), wsqdz3(NCOL, NCOL), wsqdz4(NCOL, NCOL)

        ! TT-operator components
        complex(real64) :: duy1(NCOL, NCOL), duz1(NCOL, NCOL), ddy1(NCOL, NCOL), ddz1(NCOL, NCOL), &
                            duy2(NCOL, NCOL), duz2(NCOL, NCOL), ddy2(NCOL, NCOL), ddz2(NCOL, NCOL)

        ! Plaquette operator components
        complex(real64) :: pq1(NCOL, NCOL), pq2(NCOL, NCOL), pq3(NCOL, NCOL), pq4(NCOL, NCOL), &
                            pq5(NCOL, NCOL), pq6(NCOL, NCOL), pq7(NCOL, NCOL), pq8(NCOL, NCOL)

        ! WV components
        complex(real64) :: wvuy1(NCOL, NCOL), wvuy2(NCOL, NCOL), wvdy1(NCOL, NCOL), wvdy2(NCOL, NCOL), &
                            wvuz1(NCOL, NCOL), wvuz2(NCOL, NCOL), wvdz1(NCOL, NCOL), wvdz2(NCOL, NCOL)

        ! Angular components
        complex(real64) :: wangl1(NCOL, NCOL), wangl2(NCOL, NCOL), wangl3(NCOL, NCOL), wangl4(NCOL, NCOL), &
                            wangl5(NCOL, NCOL), wangl6(NCOL, NCOL), wangl7(NCOL, NCOL), wangl8(NCOL, NCOL)

        ! ZIG components
        complex(real64) :: zig1w1(NCOL, NCOL), zig2w1(NCOL, NCOL), zig1w2(NCOL, NCOL), zig2w2(NCOL, NCOL), &
                            zig1w3(NCOL, NCOL), zig2w3(NCOL, NCOL), zig1w4(NCOL, NCOL), zig2w4(NCOL, NCOL)

        ! TIG components
        complex(real64) :: tig1w1(NCOL, NCOL), tig2w1(NCOL, NCOL), tig1w2(NCOL, NCOL), tig2w2(NCOL, NCOL), &
                            tig1w3(NCOL, NCOL), tig2w3(NCOL, NCOL), tig1w4(NCOL, NCOL), tig2w4(NCOL, NCOL)

        ! LIN components
        complex(real64) :: lin0(NCOL, NCOL), lin1(NCOL, NCOL), lin2(NCOL, NCOL), lin4(NCOL, NCOL)

        ! PLQ components
        complex(real64) :: plq1(NCOL, NCOL), plq2(NCOL, NCOL), plq3(NCOL, NCOL), plq4(NCOL, NCOL), &
                            plq5(NCOL, NCOL), plq6(NCOL, NCOL), plq7(NCOL, NCOL), plq8(NCOL, NCOL)

        ! DPLQ components
        complex(real64) :: dplq1(NCOL, NCOL), dplq2(NCOL, NCOL), dplq3(NCOL, NCOL), dplq4(NCOL, NCOL), &
                            dplq5(NCOL, NCOL), dplq6(NCOL, NCOL), dplq7(NCOL, NCOL), dplq8(NCOL, NCOL)

        ! PL
        complex(real64) :: pl(NCOL, NCOL)

        ! C-sums
        complex(real64) :: csumn, csums(4), csum2s(4), csum2ws(4), &
                            csumw(4), csum2w(4), csum3w(4), csumup(4), &
                            csumud(4)
        complex(real64) :: csumtt1(8), csumtt2(8), csumtt3(8), csumtt4(8), &
                            csumtt5(8), csumtt6(8), csumtt7(8), csumtt8(8), &
                            csumtt9(8), csumtt10(8), csumtt11(8), csumtt12(16), &
                            csumtt13(8), csumtt14(8)

        complex(real64) :: csumplq8(8, 6), csumplq16(16, 7:15) ! csumplq#

        complex(real64) :: akt1

        complex(real64) :: csumsmom(4,2), csum2smom(4,2), csum2wsmom(4,2), csumwmom(4,2), &
                            csum2wmom(4,2), csum3wmom(4,2), csumupmom(4,2), csumudmom(4,2)

        complex(real64) :: csumplqmom8(8, 2, 6), csumplqmom16(16, 2, 7:15) ! csumplqmom#

        complex(real64) :: csumttmom1(8,2), csumttmom2(8,2), csumttmom3(8,2), csumttmom4(8,2), &
                            csumttmom5(8,2), csumttmom6(8,2), csumttmom7(8,2), csumttmom8(8,2), &
                            csumttmom9(8,2), csumttmom10(8,2), csumttmom11(8,2), csumttmom12(16,2), &
                            csumttmom13(8,2), csumttmom14(8,2)

        complex(real64) :: giot, pf(2)

        !!! Setup THERML1 function
        giot = cmplx(0.0, 1.0, kind=real64)
        ls(1)=LX1
        ls(2)=LX2
        ls(3)=LX3
        lb(1)=1

        ! LB is an array containing the blocking length at each level
        do idd = 2, MAX_BLOCKING_LEVEL
            lb(idd) = 2*lb(idd-1)
        enddo

        ku = 1 ! Set directions -- X
        ju = 2 ! Y
        iu = 3 ! Z

        id = blocking_level ! Current blocking level (BL)

        ! We try and fit the blocking level lengths into LX1
        ! Fill lcnt array
        ! Example: LX1 = 20, IBL=4
        ! lcnt = 0, 0, 1, 2

        do ib = 1, id
            lcnt(ib) = 0
        enddo

        lrest = ls(ku) ! LX1

        do idg = id, 1, -1
            lcnt(idg) = lrest / lb(idg)
            lrest = lrest - lcnt(idg) * lb(idg)
        enddo

        ! Find the largest blocking level ids which fits more than
        ! one link in the lattice X direction
        ! (in which case it makes sense to define of a Polyakov loop)
        do ib = 1, id
            if (lcnt(ib) >= 1) then
                ids = ib
            end if
        end do

        ! Blocking level below ids - if it's the smallest, set it to itself
        ! Only happens if LX1 = 2 (I think)
        idsm1 = ids - 1
        if (ids == 1) idsm1 = ids

        !**********************************************************************
        ! Define our configuration at the current blocking level
        ! to be the configuration at BL corresponding to ids
        gauge_field = gauge_field_blocked(:,:,:,:,ids)

        !**********************************************************************
        ! These sums keep track of the operators, summing them to make
        ! the final operator is translation invariant
        csumn = cmplx(0.0, 0.0, kind=real64)
        ! These operators have 4 rotations
                                ! Each element will be summed over
                                ! all x
                                ! and contains one ''rotation''
        csums = cmplx(0.0, 0.0, kind=real64)
        csum2s = cmplx(0.0, 0.0, kind=real64)
        csum2ws = cmplx(0.0, 0.0, kind=real64)
        csumw = cmplx(0.0, 0.0, kind=real64)
        csum2w = cmplx(0.0, 0.0, kind=real64)
        csum3w = cmplx(0.0, 0.0, kind=real64)
        csumup = cmplx(0.0, 0.0, kind=real64)
        csumud = cmplx(0.0, 0.0, kind=real64)

        ! These operators have 4 rotations x 2 reflections
        csumtt1 = cmplx(0.0, 0.0, kind=real64)
        csumtt2 = cmplx(0.0, 0.0, kind=real64)
        csumtt3 = cmplx(0.0, 0.0, kind=real64)
        csumtt4 = cmplx(0.0, 0.0, kind=real64)
        csumtt5 = cmplx(0.0, 0.0, kind=real64)
        csumtt6 = cmplx(0.0, 0.0, kind=real64)
        csumtt7 = cmplx(0.0, 0.0, kind=real64)
        csumtt8 = cmplx(0.0, 0.0, kind=real64)
        csumtt9 = cmplx(0.0, 0.0, kind=real64)
        csumtt10 = cmplx(0.0, 0.0, kind=real64)
        csumtt11 = cmplx(0.0, 0.0, kind=real64)
        csumtt13 = cmplx(0.0, 0.0, kind=real64)
        csumtt14 = cmplx(0.0, 0.0, kind=real64)
        csumplq8 = cmplx(0.0, 0.0, kind=real64)
        ! These operators have 4 rotations x 2 reflections
        !                                  x 2 reflection
        ! from different parities
        csumtt12 = cmplx(0.0, 0.0, kind=real64)
        csumplq16 = cmplx(0.0, 0.0, kind=real64)

        ! Same as operators above but adding momentum
        csumsmom = cmplx(0.0, 0.0, kind=real64)
        csum2smom = cmplx(0.0, 0.0, kind=real64)
        csum2wsmom = cmplx(0.0, 0.0, kind=real64)
        csumwmom = cmplx(0.0, 0.0, kind=real64)
        csum2wmom = cmplx(0.0, 0.0, kind=real64)
        csum3wmom = cmplx(0.0, 0.0, kind=real64)
        csumupmom = cmplx(0.0, 0.0, kind=real64)
        csumudmom = cmplx(0.0, 0.0, kind=real64)

        csumttmom1 = cmplx(0.0, 0.0, kind=real64)
        csumttmom2 = cmplx(0.0, 0.0, kind=real64)
        csumttmom3 = cmplx(0.0, 0.0, kind=real64)
        csumttmom4 = cmplx(0.0, 0.0, kind=real64)
        csumttmom5 = cmplx(0.0, 0.0, kind=real64)
        csumttmom6 = cmplx(0.0, 0.0, kind=real64)
        csumttmom7 = cmplx(0.0, 0.0, kind=real64)
        csumttmom8 = cmplx(0.0, 0.0, kind=real64)
        csumttmom9 = cmplx(0.0, 0.0, kind=real64)
        csumttmom10 = cmplx(0.0, 0.0, kind=real64)
        csumttmom11 = cmplx(0.0, 0.0, kind=real64)
        csumttmom13 = cmplx(0.0, 0.0, kind=real64)
        csumttmom14 = cmplx(0.0, 0.0, kind=real64)

        csumplqmom8 = cmplx(0.0, 0.0, kind=real64)

        ! Note: While P_// is not well defined here due to J != 0,
        ! we double the operators for cross checking and
        ! extra statistics
        csumttmom12 = cmplx(0.0, 0.0, kind=real64)
        csumplqmom16 = cmplx(0.0, 0.0, kind=real64)

        !**********************************************************************
        ! Sum loops over lattice
        !**********************************************************************

        ! Tabulate all square pulses and plaquettes (all sites, all blocking levels) for flux along
        ! direction ku (= 1). loop_N reads them through square_pulse() and plaquette()
        call get_square_pulses(gauge_field_blocked, ku)
        call get_plaquettes(gauge_field_blocked, ku)

        ! One team processes the full slice. Each thread accumulates private partial
        ! sums over many sites; OpenMP combines them once at loop completion.
        !$omp parallel do default(shared) schedule(static) &
        !$omp&   private(site_ku, site_ju, site_iu, ix, mn, phase, pf, &
        !$omp&      iloop, i, ic, nc, ig, idg, mu, irem, li, ico, ieee, iddd, &
        !$omp&      m1, m2, m3, ml, A11, B11, C11, D11, E11, UINT11, REM11, &
        !$omp&      lin0, lin1, lin2, lin4, akt1) &
        !$omp&   reduction(+:csumn, &
        !$omp&      csums, &
        !$omp&      csum2s, &
        !$omp&      csum2ws, &
        !$omp&      csumw, &
        !$omp&      csum2w, &
        !$omp&      csum3w, &
        !$omp&      csumup, &
        !$omp&      csumud, &
        !$omp&      csumtt1, &
        !$omp&      csumtt2, &
        !$omp&      csumtt3, &
        !$omp&      csumtt4, &
        !$omp&      csumtt5, &
        !$omp&      csumtt6, &
        !$omp&      csumtt7, &
        !$omp&      csumtt8, &
        !$omp&      csumtt9, &
        !$omp&      csumtt10, &
        !$omp&      csumtt11, &
        !$omp&      csumtt12, &
        !$omp&      csumtt13, &
        !$omp&      csumtt14, &
        !$omp&      csumplq8, &
        !$omp&      csumplq16, &
        !$omp&      csumsmom, &
        !$omp&      csum2smom, &
        !$omp&      csum2wsmom, &
        !$omp&      csumwmom, &
        !$omp&      csum2wmom, &
        !$omp&      csum3wmom, &
        !$omp&      csumupmom, &
        !$omp&      csumudmom, &
        !$omp&      csumttmom1, &
        !$omp&      csumttmom2, &
        !$omp&      csumttmom3, &
        !$omp&      csumttmom4, &
        !$omp&      csumttmom5, &
        !$omp&      csumttmom6, &
        !$omp&      csumttmom7, &
        !$omp&      csumttmom8, &
        !$omp&      csumttmom9, &
        !$omp&      csumttmom10, &
        !$omp&      csumttmom11, &
        !$omp&      csumttmom12, &
        !$omp&      csumttmom13, &
        !$omp&      csumttmom14, &
        !$omp&      csumplqmom8, &
        !$omp&      csumplqmom16)
        do site = 1, SLICE_VOLUME
            ! Preserve the original traversal: x is outermost, z innermost.
            site_iu = mod(site - 1, ls(iu)) + 1
            site_ju = mod((site - 1) / ls(iu), ls(ju)) + 1
            site_ku = (site - 1) / (ls(iu) * ls(ju)) + 1
            phase = 2.0 * PI * site_ku / real(ls(ku))
            pf(1) = cmplx(cos(phase), sin(phase), kind=real64)
            phase = phase * 2.0d0
            pf(2) = cmplx(cos(phase), sin(phase), kind=real64)

            ix(iu) = site_iu
            ix(ju) = site_ju
            ix(ku) = site_ku
            mn = ix(1) + ls(1)*(ix(2)-1) + ls(1)*ls(2)*(ix(3)-1)

! TODO: Remove after cleanup
!        site = 0 ! Initial lattice point
!        do site_ku = 1, ls(ku) ! LX direction
!            ! Assign momentum phases
!            phase = 2.0 * PI * site_ku / real(ls(ku))
!            pf(1) = cmplx(cos(phase), sin(phase), kind=real64)
!            phase = phase * 2.0d0
!            PF(2)=cmplx(cos(phase), sin(phase), kind=real64)
!
!            do site_ju = 1, ls(ju) ! LY direction
!                do site_iu = 1, ls(iu) ! LZ direction
!                    ! Lexicographical lattice site definition
!                    site = site + 1
!                    ix(iu) = site_iu
!                    ix(ju) = site_ju
!                    ix(ku) = site_ku
!                    ! MN defines the current lattice site
!                    mn = ix(1) + ls(1)*(ix(2)-1) + ls(1)*ls(2)*(ix(3)-1)
!

                    ! TOdo: ANDREAS COMMENTS FROM HERE
                    !**********************************************************************
                    if (ids == id) then
                        irem = 3
                        if (ids == 1) irem = 4

                        do iloop = 1, irem
                            if (iloop == 1) li = lcnt(ids)
                            if (iloop > 1) li = lcnt(ids ) - 2**(iloop-2)
                            if (li < 0) cycle

                            m2 = mn
                            if (iloop > 1) then
                                do i = 1, 2**(iloop-2)
                                    m3 = move(m2, ku, ids)
                                    m2 = m3
                                end do
                            endif

                            E11 = cmplx(0.0, 0.0, kind=real64)
                            do ic = 1, NCOL
                                E11(ic, ic) = cmplx(1.0, 0.0, kind=real64)
                            enddo

                            do nc = 1, li
                                B11 = gauge_field(:, :, m2, ku)
                                C11 = matmul(E11, B11)
                                m3 = move(m2, ku, ids)
                                m2 = m3
                                E11 = C11
                            enddo

                            ml = m2
                            select case(iloop)
                            case(1)
                                lin0 = E11
                            case(2)
                                lin1 = E11
                            case(3)
                                lin2 = E11
                            case(4)
                                lin4 = E11
                            case default
                                error stop "iloop not a number 1,...,4"
                            end select
                        enddo
                    else
                        m2=mn
                    endif

                    !**********************************************************************
                    !*********************** REMAINING PIECE ******************************
                    !**********************************************************************

                    REM11 = cmplx(0.0, 0.0, kind=real64)
                    do ic = 1, NCOL
                        REM11(ic, ic) = cmplx(1.0, 0.0, kind=real64)
                    enddo

                    do ig=1,ids
                        idg = ids - ig + 1
                        if (idg == id) cycle
                        do nc = 1, lcnt(idg)
                            UINT11 = cmplx(0.0, 0.0, kind=real64)

                            do mu = 1, 3
                                if(mu == ku) cycle

                                B11 = gauge_field_blocked(:, :, m2, mu, idsm1)
                                m3 = move(m2, mu, idsm1)
                                C11 = gauge_field_blocked(:, :, m3, ku, idg)
                                D11 = matmul(B11, C11)
                                m1 = move(m2, ku, idg)
                                B11 = herm(gauge_field_blocked(:, :, m1, mu, idsm1))
                                C11 = matmul(D11, B11)
                                UINT11 = UINT11 + C11
                                m3 = move(m2, -mu, idsm1)
                                B11 = herm(gauge_field_blocked(:, :, m3, mu, idsm1))
                                C11 = gauge_field_blocked(:, :, m3, ku, idg)
                                D11 = matmul(B11, C11)
                                m1 = move(m3, ku, idg)
                                B11 = gauge_field_blocked(:, :, m1, mu, idsm1)
                                C11 = matmul(D11, B11)
                                UINT11 = UINT11 + C11
                            enddo
                            B11 = gauge_field_blocked(:, :, m2, ku, idg)
                            UINT11 = UINT11 + B11
                            B11 = normalise_link(UINT11)
                            C11 = matmul(REM11, B11)
                            m1 = m2
                            REM11 = C11
                            m2 = move(m1, ku, idg)
                        enddo
                    enddo

                    !! OPERATOR CONSTRUCTION
                    ! Every operator writes to the accumulators; many operators share
                    ! a rotation/reflection slot, so all sums need OpenMP reductions.
                    do iddd = 1, 337 !new!
                        ! The serial routine leaves IEEE unchanged for a short-loop
                        ! fallback. Its final successfully evaluated operator sets it to 16.
                        ! Carry that channel forward explicitly to preserve its projection.
                        if (iddd == 1) ieee = 16
                        m2 = mn
                        A11 = cmplx(0.0, 0.0, kind=real64)
                        do ic = 1,NCOL
                            A11(ic, ic) = cmplx(1.0, 0.0, kind=real64)
                        enddo

                        if (ids == id) then
                            ico=2
                            if (iddd < 5) then
                                ico=1
                            elseif (((iddd > 12).and.(iddd < 17)).and.(ids /= 1)) then
                                ico=1
                            elseif (((iddd > 16).and.(iddd < 21)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 20).and.(iddd < 25)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 24).and.(iddd < 29)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 28).and.(iddd < 33)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 40).and.(iddd < 49)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 48).and.(iddd < 57)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 56).and.(iddd < 65)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 64).and.(iddd < 69)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 68).and.(iddd < 73)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 72).and.(iddd < 81)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 80).and.(iddd < 89)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 88).and.(iddd < 97)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 96).and.(iddd < 105)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 104).and.(iddd < 113)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 112).and.(iddd < 129)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 128).and.(iddd < 137)).and.(ids == 1)) then
                                ico=4
                            elseif (((iddd > 136).and.(iddd < 145)).and.(ids /= 1)) then
                                ico=1
                            elseif (iddd == 145) then
                                ico=1
                            elseif ((iddd > 145).and.(iddd < 194)) then
                                ico=1
                            elseif ((iddd > 209).and.(iddd < 338)) then
                                ico=1 ! THIS NEEDS TO BE FIXED  !
                            else
                                continue
                            endif

                            ! Match the updated serial kernel: a short loop falls back to
                            ! the undeformed Polyakov segment; fitting operators get their
                            ! requested extra length after the operator shape is built.
                            if (lcnt(ids) < ico) then
                                A11 = LIN0
                                M2 = ML
                            else
                                ! The serial operator cases assign IEEE only when the full
                                ! operator is calculated; fallback operators retain its last value.
                                ieee = ieee_of(iddd)
                                ! Select operator to calculate
                                select case(iddd)
                                case(1)
                                    call loop_1(m2, m3, A11, ids, gauge_field_blocked)
                                case(2)
                                    call loop_2(m2, m3, A11, ids, gauge_field_blocked)
                                case(3)
                                    call loop_3(m2, m3, A11, ids, gauge_field_blocked)
                                case(4)
                                    call loop_4(m2, m3, A11, ids, gauge_field_blocked)
                                case(5)
                                    call loop_5(m2, m3, A11, ids, gauge_field_blocked)
                                case(6)
                                    call loop_6(m2, m3, A11, ids, gauge_field_blocked)
                                case(7)
                                    call loop_7(m2, m3, A11, ids, gauge_field_blocked)
                                case(8)
                                    call loop_8(m2, m3, A11, ids, gauge_field_blocked)
                                case(9)
                                    call loop_9(m2, m3, A11, ids, gauge_field_blocked)
                                case(10)
                                    call loop_10(m2, m3, A11, ids, gauge_field_blocked)
                                case(11)
                                    call loop_11(m2, m3, A11, ids, gauge_field_blocked)
                                case(12)
                                    call loop_12(m2, m3, A11, ids, gauge_field_blocked)
                                case(13)
                                    call loop_13(m2, m3, A11, ids, gauge_field_blocked)
                                case(14)
                                    call loop_14(m2, m3, A11, ids, gauge_field_blocked)
                                case(15)
                                    call loop_15(m2, m3, A11, ids, gauge_field_blocked)
                                case(16)
                                    call loop_16(m2, m3, A11, ids, gauge_field_blocked)
                                case(17)
                                    call loop_17(m2, m3, A11, ids, gauge_field_blocked)
                                case(18)
                                    call loop_18(m2, m3, A11, ids, gauge_field_blocked)
                                case(19)
                                    call loop_19(m2, m3, A11, ids, gauge_field_blocked)
                                case(20)
                                    call loop_20(m2, m3, A11, ids, gauge_field_blocked)
                                case(21)
                                    call loop_21(m2, m3, A11, ids, gauge_field_blocked)
                                case(22)
                                    call loop_22(m2, m3, A11, ids, gauge_field_blocked)
                                case(23)
                                    call loop_23(m2, m3, A11, ids, gauge_field_blocked)
                                case(24)
                                    call loop_24(m2, m3, A11, ids, gauge_field_blocked)
                                case(25)
                                    call loop_25(m2, m3, A11, ids, gauge_field_blocked)
                                case(26)
                                    call loop_26(m2, m3, A11, ids, gauge_field_blocked)
                                case(27)
                                    call loop_27(m2, m3, A11, ids, gauge_field_blocked)
                                case(28)
                                    call loop_28(m2, m3, A11, ids, gauge_field_blocked)
                                case(29)
                                    call loop_29(m2, m3, A11, ids, gauge_field_blocked)
                                case(30)
                                    call loop_30(m2, m3, A11, ids, gauge_field_blocked)
                                case(31)
                                    call loop_31(m2, m3, A11, ids, gauge_field_blocked)
                                case(32)
                                    call loop_32(m2, m3, A11, ids, gauge_field_blocked)
                                case(33)
                                    call loop_33(m2, m3, A11, ids, gauge_field_blocked)
                                case(34)
                                    call loop_34(m2, m3, A11, ids, gauge_field_blocked)
                                case(35)
                                    call loop_35(m2, m3, A11, ids, gauge_field_blocked)
                                case(36)
                                    call loop_36(m2, m3, A11, ids, gauge_field_blocked)
                                case(37)
                                    call loop_37(m2, m3, A11, ids, gauge_field_blocked)
                                case(38)
                                    call loop_38(m2, m3, A11, ids, gauge_field_blocked)
                                case(39)
                                    call loop_39(m2, m3, A11, ids, gauge_field_blocked)
                                case(40)
                                    call loop_40(m2, m3, A11, ids, gauge_field_blocked)
                                case(41)
                                    call loop_41(m2, m3, A11, ids, gauge_field_blocked)
                                case(42)
                                    call loop_42(m2, m3, A11, ids, gauge_field_blocked)
                                case(43)
                                    call loop_43(m2, m3, A11, ids, gauge_field_blocked)
                                case(44)
                                    call loop_44(m2, m3, A11, ids, gauge_field_blocked)
                                case(45)
                                    call loop_45(m2, m3, A11, ids, gauge_field_blocked)
                                case(46)
                                    call loop_46(m2, m3, A11, ids, gauge_field_blocked)
                                case(47)
                                    call loop_47(m2, m3, A11, ids, gauge_field_blocked)
                                case(48)
                                    call loop_48(m2, m3, A11, ids, gauge_field_blocked)
                                case(49)
                                    call loop_49(m2, m3, A11, ids, gauge_field_blocked)
                                case(50)
                                    call loop_50(m2, m3, A11, ids, gauge_field_blocked)
                                case(51)
                                    call loop_51(m2, m3, A11, ids, gauge_field_blocked)
                                case(52)
                                    call loop_52(m2, m3, A11, ids, gauge_field_blocked)
                                case(53)
                                    call loop_53(m2, m3, A11, ids, gauge_field_blocked)
                                case(54)
                                    call loop_54(m2, m3, A11, ids, gauge_field_blocked)
                                case(55)
                                    call loop_55(m2, m3, A11, ids, gauge_field_blocked)
                                case(56)
                                    call loop_56(m2, m3, A11, ids, gauge_field_blocked)
                                case(57)
                                    call loop_57(m2, m3, A11, ids, gauge_field_blocked)
                                case(58)
                                    call loop_58(m2, m3, A11, ids, gauge_field_blocked)
                                case(59)
                                    call loop_59(m2, m3, A11, ids, gauge_field_blocked)
                                case(60)
                                    call loop_60(m2, m3, A11, ids, gauge_field_blocked)
                                case(61)
                                    call loop_61(m2, m3, A11, ids, gauge_field_blocked)
                                case(62)
                                    call loop_62(m2, m3, A11, ids, gauge_field_blocked)
                                case(63)
                                    call loop_63(m2, m3, A11, ids, gauge_field_blocked)
                                case(64)
                                    call loop_64(m2, m3, A11, ids, gauge_field_blocked)
                                case(65)
                                    call loop_65(m2, m3, A11, ids, gauge_field_blocked)
                                case(66)
                                    call loop_66(m2, m3, A11, ids, gauge_field_blocked)
                                case(67)
                                    call loop_67(m2, m3, A11, ids, gauge_field_blocked)
                                case(68)
                                    call loop_68(m2, m3, A11, ids, gauge_field_blocked)
                                case(69)
                                    call loop_69(m2, m3, A11, ids, gauge_field_blocked)
                                case(70)
                                    call loop_70(m2, m3, A11, ids, gauge_field_blocked)
                                case(71)
                                    call loop_71(m2, m3, A11, ids, gauge_field_blocked)
                                case(72)
                                    call loop_72(m2, m3, A11, ids, gauge_field_blocked)
                                case(73)
                                    call loop_73(m2, m3, A11, ids, gauge_field_blocked)
                                case(74)
                                    call loop_74(m2, m3, A11, ids, gauge_field_blocked)
                                case(75)
                                    call loop_75(m2, m3, A11, ids, gauge_field_blocked)
                                case(76)
                                    call loop_76(m2, m3, A11, ids, gauge_field_blocked)
                                case(77)
                                    call loop_77(m2, m3, A11, ids, gauge_field_blocked)
                                case(78)
                                    call loop_78(m2, m3, A11, ids, gauge_field_blocked)
                                case(79)
                                    call loop_79(m2, m3, A11, ids, gauge_field_blocked)
                                case(80)
                                    call loop_80(m2, m3, A11, ids, gauge_field_blocked)
                                case(81)
                                    call loop_81(m2, m3, A11, ids, gauge_field_blocked)
                                case(82)
                                    call loop_82(m2, m3, A11, ids, gauge_field_blocked)
                                case(83)
                                    call loop_83(m2, m3, A11, ids, gauge_field_blocked)
                                case(84)
                                    call loop_84(m2, m3, A11, ids, gauge_field_blocked)
                                case(85)
                                    call loop_85(m2, m3, A11, ids, gauge_field_blocked)
                                case(86)
                                    call loop_86(m2, m3, A11, ids, gauge_field_blocked)
                                case(87)
                                    call loop_87(m2, m3, A11, ids, gauge_field_blocked)
                                case(88)
                                    call loop_88(m2, m3, A11, ids, gauge_field_blocked)
                                case(89)
                                    call loop_89(m2, m3, A11, ids, gauge_field_blocked)
                                case(90)
                                    call loop_90(m2, m3, A11, ids, gauge_field_blocked)
                                case(91)
                                    call loop_91(m2, m3, A11, ids, gauge_field_blocked)
                                case(92)
                                    call loop_92(m2, m3, A11, ids, gauge_field_blocked)
                                case(93)
                                    call loop_93(m2, m3, A11, ids, gauge_field_blocked)
                                case(94)
                                    call loop_94(m2, m3, A11, ids, gauge_field_blocked)
                                case(95)
                                    call loop_95(m2, m3, A11, ids, gauge_field_blocked)
                                case(96)
                                    call loop_96(m2, m3, A11, ids, gauge_field_blocked)
                                case(97)
                                    call loop_97(m2, m3, A11, ids, gauge_field_blocked)
                                case(98)
                                    call loop_98(m2, m3, A11, ids, gauge_field_blocked)
                                case(99)
                                    call loop_99(m2, m3, A11, ids, gauge_field_blocked)
                                case(100)
                                    call loop_100(m2, m3, A11, ids, gauge_field_blocked)
                                case(101)
                                    call loop_101(m2, m3, A11, ids, gauge_field_blocked)
                                case(102)
                                    call loop_102(m2, m3, A11, ids, gauge_field_blocked)
                                case(103)
                                    call loop_103(m2, m3, A11, ids, gauge_field_blocked)
                                case(104)
                                    call loop_104(m2, m3, A11, ids, gauge_field_blocked)
                                case(105)
                                    call loop_105(m2, m3, A11, ids, gauge_field_blocked)
                                case(106)
                                    call loop_106(m2, m3, A11, ids, gauge_field_blocked)
                                case(107)
                                    call loop_107(m2, m3, A11, ids, gauge_field_blocked)
                                case(108)
                                    call loop_108(m2, m3, A11, ids, gauge_field_blocked)
                                case(109)
                                    call loop_109(m2, m3, A11, ids, gauge_field_blocked)
                                case(110)
                                    call loop_110(m2, m3, A11, ids, gauge_field_blocked)
                                case(111)
                                    call loop_111(m2, m3, A11, ids, gauge_field_blocked)
                                case(112)
                                    call loop_112(m2, m3, A11, ids, gauge_field_blocked)
                                case(113)
                                    call loop_113(m2, m3, A11, ids, gauge_field_blocked)
                                case(114)
                                    call loop_114(m2, m3, A11, ids, gauge_field_blocked)
                                case(115)
                                    call loop_115(m2, m3, A11, ids, gauge_field_blocked)
                                case(116)
                                    call loop_116(m2, m3, A11, ids, gauge_field_blocked)
                                case(117)
                                    call loop_117(m2, m3, A11, ids, gauge_field_blocked)
                                case(118)
                                    call loop_118(m2, m3, A11, ids, gauge_field_blocked)
                                case(119)
                                    call loop_119(m2, m3, A11, ids, gauge_field_blocked)
                                case(120)
                                    call loop_120(m2, m3, A11, ids, gauge_field_blocked)
                                case(121)
                                    call loop_121(m2, m3, A11, ids, gauge_field_blocked)
                                case(122)
                                    call loop_122(m2, m3, A11, ids, gauge_field_blocked)
                                case(123)
                                    call loop_123(m2, m3, A11, ids, gauge_field_blocked)
                                case(124)
                                    call loop_124(m2, m3, A11, ids, gauge_field_blocked)
                                case(125)
                                    call loop_125(m2, m3, A11, ids, gauge_field_blocked)
                                case(126)
                                    call loop_126(m2, m3, A11, ids, gauge_field_blocked)
                                case(127)
                                    call loop_127(m2, m3, A11, ids, gauge_field_blocked)
                                case(128)
                                    call loop_128(m2, m3, A11, ids, gauge_field_blocked)
                                case(129)
                                    call loop_129(m2, m3, A11, ids, gauge_field_blocked)
                                case(130)
                                    call loop_130(m2, m3, A11, ids, gauge_field_blocked)
                                case(131)
                                    call loop_131(m2, m3, A11, ids, gauge_field_blocked)
                                case(132)
                                    call loop_132(m2, m3, A11, ids, gauge_field_blocked)
                                case(133)
                                    call loop_133(m2, m3, A11, ids, gauge_field_blocked)
                                case(134)
                                    call loop_134(m2, m3, A11, ids, gauge_field_blocked)
                                case(135)
                                    call loop_135(m2, m3, A11, ids, gauge_field_blocked)
                                case(136)
                                    call loop_136(m2, m3, A11, ids, gauge_field_blocked)
                                case(137)
                                    call loop_137(m2, m3, A11, ids, gauge_field_blocked)
                                case(138)
                                    call loop_138(m2, m3, A11, ids, gauge_field_blocked)
                                case(139)
                                    call loop_139(m2, m3, A11, ids, gauge_field_blocked)
                                case(140)
                                    call loop_140(m2, m3, A11, ids, gauge_field_blocked)
                                case(141)
                                    call loop_141(m2, m3, A11, ids, gauge_field_blocked)
                                case(142)
                                    call loop_142(m2, m3, A11, ids, gauge_field_blocked)
                                case(143)
                                    call loop_143(m2, m3, A11, ids, gauge_field_blocked)
                                case(144)
                                    call loop_144(m2, m3, A11, ids, gauge_field_blocked)
                                case(145)
                                    call loop_145(m2, m3, A11, ids, gauge_field_blocked)
                                case(146)
                                    call loop_146(m2, m3, A11, ids, gauge_field_blocked)
                                case(147)
                                    call loop_147(m2, m3, A11, ids, gauge_field_blocked)
                                case(148)
                                    call loop_148(m2, m3, A11, ids, gauge_field_blocked)
                                case(149)
                                    call loop_149(m2, m3, A11, ids, gauge_field_blocked)
                                case(150)
                                    call loop_150(m2, m3, A11, ids, gauge_field_blocked)
                                case(151)
                                    call loop_151(m2, m3, A11, ids, gauge_field_blocked)
                                case(152)
                                    call loop_152(m2, m3, A11, ids, gauge_field_blocked)
                                case(153)
                                    call loop_153(m2, m3, A11, ids, gauge_field_blocked)
                                case(154)
                                    call loop_154(m2, m3, A11, ids, gauge_field_blocked)
                                case(155)
                                    call loop_155(m2, m3, A11, ids, gauge_field_blocked)
                                case(156)
                                    call loop_156(m2, m3, A11, ids, gauge_field_blocked)
                                case(157)
                                    call loop_157(m2, m3, A11, ids, gauge_field_blocked)
                                case(158)
                                    call loop_158(m2, m3, A11, ids, gauge_field_blocked)
                                case(159)
                                    call loop_159(m2, m3, A11, ids, gauge_field_blocked)
                                case(160)
                                    call loop_160(m2, m3, A11, ids, gauge_field_blocked)
                                case(161)
                                    call loop_161(m2, m3, A11, ids, gauge_field_blocked)
                                case(162)
                                    call loop_162(m2, m3, A11, ids, gauge_field_blocked)
                                case(163)
                                    call loop_163(m2, m3, A11, ids, gauge_field_blocked)
                                case(164)
                                    call loop_164(m2, m3, A11, ids, gauge_field_blocked)
                                case(165)
                                    call loop_165(m2, m3, A11, ids, gauge_field_blocked)
                                case(166)
                                    call loop_166(m2, m3, A11, ids, gauge_field_blocked)
                                case(167)
                                    call loop_167(m2, m3, A11, ids, gauge_field_blocked)
                                case(168)
                                    call loop_168(m2, m3, A11, ids, gauge_field_blocked)
                                case(169)
                                    call loop_169(m2, m3, A11, ids, gauge_field_blocked)
                                case(170)
                                    call loop_170(m2, m3, A11, ids, gauge_field_blocked)
                                case(171)
                                    call loop_171(m2, m3, A11, ids, gauge_field_blocked)
                                case(172)
                                    call loop_172(m2, m3, A11, ids, gauge_field_blocked)
                                case(173)
                                    call loop_173(m2, m3, A11, ids, gauge_field_blocked)
                                case(174)
                                    call loop_174(m2, m3, A11, ids, gauge_field_blocked)
                                case(175)
                                    call loop_175(m2, m3, A11, ids, gauge_field_blocked)
                                case(176)
                                    call loop_176(m2, m3, A11, ids, gauge_field_blocked)
                                case(177)
                                    call loop_177(m2, m3, A11, ids, gauge_field_blocked)
                                case(178)
                                    call loop_178(m2, m3, A11, ids, gauge_field_blocked)
                                case(179)
                                    call loop_179(m2, m3, A11, ids, gauge_field_blocked)
                                case(180)
                                    call loop_180(m2, m3, A11, ids, gauge_field_blocked)
                                case(181)
                                    call loop_181(m2, m3, A11, ids, gauge_field_blocked)
                                case(182)
                                    call loop_182(m2, m3, A11, ids, gauge_field_blocked)
                                case(183)
                                    call loop_183(m2, m3, A11, ids, gauge_field_blocked)
                                case(184)
                                    call loop_184(m2, m3, A11, ids, gauge_field_blocked)
                                case(185)
                                    call loop_185(m2, m3, A11, ids, gauge_field_blocked)
                                case(186)
                                    call loop_186(m2, m3, A11, ids, gauge_field_blocked)
                                case(187)
                                    call loop_187(m2, m3, A11, ids, gauge_field_blocked)
                                case(188)
                                    call loop_188(m2, m3, A11, ids, gauge_field_blocked)
                                case(189)
                                    call loop_189(m2, m3, A11, ids, gauge_field_blocked)
                                case(190)
                                    call loop_190(m2, m3, A11, ids, gauge_field_blocked)
                                case(191)
                                    call loop_191(m2, m3, A11, ids, gauge_field_blocked)
                                case(192)
                                    call loop_192(m2, m3, A11, ids, gauge_field_blocked)
                                case(193)
                                    call loop_193(m2, m3, A11, ids, gauge_field_blocked)
                                case(194)
                                    call loop_194(m2, m3, A11, ids, gauge_field_blocked)
                                case(195)
                                    call loop_195(m2, m3, A11, ids, gauge_field_blocked)
                                case(196)
                                    call loop_196(m2, m3, A11, ids, gauge_field_blocked)
                                case(197)
                                    call loop_197(m2, m3, A11, ids, gauge_field_blocked)
                                case(198)
                                    call loop_198(m2, m3, A11, ids, gauge_field_blocked)
                                case(199)
                                    call loop_199(m2, m3, A11, ids, gauge_field_blocked)
                                case(200)
                                    call loop_200(m2, m3, A11, ids, gauge_field_blocked)
                                case(201)
                                    call loop_201(m2, m3, A11, ids, gauge_field_blocked)
                                case(202)
                                    call loop_202(m2, m3, A11, ids, gauge_field_blocked)
                                case(203)
                                    call loop_203(m2, m3, A11, ids, gauge_field_blocked)
                                case(204)
                                    call loop_204(m2, m3, A11, ids, gauge_field_blocked)
                                case(205)
                                    call loop_205(m2, m3, A11, ids, gauge_field_blocked)
                                case(206)
                                    call loop_206(m2, m3, A11, ids, gauge_field_blocked)
                                case(207)
                                    call loop_207(m2, m3, A11, ids, gauge_field_blocked)
                                case(208)
                                    call loop_208(m2, m3, A11, ids, gauge_field_blocked)
                                case(209)
                                    call loop_209(m2, m3, A11, ids, gauge_field_blocked)
                                case(210)
                                    call loop_210(m2, m3, A11, ids, gauge_field_blocked)
                                case(211)
                                    call loop_211(m2, m3, A11, ids, gauge_field_blocked)
                                case(212)
                                    call loop_212(m2, m3, A11, ids, gauge_field_blocked)
                                case(213)
                                    call loop_213(m2, m3, A11, ids, gauge_field_blocked)
                                case(214)
                                    call loop_214(m2, m3, A11, ids, gauge_field_blocked)
                                case(215)
                                    call loop_215(m2, m3, A11, ids, gauge_field_blocked)
                                case(216)
                                    call loop_216(m2, m3, A11, ids, gauge_field_blocked)
                                case(217)
                                    call loop_217(m2, m3, A11, ids, gauge_field_blocked)
                                case(218)
                                    call loop_218(m2, m3, A11, ids, gauge_field_blocked)
                                case(219)
                                    call loop_219(m2, m3, A11, ids, gauge_field_blocked)
                                case(220)
                                    call loop_220(m2, m3, A11, ids, gauge_field_blocked)
                                case(221)
                                    call loop_221(m2, m3, A11, ids, gauge_field_blocked)
                                case(222)
                                    call loop_222(m2, m3, A11, ids, gauge_field_blocked)
                                case(223)
                                    call loop_223(m2, m3, A11, ids, gauge_field_blocked)
                                case(224)
                                    call loop_224(m2, m3, A11, ids, gauge_field_blocked)
                                case(225)
                                    call loop_225(m2, m3, A11, ids, gauge_field_blocked)
                                case(226)
                                    call loop_226(m2, m3, A11, ids, gauge_field_blocked)
                                case(227)
                                    call loop_227(m2, m3, A11, ids, gauge_field_blocked)
                                case(228)
                                    call loop_228(m2, m3, A11, ids, gauge_field_blocked)
                                case(229)
                                    call loop_229(m2, m3, A11, ids, gauge_field_blocked)
                                case(230)
                                    call loop_230(m2, m3, A11, ids, gauge_field_blocked)
                                case(231)
                                    call loop_231(m2, m3, A11, ids, gauge_field_blocked)
                                case(232)
                                    call loop_232(m2, m3, A11, ids, gauge_field_blocked)
                                case(233)
                                    call loop_233(m2, m3, A11, ids, gauge_field_blocked)
                                case(234)
                                    call loop_234(m2, m3, A11, ids, gauge_field_blocked)
                                case(235)
                                    call loop_235(m2, m3, A11, ids, gauge_field_blocked)
                                case(236)
                                    call loop_236(m2, m3, A11, ids, gauge_field_blocked)
                                case(237)
                                    call loop_237(m2, m3, A11, ids, gauge_field_blocked)
                                case(238)
                                    call loop_238(m2, m3, A11, ids, gauge_field_blocked)
                                case(239)
                                    call loop_239(m2, m3, A11, ids, gauge_field_blocked)
                                case(240)
                                    call loop_240(m2, m3, A11, ids, gauge_field_blocked)
                                case(241)
                                    call loop_241(m2, m3, A11, ids, gauge_field_blocked)
                                case(242)
                                    call loop_242(m2, m3, A11, ids, gauge_field_blocked)
                                case(243)
                                    call loop_243(m2, m3, A11, ids, gauge_field_blocked)
                                case(244)
                                    call loop_244(m2, m3, A11, ids, gauge_field_blocked)
                                case(245)
                                    call loop_245(m2, m3, A11, ids, gauge_field_blocked)
                                case(246)
                                    call loop_246(m2, m3, A11, ids, gauge_field_blocked)
                                case(247)
                                    call loop_247(m2, m3, A11, ids, gauge_field_blocked)
                                case(248)
                                    call loop_248(m2, m3, A11, ids, gauge_field_blocked)
                                case(249)
                                    call loop_249(m2, m3, A11, ids, gauge_field_blocked)
                                case(250)
                                    call loop_250(m2, m3, A11, ids, gauge_field_blocked)
                                case(251)
                                    call loop_251(m2, m3, A11, ids, gauge_field_blocked)
                                case(252)
                                    call loop_252(m2, m3, A11, ids, gauge_field_blocked)
                                case(253)
                                    call loop_253(m2, m3, A11, ids, gauge_field_blocked)
                                case(254)
                                    call loop_254(m2, m3, A11, ids, gauge_field_blocked)
                                case(255)
                                    call loop_255(m2, m3, A11, ids, gauge_field_blocked)
                                case(256)
                                    call loop_256(m2, m3, A11, ids, gauge_field_blocked)
                                case(257)
                                    call loop_257(m2, m3, A11, ids, gauge_field_blocked)
                                case(258)
                                    call loop_258(m2, m3, A11, ids, gauge_field_blocked)
                                case(259)
                                    call loop_259(m2, m3, A11, ids, gauge_field_blocked)
                                case(260)
                                    call loop_260(m2, m3, A11, ids, gauge_field_blocked)
                                case(261)
                                    call loop_261(m2, m3, A11, ids, gauge_field_blocked)
                                case(262)
                                    call loop_262(m2, m3, A11, ids, gauge_field_blocked)
                                case(263)
                                    call loop_263(m2, m3, A11, ids, gauge_field_blocked)
                                case(264)
                                    call loop_264(m2, m3, A11, ids, gauge_field_blocked)
                                case(265)
                                    call loop_265(m2, m3, A11, ids, gauge_field_blocked)
                                case(266)
                                    call loop_266(m2, m3, A11, ids, gauge_field_blocked)
                                case(267)
                                    call loop_267(m2, m3, A11, ids, gauge_field_blocked)
                                case(268)
                                    call loop_268(m2, m3, A11, ids, gauge_field_blocked)
                                case(269)
                                    call loop_269(m2, m3, A11, ids, gauge_field_blocked)
                                case(270)
                                    call loop_270(m2, m3, A11, ids, gauge_field_blocked)
                                case(271)
                                    call loop_271(m2, m3, A11, ids, gauge_field_blocked)
                                case(272)
                                    call loop_272(m2, m3, A11, ids, gauge_field_blocked)
                                case(273)
                                    call loop_273(m2, m3, A11, ids, gauge_field_blocked)
                                case(274)
                                    call loop_274(m2, m3, A11, ids, gauge_field_blocked)
                                case(275)
                                    call loop_275(m2, m3, A11, ids, gauge_field_blocked)
                                case(276)
                                    call loop_276(m2, m3, A11, ids, gauge_field_blocked)
                                case(277)
                                    call loop_277(m2, m3, A11, ids, gauge_field_blocked)
                                case(278)
                                    call loop_278(m2, m3, A11, ids, gauge_field_blocked)
                                case(279)
                                    call loop_279(m2, m3, A11, ids, gauge_field_blocked)
                                case(280)
                                    call loop_280(m2, m3, A11, ids, gauge_field_blocked)
                                case(281)
                                    call loop_281(m2, m3, A11, ids, gauge_field_blocked)
                                case(282)
                                    call loop_282(m2, m3, A11, ids, gauge_field_blocked)
                                case(283)
                                    call loop_283(m2, m3, A11, ids, gauge_field_blocked)
                                case(284)
                                    call loop_284(m2, m3, A11, ids, gauge_field_blocked)
                                case(285)
                                    call loop_285(m2, m3, A11, ids, gauge_field_blocked)
                                case(286)
                                    call loop_286(m2, m3, A11, ids, gauge_field_blocked)
                                case(287)
                                    call loop_287(m2, m3, A11, ids, gauge_field_blocked)
                                case(288)
                                    call loop_288(m2, m3, A11, ids, gauge_field_blocked)
                                case(289)
                                    call loop_289(m2, m3, A11, ids, gauge_field_blocked)
                                case(290)
                                    call loop_290(m2, m3, A11, ids, gauge_field_blocked)
                                case(291)
                                    call loop_291(m2, m3, A11, ids, gauge_field_blocked)
                                case(292)
                                    call loop_292(m2, m3, A11, ids, gauge_field_blocked)
                                case(293)
                                    call loop_293(m2, m3, A11, ids, gauge_field_blocked)
                                case(294)
                                    call loop_294(m2, m3, A11, ids, gauge_field_blocked)
                                case(295)
                                    call loop_295(m2, m3, A11, ids, gauge_field_blocked)
                                case(296)
                                    call loop_296(m2, m3, A11, ids, gauge_field_blocked)
                                case(297)
                                    call loop_297(m2, m3, A11, ids, gauge_field_blocked)
                                case(298)
                                    call loop_298(m2, m3, A11, ids, gauge_field_blocked)
                                case(299)
                                    call loop_299(m2, m3, A11, ids, gauge_field_blocked)
                                case(300)
                                    call loop_300(m2, m3, A11, ids, gauge_field_blocked)
                                case(301)
                                    call loop_301(m2, m3, A11, ids, gauge_field_blocked)
                                case(302)
                                    call loop_302(m2, m3, A11, ids, gauge_field_blocked)
                                case(303)
                                    call loop_303(m2, m3, A11, ids, gauge_field_blocked)
                                case(304)
                                    call loop_304(m2, m3, A11, ids, gauge_field_blocked)
                                case(305)
                                    call loop_305(m2, m3, A11, ids, gauge_field_blocked)
                                case(306)
                                    call loop_306(m2, m3, A11, ids, gauge_field_blocked)
                                case(307)
                                    call loop_307(m2, m3, A11, ids, gauge_field_blocked)
                                case(308)
                                    call loop_308(m2, m3, A11, ids, gauge_field_blocked)
                                case(309)
                                    call loop_309(m2, m3, A11, ids, gauge_field_blocked)
                                case(310)
                                    call loop_310(m2, m3, A11, ids, gauge_field_blocked)
                                case(311)
                                    call loop_311(m2, m3, A11, ids, gauge_field_blocked)
                                case(312)
                                    call loop_312(m2, m3, A11, ids, gauge_field_blocked)
                                case(313)
                                    call loop_313(m2, m3, A11, ids, gauge_field_blocked)
                                case(314)
                                    call loop_314(m2, m3, A11, ids, gauge_field_blocked)
                                case(315)
                                    call loop_315(m2, m3, A11, ids, gauge_field_blocked)
                                case(316)
                                    call loop_316(m2, m3, A11, ids, gauge_field_blocked)
                                case(317)
                                    call loop_317(m2, m3, A11, ids, gauge_field_blocked)
                                case(318)
                                    call loop_318(m2, m3, A11, ids, gauge_field_blocked)
                                case(319)
                                    call loop_319(m2, m3, A11, ids, gauge_field_blocked)
                                case(320)
                                    call loop_320(m2, m3, A11, ids, gauge_field_blocked)
                                case(321)
                                    call loop_321(m2, m3, A11, ids, gauge_field_blocked)
                                case(322)
                                    call loop_322(m2, m3, A11, ids, gauge_field_blocked)
                                case(323)
                                    call loop_323(m2, m3, A11, ids, gauge_field_blocked)
                                case(324)
                                    call loop_324(m2, m3, A11, ids, gauge_field_blocked)
                                case(325)
                                    call loop_325(m2, m3, A11, ids, gauge_field_blocked)
                                case(326)
                                    call loop_326(m2, m3, A11, ids, gauge_field_blocked)
                                case(327)
                                    call loop_327(m2, m3, A11, ids, gauge_field_blocked)
                                case(328)
                                    call loop_328(m2, m3, A11, ids, gauge_field_blocked)
                                case(329)
                                    call loop_329(m2, m3, A11, ids, gauge_field_blocked)
                                case(330)
                                    call loop_330(m2, m3, A11, ids, gauge_field_blocked)
                                case(331)
                                    call loop_331(m2, m3, A11, ids, gauge_field_blocked)
                                case(332)
                                    call loop_332(m2, m3, A11, ids, gauge_field_blocked)
                                case(333)
                                    call loop_333(m2, m3, A11, ids, gauge_field_blocked)
                                case(334)
                                    call loop_334(m2, m3, A11, ids, gauge_field_blocked)
                                case(335)
                                    call loop_335(m2, m3, A11, ids, gauge_field_blocked)
                                case(336)
                                    call loop_336(m2, m3, A11, ids, gauge_field_blocked)
                                case(337)
                                    call loop_337(m2, m3, A11, ids, gauge_field_blocked)
                                end select
                            !******************************************************************C
                                if (ico == 1) then
                                    C11 = matmul(A11, LIN1)
                                    A11 = C11
                                endif

                                if (ico == 2) then
                                    C11 = matmul(A11, LIN2)
                                    A11 = C11
                                endif

                                if (ico == 4) then
                                    C11 = matmul(A11, LIN4)
                                    A11 = C11
                                endif
                                M2=ML
                            endif
                            !***********************************************************************
                        endif

                        B11 = matmul(A11, REM11)
                        AKT1 = cmplx(0.0, 0.0, kind=real64)
                        do ic = 1, NCOL
                            akt1 = akt1 + B11(ic, ic)
                        enddo
                        !************ ***********************************************************
                        if (iddd <= 4) then
                            csums(ieee) = csums(ieee) + akt1
                            csumsmom(ieee,1) = csumsmom(ieee,1) + akt1*pf(1)
                            csumsmom(ieee,2) = csumsmom(ieee,2) + akt1*pf(2)

                        elseif ((iddd > 4).and.(iddd < 9)) then
                            csum2s(ieee)=csum2s(ieee)+akt1
                            csum2smom(ieee,1)=csum2smom(ieee,1)+akt1*pf(1)
                            csum2smom(ieee,2)=csum2smom(ieee,2)+akt1*pf(2)

                        elseif ((iddd > 8).and.(iddd < 13)) then
                            csum2ws(ieee)=csum2ws(ieee)+akt1
                            csum2wsmom(ieee,1)=csum2wsmom(ieee,1)+akt1*pf(1)
                            csum2wsmom(ieee,2)=csum2wsmom(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 12).and.(iddd < 17)) then
                            csumw(ieee)=csumw(ieee)+akt1
                            csumwmom(ieee,1)=csumwmom(ieee,1)+akt1*pf(1)
                            csumwmom(ieee,2)=csumwmom(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 16).and.(iddd < 21)) then
                            csum2w(ieee)=csum2w(ieee)+akt1
                            csum2wmom(ieee,1)=csum2wmom(ieee,1)+akt1*pf(1)
                            csum2wmom(ieee,2)=csum2wmom(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 20).and.(iddd < 25)) then
                            csum3w(ieee)=csum3w(ieee)+akt1
                            csum3wmom(ieee,1)=csum3wmom(ieee,1)+akt1*pf(1)
                            csum3wmom(ieee,2)=csum3wmom(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 24).and.(iddd < 29)) then
                            csumup(ieee)=csumup(ieee)+akt1
                            csumupmom(ieee,1)=csumupmom(ieee,1)+akt1*pf(1)
                            csumupmom(ieee,2)=csumupmom(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 28).and.(iddd < 33)) then
                            csumud(ieee)=csumud(ieee)+akt1
                            csumudmom(ieee,1)=csumudmom(ieee,1)+akt1*pf(1)
                            csumudmom(ieee,2)=csumudmom(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 32).and.(iddd < 41)) then
                            csumtt1(ieee)=csumtt1(ieee)+akt1
                            csumttmom1(ieee,1)=csumttmom1(ieee,1)+akt1*pf(1)
                            csumttmom1(ieee,2)=csumttmom1(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 40).and.(iddd < 49)) then
                            csumtt2(ieee)=csumtt2(ieee)+akt1
                            csumttmom2(ieee,1)=csumttmom2(ieee,1)+akt1*pf(1)
                            csumttmom2(ieee,2)=csumttmom2(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 48).and.(iddd < 57)) then
                            csumtt3(ieee)=csumtt3(ieee)+akt1
                            csumttmom3(ieee,1)=csumttmom3(ieee,1)+akt1*pf(1)
                            csumttmom3(ieee,2)=csumttmom3(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 56).and.(iddd < 65)) then
                            csumtt4(ieee)=csumtt4(ieee)+akt1
                            csumttmom4(ieee,1)=csumttmom4(ieee,1)+akt1*pf(1)
                            csumttmom4(ieee,2)=csumttmom4(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 64).and.(iddd < 69)) then
                            csumtt5(ieee)=csumtt5(ieee)+akt1
                            csumttmom5(ieee,1)=csumttmom5(ieee,1)+akt1*pf(1)
                            csumttmom5(ieee,2)=csumttmom5(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 68).and.(iddd < 73)) then
                            csumtt6(ieee)=csumtt6(ieee)+akt1
                            csumttmom6(ieee,1)=csumttmom6(ieee,1)+akt1*pf(1)
                            csumttmom6(ieee,2)=csumttmom6(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 72).and.(iddd < 81)) then
                            csumtt7(ieee)=csumtt7(ieee)+akt1
                            csumttmom7(ieee,1)=csumttmom7(ieee,1)+akt1*pf(1)
                            csumttmom7(ieee,2)=csumttmom7(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 80).and.(iddd < 89)) then
                            csumtt8(ieee)=csumtt8(ieee)+akt1
                            csumttmom8(ieee,1)=csumttmom8(ieee,1)+akt1*pf(1)
                            csumttmom8(ieee,2)=csumttmom8(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 88).and.(iddd < 97)) then
                            csumtt9(ieee)=csumtt9(ieee)+akt1
                            csumttmom9(ieee,1)=csumttmom9(ieee,1)+akt1*pf(1)
                            csumttmom9(ieee,2)=csumttmom9(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 96).and.(iddd < 105)) then
                            csumtt10(ieee)=csumtt10(ieee)+akt1
                            csumttmom10(ieee,1)=csumttmom10(ieee,1)+akt1*pf(1)
                            csumttmom10(ieee,2)=csumttmom10(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 104).and.(iddd < 113)) then
                            csumtt11(ieee)=csumtt11(ieee)+akt1
                            csumttmom11(ieee,1)=csumttmom11(ieee,1)+akt1*pf(1)
                            csumttmom11(ieee,2)=csumttmom11(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 112).and.(iddd < 129)) then
                            csumtt12(ieee)=csumtt12(ieee)+akt1
                            csumttmom12(ieee,1)=csumttmom12(ieee,1)+akt1*pf(1)
                            csumttmom12(ieee,2)=csumttmom12(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 128).and.(iddd < 137)) then
                            csumtt13(ieee)=csumtt13(ieee)+akt1
                            csumttmom13(ieee,1)=csumttmom13(ieee,1)+akt1*pf(1)
                            csumttmom13(ieee,2)=csumttmom13(ieee,2)+akt1*pf(2)
                        

                        elseif ((iddd > 136).and.(iddd < 145)) then
                            csumtt14(ieee)=csumtt14(ieee)+akt1
                            csumttmom14(ieee,1)=csumttmom14(ieee,1)+akt1*pf(1)
                            csumttmom14(ieee,2)=csumttmom14(ieee,2)+akt1*pf(2)
                        

                        elseif (iddd == 145) then
                            csumn=csumn+akt1
                        

                        elseif ((iddd > 145).and.(iddd < 154)) then
                            csumplq8(ieee, 1) = csumplq8(ieee, 1) + akt1
                            csumplqmom8(ieee,1,1)=csumplqmom8(ieee,1,1)+akt1*pf(1)
                            csumplqmom8(ieee,2,1)=csumplqmom8(ieee,2,1)+akt1*pf(2)
                        

                        elseif ((iddd > 153).and.(iddd < 162)) then
                            csumplq8(ieee,2)=csumplq8(ieee,2)+akt1
                            csumplqmom8(ieee,1,2)=csumplqmom8(ieee,1,2)+akt1*pf(1)
                            csumplqmom8(ieee, 2, 2)=csumplqmom8(ieee, 2, 2)+akt1*pf(2)
                        

                        elseif ((iddd > 161).and.(iddd < 170)) then
                            csumplq8(ieee,3)=csumplq8(ieee,3)+akt1
                            csumplqmom8(ieee, 1, 3)=csumplqmom8(ieee, 1, 3)+akt1*pf(1)
                            csumplqmom8(ieee, 2, 3)=csumplqmom8(ieee, 2, 3)+akt1*pf(2)
                        

                        elseif ((iddd > 169).and.(iddd < 178)) then
                            csumplq8(ieee,4)=csumplq8(ieee,4)+akt1
                            csumplqmom8(ieee, 1, 4)=csumplqmom8(ieee, 1, 4)+akt1*pf(1)
                            csumplqmom8(ieee, 2, 4)=csumplqmom8(ieee, 2, 4)+akt1*pf(2)
                        

                        elseif ((iddd > 177).and.(iddd < 186)) then
                            csumplq8(ieee,5)=csumplq8(ieee,5)+akt1
                            csumplqmom8(ieee, 1, 5)=csumplqmom8(ieee, 1, 5)+akt1*pf(1)
                            csumplqmom8(ieee, 2, 5)=csumplqmom8(ieee, 2, 5)+akt1*pf(2)
                        

                        elseif ((iddd > 185).and.(iddd < 194)) then
                            csumplq8(ieee,6)=csumplq8(ieee,6)+akt1
                            csumplqmom8(ieee, 1, 6)=csumplqmom8(ieee, 1, 6)+akt1*pf(1)
                            csumplqmom8(ieee, 2, 6)=csumplqmom8(ieee, 2, 6)+akt1*pf(2)
                        

                        elseif ((iddd > 193).and.(iddd < 210)) then
                            csumplq16(ieee,7)=csumplq16(ieee,7)+akt1
                            csumplqmom16(ieee,1,7)=csumplqmom16(ieee,1,7)+akt1*pf(1)
                            csumplqmom16(ieee,2,7)=csumplqmom16(ieee,2,7)+akt1*pf(2)
                        

                        elseif ((iddd > 209).and.(iddd < 226)) then
                            csumplq16(ieee,8)=csumplq16(ieee,8)+akt1
                            csumplqmom16(ieee,1,8)=csumplqmom16(ieee,1,8)+akt1*pf(1)
                            csumplqmom16(ieee,2,8)=csumplqmom16(ieee,2,8)+akt1*pf(2)
                        

                        elseif ((iddd > 225).and.(iddd < 242)) then
                            csumplq16(ieee,9)=csumplq16(ieee,9)+akt1
                            csumplqmom16(ieee,1,9)=csumplqmom16(ieee,1,9)+akt1*pf(1)
                            csumplqmom16(ieee,2,9)=csumplqmom16(ieee,2,9)+akt1*pf(2)
                        

                        elseif ((iddd > 241).and.(iddd < 258)) then
                            csumplq16(ieee,10)=csumplq16(ieee,10)+akt1
                            csumplqmom16(ieee,1,10)=csumplqmom16(ieee,1,10)+akt1*pf(1)
                            csumplqmom16(ieee,2,10)=csumplqmom16(ieee,2,10)+akt1*pf(2)
                        

                        elseif ((iddd > 257).and.(iddd < 274)) then
                            csumplq16(ieee,11)=csumplq16(ieee,11)+akt1
                            csumplqmom16(ieee,1,11)=csumplqmom16(ieee,1,11)+akt1*pf(1)
                            csumplqmom16(ieee,2,11)=csumplqmom16(ieee,2,11)+akt1*pf(2)
                        

                        elseif ((iddd > 273).and.(iddd < 290)) then
                            csumplq16(ieee,12)=csumplq16(ieee,12)+akt1
                            csumplqmom16(ieee,1,12)=csumplqmom16(ieee,1,12)+akt1*pf(1)
                            csumplqmom16(ieee,2,12)=csumplqmom16(ieee,2,12)+akt1*pf(2)
                        

                        elseif ((iddd > 289).and.(iddd < 306)) then
                            csumplq16(ieee,13)=csumplq16(ieee,13)+akt1
                            csumplqmom16(ieee,1,13)=csumplqmom16(ieee,1,13)+akt1*pf(1)
                            csumplqmom16(ieee,2,13)=csumplqmom16(ieee,2,13)+akt1*pf(2)
                        

                        elseif ((iddd > 305).and.(iddd < 322)) then
                            csumplq16(ieee,14)=csumplq16(ieee,14)+akt1
                            csumplqmom16(ieee,1,14)=csumplqmom16(ieee,1,14)+akt1*pf(1)
                            csumplqmom16(ieee,2,14)=csumplqmom16(ieee,2,14)+akt1*pf(2)
                        

                        else
                            csumplq16(ieee,15)=csumplq16(ieee,15)+akt1
                            csumplqmom16(ieee,1,15)=csumplqmom16(ieee,1,15)+akt1*pf(1)
                            csumplqmom16(ieee,2,15)=csumplqmom16(ieee,2,15)+akt1*pf(2)
                        endif

                    enddo
            enddo
        !$omp end parallel do
    !**********************************************************************
        adiv1=1.0/(4.0*ls(ku))
        adiv2=1.0/(8.0*ls(ku))
        adiv3=1.0/(16.0*ls(ku))
        adivn=1.0/ls(ku)
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 1)=csumn*adivn
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 2)=(csums(1)+csums(2)+csums(3)+csums(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 2)=(csumsmom(1,1)+csumsmom(2,1)&
        +csumsmom(3,1)+csumsmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 2)=(csumsmom(1,2)+csumsmom(2,2)&
        +csumsmom(3,2)+csumsmom(4,2))*adiv1
    !      alinemom2(time_slice, id,3)=(csumsmom(1,3)+csumsmom(2,3)&
    !     &+csumsmom(3,3)+csumsmom(4,3))*adiv1
    !      alinemom2(time_slice, id,4)=(csumsmom(1,4)+csumsmom(2,4)&
    !     &+csumsmom(3,4)+csumsmom(4,4))*adiv1
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 3)=(csums(1)+giot*csums(2)-csums(3)-giot*csums(4))&
        *adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 3)=(csumsmom(1,1)+giot*csumsmom(2,1)&
        -csumsmom(3,1)-giot*csumsmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 3)=(csumsmom(1,2)+giot*csumsmom(2,2)&
        -csumsmom(3,2)-giot*csumsmom(4,2))*adiv1
    !      alinemom3(time_slice, id,3)=(csumsmom(1,3)+giot*csumsmom(2,3)&
    !     &-csumsmom(3,3)-giot*csumsmom(4,3))*adiv1
    !      alinemom3(time_slice, id,4)=(csumsmom(1,4)+giot*csumsmom(2,4)&
    !     &-csumsmom(3,4)-giot*csumsmom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 4)=(csums(1)-csums(2)+csums(3)-csums(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 4)=(csumsmom(1,1)-csumsmom(2,1)&
        +csumsmom(3,1)-csumsmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 4)=(csumsmom(1,2)-csumsmom(2,2)&
        +csumsmom(3,2)-csumsmom(4,2))*adiv1
    !      alinemom4(time_slice, id,3)=(csumsmom(1,3)-csumsmom(2,3)&
    !     &+csumsmom(3,3)-csumsmom(4,3))*adiv1
    !      alinemom4(time_slice, id,4)=(csumsmom(1,4)-csumsmom(2,4)&
    !     &+csumsmom(3,4)-csumsmom(4,4))*adiv1
    !**********************************************************************
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 5)=(csum2s(1)+csum2s(2)+csum2s(3)+csum2s(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 5)=(csum2smom(1,1)+csum2smom(2,1)&
        +csum2smom(3,1)+csum2smom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 5)=(csum2smom(1,2)+csum2smom(2,2)&
        +csum2smom(3,2)+csum2smom(4,2))*adiv1
    !      alinemom5(time_slice, id,3)=(csum2smom(1,3)+csum2smom(2,3)&
    !     &+csum2smom(3,3)+csum2smom(4,3))*adiv1
    !      alinemom5(time_slice, id,4)=(csum2smom(1,4)+csum2smom(2,4)&
    !     &+csum2smom(3,4)+csum2smom(4,4))*adiv1
    !***********************************************************************
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 6)=(csum2s(1)+giot*csum2s(2)-csum2s(3)&
        -giot*csum2s(4))*adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 6)=(csum2smom(1,1)+giot*csum2smom(2,1)&
        -csum2smom(3,1)-giot*csum2smom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 6)=(csum2smom(1,2)+giot*csum2smom(2,2)&
        -csum2smom(3,2)-giot*csum2smom(4,2))*adiv1
    !      alinemom6(time_slice, id,3)=(csum2smom(1,3)+giot*csum2smom(2,3)&
    !     &-csum2smom(3,3)-giot*csum2smom(4,3))*adiv1
    !      alinemom6(time_slice, id,4)=(csum2smom(1,4)+giot*csum2smom(2,4)&
    !     &-csum2smom(3,4)-giot*csum2smom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 7)=(csum2s(1)-csum2s(2)+csum2s(3)-csum2s(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 7)=(csum2smom(1,1)-csum2smom(2,1)&
        +csum2smom(3,1)-csum2smom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 7)=(csum2smom(1,2)-csum2smom(2,2)&
        +csum2smom(3,2)-csum2smom(4,2))*adiv1
    !      alinemom7(time_slice, id,3)=(csum2smom(1,3)-csum2smom(2,3)&
    !     &+csum2smom(3,3)-csum2smom(4,3))*adiv1
    !      alinemom7(time_slice, id,4)=(csum2smom(1,4)-csum2smom(2,4)&
    !     &+csum2smom(3,4)-csum2smom(4,4))*adiv1
    !**********************************************************************
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 8)=(csum2ws(1)+csum2ws(2)+csum2ws(3)+csum2ws(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 8)=(csum2wsmom(1,1)+csum2wsmom(2,1)&
        +csum2wsmom(3,1)+csum2wsmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 8)=(csum2wsmom(1,2)+csum2wsmom(2,2)&
        +csum2wsmom(3,2)+csum2wsmom(4,2))*adiv1
    !      alinemom8(time_slice, id,3)=(csum2wsmom(1,3)+csum2wsmom(2,3)&
    !     &+csum2wsmom(3,3)+csum2wsmom(4,3))*adiv1
    !      alinemom8(time_slice, id,4)=(csum2wsmom(1,4)+csum2wsmom(2,4)&
    !     &+csum2wsmom(3,4)+csum2wsmom(4,4))*adiv1
    !************************************************************************
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 9)=(csum2ws(1)+giot*csum2ws(2)-csum2ws(3)&
        -giot*csum2ws(4))*adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 9)=(csum2wsmom(1,1)+giot*csum2wsmom(2,1)&
        -csum2wsmom(3,1)-giot*csum2wsmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 9)=(csum2wsmom(1,2)+giot*csum2wsmom(2,2)&
        -csum2wsmom(3,2)-giot*csum2wsmom(4,2))*adiv1
    !      alinemom9(time_slice, id,3)=(csum2wsmom(1,3)+giot*csum2wsmom(2,3)&
    !     &-csum2wsmom(3,3)-giot*csum2wsmom(4,3))*adiv1
    !      alinemom9(time_slice, id,4)=(csum2wsmom(1,4)+giot*csum2wsmom(2,4)&
    !     &-csum2wsmom(3,4)-giot*csum2wsmom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 10)=(csum2ws(1)-csum2ws(2)+csum2ws(3)-csum2ws(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 10)=(csum2wsmom(1,1)-csum2wsmom(2,1)&
        +csum2wsmom(3,1)-csum2wsmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 10)=(csum2wsmom(1,2)-csum2wsmom(2,2)&
        +csum2wsmom(3,2)-csum2wsmom(4,2))*adiv1
    !      alinemom10(time_slice, id,3)=(csum2wsmom(1,3)-csum2wsmom(2,3)&
    !     &+csum2wsmom(3,3)-csum2wsmom(4,3))*adiv1
    !      alinemom10(time_slice, id,4)=(csum2wsmom(1,4)-csum2wsmom(2,4)&
    !     &+csum2wsmom(3,4)-csum2wsmom(4,4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 11)=(csumw(1)+csumw(2)+csumw(3)+csumw(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 11)=(csumwmom(1,1)+csumwmom(2,1)&
        +csumwmom(3,1)+csumwmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 11)=(csumwmom(1,2)+csumwmom(2,2)&
        +csumwmom(3,2)+csumwmom(4,2))*adiv1
    !      alinemom11(time_slice, id,3)=(csumwmom(1,3)+csumwmom(2,3)&
    !     &+csumwmom(3,3)+csumwmom(4,3))*adiv1
    !      alinemom11(time_slice, id,4)=(csumwmom(1,4)+csumwmom(2,4)&
    !     &+csumwmom(3,4)+csumwmom(4,4))*adiv1
    !***********************************************************************
    !     j=1, pr=-, q=0
    !***********************************************************************
        lines(time_slice, id, 12)=(csumw(1)+giot*csumw(2)-csumw(3)&
        -giot*csumw(4))*adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 12)=(csumwmom(1,1)+giot*csumwmom(2,1)&
        -csumwmom(3,1)-giot*csumwmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 12)=(csumwmom(1,2)+giot*csumwmom(2,2)&
        -csumwmom(3,2)-giot*csumwmom(4,2))*adiv1
    !      alinemom12(time_slice, id,3)=(csumwmom(1,3)+giot*csumwmom(2,3)&
    !     &-csumwmom(3,3)-giot*csumwmom(4,3))*adiv1
    !      alinemom12(time_slice, id,4)=(csumwmom(1,4)+giot*csumwmom(2,4)&
    !     &-csumwmom(3,4)-giot*csumwmom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 13)=(csumw(1)-csumw(2)+csumw(3)-csumw(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 13)=(csumwmom(1,1)-csumwmom(2,1)&
        +csumwmom(3,1)-csumwmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 13)=(csumwmom(1,2)-csumwmom(2,2)&
        +csumwmom(3,2)-csumwmom(4,2))*adiv1
    !      alinemom13(time_slice, id,3)=(csumwmom(1,3)-csumwmom(2,3)&
    !     &+csumwmom(3,3)-csumwmom(4,3))*adiv1
    !      alinemom13(time_slice, id,4)=(csumwmom(1,4)-csumwmom(2,4)&
    !     &+csumwmom(3,4)-csumwmom(4,4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 14)=(csum2w(1)+csum2w(2)+csum2w(3)+csum2w(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 14)=(csum2wmom(1,1)+csum2wmom(2,1)&
        +csum2wmom(3,1)+csum2wmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 14)=(csum2wmom(1,2)+csum2wmom(2,2)&
        +csum2wmom(3,2)+csum2wmom(4,2))*adiv1
    !      alinemom14(time_slice, id,3)=(csum2wmom(1,3)+csum2wmom(2,3)&
    !     &+csum2wmom(3,3)+csum2wmom(4,3))*adiv1
    !      alinemom14(time_slice, id,4)=(csum2wmom(1,4)+csum2wmom(2,4)&
    !     &+csum2wmom(3,4)+csum2wmom(4,4))*adiv1
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 15)=(csum2w(1)+giot*csum2w(2)-csum2w(3)&
        -giot*csum2w(4))*adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 15)=(csum2wmom(1,1)+giot*csum2wmom(2,1)&
        -csum2wmom(3,1)-giot*csum2wmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 15)=(csum2wmom(1,2)+giot*csum2wmom(2,2)&
        -csum2wmom(3,2)-giot*csum2wmom(4,2))*adiv1
    !      alinemom15(time_slice, id,3)=(csum2wmom(1,3)+giot*csum2wmom(2,3)&
    !     &-csum2wmom(3,3)-giot*csum2wmom(4,3))*adiv1
    !      alinemom15(time_slice, id,4)=(csum2wmom(1,4)+giot*csum2wmom(2,4)&
    !     &-csum2wmom(3,4)-giot*csum2wmom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 16)=(csum2w(1)-csum2w(2)+csum2w(3)-csum2w(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 16)=(csum2wmom(1,1)-csum2wmom(2,1)&
        +csum2wmom(3,1)-csum2wmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 16)=(csum2wmom(1,2)-csum2wmom(2,2)&
        +csum2wmom(3,2)-csum2wmom(4,2))*adiv1
    !      alinemom16(time_slice, id,3)=(csum2wmom(1,3)-csum2wmom(2,3)&
    !     &+csum2wmom(3,3)-csum2wmom(4,3))*adiv1
    !      alinemom16(time_slice, id,4)=(csum2wmom(1,4)-csum2wmom(2,4)&
    !     &+csum2wmom(3,4)-csum2wmom(4,4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 17)=(csum3w(1)+csum3w(2)+csum3w(3)+csum3w(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 17)=(csum3wmom(1,1)+csum3wmom(2,1)&
        +csum3wmom(3,1)+csum3wmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 17)=(csum3wmom(1,2)+csum3wmom(2,2)&
        +csum3wmom(3,2)+csum3wmom(4,2))*adiv1
    !      alinemom17(time_slice, id,3)=(csum3wmom(1,3)+csum3wmom(2,3)&
    !     &+csum3wmom(3,3)+csum3wmom(4,3))*adiv1
    !      alinemom17(time_slice, id,4)=(csum3wmom(1,4)+csum3wmom(2,4)&
    !     &+csum3wmom(3,4)+csum3wmom(4,4))*adiv1
    !***********************************************************************
    !     j=1, pr=+ q=0
    !**********************************************************************
        lines(time_slice, id, 18)=(csum3w(1)+giot*csum3w(2)-csum3w(3)-giot*csum3w(4))&
        *adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 18)=(csum3wmom(1,1)+giot*csum3wmom(2,1)&
        -csum3wmom(3,1)-giot*csum3wmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 18)=(csum3wmom(1,2)+giot*csum3wmom(2,2)&
        -csum3wmom(3,2)-giot*csum3wmom(4,2))*adiv1
    !      alinemom18(time_slice, id,3)=(csum3wmom(1,3)+giot*csum3wmom(2,3)&
    !     &-csum3wmom(3,3)-giot*csum3wmom(4,3))*adiv1
    !      alinemom18(time_slice, id,4)=(csum3wmom(1,4)+giot*csum3wmom(2,4)&
    !     &-csum3wmom(3,4)-giot*csum3wmom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 19)=(csum3w(1)-csum3w(2)+csum3w(3)-csum3w(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 19)=(csum3wmom(1,1)-csum3wmom(2,1)&
        +csum3wmom(3,1)-csum3wmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 19)=(csum3wmom(1,2)-csum3wmom(2,2)&
        +csum3wmom(3,2)-csum3wmom(4,2))*adiv1
    !      alinemom19(time_slice, id,3)=(csum3wmom(1,3)-csum3wmom(2,3)&
    !     &+csum3wmom(3,3)-csum3wmom(4,3))*adiv1
    !      alinemom19(time_slice, id,4)=(csum3wmom(1,4)-csum3wmom(2,4)&
    !     &+csum3wmom(3,4)-csum3wmom(4,4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 20)=(csumup(1)+csumup(2)+csumup(3)+csumup(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 20)=(csumupmom(1,1)+csumupmom(2,1)&
        +csumupmom(3,1)+csumupmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 20)=(csumupmom(1,2)+csumupmom(2,2)&
        +csumupmom(3,2)+csumupmom(4,2))*adiv1
    !      alinemom20(time_slice, id,3)=(csumupmom(1,3)+csumupmom(2,3)&
    !     &+csumupmom(3,3)+csumupmom(4,3))*adiv1
    !      alinemom20(time_slice, id,4)=(csumupmom(1,4)+csumupmom(2,4)&
    !     &+csumupmom(3,4)+csumupmom(4,4))*adiv1
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 21)=(csumup(1)+giot*csumup(2)-csumup(3)-giot*csumup(4))&
        *adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 21)=(csumupmom(1,1)+giot*csumupmom(2,1)&
        -csumupmom(3,1)-giot*csumupmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 21)=(csumupmom(1,2)+giot*csumupmom(2,2)&
        -csumupmom(3,2)-giot*csumupmom(4,2))*adiv1
    !      alinemom21(time_slice, id,3)=(csumupmom(1,3)+giot*csumupmom(2,3)&
    !     &-csumupmom(3,3)-giot*csumupmom(4,3))*adiv1
    !      alinemom21(time_slice, id,4)=(csumupmom(1,4)+giot*csumupmom(2,4)&
    !     &-csumupmom(3,4)-giot*csumupmom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 22)=(csumup(1)-csumup(2)+csumup(3)-csumup(4))*adiv1
    !***********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 22)=(csumupmom(1,1)-csumupmom(2,1)&
        +csumupmom(3,1)-csumupmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 22)=(csumupmom(1,2)-csumupmom(2,2)&
        +csumupmom(3,2)-csumupmom(4,2))*adiv1
    !      alinemom22(time_slice, id,3)=(csumupmom(1,3)-csumupmom(2,3)&
    !     &+csumupmom(3,3)-csumupmom(4,3))*adiv1
    !      alinemom22(time_slice, id,4)=(csumupmom(1,4)-csumupmom(2,4)&
    !     &+csumupmom(3,4)-csumupmom(4,4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 23)=(csumud(1)+csumud(2)+csumud(3)+csumud(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 23)=(csumudmom(1,1)+csumudmom(2,1)&
        +csumudmom(3,1)+csumudmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 23)=(csumudmom(1,2)+csumudmom(2,2)&
        +csumudmom(3,2)+csumudmom(4,2))*adiv1
    !      alinemom23(time_slice, id,3)=(csumudmom(1,3)+csumudmom(2,3)&
    !     &+csumudmom(3,3)+csumudmom(4,3))*adiv1
    !      alinemom23(time_slice, id,4)=(csumudmom(1,4)+csumudmom(2,4)&
    !     &+csumudmom(3,4)+csumudmom(4,4))*adiv1
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 24)=(csumud(1)+giot*csumud(2)-csumud(3)&
        -giot*csumud(4))*adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !***********************************************************************
        momentum_lines(time_slice, id, 1, 24)=(csumudmom(1,1)+giot*csumudmom(2,1)&
        -csumudmom(3,1)-giot*csumudmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 24)=(csumudmom(1,2)+giot*csumudmom(2,2)&
        -csumudmom(3,2)-giot*csumudmom(4,2))*adiv1
    !      alinemom24(time_slice, id,3)=(csumudmom(1,3)+giot*csumudmom(2,3)&
    !     &-csumudmom(3,3)-giot*csumudmom(4,3))*adiv1
    !      alinemom24(time_slice, id,4)=(csumudmom(1,4)+giot*csumudmom(2,4)&
    !     &-csumudmom(3,4)-giot*csumudmom(4,4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 25)=(csumud(1)-csumud(2)+csumud(3)-csumud(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        momentum_lines(time_slice, id, 1, 25)=(csumudmom(1,1)-csumudmom(2,1)&
        +csumudmom(3,1)-csumudmom(4,1))*adiv1
        momentum_lines(time_slice, id, 2, 25)=(csumudmom(1,2)-csumudmom(2,2)&
        +csumudmom(3,2)-csumudmom(4,2))*adiv1
    !      alinemom25(time_slice, id,3)=(csumudmom(1,3)-csumudmom(2,3)&
    !     &+csumudmom(3,3)-csumudmom(4,3))*adiv1
    !      alinemom25(time_slice, id,4)=(csumudmom(1,4)-csumudmom(2,4)&
    !     &+csumudmom(3,4)-csumudmom(4,4))*adiv1
    !**********************************************************************
    !**********************************************************************
    !**********************************************************************
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 26)=(csumtt1(1)+csumtt1(2)+csumtt1(3)+csumtt1(4)+&
        (csumtt1(7)+csumtt1(6)+csumtt1(5)+csumtt1(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 26)=(csumttmom1(1,ik)+csumttmom1(2,ik)&
        +csumttmom1(3,ik)+csumttmom1(4,ik)+csumttmom1(7,ik)&
        +csumttmom1(6,ik)+csumttmom1(5,ik)&
        +csumttmom1(8,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 27)=(csumtt1(1)+csumtt1(2)+csumtt1(3)+csumtt1(4)-&
        (csumtt1(7)+csumtt1(6)+csumtt1(5)+csumtt1(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 27)=(csumttmom1(1,ik)+csumttmom1(2,ik)&
        +csumttmom1(3,ik)+csumttmom1(4,ik)-(csumttmom1(7,ik)&
        +csumttmom1(6,ik)+csumttmom1(5,ik)&
        +csumttmom1(8,ik)))*adiv2
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 28)=(csumtt1(1)+giot*csumtt1(2)-csumtt1(3)&
        -giot*csumtt1(4)+(csumtt1(6)+giot*csumtt1(7)-csumtt1(8)&
        -giot*csumtt1(5)))*adiv2
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 28)=(csumttmom1(1,ik)+giot*csumttmom1(2,ik)&
        -csumttmom1(3,ik)-giot*csumttmom1(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 29)=(csumtt1(1)+giot*csumtt1(2)-csumtt1(3)&
        -giot*csumtt1(4)-(csumtt1(6)+giot*csumtt1(7)-csumtt1(8)&
        -giot*csumtt1(5)))*adiv2
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 29)=(csumttmom1(6,ik)+giot*csumttmom1(7,ik)&
        -csumttmom1(8,ik)-giot*csumttmom1(5,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 30)=(csumtt1(1)-csumtt1(2)+csumtt1(3)-csumtt1(4)+&
        csumtt1(7)-csumtt1(6)+csumtt1(5)-csumtt1(8))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 30)=(csumttmom1(1,ik)-csumttmom1(2,ik)&
        +csumttmom1(3,ik)-csumttmom1(4,ik)+csumttmom1(7,ik)&
        -csumttmom1(6,ik)+csumttmom1(5,ik)-csumttmom1(8,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 31)=(csumtt1(1)-csumtt1(2)+csumtt1(3)-csumtt1(4)-&
        (csumtt1(7)-csumtt1(6)+csumtt1(5)-csumtt1(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 31)=(csumttmom1(1,ik)-csumttmom1(2,ik)&
        +csumttmom1(3,ik)-csumttmom1(4,ik)-(csumttmom1(7,ik)&
        -csumttmom1(6,ik)+csumttmom1(5,ik)-csumttmom1(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 32)=(csumtt2(1)+csumtt2(2)+csumtt2(3)+csumtt2(4)+&
        csumtt2(5)+csumtt2(6)+csumtt2(7)+csumtt2(8))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 32)=(csumttmom2(1,ik)+csumttmom2(2,ik)&
        +csumttmom2(3,ik)+csumttmom2(4,ik)+csumttmom2(5,ik)&
        +csumttmom2(6,ik)+csumttmom2(7,ik)&
        +csumttmom2(8,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 33)=(csumtt2(1)+csumtt2(2)+csumtt2(3)+csumtt2(4)-&
        (csumtt2(5)+csumtt2(6)+csumtt2(7)+csumtt2(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 33)=(csumttmom2(1,ik)+csumttmom2(2,ik)&
        +csumttmom2(3,ik)+csumttmom2(4,ik)-(csumttmom2(5,ik)&
        +csumttmom2(6,ik)+csumttmom2(7,ik)&
        +csumttmom2(8,ik)))*adiv2
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 34)=(csumtt2(1)+giot*csumtt2(2)-csumtt2(3)&
        -giot*csumtt2(4)+(csumtt2(7)+giot*csumtt2(8)-csumtt2(5)&
        -giot*csumtt2(6)))*adiv2
    !***********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 34)=(csumttmom2(1,ik)+giot*csumttmom2(2,ik)&
        -csumttmom2(3,ik)-giot*csumttmom2(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 35)=(csumtt2(1)+giot*csumtt2(2)-csumtt2(3)&
        -giot*csumtt2(4)-(csumtt2(7)+giot*csumtt2(8)-csumtt2(5)&
        -giot*csumtt2(6)))*adiv2
    !***********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 35)=(csumttmom2(7,ik)+giot*csumttmom2(8,ik)&
        -csumttmom2(5,ik)-giot*csumttmom2(6,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 36)=(csumtt2(1)-csumtt2(2)+csumtt2(3)-csumtt2(4)+&
        csumtt2(6)-csumtt2(5)+csumtt2(8)-csumtt2(7))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 36)=(csumttmom2(1,ik)-csumttmom2(2,ik)&
        +csumttmom2(3,ik)-csumttmom2(4,ik)+csumttmom2(6,ik)&
        -csumttmom2(5,ik)+csumttmom2(8,ik)&
        -csumttmom2(7,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 37)=(csumtt2(1)-csumtt2(2)+csumtt2(3)-csumtt2(4)-&
        (csumtt2(6)-csumtt2(5)+csumtt2(8)-csumtt2(7)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 37)=(csumttmom2(1,ik)-csumttmom2(2,ik)&
        +csumttmom2(3,ik)-csumttmom2(4,ik)-(csumttmom2(6,ik)&
        -csumttmom2(5,ik)+csumttmom2(8,ik)&
        -csumttmom2(7,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 38)=(csumtt3(1)+csumtt3(2)+csumtt3(3)+csumtt3(4)+&
        csumtt3(5)+csumtt3(6)+csumtt3(7)+csumtt3(8))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 38)=(csumttmom3(1,ik)+csumttmom3(2,ik)&
        +csumttmom3(3,ik)+csumttmom3(4,ik)+csumttmom3(5,ik)&
        +csumttmom3(6,ik)+csumttmom3(7,ik)&
        +csumttmom3(8,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 39)=(csumtt3(1)+csumtt3(2)+csumtt3(3)+csumtt3(4)-&
        (csumtt3(5)+csumtt3(6)+csumtt3(7)+csumtt3(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 39)=(csumttmom3(1,ik)+csumttmom3(2,ik)&
        +csumttmom3(3,ik)+csumttmom3(4,ik)-(csumttmom3(5,ik)&
        +csumttmom3(6,ik)+csumttmom3(7,ik)&
        +csumttmom3(8,ik)))*adiv2
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 40)=(csumtt3(1)+giot*csumtt3(2)-csumtt3(3)&
        -giot*csumtt3(4)+(csumtt3(6)+giot*csumtt3(7)-csumtt3(8)&
        -giot*csumtt3(5)))*adiv2
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 40)=(csumttmom3(1,ik)+giot*csumttmom3(2,ik)&
        -csumttmom3(3,ik)-giot*csumttmom3(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 41)=(csumtt3(1)+giot*csumtt3(2)-csumtt3(3)&
        -giot*csumtt3(4)-(csumtt3(6)+giot*csumtt3(7)-csumtt3(8)&
        -giot*csumtt3(5)))*adiv2
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 41)=(csumttmom3(6,ik)+giot*csumttmom3(7,ik)&
        -csumttmom3(8,ik)-giot*csumttmom3(5,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 42)=(csumtt3(1)-csumtt3(2)+csumtt3(3)-csumtt3(4)+&
        (csumtt3(7)-csumtt3(6)+csumtt3(5)-csumtt3(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 42)=(csumttmom3(1,ik)-csumttmom3(2,ik)&
        +csumttmom3(3,ik)-csumttmom3(4,ik)+(csumttmom3(7,ik)&
        -csumttmom3(6,ik)+csumttmom3(5,ik)-csumttmom3(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 43)=(csumtt3(1)-csumtt3(2)+csumtt3(3)-csumtt3(4)-&
        (csumtt3(7)-csumtt3(6)+csumtt3(5)-csumtt3(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 43)=(csumttmom3(1,ik)-csumttmom3(2,ik)&
        +csumttmom3(3,ik)-csumttmom3(4,ik)-(csumttmom3(7,ik)&
        -csumttmom3(6,ik)+csumttmom3(5,ik)-csumttmom3(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 44)=(csumtt4(1)+csumtt4(2)+csumtt4(3)+csumtt4(4)+&
        csumtt4(5)+csumtt4(6)+csumtt4(7)+csumtt4(8))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 44)=(csumttmom4(1,ik)+csumttmom4(2,ik)&
        +csumttmom4(3,ik)+csumttmom4(4,ik)+csumttmom4(5,ik)&
        +csumttmom4(6,ik)+csumttmom4(7,ik)&
        +csumttmom4(8,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 45)=(csumtt4(1)+csumtt4(2)+csumtt4(3)+csumtt4(4)-&
        (csumtt4(5)+csumtt4(6)+csumtt4(7)+csumtt4(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 45)=(csumttmom4(1,ik)+csumttmom4(2,ik)&
        +csumttmom4(3,ik)+csumttmom4(4,ik)-(csumttmom4(5,ik)&
        +csumttmom4(6,ik)+csumttmom4(7,ik)&
        +csumttmom4(8,ik)))*adiv2
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 46)=(csumtt4(1)+giot*csumtt4(2)-csumtt4(3)&
        -giot*csumtt4(4)+(csumtt4(7)+giot*csumtt4(8)-csumtt4(5)&
        -giot*csumtt4(6)))*adiv2
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 46)=(csumttmom4(1,ik)+giot*csumttmom4(2,ik)&
        -csumttmom4(3,ik)-giot*csumttmom4(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 47)=(csumtt4(1)+giot*csumtt4(2)-csumtt4(3)&
        -giot*csumtt4(4)-(csumtt4(7)+giot*csumtt4(8)-csumtt4(5)&
        -giot*csumtt4(6)))*adiv2
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 47)=(csumttmom4(7,ik)+giot*csumttmom4(8,ik)&
        -csumttmom4(5,ik)-giot*csumttmom4(6,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 48)=(csumtt4(1)-csumtt4(2)+csumtt4(3)-csumtt4(4)+&
        csumtt4(8)-csumtt4(7)+csumtt4(6)-csumtt4(5))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 48)=(csumttmom4(1,ik)-csumttmom4(2,ik)&
        +csumttmom4(3,ik)-csumttmom4(4,ik)+csumttmom4(8,ik)&
        -csumttmom4(7,ik)+csumttmom4(6,ik)-csumttmom4(5,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 49)=(csumtt4(1)-csumtt4(2)+csumtt4(3)-csumtt4(4)-&
        (csumtt4(8)-csumtt4(7)+csumtt4(6)-csumtt4(5)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 49)=(csumttmom4(1,ik)-csumttmom4(2,ik)&
        +csumttmom4(3,ik)-csumttmom4(4,ik)-(csumttmom4(8,ik)&
        -csumttmom4(7,ik)+csumttmom4(6,ik)-csumttmom4(5,ik)))*adiv2
        enddo
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-5 operators
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 50)=(csumtt5(1)+csumtt5(2)+csumtt5(3)+csumtt5(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 50)=(csumttmom5(1,ik)+csumttmom5(2,ik)&
        +csumttmom5(3,ik)+csumttmom5(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 51)=(csumtt5(1)+giot*csumtt5(2)-csumtt5(3)&
        -giot*csumtt5(4))*adiv1
    !***********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 51)=(csumttmom5(1,ik)+giot*csumttmom5(2,ik)&
        -csumttmom5(3,ik)-giot*csumttmom5(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 52)=(csumtt5(1)-csumtt5(2)+csumtt5(3)-csumtt5(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 52)=(csumttmom5(1,ik)-csumttmom5(2,ik)&
        +csumttmom5(3,ik)-csumttmom5(4,ik))*adiv1
        enddo
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! 6
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-6 operators
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 53)=(csumtt6(1)+csumtt6(2)+csumtt6(3)+csumtt6(4))*adiv1
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 53)=(csumttmom6(1,ik)+csumttmom6(2,ik)&
        +csumttmom6(3,ik)+csumttmom6(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 54)=(csumtt6(1)+giot*csumtt6(2)-csumtt6(3)&
        -giot*csumtt6(4))*adiv1
    !***********************************************************************
    !     j=1, q
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 54)=(csumttmom6(1,ik)+giot*csumttmom6(2,ik)&
        -csumttmom6(3,ik)-giot*csumttmom6(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 55)=(csumtt6(1)-csumtt6(2)+csumtt6(3)-csumtt6(4))*adiv1
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 55)=(csumttmom6(1,ik)-csumttmom6(2,ik)&
        +csumttmom6(3,ik)-csumttmom6(4,ik))*adiv1
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! 7
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-7 operators
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 56)=(csumtt7(1)+csumtt7(2)+csumtt7(3)+csumtt7(4)+&
        (csumtt7(5)+csumtt7(6)+csumtt7(7)+csumtt7(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 56)=(csumttmom7(1,ik)+csumttmom7(2,ik)&
        +csumttmom7(3,ik)+csumttmom7(4,ik)+(csumttmom7(5,ik)&
        +csumttmom7(6,ik)+csumttmom7(7,ik)&
        +csumttmom7(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 57)=(csumtt7(1)+csumtt7(2)+csumtt7(3)+csumtt7(4)-&
        (csumtt7(5)+csumtt7(6)+csumtt7(7)+csumtt7(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 57)=(csumttmom7(1,ik)+csumttmom7(2,ik)&
        +csumttmom7(3,ik)+csumttmom7(4,ik)-(csumttmom7(5,ik)&
        +csumttmom7(6,ik)+csumttmom7(7,ik)&
        +csumttmom7(8,ik)))*adiv2
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 58)=(csumtt7(1)+giot*csumtt7(2)-csumtt7(3)&
        -giot*csumtt7(4)+&
        (csumtt7(7)+giot*csumtt7(8)-csumtt7(5)-giot*csumtt7(6)))*adiv2
    !***********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 58)=(csumttmom7(1,ik)+giot*csumttmom7(2,ik)&
        -csumttmom7(3,ik)-giot*csumttmom7(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 59)=(csumtt7(1)+giot*csumtt7(2)-csumtt7(3)&
        -giot*csumtt7(4)-&
        (csumtt7(7)+giot*csumtt7(8)-csumtt7(5)-giot*csumtt7(6)))*adiv2
    !***********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 59)=(csumttmom7(7,ik)+giot*csumttmom7(8,ik)&
        -csumttmom7(5,ik)-giot*csumttmom7(6,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 60)=(csumtt7(1)-csumtt7(2)+csumtt7(3)-csumtt7(4)+&
        (csumtt7(8)-csumtt7(7)+csumtt7(6)-csumtt7(5)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 60)=(csumttmom7(1,ik)-csumttmom7(2,ik)&
        +csumttmom7(3,ik)-csumttmom7(4,ik)+(csumttmom7(8,ik)&
        -csumttmom7(7,ik)+csumttmom7(6,ik)&
        -csumttmom7(5,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 61)=(csumtt7(1)-csumtt7(2)+csumtt7(3)-csumtt7(4)-&
        (csumtt7(8)-csumtt7(7)+csumtt7(6)-csumtt7(5)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 61)=(csumttmom7(1,ik)-csumttmom7(2,ik)&
        +csumttmom7(3,ik)-csumttmom7(4,ik)-(csumttmom7(8,ik)&
        -csumttmom7(7,ik)+csumttmom7(6,ik)&
        -csumttmom7(5,ik)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! 8
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-8 operators cp=+,pr=+
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 62)=(csumtt8(1)+csumtt8(2)+csumtt8(3)+csumtt8(4)&
        +(csumtt8(5)+csumtt8(6)+csumtt8(7)+csumtt8(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 62)=(csumttmom8(1,ik)+csumttmom8(2,ik)&
        +csumttmom8(3,ik)+csumttmom8(4,ik)+(csumttmom8(5,ik)&
        +csumttmom8(6,ik)+csumttmom8(7,ik)&
        +csumttmom8(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 63)=(csumtt8(1)+csumtt8(2)+csumtt8(3)+csumtt8(4)&
        -(csumtt8(5)+csumtt8(6)+csumtt8(7)+csumtt8(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 63)=(csumttmom8(1,ik)+csumttmom8(2,ik)&
        +csumttmom8(3,ik)+csumttmom8(4,ik)-(csumttmom8(5,ik)&
        +csumttmom8(6,ik)+csumttmom8(7,ik)&
        +csumttmom8(8,ik)))*adiv2
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 64)=(csumtt8(1)+giot*csumtt8(2)-csumtt8(3)&
        -giot*csumtt8(4))*adiv1
    !***********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 64)=(csumttmom8(1,ik)+giot*csumttmom8(2,ik)&
        -csumttmom8(3,ik)-giot*csumttmom8(4,ik))*adiv1
        enddo
    !***********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 65)=(csumtt8(7)+giot*csumtt8(6)-csumtt8(5)&
        -giot*csumtt8(8))*adiv1
    !***********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 65)=(csumttmom8(7,ik)+giot*csumttmom8(6,ik)&
        -csumttmom8(5,ik)-giot*csumttmom8(8,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 66)=(csumtt8(1)-csumtt8(2)+csumtt8(3)-csumtt8(4)&
        +(csumtt8(5)-csumtt8(6)+csumtt8(7)-csumtt8(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 66)=(csumttmom8(1,ik)-csumttmom8(2,ik)&
        +csumttmom8(3,ik)-csumttmom8(4,ik)+(csumttmom8(5,ik)&
        -csumttmom8(6,ik)+csumttmom8(7,ik)&
        -csumttmom8(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 67)=(csumtt8(1)-csumtt8(2)+csumtt8(3)-csumtt8(4)&
        -(csumtt8(5)-csumtt8(6)+csumtt8(7)-csumtt8(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 67)=(csumttmom8(1,ik)-csumttmom8(2,ik)&
        +csumttmom8(3,ik)-csumttmom8(4,ik)-(csumttmom8(5,ik)&
        -csumttmom8(6,ik)+csumttmom8(7,ik)&
        -csumttmom8(8,ik)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! 9
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-9 operators
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 68)=(csumtt9(1)+csumtt9(2)+csumtt9(3)+csumtt9(4)+&
        (csumtt9(5)+csumtt9(6)+csumtt9(7)+csumtt9(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 68)=(csumttmom9(1,ik)+csumttmom9(2,ik)+&
        csumttmom9(3,ik)+csumttmom9(4,ik)&
        +(csumttmom9(5,ik)+csumttmom9(6,ik)+csumttmom9(7,ik)&
        +csumttmom9(8,ik)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-9 operators cp=-,pz=-,j=0
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 69)=(csumtt9(1)+csumtt9(2)+csumtt9(3)+csumtt9(4)&
        -(csumtt9(5)+csumtt9(6)+csumtt9(7)+csumtt9(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 69)=(csumttmom9(1,ik)+csumttmom9(2,ik)+&
        csumttmom9(3,ik)+csumttmom9(4,ik)&
        -(csumttmom9(5,ik)+csumttmom9(6,ik)+csumttmom9(7,ik)&
        +csumttmom9(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 70)=(csumtt9(1)+giot*csumtt9(2)-csumtt9(3)&
        -giot*csumtt9(4)+(csumtt9(7)+giot*csumtt9(8)-csumtt9(5)&
        -giot*csumtt9(6)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 70)=(csumttmom9(1,ik)+giot*csumttmom9(2,ik)&
        -csumttmom9(3,ik)-giot*csumttmom9(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 71)=(csumtt9(1)+giot*csumtt9(2)-csumtt9(3)&
        -giot*csumtt9(4)-(csumtt9(7)+giot*csumtt9(8)-csumtt9(5)&
        -giot*csumtt9(6)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 71)=(csumttmom9(7,ik)+giot*csumttmom9(8,ik)&
        -csumttmom9(5,ik)-giot*csumttmom9(6,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 72)=(csumtt9(1)-csumtt9(2)+csumtt9(3)-csumtt9(4)&
        +(csumtt9(7)-csumtt9(6)+csumtt9(5)-csumtt9(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 72)=(csumttmom9(1,ik)-csumttmom9(2,ik)&
        +csumttmom9(3,ik)-csumttmom9(4,ik)&
        +(csumttmom9(7,ik)-csumttmom9(6,ik)+csumttmom9(5,ik)&
        -csumttmom9(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0 !here!
    !**********************************************************************
        lines(time_slice, id, 73)=(csumtt9(1)-csumtt9(2)+csumtt9(3)-csumtt9(4)&
        -(csumtt9(7)-csumtt9(6)+csumtt9(5)-csumtt9(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 73)=(csumttmom9(1,ik)-csumttmom9(2,ik)&
        +csumttmom9(3,ik)-csumttmom9(4,ik)&
        -(csumttmom9(7,ik)-csumttmom9(6,ik)+csumttmom9(5,ik)&
        -csumttmom9(8,ik)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-10 operators
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 74)=(csumtt10(1)+csumtt10(2)+csumtt10(3)+csumtt10(4)&
        +(csumtt10(5)+csumtt10(6)+csumtt10(7)+csumtt10(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 74)=(csumttmom10(1,ik)+csumttmom10(2,ik)&
        +csumttmom10(3,ik)+csumttmom10(4,ik)&
        +(csumttmom10(5,ik)+csumttmom10(6,ik)+csumttmom10(7,ik)&
        +csumttmom10(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 75)=(csumtt10(1)+csumtt10(2)+csumtt10(3)+csumtt10(4)&
        -(csumtt10(5)+csumtt10(6)+csumtt10(7)+csumtt10(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 75)=(csumttmom10(1,ik)+csumttmom10(2,ik)&
        +csumttmom10(3,ik)+csumttmom10(4,ik)&
        -(csumttmom10(5,ik)+csumttmom10(6,ik)+csumttmom10(7,ik)&
        +csumttmom10(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 76)=(csumtt10(1)+giot*csumtt10(2)-csumtt10(3)&
        -giot*csumtt10(4)+(csumtt10(7)+giot*csumtt10(6)-csumtt10(5)&
        -giot*csumtt10(8)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 76)=(csumttmom10(1,ik)+giot*csumttmom10(2,ik)&
        -csumttmom10(3,ik)-giot*csumttmom10(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 77)=(csumtt10(1)+giot*csumtt10(2)-csumtt10(3)&
        -giot*csumtt10(4)-(csumtt10(7)+giot*csumtt10(6)-csumtt10(5)&
        -giot*csumtt10(8)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 77)=(csumttmom10(7,ik)+giot*csumttmom10(6,ik)&
        -csumttmom10(5,ik)-giot*csumttmom10(8,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 78)=(csumtt10(1)-csumtt10(2)+csumtt10(3)-csumtt10(4)&
        +(csumtt10(5)-csumtt10(6)+csumtt10(7)-csumtt10(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 78)=(csumttmom10(1,ik)-csumttmom10(2,ik)&
        +csumttmom10(3,ik)-csumttmom10(4,ik)&
        +(csumttmom10(5,ik)-csumttmom10(6,ik)+csumttmom10(7,ik)&
        -csumttmom10(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 79)=(csumtt10(1)-csumtt10(2)+csumtt10(3)-csumtt10(4)&
        -(csumtt10(5)-csumtt10(6)+csumtt10(7)-csumtt10(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 79)=(csumttmom10(1,ik)-csumttmom10(2,ik)&
        +csumttmom10(3,ik)-csumttmom10(4,ik)&
        -(csumttmom10(5,ik)-csumttmom10(6,ik)+csumttmom10(7,ik)&
        -csumttmom10(8,ik)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-11 operators
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 80)=(csumtt11(1)+csumtt11(2)+csumtt11(3)+csumtt11(4)&
        +(csumtt11(5)+csumtt11(6)+csumtt11(7)+csumtt11(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 80)=(csumttmom11(1,ik)+csumttmom11(2,ik)&
        +csumttmom11(3,ik)+csumttmom11(4,ik)&
        +(csumttmom11(5,ik)+csumttmom11(6,ik)+csumttmom11(7,ik)&
        +csumttmom11(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 81)=(csumtt11(1)+csumtt11(2)+csumtt11(3)+csumtt11(4)&
        -(csumtt11(5)+csumtt11(6)+csumtt11(7)+csumtt11(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 81)=(csumttmom11(1,ik)+csumttmom11(2,ik)&
        +csumttmom11(3,ik)+csumttmom11(4,ik)&
        -(csumttmom11(5,ik)+csumttmom11(6,ik)+csumttmom11(7,ik)&
        +csumttmom11(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 82)=(csumtt11(1)+giot*csumtt11(2)-csumtt11(3)&
        -giot*csumtt11(4))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 82)=(csumttmom11(1,ik)+giot*csumttmom11(2,ik)&
        -csumttmom11(3,ik)-giot*csumttmom11(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 83)=(csumtt11(7)+giot*csumtt11(6)-csumtt11(5)&
        -giot*csumtt11(8))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 83)=(csumttmom11(7,ik)+giot*csumttmom11(6,ik)&
        -csumttmom11(5,ik)-giot*csumttmom11(8,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 84)=(csumtt11(1)-csumtt11(2)+csumtt11(3)-csumtt11(4)&
        +(csumtt11(5)-csumtt11(6)+csumtt11(7)-csumtt11(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 84)=(csumttmom11(1,ik)-csumttmom11(2,ik)&
        +csumttmom11(3,ik)-csumttmom11(4,ik)&
        +(csumttmom11(5,ik)-csumttmom11(6,ik)+csumttmom11(7,ik)&
        -csumttmom11(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 85)=(csumtt11(1)-csumtt11(2)+csumtt11(3)-csumtt11(4)&
        -(csumtt11(5)-csumtt11(6)+csumtt11(7)-csumtt11(8)))*adiv2
    !c**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 85)=(csumttmom11(1,ik)-csumttmom11(2,ik)&
        +csumttmom11(3,ik)-csumttmom11(4,ik)&
        -(csumttmom11(5,ik)-csumttmom11(6,ik)+csumttmom11(7,ik)&
        -csumttmom11(8,ik)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     12 this!!!!!
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-12 operators
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 86)=(csumtt12(1)+csumtt12(2)+csumtt12(3)+csumtt12(4)&
        +(csumtt12(5)+csumtt12(6)+csumtt12(7)+csumtt12(8))&
        +(csumtt12(9)+csumtt12(10)+csumtt12(11)+csumtt12(12))&
        +(csumtt12(13)+csumtt12(14)+csumtt12(15)+csumtt12(16)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 86)=(csumttmom12(1,ik)+csumttmom12(2,ik)&
        +csumttmom12(3,ik)+csumttmom12(4,ik)&
        +(csumttmom12(9,ik)+csumttmom12(10,ik)+csumttmom12(11,ik)&
        +csumttmom12(12,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 87)=(csumtt12(1)+csumtt12(2)+csumtt12(3)+csumtt12(4)&
        +(csumtt12(9)+csumtt12(10)+csumtt12(11)+csumtt12(12))&
        -(csumtt12(7)+csumtt12(8)+csumtt12(5)+csumtt12(6))&
        -(csumtt12(15)+csumtt12(16)+csumtt12(13)+csumtt12(14)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 87)=(csumttmom12(7,ik)+csumttmom12(8,ik)&
        +csumttmom12(5,ik)+csumttmom12(6,ik)+(csumttmom12(15,ik)&
        +csumttmom12(16,ik)+csumttmom12(13,ik)+csumttmom12(14,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 88)=(csumtt12(1)+csumtt12(2)+csumtt12(3)+csumtt12(4)&
        -(csumtt12(9)+csumtt12(10)+csumtt12(11)+csumtt12(12))&
        +(csumtt12(7)+csumtt12(8)+csumtt12(5)+csumtt12(6))&
        -(csumtt12(15)+csumtt12(16)+csumtt12(13)+csumtt12(14)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 88)=(csumttmom12(1,ik)+csumttmom12(2,ik)&
        +csumttmom12(3,ik)+csumttmom12(4,ik)-(csumttmom12(9,ik)&
        +csumttmom12(10,ik)+csumttmom12(11,ik)+csumttmom12(12,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 89)=(csumtt12(1)+csumtt12(2)+csumtt12(3)+csumtt12(4)&
        -(csumtt12(9)+csumtt12(10)+csumtt12(11)+csumtt12(12))&
        -(csumtt12(7)+csumtt12(8)+csumtt12(5)+csumtt12(6))&
        +(csumtt12(15)+csumtt12(16)+csumtt12(13)+csumtt12(14)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 89)=(csumttmom12(7,ik)+csumttmom12(8,ik)&
        +csumttmom12(5,ik)+csumttmom12(6,ik)-(csumttmom12(15,ik)&
        +csumttmom12(16,ik)+csumttmom12(13,ik)+csumttmom12(14,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 90)=(csumtt12(1)+giot*csumtt12(2)-csumtt12(3)&
        -giot*csumtt12(4)&
        +(csumtt12(7)+giot*csumtt12(8)-csumtt12(5)-giot*csumtt12(6))&
        )*adiv2
    !**********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 90)=(csumttmom12(1,ik)+giot*csumttmom12(2,ik)&
        -csumttmom12(3,ik)-giot*csumttmom12(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 91)=(csumtt12(9)+giot*csumtt12(10)-csumtt12(11)&
        -giot*csumtt12(12)&
        +(csumtt12(15)+giot*csumtt12(16)-csumtt12(13)-giot*csumtt12(14))&
        )*adiv2
    !**********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 91)=(csumttmom12(9,ik)+giot*csumttmom12(10,ik)&
        -csumttmom12(11,ik)-giot*csumttmom12(12,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 92)=(csumtt12(1)+giot*csumtt12(2)-csumtt12(3)&
        -giot*csumtt12(4)&
        -(csumtt12(7)+giot*csumtt12(8)-csumtt12(5)-giot*csumtt12(6))&
        )*adiv2
    !**********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 92)=(csumttmom12(7,ik)+giot*csumttmom12(8,ik)&
        -csumttmom12(5,ik)-giot*csumttmom12(6,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 93)=(csumtt12(9)+giot*csumtt12(10)-csumtt12(11)&
        -giot*csumtt12(12)&
        -(csumtt12(15)+giot*csumtt12(16)-csumtt12(13)-giot*csumtt12(14))&
        )*adiv2
    !**********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 93)=(csumttmom12(15,ik)&
        +giot*csumttmom12(16,ik)&
        -csumttmom12(13,ik)-giot*csumttmom12(14,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 94)=(csumtt12(1)-csumtt12(2)+csumtt12(3)-csumtt12(4)&
        +(csumtt12(9)-csumtt12(10)+csumtt12(11)-csumtt12(12))&
        +(csumtt12(7)-csumtt12(8)+csumtt12(5)-csumtt12(6))&
        +(csumtt12(15)-csumtt12(16)+csumtt12(13)-csumtt12(14)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 94)=(csumttmom12(1,ik)-csumttmom12(2,ik)&
        +csumttmom12(3,ik)-csumttmom12(4,ik)&
        +(csumttmom12(9,ik)-csumttmom12(10,ik)+csumttmom12(11,ik)&
        -csumttmom12(12,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 95)=(csumtt12(1)-csumtt12(2)+csumtt12(3)-csumtt12(4)&
        +(csumtt12(9)-csumtt12(10)+csumtt12(11)-csumtt12(12))&
        -(csumtt12(7)-csumtt12(8)+csumtt12(5)-csumtt12(6))&
        -(csumtt12(15)-csumtt12(16)+csumtt12(13)-csumtt12(14)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 95)=(csumttmom12(7,ik)-csumttmom12(8,ik)&
        +csumttmom12(5,ik)-csumttmom12(6,ik)+(csumttmom12(15,ik)&
        -csumttmom12(16,ik)+csumttmom12(13,ik)-csumttmom12(14,ik))&
        )*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 96)=(csumtt12(1)-csumtt12(2)+csumtt12(3)-csumtt12(4)&
        -(csumtt12(9)-csumtt12(10)+csumtt12(11)-csumtt12(12))&
        +(csumtt12(7)-csumtt12(8)+csumtt12(5)-csumtt12(6))&
        -(csumtt12(15)-csumtt12(16)+csumtt12(13)-csumtt12(14)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 96)=(csumttmom12(1,ik)-csumttmom12(2,ik)&
        +csumttmom12(3,ik)-csumttmom12(4,ik)-(csumttmom12(9,ik)&
        -csumttmom12(10,ik)+csumttmom12(11,ik)-csumttmom12(12,ik))&
        )*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 97)=(csumtt12(1)-csumtt12(2)+csumtt12(3)-csumtt12(4)&
        -(csumtt12(9)-csumtt12(10)+csumtt12(11)-csumtt12(12))&
        -(csumtt12(7)-csumtt12(8)+csumtt12(5)-csumtt12(6))&
        +(csumtt12(15)-csumtt12(16)+csumtt12(13)-csumtt12(14)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 97)=(csumttmom12(7,ik)-csumttmom12(8,ik)&
        +csumttmom12(5,ik)-csumttmom12(6,ik)-(csumttmom12(15,ik)&
        -csumttmom12(16,ik)+csumttmom12(13,ik)-csumttmom12(14,ik))&
        )*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     13
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     tt-13 operators cp=+,pz=+,j=0
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 98)=(csumtt13(1)+csumtt13(2)+csumtt13(3)+csumtt13(4)+&
        csumtt13(5)+csumtt13(6)+csumtt13(7)+csumtt13(8))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 98)=(csumttmom13(1,ik)+csumttmom13(2,ik)&
        +csumttmom13(3,ik)+csumttmom13(4,ik)+csumttmom13(5,ik)&
        +csumttmom13(6,ik)+csumttmom13(7,ik)&
        +csumttmom13(8,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 99)=(csumtt13(1)+csumtt13(2)+csumtt13(3)+csumtt13(4)&
        -(csumtt13(5)+csumtt13(6)+csumtt13(7)+csumtt13(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 99)=(csumttmom13(1,ik)+csumttmom13(2,ik)&
        +csumttmom13(3,ik)+csumttmom13(4,ik)-(csumttmom13(5,ik)&
        +csumttmom13(6,ik)+csumttmom13(7,ik)&
        +csumttmom13(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 100)=(csumtt13(1)+giot*csumtt13(2)-csumtt13(3)&
        -giot*csumtt13(4)+(csumtt13(7)+giot*csumtt13(8)-csumtt13(5)&
        -giot*csumtt13(6)))*adiv2
    !**********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 100)=(csumttmom13(1,ik)+giot*csumttmom13(2,ik)&
        -csumttmom13(3,ik)-giot*csumttmom13(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 101)=(csumtt13(1)+giot*csumtt13(2)-csumtt13(3)&
        -giot*csumtt13(4)&
        -(csumtt13(7)+giot*csumtt13(8)-csumtt13(5)&
        -giot*csumtt13(6)))*adiv2
    !**********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 101)=(csumttmom13(7,ik)+giot*csumttmom13(8,ik)&
        -csumttmom13(5,ik)-giot*csumttmom13(6,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 102)=(csumtt13(1)-csumtt13(2)+csumtt13(3)-csumtt13(4)+&
        csumtt13(6)-csumtt13(5)+csumtt13(8)-csumtt13(7))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 102)=(csumttmom13(1,ik)-csumttmom13(2,ik)&
        +csumttmom13(3,ik)-csumttmom13(4,ik)+csumttmom13(6,ik)&
        -csumttmom13(5,ik)+csumttmom13(8,ik)&
        -csumttmom13(7,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 103)=(csumtt13(1)-csumtt13(2)+csumtt13(3)-csumtt13(4)&
        -(csumtt13(6)-csumtt13(5)+csumtt13(8)-csumtt13(7)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 103)=(csumttmom13(1,ik)-csumttmom13(2,ik)&
        +csumttmom13(3,ik)-csumttmom13(4,ik)-(csumttmom13(6,ik)&
        -csumttmom13(5,ik)+csumttmom13(8,ik)&
        -csumttmom13(7,ik)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     14 thissssssss!!!!!!!
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 104)=(csumtt14(1)+csumtt14(2)+csumtt14(3)+csumtt14(4)+&
        csumtt14(5)+csumtt14(6)+csumtt14(7)+csumtt14(8))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=0
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 104)=(csumttmom14(1,ik)+csumttmom14(2,ik)&
        +csumttmom14(3,ik)+csumttmom14(4,ik)+csumttmom14(5,ik)&
        +csumttmom14(6,ik)+csumttmom14(7,ik)&
        +csumttmom14(8,ik))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 105)=(csumtt14(1)+csumtt14(2)+csumtt14(3)+csumtt14(4)-&
        (csumtt14(5)+csumtt14(6)+csumtt14(7)+csumtt14(8)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 105)=(csumttmom14(1,ik)+csumttmom14(2,ik)&
        +csumttmom14(3,ik)+csumttmom14(4,ik)-(csumttmom14(5,ik)&
        +csumttmom14(6,ik)+csumttmom14(7,ik)&
        +csumttmom14(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 106)=(csumtt14(1)+giot*csumtt14(2)-csumtt14(3)&
        -giot*csumtt14(4)+(csumtt14(6)+giot*csumtt14(7)-csumtt14(8)&
        -giot*csumtt14(5)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 106)=(csumttmom14(1,ik)+giot*csumttmom14(2,ik)&
        -csumttmom14(3,ik)-giot*csumttmom14(4,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 107)=(csumtt14(1)+giot*csumtt14(2)-csumtt14(3)&
        -giot*csumtt14(4)-(csumtt14(6)+giot*csumtt14(7)-csumtt14(8)&
        -giot*csumtt14(5)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 107)=(csumttmom14(6,ik)+giot*csumttmom14(7,ik)&
        -csumttmom14(8,ik)-giot*csumttmom14(5,ik))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 108)=(csumtt14(1)-csumtt14(2)+csumtt14(3)-csumtt14(4)+&
        (csumtt14(7)-csumtt14(6)+csumtt14(5)-csumtt14(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
            momentum_lines(time_slice, id, ik, 108)=(csumttmom14(1,ik)-csumttmom14(2,ik)&
        +csumttmom14(3,ik)-csumttmom14(4,ik)+(csumttmom14(7,ik)&
        -csumttmom14(6,ik)+csumttmom14(5,ik)-csumttmom14(8,ik)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 109)=(csumtt14(1)-csumtt14(2)+csumtt14(3)-csumtt14(4)-&
        (csumtt14(7)-csumtt14(6)+csumtt14(5)-csumtt14(8)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2,3,4
    !**********************************************************************
        do ik=1,2
        momentum_lines(time_slice, id, ik, 109)=(csumttmom14(1,ik)-csumttmom14(2,ik)&
        +csumttmom14(3,ik)-csumttmom14(4,ik)-(csumttmom14(7,ik)&
        -csumttmom14(6,ik)+csumttmom14(5,ik)-csumttmom14(8,ik)))*adiv2
        enddo
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operator 1
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 110)=sum(csumplq8(:,1))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 110)=(csumplqmom8(1,ik,1)+csumplqmom8(2,ik,1)&
        +csumplqmom8(3,ik,1)+csumplqmom8(4,ik,1)+csumplqmom8(5,ik,1)&
        +csumplqmom8(6,ik,1)+csumplqmom8(7,ik,1)+csumplqmom8(8,ik,1))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 111)=(csumplq8(1,1)+csumplq8(2,1)+csumplq8(3,1)+csumplq8(4,1)&
        -csumplq8(5,1)-csumplq8(6,1)-csumplq8(7,1)-csumplq8(8,1))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 111)=(csumplqmom8(1,ik,1)+csumplqmom8(2,ik,1)&
        +csumplqmom8(3,ik,1)+csumplqmom8(4,ik,1)-csumplqmom8(5,ik,1)&
        -csumplqmom8(6,ik,1)-csumplqmom8(7,ik,1)-csumplqmom8(8,ik,1))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 112)=(csumplq8(1,1)+giot*csumplq8(2,1)-csumplq8(3,1)&
        -giot*csumplq8(4,1)&
        +(csumplq8(6,1)+giot*csumplq8(5,1)-csumplq8(8,1)-giot*csumplq8(7,1)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 112)=(csumplqmom8(1,ik,1)+giot*csumplqmom8(2,ik,1)&
        -csumplqmom8(3,ik,1)-giot*csumplqmom8(4,ik,1))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 113)=(csumplq8(1,1)+giot*csumplq8(2,1)-csumplq8(3,1)&
        -giot*csumplq8(4,1)&
        -(csumplq8(6,1)+giot*csumplq8(5,1)-csumplq8(8,1)-giot*csumplq8(7,1)))*adiv2
    !**********************************************************************
    !     j=1, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 113)=(csumplqmom8(6,ik,1)+giot*csumplqmom8(5,ik,1)&
        -csumplqmom8(8,ik,1)-giot*csumplqmom8(7,ik,1))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 114)=(csumplq8(1,1)-csumplq8(2,1)+csumplq8(3,1)-csumplq8(4,1)&
        +csumplq8(5,1)-csumplq8(6,1)+csumplq8(7,1)-csumplq8(8,1))*adiv2
    !**********************************************************************
    !     j=2, pr=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 114)=(csumplqmom8(1,ik,1)-csumplqmom8(2,ik,1)&
        +csumplqmom8(3,ik,1)-csumplqmom8(4,ik,1)+csumplqmom8(5,ik,1)&
        -csumplqmom8(6,ik,1)+csumplqmom8(7,ik,1)-csumplqmom8(8,ik,1))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 115)=(csumplq8(1,1)-csumplq8(2,1)+csumplq8(3,1)-csumplq8(4,1)&
        -(csumplq8(5,1)-csumplq8(6,1)+csumplq8(7,1)-csumplq8(8,1)))*adiv2
    !**********************************************************************
    !     j=2, pr=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 115)=(csumplqmom8(1,ik,1)-csumplqmom8(2,ik,1)&
        +csumplqmom8(3,ik,1)-csumplqmom8(4,ik,1)-(csumplqmom8(5,ik,1)&
        -csumplqmom8(6,ik,1)+csumplqmom8(7,ik,1)-csumplqmom8(8,ik,1)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operator 2
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 116)=(csumplq8(1,2)+csumplq8(2,2)+csumplq8(3,2)+csumplq8(4,2)&
        +csumplq8(5,2)+csumplq8(6,2)+csumplq8(7,2)+csumplq8(8,2))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 116)=(csumplqmom8(1,ik,2)+csumplqmom8(2,ik,2)&
        +csumplqmom8(3,ik,2)+csumplqmom8(4,ik,2)+csumplqmom8(5,ik,2)&
        +csumplqmom8(6,ik,2)+csumplqmom8(7,ik,2)+csumplqmom8(8,ik,2))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 117)=(csumplq8(1,2)+csumplq8(2,2)+csumplq8(3,2)+csumplq8(4,2)&
        -(csumplq8(5,2)+csumplq8(6,2)+csumplq8(7,2)+csumplq8(8,2)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 117)=(csumplqmom8(1,ik,2)+csumplqmom8(2,ik,2)&
        +csumplqmom8(3,ik,2)+csumplqmom8(4,ik,2)-(csumplqmom8(5,ik,2)&
        +csumplqmom8(6,ik,2)+csumplqmom8(7,ik,2)+csumplqmom8(8,ik,2)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 118)=(csumplq8(1,2)+giot*csumplq8(2,2)-csumplq8(3,2)&
        -giot*csumplq8(4,2)&
        +csumplq8(6,2)+giot*csumplq8(5,2)-csumplq8(8,2)-giot*csumplq8(7,2))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 118)=(csumplqmom8(1,ik,2)+giot*csumplqmom8(2,ik,2)&
        -csumplqmom8(3,ik,2)-giot*csumplqmom8(4,ik,2))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 119)=(csumplq8(1,2)+giot*csumplq8(2,2)-csumplq8(3,2)&
        -giot*csumplq8(4,2)&
        -(csumplq8(6,2)+giot*csumplq8(5,2)-csumplq8(8,2)-giot*csumplq8(7,2)))&
        *adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 119)=(csumplqmom8(6,ik,2)+giot*csumplqmom8(5,ik,2)&
        -csumplqmom8(8,ik,2)-giot*csumplqmom8(7,ik,2))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 120)=(csumplq8(1,2)-csumplq8(2,2)+csumplq8(3,2)-csumplq8(4,2)&
        +csumplq8(5,2)-csumplq8(6,2)+csumplq8(7,2)-csumplq8(8,2))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 120)=(csumplqmom8(1,ik,2)-csumplqmom8(2,ik,2)&
        +csumplqmom8(3,ik,2)-csumplqmom8(4,ik,2)+csumplqmom8(5,ik,2)&
        -csumplqmom8(6,ik,2)+csumplqmom8(7,ik,2)-csumplqmom8(8,ik,2))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 121)=(csumplq8(1,2)-csumplq8(2,2)+csumplq8(3,2)-csumplq8(4,2)&
        -(csumplq8(5,2)-csumplq8(6,2)+csumplq8(7,2)-csumplq8(8,2)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 121)=(csumplqmom8(1,ik,2)-csumplqmom8(2,ik,2)&
        +csumplqmom8(3,ik,2)-csumplqmom8(4,ik,2)-(csumplqmom8(5,ik,2)&
        -csumplqmom8(6,ik,2)+csumplqmom8(7,ik,2)-csumplqmom8(8,ik,2)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! plaquette operators 3
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 122)=(csumplq8(1,3)+csumplq8(2,3)+csumplq8(3,3)&
        +csumplq8(4,3)+(csumplq8(5,3)+csumplq8(6,3)+csumplq8(7,3)&
        +csumplq8(8,3)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 122)=(csumplqmom8(1,ik,3)+csumplqmom8(2,ik,3)&
        +csumplqmom8(3,ik,3)+csumplqmom8(4,ik,3)&
        +(csumplqmom8(5,ik,3)+csumplqmom8(6,ik,3)+csumplqmom8(7,ik,3)&
        +csumplqmom8(8,ik,3)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 123)=(csumplq8(1,3)+csumplq8(2,3)+csumplq8(3,3)&
        +csumplq8(4,3)-(csumplq8(5,3)+csumplq8(6,3)+csumplq8(7,3)&
        +csumplq8(8,3)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 123)=(csumplqmom8(1,ik,3)+csumplqmom8(2,ik,3)&
        +csumplqmom8(3,ik,3)+csumplqmom8(4,ik,3)&
        -(csumplqmom8(5,ik,3)+csumplqmom8(6,ik,3)+csumplqmom8(7,ik,3)&
        +csumplqmom8(8,ik,3)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 124)=(csumplq8(1,3)+giot*csumplq8(2,3)-csumplq8(3,3)&
        -giot*csumplq8(4,3))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 124)=(csumplqmom8(1,ik,3)+giot*csumplqmom8(2,ik,3)&
        -csumplqmom8(3,ik,3)-giot*csumplqmom8(4,ik,3))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 125)=(csumplq8(7,3)+giot*csumplq8(6,3)-csumplq8(5,3)&
        -giot*csumplq8(8,3))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 125)=(csumplqmom8(7,ik,3)+giot*csumplqmom8(6,ik,3)&
        -csumplqmom8(5,ik,3)-giot*csumplqmom8(8,ik,3))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 126)=(csumplq8(1,3)-csumplq8(2,3)+csumplq8(3,3)&
        -csumplq8(4,3)+(csumplq8(5,3)-csumplq8(6,3)+csumplq8(7,3)&
        -csumplq8(8,3)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 126)=(csumplqmom8(1,ik,3)-csumplqmom8(2,ik,3)&
        +csumplqmom8(3,ik,3)-csumplqmom8(4,ik,3)&
        +(csumplqmom8(5,ik,3)-csumplqmom8(6,ik,3)+csumplqmom8(7,ik,3)&
        -csumplqmom8(8,ik,3)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 127)=(csumplq8(1,3)-csumplq8(2,3)+csumplq8(3,3)&
        -csumplq8(4,3)-(csumplq8(5,3)-csumplq8(6,3)+csumplq8(7,3)&
        -csumplq8(8,3)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 127)=(csumplqmom8(1,ik,3)-csumplqmom8(2,ik,3)&
        +csumplqmom8(3,ik,3)-csumplqmom8(4,ik,3)&
        -(csumplqmom8(5,ik,3)-csumplqmom8(6,ik,3)+csumplqmom8(7,ik,3)&
        -csumplqmom8(8,ik,3)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 4
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 128)=(csumplq8(1,4)+csumplq8(2,4)+csumplq8(3,4)&
        +csumplq8(4,4)+(csumplq8(5,4)+csumplq8(6,4)+csumplq8(7,4)&
        +csumplq8(8,4)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 128)=(csumplqmom8(1,ik,4)+csumplqmom8(2,ik,4)&
        +csumplqmom8(3,ik,4)+csumplqmom8(4,ik,4)&
        +(csumplqmom8(5,ik,4)+csumplqmom8(6,ik,4)+csumplqmom8(7,ik,4)&
        +csumplqmom8(8,ik,4)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 129)=(csumplq8(1,4)+csumplq8(2,4)+csumplq8(3,4)&
        +csumplq8(4,4)-(csumplq8(5,4)+csumplq8(6,4)+csumplq8(7,4)&
        +csumplq8(8,4)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 129)=(csumplqmom8(1,ik,4)+csumplqmom8(2,ik,4)&
        +csumplqmom8(3,ik,4)+csumplqmom8(4,ik,4)&
        -(csumplqmom8(5,ik,4)+csumplqmom8(6,ik,4)+csumplqmom8(7,ik,4)&
        +csumplqmom8(8,ik,4)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 130)=(csumplq8(1,4)+giot*csumplq8(2,4)-csumplq8(3,4)&
        -giot*csumplq8(4,4))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 130)=(csumplqmom8(1,ik,4)+giot*csumplqmom8(2,ik,4)&
        -csumplqmom8(3,ik,4)-giot*csumplqmom8(4,ik,4))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 131)=(csumplq8(7,4)+giot*csumplq8(6,4)-csumplq8(5,4)&
        -giot*csumplq8(8,4))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 131)=(csumplqmom8(7,ik,4)+giot*csumplqmom8(6,ik,4)&
        -csumplqmom8(5,ik,4)-giot*csumplqmom8(8,ik,4))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 132)=(csumplq8(1,4)-csumplq8(2,4)+csumplq8(3,4)&
        -csumplq8(4,4)+(csumplq8(5,4)-csumplq8(6,4)+csumplq8(7,4)&
        -csumplq8(8,4)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 132)=(csumplqmom8(1,ik,4)-csumplqmom8(2,ik,4)&
        +csumplqmom8(3,ik,4)-csumplqmom8(4,ik,4)&
        +(csumplqmom8(5,ik,4)-csumplqmom8(6,ik,4)+csumplqmom8(7,ik,4)&
        -csumplqmom8(8,ik,4)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 133)=(csumplq8(1,4)-csumplq8(2,4)+csumplq8(3,4)&
        -csumplq8(4,4)-(csumplq8(5,4)-csumplq8(6,4)+csumplq8(7,4)&
        -csumplq8(8,4)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 133)=(csumplqmom8(1,ik,4)-csumplqmom8(2,ik,4)&
        +csumplqmom8(3,ik,4)-csumplqmom8(4,ik,4)&
        -(csumplqmom8(5,ik,4)-csumplqmom8(6,ik,4)+csumplqmom8(7,ik,4)&
        -csumplqmom8(8,ik,4)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! plaquette operator 5
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 134)=(csumplq8(1,5)+csumplq8(2,5)+csumplq8(3,5)&
        +csumplq8(4,5)+(csumplq8(5,5)+csumplq8(6,5)+csumplq8(7,5)&
        +csumplq8(8,5)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 134)=(csumplqmom8(1,ik,5)+csumplqmom8(2,ik,5)&
        +csumplqmom8(3,ik,5)+csumplqmom8(4,ik,5)&
        +(csumplqmom8(5,ik,5)+csumplqmom8(6,ik,5)+csumplqmom8(7,ik,5)&
        +csumplqmom8(8,ik,5)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 135)=(csumplq8(1,5)+csumplq8(2,5)+csumplq8(3,5)&
        +csumplq8(4,5)-(csumplq8(5,5)+csumplq8(6,5)+csumplq8(7,5)&
        +csumplq8(8,5)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 135)=(csumplqmom8(1,ik,5)+csumplqmom8(2,ik,5)&
        +csumplqmom8(3,ik,5)+csumplqmom8(4,ik,5)&
        -(csumplqmom8(5,ik,5)+csumplqmom8(6,ik,5)+csumplqmom8(7,ik,5)&
        +csumplqmom8(8,ik,5)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 136)=(csumplq8(1,5)+giot*csumplq8(2,5)-csumplq8(3,5)&
        -giot*csumplq8(4,5))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 136)=(csumplqmom8(1,ik,5)+giot*csumplqmom8(2,ik,5)&
        -csumplqmom8(3,ik,5)-giot*csumplqmom8(4,ik,5))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 137)=(csumplq8(7,5)+giot*csumplq8(6,5)-csumplq8(5,5)&
        -giot*csumplq8(8,5))*adiv1
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 137)=(csumplqmom8(7,ik,5)+giot*csumplqmom8(6,ik,5)&
        -csumplqmom8(5,ik,5)-giot*csumplqmom8(8,ik,5))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 138)=(csumplq8(1,5)-csumplq8(2,5)+csumplq8(3,5)&
        -csumplq8(4,5)+(csumplq8(5,5)-csumplq8(6,5)+csumplq8(7,5)&
        -csumplq8(8,5)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 138)=(csumplqmom8(1,ik,5)-csumplqmom8(2,ik,5)&
        +csumplqmom8(3,ik,5)-csumplqmom8(4,ik,5)&
        +(csumplqmom8(5,ik,5)-csumplqmom8(6,ik,5)+csumplqmom8(7,ik,5)&
        -csumplqmom8(8,ik,5)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 139)=(csumplq8(1,5)-csumplq8(2,5)+csumplq8(3,5)&
        -csumplq8(4,5)-(csumplq8(5,5)-csumplq8(6,5)+csumplq8(7,5)&
        -csumplq8(8,5)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 139)=(csumplqmom8(1,ik,5)-csumplqmom8(2,ik,5)&
        +csumplqmom8(3,ik,5)-csumplqmom8(4,ik,5)&
        -(csumplqmom8(5,ik,5)-csumplqmom8(6,ik,5)+csumplqmom8(7,ik,5)&
        -csumplqmom8(8,ik,5)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 6
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 140)=(csumplq8(1,6)+csumplq8(2,6)+csumplq8(3,6)&
        +csumplq8(4,6)+(csumplq8(5,6)+csumplq8(6,6)+csumplq8(7,6)&
        +csumplq8(8,6)))*adiv2
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 140)=(csumplqmom8(1,ik,6)+csumplqmom8(2,ik,6)&
        +csumplqmom8(3,ik,6)+csumplqmom8(4,ik,6)&
        +(csumplqmom8(5,ik,6)+csumplqmom8(6,ik,6)+csumplqmom8(7,ik,6)&
        +csumplqmom8(8,ik,6)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 141)=(csumplq8(1,6)+csumplq8(2,6)+csumplq8(3,6)&
        +csumplq8(4,6)-(csumplq8(5,6)+csumplq8(6,6)+csumplq8(7,6)&
            +csumplq8(8,6)))*adiv2
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 141)=(csumplqmom8(1,ik,6)+csumplqmom8(2,ik,6)&
        +csumplqmom8(3,ik,6)+csumplqmom8(4,ik,6)&
        -(csumplqmom8(5,ik,6)+csumplqmom8(6,ik,6)+csumplqmom8(7,ik,6)&
        +csumplqmom8(8,ik,6)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 142)=(csumplq8(1,6)+giot*csumplq8(2,6)-csumplq8(3,6)&
        -giot*csumplq8(4,6)+(csumplq8(8,6)+giot*csumplq8(7,6)-csumplq8(6,6)&
        -giot*csumplq8(5,6)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 142)=(csumplqmom8(1,ik,6)+giot*csumplqmom8(2,ik,6)&
        -csumplqmom8(3,ik,6)-giot*csumplqmom8(4,ik,6))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 143)=(csumplq8(1,6)+giot*csumplq8(2,6)-csumplq8(3,6)&
        -giot*csumplq8(4,6)-(csumplq8(8,6)+giot*csumplq8(7,6)-csumplq8(6,6)&
        -giot*csumplq8(5,6)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 143)=(csumplqmom8(8,ik,6)+giot*csumplqmom8(7,ik,6)&
        -csumplqmom8(6,ik,6)-giot*csumplqmom8(5,ik,6))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 144)=(csumplq8(1,6)-csumplq8(2,6)+csumplq8(3,6)&
        -csumplq8(4,6)+(csumplq8(5,6)-csumplq8(6,6)+csumplq8(7,6)&
        -csumplq8(8,6)))*adiv2
    !**********************************************************************
    !     j=2, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 144)=(csumplqmom8(1,ik,6)-csumplqmom8(2,ik,6)&
        +csumplqmom8(3,ik,6)-csumplqmom8(4,ik,6)&
        +(csumplqmom8(5,ik,6)-csumplqmom8(6,ik,6)+csumplqmom8(7,ik,6)&
        -csumplqmom8(8,ik,6)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 145)=(csumplq8(1,6)-csumplq8(2,6)+csumplq8(3,6)&
        -csumplq8(4,6)-(csumplq8(5,6)-csumplq8(6,6)+csumplq8(7,6)&
        -csumplq8(8,6)))*adiv2
    !**********************************************************************
    !     j=2, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 145)=(csumplqmom8(1,ik,6)-csumplqmom8(2,ik,6)&
        +csumplqmom8(3,ik,6)-csumplqmom8(4,ik,6)&
        -(csumplqmom8(5,ik,6)-csumplqmom8(6,ik,6)+csumplqmom8(7,ik,6)&
        -csumplqmom8(8,ik,6)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 7
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 146)=(csumplq16(1,7)+csumplq16(2,7)+csumplq16(3,7)&
        +csumplq16(4,7)&
        +(csumplq16(5,7)+csumplq16(6,7)+csumplq16(7,7)+csumplq16(8,7))&
        +(csumplq16(9,7)+csumplq16(10,7)+csumplq16(11,7)+csumplq16(12,7))&
        +(csumplq16(13,7)+csumplq16(14,7)+csumplq16(15,7)+csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 146)=(csumplqmom16(1,ik,7)+csumplqmom16(2,ik,7)&
        +csumplqmom16(3,ik,7)+csumplqmom16(4,ik,7)&
        +(csumplqmom16(9,ik,7)+csumplqmom16(10,ik,7)+csumplqmom16(11,ik,7)&
        +csumplqmom16(12,ik,7)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 147)=(csumplq16(1,7)+csumplq16(2,7)+csumplq16(3,7)&
        +csumplq16(4,7)&
        -(csumplq16(5,7)+csumplq16(6,7)+csumplq16(7,7)+csumplq16(8,7))&
        +(csumplq16(9,7)+csumplq16(10,7)+csumplq16(11,7)+csumplq16(12,7))&
        -(csumplq16(13,7)+csumplq16(14,7)+csumplq16(15,7)+csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 147)=(csumplqmom16(5,ik,7)+csumplqmom16(6,ik,7)&
        +csumplqmom16(7,ik,7)+csumplqmom16(8,ik,7)&
        +(csumplqmom16(13,ik,7)+csumplqmom16(14,ik,7)+csumplqmom16(15,ik,7)&
        +csumplqmom16(16,ik,7)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 148)=(csumplq16(1,7)+csumplq16(2,7)+csumplq16(3,7)&
        +csumplq16(4,7)&
        +(csumplq16(5,7)+csumplq16(6,7)+csumplq16(7,7)+csumplq16(8,7))&
        -(csumplq16(9,7)+csumplq16(10,7)+csumplq16(11,7)+csumplq16(12,7))&
        -(csumplq16(13,7)+csumplq16(14,7)+csumplq16(15,7)+csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 148)=(csumplqmom16(1,ik,7)+csumplqmom16(2,ik,7)&
        +csumplqmom16(3,ik,7)+csumplqmom16(4,ik,7)&
        -(csumplqmom16(9,ik,7)+csumplqmom16(10,ik,7)+csumplqmom16(11,ik,7)&
        +csumplqmom16(12,ik,7)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 149)=(csumplq16(1,7)+csumplq16(2,7)+csumplq16(3,7)&
        +csumplq16(4,7)&
        -(csumplq16(5,7)+csumplq16(6,7)+csumplq16(7,7)+csumplq16(8,7))&
        -(csumplq16(9,7)+csumplq16(10,7)+csumplq16(11,7)+csumplq16(12,7))&
        +(csumplq16(13,7)+csumplq16(14,7)+csumplq16(15,7)+csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 149)=(csumplqmom16(5,ik,7)+csumplqmom16(6,ik,7)&
        +csumplqmom16(7,ik,7)+csumplqmom16(8,ik,7)&
        -(csumplqmom16(13,ik,7)+csumplqmom16(14,ik,7)+csumplqmom16(15,ik,7)&
        +csumplqmom16(16,ik,7)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 150)=(csumplq16(1,7)+giot*csumplq16(2,7)-csumplq16(3,7)&
        -giot*csumplq16(4,7)+(csumplq16(5,7)+giot*csumplq16(6,7)-csumplq16(7,7)&
        -giot*csumplq16(8,7)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 150)=(csumplqmom16(1,ik,7)+giot*csumplqmom16(2,ik,7)&
        -csumplqmom16(3,ik,7)-giot*csumplqmom16(4,ik,7))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 151)=(csumplq16(1,7)+giot*csumplq16(2,7)-csumplq16(3,7)&
        -giot*csumplq16(4,7)-(csumplq16(5,7)+giot*csumplq16(6,7)-csumplq16(7,7)&
        -giot*csumplq16(8,7)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 151)=(csumplqmom16(5,ik,7)+giot*csumplqmom16(6,ik,7)&
        -csumplqmom16(7,ik,7)-giot*csumplqmom16(8,ik,7))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 152)=(csumplq16(1,7)-csumplq16(2,7)+csumplq16(3,7)&
        -csumplq16(4,7)&
        +(csumplq16(5,7)-csumplq16(6,7)+csumplq16(7,7)-csumplq16(8,7))&
        +(csumplq16(9,7)-csumplq16(10,7)+csumplq16(11,7)-csumplq16(12,7))&
        +(csumplq16(13,7)-csumplq16(14,7)+csumplq16(15,7)-csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 152)=(csumplqmom16(1,ik,7)-csumplqmom16(2,ik,7)&
        +csumplqmom16(3,ik,7)-csumplqmom16(4,ik,7)&
        +(csumplqmom16(9,ik,7)-csumplqmom16(10,ik,7)+csumplqmom16(11,ik,7)&
        -csumplqmom16(12,ik,7)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 153)=(csumplq16(1,7)-csumplq16(2,7)+csumplq16(3,7)&
        -csumplq16(4,7)&
        -(csumplq16(5,7)-csumplq16(6,7)+csumplq16(7,7)-csumplq16(8,7))&
        +(csumplq16(9,7)-csumplq16(10,7)+csumplq16(11,7)-csumplq16(12,7))&
        -(csumplq16(13,7)-csumplq16(14,7)+csumplq16(15,7)-csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 153)=(csumplqmom16(5,ik,7)-csumplqmom16(6,ik,7)&
        +csumplqmom16(7,ik,7)-csumplqmom16(8,ik,7)&
        +(csumplqmom16(13,ik,7)-csumplqmom16(14,ik,7)+csumplqmom16(15,ik,7)&
        -csumplqmom16(16,ik,7)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 154)=(csumplq16(1,7)-csumplq16(2,7)+csumplq16(3,7)&
        -csumplq16(4,7)&
        +(csumplq16(5,7)-csumplq16(6,7)+csumplq16(7,7)-csumplq16(8,7))&
        -(csumplq16(9,7)-csumplq16(10,7)+csumplq16(11,7)-csumplq16(12,7))&
        -(csumplq16(13,7)-csumplq16(14,7)+csumplq16(15,7)-csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 154)=(csumplqmom16(1,ik,7)-csumplqmom16(2,ik,7)&
        +csumplqmom16(3,ik,7)-csumplqmom16(4,ik,7)&
        -(csumplqmom16(9,ik,7)-csumplqmom16(10,ik,7)+csumplqmom16(11,ik,7)&
        -csumplqmom16(12,ik,7)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 155)=(csumplq16(1,7)-csumplq16(2,7)+csumplq16(3,7)&
        -csumplq16(4,7)&
        -(csumplq16(5,7)-csumplq16(6,7)+csumplq16(7,7)-csumplq16(8,7))&
        -(csumplq16(9,7)-csumplq16(10,7)+csumplq16(11,7)-csumplq16(12,7))&
        +(csumplq16(13,7)-csumplq16(14,7)+csumplq16(15,7)-csumplq16(16,7)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 155)=(csumplqmom16(5,ik,7)-csumplqmom16(6,ik,7)&
        +csumplqmom16(7,ik,7)-csumplqmom16(8,ik,7)&
        -(csumplqmom16(13,ik,7)-csumplqmom16(14,ik,7)+csumplqmom16(15,ik,7)&
        -csumplqmom16(16,ik,7)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 8
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 156)=(csumplq16(1,8)+csumplq16(2,8)+csumplq16(3,8)&
        +csumplq16(4,8)&
        +(csumplq16(5,8)+csumplq16(6,8)+csumplq16(7,8)+csumplq16(8,8))&
        +(csumplq16(9,8)+csumplq16(10,8)+csumplq16(11,8)+csumplq16(12,8))&
        +(csumplq16(13,8)+csumplq16(14,8)+csumplq16(15,8)+csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 156)=(csumplqmom16(1,ik,8)+csumplqmom16(2,ik,8)&
        +csumplqmom16(3,ik,8)+csumplqmom16(4,ik,8)&
        +(csumplqmom16(9,ik,8)+csumplqmom16(10,ik,8)+csumplqmom16(11,ik,8)&
        +csumplqmom16(12,ik,8)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 157)=(csumplq16(1,8)+csumplq16(2,8)+csumplq16(3,8)&
        +csumplq16(4,8)&
        -(csumplq16(5,8)+csumplq16(6,8)+csumplq16(7,8)+csumplq16(8,8))&
        +(csumplq16(9,8)+csumplq16(10,8)+csumplq16(11,8)+csumplq16(12,8))&
        -(csumplq16(13,8)+csumplq16(14,8)+csumplq16(15,8)+csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 157)=(csumplqmom16(5,ik,8)+csumplqmom16(6,ik,8)&
        +csumplqmom16(7,ik,8)+csumplqmom16(8,ik,8)&
        +(csumplqmom16(13,ik,8)+csumplqmom16(14,ik,8)+csumplqmom16(15,ik,8)&
        +csumplqmom16(16,ik,8)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 158)=(csumplq16(1,8)+csumplq16(2,8)+csumplq16(3,8)&
        +csumplq16(4,8)&
        +(csumplq16(5,8)+csumplq16(6,8)+csumplq16(7,8)+csumplq16(8,8))&
        -(csumplq16(9,8)+csumplq16(10,8)+csumplq16(11,8)+csumplq16(12,8))&
        -(csumplq16(13,8)+csumplq16(14,8)+csumplq16(15,8)+csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 158)=(csumplqmom16(1,ik,8)+csumplqmom16(2,ik,8)&
        +csumplqmom16(3,ik,8)+csumplqmom16(4,ik,8)&
        -(csumplqmom16(9,ik,8)+csumplqmom16(10,ik,8)+csumplqmom16(11,ik,8)&
        +csumplqmom16(12,ik,8)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 159)=(csumplq16(1,8)+csumplq16(2,8)+csumplq16(3,8)&
        +csumplq16(4,8)&
        -(csumplq16(5,8)+csumplq16(6,8)+csumplq16(7,8)+csumplq16(8,8))&
        -(csumplq16(9,8)+csumplq16(10,8)+csumplq16(11,8)+csumplq16(12,8))&
        +(csumplq16(13,8)+csumplq16(14,8)+csumplq16(15,8)+csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 159)=(csumplqmom16(5,ik,8)+csumplqmom16(6,ik,8)&
        +csumplqmom16(7,ik,8)+csumplqmom16(8,ik,8)&
        -(csumplqmom16(13,ik,8)+csumplqmom16(14,ik,8)+csumplqmom16(15,ik,8)&
        +csumplqmom16(16,ik,8)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 160)=(csumplq16(1,8)+giot*csumplq16(2,8)-csumplq16(3,8)&
        -giot*csumplq16(4,8)+(csumplq16(5,8)+giot*csumplq16(6,8)-csumplq16(7,8)&
        -giot*csumplq16(8,8)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 160)=(csumplqmom16(1,ik,8)+giot*csumplqmom16(2,ik,8)&
        -csumplqmom16(3,ik,8)-giot*csumplqmom16(4,ik,8))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 161)=(csumplq16(1,8)+giot*csumplq16(2,8)-csumplq16(3,8)&
        -giot*csumplq16(4,8)-(csumplq16(5,8)+giot*csumplq16(6,8)-csumplq16(7,8)&
        -giot*csumplq16(8,8)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 161)=(csumplqmom16(5,ik,8)+giot*csumplqmom16(6,ik,8)&
        -csumplqmom16(7,ik,8)-giot*csumplqmom16(8,ik,8))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 162)=(csumplq16(1,8)-csumplq16(2,8)+csumplq16(3,8)&
        -csumplq16(4,8)&
        +(csumplq16(5,8)-csumplq16(6,8)+csumplq16(7,8)-csumplq16(8,8))&
        +(csumplq16(9,8)-csumplq16(10,8)+csumplq16(11,8)-csumplq16(12,8))&
        +(csumplq16(13,8)-csumplq16(14,8)+csumplq16(15,8)-csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 162)=(csumplqmom16(1,ik,8)-csumplqmom16(2,ik,8)&
        +csumplqmom16(3,ik,8)-csumplqmom16(4,ik,8)&
        +(csumplqmom16(9,ik,8)-csumplqmom16(10,ik,8)+csumplqmom16(11,ik,8)&
        -csumplqmom16(12,ik,8)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 163)=(csumplq16(1,8)-csumplq16(2,8)+csumplq16(3,8)&
        -csumplq16(4,8)&
        -(csumplq16(5,8)-csumplq16(6,8)+csumplq16(7,8)-csumplq16(8,8))&
        +(csumplq16(9,8)-csumplq16(10,8)+csumplq16(11,8)-csumplq16(12,8))&
        -(csumplq16(13,8)-csumplq16(14,8)+csumplq16(15,8)-csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 163)=(csumplqmom16(5,ik,8)-csumplqmom16(6,ik,8)&
        +csumplqmom16(7,ik,8)-csumplqmom16(8,ik,8)&
        +(csumplqmom16(13,ik,8)-csumplqmom16(14,ik,8)+csumplqmom16(15,ik,8)&
        -csumplqmom16(16,ik,8)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 164)=(csumplq16(1,8)-csumplq16(2,8)+csumplq16(3,8)&
        -csumplq16(4,8)&
        +(csumplq16(5,8)-csumplq16(6,8)+csumplq16(7,8)-csumplq16(8,8))&
        -(csumplq16(9,8)-csumplq16(10,8)+csumplq16(11,8)-csumplq16(12,8))&
        -(csumplq16(13,8)-csumplq16(14,8)+csumplq16(15,8)-csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 164)=(csumplqmom16(1,ik,8)-csumplqmom16(2,ik,8)&
        +csumplqmom16(3,ik,8)-csumplqmom16(4,ik,8)&
        -(csumplqmom16(9,ik,8)-csumplqmom16(10,ik,8)+csumplqmom16(11,ik,8)&
        -csumplqmom16(12,ik,8)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 165)=(csumplq16(1,8)-csumplq16(2,8)+csumplq16(3,8)&
        -csumplq16(4,8)&
        -(csumplq16(5,8)-csumplq16(6,8)+csumplq16(7,8)-csumplq16(8,8))&
        -(csumplq16(9,8)-csumplq16(10,8)+csumplq16(11,8)-csumplq16(12,8))&
        +(csumplq16(13,8)-csumplq16(14,8)+csumplq16(15,8)-csumplq16(16,8)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 165)=(csumplqmom16(5,ik,8)-csumplqmom16(6,ik,8)&
        +csumplqmom16(7,ik,8)-csumplqmom16(8,ik,8)&
        -(csumplqmom16(13,ik,8)-csumplqmom16(14,ik,8)+csumplqmom16(15,ik,8)&
        -csumplqmom16(16,ik,8)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 9
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 166)=(csumplq16(1,9)+csumplq16(2,9)+csumplq16(3,9)&
        +csumplq16(4,9)&
        +(csumplq16(5,9)+csumplq16(6,9)+csumplq16(7,9)+csumplq16(8,9))&
        +(csumplq16(9,9)+csumplq16(10,9)+csumplq16(11,9)+csumplq16(12,9))&
        +(csumplq16(13,9)+csumplq16(14,9)+csumplq16(15,9)+csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 166)=(csumplqmom16(1,ik,9)+csumplqmom16(2,ik,9)&
        +csumplqmom16(3,ik,9)+csumplqmom16(4,ik,9)&
        +(csumplqmom16(9,ik,9)+csumplqmom16(10,ik,9)+csumplqmom16(11,ik,9)&
        +csumplqmom16(12,ik,9)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 167)=(csumplq16(1,9)+csumplq16(2,9)+csumplq16(3,9)&
        +csumplq16(4,9)&
        -(csumplq16(5,9)+csumplq16(6,9)+csumplq16(7,9)+csumplq16(8,9))&
        +(csumplq16(9,9)+csumplq16(10,9)+csumplq16(11,9)+csumplq16(12,9))&
        -(csumplq16(13,9)+csumplq16(14,9)+csumplq16(15,9)+csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 167)=(csumplqmom16(5,ik,9)+csumplqmom16(6,ik,9)&
        +csumplqmom16(7,ik,9)+csumplqmom16(8,ik,9)&
        +(csumplqmom16(13,ik,9)+csumplqmom16(14,ik,9)+csumplqmom16(15,ik,9)&
        +csumplqmom16(16,ik,9)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 168)=(csumplq16(1,9)+csumplq16(2,9)+csumplq16(3,9)&
        +csumplq16(4,9)&
        +(csumplq16(5,9)+csumplq16(6,9)+csumplq16(7,9)+csumplq16(8,9))&
        -(csumplq16(9,9)+csumplq16(10,9)+csumplq16(11,9)+csumplq16(12,9))&
        -(csumplq16(13,9)+csumplq16(14,9)+csumplq16(15,9)+csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 168)=(csumplqmom16(1,ik,9)+csumplqmom16(2,ik,9)&
        +csumplqmom16(3,ik,9)+csumplqmom16(4,ik,9)&
        -(csumplqmom16(9,ik,9)+csumplqmom16(10,ik,9)+csumplqmom16(11,ik,9)&
        +csumplqmom16(12,ik,9)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 169)=(csumplq16(1,9)+csumplq16(2,9)+csumplq16(3,9)&
        +csumplq16(4,9)&
        -(csumplq16(5,9)+csumplq16(6,9)+csumplq16(7,9)+csumplq16(8,9))&
        -(csumplq16(9,9)+csumplq16(10,9)+csumplq16(11,9)+csumplq16(12,9))&
        +(csumplq16(13,9)+csumplq16(14,9)+csumplq16(15,9)+csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 169)=(csumplqmom16(5,ik,9)+csumplqmom16(6,ik,9)&
        +csumplqmom16(7,ik,9)+csumplqmom16(8,ik,9)&
        -(csumplqmom16(13,ik,9)+csumplqmom16(14,ik,9)+csumplqmom16(15,ik,9)&
        +csumplqmom16(16,ik,9)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 170)=(csumplq16(1,9)+giot*csumplq16(2,9)-csumplq16(3,9)&
        -giot*csumplq16(4,9)+(csumplq16(5,9)+giot*csumplq16(6,9)-csumplq16(7,9)&
        -giot*csumplq16(8,9)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 170)=(csumplqmom16(1,ik,9)+giot*csumplqmom16(2,ik,9)&
        -csumplqmom16(3,ik,9)-giot*csumplqmom16(4,ik,9))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 171)=(csumplq16(1,9)+giot*csumplq16(2,9)-csumplq16(3,9)&
        -giot*csumplq16(4,9)-(csumplq16(5,9)+giot*csumplq16(6,9)-csumplq16(7,9)&
        -giot*csumplq16(8,9)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 171)=(csumplqmom16(5,ik,9)+giot*csumplqmom16(6,ik,9)&
        -csumplqmom16(7,ik,9)-giot*csumplqmom16(8,ik,9))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 172)=(csumplq16(1,9)-csumplq16(2,9)+csumplq16(3,9)&
        -csumplq16(4,9)&
        +(csumplq16(5,9)-csumplq16(6,9)+csumplq16(7,9)-csumplq16(8,9))&
        +(csumplq16(9,9)-csumplq16(10,9)+csumplq16(11,9)-csumplq16(12,9))&
        +(csumplq16(13,9)-csumplq16(14,9)+csumplq16(15,9)-csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 172)=(csumplqmom16(1,ik,9)-csumplqmom16(2,ik,9)&
        +csumplqmom16(3,ik,9)-csumplqmom16(4,ik,9)&
        +(csumplqmom16(9,ik,9)-csumplqmom16(10,ik,9)+csumplqmom16(11,ik,9)&
        -csumplqmom16(12,ik,9)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 173)=(csumplq16(1,9)-csumplq16(2,9)+csumplq16(3,9)&
        -csumplq16(4,9)&
        -(csumplq16(5,9)-csumplq16(6,9)+csumplq16(7,9)-csumplq16(8,9))&
        +(csumplq16(9,9)-csumplq16(10,9)+csumplq16(11,9)-csumplq16(12,9))&
        -(csumplq16(13,9)-csumplq16(14,9)+csumplq16(15,9)-csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 173)=(csumplqmom16(5,ik,9)-csumplqmom16(6,ik,9)&
        +csumplqmom16(7,ik,9)-csumplqmom16(8,ik,9)&
        +(csumplqmom16(13,ik,9)-csumplqmom16(14,ik,9)+csumplqmom16(15,ik,9)&
        -csumplqmom16(16,ik,9)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 174)=(csumplq16(1,9)-csumplq16(2,9)+csumplq16(3,9)&
        -csumplq16(4,9)&
        +(csumplq16(5,9)-csumplq16(6,9)+csumplq16(7,9)-csumplq16(8,9))&
        -(csumplq16(9,9)-csumplq16(10,9)+csumplq16(11,9)-csumplq16(12,9))&
        -(csumplq16(13,9)-csumplq16(14,9)+csumplq16(15,9)-csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 174)=(csumplqmom16(1,ik,9)-csumplqmom16(2,ik,9)&
        +csumplqmom16(3,ik,9)-csumplqmom16(4,ik,9)&
        -(csumplqmom16(9,ik,9)-csumplqmom16(10,ik,9)+csumplqmom16(11,ik,9)&
        -csumplqmom16(12,ik,9)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 175)=(csumplq16(1,9)-csumplq16(2,9)+csumplq16(3,9)&
        -csumplq16(4,9)&
        -(csumplq16(5,9)-csumplq16(6,9)+csumplq16(7,9)-csumplq16(8,9))&
        -(csumplq16(9,9)-csumplq16(10,9)+csumplq16(11,9)-csumplq16(12,9))&
        +(csumplq16(13,9)-csumplq16(14,9)+csumplq16(15,9)-csumplq16(16,9)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 175)=(csumplqmom16(5,ik,9)-csumplqmom16(6,ik,9)&
        +csumplqmom16(7,ik,9)-csumplqmom16(8,ik,9)&
        -(csumplqmom16(13,ik,9)-csumplqmom16(14,ik,9)+csumplqmom16(15,ik,9)&
        -csumplqmom16(16,ik,9)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 10
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 176)=(csumplq16(1,10)+csumplq16(2,10)+csumplq16(3,10)&
        +csumplq16(4,10)&
        +(csumplq16(5,10)+csumplq16(6,10)+csumplq16(7,10)+csumplq16(8,10))&
        +(csumplq16(9,10)+csumplq16(10,10)+csumplq16(11,10)+csumplq16(12,10))&
        +(csumplq16(13,10)+csumplq16(14,10)+csumplq16(15,10)+csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 176)=(csumplqmom16(1,ik,10)+csumplqmom16(2,ik,10)&
        +csumplqmom16(3,ik,10)+csumplqmom16(4,ik,10)&
        +(csumplqmom16(9,ik,10)+csumplqmom16(10,ik,10)+csumplqmom16(11,ik,10)&
        +csumplqmom16(12,ik,10)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 177)=(csumplq16(1,10)+csumplq16(2,10)+csumplq16(3,10)&
        +csumplq16(4,10)&
        -(csumplq16(5,10)+csumplq16(6,10)+csumplq16(7,10)+csumplq16(8,10))&
        +(csumplq16(9,10)+csumplq16(10,10)+csumplq16(11,10)+csumplq16(12,10))&
        -(csumplq16(13,10)+csumplq16(14,10)+csumplq16(15,10)+csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 177)=(csumplqmom16(5,ik,10)+csumplqmom16(6,ik,10)&
        +csumplqmom16(7,ik,10)+csumplqmom16(8,ik,10)&
        +(csumplqmom16(13,ik,10)+csumplqmom16(14,ik,10)+csumplqmom16(15,ik,10)&
        +csumplqmom16(16,ik,10)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 178)=(csumplq16(1,10)+csumplq16(2,10)+csumplq16(3,10)&
        +csumplq16(4,10)&
        +(csumplq16(5,10)+csumplq16(6,10)+csumplq16(7,10)+csumplq16(8,10))&
        -(csumplq16(9,10)+csumplq16(10,10)+csumplq16(11,10)+csumplq16(12,10))&
        -(csumplq16(13,10)+csumplq16(14,10)+csumplq16(15,10)+csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 178)=(csumplqmom16(1,ik,10)+csumplqmom16(2,ik,10)&
        +csumplqmom16(3,ik,10)+csumplqmom16(4,ik,10)&
        -(csumplqmom16(9,ik,10)+csumplqmom16(10,ik,10)+csumplqmom16(11,ik,10)&
        +csumplqmom16(12,ik,10)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 179)=(csumplq16(1,10)+csumplq16(2,10)+csumplq16(3,10)&
        +csumplq16(4,10)&
        -(csumplq16(5,10)+csumplq16(6,10)+csumplq16(7,10)+csumplq16(8,10))&
        -(csumplq16(9,10)+csumplq16(10,10)+csumplq16(11,10)+csumplq16(12,10))&
        +(csumplq16(13,10)+csumplq16(14,10)+csumplq16(15,10)+csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 179)=(csumplqmom16(5,ik,10)+csumplqmom16(6,ik,10)&
        +csumplqmom16(7,ik,10)+csumplqmom16(8,ik,10)&
        -(csumplqmom16(13,ik,10)+csumplqmom16(14,ik,10)+csumplqmom16(15,ik,10)&
        +csumplqmom16(16,ik,10)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 180)=(csumplq16(1,10)+giot*csumplq16(2,10)-csumplq16(3,10)&
        -giot*csumplq16(4,10)+(csumplq16(5,10)+giot*csumplq16(6,10)-csumplq16(7,10)&
        -giot*csumplq16(8,10)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 180)=(csumplqmom16(1,ik,10)&
        +giot*csumplqmom16(2,ik,10)&
        -csumplqmom16(3,ik,10)-giot*csumplqmom16(4,ik,10))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 181)=(csumplq16(1,10)+giot*csumplq16(2,10)-csumplq16(3,10)&
        -giot*csumplq16(4,10)-(csumplq16(5,10)+giot*csumplq16(6,10)-csumplq16(7,10)&
        -giot*csumplq16(8,10)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 181)=(csumplqmom16(5,ik,10)&
        +giot*csumplqmom16(6,ik,10)&
        -csumplqmom16(7,ik,10)-giot*csumplqmom16(8,ik,10))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 182)=(csumplq16(1,10)-csumplq16(2,10)+csumplq16(3,10)&
        -csumplq16(4,10)&
        +(csumplq16(5,10)-csumplq16(6,10)+csumplq16(7,10)-csumplq16(8,10))&
        +(csumplq16(9,10)-csumplq16(10,10)+csumplq16(11,10)-csumplq16(12,10))&
        +(csumplq16(13,10)-csumplq16(14,10)+csumplq16(15,10)-csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 182)=(csumplqmom16(1,ik,10)-csumplqmom16(2,ik,10)&
        +csumplqmom16(3,ik,10)-csumplqmom16(4,ik,10)&
        +(csumplqmom16(9,ik,10)-csumplqmom16(10,ik,10)+csumplqmom16(11,ik,10)&
        -csumplqmom16(12,ik,10)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 183)=(csumplq16(1,10)-csumplq16(2,10)+csumplq16(3,10)&
        -csumplq16(4,10)&
        -(csumplq16(5,10)-csumplq16(6,10)+csumplq16(7,10)-csumplq16(8,10))&
        +(csumplq16(9,10)-csumplq16(10,10)+csumplq16(11,10)-csumplq16(12,10))&
        -(csumplq16(13,10)-csumplq16(14,10)+csumplq16(15,10)-csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 183)=(csumplqmom16(5,ik,10)-csumplqmom16(6,ik,10)&
        +csumplqmom16(7,ik,10)-csumplqmom16(8,ik,10)&
        +(csumplqmom16(13,ik,10)-csumplqmom16(14,ik,10)+csumplqmom16(15,ik,10)&
        -csumplqmom16(16,ik,10)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 184)=(csumplq16(1,10)-csumplq16(2,10)+csumplq16(3,10)&
        -csumplq16(4,10)&
        +(csumplq16(5,10)-csumplq16(6,10)+csumplq16(7,10)-csumplq16(8,10))&
        -(csumplq16(9,10)-csumplq16(10,10)+csumplq16(11,10)-csumplq16(12,10))&
        -(csumplq16(13,10)-csumplq16(14,10)+csumplq16(15,10)-csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 184)=(csumplqmom16(1,ik,10)-csumplqmom16(2,ik,10)&
        +csumplqmom16(3,ik,10)-csumplqmom16(4,ik,10)&
        -(csumplqmom16(9,ik,10)-csumplqmom16(10,ik,10)+csumplqmom16(11,ik,10)&
        -csumplqmom16(12,ik,10)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 185)=(csumplq16(1,10)-csumplq16(2,10)+csumplq16(3,10)&
        -csumplq16(4,10)&
        -(csumplq16(5,10)-csumplq16(6,10)+csumplq16(7,10)-csumplq16(8,10))&
        -(csumplq16(9,10)-csumplq16(10,10)+csumplq16(11,10)-csumplq16(12,10))&
        +(csumplq16(13,10)-csumplq16(14,10)+csumplq16(15,10)-csumplq16(16,10)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 185)=(csumplqmom16(5,ik,10)-csumplqmom16(6,ik,10)&
        +csumplqmom16(7,ik,10)-csumplqmom16(8,ik,10)&
        -(csumplqmom16(13,ik,10)-csumplqmom16(14,ik,10)+csumplqmom16(15,ik,10)&
        -csumplqmom16(16,ik,10)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 11
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 186)=(csumplq16(1,11)+csumplq16(2,11)+csumplq16(3,11)&
        +csumplq16(4,11)&
        +(csumplq16(5,11)+csumplq16(6,11)+csumplq16(7,11)+csumplq16(8,11))&
        +(csumplq16(9,11)+csumplq16(10,11)+csumplq16(11,11)+csumplq16(12,11))&
        +(csumplq16(13,11)+csumplq16(14,11)+csumplq16(15,11)+csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 186)=(csumplqmom16(1,ik,11)+csumplqmom16(2,ik,11)&
        +csumplqmom16(3,ik,11)+csumplqmom16(4,ik,11)&
        +(csumplqmom16(9,ik,11)+csumplqmom16(10,ik,11)+csumplqmom16(11,ik,11)&
        +csumplqmom16(12,ik,11)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 187)=(csumplq16(1,11)+csumplq16(2,11)+csumplq16(3,11)&
        +csumplq16(4,11)&
        -(csumplq16(5,11)+csumplq16(6,11)+csumplq16(7,11)+csumplq16(8,11))&
        +(csumplq16(9,11)+csumplq16(10,11)+csumplq16(11,11)+csumplq16(12,11))&
        -(csumplq16(13,11)+csumplq16(14,11)+csumplq16(15,11)+csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 187)=(csumplqmom16(5,ik,11)+csumplqmom16(6,ik,11)&
        +csumplqmom16(7,ik,11)+csumplqmom16(8,ik,11)&
        +(csumplqmom16(13,ik,11)+csumplqmom16(14,ik,11)+csumplqmom16(15,ik,11)&
        +csumplqmom16(16,ik,11)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 188)=(csumplq16(1,11)+csumplq16(2,11)+csumplq16(3,11)&
        +csumplq16(4,11)&
        +(csumplq16(5,11)+csumplq16(6,11)+csumplq16(7,11)+csumplq16(8,11))&
        -(csumplq16(9,11)+csumplq16(10,11)+csumplq16(11,11)+csumplq16(12,11))&
        -(csumplq16(13,11)+csumplq16(14,11)+csumplq16(15,11)+csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 188)=(csumplqmom16(1,ik,11)+csumplqmom16(2,ik,11)&
        +csumplqmom16(3,ik,11)+csumplqmom16(4,ik,11)&
        -(csumplqmom16(9,ik,11)+csumplqmom16(10,ik,11)+csumplqmom16(11,ik,11)&
        +csumplqmom16(12,ik,11)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 189)=(csumplq16(1,11)+csumplq16(2,11)+csumplq16(3,11)&
        +csumplq16(4,11)&
        -(csumplq16(5,11)+csumplq16(6,11)+csumplq16(7,11)+csumplq16(8,11))&
        -(csumplq16(9,11)+csumplq16(10,11)+csumplq16(11,11)+csumplq16(12,11))&
        +(csumplq16(13,11)+csumplq16(14,11)+csumplq16(15,11)+csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 189)=(csumplqmom16(5,ik,11)+csumplqmom16(6,ik,11)&
        +csumplqmom16(7,ik,11)+csumplqmom16(8,ik,11)&
        -(csumplqmom16(13,ik,11)+csumplqmom16(14,ik,11)+csumplqmom16(15,ik,11)&
        +csumplqmom16(16,ik,11)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 190)=(csumplq16(1,11)+giot*csumplq16(2,11)-csumplq16(3,11)&
        -giot*csumplq16(4,11)+(csumplq16(5,11)+giot*csumplq16(6,11)-csumplq16(7,11)&
        -giot*csumplq16(8,11)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 190)=(csumplqmom16(1,ik,11)&
        +giot*csumplqmom16(2,ik,11)&
        -csumplqmom16(3,ik,11)-giot*csumplqmom16(4,ik,11))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 191)=(csumplq16(1,11)+giot*csumplq16(2,11)-csumplq16(3,11)&
        -giot*csumplq16(4,11)-(csumplq16(5,11)+giot*csumplq16(6,11)-csumplq16(7,11)&
        -giot*csumplq16(8,11)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 191)=(csumplqmom16(5,ik,11)&
        +giot*csumplqmom16(6,ik,11)&
        -csumplqmom16(7,ik,11)-giot*csumplqmom16(8,ik,11))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 192)=(csumplq16(1,11)-csumplq16(2,11)+csumplq16(3,11)&
        -csumplq16(4,11)&
        +(csumplq16(5,11)-csumplq16(6,11)+csumplq16(7,11)-csumplq16(8,11))&
        +(csumplq16(9,11)-csumplq16(10,11)+csumplq16(11,11)-csumplq16(12,11))&
        +(csumplq16(13,11)-csumplq16(14,11)+csumplq16(15,11)-csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 192)=(csumplqmom16(1,ik,11)-csumplqmom16(2,ik,11)&
        +csumplqmom16(3,ik,11)-csumplqmom16(4,ik,11)&
        +(csumplqmom16(9,ik,11)-csumplqmom16(10,ik,11)+csumplqmom16(11,ik,11)&
        -csumplqmom16(12,ik,11)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 193)=(csumplq16(1,11)-csumplq16(2,11)+csumplq16(3,11)&
        -csumplq16(4,11)&
        -(csumplq16(5,11)-csumplq16(6,11)+csumplq16(7,11)-csumplq16(8,11))&
        +(csumplq16(9,11)-csumplq16(10,11)+csumplq16(11,11)-csumplq16(12,11))&
        -(csumplq16(13,11)-csumplq16(14,11)+csumplq16(15,11)-csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 193)=(csumplqmom16(5,ik,11)-csumplqmom16(6,ik,11)&
        +csumplqmom16(7,ik,11)-csumplqmom16(8,ik,11)&
        +(csumplqmom16(13,ik,11)-csumplqmom16(14,ik,11)+csumplqmom16(15,ik,11)&
        -csumplqmom16(16,ik,11)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 194)=(csumplq16(1,11)-csumplq16(2,11)+csumplq16(3,11)&
        -csumplq16(4,11)&
        +(csumplq16(5,11)-csumplq16(6,11)+csumplq16(7,11)-csumplq16(8,11))&
        -(csumplq16(9,11)-csumplq16(10,11)+csumplq16(11,11)-csumplq16(12,11))&
        -(csumplq16(13,11)-csumplq16(14,11)+csumplq16(15,11)-csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 194)=(csumplqmom16(1,ik,11)-csumplqmom16(2,ik,11)&
        +csumplqmom16(3,ik,11)-csumplqmom16(4,ik,11)&
        -(csumplqmom16(9,ik,11)-csumplqmom16(10,ik,11)+csumplqmom16(11,ik,11)&
        -csumplqmom16(12,ik,11)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 195)=(csumplq16(1,11)-csumplq16(2,11)+csumplq16(3,11)&
        -csumplq16(4,11)&
        -(csumplq16(5,11)-csumplq16(6,11)+csumplq16(7,11)-csumplq16(8,11))&
        -(csumplq16(9,11)-csumplq16(10,11)+csumplq16(11,11)-csumplq16(12,11))&
        +(csumplq16(13,11)-csumplq16(14,11)+csumplq16(15,11)-csumplq16(16,11)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 195)=(csumplqmom16(5,ik,11)-csumplqmom16(6,ik,11)&
        +csumplqmom16(7,ik,11)-csumplqmom16(8,ik,11)&
        -(csumplqmom16(13,ik,11)-csumplqmom16(14,ik,11)+csumplqmom16(15,ik,11)&
        -csumplqmom16(16,ik,11)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 12
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 196)=(csumplq16(1,12)+csumplq16(2,12)+csumplq16(3,12)&
        +csumplq16(4,12)&
        +(csumplq16(5,12)+csumplq16(6,12)+csumplq16(7,12)+csumplq16(8,12))&
        +(csumplq16(9,12)+csumplq16(10,12)+csumplq16(11,12)+csumplq16(12,12))&
        +(csumplq16(13,12)+csumplq16(14,12)+csumplq16(15,12)+csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 196)=(csumplqmom16(1,ik,12)+csumplqmom16(2,ik,12)&
        +csumplqmom16(3,ik,12)+csumplqmom16(4,ik,12)&
        +(csumplqmom16(9,ik,12)+csumplqmom16(10,ik,12)+csumplqmom16(11,ik,12)&
        +csumplqmom16(12,ik,12)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 197)=(csumplq16(1,12)+csumplq16(2,12)+csumplq16(3,12)&
        +csumplq16(4,12)&
        -(csumplq16(5,12)+csumplq16(6,12)+csumplq16(7,12)+csumplq16(8,12))&
        +(csumplq16(9,12)+csumplq16(10,12)+csumplq16(11,12)+csumplq16(12,12))&
        -(csumplq16(13,12)+csumplq16(14,12)+csumplq16(15,12)+csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 197)=(csumplqmom16(5,ik,12)+csumplqmom16(6,ik,12)&
        +csumplqmom16(7,ik,12)+csumplqmom16(8,ik,12)&
        +(csumplqmom16(13,ik,12)+csumplqmom16(14,ik,12)+csumplqmom16(15,ik,12)&
        +csumplqmom16(16,ik,12)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 198)=(csumplq16(1,12)+csumplq16(2,12)+csumplq16(3,12)&
        +csumplq16(4,12)&
        +(csumplq16(5,12)+csumplq16(6,12)+csumplq16(7,12)+csumplq16(8,12))&
        -(csumplq16(9,12)+csumplq16(10,12)+csumplq16(11,12)+csumplq16(12,12))&
        -(csumplq16(13,12)+csumplq16(14,12)+csumplq16(15,12)+csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 198)=(csumplqmom16(1,ik,12)+csumplqmom16(2,ik,12)&
        +csumplqmom16(3,ik,12)+csumplqmom16(4,ik,12)&
        -(csumplqmom16(9,ik,12)+csumplqmom16(10,ik,12)+csumplqmom16(11,ik,12)&
        +csumplqmom16(12,ik,12)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 199)=(csumplq16(1,12)+csumplq16(2,12)+csumplq16(3,12)&
        +csumplq16(4,12)&
        -(csumplq16(5,12)+csumplq16(6,12)+csumplq16(7,12)+csumplq16(8,12))&
        -(csumplq16(9,12)+csumplq16(10,12)+csumplq16(11,12)+csumplq16(12,12))&
        +(csumplq16(13,12)+csumplq16(14,12)+csumplq16(15,12)+csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 199)=(csumplqmom16(5,ik,12)+csumplqmom16(6,ik,12)&
        +csumplqmom16(7,ik,12)+csumplqmom16(8,ik,12)&
        -(csumplqmom16(13,ik,12)+csumplqmom16(14,ik,12)+csumplqmom16(15,ik,12)&
        +csumplqmom16(16,ik,12)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 200)=(csumplq16(1,12)+giot*csumplq16(2,12)-csumplq16(3,12)&
        -giot*csumplq16(4,12)+(csumplq16(5,12)+giot*csumplq16(6,12)-csumplq16(7,12)&
        -giot*csumplq16(8,12)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 200)=(csumplqmom16(1,ik,12)&
        +giot*csumplqmom16(2,ik,12)&
        -csumplqmom16(3,ik,12)-giot*csumplqmom16(4,ik,12))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 201)=(csumplq16(1,12)+giot*csumplq16(2,12)-csumplq16(3,12)&
        -giot*csumplq16(4,12)-(csumplq16(5,12)+giot*csumplq16(6,12)-csumplq16(7,12)&
        -giot*csumplq16(8,12)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 201)=(csumplqmom16(5,ik,12)&
        +giot*csumplqmom16(6,ik,12)&
        -csumplqmom16(7,ik,12)-giot*csumplqmom16(8,ik,12))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 202)=(csumplq16(1,12)-csumplq16(2,12)+csumplq16(3,12)&
        -csumplq16(4,12)&
        +(csumplq16(5,12)-csumplq16(6,12)+csumplq16(7,12)-csumplq16(8,12))&
        +(csumplq16(9,12)-csumplq16(10,12)+csumplq16(11,12)-csumplq16(12,12))&
        +(csumplq16(13,12)-csumplq16(14,12)+csumplq16(15,12)-csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 202)=(csumplqmom16(1,ik,12)-csumplqmom16(2,ik,12)&
        +csumplqmom16(3,ik,12)-csumplqmom16(4,ik,12)&
        +(csumplqmom16(9,ik,12)-csumplqmom16(10,ik,12)+csumplqmom16(11,ik,12)&
        -csumplqmom16(12,ik,12)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 203)=(csumplq16(1,12)-csumplq16(2,12)+csumplq16(3,12)&
        -csumplq16(4,12)&
        -(csumplq16(5,12)-csumplq16(6,12)+csumplq16(7,12)-csumplq16(8,12))&
        +(csumplq16(9,12)-csumplq16(10,12)+csumplq16(11,12)-csumplq16(12,12))&
        -(csumplq16(13,12)-csumplq16(14,12)+csumplq16(15,12)-csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 203)=(csumplqmom16(5,ik,12)-csumplqmom16(6,ik,12)&
        +csumplqmom16(7,ik,12)-csumplqmom16(8,ik,12)&
        +(csumplqmom16(13,ik,12)-csumplqmom16(14,ik,12)+csumplqmom16(15,ik,12)&
        -csumplqmom16(16,ik,12)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 204)=(csumplq16(1,12)-csumplq16(2,12)+csumplq16(3,12)&
        -csumplq16(4,12)&
        +(csumplq16(5,12)-csumplq16(6,12)+csumplq16(7,12)-csumplq16(8,12))&
        -(csumplq16(9,12)-csumplq16(10,12)+csumplq16(11,12)-csumplq16(12,12))&
        -(csumplq16(13,12)-csumplq16(14,12)+csumplq16(15,12)-csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 204)=(csumplqmom16(1,ik,12)-csumplqmom16(2,ik,12)&
        +csumplqmom16(3,ik,12)-csumplqmom16(4,ik,12)&
        -(csumplqmom16(9,ik,12)-csumplqmom16(10,ik,12)+csumplqmom16(11,ik,12)&
        -csumplqmom16(12,ik,12)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 205)=(csumplq16(1,12)-csumplq16(2,12)+csumplq16(3,12)&
        -csumplq16(4,12)&
        -(csumplq16(5,12)-csumplq16(6,12)+csumplq16(7,12)-csumplq16(8,12))&
        -(csumplq16(9,12)-csumplq16(10,12)+csumplq16(11,12)-csumplq16(12,12))&
        +(csumplq16(13,12)-csumplq16(14,12)+csumplq16(15,12)-csumplq16(16,12)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 205)=(csumplqmom16(5,ik,12)-csumplqmom16(6,ik,12)&
        +csumplqmom16(7,ik,12)-csumplqmom16(8,ik,12)&
        -(csumplqmom16(13,ik,12)-csumplqmom16(14,ik,12)+csumplqmom16(15,ik,12)&
        -csumplqmom16(16,ik,12)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 13
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 206)=(csumplq16(1,13)+csumplq16(2,13)+csumplq16(3,13)&
        +csumplq16(4,13)&
        +(csumplq16(5,13)+csumplq16(6,13)+csumplq16(7,13)+csumplq16(8,13))&
        +(csumplq16(9,13)+csumplq16(10,13)+csumplq16(11,13)+csumplq16(12,13))&
        +(csumplq16(13,13)+csumplq16(14,13)+csumplq16(15,13)+csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 206)=(csumplqmom16(1,ik,13)+csumplqmom16(2,ik,13)&
        +csumplqmom16(3,ik,13)+csumplqmom16(4,ik,13)&
        +(csumplqmom16(9,ik,13)+csumplqmom16(10,ik,13)+csumplqmom16(11,ik,13)&
        +csumplqmom16(12,ik,13)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 207)=(csumplq16(1,13)+csumplq16(2,13)+csumplq16(3,13)&
        +csumplq16(4,13)&
        -(csumplq16(5,13)+csumplq16(6,13)+csumplq16(7,13)+csumplq16(8,13))&
        +(csumplq16(9,13)+csumplq16(10,13)+csumplq16(11,13)+csumplq16(12,13))&
        -(csumplq16(13,13)+csumplq16(14,13)+csumplq16(15,13)+csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 207)=(csumplqmom16(5,ik,13)+csumplqmom16(6,ik,13)&
        +csumplqmom16(7,ik,13)+csumplqmom16(8,ik,13)&
        +(csumplqmom16(13,ik,13)+csumplqmom16(14,ik,13)+csumplqmom16(15,ik,13)&
        +csumplqmom16(16,ik,13)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 208)=(csumplq16(1,13)+csumplq16(2,13)+csumplq16(3,13)&
        +csumplq16(4,13)&
        +(csumplq16(5,13)+csumplq16(6,13)+csumplq16(7,13)+csumplq16(8,13))&
        -(csumplq16(9,13)+csumplq16(10,13)+csumplq16(11,13)+csumplq16(12,13))&
        -(csumplq16(13,13)+csumplq16(14,13)+csumplq16(15,13)+csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 208)=(csumplqmom16(1,ik,13)+csumplqmom16(2,ik,13)&
        +csumplqmom16(3,ik,13)+csumplqmom16(4,ik,13)&
        -(csumplqmom16(9,ik,13)+csumplqmom16(10,ik,13)+csumplqmom16(11,ik,13)&
        +csumplqmom16(12,ik,13)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 209)=(csumplq16(1,13)+csumplq16(2,13)+csumplq16(3,13)&
        +csumplq16(4,13)&
        -(csumplq16(5,13)+csumplq16(6,13)+csumplq16(7,13)+csumplq16(8,13))&
        -(csumplq16(9,13)+csumplq16(10,13)+csumplq16(11,13)+csumplq16(12,13))&
        +(csumplq16(13,13)+csumplq16(14,13)+csumplq16(15,13)+csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 209)=(csumplqmom16(5,ik,13)+csumplqmom16(6,ik,13)&
        +csumplqmom16(7,ik,13)+csumplqmom16(8,ik,13)&
        -(csumplqmom16(13,ik,13)+csumplqmom16(14,ik,13)+csumplqmom16(15,ik,13)&
        +csumplqmom16(16,ik,13)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 210)=(csumplq16(1,13)+giot*csumplq16(2,13)-csumplq16(3,13)&
        -giot*csumplq16(4,13)+(csumplq16(5,13)+giot*csumplq16(6,13)-csumplq16(7,13)&
        -giot*csumplq16(8,13)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 210)=(csumplqmom16(1,ik,13)&
        +giot*csumplqmom16(2,ik,13)&
        -csumplqmom16(3,ik,13)-giot*csumplqmom16(4,ik,13))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 211)=(csumplq16(1,13)+giot*csumplq16(2,13)-csumplq16(3,13)&
        -giot*csumplq16(4,13)-(csumplq16(5,13)+giot*csumplq16(6,13)-csumplq16(7,13)&
        -giot*csumplq16(8,13)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 211)=(csumplqmom16(5,ik,13)&
        +giot*csumplqmom16(6,ik,13)&
        -csumplqmom16(7,ik,13)-giot*csumplqmom16(8,ik,13))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 212)=(csumplq16(1,13)-csumplq16(2,13)+csumplq16(3,13)&
        -csumplq16(4,13)&
        +(csumplq16(5,13)-csumplq16(6,13)+csumplq16(7,13)-csumplq16(8,13))&
        +(csumplq16(9,13)-csumplq16(10,13)+csumplq16(11,13)-csumplq16(12,13))&
        +(csumplq16(13,13)-csumplq16(14,13)+csumplq16(15,13)-csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 212)=(csumplqmom16(1,ik,13)-csumplqmom16(2,ik,13)&
        +csumplqmom16(3,ik,13)-csumplqmom16(4,ik,13)&
        +(csumplqmom16(9,ik,13)-csumplqmom16(10,ik,13)+csumplqmom16(11,ik,13)&
        -csumplqmom16(12,ik,13)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 213)=(csumplq16(1,13)-csumplq16(2,13)+csumplq16(3,13)&
        -csumplq16(4,13)&
        -(csumplq16(5,13)-csumplq16(6,13)+csumplq16(7,13)-csumplq16(8,13))&
        +(csumplq16(9,13)-csumplq16(10,13)+csumplq16(11,13)-csumplq16(12,13))&
        -(csumplq16(13,13)-csumplq16(14,13)+csumplq16(15,13)-csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 213)=(csumplqmom16(5,ik,13)-csumplqmom16(6,ik,13)&
        +csumplqmom16(7,ik,13)-csumplqmom16(8,ik,13)&
        +(csumplqmom16(13,ik,13)-csumplqmom16(14,ik,13)+csumplqmom16(15,ik,13)&
        -csumplqmom16(16,ik,13)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 214)=(csumplq16(1,13)-csumplq16(2,13)+csumplq16(3,13)&
        -csumplq16(4,13)&
        +(csumplq16(5,13)-csumplq16(6,13)+csumplq16(7,13)-csumplq16(8,13))&
        -(csumplq16(9,13)-csumplq16(10,13)+csumplq16(11,13)-csumplq16(12,13))&
        -(csumplq16(13,13)-csumplq16(14,13)+csumplq16(15,13)-csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 214)=(csumplqmom16(1,ik,13)-csumplqmom16(2,ik,13)&
        +csumplqmom16(3,ik,13)-csumplqmom16(4,ik,13)&
        -(csumplqmom16(9,ik,13)-csumplqmom16(10,ik,13)+csumplqmom16(11,ik,13)&
        -csumplqmom16(12,ik,13)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 215)=(csumplq16(1,13)-csumplq16(2,13)+csumplq16(3,13)&
        -csumplq16(4,13)&
        -(csumplq16(5,13)-csumplq16(6,13)+csumplq16(7,13)-csumplq16(8,13))&
        -(csumplq16(9,13)-csumplq16(10,13)+csumplq16(11,13)-csumplq16(12,13))&
        +(csumplq16(13,13)-csumplq16(14,13)+csumplq16(15,13)-csumplq16(16,13)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 215)=(csumplqmom16(5,ik,13)-csumplqmom16(6,ik,13)&
        +csumplqmom16(7,ik,13)-csumplqmom16(8,ik,13)&
        -(csumplqmom16(13,ik,13)-csumplqmom16(14,ik,13)+csumplqmom16(15,ik,13)&
        -csumplqmom16(16,ik,13)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 14
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 216)=(csumplq16(1,14)+csumplq16(2,14)+csumplq16(3,14)&
        +csumplq16(4,14)&
        +(csumplq16(5,14)+csumplq16(6,14)+csumplq16(7,14)+csumplq16(8,14))&
        +(csumplq16(9,14)+csumplq16(10,14)+csumplq16(11,14)+csumplq16(12,14))&
        +(csumplq16(13,14)+csumplq16(14,14)+csumplq16(15,14)+csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 216)=(csumplqmom16(1,ik,14)+csumplqmom16(2,ik,14)&
        +csumplqmom16(3,ik,14)+csumplqmom16(4,ik,14)&
        +(csumplqmom16(9,ik,14)+csumplqmom16(10,ik,14)+csumplqmom16(11,ik,14)&
        +csumplqmom16(12,ik,14)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 217)=(csumplq16(1,14)+csumplq16(2,14)+csumplq16(3,14)&
        +csumplq16(4,14)&
        -(csumplq16(5,14)+csumplq16(6,14)+csumplq16(7,14)+csumplq16(8,14))&
        +(csumplq16(9,14)+csumplq16(10,14)+csumplq16(11,14)+csumplq16(12,14))&
        -(csumplq16(13,14)+csumplq16(14,14)+csumplq16(15,14)+csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 217)=(csumplqmom16(5,ik,14)+csumplqmom16(6,ik,14)&
        +csumplqmom16(7,ik,14)+csumplqmom16(8,ik,14)&
        +(csumplqmom16(13,ik,14)+csumplqmom16(14,ik,14)+csumplqmom16(15,ik,14)&
        +csumplqmom16(16,ik,14)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 218)=(csumplq16(1,14)+csumplq16(2,14)+csumplq16(3,14)&
        +csumplq16(4,14)&
        +(csumplq16(5,14)+csumplq16(6,14)+csumplq16(7,14)+csumplq16(8,14))&
        -(csumplq16(9,14)+csumplq16(10,14)+csumplq16(11,14)+csumplq16(12,14))&
        -(csumplq16(13,14)+csumplq16(14,14)+csumplq16(15,14)+csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 218)=(csumplqmom16(1,ik,14)+csumplqmom16(2,ik,14)&
        +csumplqmom16(3,ik,14)+csumplqmom16(4,ik,14)&
        -(csumplqmom16(9,ik,14)+csumplqmom16(10,ik,14)+csumplqmom16(11,ik,14)&
        +csumplqmom16(12,ik,14)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 219)=(csumplq16(1,14)+csumplq16(2,14)+csumplq16(3,14)&
        +csumplq16(4,14)&
        -(csumplq16(5,14)+csumplq16(6,14)+csumplq16(7,14)+csumplq16(8,14))&
        -(csumplq16(9,14)+csumplq16(10,14)+csumplq16(11,14)+csumplq16(12,14))&
        +(csumplq16(13,14)+csumplq16(14,14)+csumplq16(15,14)+csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 219)=(csumplqmom16(5,ik,14)+csumplqmom16(6,ik,14)&
        +csumplqmom16(7,ik,14)+csumplqmom16(8,ik,14)&
        -(csumplqmom16(13,ik,14)+csumplqmom16(14,ik,14)+csumplqmom16(15,ik,14)&
        +csumplqmom16(16,ik,14)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 220)=(csumplq16(1,14)+giot*csumplq16(2,14)-csumplq16(3,14)&
        -giot*csumplq16(4,14)+(csumplq16(5,14)+giot*csumplq16(6,14)-csumplq16(7,14)&
        -giot*csumplq16(8,14)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 220)=(csumplqmom16(1,ik,14)&
        +giot*csumplqmom16(2,ik,14)&
        -csumplqmom16(3,ik,14)-giot*csumplqmom16(4,ik,14))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 221)=(csumplq16(1,14)+giot*csumplq16(2,14)-csumplq16(3,14)&
        -giot*csumplq16(4,14)-(csumplq16(5,14)+giot*csumplq16(6,14)-csumplq16(7,14)&
        -giot*csumplq16(8,14)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 221)=(csumplqmom16(5,ik,14)&
        +giot*csumplqmom16(6,ik,14)&
        -csumplqmom16(7,ik,14)-giot*csumplqmom16(8,ik,14))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 222)=(csumplq16(1,14)-csumplq16(2,14)+csumplq16(3,14)&
        -csumplq16(4,14)&
        +(csumplq16(5,14)-csumplq16(6,14)+csumplq16(7,14)-csumplq16(8,14))&
        +(csumplq16(9,14)-csumplq16(10,14)+csumplq16(11,14)-csumplq16(12,14))&
        +(csumplq16(13,14)-csumplq16(14,14)+csumplq16(15,14)-csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 222)=(csumplqmom16(1,ik,14)-csumplqmom16(2,ik,14)&
        +csumplqmom16(3,ik,14)-csumplqmom16(4,ik,14)&
        +(csumplqmom16(9,ik,14)-csumplqmom16(10,ik,14)+csumplqmom16(11,ik,14)&
        -csumplqmom16(12,ik,14)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 223)=(csumplq16(1,14)-csumplq16(2,14)+csumplq16(3,14)&
        -csumplq16(4,14)&
        -(csumplq16(5,14)-csumplq16(6,14)+csumplq16(7,14)-csumplq16(8,14))&
        +(csumplq16(9,14)-csumplq16(10,14)+csumplq16(11,14)-csumplq16(12,14))&
        -(csumplq16(13,14)-csumplq16(14,14)+csumplq16(15,14)-csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 223)=(csumplqmom16(5,ik,14)-csumplqmom16(6,ik,14)&
        +csumplqmom16(7,ik,14)-csumplqmom16(8,ik,14)&
        +(csumplqmom16(13,ik,14)-csumplqmom16(14,ik,14)+csumplqmom16(15,ik,14)&
        -csumplqmom16(16,ik,14)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 224)=(csumplq16(1,14)-csumplq16(2,14)+csumplq16(3,14)&
        -csumplq16(4,14)&
        +(csumplq16(5,14)-csumplq16(6,14)+csumplq16(7,14)-csumplq16(8,14))&
        -(csumplq16(9,14)-csumplq16(10,14)+csumplq16(11,14)-csumplq16(12,14))&
        -(csumplq16(13,14)-csumplq16(14,14)+csumplq16(15,14)-csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 224)=(csumplqmom16(1,ik,14)-csumplqmom16(2,ik,14)&
        +csumplqmom16(3,ik,14)-csumplqmom16(4,ik,14)&
        -(csumplqmom16(9,ik,14)-csumplqmom16(10,ik,14)+csumplqmom16(11,ik,14)&
        -csumplqmom16(12,ik,14)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 225)=(csumplq16(1,14)-csumplq16(2,14)+csumplq16(3,14)&
        -csumplq16(4,14)&
        -(csumplq16(5,14)-csumplq16(6,14)+csumplq16(7,14)-csumplq16(8,14))&
        -(csumplq16(9,14)-csumplq16(10,14)+csumplq16(11,14)-csumplq16(12,14))&
        +(csumplq16(13,14)-csumplq16(14,14)+csumplq16(15,14)-csumplq16(16,14)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 225)=(csumplqmom16(5,ik,14)-csumplqmom16(6,ik,14)&
        +csumplqmom16(7,ik,14)-csumplqmom16(8,ik,14)&
        -(csumplqmom16(13,ik,14)-csumplqmom16(14,ik,14)+csumplqmom16(15,ik,14)&
        -csumplqmom16(16,ik,14)))*adiv2
        enddo
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     plaquette operators 15
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !**********************************************************************
    !     j=0, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 226)=(csumplq16(1,15)+csumplq16(2,15)+csumplq16(3,15)&
        +csumplq16(4,15)&
        +(csumplq16(5,15)+csumplq16(6,15)+csumplq16(7,15)+csumplq16(8,15))&
        +(csumplq16(9,15)+csumplq16(10,15)+csumplq16(11,15)+csumplq16(12,15))&
        +(csumplq16(13,15)+csumplq16(14,15)+csumplq16(15,15)+csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 226)=(csumplqmom16(1,ik,15)+csumplqmom16(2,ik,15)&
        +csumplqmom16(3,ik,15)+csumplqmom16(4,ik,15)&
        +(csumplqmom16(9,ik,15)+csumplqmom16(10,ik,15)+csumplqmom16(11,ik,15)&
        +csumplqmom16(12,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 227)=(csumplq16(1,15)+csumplq16(2,15)+csumplq16(3,15)&
        +csumplq16(4,15)&
        -(csumplq16(5,15)+csumplq16(6,15)+csumplq16(7,15)+csumplq16(8,15))&
        +(csumplq16(9,15)+csumplq16(10,15)+csumplq16(11,15)+csumplq16(12,15))&
        -(csumplq16(13,15)+csumplq16(14,15)+csumplq16(15,15)+csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=0, pp=+, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 227)=(csumplqmom16(5,ik,15)+csumplqmom16(6,ik,15)&
        +csumplqmom16(7,ik,15)+csumplqmom16(8,ik,15)&
        +(csumplqmom16(13,ik,15)+csumplqmom16(14,ik,15)+csumplqmom16(15,ik,15)&
        +csumplqmom16(16,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 228)=(csumplq16(1,15)+csumplq16(2,15)+csumplq16(3,15)&
        +csumplq16(4,15)&
        +(csumplq16(5,15)+csumplq16(6,15)+csumplq16(7,15)+csumplq16(8,15))&
        -(csumplq16(9,15)+csumplq16(10,15)+csumplq16(11,15)+csumplq16(12,15))&
        -(csumplq16(13,15)+csumplq16(14,15)+csumplq16(15,15)+csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 228)=(csumplqmom16(1,ik,15)+csumplqmom16(2,ik,15)&
        +csumplqmom16(3,ik,15)+csumplqmom16(4,ik,15)&
        -(csumplqmom16(9,ik,15)+csumplqmom16(10,ik,15)+csumplqmom16(11,ik,15)&
        +csumplqmom16(12,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !     j=0, pp=-, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 229)=(csumplq16(1,15)+csumplq16(2,15)+csumplq16(3,15)&
        +csumplq16(4,15)&
        -(csumplq16(5,15)+csumplq16(6,15)+csumplq16(7,15)+csumplq16(8,15))&
        -(csumplq16(9,15)+csumplq16(10,15)+csumplq16(11,15)+csumplq16(12,15))&
        +(csumplq16(13,15)+csumplq16(14,15)+csumplq16(15,15)+csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=0, pp=-, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 229)=(csumplqmom16(5,ik,15)+csumplqmom16(6,ik,15)&
        +csumplqmom16(7,ik,15)+csumplqmom16(8,ik,15)&
        -(csumplqmom16(13,ik,15)+csumplqmom16(14,ik,15)+csumplqmom16(15,ik,15)&
        +csumplqmom16(16,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !     j=1, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 230)=(csumplq16(1,15)+giot*csumplq16(2,15)-csumplq16(3,15)&
        -giot*csumplq16(4,15)+(csumplq16(5,15)+giot*csumplq16(6,15)-csumplq16(7,15)&
        -giot*csumplq16(8,15)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 230)=(csumplqmom16(1,ik,15)&
        +giot*csumplqmom16(2,ik,15)&
        -csumplqmom16(3,ik,15)-giot*csumplqmom16(4,ik,15))*adiv1
        enddo
    !**********************************************************************
    !     j=1, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 231)=(csumplq16(1,15)+giot*csumplq16(2,15)-csumplq16(3,15)&
        -giot*csumplq16(4,15)-(csumplq16(5,15)+giot*csumplq16(6,15)-csumplq16(7,15)&
        -giot*csumplq16(8,15)))*adiv2
    !**********************************************************************
    !     j=1, q=1,2  here!
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 231)=(csumplqmom16(5,ik,15)&
        +giot*csumplqmom16(6,ik,15)&
        -csumplqmom16(7,ik,15)-giot*csumplqmom16(8,ik,15))*adiv1
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 232)=(csumplq16(1,15)-csumplq16(2,15)+csumplq16(3,15)&
        -csumplq16(4,15)&
        +(csumplq16(5,15)-csumplq16(6,15)+csumplq16(7,15)-csumplq16(8,15))&
        +(csumplq16(9,15)-csumplq16(10,15)+csumplq16(11,15)-csumplq16(12,15))&
        +(csumplq16(13,15)-csumplq16(14,15)+csumplq16(15,15)-csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 232)=(csumplqmom16(1,ik,15)-csumplqmom16(2,ik,15)&
        +csumplqmom16(3,ik,15)-csumplqmom16(4,ik,15)&
        +(csumplqmom16(9,ik,15)-csumplqmom16(10,ik,15)+csumplqmom16(11,ik,15)&
        -csumplqmom16(12,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=+, pr=-, q=0
    !**********************************************************************
        lines(time_slice, id, 233)=(csumplq16(1,15)-csumplq16(2,15)+csumplq16(3,15)&
        -csumplq16(4,15)&
        -(csumplq16(5,15)-csumplq16(6,15)+csumplq16(7,15)-csumplq16(8,15))&
        +(csumplq16(9,15)-csumplq16(10,15)+csumplq16(11,15)-csumplq16(12,15))&
        -(csumplq16(13,15)-csumplq16(14,15)+csumplq16(15,15)-csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=2, pp=+, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 233)=(csumplqmom16(5,ik,15)-csumplqmom16(6,ik,15)&
        +csumplqmom16(7,ik,15)-csumplqmom16(8,ik,15)&
        +(csumplqmom16(13,ik,15)-csumplqmom16(14,ik,15)+csumplqmom16(15,ik,15)&
        -csumplqmom16(16,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=+, q=0
    !**********************************************************************
        lines(time_slice, id, 234)=(csumplq16(1,15)-csumplq16(2,15)+csumplq16(3,15)&
        -csumplq16(4,15)&
        +(csumplq16(5,15)-csumplq16(6,15)+csumplq16(7,15)-csumplq16(8,15))&
        -(csumplq16(9,15)-csumplq16(10,15)+csumplq16(11,15)-csumplq16(12,15))&
        -(csumplq16(13,15)-csumplq16(14,15)+csumplq16(15,15)-csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 234)=(csumplqmom16(1,ik,15)-csumplqmom16(2,ik,15)&
        +csumplqmom16(3,ik,15)-csumplqmom16(4,ik,15)&
        -(csumplqmom16(9,ik,15)-csumplqmom16(10,ik,15)+csumplqmom16(11,ik,15)&
        -csumplqmom16(12,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !     j=2, pp=-, pr=-, q=0  test
    !**********************************************************************
        lines(time_slice, id, 235)=(csumplq16(1,15)-csumplq16(2,15)+csumplq16(3,15)&
        -csumplq16(4,15)&
        -(csumplq16(5,15)-csumplq16(6,15)+csumplq16(7,15)-csumplq16(8,15))&
        -(csumplq16(9,15)-csumplq16(10,15)+csumplq16(11,15)-csumplq16(12,15))&
        +(csumplq16(13,15)-csumplq16(14,15)+csumplq16(15,15)-csumplq16(16,15)))*adiv3
    !**********************************************************************
    !     j=2, pp=-, q=0
    !**********************************************************************
        do ik=1, 2
            momentum_lines(time_slice, id, ik, 235)=(csumplqmom16(5,ik,15)-csumplqmom16(6,ik,15)&
        +csumplqmom16(7,ik,15)-csumplqmom16(8,ik,15)&
        -(csumplqmom16(13,ik,15)-csumplqmom16(14,ik,15)+csumplqmom16(15,ik,15)&
        -csumplqmom16(16,ik,15)))*adiv2
        enddo
    !**********************************************************************
    !**********************************************************************
        return
        end
    !**********************************************************************
end module
