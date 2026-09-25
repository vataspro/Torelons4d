module here_be_dragons
    use parameters
    use lattice
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
        complex(real64) :: act(3), ast(3)
        complex(real64) :: A11(NCOL, NCOL), B11(NCOL, NCOL), C11(NCOL, NCOL), &
        D11(NCOL, NCOL), E11(NCOL, NCOL), F11(NCOL, NCOL), UINT11(NCOL, NCOL), REM11(NCOL, NCOL)

        ! Loop variables
        integer :: idd, ib, idg, ijn, iloop, i, ic, nc, ig, mu

        ! Site indices
        integer :: site, site_ku, site_ju, site_iu, mn, ml, m1, m2, m3, m4, m5

        ! Direction indices
        integer :: iu, ju, ku

        ! Blocking level variables
        integer :: id, lrest, ids, idsm1, irem, li, ico

        ! Operator variables
        integer :: iddd, ieee

        ! Momentum variables
        real(real64) :: phase

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

        ! G11
        complex(real64) :: g11(NCOL, NCOL)

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
                            csumud(4), csumrc(4), csumrcw(4)
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
    
        complex(real64) :: giot, pf(2), uren11(NCOL, NCOL), aa(NCOL, NCOL), det, cnorm, cdet, dncol

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

        site = 0 ! Initial lattice point
        do site_ku = 1, ls(ku) ! LX direction
            ! Assign momentum phases
            phase = 2.0 * PI * site_ku / real(ls(ku))
            pf(1) = cmplx(cos(phase), sin(phase), kind=real64)
            phase = phase * 2.0d0
            PF(2)=cmplx(cos(phase), sin(phase), kind=real64)

            do site_ju = 1, ls(ju) ! LY direction
                do site_iu = 1, ls(iu) ! LZ direction 
                ! Lexicographical lattice site definition
                site = site + 1
                ix(iu) = site_iu
                ix(ju) = site_ju
                ix(ku) = site_ku
                ! MN defines the current lattice site
                mn = ix(1) + ls(1)*(ix(2)-1) + ls(1)*ls(2)*(ix(3)-1)


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
                do 9 iddd = 1, 337 !new!
                    m2 = mn
                    do ic = 1,NCOL
                        A11(ic, ic) = cmplx(1.0, 0.0, kind=real64)
                    enddo

                    if (ids == id) then
                        ico=2
                        if (iddd < 5) ico=1
                        if (((iddd > 12).and.(iddd < 17)).and.(ids /= 1)) ico=1
                        if (((iddd > 16).and.(iddd < 21)).and.(ids == 1)) ico=4
                        if (((iddd > 20).and.(iddd < 25)).and.(ids == 1)) ico=4
                        if (((iddd > 24).and.(iddd < 29)).and.(ids == 1)) ico=4
                        if (((iddd > 28).and.(iddd < 33)).and.(ids == 1)) ico=4
                        if (((iddd > 40).and.(iddd < 49)).and.(ids == 1)) ico=4
                        if (((iddd > 48).and.(iddd < 57)).and.(ids == 1)) ico=4
                        if (((iddd > 56).and.(iddd < 65)).and.(ids == 1)) ico=4
                        if (((iddd > 64).and.(iddd < 69)).and.(ids == 1)) ico=4
                        if (((iddd > 68).and.(iddd < 73)).and.(ids == 1)) ico=4
                        if (((iddd > 72).and.(iddd < 81)).and.(ids == 1)) ico=4
                        if (((iddd > 80).and.(iddd < 89)).and.(ids == 1)) ico=4
                        if (((iddd > 88).and.(iddd < 97)).and.(ids == 1)) ico=4
                        if (((iddd > 96).and.(iddd < 105)).and.(ids == 1)) ico=4
                        if (((iddd > 104).and.(iddd < 113)).and.(ids == 1)) ico=4
                        if (((iddd > 112).and.(iddd < 129)).and.(ids == 1)) ico=4
                        if (((iddd > 128).and.(iddd < 137)).and.(ids == 1)) ico=4
                        if (((iddd > 136).and.(iddd < 145)).and.(ids /= 1)) ico=1
                        if (iddd == 145) ico=1
                        if ((iddd > 145).and.(iddd < 194)) ico=1
                        if ((iddd > 209).and.(iddd < 338)) ico=1 ! THIS NEEDS TO BE FIXED  !

                        if (lcnt(ids) >= ico) then
                            !**********************************************************************C
                            !                     UP SQUARE PULSE                                  C 
                            !**********************************************************************C
                            !                     UP Y
                            !**********************************************************************
                            if (iddd == 1) then
                                do IC=1,NCOL2
                                B11(IC)=UC11(IC,M2,JU)
                                enddo
                                M3=IUP(M2,JU)
                                do IC=1,NCOL2
                                C11(IC)=UC11(IC,M3,KU)
                                enddo
                                CALL VMX(1,B11,C11,D11,1)
                                M3=IUP(M2,KU)
                                M2=M3
                                do IC=1,NCOL2
                                C11(IC)=UC11(IC,M2,JU)
                                enddo
                                CALL HERM(1,C11,DUM11,1)
                                CALL VMX(1,D11,C11,SQUY1,1)
                                CALL VMX(1,A11,SQUY1,C11,1)
       
                                do IC=1,NCOL2
                                A11(IC)=C11(IC)
                                enddo
         
                                IEEE=1
         
                            endif 
                            !**********************************************************************C
                            !     UP Z                                                             C 
                            !**********************************************************************C
                            if (iddd == 2) then
         
                            do 16 IC=1,NCOL2
                                B11(IC)=UC11(IC,M2,IU)
    16                        continue
                            M3=IUP(M2,IU)
                            do 17 IC=1,NCOL2
                                C11(IC)=UC11(IC,M3,KU)
    17                        continue
                            CALL VMX(1,B11,C11,D11,1)
                            M3=IUP(M2,KU)
                            M2=M3
                            do 18 IC=1,NCOL2
                                C11(IC)=UC11(IC,M2,IU)
    18                        continue
                            CALL HERM(1,C11,DUM11,1)
                            CALL VMX(1,D11,C11,SQUZ1,1)
                            CALL VMX(1,A11,SQUZ1,C11,1)
    !     
                            do 19 IC=1,NCOL2
                                A11(IC)=C11(IC)
    19                        continue
    !     
                            IEEE=2
    !     
                            endif 
    !**********************************************************************C 
    !                       doWN Y                                         C
    !**********************************************************************C
                            if (iddd == 3) then
    !     
                            M3=IDN(M2,JU)
                            do 20 IC=1,NCOL2
                                D11(IC)=UC11(IC,M3,JU)
    20                        continue
                            CALL HERM(1,D11,DUM11,1)
                            do 21 IC=1,NCOL2
                                C11(IC)=UC11(IC,M3,KU)
    21                        continue
                            CALL VMX(1,D11,C11,B11,1)
                            M4=IUP(M3,KU)
                            do 22 IC=1,NCOL2
                                D11(IC)=UC11(IC,M4,JU)
    22                        continue
                            CALL VMX(1,B11,D11,SQDY1,1)
                            CALL VMX(1,A11,SQDY1,C11,1)
    !     
                            do 23 IC=1,NCOL2
                                A11(IC)=C11(IC)
    23                        continue
    !     
                            IEEE=3
    !     
                            endif 
    !**********************************************************************C 
    !                       doWN Z                                 C
    !**********************************************************************C
                            if (iddd == 4) then
    !     
                            M3=IDN(M2,IU)
                            do 24 IC=1,NCOL2
                                D11(IC)=UC11(IC,M3,IU)
    24                        continue
                            CALL HERM(1,D11,DUM11,1)
                            do 25 IC=1,NCOL2
                                C11(IC)=UC11(IC,M3,KU)
    25                        continue
                            CALL VMX(1,D11,C11,B11,1)
                            M4=IUP(M3,KU)
                            do 26 IC=1,NCOL2
                                D11(IC)=UC11(IC,M4,IU)
    26                        continue
                            CALL VMX(1,B11,D11,SQDZ1,1)
                            CALL VMX(1,A11,SQDZ1,C11,1)
    !     
                            do 27 IC=1,NCOL2
                                A11(IC)=C11(IC)
    27                        continue
    !     
                            IEEE=4
    !     
                            endif 
    !**********************************************************************C
    !                       UP - UP SQUARE PULSES                         *C 
    !**********************************************************************C
    !                       UP Y
    !**********************************************************************C
                            if (iddd == 5) then
    !     
                            M3=IUP(M2,KU)
                            do 28 IC=1,NCOL2
                                D11(IC)=UC11(IC,M3,JU)
    28                        continue
                            M4=IUP(M3,JU)
                            do 29 IC=1,NCOL2
                                C11(IC)=UC11(IC,M4,KU)
    29                        continue
                            CALL VMX(1,D11,C11,B11,1)
                            M2=IUP(M3,KU)
                            do 30 IC=1,NCOL2
                                C11(IC)=UC11(IC,M2,JU)
    30                        continue
                            CALL HERM(1,C11,DUM11,1)
                            CALL VMX(1,B11,C11,SQUY2,1)
                            CALL VMX(1,SQUY1,SQUY2,A11,1)
    !     
                            IEEE=1
    !     
                            endif 
    !**********************************************************************C
    !              UP Z
    !**********************************************************************C
                            if (iddd == 6) then
    !     
                            M3=IUP(M2,KU)
                            do 31 IC=1,NCOL2
                                D11(IC)=UC11(IC,M3,IU)
    31                        continue
                            M4=IUP(M3,IU)
                            do 32 IC=1,NCOL2
                                C11(IC)=UC11(IC,M4,KU)
    32                        continue
                            CALL VMX(1,D11,C11,B11,1)
                            M2=IUP(M3,KU)
                            do 33 IC=1,NCOL2
                                C11(IC)=UC11(IC,M2,IU)
    33                        continue
                            CALL HERM(1,C11,DUM11,1)
                            CALL VMX(1,B11,C11,SQUZ2,1)
                            CALL VMX(1,SQUZ1,SQUZ2,A11,1)
    !     
                            IEEE=2
    !     
                            endif 
    !**********************************************************************C 
    !              doWN Y                                              C
    !**********************************************************************C
                            if (iddd == 7) then
    !     
                            M3=IUP(M2,KU)
                            M4=IDN(M3,JU)
                            do 34 IC=1,NCOL2
                                D11(IC)=UC11(IC,M4,JU)
    34                        continue
                            CALL HERM(1,D11,DUM11,1)
                            do 35 IC=1,NCOL2
                                C11(IC)=UC11(IC,M4,KU)
    35                        continue
                            CALL VMX(1,D11,C11,B11,1)
                            M1=IUP(M4,KU)
                            do 36 IC=1,NCOL2
                                C11(IC)=UC11(IC,M1,JU)
    36                        continue
                            CALL VMX(1,B11,C11,SQDY2,1)
                            CALL VMX(1,SQDY1,SQDY2,A11,1)
                            M2=IUP(M3,KU)
    !     
                            IEEE=3
    !     
                            endif 
    !**********************************************************************C 
    !              doWN Z                                              C
    !**********************************************************************C
                            if (iddd == 8) then
    !     
                            M3=IUP(M2,KU)
                            M4=IDN(M3,IU)
                            do 37 IC=1,NCOL2
                                D11(IC)=UC11(IC,M4,IU)
    37                        continue
                            CALL HERM(1,D11,DUM11,1)
                            do 38 IC=1,NCOL2
                                C11(IC)=UC11(IC,M4,KU)
    38                        continue
                            CALL VMX(1,D11,C11,B11,1)
                            M1=IUP(M4,KU)
                            do 39 IC=1,NCOL2
                                C11(IC)=UC11(IC,M1,IU)
    39                        continue
                            CALL VMX(1,B11,C11,SQDZ2,1)
                            CALL VMX(1,SQDZ1,SQDZ2,A11,1)
                            M2=IUP(M3,KU)
    !     
                            IEEE=4
    !     
                            endif 
    !**********************************************************************C
    !**********************************************************************C
    !              UP - doWN SQUARE PULSES                                 C 
    !**********************************************************************C
    !**********************************************************************C
    !                 UP Y 
    !**********************************************************************C
                            if (iddd == 9) then
                            CALL VMX(1,SQUY1,SQDY2,SQUDY1,1)
                            do IC=1, NCOL2
                                A11(IC)=SQUDY1(IC) !NEW!
                            enddo
                            M3=IUP(M2,KU)
                            M2=IUP(M3,KU)
    !     
                            IEEE=1
    !     
                            endif 
    !**********************************************************************C
    !                 UP Z
    !**********************************************************************C
                            if (iddd == 10) then    
                            CALL VMX(1,SQUZ1,SQDZ2,SQUDZ1,1)
                            do IC=1, NCOL2
                                A11(IC)=SQUDZ1(IC) !NEW!
                            enddo
                            M3=IUP(M2,KU)
                            M2=IUP(M3,KU)
    !     
                            IEEE=2
    !     
                            endif 
    !**********************************************************************C 
    !                 doWN Y
    !**********************************************************************C
                            if (iddd == 11) then
                            CALL VMX(1,SQDY1,SQUY2,SQDUY1,1)
                            do IC=1, NCOL2
                                A11(IC)=SQDUY1(IC) !NEW!
                            enddo
                            M3=IUP(M2,KU)
                            M2=IUP(M3,KU)
    !     
                            IEEE=3
    !     
                            endif 
    !**********************************************************************C 
    !                 doWN Z
    !**********************************************************************C
                            if (iddd == 12) then
                            CALL VMX(1,SQDZ1,SQUZ2,SQDUZ1,1)
                            do IC=1, NCOL2
                                A11(IC)=SQDUZ1(IC) !NEW!
                            enddo
                            M3=IUP(M2,KU)
                            M2=IUP(M3,KU)
    !
                            IEEE=4
    !     
                            endif 
!**********************************************************************C
!**********************************************************************C
!                 UP WAVE - LIKE PULSE                                *C
!**********************************************************************C
!**********************************************************************C
!                 UP Y
!**********************************************************************C
                if (iddd == 13) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!     
                    do 40 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,JU,idsW)
40                  continue
                    M3=IUPB(M2,JU,idsW)
                    do 41 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
41                  continue
                    CALL VMX(1,D11,C11,B11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 42 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,JU,idsW)
42                  continue
                    CALL HERM(1,C11,DUM11,1)
                    CALL VMX(1,B11,C11,WSQUY1,1)
!                     
                    M3=IDNB(M2,JU,idsW)
                    do 43 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,JU,idsW)
43                  continue
                    CALL HERM(1,B11,DUM11,1)
                    do 44 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
44                  continue
                    CALL VMX(1,B11,C11,D11,1)
                    M4=IUPB(M3,KU,idsW)
                    do 45 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,JU,idsW)
45                  continue
                    CALL VMX(1,D11,B11,WSQDY2,1)
                    CALL VMX(1,WSQUY1,WSQDY2,WVUY1,1)
                    do 46 IC=1,NCOL2
                        A11(IC)=WVUY1(IC)
46                  continue
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=1
!
                endif 
!**********************************************************************C
!                 UP Z
!**********************************************************************C
                if (iddd == 14) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!     
                    do 47 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,IU,idsW)
47                  continue
                    M3=IUPB(M2,IU,idsW)
                    do 48 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
48                  continue
                    CALL VMX(1,D11,C11,B11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 49 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,IU,idsW)
49                  continue
                    CALL HERM(1,C11,DUM11,1)
                    CALL VMX(1,B11,C11,WSQUZ1,1)
!                     
                    M3=IDNB(M2,IU,idsW)
                    do 50 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,IU,idsW)
50                  continue
                    CALL HERM(1,B11,DUM11,1)
                    do 51 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
51                  continue
                    CALL VMX(1,B11,C11,D11,1)
                    M4=IUPB(M3,KU,idsW)
                    do 52 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,IU,idsW)
52                  continue
                    CALL VMX(1,D11,B11,WSQDZ2,1)
                    CALL VMX(1,WSQUZ1,WSQDZ2,WVUZ1,1)
                    do 53 IC=1,NCOL2
                        A11(IC)=WVUZ1(IC)
53                  continue
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=2
!
                endif 
!**********************************************************************C
!                 doWN Y                                              *C 
!**********************************************************************C
                if (iddd == 15) then
!     
                        idsW=ids-1
!
                        if (idsW == 0) then
                        idsW=1
                        endif 
!     
                        M3=IDNB(M2,JU,idsW)
                        do 54 IC=1,NCOL2
                        D11(IC)=UB11(IC,M3,JU,idsW)
54                     continue
                        CALL HERM(1,D11,DUM11,1)
                        do 55 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
55                     continue
                        CALL VMX(1,D11,C11,B11,1)
                        M4=IUPB(M3,KU,idsW)
                        do 56 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,JU,idsW)
56                     continue
                        CALL VMX(1,B11,C11,WSQDY1,1)
                        M3=IUPB(M2,KU,idsW)
                        M2=M3
!     
                        do 57 IC=1,NCOL2
                        B11(IC)=UB11(IC,M2,JU,idsW)
57                     continue
                        M3=IUPB(M2,JU,idsW)
                        do 58 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
58                     continue
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IUPB(M2,KU,idsW)
                        M2=M3
                        do 59 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,JU,idsW)
59                     continue
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,D11,C11,WSQUY2,1)
                        CALL VMX(1,WSQDY1,WSQUY2,WVDY1,1)
!     
                        do 60 IC=1, NCOL2
                        A11(IC)=WVDY1(IC)
60                     continue
!
                        IEEE=3
! 
                endif 
!**********************************************************************C
!                 doWN Z                                              *C 
!**********************************************************************C
                if (iddd == 16) then
!     
                        idsW=ids-1
!
                        if (idsW == 0) then
                        idsW=1
                        endif 
!     
                        M3=IDNB(M2,IU,idsW)
                        do 61 IC=1,NCOL2
                        D11(IC)=UB11(IC,M3,IU,idsW)
61                     continue
                        CALL HERM(1,D11,DUM11,1)
                        do 62 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
62                     continue
                        CALL VMX(1,D11,C11,B11,1)
                        M4=IUPB(M3,KU,idsW)
                        do 63 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,IU,idsW)
63                     continue
                        CALL VMX(1,B11,C11,WSQDZ1,1)
                        M3=IUPB(M2,KU,idsW)
                        M2=M3
!     
                        do 64 IC=1,NCOL2
                        B11(IC)=UB11(IC,M2,IU,idsW)
64                     continue
                        M3=IUPB(M2,IU,idsW)
                        do 65 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
65                     continue
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IUPB(M2,KU,idsW)
                        M2=M3
                        do 66 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,IU,idsW)
66                     continue
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,D11,C11,WSQUZ2,1)
                        CALL VMX(1,WSQDZ1,WSQUZ2,WVDZ1,1)
!     
                        do 67 IC=1, NCOL2
                        A11(IC)=WVDZ1(IC)
67                     continue
!
                        IEEE=4
!     
                endif 
!**********************************************************************C
!**********************************************************************C
!              UP - UP WAVE-LIKE PULSE                                 C
!**********************************************************************C
!**********************************************************************C
!              UP Y
!**********************************************************************C
                if (iddd == 17) then
!
                    idsW=ids-1

                    if (idsW == 0) then 
                        idsW=1
                    endif 
! 
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
! 
                    do 68 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,JU,idsW)
68                  continue 
                    M3=IUPB(M2,JU,idsW)
                    do 69 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
69                  continue
                    CALL VMX(1,D11,C11,B11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 70 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,JU,idsW)
70                  continue
                    CALL HERM(1,C11,DUM11,1)
                    CALL VMX(1,B11,C11,WSQUY3,1)
!
                    M3=IDNB(M2,JU,idsW)
                    do 71 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,JU,idsW)
71                  continue
                    CALL HERM(1,B11,DUM11,1)
                    do 72 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
72                  continue
                    CALL VMX(1,B11,C11,D11,1)
                    M4=IUPB(M3,KU,idsW)
                    do 73 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,JU,idsW)
73                  continue
                    CALL VMX(1,D11,B11,WSQDY4,1)
                    CALL VMX(1,WSQUY3,WSQDY4,WVUY2,1)
                    CALL VMX(1,WVUY1,WVUY2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=1
!     
                endif 
!**********************************************************************C
!                 UP Z
!**********************************************************************C
                if (iddd == 18) then
!
                    idsW=ids-1

                    if (idsW == 0) then 
                        idsW=1
                    endif 
! 
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
! 
                    do 74 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,IU,idsW)
74                  continue 
                    M3=IUPB(M2,IU,idsW)
                    do 75 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
75                  continue
                    CALL VMX(1,D11,C11,B11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 76 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,IU,idsW)
76                  continue
                    CALL HERM(1,C11,DUM11,1)
                    CALL VMX(1,B11,C11,WSQUZ3,1)
!
                    M3=IDNB(M2,IU,idsW)
                    do 77 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,IU,idsW)
77                  continue
                    CALL HERM(1,B11,DUM11,1)
                    do 78 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
78                  continue
                    CALL VMX(1,B11,C11,D11,1)
                    M4=IUPB(M3,KU,idsW)
                    do 79 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,IU,idsW)
79                  continue
                    CALL VMX(1,D11,B11,WSQDZ4,1)
                    CALL VMX(1,WSQUZ3,WSQDZ4,WVUZ2,1)
                    CALL VMX(1,WVUZ1,WVUZ2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=2
!     
                endif 
!**********************************************************************C
!                 doWN Y                                               C
!**********************************************************************C
                if (iddd == 19) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!     
                    M3=IDNB(M2,JU,idsW)
                    do 80 IC=1,NCOL2
                        D11(IC)=UB11(IC,M3,JU,idsW)
80                  continue
                    CALL HERM(1,D11,DUM11,1)
                    do 81 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
81                  continue
                    CALL VMX(1,D11,C11,B11,1)
                    M4=IUPB(M3,KU,idsW)
                    do 82 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,JU,idsW)
82                  continue
                    CALL VMX(1,B11,C11,WSQDY3,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 83 IC=1,NCOL2
                        B11(IC)=UB11(IC,M2,JU,idsW)
83                  continue
                    M3=IUPB(M2,JU,idsW)
                    do 84 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
84                  continue
                    CALL VMX(1,B11,C11,D11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 85 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,JU,idsW)
85                  continue
                    CALL HERM(1,C11,DUM11,1)
                    CALL VMX(1,D11,C11,WSQUY4,1)
                    CALL VMX(1,WSQDY3,WSQUY4,WVDY2,1)
                    CALL VMX(1,WVDY1,WVDY2,A11,1)
!
                        IEEE=3
!     
                endif 
!**********************************************************************C
!                 doWN Y                                               C
!**********************************************************************C
                if (iddd == 20) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!     
                    M3=IDNB(M2,IU,idsW)
                    do 86 IC=1,NCOL2
                        D11(IC)=UB11(IC,M3,IU,idsW)
86                  continue
                    CALL HERM(1,D11,DUM11,1)
                    do 87 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
87                  continue
                    CALL VMX(1,D11,C11,B11,1)
                    M4=IUPB(M3,KU,idsW)
                    do 88 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,IU,idsW)
88                  continue
                    CALL VMX(1,B11,C11,WSQDZ3,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 89 IC=1,NCOL2
                        B11(IC)=UB11(IC,M2,IU,idsW)
89                  continue
                    M3=IUPB(M2,IU,idsW)
                    do 90 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,idsW)
90                  continue
                    CALL VMX(1,B11,C11,D11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    do 91 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,IU,idsW)
91                  continue
                    CALL HERM(1,C11,DUM11,1)
                    CALL VMX(1,D11,C11,WSQUZ4,1)
                    CALL VMX(1,WSQDZ3,WSQUZ4,WVDZ2,1)
                    CALL VMX(1,WVDZ1,WVDZ2,A11,1)
!
                        IEEE=4
!     
                endif 
!**********************************************************************C
!**********************************************************************C
!                 UP - doWN WAVE-LIKE PULSE                           *C
!**********************************************************************C
!**********************************************************************C
!                 UP Y                                                 C
!**********************************************************************C
                if (iddd == 21) then
!     
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!     
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVUY1,WVDY2,A11,1)
!
                        IEEE=1
!
                endif 
!**********************************************************************C
!                 UP Z                                                 C
!**********************************************************************C
                if (iddd == 22) then
!     
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!     
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVUZ1,WVDZ2,A11,1)
!
                        IEEE=2
!
                endif 
!**********************************************************************C
!                 doWN Y                                              *C
!**********************************************************************C
                if (iddd == 23) then
!     
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!     
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVDY1,WVUY2,A11,1)
!
                        IEEE=3
!     
                endif 
!**********************************************************************C
!                 doWN Z                                              *C
!**********************************************************************C
                if (iddd == 24) then
!     
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!     
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVDZ1,WVUZ2,A11,1)
!
                        IEEE=4
!     
                endif 
!**********************************************************************C
!**********************************************************************C
!                   /\_____/\  UP PULSE                                C
!**********************************************************************C
!**********************************************************************C
!                  UP Y
!**********************************************************************C
                if (iddd == 25) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 139 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,KU,idsW)
139                    continue
                        M4=IUPB(M3,KU,idsW)
                        do 140 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,KU,idsW)
140                    continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 141 IC=1, NCOL2
                        D11(IC)=UB11(IC,M3,KU,idsW+1)
141                    continue
                    endif 
!
                    CALL VMX(1,WSQUY1,D11,C11,1)
                    CALL VMX(1,C11,WSQUY4,A11,1)
!
!
                        IEEE=1
!
                endif 
!**********************************************************************C
!                  UP Z
!**********************************************************************C
                if (iddd == 26) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 142 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,KU,idsW)
142                    continue
                        M4=IUPB(M3,KU,idsW)
                        do 143 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,KU,idsW)
143                    continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 144 IC=1, NCOL2
                        D11(IC)=UB11(IC,M3,KU,idsW+1)
144                    continue
                    endif 
!
                    CALL VMX(1,WSQUZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQUZ4,A11,1)
!
                        IEEE=2
!
                endif 
!**********************************************************************C
!                   doWN Y                                             C
!**********************************************************************C
                if (iddd == 27) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 145 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
145                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 146 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
146                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 147 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
147                     continue
                    endif 
!     
                    CALL VMX(1,WSQDY1,D11,C11,1)
                    CALL VMX(1,C11,WSQDY4,A11,1)
!
                        IEEE=3
!
                endif 
!**********************************************************************C
!                   doWN Z                                             C
!**********************************************************************C
                if (iddd == 28) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 149 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
149                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 150 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
150                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 151 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
151                     continue
                    endif 
!     
                    CALL VMX(1,WSQDZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQDZ4,A11,1)
!
                        IEEE=4
!
                endif 
!**********************************************************************C
!**********************************************************************C
!                   /\-----\/  PULSE                                   C
!**********************************************************************C
!**********************************************************************C
!                  UP Y
!**********************************************************************C
                if (iddd == 29) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 152 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
152                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 153 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
153                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 154 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
154                     continue
                    endif 
!
                    CALL VMX(1,WSQUY1,D11,C11,1)
                    CALL VMX(1,C11,WSQDY4,A11,1)
!
                        IEEE=1
!
                endif 
!**********************************************************************C
!                  UP Z
!**********************************************************************C
                if (iddd == 30) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 155 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
155                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 156 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
156                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 157 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
157                     continue
                    endif 
!
                    CALL VMX(1,WSQUZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQDZ4,A11,1)
!
                        IEEE=2
!
                endif 
!**********************************************************************C
!                 doWN Y                                               C
!**********************************************************************C
                if (iddd == 31) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 158 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,KU,idsW)
158                    continue
                        M4=IUPB(M3,KU,idsW)
                        do 159 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,KU,idsW)  
159                    continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 160 IC=1, NCOL2
                        D11(IC)=UB11(IC,M3,KU,idsW+1)
160                    continue
                    endif 
!     
                    CALL VMX(1,WSQDY1,D11,C11,1)
                    CALL VMX(1,C11,WSQUY4,A11,1)
!
                        IEEE=3
!     
                endif 
!**********************************************************************C
!                 doWN Z                                               C
!**********************************************************************C
                if (iddd == 32) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 161 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,KU,idsW)
161                    continue
                        M4=IUPB(M3,KU,idsW)
                        do 162 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,KU,idsW)  
162                    continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 262 IC=1, NCOL2
                        D11(IC)=UB11(IC,M3,KU,idsW+1)
262                    continue
                    endif 
!     
                    CALL VMX(1,WSQDZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQUZ4,A11,1)
!
                        IEEE=4
!     
                endif 
!**********************************************************************
!     TT-1 OPERATORS
!**********************************************************************
                        if (iddd == 33)then
!     
                        CALL VMX(1,SQUY1,SQUZ2,A11,1)
!
                        do IC=1, NCOL2
                            WANGL1(IC)=A11(IC) !NEW!
                        enddo
!                           
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        IEEE=1
!     
                        endif 
!**********************************************************************
                        if (iddd == 34) then
!     
                        CALL VMX(1,SQUZ1,SQDY2,A11,1)
!
                        do IC=1, NCOL2
                            WANGL2(IC)=A11(IC) !NEW!
                        enddo
!                                                      
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        IEEE=2
!     
                        endif 
!**********************************************************************
                        if (iddd == 35) then
!     
                        CALL VMX(1,SQDY1,SQDZ2,A11,1)
!    
                        do IC=1, NCOL2
                            WANGL3(IC)=A11(IC) !NEW!
                        enddo
!                           
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        IEEE=3
!     
                        endif 
!**********************************************************************
                        if (iddd == 36) then
!     
                        CALL VMX(1,SQDZ1,SQUY2,A11,1)
!                           
                        do IC=1, NCOL2
                            WANGL4(IC)=A11(IC) !NEW!
                        enddo                           
!
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        IEEE=4
!     
                        endif 
!**********************************************************************
                        if (iddd == 37) then
!     
                        CALL VMX(1,SQUY1,SQDZ2,A11,1)
!                           
                        do IC=1, NCOL2
                            WANGL5(IC)=A11(IC) !NEW!
                        enddo                           
!
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        IEEE=5
!     
                        endif 
!**********************************************************************
                        if (iddd == 38) then
!     
                        CALL VMX(1,SQUZ1,SQUY2,A11,1)
!                           
                        do IC=1, NCOL2
                            WANGL6(IC)=A11(IC) !NEW!
                        enddo                           
!
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        IEEE=6
!     
                        endif 
!**********************************************************************
                        if (iddd == 39) then
!     
                        CALL VMX(1,SQDY1,SQUZ2,A11,1)
!                           
                        do IC=1, NCOL2
                            WANGL7(IC)=A11(IC) !NEW!
                        enddo                           
!                           
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        IEEE=7
!     
                        endif 
!**********************************************************************
                        if (iddd == 40) then
                        M3=IUP(M2,KU)
                        M2=IUP(M3,KU)
!     
                        CALL VMX(1,SQDZ1,SQDY2,A11,1)
!                           
                        do IC=1, NCOL2
                            WANGL8(IC)=A11(IC) !NEW!
                        enddo                                                      
!     
                        IEEE=8
!     
                        endif 
!**********************************************************************
!     TT-2 OPERATORS
!**********************************************************************
                if (iddd == 41) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVUY1,WVUZ2,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************
                if (iddd == 42) then
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVUZ1,WVDY2,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************
                if (iddd == 43) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVDY1,WVDZ2,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************
                if (iddd == 44) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVDZ1,WVUY2,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************
                if (iddd == 45) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVUZ1,WVUY2,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************
                if (iddd == 46) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVDY1,WVUZ2,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************
                if (iddd == 47) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVDZ1,WVDY2,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************
                if (iddd == 48) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                    endif 
!
                    M3=IUPB(M2,KU,idsW)
                    M4=IUPB(M3,KU,idsW)
                    M1=IUPB(M4,KU,idsW)
                    M2=IUPB(M1,KU,idsW)
!
                    CALL VMX(1,WVUY1,WVDZ2,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
!**********************************************************************C
!                   /\----/_/  UP PULSE-TT3                            C
!**********************************************************************C
!**********************************************************************C
!                  1
!**********************************************************************C
                if (iddd == 49) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 301 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,KU,idsW)
301                    continue
                        M4=IUPB(M3,KU,idsW)
                        do 302 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,KU,idsW)
302                    continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 303 IC=1, NCOL2
                        D11(IC)=UB11(IC,M3,KU,idsW+1)
303                    continue
                    endif 
!                    !   WARNING WE CAN SPEED UP THE COMPUTATION  !
                    CALL VMX(1,WSQUY1,D11,C11,1)
                    CALL VMX(1,C11,WSQUZ4,A11,1)
!
!
                        IEEE=1
!
                endif 
!**********************************************************************C
!                  2
!**********************************************************************C
                if (iddd == 50) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 304 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,KU,idsW)
304                    continue
                        M4=IUPB(M3,KU,idsW)
                        do 305 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,KU,idsW)
305                    continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 306 IC=1, NCOL2
                        D11(IC)=UB11(IC,M3,KU,idsW+1)
306                    continue
                    endif 
!
                    CALL VMX(1,WSQUZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQDY4,A11,1)
!
                        IEEE=2
!
                endif 
!**********************************************************************C
!                   3                                             C
!**********************************************************************C
                if (iddd == 51) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 307 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
307                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 308 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
308                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 309 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
309                     continue
                    endif 
!     
                    CALL VMX(1,WSQDY1,D11,C11,1)
                    CALL VMX(1,C11,WSQDZ4,A11,1)
!
                        IEEE=3
!
                endif 
!**********************************************************************C
!                 4                                                    C
!**********************************************************************C
                if (iddd == 52) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 310 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
310                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 311 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
311                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 312 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
312                     continue
                    endif 
!     
                    CALL VMX(1,WSQDZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQUY4,A11,1)
!
                        IEEE=4
!
                endif 
!**********************************************************************C
!                 5                                                    C
!**********************************************************************C
                if (iddd == 53) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 313 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
313                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 314 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
314                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 315 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
315                     continue
                    endif 
!
                    CALL VMX(1,WSQUY1,D11,C11,1)
                    CALL VMX(1,C11,WSQDZ4,A11,1)     
!
                        IEEE=5
!
                endif 
!**********************************************************************C
!                 6                                                    C
!**********************************************************************C
                if (iddd == 54) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 316 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
316                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 317 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
317                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 318 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
318                     continue
                    endif 
!
                    CALL VMX(1,WSQUZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQUY4,A11,1)     
!
                        IEEE=6
!
                endif 
!**********************************************************************C
!                 7                                                    C
!**********************************************************************C
                if (iddd == 55) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 319 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
319                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 320 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
320                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 321 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
321                     continue
                    endif 
!
                    CALL VMX(1,WSQDY1,D11,C11,1)
                    CALL VMX(1,C11,WSQUZ4,A11,1)     
!
                        IEEE=7
!
                endif 
!**********************************************************************C
!                 8                                                    C
!**********************************************************************C
                if (iddd == 56) then
!
                    idsW=ids-1
!
                    if (idsW == 0) then
                        idsW=1
                        M3=IUPB(M2,KU,idsW)
                        do 322 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,idsW)
322                     continue
                        M4=IUPB(M3,KU,idsW)
                        do 323 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,idsW)  
323                     continue
                        CALL VMX(1,B11,C11,D11,1)
                    ELSE
                        M3=IUPB(M2,KU,idsW)
                        do 324 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,idsW+1)
324                     continue
                    endif 
!
                    CALL VMX(1,WSQDZ1,D11,C11,1)
                    CALL VMX(1,C11,WSQDY4,A11,1)     
!
                        IEEE=8
!
                endif 
!**********************************************************************
!                 T-T4 OPERATORS
!**********************************************************************
                if (iddd == 57) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,ZIG1W1,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,ZIG2W1,1)
                    CALL VMX(1,ZIG1W1,ZIG2W1,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************
                if (iddd == 58) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,ZIG1W2,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,ZIG2W2,1)
                    CALL VMX(1,ZIG1W2,ZIG2W2,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************
                if (iddd == 59) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,ZIG1W3,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,ZIG2W3,1)
                    CALL VMX(1,ZIG1W3,ZIG2W3,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************
                if (iddd == 60) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,ZIG1W4,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,ZIG2W4,1)
                    CALL VMX(1,ZIG1W4,ZIG2W4,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************
                if (iddd == 61) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,TIG1W4,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,TIG2W4,1)
                    CALL VMX(1,TIG1W4,TIG2W4,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************
                if (iddd == 62) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,TIG1W1,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,TIG2W1,1)
                    CALL VMX(1,TIG1W1,TIG2W1,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************
                if (iddd == 63) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,TIG1W2,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,TIG2W2,1)
                    CALL VMX(1,TIG1W2,TIG2W2,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************
                if (iddd == 64) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,TIG1W3,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,TIG2W3,1)
                    CALL VMX(1,TIG1W3,TIG2W3,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************
!                 T-T5 OPERATORS
!**********************************************************************
                if (iddd == 65) then
!
                    CALL VMX(1,WSQUY1,WSQUY2,DUY1,1)
                    CALL VMX(1,WSQUY3,WSQUY4,DUY2,1)
                    CALL VMX(1,DUY1,DUY2,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************C
                if (iddd == 66) then
!
                    CALL VMX(1,WSQUZ1,WSQUZ2,DUZ1,1)
                    CALL VMX(1,WSQUZ3,WSQUZ4,DUZ2,1)
                    CALL VMX(1,DUZ1,DUZ2,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************C
                if (iddd == 67) then
!
                    CALL VMX(1,WSQDY1,WSQDY2,DDY1,1)
                    CALL VMX(1,WSQDY3,WSQDY4,DDY2,1)
                    CALL VMX(1,DDY1,DDY2,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************C
                if (iddd == 68) then
!
                    CALL VMX(1,WSQDZ1,WSQDZ2,DDZ1,1)
                    CALL VMX(1,WSQDZ3,WSQDZ4,DDZ2,1)
                    CALL VMX(1,DDZ1,DDZ2,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
!                 TT-6 OPERATORS
!**********************************************************************C
                if (iddd == 69) then
!
                    CALL VMX(1,DUY1,DDY2,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************C
                if (iddd == 70) then
!
                    CALL VMX(1,DUZ1,DDZ2,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************C
                if (iddd == 71) then
!
                    CALL VMX(1,DDY1,DUY2,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************C
                if (iddd == 72) then
!
                    CALL VMX(1,DDZ1,DUZ2,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
!                 T-T 7 OPERATORS
!**********************************************************************C
                if (iddd == 73) then
!
                    CALL VMX(1,DUY1,DUZ2,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************C
                if (iddd == 74) then
!
                    CALL VMX(1,DUZ1,DDY2,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************C
                if (iddd == 75) then
!
                    CALL VMX(1,DDY1,DDZ2,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************C
                if (iddd == 76) then
!
                    CALL VMX(1,DDZ1,DUY2,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 77) then
!
                    CALL VMX(1,DDZ1,DDY2,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************C
                if (iddd == 78) then
!
                    CALL VMX(1,DUY1,DDZ2,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************C
                if (iddd == 79) then
!
                    CALL VMX(1,DUZ1,DUY2,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************C
                if (iddd == 80) then
!
                    CALL VMX(1,DDY1,DUZ2,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
!                 T-T 8 OPERATORS
!**********************************************************************C
                if (iddd == 81) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************c
                if (iddd == 82) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************c
                if (iddd == 83) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************c
                if (iddd == 84) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 85) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************c
                if (iddd == 86) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************c
                if (iddd == 87) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************c
                if (iddd == 88) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
!                 T-T 9 OPERATORS
!**********************************************************************C
                if (iddd == 89) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************c
                if (iddd == 90) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************c
                if (iddd == 91) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************c
                if (iddd == 92) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 93) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************c
                if (iddd == 94) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************c
                if (iddd == 95) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************c
                if (iddd == 96) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
!                 T-T 10 OPERATORS
!**********************************************************************C
                if (iddd == 97) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************c
                if (iddd == 98) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************c
                if (iddd == 99) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************c
                if (iddd == 100) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 101) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************c
                if (iddd == 102) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************c
                if (iddd == 103) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************c
                if (iddd == 104) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
!                 T-T 11 OPERATORS
!**********************************************************************C
                if (iddd == 105) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************c
                if (iddd == 106) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************c
                if (iddd == 107) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************c
                if (iddd == 108) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 109) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************c
                if (iddd == 110) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************c
                if (iddd == 111) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************c
                if (iddd == 112) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
!                 T-T 12 OPERATORS
!**********************************************************************C
                if (iddd == 113) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************c
                if (iddd == 114) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************c
                if (iddd == 115) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************c
                if (iddd == 116) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 117) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************c
                if (iddd == 118) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************c
                if (iddd == 119) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************c
                if (iddd == 120) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
                if (iddd == 121) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=9
!
                endif 
!**********************************************************************c
                if (iddd == 122) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=10
!
                endif 
!**********************************************************************c
                if (iddd == 123) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=11
!
                endif 
!**********************************************************************c
                if (iddd == 124) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=12
!
                endif 
!**********************************************************************C
                if (iddd == 125) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=13
!
                endif 
!**********************************************************************c
                if (iddd == 126) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=14
!
                endif 
!**********************************************************************c
                if (iddd == 127) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=15
!
                endif 
!**********************************************************************c
                if (iddd == 128) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=16
!
                endif 
!**********************************************************************C
!                 T-T 13 OPERATORS
!**********************************************************************C
                if (iddd == 129) then
!
                    CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=1
!
                endif 
!**********************************************************************c
                if (iddd == 130) then
!
                    CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=2
!
                endif 
!**********************************************************************c
                if (iddd == 131) then
!
                    CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=3
!
                endif 
!**********************************************************************c
                if (iddd == 132) then
!
                    CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 133) then
!
                    CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                    CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=5
!
                endif 
!**********************************************************************c
                if (iddd == 134) then
!
                    CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                    CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=6
!
                endif 
!**********************************************************************c
                if (iddd == 135) then
!
                    CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                    CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=7
!
                endif 
!**********************************************************************c
                if (iddd == 136) then
!
                    CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                    CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                    CALL VMX(1,B11,C11,A11,1)
!
                    IEEE=8
!
                endif 
!**********************************************************************C
!**********************************************************************C
!                 T-T 14 OPERATORS
!**********************************************************************C
                if (iddd == 137) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQUY1,WSQUZ2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=1
!
                endif 
!**********************************************************************C
                if (iddd == 138) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQUZ1,WSQDY2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=2
!
                endif 
!**********************************************************************C
                if (iddd == 139) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQDY1,WSQDZ2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=3
!
                endif 
!**********************************************************************C
                if (iddd == 140) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQDZ1,WSQUY2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=4
!
                endif 
!**********************************************************************C
                if (iddd == 141) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQUY1,WSQDZ2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=5
!
                endif 
!**********************************************************************C
                if (iddd == 142) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQUZ1,WSQUY2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=6
!
                endif 
!**********************************************************************C
                if (iddd == 143) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQDY1,WSQUZ2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=7
!
                endif 
!**********************************************************************C
                if (iddd == 144) then
!     
                    idsW=ids-1
!     
                    if (idsW == 0) then
                        idsW=1
                    endif 
                    CALL VMX(1,WSQDZ1,WSQDY2,A11,1)
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
                    M3=IUPB(M2,KU,idsW)
                    M2=M3
!
                        IEEE=8
!
                endif 
!**********************************************************************C
!                         NORMAL POLYAKOV LOOP                         C
!**********************************************************************C
                        if  (iddd == 145) then
                        do 167 IC=1, NCOL2
                            A11(IC)=UC11(IC,M2,KU)
                            PL(IC)=UC11(IC,M2,KU)
167                       continue
                        endif 
!**********************************************************************C
!***                      PLAQUETTE OPERATOR                        ***C
!**********************************************************************C
!***                      UP Y OPERATOR
!**********************************************************************C
                        if  (iddd == 146) then
                        do 401 IC=1, NCOL2
                            B11(IC)=UC11(IC,M2,JU)
401                       continue
                        M3=IUP(M2,JU)
                        do 402 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,IU)
402                       continue
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IUP(M2,IU)
                        do 403 IC=1, NCOL2
                            B11(IC)=UC11(IC,M3,JU)
403                       continue
                        CALL HERM(1,B11,DUM11,1)
                        CALL VMX(1,D11,B11,C11,1)
                        do 404 IC=1,NCOL2
                            D11(IC)=UC11(IC,M2,IU)
404                       continue
                        CALL HERM(1,D11,DUM11,1)
                        CALL VMX(1,C11,D11,PLQ1,1)
!     
                        CALL VMX(1,PLQ1,PL,PQ1,1)
!
                        do 505 IC=1,NCOL2
                            A11(IC)=PQ1(IC)
505                       continue
!
                        IEEE=1
!     
                        endif 
!**********************************************************************C
!***                      UP Z OPERATOR
!**********************************************************************C
                        if (iddd == 147) then
!     
                        do 405 IC=1, NCOL2
                            B11(IC)=UC11(IC,M2,IU)
405                       continue
                        M3=IUP(M2,IU)
                        M1=IDN(M3,JU)
                        do 406 IC=1, NCOL2
                            C11(IC)=UC11(IC,M1,JU)
406                       continue
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IDN(M2,JU)
                        do 407 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,IU)
407                       continue
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,D11,C11,B11,1)
                        do 408 IC=1, NCOL2
                            D11(IC)=UC11(IC,M3,JU)
408                       continue
                        CALL VMX(1,B11,D11,PLQ2,1)
!     
                        CALL VMX(1,PLQ2,PL,PQ2,1)
!     
                        do 509 IC=1,NCOL2
                            A11(IC)=PQ2(IC)
509                       continue
!
                        IEEE=2
!     
                        endif 
!**********************************************************************C
!***                      doWN Y OPERATOR
!**********************************************************************C
                        if (iddd == 148) then
!     
                        M3=IDN(M2,JU)
                        do 409 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,JU)
409                       continue
                        CALL HERM(1,C11,DUM11,1)
                        M1=IDN(M3,IU)
                        do 410 IC=1, NCOL2
                            B11(IC)=UC11(IC,M1,IU)
410                       continue
                        CALL HERM(1,B11,DUM11,1)
                        CALL VMX(1,C11,B11,D11,1)
                        do 411 IC=1, NCOL2
                            C11(IC)=UC11(IC,M1,JU)
411                       continue
                        CALL VMX(1,D11,C11,B11,1)
                        M3=IDN(M2,IU)
                        do 412 IC=1, NCOL2
                            D11(IC)=UC11(IC,M3,IU)
412                       continue
                        CALL VMX(1,B11,D11,PLQ3,1)
!     
                        CALL VMX(1,PLQ3,PL,PQ3,1)
!
                        do 513 IC=1,NCOL2
                            A11(IC)=PQ3(IC)
513                       continue
!     
                        IEEE=3
!     
                        endif 
!**********************************************************************C
!***                      doWN Z OPERATOR
!**********************************************************************C
                        if (iddd == 149) then
!     
                        M3=IDN(M2,IU)
                        do 413 IC=1, NCOL2
                            B11(IC)=UC11(IC,M3,IU)
413                       continue
                        CALL HERM(1,B11,DUM11,1)
                        do 414 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,JU)
414                       continue
                        CALL VMX(1,B11,C11,D11,1)
                        M1=IUP(M3,JU)
                        do 415 IC=1, NCOL2
                            B11(IC)=UC11(IC,M1,IU)
415                       continue
                        CALL VMX(1,D11,B11,C11,1)
                        do 416 IC=1, NCOL2
                            B11(IC)=UC11(IC,M2,JU)
416                       continue
                        CALL HERM(1,B11,DUM11,1)
                        CALL VMX(1,C11,B11,PLQ4,1)
!     
                        CALL VMX(1,PLQ4,PL,PQ4,1)
!
                        do 517 IC=1,NCOL2
                            A11(IC)=PQ4(IC)
517                       continue
!
                        IEEE=4
!     
                        endif 
!******************************************************************C
!***                       UP Y OPERATOR +
!******************************************************************C
                        if (iddd == 150)then
!
                        do 801 IC=1, NCOL2
                            PLQ5(IC)=PLQ2(IC)
801                       continue
!
                        CALL HERM(1,PLQ5,DUM11,1)
                        CALL VMX(1,PLQ5,PL,PQ5,1)
!
                        do 518 IC=1,NCOL2
                            A11(IC)=PQ5(IC)
518                       continue
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
!***                       UP Z OPERATOR +
!******************************************************************C
                        if (iddd == 151)then
!     
                        do 802 IC=1, NCOL2
                            PLQ6(IC)=PLQ1(IC)
802                       continue
!
                        CALL HERM(1,PLQ6,DUM11,1)
                        CALL VMX(1,PLQ6,PL,PQ6,1)
!
                        do 519 IC=1,NCOL2
                            A11(IC)=PQ6(IC)
519                       continue
!
                        IEEE=6
!     
                        endif 
!******************************************************************C
!***                       doWN Y OPERATOR +
!******************************************************************C
                        if (iddd == 152)then
!     
                        do 803 IC=1, NCOL2
                            PLQ7(IC)=PLQ4(IC)
803                       continue
!
                        CALL HERM(1,PLQ7,DUM11,1)
                        CALL VMX(1,PLQ7,PL,PQ7,1)
!     
                        do 520 IC=1,NCOL2
                            A11(IC)=PQ7(IC)
520                       continue
!
                        IEEE=7
!     
                        endif 
!******************************************************************C
!***                       doWN Z OPERATOR +
!******************************************************************C
                        if (iddd == 153)then
!     
                        do 804 IC=1, NCOL2
                            PLQ8(IC)=PLQ3(IC)
804                       continue
!
                        CALL HERM(1,PLQ8,DUM11,1)
                        CALL VMX(1,PLQ8,PL,PQ8,1)
!
                        do 521 IC=1,NCOL2
                            A11(IC)=PQ8(IC)
521                       continue
!
                        IEEE=8
!
                        endif 
!**********************************************************************C
!***                      PLAQUETTE OPERATOR 2                      ***C
!**********************************************************************C
!***                      UP Y OPERATOR
!**********************************************************************C
                        if  (iddd == 154) then
                        M5=IUP(M2,KU)
                        do 601 IC=1, NCOL2
                            B11(IC)=UC11(IC,M5,JU)
601                       continue
                        M3=IUP(M5,JU)
                        do 602 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,IU)
602                       continue
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IUP(M5,IU)
                        do 603 IC=1, NCOL2
                            B11(IC)=UC11(IC,M3,JU)
603                       continue
                        CALL HERM(1,B11,DUM11,1)
                        CALL VMX(1,D11,B11,C11,1)
                        do 604 IC=1,NCOL2
                            D11(IC)=UC11(IC,M5,IU)
604                       continue
                        CALL HERM(1,D11,DUM11,1)
                        CALL VMX(1,C11,D11,DPLQ1,1)
!     
                        CALL VMX(1,PQ1,DPLQ1,A11,1)
!
                        IEEE=1
!     
                        endif 
!**********************************************************************C
!***                      UP Z OPERATOR
!**********************************************************************C
                        if (iddd == 155) then
!     
                        M5=IUP(M2,KU)
                        do 605 IC=1, NCOL2
                            B11(IC)=UC11(IC,M5,IU)
605                       continue
                        M3=IUP(M5,IU)
                        M1=IDN(M3,JU)
                        do 606 IC=1, NCOL2
                            C11(IC)=UC11(IC,M1,JU)
606                       continue
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IDN(M5,JU)
                        do 607 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,IU)
607                       continue
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,D11,C11,B11,1)
                        do 608 IC=1, NCOL2
                            D11(IC)=UC11(IC,M3,JU)
608                       continue
                        CALL VMX(1,B11,D11,DPLQ2,1)
!     
                        CALL VMX(1,PQ2,DPLQ2,A11,1)
!     
                        IEEE=2
!     
                        endif 
!**********************************************************************C
!***                      doWN Y OPERATOR
!**********************************************************************C
                        if (iddd == 156) then
!     
                        M5=IUP(M2,KU)
                        M3=IDN(M5,JU)
                        do 609 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,JU)
609                       continue
                        CALL HERM(1,C11,DUM11,1)
                        M1=IDN(M3,IU)
                        do 610 IC=1, NCOL2
                            B11(IC)=UC11(IC,M1,IU)
610                       continue
                        CALL HERM(1,B11,DUM11,1)
                        CALL VMX(1,C11,B11,D11,1)
                        do 611 IC=1, NCOL2
                            C11(IC)=UC11(IC,M1,JU)
611                       continue
                        CALL VMX(1,D11,C11,B11,1)
                        M3=IDN(M5,IU)
                        do 612 IC=1, NCOL2
                            D11(IC)=UC11(IC,M3,IU)
612                       continue
                        CALL VMX(1,B11,D11,DPLQ3,1)
!     
                        CALL VMX(1,PQ3,DPLQ3,A11,1)
!     
                        IEEE=3
!     
                        endif 
!**********************************************************************C
!***                      doWN Z OPERATOR
!**********************************************************************C
                        if (iddd == 157) then
!     
                        M5=IUP(M2,KU)
                        M3=IDN(M5,IU)
                        do 613 IC=1, NCOL2
                            B11(IC)=UC11(IC,M3,IU)
613                       continue
                        CALL HERM(1,B11,DUM11,1)
                        do 614 IC=1, NCOL2
                            C11(IC)=UC11(IC,M3,JU)
614                       continue
                        CALL VMX(1,B11,C11,D11,1)
                        M1=IUP(M3,JU)
                        do 615 IC=1, NCOL2
                            B11(IC)=UC11(IC,M1,IU)
615                       continue
                        CALL VMX(1,D11,B11,C11,1)
                        do 616 IC=1, NCOL2
                            B11(IC)=UC11(IC,M5,JU)
616                       continue
                        CALL HERM(1,B11,DUM11,1)
                        CALL VMX(1,C11,B11,DPLQ4,1)
!     
                        CALL VMX(1,PQ4,DPLQ4,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
!***                       UP Y OPERATOR +
!******************************************************************C
                        if (iddd == 158)then
!
                        do 701 IC=1, NCOL2
                            DPLQ5(IC)=DPLQ2(IC)
701                       continue
!     
                        CALL HERM(1,DPLQ5,DUM11,1)
                        CALL VMX(1,PQ5,DPLQ5,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
!***                       UP Z OPERATOR +
!******************************************************************C
                        if (iddd == 159)then
!     
                        do 702 IC=1, NCOL2
                            DPLQ6(IC)=DPLQ1(IC)
702                       continue
!
                        CALL HERM(1,DPLQ6,DUM11,1)
                        CALL VMX(1,PQ6,DPLQ6,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
!***  doWN Y OPERATOR +
!******************************************************************C
                        if (iddd == 160)then
!     
                        do 703 IC=1, NCOL2
                            DPLQ7(IC)=DPLQ4(IC)
703                       continue
!
                        CALL HERM(1,DPLQ7,DUM11,1)
                        CALL VMX(1,PQ7,DPLQ7,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
!***                       doWN Z OPERATOR +
!******************************************************************C
                        if (iddd == 161)then
!     
                        do 704 IC=1, NCOL2
                            DPLQ8(IC)=DPLQ3(IC)
704                       continue
!
                        CALL HERM(1,DPLQ8,DUM11,1)
                        CALL VMX(1,PQ8,DPLQ8,A11,1)
!
                        IEEE=8
!
                        endif 
!******************************************************************C
!***             PLAQUETTE OPERATPORS 3
!******************************************************************C
                        if (iddd == 162)then
!
!                           do 705 IC=1, NCOL2
!                              B11(IC)=DPLQ1(IC)
! 705                       continue
!
!                           CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ1,B11,A11,1)
                            CALL VMX(1,PQ1,DPLQ6,A11,1) ! New altered !                           
!     
                        IEEE=1
!
                        endif 
!******************************************************************C
                        if (iddd == 163)then
!
!                           do 706 IC=1, NCOL2
!                              B11(IC)=DPLQ2(IC)
! 706                       continue
!C
!                          CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ2,B11,A11,1)                           
                        CALL VMX(1,PQ2,DPLQ5,A11,1) ! New altered !                           
!     
                        IEEE=2
!
                        endif 
!******************************************************************C
                        if (iddd == 164)then
!
!                           do 707 IC=1, NCOL2
!                              B11(IC)=DPLQ3(IC)
! 707                       continue
!C
!                           CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ3,B11,A11,1)
                        
                            CALL VMX(1,PQ3,DPLQ8,A11,1) ! New altered !
!     
                        IEEE=3
!
                        endif 
!******************************************************************C
                        if (iddd == 165)then
!
!                           do 708 IC=1, NCOL2
!                              B11(IC)=DPLQ4(IC)
! 708                       continue
!C
!                          CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ4,B11,A11,1)                           
                        CALL VMX(1,PQ4,DPLQ7,A11,1) ! New altered ! 
!     
                        IEEE=4
!
                        endif 
!******************************************************************C
                        if (iddd == 166)then
!
                        CALL VMX(1,PQ5,DPLQ2,A11,1)
!     
                        IEEE=5
!
                        endif 
!******************************************************************C
                        if (iddd == 167)then
!
                        CALL VMX(1,PQ6,DPLQ1,A11,1)
!     
                        IEEE=6
!
                        endif 
!******************************************************************C
                        if (iddd == 168)then
!
                        CALL VMX(1,PQ7,DPLQ4,A11,1)
!     
                        IEEE=7
!
                        endif 
!******************************************************************C
                        if (iddd == 169)then
!
                        CALL VMX(1,PQ8,DPLQ3,A11,1)
!     
                        IEEE=8
!
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 4 
!******************************************************************C
                        if (iddd == 170)then
!
!                           do 709 IC=1, NCOL2
!                              C11(IC)=DPLQ1(IC)
! 709                       continue
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ1,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ6,A11,1)
!
                        IEEE=1
!
                        endif 
!******************************************************************C
                        if (iddd == 171)then
!
!                           do 710 IC=1, NCOL2
!                              C11(IC)=DPLQ2(IC)
! 710                       continue
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ2,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ5,A11,1)
!
                        IEEE=2
!
                        endif 
!******************************************************************C
                        if (iddd == 172)then
!
!                           do 711 IC=1, NCOL2
!                              C11(IC)=DPLQ3(IC)
! 711                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ3,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ8,A11,1)
!
                        IEEE=3
!
                        endif 
!******************************************************************C
                        if (iddd == 173)then
!
!                           do 712 IC=1, NCOL2
!                              C11(IC)=DPLQ4(IC)
! 712                       continue
!
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ4,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ7,A11,1)
!
                        IEEE=4
!
                        endif 
!******************************************************************C
                        if (iddd == 174)then
!
!                           do 713 IC=1, NCOL2
!                              C11(IC)=PLQ2(IC)
! 713                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ5,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ2,A11,1)
!
                        IEEE=5
!
                        endif 
!******************************************************************C
                        if (iddd == 175)then
!
!                           do 714 IC=1, NCOL2
!                              C11(IC)=PLQ1(IC)
! 714                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ6,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ1,A11,1)
!
                        IEEE=6
!
                        endif 
!******************************************************************C
                        if (iddd == 176)then
!
!                           do 715 IC=1, NCOL2
!                              C11(IC)=PLQ4(IC)
! 715                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ7,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ4,A11,1)
!
                        IEEE=7
!
                        endif 
!******************************************************************C
                        if (iddd == 177)then
!
!                           do 716 IC=1, NCOL2
!                              C11(IC)=PLQ3(IC)
! 716                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PLQ8,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ3,A11,1)
!
                        IEEE=8
!
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 5
!******************************************************************C
                        if (iddd == 178)then
!
!                           do 717 IC=1, NCOL2
!                              C11(IC)=DPLQ3(IC)
! 717                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ1,DPLQ8,A11,1)
!
                        IEEE=1
!
                        endif 
!******************************************************************C
                        if (iddd == 179)then
!
!                           do 718 IC=1, NCOL2
!                              C11(IC)=DPLQ4(IC)
! 718                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ2,DPLQ7,A11,1)
!
                        IEEE=2
!
                        endif 
!******************************************************************C
                        if (iddd == 180)then
!
!                           do 719 IC=1, NCOL2
!                              C11(IC)=DPLQ1(IC)
! 719                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ3,DPLQ6,A11,1)
!
                        IEEE=3
!
                        endif 
!******************************************************************C
                        if (iddd == 181)then
!
!                           do 720 IC=1, NCOL2
!                              C11(IC)=DPLQ2(IC)
! 720                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ4,DPLQ5,A11,1)
!
                        IEEE=4
!
                        endif 
!******************************************************************C
                        if (iddd == 182)then
!
                        CALL VMX(1,PQ5,DPLQ4,A11,1)
!
                        IEEE=5
!
                        endif 
!******************************************************************C
                        if (iddd == 183)then
!
                        CALL VMX(1,PQ6,DPLQ3,A11,1)
!
                        IEEE=6
!
                        endif 
!******************************************************************C
                        if (iddd == 184)then
!
                        CALL VMX(1,PQ7,DPLQ2,A11,1)
!
                        IEEE=7
!
                        endif 
!******************************************************************C
                        if (iddd == 185)then
!
                        CALL VMX(1,PQ8,DPLQ1,A11,1)
!
                        IEEE=8
!
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 6
!******************************************************************C
                        if (iddd == 186)then
!
                        CALL VMX(1,PQ1,DPLQ3,A11,1)
!     
                        IEEE=1
!
                        endif 
!******************************************************************C
                        if (iddd == 187)then
!
                        CALL VMX(1,PQ2,DPLQ4,A11,1)
!     
                        IEEE=2
!
                        endif 
!******************************************************************C
                        if (iddd == 188)then
!
                        CALL VMX(1,PQ3,DPLQ1,A11,1)
!     
                        IEEE=3
!
                        endif 
!******************************************************************C
                        if (iddd == 189)then
!
                        CALL VMX(1,PQ4,DPLQ2,A11,1)
!     
                        IEEE=4
!
                        endif 
!******************************************************************C
                        if (iddd == 190)then
!
!                           do 721 IC=1, NCOL2
!                              C11(IC)=DPLQ4(IC)
! 721                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ5,DPLQ7,A11,1)
!     
                        IEEE=5
!
                        endif 
!******************************************************************C
                        if (iddd == 191)then
!
!                           do 722 IC=1, NCOL2
!                              C11(IC)=DPLQ3(IC)
! 722                       continue
!
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ6,DPLQ8,A11,1)
!     
                        IEEE=6
!
                        endif 
!******************************************************************C
                        if (iddd == 192)then
!
!                           do 723 IC=1, NCOL2
!                              C11(IC)=DPLQ2(IC)
! 723                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ7,DPLQ5,A11,1)
!     
                        IEEE=7
!
                        endif 
!******************************************************************C
                        if (iddd == 193)then
!
!                           do 724 IC=1, NCOL2
!                              C11(IC)=DPLQ1(IC)
! 724                       continue
!C
!                           CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,PQ8,DPLQ6,A11,1)
!     
                        IEEE=8
!
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 7
!******************************************************************C
                        if (iddd == 194)then
!
                        CALL VMX(1,SQUY1,DPLQ1,B11,1)
                        CALL VMX(1,B11,SQDZ2,A11,1)
!     
                        IEEE=1
!
                        endif 
!******************************************************************C                        
                        if (iddd == 195)then
!
                        CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                        CALL VMX(1,B11,SQUY2,A11,1)
!     
                        IEEE=2
!
                        endif 
!******************************************************************C                        
                        if (iddd == 196)then
!
                        CALL VMX(1,SQDY1,DPLQ3,B11,1)
                        CALL VMX(1,B11,SQUZ2,A11,1)
!     
                        IEEE=3
!
                        endif 
!******************************************************************C                        
                        if (iddd == 197)then
!
                        CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                        CALL VMX(1,B11,SQDY2,A11,1)
!     
                        IEEE=4
!
                        endif 
!******************************************************************C                        
                        if (iddd == 198)then
!
                        CALL VMX(1,SQDZ1,DPLQ6,B11,1)
                        CALL VMX(1,B11,SQUY2,A11,1)
!     
                        IEEE=5
!
                        endif 
!******************************************************************C                        
                        if (iddd == 199)then
!
                        CALL VMX(1,SQUY1,DPLQ5,B11,1)
                        CALL VMX(1,B11,SQUZ2,A11,1)
!     
                        IEEE=6
!
                        endif 
!******************************************************************C                        
                        if (iddd == 200)then
!
                        CALL VMX(1,SQUZ1,DPLQ8,B11,1)
                        CALL VMX(1,B11,SQDY2,A11,1)
!     
                        IEEE=7
!
                        endif 
!******************************************************************C                        
                        if (iddd == 201)then
!
                        CALL VMX(1,SQDY1,DPLQ7,B11,1)
                        CALL VMX(1,B11,SQDZ2,A11,1)
!     
                        IEEE=8
!
                        endif 
!******************************************************************C                        
                        if (iddd == 202)then
!
                        CALL VMX(1,SQDY1,DPLQ5,B11,1)
                        CALL VMX(1,B11,SQDZ2,A11,1)
!     
                        IEEE=9
!
                        endif 
!******************************************************************C                        
                        if (iddd == 203)then
!
                        CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                        CALL VMX(1,B11,SQDY2,A11,1)
!     
                        IEEE=10
!
                        endif 
!******************************************************************C                        
                        if (iddd == 204)then
!
                        CALL VMX(1,SQUY1,DPLQ7,B11,1)
                        CALL VMX(1,B11,SQUZ2,A11,1)
!     
                        IEEE=11
!
                        endif 
!******************************************************************C                        
                        if (iddd == 205)then
!
                        CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                        CALL VMX(1,B11,SQUY2,A11,1)
!     
                        IEEE=12
!
                        endif 
!******************************************************************C                        
                        if (iddd == 206)then
!
                        CALL VMX(1,SQDZ1,DPLQ2,B11,1)
                        CALL VMX(1,B11,SQDY2,A11,1)
!     
                        IEEE=13
!
                        endif 
!******************************************************************C                        
                        if (iddd == 207)then
!
                        CALL VMX(1,SQDY1,DPLQ1,B11,1)
                        CALL VMX(1,B11,SQUZ2,A11,1)
!     
                        IEEE=14
!
                        endif 
!******************************************************************C                        
                        if (iddd == 208)then
!
                        CALL VMX(1,SQUZ1,DPLQ4,B11,1)
                        CALL VMX(1,B11,SQUY2,A11,1)
!     
                        IEEE=15
!
                        endif 
!******************************************************************C                        
                        if (iddd == 209)then
!
                        CALL VMX(1,SQUY1,DPLQ3,B11,1)
                        CALL VMX(1,B11,SQDZ2,A11,1)
!     
                        IEEE=16
!
                        endif                             
!******************************************************************C
!     PLAQUETTE OPERATORS 8
!******************************************************************C
                        if (iddd == 210)then
!
                        CALL VMX(1,SQUY1,DPLQ1,B11,1)
                        CALL VMX(1,B11,DPLQ2,A11,1)
!     
                        IEEE=1
!     
                        endif 
!******************************************************************c
                        if (iddd == 211)then
!
                        CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                        CALL VMX(1,B11,DPLQ3,A11,1)
!     
                        IEEE=2
!     
                        endif  
!******************************************************************c
                        if (iddd == 212)then
!
                        CALL VMX(1,SQDY1,DPLQ3,B11,1)
                        CALL VMX(1,B11,DPLQ4,A11,1)
!     
                        IEEE=3
!     
                        endif  
!******************************************************************c
                        if (iddd == 213)then
!
                        CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                        CALL VMX(1,B11,DPLQ1,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************c
                        if (iddd == 214)then
!
                        CALL VMX(1,PLQ5,PLQ6,B11,1)
                        CALL VMX(1,B11,SQUY1,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************c
                        if (iddd == 215)then
!
                        CALL VMX(1,PLQ8,PLQ5,B11,1)
                        CALL VMX(1,B11,SQUZ1,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************c
                        if (iddd == 216)then
!
                        CALL VMX(1,PLQ7,PLQ8,B11,1)
                        CALL VMX(1,B11,SQDY1,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************c
                        if (iddd == 217)then
!
                        CALL VMX(1,PLQ6,PLQ7,B11,1)
                        CALL VMX(1,B11,SQDZ1,A11,1)
!     
                        IEEE=8
!     
                        endif 
!******************************************************************c
                        if (iddd == 218)then
!
                        CALL VMX(1,SQDY1,DPLQ5,B11,1)
                        CALL VMX(1,B11,DPLQ6,A11,1)
!     
                        IEEE=9
!     
                        endif 
!******************************************************************c
                        if (iddd == 219)then
!
                        CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                        CALL VMX(1,B11,DPLQ7,A11,1)
!     
                        IEEE=10
!     
                        endif 
!******************************************************************c
                        if (iddd == 220)then
!
                        CALL VMX(1,SQUY1,DPLQ7,B11,1)
                        CALL VMX(1,B11,DPLQ8,A11,1)
!     
                        IEEE=11
!     
                        endif 
!******************************************************************c
                        if (iddd == 221)then
!
                        CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                        CALL VMX(1,B11,DPLQ5,A11,1)
!     
                        IEEE=12
!     
                        endif 
!******************************************************************c
                        if (iddd == 222)then
!
                        CALL VMX(1,PLQ1,PLQ2,B11,1)
                        CALL VMX(1,B11,SQDY1,A11,1)
!     
                        IEEE=13
!     
                        endif 
!******************************************************************c
                        if (iddd == 223)then
!
                        CALL VMX(1,PLQ4,PLQ1,B11,1)
                        CALL VMX(1,B11,SQUZ1,A11,1)
!     
                        IEEE=14
!     
                        endif 
!******************************************************************c
                        if (iddd == 224)then
!
                        CALL VMX(1,PLQ3,PLQ4,B11,1)
                        CALL VMX(1,B11,SQUY1,A11,1)
!     
                        IEEE=15
!     
                        endif 
!******************************************************************c
                        if (iddd == 225)then
!
                        CALL VMX(1,PLQ2,PLQ3,B11,1)
                        CALL VMX(1,B11,SQDZ1,A11,1)
!     
                        IEEE=16
!     
                        endif                         
!******************************************************************C
!     PLAQUETTE OPERATORS 9
!******************************************************************C
                        if (iddd == 226)then
!     
                        CALL VMX(1,SQUY1,DPLQ1,B11,1)
                        CALL VMX(1,B11,DPLQ2,C11,1)
                        CALL VMX(1,C11,DPLQ3,A11,1)                           
!     
                        IEEE=1
!     
                        endif 
!******************************************************************C
                        if (iddd == 227)then
!     
                        CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                        CALL VMX(1,B11,DPLQ3,C11,1)
                        CALL VMX(1,C11,DPLQ4,A11,1)                           
!     
                        IEEE=2
!     
                        endif 
!******************************************************************C
                        if (iddd == 228)then
!     
                        CALL VMX(1,SQDY1,DPLQ3,B11,1)
                        CALL VMX(1,B11,DPLQ4,C11,1)
                        CALL VMX(1,C11,DPLQ1,A11,1)                           
!     
                        IEEE=3
!     
                        endif 
!******************************************************************C
                        if (iddd == 229)then
!     
                        CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                        CALL VMX(1,B11,DPLQ1,C11,1)
                        CALL VMX(1,C11,DPLQ2,A11,1)                           
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
                        if (iddd == 230)then
!     
                        CALL VMX(1,PLQ8,PLQ5,B11,1)
                        CALL VMX(1,B11,PLQ6,C11,1)
                        CALL VMX(1,C11,SQUY1,A11,1)                           
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
                        if (iddd == 231)then
!     
                        CALL VMX(1,PLQ7,PLQ8,B11,1)
                        CALL VMX(1,B11,PLQ5,C11,1)
                        CALL VMX(1,C11,SQUZ1,A11,1)                           
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
                        if (iddd == 232)then
!     
                        CALL VMX(1,PLQ6,PLQ7,B11,1)
                        CALL VMX(1,B11,PLQ8,C11,1)
                        CALL VMX(1,C11,SQDY1,A11,1)                           
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
                        if (iddd == 233)then
!     
                        CALL VMX(1,PLQ5,PLQ6,B11,1)
                        CALL VMX(1,B11,PLQ7,C11,1)
                        CALL VMX(1,C11,SQDZ1,A11,1)                           
!     
                        IEEE=8
!     
                        endif 
!******************************************************************C
                        if (iddd == 234)then
!
                        CALL VMX(1,SQDY1,DPLQ5,B11,1)
                        CALL VMX(1,B11,DPLQ6,C11,1)
                        CALL VMX(1,C11,DPLQ7,A11,1)    
!     
                        IEEE=9
!     
                        endif 
!******************************************************************C
                        if (iddd == 235)then
!
                        CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                        CALL VMX(1,B11,DPLQ7,C11,1)
                        CALL VMX(1,C11,DPLQ8,A11,1)    
!     
                        IEEE=10
!     
                        endif 
!******************************************************************C
                        if (iddd == 236)then
!
                        CALL VMX(1,SQUY1,DPLQ7,B11,1)
                        CALL VMX(1,B11,DPLQ8,C11,1)
                        CALL VMX(1,C11,DPLQ5,A11,1)    
!     
                        IEEE=11
!     
                        endif 
!******************************************************************C
                        if (iddd == 237)then
!
                        CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                        CALL VMX(1,B11,DPLQ5,C11,1)
                        CALL VMX(1,C11,DPLQ6,A11,1)    
!     
                        IEEE=12
!     
                        endif 
!******************************************************************C
                        if (iddd == 238)then
!
                        CALL VMX(1,PLQ4,PLQ1,B11,1)
                        CALL VMX(1,B11,PLQ2,C11,1)
                        CALL VMX(1,C11,SQDY1,A11,1)    
!     
                        IEEE=13
!     
                        endif 
!******************************************************************C
                        if (iddd == 239)then
!
                        CALL VMX(1,PLQ3,PLQ4,B11,1)
                        CALL VMX(1,B11,PLQ1,C11,1)
                        CALL VMX(1,C11,SQUZ1,A11,1)    
!     
                        IEEE=14
!     
                        endif 
!******************************************************************C
                        if (iddd == 240)then
!
                        CALL VMX(1,PLQ2,PLQ3,B11,1)
                        CALL VMX(1,B11,PLQ4,C11,1)
                        CALL VMX(1,C11,SQUY1,A11,1)    
!     
                        IEEE=15
!     
                        endif 
!******************************************************************C
                        if (iddd == 241)then
!
                        CALL VMX(1,PLQ1,PLQ2,B11,1)
                        CALL VMX(1,B11,PLQ3,C11,1)
                        CALL VMX(1,C11,SQDZ1,A11,1)    
!     
                        IEEE=16
!     
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 10
!******************************************************************C
                        if (iddd == 242)then
!     
                        CALL VMX(1,PLQ2,PLQ7,B11,1)
                        CALL VMX(1,B11,SQDZ1,A11,1)
!     
                        IEEE=1
!     
                        endif 
!******************************************************************C
                        if (iddd == 243)then
!     
                        CALL VMX(1,PLQ3,PLQ6,B11,1)
                        CALL VMX(1,B11,SQUY1,A11,1)
!     
                        IEEE=2
!     
                        endif 
!******************************************************************C
                        if (iddd == 244)then
!     
                        CALL VMX(1,PLQ4,PLQ5,B11,1)
                        CALL VMX(1,B11,SQUZ1,A11,1)
!     
                        IEEE=3
!     
                        endif 
!******************************************************************C
                        if (iddd == 245)then
!     
                        CALL VMX(1,PLQ1,PLQ8,B11,1)
                        CALL VMX(1,B11,SQDY1,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
                        if (iddd == 246)then
!     
                        CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                        CALL VMX(1,B11,DPLQ5,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
                        if (iddd == 247)then
!     
                        CALL VMX(1,SQUY1,DPLQ1,B11,1)
                        CALL VMX(1,B11,DPLQ8,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
                        if (iddd == 248)then
!     
                        CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                        CALL VMX(1,B11,DPLQ7,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
                        if (iddd == 249)then
!     
                        CALL VMX(1,SQDY1,DPLQ3,B11,1)
                        CALL VMX(1,B11,DPLQ6,A11,1)
!     
                        IEEE=8
!     
                        endif 
!******************************************************************C
                        if (iddd == 250)then
!     
                        CALL VMX(1,PLQ6,PLQ3,B11,1)
                        CALL VMX(1,B11,SQDZ1,A11,1)
!     
                        IEEE=9
!     
                        endif 
!******************************************************************C
                        if (iddd == 251)then
!     
                        CALL VMX(1,PLQ7,PLQ2,B11,1)
                        CALL VMX(1,B11,SQDY1,A11,1)
!     
                        IEEE=10
!     
                        endif 
!******************************************************************C
                        if (iddd == 252)then
!     
                        CALL VMX(1,PLQ8,PLQ1,B11,1)
                        CALL VMX(1,B11,SQUZ1,A11,1)
!     
                        IEEE=11
!     
                        endif 
!******************************************************************C
                        if (iddd == 253)then
!     
                        CALL VMX(1,PLQ5,PLQ4,B11,1)
                        CALL VMX(1,B11,SQUY1,A11,1)
!     
                        IEEE=12
!     
                        endif 
!******************************************************************C
                        if (iddd == 254)then
!     
                        CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                        CALL VMX(1,B11,DPLQ1,A11,1)
!     
                        IEEE=13
!     
                        endif 
!******************************************************************C
                        if (iddd == 255)then
!     
                        CALL VMX(1,SQDY1,DPLQ5,B11,1)
                        CALL VMX(1,B11,DPLQ4,A11,1)
!     
                        IEEE=14
!     
                        endif 
!******************************************************************C
                        if (iddd == 256)then
!     
                        CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                        CALL VMX(1,B11,DPLQ3,A11,1)
!     
                        IEEE=15
!     
                        endif 
!******************************************************************C
                        if (iddd == 257)then
!     
                        CALL VMX(1,SQUY1,DPLQ7,B11,1)
                        CALL VMX(1,B11,DPLQ2,A11,1)
!     
                        IEEE=16
!     
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 11
!******************************************************************C
                        if (iddd == 258)then
!     
                        CALL VMX(1,PLQ2,PLQ4,B11,1)
                        CALL VMX(1,B11,SQUY1,A11,1)
!     
                        IEEE=1
!     
                        endif 
!******************************************************************C
                        if (iddd == 259)then
!     
                        CALL VMX(1,PLQ3,PLQ1,B11,1)
                        CALL VMX(1,B11,SQUZ1,A11,1)
!     
                        IEEE=2
!     
                        endif 
!******************************************************************C
                        if (iddd == 260)then
!     
                        CALL VMX(1,PLQ4,PLQ2,B11,1)
                        CALL VMX(1,B11,SQDY1,A11,1)
!     
                        IEEE=3
!     
                        endif 
!******************************************************************C
                        if (iddd == 261)then
!     
                        CALL VMX(1,PLQ1,PLQ3,B11,1)
                        CALL VMX(1,B11,SQDZ1,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
                        if (iddd == 262)then
!
                        CALL VMX(1,SQUY1,DPLQ7,B11,1)
                        CALL VMX(1,B11,DPLQ5,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
                        if (iddd == 263)then
!
                        CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                        CALL VMX(1,B11,DPLQ8,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
                        if (iddd == 264)then
!
                        CALL VMX(1,SQDY1,DPLQ5,B11,1)
                        CALL VMX(1,B11,DPLQ7,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
                        if (iddd == 265)then
!
                        CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                        CALL VMX(1,B11,DPLQ6,A11,1)
!     
                        IEEE=8
!     
                        endif 
!******************************************************************C
                        if (iddd == 266)then
!     
                        CALL VMX(1,PLQ6,PLQ8,B11,1)
                        CALL VMX(1,B11,SQDY1,A11,1)
!     
                        IEEE=9
!     
                        endif 
!******************************************************************C
                        if (iddd == 267)then
!     
                        CALL VMX(1,PLQ7,PLQ5,B11,1)
                        CALL VMX(1,B11,SQUZ1,A11,1)
!     
                        IEEE=10
!     
                        endif 
!******************************************************************C
                        if (iddd == 268)then
!     
                        CALL VMX(1,PLQ8,PLQ6,B11,1)
                        CALL VMX(1,B11,SQUY1,A11,1)
!     
                        IEEE=11
!     
                        endif 
!******************************************************************C
                        if (iddd == 269)then
!     
                        CALL VMX(1,PLQ5,PLQ7,B11,1)
                        CALL VMX(1,B11,SQDZ1,A11,1)
!     
                        IEEE=12
!     
                        endif 
!******************************************************************C
                        if (iddd == 270)then
!
                        CALL VMX(1,SQDY1,DPLQ3,B11,1)
                        CALL VMX(1,B11,DPLQ1,A11,1)
!     
                        IEEE=13
!     
                        endif 
!******************************************************************C
                        if (iddd == 271)then
!
                        CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                        CALL VMX(1,B11,DPLQ4,A11,1)
!     
                        IEEE=14
!     
                        endif                         
!******************************************************************C
                        if (iddd == 272)then
!
                        CALL VMX(1,SQUY1,DPLQ1,B11,1)
                        CALL VMX(1,B11,DPLQ3,A11,1)
!     
                        IEEE=15
!     
                        endif 
!******************************************************************C
                        if (iddd == 273)then
!
                        CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                        CALL VMX(1,B11,DPLQ2,A11,1)
!     
                        IEEE=16
!     
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 12
!******************************************************************C
                        if (iddd == 274)then
!     
                        CALL VMX(1,PLQ2,PLQ3,B11,1)
                        CALL VMX(1,B11,SQUY1,C11,1)
                        CALL VMX(1,C11,DPLQ7,A11,1)
!     
                        IEEE=1
!     
                        endif 
!******************************************************************C
                        if (iddd == 275)then
!     
                        CALL VMX(1,PLQ3,PLQ4,B11,1)
                        CALL VMX(1,B11,SQUZ1,C11,1)
                        CALL VMX(1,C11,DPLQ6,A11,1)
!     
                        IEEE=2
!     
                        endif 
!******************************************************************C
                        if (iddd == 276)then
!     
                        CALL VMX(1,PLQ4,PLQ1,B11,1)
                        CALL VMX(1,B11,SQDY1,C11,1)
                        CALL VMX(1,C11,DPLQ5,A11,1)
!     
                        IEEE=3
!     
                        endif 
!******************************************************************C
                        if (iddd == 277)then
!     
                        CALL VMX(1,PLQ1,PLQ2,B11,1)
                        CALL VMX(1,B11,SQDZ1,C11,1)
                        CALL VMX(1,C11,DPLQ8,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
                        if (iddd == 278)then
!     
                        CALL VMX(1,PLQ4,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ8,C11,1)
                        CALL VMX(1,C11,DPLQ5,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
                        if (iddd == 279)then
!     
                        CALL VMX(1,PLQ1,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ7,C11,1)
                        CALL VMX(1,C11,DPLQ8,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
                        if (iddd == 280)then
!     
                        CALL VMX(1,PLQ2,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ6,C11,1)
                        CALL VMX(1,C11,DPLQ7,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
                        if (iddd == 281)then
!     
                        CALL VMX(1,PLQ3,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ5,C11,1)
                        CALL VMX(1,C11,DPLQ6,A11,1)
!     
                        IEEE=8
!     
                        endif 
!******************************************************************C
                        if (iddd == 282)then
!     
                        CALL VMX(1,PLQ6,PLQ7,B11,1)
                        CALL VMX(1,B11,SQDY1,C11,1)
                        CALL VMX(1,C11,DPLQ3,A11,1)
!     
                        IEEE=9
!     
                        endif 
!******************************************************************C
                        if (iddd == 283)then
!     
                        CALL VMX(1,PLQ7,PLQ8,B11,1)
                        CALL VMX(1,B11,SQUZ1,C11,1)
                        CALL VMX(1,C11,DPLQ2,A11,1)
!     
                        IEEE=10
!     
                        endif 
!******************************************************************C
                        if (iddd == 284)then
!     
                        CALL VMX(1,PLQ8,PLQ5,B11,1)
                        CALL VMX(1,B11,SQUY1,C11,1)
                        CALL VMX(1,C11,DPLQ1,A11,1)
!     
                        IEEE=11
!     
                        endif 
!******************************************************************C
                        if (iddd == 285)then
!     
                        CALL VMX(1,PLQ5,PLQ6,B11,1)
                        CALL VMX(1,B11,SQDZ1,C11,1)
                        CALL VMX(1,C11,DPLQ4,A11,1)
!     
                        IEEE=12
!     
                        endif 
!******************************************************************C
                        if (iddd == 286)then
!     
                        CALL VMX(1,PLQ8,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ4,C11,1)
                        CALL VMX(1,C11,DPLQ1,A11,1)
!     
                        IEEE=13
!     
                        endif 
!******************************************************************C
                        if (iddd == 287)then
!     
                        CALL VMX(1,PLQ5,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ3,C11,1)
                        CALL VMX(1,C11,DPLQ4,A11,1)
!     
                        IEEE=14
!     
                        endif 
!******************************************************************C
                        if (iddd == 288)then
!     
                        CALL VMX(1,PLQ6,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ2,C11,1)
                        CALL VMX(1,C11,DPLQ3,A11,1)
!     
                        IEEE=15
!     
                        endif 
!******************************************************************C
                        if (iddd == 289)then
!     
                        CALL VMX(1,PLQ7,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ1,C11,1)
                        CALL VMX(1,C11,DPLQ2,A11,1)
!     
                        IEEE=16
!     
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 13
!******************************************************************C
                        if (iddd == 290)then
!     
                        CALL VMX(1,PLQ2,PLQ3,B11,1)
                        CALL VMX(1,B11,PLQ4,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ1,C11,1)
                        CALL VMX(1,C11,DPLQ2,A11,1)
!     
                        IEEE=1
!     
                        endif 
!******************************************************************C
                        if (iddd == 291)then
!     
                        CALL VMX(1,PLQ3,PLQ4,B11,1)
                        CALL VMX(1,B11,PLQ1,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ2,C11,1)
                        CALL VMX(1,C11,DPLQ3,A11,1)
!     
                        IEEE=2
!     
                        endif 
!******************************************************************C
                        if (iddd == 292)then
!     
                        CALL VMX(1,PLQ4,PLQ1,B11,1)
                        CALL VMX(1,B11,PLQ2,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ3,C11,1)
                        CALL VMX(1,C11,DPLQ4,A11,1)
!     
                        IEEE=3
!     
                        endif 
!******************************************************************C
                        if (iddd == 293)then
!     
                        CALL VMX(1,PLQ1,PLQ2,B11,1)
                        CALL VMX(1,B11,PLQ3,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ4,C11,1)
                        CALL VMX(1,C11,DPLQ1,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
                        if (iddd == 294)then
!     
                        CALL VMX(1,PLQ5,PLQ6,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ7,B11,1)
                        CALL VMX(1,B11,DPLQ8,C11,1)
                        CALL VMX(1,C11,DPLQ5,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
                        if (iddd == 295)then
!     
                        CALL VMX(1,PLQ8,PLQ5,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ6,B11,1)
                        CALL VMX(1,B11,DPLQ7,C11,1)
                        CALL VMX(1,C11,DPLQ8,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
                        if (iddd == 296)then
!     
                        CALL VMX(1,PLQ7,PLQ8,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ5,B11,1)
                        CALL VMX(1,B11,DPLQ6,C11,1)
                        CALL VMX(1,C11,DPLQ7,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
                        if (iddd == 297)then
!     
                        CALL VMX(1,PLQ6,PLQ7,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ8,B11,1)
                        CALL VMX(1,B11,DPLQ5,C11,1)
                        CALL VMX(1,C11,DPLQ6,A11,1)
!     
                        IEEE=8
!     
                        endif 
!******************************************************************C
                        if (iddd == 298)then
!     
                        CALL VMX(1,PLQ6,PLQ7,B11,1)
                        CALL VMX(1,B11,PLQ8,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ5,C11,1)
                        CALL VMX(1,C11,DPLQ6,A11,1)
!     
                        IEEE=9
!     
                        endif 
!******************************************************************C
                        if (iddd == 299)then
!     
                        CALL VMX(1,PLQ7,PLQ8,B11,1)
                        CALL VMX(1,B11,PLQ5,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ6,C11,1)
                        CALL VMX(1,C11,DPLQ7,A11,1)
!     
                        IEEE=10
!     
                        endif 
!******************************************************************C
                        if (iddd == 300)then
!     
                        CALL VMX(1,PLQ8,PLQ5,B11,1)
                        CALL VMX(1,B11,PLQ6,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ7,C11,1)
                        CALL VMX(1,C11,DPLQ8,A11,1)
!     
                        IEEE=11
!     
                        endif 
!******************************************************************C
                        if (iddd == 301)then
!     
                        CALL VMX(1,PLQ5,PLQ6,B11,1)
                        CALL VMX(1,B11,PLQ7,C11,1)
                        CALL VMX(1,C11,PL,B11,1)
                        CALL VMX(1,B11,DPLQ8,C11,1)
                        CALL VMX(1,C11,DPLQ5,A11,1)
!     
                        IEEE=12
!     
                        endif 
!******************************************************************C
                        if (iddd == 302)then
!     
                        CALL VMX(1,PLQ1,PLQ2,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ3,B11,1)
                        CALL VMX(1,B11,DPLQ4,C11,1)
                        CALL VMX(1,C11,DPLQ1,A11,1)
!     
                        IEEE=13
!     
                        endif 
!******************************************************************C
                        if (iddd == 303)then
!     
                        CALL VMX(1,PLQ4,PLQ1,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ2,B11,1)
                        CALL VMX(1,B11,DPLQ3,C11,1)
                        CALL VMX(1,C11,DPLQ4,A11,1)
!     
                        IEEE=14
!     
                        endif 
!******************************************************************C
                        if (iddd == 304)then
!     
                        CALL VMX(1,PLQ3,PLQ4,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ1,B11,1)
                        CALL VMX(1,B11,DPLQ2,C11,1)
                        CALL VMX(1,C11,DPLQ3,A11,1)
!     
                        IEEE=15
!     
                        endif 
!******************************************************************C
                        if (iddd == 305)then
!     
                        CALL VMX(1,PLQ2,PLQ3,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ4,B11,1)
                        CALL VMX(1,B11,DPLQ1,C11,1)
                        CALL VMX(1,C11,DPLQ2,A11,1)
!     
                        IEEE=16
!     
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 14
!******************************************************************C
                        if (iddd == 306)then
!     
                        CALL VMX(1,PLQ2,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ7,A11,1)
!     
                        IEEE=1
!     
                        endif 
!******************************************************************C
                        if (iddd == 307)then
!     
                        CALL VMX(1,PLQ3,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ6,A11,1)
!     
                        IEEE=2
!     
                        endif 
!******************************************************************C
                        if (iddd == 308)then
!     
                        CALL VMX(1,PLQ4,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ5,A11,1)
!     
                        IEEE=3
!     
                        endif 
!******************************************************************C
                        if (iddd == 309)then
!     
                        CALL VMX(1,PLQ1,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ8,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
                        if (iddd == 310)then
!     
                        CALL VMX(1,PLQ4,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ5,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
                        if (iddd == 311)then
!     
                        CALL VMX(1,PLQ1,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ8,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
                        if (iddd == 312)then
!     
                        CALL VMX(1,PLQ2,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ7,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
                        if (iddd == 313)then
!     
                        CALL VMX(1,PLQ3,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ6,A11,1)
!     
                        IEEE=8
!     
                        endif 
!******************************************************************C
                        if (iddd == 314)then
!     
                        CALL VMX(1,PLQ6,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ3,A11,1)
!     
                        IEEE=9
!     
                        endif 
!******************************************************************C
                        if (iddd == 315)then
!     
                        CALL VMX(1,PLQ7,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ2,A11,1)
!     
                        IEEE=10
!     
                        endif 
!******************************************************************C
                        if (iddd == 316)then
!     
                        CALL VMX(1,PLQ8,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ1,A11,1)
!     
                        IEEE=11
!     
                        endif 
!******************************************************************C
                        if (iddd == 317)then
!     
                        CALL VMX(1,PLQ5,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ4,A11,1)
!     
                        IEEE=12
!     
                        endif 
!******************************************************************C
                        if (iddd == 318)then
!     
                        CALL VMX(1,PLQ8,SQDY1,B11,1)
                        CALL VMX(1,B11,DPLQ1,A11,1)
!     
                        IEEE=13
!     
                        endif 
!******************************************************************C
                        if (iddd == 319)then
!     
                        CALL VMX(1,PLQ5,SQUZ1,B11,1)
                        CALL VMX(1,B11,DPLQ4,A11,1)
!     
                        IEEE=14
!     
                        endif 
!******************************************************************C
                        if (iddd == 320)then
!     
                        CALL VMX(1,PLQ6,SQUY1,B11,1)
                        CALL VMX(1,B11,DPLQ3,A11,1)
!     
                        IEEE=15
!     
                        endif 
!******************************************************************C
                        if (iddd == 321)then
!     
                        CALL VMX(1,PLQ7,SQDZ1,B11,1)
                        CALL VMX(1,B11,DPLQ2,A11,1)
!     
                        IEEE=16
!     
                        endif 
!******************************************************************C
!     PLAQUETTE OPERATORS 15
!******************************************************************C
                        if (iddd == 322)then
!     
                        CALL VMX(1,PLQ2,PLQ7,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ6,A11,1)
!     
                        IEEE=1
!     
                        endif 
!******************************************************************C
                        if (iddd == 323)then
!     
                        CALL VMX(1,PLQ3,PLQ6,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ5,A11,1)
!     
                        IEEE=2
!     
                        endif 
!******************************************************************C
                        if (iddd == 324)then
!     
                        CALL VMX(1,PLQ4,PLQ5,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ8,A11,1)
!     
                        IEEE=3
!     
                        endif 
!******************************************************************C
                        if (iddd == 325)then
!     
                        CALL VMX(1,PLQ1,PLQ8,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ7,A11,1)
!     
                        IEEE=4
!     
                        endif 
!******************************************************************C
                        if (iddd == 326)then
!     
                        CALL VMX(1,PLQ1,PL,B11,1)
                        CALL VMX(1,B11,DPLQ4,C11,1)
                        CALL VMX(1,C11,DPLQ5,A11,1)
!     
                        IEEE=5
!     
                        endif 
!******************************************************************C
                        if (iddd == 327)then
!     
                        CALL VMX(1,PLQ2,PL,B11,1)
                        CALL VMX(1,B11,DPLQ1,C11,1)
                        CALL VMX(1,C11,DPLQ8,A11,1)
!     
                        IEEE=6
!     
                        endif 
!******************************************************************C
                        if (iddd == 328)then
!     
                        CALL VMX(1,PLQ3,PL,B11,1)
                        CALL VMX(1,B11,DPLQ2,C11,1)
                        CALL VMX(1,C11,DPLQ7,A11,1)
!     
                        IEEE=7
!     
                        endif 
!******************************************************************C
                        if (iddd == 329)then
!     
                        CALL VMX(1,PLQ4,PL,B11,1)
                        CALL VMX(1,B11,DPLQ3,C11,1)
                        CALL VMX(1,C11,DPLQ6,A11,1)
!     
                        IEEE=8
!     
                        endif 
!******************************************************************C
                        if (iddd == 330)then
!     
                        CALL VMX(1,PLQ6,PLQ3,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ2,A11,1)
!     
                        IEEE=9
!     
                        endif 
!******************************************************************C
                        if (iddd == 331)then
!     
                        CALL VMX(1,PLQ7,PLQ2,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ1,A11,1)
!     
                        IEEE=10
!     
                        endif 
!******************************************************************C
                        if (iddd == 332)then
!     
                        CALL VMX(1,PLQ8,PLQ1,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ4,A11,1)
!     
                        IEEE=11
!     
                        endif 
!******************************************************************C
                        if (iddd == 333)then
!     
                        CALL VMX(1,PLQ5,PLQ4,B11,1)
                        CALL VMX(1,B11,PL,C11,1)
                        CALL VMX(1,C11,DPLQ3,A11,1)
!     
                        IEEE=12
!     
                        endif 
!******************************************************************C
                        if (iddd == 334)then
!     
                        CALL VMX(1,PLQ5,PL,B11,1)
                        CALL VMX(1,B11,DPLQ8,C11,1)
                        CALL VMX(1,C11,DPLQ1,A11,1)
!     
                        IEEE=13
!     
                        endif 
!******************************************************************C
                        if (iddd == 335)then
!     
                        CALL VMX(1,PLQ6,PL,B11,1)
                        CALL VMX(1,B11,DPLQ5,C11,1)
                        CALL VMX(1,C11,DPLQ4,A11,1)
!     
                        IEEE=14
!     
                        endif 
!******************************************************************C
                        if (iddd == 336)then
!     
                        CALL VMX(1,PLQ7,PL,B11,1)
                        CALL VMX(1,B11,DPLQ6,C11,1)
                        CALL VMX(1,C11,DPLQ3,A11,1)
!     
                        IEEE=15
!     
                        endif 
!******************************************************************C
                        if (iddd == 337)then
!     
                        CALL VMX(1,PLQ8,PL,B11,1)
                        CALL VMX(1,B11,DPLQ7,C11,1)
                        CALL VMX(1,C11,DPLQ2,A11,1)
!     
                        IEEE=16
!     
                        endif                         
!******************************************************************c
        endif 
!******************************************************************C
! If the operator does not fit in this blocking level (more than once)
            if (lcnt(ids) < ico) then
                do 168 IC=1, NCOL2
                    A11(IC)=LIN0(IC)
168              continue
                M2=ML
            ELSE
                if (ico == 1) then
                    CALL VMX(1,A11,LIN1,C11,1)
                    do 169 IC=1,NCOL2
                        A11(IC)=C11(IC)
169                 continue
                endif 
!
                if (ico == 2) then
                    CALL VMX(1,A11,LIN2,C11,1)
                    do 170 IC=1,NCOL2
                        A11(IC)=C11(IC)
170                 continue
                endif 
!
                if (ico == 4) then
                    CALL VMX(1,A11,LIN4,C11,1)
                    do 171 IC=1,NCOL2
                        A11(IC)=C11(IC)
171                 continue
                endif 
                M2=ML
            endif 
!***********************************************************************
            endif 
            CALL VMX(1,A11,REM11,B11,1)
!***********************************************************************
!            CALL RENORMBS(B11,UREN11)
!            do IJ=1,NCOL2
!               AA(IJ)=UREN11(IJ)
!            enddo
!            JMAT=NCOL
!            CALL DETNANT(JMAT,DET,AA)
!c     
!            DNCOL=CMPLX(1.0/NCOL)
!            CDET=DET**DNCOL
!            CNORM=1.0/CDET
!c     
!            do IC=1,NCOL2
!               B11(IC)=UREN11(IC)*CNORM
!            enddo
!***********************************************************************
            AKT1 = cmplx(0.0, 0.0, kind=real64)
!***********************************************************************
            do 192 N1=1,NCOL
            IJ=N1+NCOL*(N1-1)
            AKT1=AKT1+B11(IJ)
192        continue
!************ ***********************************************************     
            if (iddd <= 4) then
            CSUMS(IEEE)=CSUMS(IEEE)+AKT1
            CSUMSMOM(IEEE,1)=CSUMSMOM(IEEE,1)+AKT1*PF(1)
            CSUMSMOM(IEEE,2)=CSUMSMOM(IEEE,2)+AKT1*PF(2)
!               CSUMSMOM(IEEE,3)=CSUMSMOM(IEEE,3)+AKT1*PF(3)
!               CSUMSMOM(IEEE,4)=CSUMSMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 4).and.(iddd < 9)) then
            CSUM2S(IEEE)=CSUM2S(IEEE)+AKT1
            CSUM2SMOM(IEEE,1)=CSUM2SMOM(IEEE,1)+AKT1*PF(1)
            CSUM2SMOM(IEEE,2)=CSUM2SMOM(IEEE,2)+AKT1*PF(2)
!               CSUM2SMOM(IEEE,3)=CSUM2SMOM(IEEE,3)+AKT1*PF(3)
!               CSUM2SMOM(IEEE,4)=CSUM2SMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 8).and.(iddd < 13)) then
            CSUM2WS(IEEE)=CSUM2WS(IEEE)+AKT1
            CSUM2WSMOM(IEEE,1)=CSUM2WSMOM(IEEE,1)+AKT1*PF(1)
            CSUM2WSMOM(IEEE,2)=CSUM2WSMOM(IEEE,2)+AKT1*PF(2)
!               CSUM2WSMOM(IEEE,3)=CSUM2WSMOM(IEEE,3)+AKT1*PF(3)
!               CSUM2WSMOM(IEEE,4)=CSUM2WSMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 12).and.(iddd < 17)) then
            CSUMW(IEEE)=CSUMW(IEEE)+AKT1
            CSUMWMOM(IEEE,1)=CSUMWMOM(IEEE,1)+AKT1*PF(1)
            CSUMWMOM(IEEE,2)=CSUMWMOM(IEEE,2)+AKT1*PF(2)
!               CSUMWMOM(IEEE,3)=CSUMWMOM(IEEE,3)+AKT1*PF(3)
!               CSUMWMOM(IEEE,4)=CSUMWMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 16).and.(iddd < 21)) then
            CSUM2W(IEEE)=CSUM2W(IEEE)+AKT1
            CSUM2WMOM(IEEE,1)=CSUM2WMOM(IEEE,1)+AKT1*PF(1)
            CSUM2WMOM(IEEE,2)=CSUM2WMOM(IEEE,2)+AKT1*PF(2)
!               CSUM2WMOM(IEEE,3)=CSUM2WMOM(IEEE,3)+AKT1*PF(3)
!               CSUM2WMOM(IEEE,4)=CSUM2WMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 20).and.(iddd < 25)) then
            CSUM3W(IEEE)=CSUM3W(IEEE)+AKT1
            CSUM3WMOM(IEEE,1)=CSUM3WMOM(IEEE,1)+AKT1*PF(1)
            CSUM3WMOM(IEEE,2)=CSUM3WMOM(IEEE,2)+AKT1*PF(2)
!               CSUM3WMOM(IEEE,3)=CSUM3WMOM(IEEE,3)+AKT1*PF(3)
!               CSUM3WMOM(IEEE,4)=CSUM3WMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 24).and.(iddd < 29)) then
            CSUMUP(IEEE)=CSUMUP(IEEE)+AKT1
            CSUMUPMOM(IEEE,1)=CSUMUPMOM(IEEE,1)+AKT1*PF(1)
            CSUMUPMOM(IEEE,2)=CSUMUPMOM(IEEE,2)+AKT1*PF(2)
!               CSUMUPMOM(IEEE,3)=CSUMUPMOM(IEEE,3)+AKT1*PF(3)
!               CSUMUPMOM(IEEE,4)=CSUMUPMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 28).and.(iddd < 33)) then
            CSUMUD(IEEE)=CSUMUD(IEEE)+AKT1
            CSUMUDMOM(IEEE,1)=CSUMUDMOM(IEEE,1)+AKT1*PF(1)
            CSUMUDMOM(IEEE,2)=CSUMUDMOM(IEEE,2)+AKT1*PF(2)
!               CSUMUDMOM(IEEE,3)=CSUMUDMOM(IEEE,3)+AKT1*PF(3)
!               CSUMUDMOM(IEEE,4)=CSUMUDMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 32).and.(iddd < 41)) then
            CSUMTT1(IEEE)=CSUMTT1(IEEE)+AKT1
            CSUMTTMOM1(IEEE,1)=CSUMTTMOM1(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM1(IEEE,2)=CSUMTTMOM1(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM1(IEEE,3)=CSUMTTMOM1(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM1(IEEE,4)=CSUMTTMOM1(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 40).and.(iddd < 49)) then
            CSUMTT2(IEEE)=CSUMTT2(IEEE)+AKT1
            CSUMTTMOM2(IEEE,1)=CSUMTTMOM2(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM2(IEEE,2)=CSUMTTMOM2(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM2(IEEE,3)=CSUMTTMOM2(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM2(IEEE,4)=CSUMTTMOM2(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 48).and.(iddd < 57)) then
            CSUMTT3(IEEE)=CSUMTT3(IEEE)+AKT1
            CSUMTTMOM3(IEEE,1)=CSUMTTMOM3(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM3(IEEE,2)=CSUMTTMOM3(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM3(IEEE,3)=CSUMTTMOM3(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM3(IEEE,4)=CSUMTTMOM3(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 56).and.(iddd < 65)) then
            CSUMTT4(IEEE)=CSUMTT4(IEEE)+AKT1
            CSUMTTMOM4(IEEE,1)=CSUMTTMOM4(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM4(IEEE,2)=CSUMTTMOM4(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM4(IEEE,3)=CSUMTTMOM4(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM4(IEEE,4)=CSUMTTMOM4(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 64).and.(iddd < 69)) then
            CSUMTT5(IEEE)=CSUMTT5(IEEE)+AKT1
            CSUMTTMOM5(IEEE,1)=CSUMTTMOM5(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM5(IEEE,2)=CSUMTTMOM5(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM5(IEEE,3)=CSUMTTMOM5(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM5(IEEE,4)=CSUMTTMOM5(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 68).and.(iddd < 73)) then
            CSUMTT6(IEEE)=CSUMTT6(IEEE)+AKT1
            CSUMTTMOM6(IEEE,1)=CSUMTTMOM6(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM6(IEEE,2)=CSUMTTMOM6(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM6(IEEE,3)=CSUMTTMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM6(IEEE,4)=CSUMTTMOM6(IEEE,4)+AKT1*PF(4)
            endif 
!            
            if ((iddd > 72).and.(iddd < 81)) then
            CSUMTT7(IEEE)=CSUMTT7(IEEE)+AKT1
            CSUMTTMOM7(IEEE,1)=CSUMTTMOM7(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM7(IEEE,2)=CSUMTTMOM7(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM7(IEEE,3)=CSUMTTMOM7(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM7(IEEE,4)=CSUMTTMOM7(IEEE,4)+AKT1*PF(4)
            endif 
!            
            if ((iddd > 80).and.(iddd < 89)) then
            CSUMTT8(IEEE)=CSUMTT8(IEEE)+AKT1
            CSUMTTMOM8(IEEE,1)=CSUMTTMOM8(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM8(IEEE,2)=CSUMTTMOM8(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM8(IEEE,3)=CSUMTTMOM8(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM8(IEEE,4)=CSUMTTMOM8(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 88).and.(iddd < 97)) then
            CSUMTT9(IEEE)=CSUMTT9(IEEE)+AKT1
            CSUMTTMOM9(IEEE,1)=CSUMTTMOM9(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM9(IEEE,2)=CSUMTTMOM9(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM9(IEEE,3)=CSUMTTMOM9(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM9(IEEE,4)=CSUMTTMOM9(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 96).and.(iddd < 105)) then
            CSUMTT10(IEEE)=CSUMTT10(IEEE)+AKT1
            CSUMTTMOM10(IEEE,1)=CSUMTTMOM10(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM10(IEEE,2)=CSUMTTMOM10(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM10(IEEE,3)=CSUMTTMOM10(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM10(IEEE,4)=CSUMTTMOM10(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 104).and.(iddd < 113)) then
            CSUMTT11(IEEE)=CSUMTT11(IEEE)+AKT1
            CSUMTTMOM11(IEEE,1)=CSUMTTMOM11(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM11(IEEE,2)=CSUMTTMOM11(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM11(IEEE,3)=CSUMTTMOM11(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM11(IEEE,4)=CSUMTTMOM11(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 112).and.(iddd < 129)) then
            CSUMTT12(IEEE)=CSUMTT12(IEEE)+AKT1
            CSUMTTMOM12(IEEE,1)=CSUMTTMOM12(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM12(IEEE,2)=CSUMTTMOM12(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM12(IEEE,3)=CSUMTTMOM12(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM12(IEEE,4)=CSUMTTMOM12(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 128).and.(iddd < 137)) then
            CSUMTT13(IEEE)=CSUMTT13(IEEE)+AKT1
            CSUMTTMOM13(IEEE,1)=CSUMTTMOM13(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM13(IEEE,2)=CSUMTTMOM13(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM13(IEEE,3)=CSUMTTMOM13(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM13(IEEE,4)=CSUMTTMOM13(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 136).and.(iddd < 145)) then
            CSUMTT14(IEEE)=CSUMTT14(IEEE)+AKT1
            CSUMTTMOM14(IEEE,1)=CSUMTTMOM14(IEEE,1)+AKT1*PF(1)
            CSUMTTMOM14(IEEE,2)=CSUMTTMOM14(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM14(IEEE,3)=CSUMTTMOM14(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM14(IEEE,4)=CSUMTTMOM14(IEEE,4)+AKT1*PF(4)
            endif 
!
            if (iddd == 145) then
            CSUMN=CSUMN+AKT1
            endif 
!
            if ((iddd > 145).and.(iddd < 154)) then
            CSUMPLQ(IEEE)=CSUMPLQ(IEEE)+AKT1
            CSUMPLQMOM(IEEE,1)=CSUMPLQMOM(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM(IEEE,2)=CSUMPLQMOM(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM(IEEE,3)=CSUMPLQMOM(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM(IEEE,4)=CSUMPLQMOM(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 153).and.(iddd < 162)) then
            CSUMPLQ2(IEEE)=CSUMPLQ2(IEEE)+AKT1
            CSUMPLQMOM2(IEEE,1)=CSUMPLQMOM2(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM2(IEEE,2)=CSUMPLQMOM2(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM2(IEEE,3)=CSUMPLQMOM2(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM2(IEEE,4)=CSUMPLQMOM2(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 161).and.(iddd < 170)) then
            CSUMPLQ3(IEEE)=CSUMPLQ3(IEEE)+AKT1
            CSUMPLQMOM3(IEEE,1)=CSUMPLQMOM3(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM3(IEEE,2)=CSUMPLQMOM3(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM3(IEEE,3)=CSUMPLQMOM3(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM3(IEEE,4)=CSUMPLQMOM3(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 169).and.(iddd < 178)) then
            CSUMPLQ4(IEEE)=CSUMPLQ4(IEEE)+AKT1
            CSUMPLQMOM4(IEEE,1)=CSUMPLQMOM4(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM4(IEEE,2)=CSUMPLQMOM4(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM4(IEEE,3)=CSUMPLQMOM4(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM4(IEEE,4)=CSUMPLQMOM4(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 177).and.(iddd < 186)) then
            CSUMPLQ5(IEEE)=CSUMPLQ5(IEEE)+AKT1
            CSUMPLQMOM5(IEEE,1)=CSUMPLQMOM5(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM5(IEEE,2)=CSUMPLQMOM5(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM5(IEEE,3)=CSUMPLQMOM5(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM5(IEEE,4)=CSUMPLQMOM5(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 185).and.(iddd < 194)) then
            CSUMPLQ6(IEEE)=CSUMPLQ6(IEEE)+AKT1
            CSUMPLQMOM6(IEEE,1)=CSUMPLQMOM6(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM6(IEEE,2)=CSUMPLQMOM6(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 193).and.(iddd < 210)) then
            CSUMPLQ7(IEEE)=CSUMPLQ7(IEEE)+AKT1
            CSUMPLQMOM7(IEEE,1)=CSUMPLQMOM7(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM7(IEEE,2)=CSUMPLQMOM7(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 209).and.(iddd < 226)) then 
            CSUMPLQ8(IEEE)=CSUMPLQ8(IEEE)+AKT1
            CSUMPLQMOM8(IEEE,1)=CSUMPLQMOM8(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM8(IEEE,2)=CSUMPLQMOM8(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 225).and.(iddd < 242)) then 
            CSUMPLQ9(IEEE)=CSUMPLQ9(IEEE)+AKT1
            CSUMPLQMOM9(IEEE,1)=CSUMPLQMOM9(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM9(IEEE,2)=CSUMPLQMOM9(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 241).and.(iddd < 258)) then 
            CSUMPLQ10(IEEE)=CSUMPLQ10(IEEE)+AKT1
            CSUMPLQMOM10(IEEE,1)=CSUMPLQMOM10(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM10(IEEE,2)=CSUMPLQMOM10(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif    
!
            if ((iddd > 257).and.(iddd < 274)) then 
            CSUMPLQ11(IEEE)=CSUMPLQ11(IEEE)+AKT1
            CSUMPLQMOM11(IEEE,1)=CSUMPLQMOM11(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM11(IEEE,2)=CSUMPLQMOM11(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif 
!
            if ((iddd > 273).and.(iddd < 290)) then 
            CSUMPLQ12(IEEE)=CSUMPLQ12(IEEE)+AKT1
            CSUMPLQMOM12(IEEE,1)=CSUMPLQMOM12(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM12(IEEE,2)=CSUMPLQMOM12(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif   
!
            if ((iddd > 289).and.(iddd < 306)) then 
            CSUMPLQ13(IEEE)=CSUMPLQ13(IEEE)+AKT1
            CSUMPLQMOM13(IEEE,1)=CSUMPLQMOM13(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM13(IEEE,2)=CSUMPLQMOM13(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif   
!
            if ((iddd > 305).and.(iddd < 322)) then 
            CSUMPLQ14(IEEE)=CSUMPLQ14(IEEE)+AKT1
            CSUMPLQMOM14(IEEE,1)=CSUMPLQMOM14(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM14(IEEE,2)=CSUMPLQMOM14(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif   
!
            if ((iddd > 321).and.(iddd < 338)) then 
            CSUMPLQ15(IEEE)=CSUMPLQ15(IEEE)+AKT1
            CSUMPLQMOM15(IEEE,1)=CSUMPLQMOM15(IEEE,1)+AKT1*PF(1)
            CSUMPLQMOM15(IEEE,2)=CSUMPLQMOM15(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            endif   
    !     
    9       continue ! OPERATOR CONSTRUCTION LOOP end
    !
        enddo
        enddo
        enddo
    !**********************************************************************     
        ADIV1=1.0/(4.0*LS(KU))
        ADIV2=1.0/(8.0*LS(KU))
        ADIV3=1.0/(16.0*LS(KU))
        ADIVN=1.0/LS(KU)
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE1(N4,ID)=CSUMN*ADIVN
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE2(N4,ID)=(CSUMS(1)+CSUMS(2)+CSUMS(3)+CSUMS(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM2(N4,ID,1)=(CSUMSMOM(1,1)+CSUMSMOM(2,1)&
        +CSUMSMOM(3,1)+CSUMSMOM(4,1))*ADIV1
        ALINEMOM2(N4,ID,2)=(CSUMSMOM(1,2)+CSUMSMOM(2,2)&
        +CSUMSMOM(3,2)+CSUMSMOM(4,2))*ADIV1
    !      ALINEMOM2(N4,ID,3)=(CSUMSMOM(1,3)+CSUMSMOM(2,3)&
    !     &+CSUMSMOM(3,3)+CSUMSMOM(4,3))*ADIV1
    !      ALINEMOM2(N4,ID,4)=(CSUMSMOM(1,4)+CSUMSMOM(2,4)&
    !     &+CSUMSMOM(3,4)+CSUMSMOM(4,4))*ADIV1
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE3(N4,ID)=(CSUMS(1)+GIOT*CSUMS(2)-CSUMS(3)-GIOT*CSUMS(4))&
        *ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM3(N4,ID,1)=(CSUMSMOM(1,1)+GIOT*CSUMSMOM(2,1)&
        -CSUMSMOM(3,1)-GIOT*CSUMSMOM(4,1))*ADIV1
        ALINEMOM3(N4,ID,2)=(CSUMSMOM(1,2)+GIOT*CSUMSMOM(2,2)&
        -CSUMSMOM(3,2)-GIOT*CSUMSMOM(4,2))*ADIV1
    !      ALINEMOM3(N4,ID,3)=(CSUMSMOM(1,3)+GIOT*CSUMSMOM(2,3)&
    !     &-CSUMSMOM(3,3)-GIOT*CSUMSMOM(4,3))*ADIV1
    !      ALINEMOM3(N4,ID,4)=(CSUMSMOM(1,4)+GIOT*CSUMSMOM(2,4)&
    !     &-CSUMSMOM(3,4)-GIOT*CSUMSMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE4(N4,ID)=(CSUMS(1)-CSUMS(2)+CSUMS(3)-CSUMS(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM4(N4,ID,1)=(CSUMSMOM(1,1)-CSUMSMOM(2,1)&
        +CSUMSMOM(3,1)-CSUMSMOM(4,1))*ADIV1
        ALINEMOM4(N4,ID,2)=(CSUMSMOM(1,2)-CSUMSMOM(2,2)&
        +CSUMSMOM(3,2)-CSUMSMOM(4,2))*ADIV1
    !      ALINEMOM4(N4,ID,3)=(CSUMSMOM(1,3)-CSUMSMOM(2,3)&
    !     &+CSUMSMOM(3,3)-CSUMSMOM(4,3))*ADIV1
    !      ALINEMOM4(N4,ID,4)=(CSUMSMOM(1,4)-CSUMSMOM(2,4)&
    !     &+CSUMSMOM(3,4)-CSUMSMOM(4,4))*ADIV1
    !**********************************************************************
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE5(N4,ID)=(CSUM2S(1)+CSUM2S(2)+CSUM2S(3)+CSUM2S(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM5(N4,ID,1)=(CSUM2SMOM(1,1)+CSUM2SMOM(2,1)&
        +CSUM2SMOM(3,1)+CSUM2SMOM(4,1))*ADIV1
        ALINEMOM5(N4,ID,2)=(CSUM2SMOM(1,2)+CSUM2SMOM(2,2)&
        +CSUM2SMOM(3,2)+CSUM2SMOM(4,2))*ADIV1
    !      ALINEMOM5(N4,ID,3)=(CSUM2SMOM(1,3)+CSUM2SMOM(2,3)&
    !     &+CSUM2SMOM(3,3)+CSUM2SMOM(4,3))*ADIV1
    !      ALINEMOM5(N4,ID,4)=(CSUM2SMOM(1,4)+CSUM2SMOM(2,4)&
    !     &+CSUM2SMOM(3,4)+CSUM2SMOM(4,4))*ADIV1
    !***********************************************************************
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE6(N4,ID)=(CSUM2S(1)+GIOT*CSUM2S(2)-CSUM2S(3)&
        -GIOT*CSUM2S(4))*ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM6(N4,ID,1)=(CSUM2SMOM(1,1)+GIOT*CSUM2SMOM(2,1)&
        -CSUM2SMOM(3,1)-GIOT*CSUM2SMOM(4,1))*ADIV1
        ALINEMOM6(N4,ID,2)=(CSUM2SMOM(1,2)+GIOT*CSUM2SMOM(2,2)&
        -CSUM2SMOM(3,2)-GIOT*CSUM2SMOM(4,2))*ADIV1
    !      ALINEMOM6(N4,ID,3)=(CSUM2SMOM(1,3)+GIOT*CSUM2SMOM(2,3)&
    !     &-CSUM2SMOM(3,3)-GIOT*CSUM2SMOM(4,3))*ADIV1
    !      ALINEMOM6(N4,ID,4)=(CSUM2SMOM(1,4)+GIOT*CSUM2SMOM(2,4)&
    !     &-CSUM2SMOM(3,4)-GIOT*CSUM2SMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE7(N4,ID)=(CSUM2S(1)-CSUM2S(2)+CSUM2S(3)-CSUM2S(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM7(N4,ID,1)=(CSUM2SMOM(1,1)-CSUM2SMOM(2,1)&
        +CSUM2SMOM(3,1)-CSUM2SMOM(4,1))*ADIV1
        ALINEMOM7(N4,ID,2)=(CSUM2SMOM(1,2)-CSUM2SMOM(2,2)&
        +CSUM2SMOM(3,2)-CSUM2SMOM(4,2))*ADIV1
    !      ALINEMOM7(N4,ID,3)=(CSUM2SMOM(1,3)-CSUM2SMOM(2,3)&
    !     &+CSUM2SMOM(3,3)-CSUM2SMOM(4,3))*ADIV1
    !      ALINEMOM7(N4,ID,4)=(CSUM2SMOM(1,4)-CSUM2SMOM(2,4)&
    !     &+CSUM2SMOM(3,4)-CSUM2SMOM(4,4))*ADIV1
    !**********************************************************************
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE8(N4,ID)=(CSUM2WS(1)+CSUM2WS(2)+CSUM2WS(3)+CSUM2WS(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM8(N4,ID,1)=(CSUM2WSMOM(1,1)+CSUM2WSMOM(2,1)&
        +CSUM2WSMOM(3,1)+CSUM2WSMOM(4,1))*ADIV1
        ALINEMOM8(N4,ID,2)=(CSUM2WSMOM(1,2)+CSUM2WSMOM(2,2)&
        +CSUM2WSMOM(3,2)+CSUM2WSMOM(4,2))*ADIV1
    !      ALINEMOM8(N4,ID,3)=(CSUM2WSMOM(1,3)+CSUM2WSMOM(2,3)&
    !     &+CSUM2WSMOM(3,3)+CSUM2WSMOM(4,3))*ADIV1
    !      ALINEMOM8(N4,ID,4)=(CSUM2WSMOM(1,4)+CSUM2WSMOM(2,4)&
    !     &+CSUM2WSMOM(3,4)+CSUM2WSMOM(4,4))*ADIV1
    !************************************************************************
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE9(N4,ID)=(CSUM2WS(1)+GIOT*CSUM2WS(2)-CSUM2WS(3)&
        -GIOT*CSUM2WS(4))*ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM9(N4,ID,1)=(CSUM2WSMOM(1,1)+GIOT*CSUM2WSMOM(2,1)&
        -CSUM2WSMOM(3,1)-GIOT*CSUM2WSMOM(4,1))*ADIV1
        ALINEMOM9(N4,ID,2)=(CSUM2WSMOM(1,2)+GIOT*CSUM2WSMOM(2,2)&
        -CSUM2WSMOM(3,2)-GIOT*CSUM2WSMOM(4,2))*ADIV1
    !      ALINEMOM9(N4,ID,3)=(CSUM2WSMOM(1,3)+GIOT*CSUM2WSMOM(2,3)&
    !     &-CSUM2WSMOM(3,3)-GIOT*CSUM2WSMOM(4,3))*ADIV1
    !      ALINEMOM9(N4,ID,4)=(CSUM2WSMOM(1,4)+GIOT*CSUM2WSMOM(2,4)&
    !     &-CSUM2WSMOM(3,4)-GIOT*CSUM2WSMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE10(N4,ID)=(CSUM2WS(1)-CSUM2WS(2)+CSUM2WS(3)-CSUM2WS(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM10(N4,ID,1)=(CSUM2WSMOM(1,1)-CSUM2WSMOM(2,1)&
        +CSUM2WSMOM(3,1)-CSUM2WSMOM(4,1))*ADIV1
        ALINEMOM10(N4,ID,2)=(CSUM2WSMOM(1,2)-CSUM2WSMOM(2,2)&
        +CSUM2WSMOM(3,2)-CSUM2WSMOM(4,2))*ADIV1
    !      ALINEMOM10(N4,ID,3)=(CSUM2WSMOM(1,3)-CSUM2WSMOM(2,3)&
    !     &+CSUM2WSMOM(3,3)-CSUM2WSMOM(4,3))*ADIV1
    !      ALINEMOM10(N4,ID,4)=(CSUM2WSMOM(1,4)-CSUM2WSMOM(2,4)&
    !     &+CSUM2WSMOM(3,4)-CSUM2WSMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE11(N4,ID)=(CSUMW(1)+CSUMW(2)+CSUMW(3)+CSUMW(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM11(N4,ID,1)=(CSUMWMOM(1,1)+CSUMWMOM(2,1)&
        +CSUMWMOM(3,1)+CSUMWMOM(4,1))*ADIV1
        ALINEMOM11(N4,ID,2)=(CSUMWMOM(1,2)+CSUMWMOM(2,2)&
        +CSUMWMOM(3,2)+CSUMWMOM(4,2))*ADIV1
    !      ALINEMOM11(N4,ID,3)=(CSUMWMOM(1,3)+CSUMWMOM(2,3)&
    !     &+CSUMWMOM(3,3)+CSUMWMOM(4,3))*ADIV1
    !      ALINEMOM11(N4,ID,4)=(CSUMWMOM(1,4)+CSUMWMOM(2,4)&
    !     &+CSUMWMOM(3,4)+CSUMWMOM(4,4))*ADIV1
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !***********************************************************************
        ALINE12(N4,ID)=(CSUMW(1)+GIOT*CSUMW(2)-CSUMW(3)&
        -GIOT*CSUMW(4))*ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM12(N4,ID,1)=(CSUMWMOM(1,1)+GIOT*CSUMWMOM(2,1)&
        -CSUMWMOM(3,1)-GIOT*CSUMWMOM(4,1))*ADIV1
        ALINEMOM12(N4,ID,2)=(CSUMWMOM(1,2)+GIOT*CSUMWMOM(2,2)&
        -CSUMWMOM(3,2)-GIOT*CSUMWMOM(4,2))*ADIV1
    !      ALINEMOM12(N4,ID,3)=(CSUMWMOM(1,3)+GIOT*CSUMWMOM(2,3)&
    !     &-CSUMWMOM(3,3)-GIOT*CSUMWMOM(4,3))*ADIV1
    !      ALINEMOM12(N4,ID,4)=(CSUMWMOM(1,4)+GIOT*CSUMWMOM(2,4)&
    !     &-CSUMWMOM(3,4)-GIOT*CSUMWMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE13(N4,ID)=(CSUMW(1)-CSUMW(2)+CSUMW(3)-CSUMW(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM13(N4,ID,1)=(CSUMWMOM(1,1)-CSUMWMOM(2,1)&
        +CSUMWMOM(3,1)-CSUMWMOM(4,1))*ADIV1
        ALINEMOM13(N4,ID,2)=(CSUMWMOM(1,2)-CSUMWMOM(2,2)&
        +CSUMWMOM(3,2)-CSUMWMOM(4,2))*ADIV1
    !      ALINEMOM13(N4,ID,3)=(CSUMWMOM(1,3)-CSUMWMOM(2,3)&
    !     &+CSUMWMOM(3,3)-CSUMWMOM(4,3))*ADIV1
    !      ALINEMOM13(N4,ID,4)=(CSUMWMOM(1,4)-CSUMWMOM(2,4)&
    !     &+CSUMWMOM(3,4)-CSUMWMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE14(N4,ID)=(CSUM2W(1)+CSUM2W(2)+CSUM2W(3)+CSUM2W(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM14(N4,ID,1)=(CSUM2WMOM(1,1)+CSUM2WMOM(2,1)&
        +CSUM2WMOM(3,1)+CSUM2WMOM(4,1))*ADIV1
        ALINEMOM14(N4,ID,2)=(CSUM2WMOM(1,2)+CSUM2WMOM(2,2)&
        +CSUM2WMOM(3,2)+CSUM2WMOM(4,2))*ADIV1
    !      ALINEMOM14(N4,ID,3)=(CSUM2WMOM(1,3)+CSUM2WMOM(2,3)&
    !     &+CSUM2WMOM(3,3)+CSUM2WMOM(4,3))*ADIV1
    !      ALINEMOM14(N4,ID,4)=(CSUM2WMOM(1,4)+CSUM2WMOM(2,4)&
    !     &+CSUM2WMOM(3,4)+CSUM2WMOM(4,4))*ADIV1
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE15(N4,ID)=(CSUM2W(1)+GIOT*CSUM2W(2)-CSUM2W(3)&
        -GIOT*CSUM2W(4))*ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM15(N4,ID,1)=(CSUM2WMOM(1,1)+GIOT*CSUM2WMOM(2,1)&
        -CSUM2WMOM(3,1)-GIOT*CSUM2WMOM(4,1))*ADIV1
        ALINEMOM15(N4,ID,2)=(CSUM2WMOM(1,2)+GIOT*CSUM2WMOM(2,2)&
        -CSUM2WMOM(3,2)-GIOT*CSUM2WMOM(4,2))*ADIV1
    !      ALINEMOM15(N4,ID,3)=(CSUM2WMOM(1,3)+GIOT*CSUM2WMOM(2,3)&
    !     &-CSUM2WMOM(3,3)-GIOT*CSUM2WMOM(4,3))*ADIV1
    !      ALINEMOM15(N4,ID,4)=(CSUM2WMOM(1,4)+GIOT*CSUM2WMOM(2,4)&
    !     &-CSUM2WMOM(3,4)-GIOT*CSUM2WMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE16(N4,ID)=(CSUM2W(1)-CSUM2W(2)+CSUM2W(3)-CSUM2W(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM16(N4,ID,1)=(CSUM2WMOM(1,1)-CSUM2WMOM(2,1)&
        +CSUM2WMOM(3,1)-CSUM2WMOM(4,1))*ADIV1
        ALINEMOM16(N4,ID,2)=(CSUM2WMOM(1,2)-CSUM2WMOM(2,2)&
        +CSUM2WMOM(3,2)-CSUM2WMOM(4,2))*ADIV1
    !      ALINEMOM16(N4,ID,3)=(CSUM2WMOM(1,3)-CSUM2WMOM(2,3)&
    !     &+CSUM2WMOM(3,3)-CSUM2WMOM(4,3))*ADIV1
    !      ALINEMOM16(N4,ID,4)=(CSUM2WMOM(1,4)-CSUM2WMOM(2,4)&
    !     &+CSUM2WMOM(3,4)-CSUM2WMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE17(N4,ID)=(CSUM3W(1)+CSUM3W(2)+CSUM3W(3)+CSUM3W(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM17(N4,ID,1)=(CSUM3WMOM(1,1)+CSUM3WMOM(2,1)&
        +CSUM3WMOM(3,1)+CSUM3WMOM(4,1))*ADIV1
        ALINEMOM17(N4,ID,2)=(CSUM3WMOM(1,2)+CSUM3WMOM(2,2)&
        +CSUM3WMOM(3,2)+CSUM3WMOM(4,2))*ADIV1
    !      ALINEMOM17(N4,ID,3)=(CSUM3WMOM(1,3)+CSUM3WMOM(2,3)&
    !     &+CSUM3WMOM(3,3)+CSUM3WMOM(4,3))*ADIV1
    !      ALINEMOM17(N4,ID,4)=(CSUM3WMOM(1,4)+CSUM3WMOM(2,4)&
    !     &+CSUM3WMOM(3,4)+CSUM3WMOM(4,4))*ADIV1
    !***********************************************************************
    !     J=1, Pr=+ q=0
    !**********************************************************************
        ALINE18(N4,ID)=(CSUM3W(1)+GIOT*CSUM3W(2)-CSUM3W(3)-GIOT*CSUM3W(4))&
        *ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM18(N4,ID,1)=(CSUM3WMOM(1,1)+GIOT*CSUM3WMOM(2,1)&
        -CSUM3WMOM(3,1)-GIOT*CSUM3WMOM(4,1))*ADIV1
        ALINEMOM18(N4,ID,2)=(CSUM3WMOM(1,2)+GIOT*CSUM3WMOM(2,2)&
        -CSUM3WMOM(3,2)-GIOT*CSUM3WMOM(4,2))*ADIV1
    !      ALINEMOM18(N4,ID,3)=(CSUM3WMOM(1,3)+GIOT*CSUM3WMOM(2,3)&
    !     &-CSUM3WMOM(3,3)-GIOT*CSUM3WMOM(4,3))*ADIV1
    !      ALINEMOM18(N4,ID,4)=(CSUM3WMOM(1,4)+GIOT*CSUM3WMOM(2,4)&
    !     &-CSUM3WMOM(3,4)-GIOT*CSUM3WMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE19(N4,ID)=(CSUM3W(1)-CSUM3W(2)+CSUM3W(3)-CSUM3W(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM19(N4,ID,1)=(CSUM3WMOM(1,1)-CSUM3WMOM(2,1)&
        +CSUM3WMOM(3,1)-CSUM3WMOM(4,1))*ADIV1
        ALINEMOM19(N4,ID,2)=(CSUM3WMOM(1,2)-CSUM3WMOM(2,2)&
        +CSUM3WMOM(3,2)-CSUM3WMOM(4,2))*ADIV1
    !      ALINEMOM19(N4,ID,3)=(CSUM3WMOM(1,3)-CSUM3WMOM(2,3)&
    !     &+CSUM3WMOM(3,3)-CSUM3WMOM(4,3))*ADIV1
    !      ALINEMOM19(N4,ID,4)=(CSUM3WMOM(1,4)-CSUM3WMOM(2,4)&
    !     &+CSUM3WMOM(3,4)-CSUM3WMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE20(N4,ID)=(CSUMUP(1)+CSUMUP(2)+CSUMUP(3)+CSUMUP(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM20(N4,ID,1)=(CSUMUPMOM(1,1)+CSUMUPMOM(2,1)&
        +CSUMUPMOM(3,1)+CSUMUPMOM(4,1))*ADIV1
        ALINEMOM20(N4,ID,2)=(CSUMUPMOM(1,2)+CSUMUPMOM(2,2)&
        +CSUMUPMOM(3,2)+CSUMUPMOM(4,2))*ADIV1
    !      ALINEMOM20(N4,ID,3)=(CSUMUPMOM(1,3)+CSUMUPMOM(2,3)&
    !     &+CSUMUPMOM(3,3)+CSUMUPMOM(4,3))*ADIV1
    !      ALINEMOM20(N4,ID,4)=(CSUMUPMOM(1,4)+CSUMUPMOM(2,4)&
    !     &+CSUMUPMOM(3,4)+CSUMUPMOM(4,4))*ADIV1
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE21(N4,ID)=(CSUMUP(1)+GIOT*CSUMUP(2)-CSUMUP(3)-GIOT*CSUMUP(4))&
        *ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM21(N4,ID,1)=(CSUMUPMOM(1,1)+GIOT*CSUMUPMOM(2,1)&
        -CSUMUPMOM(3,1)-GIOT*CSUMUPMOM(4,1))*ADIV1
        ALINEMOM21(N4,ID,2)=(CSUMUPMOM(1,2)+GIOT*CSUMUPMOM(2,2)&
        -CSUMUPMOM(3,2)-GIOT*CSUMUPMOM(4,2))*ADIV1
    !      ALINEMOM21(N4,ID,3)=(CSUMUPMOM(1,3)+GIOT*CSUMUPMOM(2,3)&
    !     &-CSUMUPMOM(3,3)-GIOT*CSUMUPMOM(4,3))*ADIV1
    !      ALINEMOM21(N4,ID,4)=(CSUMUPMOM(1,4)+GIOT*CSUMUPMOM(2,4)&
    !     &-CSUMUPMOM(3,4)-GIOT*CSUMUPMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE22(N4,ID)=(CSUMUP(1)-CSUMUP(2)+CSUMUP(3)-CSUMUP(4))*ADIV1
    !***********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !***********************************************************************
        ALINEMOM22(N4,ID,1)=(CSUMUPMOM(1,1)-CSUMUPMOM(2,1)&
        +CSUMUPMOM(3,1)-CSUMUPMOM(4,1))*ADIV1
        ALINEMOM22(N4,ID,2)=(CSUMUPMOM(1,2)-CSUMUPMOM(2,2)&
        +CSUMUPMOM(3,2)-CSUMUPMOM(4,2))*ADIV1
    !      ALINEMOM22(N4,ID,3)=(CSUMUPMOM(1,3)-CSUMUPMOM(2,3)&
    !     &+CSUMUPMOM(3,3)-CSUMUPMOM(4,3))*ADIV1
    !      ALINEMOM22(N4,ID,4)=(CSUMUPMOM(1,4)-CSUMUPMOM(2,4)&
    !     &+CSUMUPMOM(3,4)-CSUMUPMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE23(N4,ID)=(CSUMUD(1)+CSUMUD(2)+CSUMUD(3)+CSUMUD(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM23(N4,ID,1)=(CSUMUDMOM(1,1)+CSUMUDMOM(2,1)&
        +CSUMUDMOM(3,1)+CSUMUDMOM(4,1))*ADIV1
        ALINEMOM23(N4,ID,2)=(CSUMUDMOM(1,2)+CSUMUDMOM(2,2)&
        +CSUMUDMOM(3,2)+CSUMUDMOM(4,2))*ADIV1
    !      ALINEMOM23(N4,ID,3)=(CSUMUDMOM(1,3)+CSUMUDMOM(2,3)&
    !     &+CSUMUDMOM(3,3)+CSUMUDMOM(4,3))*ADIV1
    !      ALINEMOM23(N4,ID,4)=(CSUMUDMOM(1,4)+CSUMUDMOM(2,4)&
    !     &+CSUMUDMOM(3,4)+CSUMUDMOM(4,4))*ADIV1
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE24(N4,ID)=(CSUMUD(1)+GIOT*CSUMUD(2)-CSUMUD(3)&
        -GIOT*CSUMUD(4))*ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !***********************************************************************
        ALINEMOM24(N4,ID,1)=(CSUMUDMOM(1,1)+GIOT*CSUMUDMOM(2,1)&
        -CSUMUDMOM(3,1)-GIOT*CSUMUDMOM(4,1))*ADIV1
        ALINEMOM24(N4,ID,2)=(CSUMUDMOM(1,2)+GIOT*CSUMUDMOM(2,2)&
        -CSUMUDMOM(3,2)-GIOT*CSUMUDMOM(4,2))*ADIV1
    !      ALINEMOM24(N4,ID,3)=(CSUMUDMOM(1,3)+GIOT*CSUMUDMOM(2,3)&
    !     &-CSUMUDMOM(3,3)-GIOT*CSUMUDMOM(4,3))*ADIV1
    !      ALINEMOM24(N4,ID,4)=(CSUMUDMOM(1,4)+GIOT*CSUMUDMOM(2,4)&
    !     &-CSUMUDMOM(3,4)-GIOT*CSUMUDMOM(4,4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE25(N4,ID)=(CSUMUD(1)-CSUMUD(2)+CSUMUD(3)-CSUMUD(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        ALINEMOM25(N4,ID,1)=(CSUMUDMOM(1,1)-CSUMUDMOM(2,1)&
        +CSUMUDMOM(3,1)-CSUMUDMOM(4,1))*ADIV1
        ALINEMOM25(N4,ID,2)=(CSUMUDMOM(1,2)-CSUMUDMOM(2,2)&
        +CSUMUDMOM(3,2)-CSUMUDMOM(4,2))*ADIV1
    !      ALINEMOM25(N4,ID,3)=(CSUMUDMOM(1,3)-CSUMUDMOM(2,3)&
    !     &+CSUMUDMOM(3,3)-CSUMUDMOM(4,3))*ADIV1
    !      ALINEMOM25(N4,ID,4)=(CSUMUDMOM(1,4)-CSUMUDMOM(2,4)&
    !     &+CSUMUDMOM(3,4)-CSUMUDMOM(4,4))*ADIV1
    !**********************************************************************
    !**********************************************************************
    !**********************************************************************
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE26(N4,ID)=(CSUMTT1(1)+CSUMTT1(2)+CSUMTT1(3)+CSUMTT1(4)+&
        (CSUMTT1(7)+CSUMTT1(6)+CSUMTT1(5)+CSUMTT1(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM26(N4,ID,IK)=(CSUMTTMOM1(1,IK)+CSUMTTMOM1(2,IK)&
        +CSUMTTMOM1(3,IK)+CSUMTTMOM1(4,IK)+CSUMTTMOM1(7,IK)&
        +CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)&
        +CSUMTTMOM1(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE27(N4,ID)=(CSUMTT1(1)+CSUMTT1(2)+CSUMTT1(3)+CSUMTT1(4)-&
        (CSUMTT1(7)+CSUMTT1(6)+CSUMTT1(5)+CSUMTT1(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM27(N4,ID,IK)=(CSUMTTMOM1(1,IK)+CSUMTTMOM1(2,IK)&
        +CSUMTTMOM1(3,IK)+CSUMTTMOM1(4,IK)-(CSUMTTMOM1(7,IK)&
        +CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)&
        +CSUMTTMOM1(8,IK)))*ADIV2
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE28(N4,ID)=(CSUMTT1(1)+GIOT*CSUMTT1(2)-CSUMTT1(3)&
        -GIOT*CSUMTT1(4)+(CSUMTT1(6)+GIOT*CSUMTT1(7)-CSUMTT1(8)&
        -GIOT*CSUMTT1(5)))*ADIV2
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM28(N4,ID,IK)=(CSUMTTMOM1(1,IK)+GIOT*CSUMTTMOM1(2,IK)&
        -CSUMTTMOM1(3,IK)-GIOT*CSUMTTMOM1(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE29(N4,ID)=(CSUMTT1(1)+GIOT*CSUMTT1(2)-CSUMTT1(3)&
        -GIOT*CSUMTT1(4)-(CSUMTT1(6)+GIOT*CSUMTT1(7)-CSUMTT1(8)&
        -GIOT*CSUMTT1(5)))*ADIV2
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM29(N4,ID,IK)=(CSUMTTMOM1(6,IK)+GIOT*CSUMTTMOM1(7,IK)&
        -CSUMTTMOM1(8,IK)-GIOT*CSUMTTMOM1(5,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE30(N4,ID)=(CSUMTT1(1)-CSUMTT1(2)+CSUMTT1(3)-CSUMTT1(4)+&
        CSUMTT1(7)-CSUMTT1(6)+CSUMTT1(5)-CSUMTT1(8))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM30(N4,ID,IK)=(CSUMTTMOM1(1,IK)-CSUMTTMOM1(2,IK)&
        +CSUMTTMOM1(3,IK)-CSUMTTMOM1(4,IK)+CSUMTTMOM1(7,IK)&
        -CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)-CSUMTTMOM1(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE31(N4,ID)=(CSUMTT1(1)-CSUMTT1(2)+CSUMTT1(3)-CSUMTT1(4)-&
        (CSUMTT1(7)-CSUMTT1(6)+CSUMTT1(5)-CSUMTT1(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM31(N4,ID,IK)=(CSUMTTMOM1(1,IK)-CSUMTTMOM1(2,IK)&
        +CSUMTTMOM1(3,IK)-CSUMTTMOM1(4,IK)-(CSUMTTMOM1(7,IK)&
        -CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)-CSUMTTMOM1(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE32(N4,ID)=(CSUMTT2(1)+CSUMTT2(2)+CSUMTT2(3)+CSUMTT2(4)+&
        CSUMTT2(5)+CSUMTT2(6)+CSUMTT2(7)+CSUMTT2(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM32(N4,ID,IK)=(CSUMTTMOM2(1,IK)+CSUMTTMOM2(2,IK)&
        +CSUMTTMOM2(3,IK)+CSUMTTMOM2(4,IK)+CSUMTTMOM2(5,IK)&
        +CSUMTTMOM2(6,IK)+CSUMTTMOM2(7,IK)&
        +CSUMTTMOM2(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE33(N4,ID)=(CSUMTT2(1)+CSUMTT2(2)+CSUMTT2(3)+CSUMTT2(4)-&
        (CSUMTT2(5)+CSUMTT2(6)+CSUMTT2(7)+CSUMTT2(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q
    !**********************************************************************
        do IK=1,2
        ALINEMOM33(N4,ID,IK)=(CSUMTTMOM2(1,IK)+CSUMTTMOM2(2,IK)&
        +CSUMTTMOM2(3,IK)+CSUMTTMOM2(4,IK)-(CSUMTTMOM2(5,IK)&
        +CSUMTTMOM2(6,IK)+CSUMTTMOM2(7,IK)&
        +CSUMTTMOM2(8,IK)))*ADIV2
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE34(N4,ID)=(CSUMTT2(1)+GIOT*CSUMTT2(2)-CSUMTT2(3)&
        -GIOT*CSUMTT2(4)+(CSUMTT2(7)+GIOT*CSUMTT2(8)-CSUMTT2(5)&
        -GIOT*CSUMTT2(6)))*ADIV2
    !***********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM34(N4,ID,IK)=(CSUMTTMOM2(1,IK)+GIOT*CSUMTTMOM2(2,IK)&
        -CSUMTTMOM2(3,IK)-GIOT*CSUMTTMOM2(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE35(N4,ID)=(CSUMTT2(1)+GIOT*CSUMTT2(2)-CSUMTT2(3)&
        -GIOT*CSUMTT2(4)-(CSUMTT2(7)+GIOT*CSUMTT2(8)-CSUMTT2(5)&
        -GIOT*CSUMTT2(6)))*ADIV2
    !***********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM35(N4,ID,IK)=(CSUMTTMOM2(7,IK)+GIOT*CSUMTTMOM2(8,IK)&
        -CSUMTTMOM2(5,IK)-GIOT*CSUMTTMOM2(6,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE36(N4,ID)=(CSUMTT2(1)-CSUMTT2(2)+CSUMTT2(3)-CSUMTT2(4)+&
        CSUMTT2(6)-CSUMTT2(5)+CSUMTT2(8)-CSUMTT2(7))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM36(N4,ID,IK)=(CSUMTTMOM2(1,IK)-CSUMTTMOM2(2,IK)&
        +CSUMTTMOM2(3,IK)-CSUMTTMOM2(4,IK)+CSUMTTMOM2(6,IK)&
        -CSUMTTMOM2(5,IK)+CSUMTTMOM2(8,IK)&
        -CSUMTTMOM2(7,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE37(N4,ID)=(CSUMTT2(1)-CSUMTT2(2)+CSUMTT2(3)-CSUMTT2(4)-&
        (CSUMTT2(6)-CSUMTT2(5)+CSUMTT2(8)-CSUMTT2(7)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM37(N4,ID,IK)=(CSUMTTMOM2(1,IK)-CSUMTTMOM2(2,IK)&
        +CSUMTTMOM2(3,IK)-CSUMTTMOM2(4,IK)-(CSUMTTMOM2(6,IK)&
        -CSUMTTMOM2(5,IK)+CSUMTTMOM2(8,IK)&
        -CSUMTTMOM2(7,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE38(N4,ID)=(CSUMTT3(1)+CSUMTT3(2)+CSUMTT3(3)+CSUMTT3(4)+&
        CSUMTT3(5)+CSUMTT3(6)+CSUMTT3(7)+CSUMTT3(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM38(N4,ID,IK)=(CSUMTTMOM3(1,IK)+CSUMTTMOM3(2,IK)&
        +CSUMTTMOM3(3,IK)+CSUMTTMOM3(4,IK)+CSUMTTMOM3(5,IK)&
        +CSUMTTMOM3(6,IK)+CSUMTTMOM3(7,IK)&
        +CSUMTTMOM3(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE39(N4,ID)=(CSUMTT3(1)+CSUMTT3(2)+CSUMTT3(3)+CSUMTT3(4)-&
        (CSUMTT3(5)+CSUMTT3(6)+CSUMTT3(7)+CSUMTT3(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM39(N4,ID,IK)=(CSUMTTMOM3(1,IK)+CSUMTTMOM3(2,IK)&
        +CSUMTTMOM3(3,IK)+CSUMTTMOM3(4,IK)-(CSUMTTMOM3(5,IK)&
        +CSUMTTMOM3(6,IK)+CSUMTTMOM3(7,IK)&
        +CSUMTTMOM3(8,IK)))*ADIV2
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE40(N4,ID)=(CSUMTT3(1)+GIOT*CSUMTT3(2)-CSUMTT3(3)&
        -GIOT*CSUMTT3(4)+(CSUMTT3(6)+GIOT*CSUMTT3(7)-CSUMTT3(8)&
        -GIOT*CSUMTT3(5)))*ADIV2
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM40(N4,ID,IK)=(CSUMTTMOM3(1,IK)+GIOT*CSUMTTMOM3(2,IK)&
        -CSUMTTMOM3(3,IK)-GIOT*CSUMTTMOM3(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE41(N4,ID)=(CSUMTT3(1)+GIOT*CSUMTT3(2)-CSUMTT3(3)&
        -GIOT*CSUMTT3(4)-(CSUMTT3(6)+GIOT*CSUMTT3(7)-CSUMTT3(8)&
        -GIOT*CSUMTT3(5)))*ADIV2
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM41(N4,ID,IK)=(CSUMTTMOM3(6,IK)+GIOT*CSUMTTMOM3(7,IK)&
        -CSUMTTMOM3(8,IK)-GIOT*CSUMTTMOM3(5,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE42(N4,ID)=(CSUMTT3(1)-CSUMTT3(2)+CSUMTT3(3)-CSUMTT3(4)+&
        (CSUMTT3(7)-CSUMTT3(6)+CSUMTT3(5)-CSUMTT3(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM42(N4,ID,IK)=(CSUMTTMOM3(1,IK)-CSUMTTMOM3(2,IK)&
        +CSUMTTMOM3(3,IK)-CSUMTTMOM3(4,IK)+(CSUMTTMOM3(7,IK)&
        -CSUMTTMOM3(6,IK)+CSUMTTMOM3(5,IK)-CSUMTTMOM3(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE43(N4,ID)=(CSUMTT3(1)-CSUMTT3(2)+CSUMTT3(3)-CSUMTT3(4)-&
        (CSUMTT3(7)-CSUMTT3(6)+CSUMTT3(5)-CSUMTT3(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM43(N4,ID,IK)=(CSUMTTMOM3(1,IK)-CSUMTTMOM3(2,IK)&
        +CSUMTTMOM3(3,IK)-CSUMTTMOM3(4,IK)-(CSUMTTMOM3(7,IK)&
        -CSUMTTMOM3(6,IK)+CSUMTTMOM3(5,IK)-CSUMTTMOM3(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE44(N4,ID)=(CSUMTT4(1)+CSUMTT4(2)+CSUMTT4(3)+CSUMTT4(4)+&
        CSUMTT4(5)+CSUMTT4(6)+CSUMTT4(7)+CSUMTT4(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM44(N4,ID,IK)=(CSUMTTMOM4(1,IK)+CSUMTTMOM4(2,IK)&
        +CSUMTTMOM4(3,IK)+CSUMTTMOM4(4,IK)+CSUMTTMOM4(5,IK)&
        +CSUMTTMOM4(6,IK)+CSUMTTMOM4(7,IK)&
        +CSUMTTMOM4(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE45(N4,ID)=(CSUMTT4(1)+CSUMTT4(2)+CSUMTT4(3)+CSUMTT4(4)-&
        (CSUMTT4(5)+CSUMTT4(6)+CSUMTT4(7)+CSUMTT4(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM45(N4,ID,IK)=(CSUMTTMOM4(1,IK)+CSUMTTMOM4(2,IK)&
        +CSUMTTMOM4(3,IK)+CSUMTTMOM4(4,IK)-(CSUMTTMOM4(5,IK)&
        +CSUMTTMOM4(6,IK)+CSUMTTMOM4(7,IK)&
        +CSUMTTMOM4(8,IK)))*ADIV2
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE46(N4,ID)=(CSUMTT4(1)+GIOT*CSUMTT4(2)-CSUMTT4(3)&
        -GIOT*CSUMTT4(4)+(CSUMTT4(7)+GIOT*CSUMTT4(8)-CSUMTT4(5)&
        -GIOT*CSUMTT4(6)))*ADIV2
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM46(N4,ID,IK)=(CSUMTTMOM4(1,IK)+GIOT*CSUMTTMOM4(2,IK)&
        -CSUMTTMOM4(3,IK)-GIOT*CSUMTTMOM4(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE47(N4,ID)=(CSUMTT4(1)+GIOT*CSUMTT4(2)-CSUMTT4(3)&
        -GIOT*CSUMTT4(4)-(CSUMTT4(7)+GIOT*CSUMTT4(8)-CSUMTT4(5)&
        -GIOT*CSUMTT4(6)))*ADIV2
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM47(N4,ID,IK)=(CSUMTTMOM4(7,IK)+GIOT*CSUMTTMOM4(8,IK)&
        -CSUMTTMOM4(5,IK)-GIOT*CSUMTTMOM4(6,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE48(N4,ID)=(CSUMTT4(1)-CSUMTT4(2)+CSUMTT4(3)-CSUMTT4(4)+&
        CSUMTT4(8)-CSUMTT4(7)+CSUMTT4(6)-CSUMTT4(5))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM48(N4,ID,IK)=(CSUMTTMOM4(1,IK)-CSUMTTMOM4(2,IK)&
        +CSUMTTMOM4(3,IK)-CSUMTTMOM4(4,IK)+CSUMTTMOM4(8,IK)&
        -CSUMTTMOM4(7,IK)+CSUMTTMOM4(6,IK)-CSUMTTMOM4(5,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE49(N4,ID)=(CSUMTT4(1)-CSUMTT4(2)+CSUMTT4(3)-CSUMTT4(4)-&
        (CSUMTT4(8)-CSUMTT4(7)+CSUMTT4(6)-CSUMTT4(5)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM49(N4,ID,IK)=(CSUMTTMOM4(1,IK)-CSUMTTMOM4(2,IK)&
        +CSUMTTMOM4(3,IK)-CSUMTTMOM4(4,IK)-(CSUMTTMOM4(8,IK)&
        -CSUMTTMOM4(7,IK)+CSUMTTMOM4(6,IK)-CSUMTTMOM4(5,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-5 OPERATORS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE50(N4,ID)=(CSUMTT5(1)+CSUMTT5(2)+CSUMTT5(3)+CSUMTT5(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM50(N4,ID,IK)=(CSUMTTMOM5(1,IK)+CSUMTTMOM5(2,IK)&
        +CSUMTTMOM5(3,IK)+CSUMTTMOM5(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE51(N4,ID)=(CSUMTT5(1)+GIOT*CSUMTT5(2)-CSUMTT5(3)&
        -GIOT*CSUMTT5(4))*ADIV1
    !***********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM51(N4,ID,IK)=(CSUMTTMOM5(1,IK)+GIOT*CSUMTTMOM5(2,IK)&
        -CSUMTTMOM5(3,IK)-GIOT*CSUMTTMOM5(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE52(N4,ID)=(CSUMTT5(1)-CSUMTT5(2)+CSUMTT5(3)-CSUMTT5(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM52(N4,ID,IK)=(CSUMTTMOM5(1,IK)-CSUMTTMOM5(2,IK)&
        +CSUMTTMOM5(3,IK)-CSUMTTMOM5(4,IK))*ADIV1
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! 6
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-6 OPERATORS 
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE53(N4,ID)=(CSUMTT6(1)+CSUMTT6(2)+CSUMTT6(3)+CSUMTT6(4))*ADIV1
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM53(N4,ID,IK)=(CSUMTTMOM6(1,IK)+CSUMTTMOM6(2,IK)&
        +CSUMTTMOM6(3,IK)+CSUMTTMOM6(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE54(N4,ID)=(CSUMTT6(1)+GIOT*CSUMTT6(2)-CSUMTT6(3)&
        -GIOT*CSUMTT6(4))*ADIV1
    !***********************************************************************
    !     J=1, q
    !**********************************************************************
        do IK=1,2
        ALINEMOM54(N4,ID,IK)=(CSUMTTMOM6(1,IK)+GIOT*CSUMTTMOM6(2,IK)&
        -CSUMTTMOM6(3,IK)-GIOT*CSUMTTMOM6(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE55(N4,ID)=(CSUMTT6(1)-CSUMTT6(2)+CSUMTT6(3)-CSUMTT6(4))*ADIV1
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM55(N4,ID,IK)=(CSUMTTMOM6(1,IK)-CSUMTTMOM6(2,IK)&
        +CSUMTTMOM6(3,IK)-CSUMTTMOM6(4,IK))*ADIV1
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! 7
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-7 OPERATORS 
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE56(N4,ID)=(CSUMTT7(1)+CSUMTT7(2)+CSUMTT7(3)+CSUMTT7(4)+&
        (CSUMTT7(5)+CSUMTT7(6)+CSUMTT7(7)+CSUMTT7(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM56(N4,ID,IK)=(CSUMTTMOM7(1,IK)+CSUMTTMOM7(2,IK)&
        +CSUMTTMOM7(3,IK)+CSUMTTMOM7(4,IK)+(CSUMTTMOM7(5,IK)&
        +CSUMTTMOM7(6,IK)+CSUMTTMOM7(7,IK)&
        +CSUMTTMOM7(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE57(N4,ID)=(CSUMTT7(1)+CSUMTT7(2)+CSUMTT7(3)+CSUMTT7(4)-&
        (CSUMTT7(5)+CSUMTT7(6)+CSUMTT7(7)+CSUMTT7(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM57(N4,ID,IK)=(CSUMTTMOM7(1,IK)+CSUMTTMOM7(2,IK)&
        +CSUMTTMOM7(3,IK)+CSUMTTMOM7(4,IK)-(CSUMTTMOM7(5,IK)&
        +CSUMTTMOM7(6,IK)+CSUMTTMOM7(7,IK)&
        +CSUMTTMOM7(8,IK)))*ADIV2
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE58(N4,ID)=(CSUMTT7(1)+GIOT*CSUMTT7(2)-CSUMTT7(3)&
        -GIOT*CSUMTT7(4)+&
        (CSUMTT7(7)+GIOT*CSUMTT7(8)-CSUMTT7(5)-GIOT*CSUMTT7(6)))*ADIV2
    !***********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM58(N4,ID,IK)=(CSUMTTMOM7(1,IK)+GIOT*CSUMTTMOM7(2,IK)&
        -CSUMTTMOM7(3,IK)-GIOT*CSUMTTMOM7(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE59(N4,ID)=(CSUMTT7(1)+GIOT*CSUMTT7(2)-CSUMTT7(3)&
        -GIOT*CSUMTT7(4)-&
        (CSUMTT7(7)+GIOT*CSUMTT7(8)-CSUMTT7(5)-GIOT*CSUMTT7(6)))*ADIV2
    !***********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM59(N4,ID,IK)=(CSUMTTMOM7(7,IK)+GIOT*CSUMTTMOM7(8,IK)&
        -CSUMTTMOM7(5,IK)-GIOT*CSUMTTMOM7(6,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE60(N4,ID)=(CSUMTT7(1)-CSUMTT7(2)+CSUMTT7(3)-CSUMTT7(4)+&
        (CSUMTT7(8)-CSUMTT7(7)+CSUMTT7(6)-CSUMTT7(5)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM60(N4,ID,IK)=(CSUMTTMOM7(1,IK)-CSUMTTMOM7(2,IK)&
        +CSUMTTMOM7(3,IK)-CSUMTTMOM7(4,IK)+(CSUMTTMOM7(8,IK)&
        -CSUMTTMOM7(7,IK)+CSUMTTMOM7(6,IK)&
        -CSUMTTMOM7(5,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE61(N4,ID)=(CSUMTT7(1)-CSUMTT7(2)+CSUMTT7(3)-CSUMTT7(4)-&
        (CSUMTT7(8)-CSUMTT7(7)+CSUMTT7(6)-CSUMTT7(5)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM61(N4,ID,IK)=(CSUMTTMOM7(1,IK)-CSUMTTMOM7(2,IK)&
        +CSUMTTMOM7(3,IK)-CSUMTTMOM7(4,IK)-(CSUMTTMOM7(8,IK)&
        -CSUMTTMOM7(7,IK)+CSUMTTMOM7(6,IK)&
        -CSUMTTMOM7(5,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! 8
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-8 OPERATORS CP=+,Pr=+
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE62(N4,ID)=(CSUMTT8(1)+CSUMTT8(2)+CSUMTT8(3)+CSUMTT8(4)&
        +(CSUMTT8(5)+CSUMTT8(6)+CSUMTT8(7)+CSUMTT8(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM62(N4,ID,IK)=(CSUMTTMOM8(1,IK)+CSUMTTMOM8(2,IK)&
        +CSUMTTMOM8(3,IK)+CSUMTTMOM8(4,IK)+(CSUMTTMOM8(5,IK)&
        +CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)&
        +CSUMTTMOM8(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE63(N4,ID)=(CSUMTT8(1)+CSUMTT8(2)+CSUMTT8(3)+CSUMTT8(4)&
        -(CSUMTT8(5)+CSUMTT8(6)+CSUMTT8(7)+CSUMTT8(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM63(N4,ID,IK)=(CSUMTTMOM8(1,IK)+CSUMTTMOM8(2,IK)&
        +CSUMTTMOM8(3,IK)+CSUMTTMOM8(4,IK)-(CSUMTTMOM8(5,IK)&
        +CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)&
        +CSUMTTMOM8(8,IK)))*ADIV2
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE64(N4,ID)=(CSUMTT8(1)+GIOT*CSUMTT8(2)-CSUMTT8(3)&
        -GIOT*CSUMTT8(4))*ADIV1
    !***********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
            ALINEMOM64(N4,ID,IK)=(CSUMTTMOM8(1,IK)+GIOT*CSUMTTMOM8(2,IK)&
        -CSUMTTMOM8(3,IK)-GIOT*CSUMTTMOM8(4,IK))*ADIV1
        enddo
    !***********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE65(N4,ID)=(CSUMTT8(7)+GIOT*CSUMTT8(6)-CSUMTT8(5)&
        -GIOT*CSUMTT8(8))*ADIV1
    !***********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
            ALINEMOM65(N4,ID,IK)=(CSUMTTMOM8(7,IK)+GIOT*CSUMTTMOM8(6,IK)&
        -CSUMTTMOM8(5,IK)-GIOT*CSUMTTMOM8(8,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE66(N4,ID)=(CSUMTT8(1)-CSUMTT8(2)+CSUMTT8(3)-CSUMTT8(4)&
        +(CSUMTT8(5)-CSUMTT8(6)+CSUMTT8(7)-CSUMTT8(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM66(N4,ID,IK)=(CSUMTTMOM8(1,IK)-CSUMTTMOM8(2,IK)&
        +CSUMTTMOM8(3,IK)-CSUMTTMOM8(4,IK)+(CSUMTTMOM8(5,IK)&
        -CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)&
        -CSUMTTMOM8(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE67(N4,ID)=(CSUMTT8(1)-CSUMTT8(2)+CSUMTT8(3)-CSUMTT8(4)&
        -(CSUMTT8(5)-CSUMTT8(6)+CSUMTT8(7)-CSUMTT8(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM67(N4,ID,IK)=(CSUMTTMOM8(1,IK)-CSUMTTMOM8(2,IK)&
        +CSUMTTMOM8(3,IK)-CSUMTTMOM8(4,IK)-(CSUMTTMOM8(5,IK)&
        -CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)&
        -CSUMTTMOM8(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! 9
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-9 OPERATORS 
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE68(N4,ID)=(CSUMTT9(1)+CSUMTT9(2)+CSUMTT9(3)+CSUMTT9(4)+&
        (CSUMTT9(5)+CSUMTT9(6)+CSUMTT9(7)+CSUMTT9(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1,2
            ALINEMOM68(N4,ID,IK)=(CSUMTTMOM9(1,IK)+CSUMTTMOM9(2,IK)+&
        CSUMTTMOM9(3,IK)+CSUMTTMOM9(4,IK)&
        +(CSUMTTMOM9(5,IK)+CSUMTTMOM9(6,IK)+CSUMTTMOM9(7,IK)&
        +CSUMTTMOM9(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-9 OPERATORS CP=-,Pz=-,J=0
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE69(N4,ID)=(CSUMTT9(1)+CSUMTT9(2)+CSUMTT9(3)+CSUMTT9(4)&
        -(CSUMTT9(5)+CSUMTT9(6)+CSUMTT9(7)+CSUMTT9(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1,2
            ALINEMOM69(N4,ID,IK)=(CSUMTTMOM9(1,IK)+CSUMTTMOM9(2,IK)+&
        CSUMTTMOM9(3,IK)+CSUMTTMOM9(4,IK)&
        -(CSUMTTMOM9(5,IK)+CSUMTTMOM9(6,IK)+CSUMTTMOM9(7,IK)&
        +CSUMTTMOM9(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************      
        ALINE70(N4,ID)=(CSUMTT9(1)+GIOT*CSUMTT9(2)-CSUMTT9(3)&
        -GIOT*CSUMTT9(4)+(CSUMTT9(7)+GIOT*CSUMTT9(8)-CSUMTT9(5)&
        -GIOT*CSUMTT9(6)))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************      
        do IK=1,2
            ALINEMOM70(N4,ID,IK)=(CSUMTTMOM9(1,IK)+GIOT*CSUMTTMOM9(2,IK)&
        -CSUMTTMOM9(3,IK)-GIOT*CSUMTTMOM9(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************            
        ALINE71(N4,ID)=(CSUMTT9(1)+GIOT*CSUMTT9(2)-CSUMTT9(3)&
        -GIOT*CSUMTT9(4)-(CSUMTT9(7)+GIOT*CSUMTT9(8)-CSUMTT9(5)&
        -GIOT*CSUMTT9(6)))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************            
        do IK=1,2
            ALINEMOM71(N4,ID,IK)=(CSUMTTMOM9(7,IK)+GIOT*CSUMTTMOM9(8,IK)&
        -CSUMTTMOM9(5,IK)-GIOT*CSUMTTMOM9(6,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************            
        ALINE72(N4,ID)=(CSUMTT9(1)-CSUMTT9(2)+CSUMTT9(3)-CSUMTT9(4)&
        +(CSUMTT9(7)-CSUMTT9(6)+CSUMTT9(5)-CSUMTT9(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2
    !**********************************************************************            
        do IK=1,2
            ALINEMOM72(N4,ID,IK)=(CSUMTTMOM9(1,IK)-CSUMTTMOM9(2,IK)&
        +CSUMTTMOM9(3,IK)-CSUMTTMOM9(4,IK)&
        +(CSUMTTMOM9(7,IK)-CSUMTTMOM9(6,IK)+CSUMTTMOM9(5,IK)&
        -CSUMTTMOM9(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=-, q=0 !Here!
    !**********************************************************************                  
        ALINE73(N4,ID)=(CSUMTT9(1)-CSUMTT9(2)+CSUMTT9(3)-CSUMTT9(4)&
        -(CSUMTT9(7)-CSUMTT9(6)+CSUMTT9(5)-CSUMTT9(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2
    !**********************************************************************                  
        do IK=1,2
            ALINEMOM73(N4,ID,IK)=(CSUMTTMOM9(1,IK)-CSUMTTMOM9(2,IK)&
        +CSUMTTMOM9(3,IK)-CSUMTTMOM9(4,IK)&
        -(CSUMTTMOM9(7,IK)-CSUMTTMOM9(6,IK)+CSUMTTMOM9(5,IK)&
        -CSUMTTMOM9(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-10 OPERATORS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************                  
        ALINE74(N4,ID)=(CSUMTT10(1)+CSUMTT10(2)+CSUMTT10(3)+CSUMTT10(4)&
        +(CSUMTT10(5)+CSUMTT10(6)+CSUMTT10(7)+CSUMTT10(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2
    !**********************************************************************                  
        do IK=1,2
            ALINEMOM74(N4,ID,IK)=(CSUMTTMOM10(1,IK)+CSUMTTMOM10(2,IK)&
        +CSUMTTMOM10(3,IK)+CSUMTTMOM10(4,IK)&
        +(CSUMTTMOM10(5,IK)+CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)&
        +CSUMTTMOM10(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************                        
        ALINE75(N4,ID)=(CSUMTT10(1)+CSUMTT10(2)+CSUMTT10(3)+CSUMTT10(4)&
        -(CSUMTT10(5)+CSUMTT10(6)+CSUMTT10(7)+CSUMTT10(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2
    !**********************************************************************                        
        do IK=1,2
            ALINEMOM75(N4,ID,IK)=(CSUMTTMOM10(1,IK)+CSUMTTMOM10(2,IK)&
        +CSUMTTMOM10(3,IK)+CSUMTTMOM10(4,IK)&
        -(CSUMTTMOM10(5,IK)+CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)&
        +CSUMTTMOM10(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************                  
        ALINE76(N4,ID)=(CSUMTT10(1)+GIOT*CSUMTT10(2)-CSUMTT10(3)&
        -GIOT*CSUMTT10(4)+(CSUMTT10(7)+GIOT*CSUMTT10(6)-CSUMTT10(5)&
        -GIOT*CSUMTT10(8)))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1,2
            ALINEMOM76(N4,ID,IK)=(CSUMTTMOM10(1,IK)+GIOT*CSUMTTMOM10(2,IK)&
        -CSUMTTMOM10(3,IK)-GIOT*CSUMTTMOM10(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE77(N4,ID)=(CSUMTT10(1)+GIOT*CSUMTT10(2)-CSUMTT10(3)&
        -GIOT*CSUMTT10(4)-(CSUMTT10(7)+GIOT*CSUMTT10(6)-CSUMTT10(5)&
        -GIOT*CSUMTT10(8)))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************
        do IK=1,2
            ALINEMOM77(N4,ID,IK)=(CSUMTTMOM10(7,IK)+GIOT*CSUMTTMOM10(6,IK)&
        -CSUMTTMOM10(5,IK)-GIOT*CSUMTTMOM10(8,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE78(N4,ID)=(CSUMTT10(1)-CSUMTT10(2)+CSUMTT10(3)-CSUMTT10(4)&
        +(CSUMTT10(5)-CSUMTT10(6)+CSUMTT10(7)-CSUMTT10(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2
    !**********************************************************************
        do IK=1,2
            ALINEMOM78(N4,ID,IK)=(CSUMTTMOM10(1,IK)-CSUMTTMOM10(2,IK)&
        +CSUMTTMOM10(3,IK)-CSUMTTMOM10(4,IK)&
        +(CSUMTTMOM10(5,IK)-CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)&
        -CSUMTTMOM10(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************      
        ALINE79(N4,ID)=(CSUMTT10(1)-CSUMTT10(2)+CSUMTT10(3)-CSUMTT10(4)&
        -(CSUMTT10(5)-CSUMTT10(6)+CSUMTT10(7)-CSUMTT10(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2
    !**********************************************************************      
        do IK=1,2
            ALINEMOM79(N4,ID,IK)=(CSUMTTMOM10(1,IK)-CSUMTTMOM10(2,IK)&
        +CSUMTTMOM10(3,IK)-CSUMTTMOM10(4,IK)&
        -(CSUMTTMOM10(5,IK)-CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)&
        -CSUMTTMOM10(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-11 OPERATORS
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************            
        ALINE80(N4,ID)=(CSUMTT11(1)+CSUMTT11(2)+CSUMTT11(3)+CSUMTT11(4)&
        +(CSUMTT11(5)+CSUMTT11(6)+CSUMTT11(7)+CSUMTT11(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2
    !**********************************************************************            
        do IK=1,2
            ALINEMOM80(N4,ID,IK)=(CSUMTTMOM11(1,IK)+CSUMTTMOM11(2,IK)&
        +CSUMTTMOM11(3,IK)+CSUMTTMOM11(4,IK)&
        +(CSUMTTMOM11(5,IK)+CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)&
        +CSUMTTMOM11(8,IK)))*ADIV2
        enddo      
    !**********************************************************************
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************            
        ALINE81(N4,ID)=(CSUMTT11(1)+CSUMTT11(2)+CSUMTT11(3)+CSUMTT11(4)&
        -(CSUMTT11(5)+CSUMTT11(6)+CSUMTT11(7)+CSUMTT11(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2
    !**********************************************************************            
        do IK=1,2
            ALINEMOM81(N4,ID,IK)=(CSUMTTMOM11(1,IK)+CSUMTTMOM11(2,IK)&
        +CSUMTTMOM11(3,IK)+CSUMTTMOM11(4,IK)&
        -(CSUMTTMOM11(5,IK)+CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)&
        +CSUMTTMOM11(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE82(N4,ID)=(CSUMTT11(1)+GIOT*CSUMTT11(2)-CSUMTT11(3)&
        -GIOT*CSUMTT11(4))*ADIV1
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1,2
            ALINEMOM82(N4,ID,IK)=(CSUMTTMOM11(1,IK)+GIOT*CSUMTTMOM11(2,IK)&
        -CSUMTTMOM11(3,IK)-GIOT*CSUMTTMOM11(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************                        
        ALINE83(N4,ID)=(CSUMTT11(7)+GIOT*CSUMTT11(6)-CSUMTT11(5)&
        -GIOT*CSUMTT11(8))*ADIV1
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1,2
            ALINEMOM83(N4,ID,IK)=(CSUMTTMOM11(7,IK)+GIOT*CSUMTTMOM11(6,IK)&
        -CSUMTTMOM11(5,IK)-GIOT*CSUMTTMOM11(8,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE84(N4,ID)=(CSUMTT11(1)-CSUMTT11(2)+CSUMTT11(3)-CSUMTT11(4)&
        +(CSUMTT11(5)-CSUMTT11(6)+CSUMTT11(7)-CSUMTT11(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2
    !**********************************************************************                  
        do IK=1,2
            ALINEMOM84(N4,ID,IK)=(CSUMTTMOM11(1,IK)-CSUMTTMOM11(2,IK)&
        +CSUMTTMOM11(3,IK)-CSUMTTMOM11(4,IK)&
        +(CSUMTTMOM11(5,IK)-CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)&
        -CSUMTTMOM11(8,IK)))*ADIV2
        enddo      
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************                        
        ALINE85(N4,ID)=(CSUMTT11(1)-CSUMTT11(2)+CSUMTT11(3)-CSUMTT11(4)&
        -(CSUMTT11(5)-CSUMTT11(6)+CSUMTT11(7)-CSUMTT11(8)))*ADIV2
    !c**********************************************************************
    !     J=2, Pp=-, q=1,2
    !**********************************************************************                        
        do IK=1,2
            ALINEMOM85(N4,ID,IK)=(CSUMTTMOM11(1,IK)-CSUMTTMOM11(2,IK)&
        +CSUMTTMOM11(3,IK)-CSUMTTMOM11(4,IK)&
        -(CSUMTTMOM11(5,IK)-CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)&
        -CSUMTTMOM11(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     12 THIS!!!!!
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-12 OPERATORS 
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE86(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)&
        +(CSUMTT12(5)+CSUMTT12(6)+CSUMTT12(7)+CSUMTT12(8))&
        +(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))&
        +(CSUMTT12(13)+CSUMTT12(14)+CSUMTT12(15)+CSUMTT12(16)))*ADIV3
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM86(N4,ID,IK)=(CSUMTTMOM12(1,IK)+CSUMTTMOM12(2,IK)&
        +CSUMTTMOM12(3,IK)+CSUMTTMOM12(4,IK)&
        +(CSUMTTMOM12(9,IK)+CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)&
        +CSUMTTMOM12(12,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE87(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)&
        +(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))&
        -(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6))&
        -(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM87(N4,ID,IK)=(CSUMTTMOM12(7,IK)+CSUMTTMOM12(8,IK)&
        +CSUMTTMOM12(5,IK)+CSUMTTMOM12(6,IK)+(CSUMTTMOM12(15,IK)&
        +CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)+CSUMTTMOM12(14,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE88(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)&
        -(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))&
        +(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6))&
        -(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM88(N4,ID,IK)=(CSUMTTMOM12(1,IK)+CSUMTTMOM12(2,IK)&
        +CSUMTTMOM12(3,IK)+CSUMTTMOM12(4,IK)-(CSUMTTMOM12(9,IK)&
        +CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)+CSUMTTMOM12(12,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE89(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)&
        -(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))&
        -(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6))&
        +(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM89(N4,ID,IK)=(CSUMTTMOM12(7,IK)+CSUMTTMOM12(8,IK)&
        +CSUMTTMOM12(5,IK)+CSUMTTMOM12(6,IK)-(CSUMTTMOM12(15,IK)&
        +CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)+CSUMTTMOM12(14,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE90(N4,ID)=(CSUMTT12(1)+GIOT*CSUMTT12(2)-CSUMTT12(3)&
        -GIOT*CSUMTT12(4)&
        +(CSUMTT12(7)+GIOT*CSUMTT12(8)-CSUMTT12(5)-GIOT*CSUMTT12(6))&
        )*ADIV2
    !**********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM90(N4,ID,IK)=(CSUMTTMOM12(1,IK)+GIOT*CSUMTTMOM12(2,IK)&
        -CSUMTTMOM12(3,IK)-GIOT*CSUMTTMOM12(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE91(N4,ID)=(CSUMTT12(9)+GIOT*CSUMTT12(10)-CSUMTT12(11)&
        -GIOT*CSUMTT12(12)&
        +(CSUMTT12(15)+GIOT*CSUMTT12(16)-CSUMTT12(13)-GIOT*CSUMTT12(14))&
        )*ADIV2
    !**********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM91(N4,ID,IK)=(CSUMTTMOM12(9,IK)+GIOT*CSUMTTMOM12(10,IK)&
        -CSUMTTMOM12(11,IK)-GIOT*CSUMTTMOM12(12,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE92(N4,ID)=(CSUMTT12(1)+GIOT*CSUMTT12(2)-CSUMTT12(3)&
        -GIOT*CSUMTT12(4)&
        -(CSUMTT12(7)+GIOT*CSUMTT12(8)-CSUMTT12(5)-GIOT*CSUMTT12(6))&
        )*ADIV2
    !**********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM92(N4,ID,IK)=(CSUMTTMOM12(7,IK)+GIOT*CSUMTTMOM12(8,IK)&
        -CSUMTTMOM12(5,IK)-GIOT*CSUMTTMOM12(6,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE93(N4,ID)=(CSUMTT12(9)+GIOT*CSUMTT12(10)-CSUMTT12(11)&
        -GIOT*CSUMTT12(12)&
        -(CSUMTT12(15)+GIOT*CSUMTT12(16)-CSUMTT12(13)-GIOT*CSUMTT12(14))&
        )*ADIV2
    !**********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM93(N4,ID,IK)=(CSUMTTMOM12(15,IK)&
        +GIOT*CSUMTTMOM12(16,IK)&
        -CSUMTTMOM12(13,IK)-GIOT*CSUMTTMOM12(14,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE94(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)&
        +(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))&
        +(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))&
        +(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM94(N4,ID,IK)=(CSUMTTMOM12(1,IK)-CSUMTTMOM12(2,IK)&
        +CSUMTTMOM12(3,IK)-CSUMTTMOM12(4,IK)&
        +(CSUMTTMOM12(9,IK)-CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)&
        -CSUMTTMOM12(12,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE95(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)&
        +(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))&
        -(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))&
        -(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM95(N4,ID,IK)=(CSUMTTMOM12(7,IK)-CSUMTTMOM12(8,IK)&
        +CSUMTTMOM12(5,IK)-CSUMTTMOM12(6,IK)+(CSUMTTMOM12(15,IK)&
        -CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)-CSUMTTMOM12(14,IK))&
        )*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE96(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)&
        -(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))&
        +(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))&
        -(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM96(N4,ID,IK)=(CSUMTTMOM12(1,IK)-CSUMTTMOM12(2,IK)&
        +CSUMTTMOM12(3,IK)-CSUMTTMOM12(4,IK)-(CSUMTTMOM12(9,IK)&
        -CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)-CSUMTTMOM12(12,IK))&
        )*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE97(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)&
        -(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))&
        -(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))&
        +(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM97(N4,ID,IK)=(CSUMTTMOM12(7,IK)-CSUMTTMOM12(8,IK)&
        +CSUMTTMOM12(5,IK)-CSUMTTMOM12(6,IK)-(CSUMTTMOM12(15,IK)&
        -CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)-CSUMTTMOM12(14,IK))&
        )*ADIV2
        enddo  
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     13
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     TT-13 OPERATORS CP=+,Pz=+,J=0
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE98(N4,ID)=(CSUMTT13(1)+CSUMTT13(2)+CSUMTT13(3)+CSUMTT13(4)+&
        CSUMTT13(5)+CSUMTT13(6)+CSUMTT13(7)+CSUMTT13(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM98(N4,ID,IK)=(CSUMTTMOM13(1,IK)+CSUMTTMOM13(2,IK)&
        +CSUMTTMOM13(3,IK)+CSUMTTMOM13(4,IK)+CSUMTTMOM13(5,IK)&
        +CSUMTTMOM13(6,IK)+CSUMTTMOM13(7,IK)&
        +CSUMTTMOM13(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE99(N4,ID)=(CSUMTT13(1)+CSUMTT13(2)+CSUMTT13(3)+CSUMTT13(4)&
        -(CSUMTT13(5)+CSUMTT13(6)+CSUMTT13(7)+CSUMTT13(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM99(N4,ID,IK)=(CSUMTTMOM13(1,IK)+CSUMTTMOM13(2,IK)&
        +CSUMTTMOM13(3,IK)+CSUMTTMOM13(4,IK)-(CSUMTTMOM13(5,IK)&
        +CSUMTTMOM13(6,IK)+CSUMTTMOM13(7,IK)&
        +CSUMTTMOM13(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE100(N4,ID)=(CSUMTT13(1)+GIOT*CSUMTT13(2)-CSUMTT13(3)&
        -GIOT*CSUMTT13(4)+(CSUMTT13(7)+GIOT*CSUMTT13(8)-CSUMTT13(5)&
        -GIOT*CSUMTT13(6)))*ADIV2
    !**********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM100(N4,ID,IK)=(CSUMTTMOM13(1,IK)+GIOT*CSUMTTMOM13(2,IK)&
        -CSUMTTMOM13(3,IK)-GIOT*CSUMTTMOM13(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE101(N4,ID)=(CSUMTT13(1)+GIOT*CSUMTT13(2)-CSUMTT13(3)&
        -GIOT*CSUMTT13(4)&
        -(CSUMTT13(7)+GIOT*CSUMTT13(8)-CSUMTT13(5)&
        -GIOT*CSUMTT13(6)))*ADIV2
    !**********************************************************************
    !     J=1, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM101(N4,ID,IK)=(CSUMTTMOM13(7,IK)+GIOT*CSUMTTMOM13(8,IK)&
        -CSUMTTMOM13(5,IK)-GIOT*CSUMTTMOM13(6,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE102(N4,ID)=(CSUMTT13(1)-CSUMTT13(2)+CSUMTT13(3)-CSUMTT13(4)+&
        CSUMTT13(6)-CSUMTT13(5)+CSUMTT13(8)-CSUMTT13(7))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM102(N4,ID,IK)=(CSUMTTMOM13(1,IK)-CSUMTTMOM13(2,IK)&
        +CSUMTTMOM13(3,IK)-CSUMTTMOM13(4,IK)+CSUMTTMOM13(6,IK)&
        -CSUMTTMOM13(5,IK)+CSUMTTMOM13(8,IK)&
        -CSUMTTMOM13(7,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE103(N4,ID)=(CSUMTT13(1)-CSUMTT13(2)+CSUMTT13(3)-CSUMTT13(4)&
        -(CSUMTT13(6)-CSUMTT13(5)+CSUMTT13(8)-CSUMTT13(7)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM103(N4,ID,IK)=(CSUMTTMOM13(1,IK)-CSUMTTMOM13(2,IK)&
        +CSUMTTMOM13(3,IK)-CSUMTTMOM13(4,IK)-(CSUMTTMOM13(6,IK)&
        -CSUMTTMOM13(5,IK)+CSUMTTMOM13(8,IK)&
        -CSUMTTMOM13(7,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     14 THISSSSSSSS!!!!!!!
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE104(N4,ID)=(CSUMTT14(1)+CSUMTT14(2)+CSUMTT14(3)+CSUMTT14(4)+&
        CSUMTT14(5)+CSUMTT14(6)+CSUMTT14(7)+CSUMTT14(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=0
    !**********************************************************************
        do IK=1,2
        ALINEMOM104(N4,ID,IK)=(CSUMTTMOM14(1,IK)+CSUMTTMOM14(2,IK)&
        +CSUMTTMOM14(3,IK)+CSUMTTMOM14(4,IK)+CSUMTTMOM14(5,IK)&
        +CSUMTTMOM14(6,IK)+CSUMTTMOM14(7,IK)&
        +CSUMTTMOM14(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE105(N4,ID)=(CSUMTT14(1)+CSUMTT14(2)+CSUMTT14(3)+CSUMTT14(4)-&
        (CSUMTT14(5)+CSUMTT14(6)+CSUMTT14(7)+CSUMTT14(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM105(N4,ID,IK)=(CSUMTTMOM14(1,IK)+CSUMTTMOM14(2,IK)&
        +CSUMTTMOM14(3,IK)+CSUMTTMOM14(4,IK)-(CSUMTTMOM14(5,IK)&
        +CSUMTTMOM14(6,IK)+CSUMTTMOM14(7,IK)&
        +CSUMTTMOM14(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************
        ALINE106(N4,ID)=(CSUMTT14(1)+GIOT*CSUMTT14(2)-CSUMTT14(3)&
        -GIOT*CSUMTT14(4)+(CSUMTT14(6)+GIOT*CSUMTT14(7)-CSUMTT14(8)&
        -GIOT*CSUMTT14(5)))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
        ALINEMOM106(N4,ID,IK)=(CSUMTTMOM14(1,IK)+GIOT*CSUMTTMOM14(2,IK)&
        -CSUMTTMOM14(3,IK)-GIOT*CSUMTTMOM14(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************
        ALINE107(N4,ID)=(CSUMTT14(1)+GIOT*CSUMTT14(2)-CSUMTT14(3)&
        -GIOT*CSUMTT14(4)-(CSUMTT14(6)+GIOT*CSUMTT14(7)-CSUMTT14(8)&
        -GIOT*CSUMTT14(5)))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM107(N4,ID,IK)=(CSUMTTMOM14(6,IK)+GIOT*CSUMTTMOM14(7,IK)&
        -CSUMTTMOM14(8,IK)-GIOT*CSUMTTMOM14(5,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE108(N4,ID)=(CSUMTT14(1)-CSUMTT14(2)+CSUMTT14(3)-CSUMTT14(4)+&
        (CSUMTT14(7)-CSUMTT14(6)+CSUMTT14(5)-CSUMTT14(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2,3,4
    !**********************************************************************
        do IK=1,2
            ALINEMOM108(N4,ID,IK)=(CSUMTTMOM14(1,IK)-CSUMTTMOM14(2,IK)&
        +CSUMTTMOM14(3,IK)-CSUMTTMOM14(4,IK)+(CSUMTTMOM14(7,IK)&
        -CSUMTTMOM14(6,IK)+CSUMTTMOM14(5,IK)-CSUMTTMOM14(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE109(N4,ID)=(CSUMTT14(1)-CSUMTT14(2)+CSUMTT14(3)-CSUMTT14(4)-&
        (CSUMTT14(7)-CSUMTT14(6)+CSUMTT14(5)-CSUMTT14(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2,3,4
    !**********************************************************************
        do IK=1,2 
        ALINEMOM109(N4,ID,IK)=(CSUMTTMOM14(1,IK)-CSUMTTMOM14(2,IK)&
        +CSUMTTMOM14(3,IK)-CSUMTTMOM14(4,IK)-(CSUMTTMOM14(7,IK)&
        -CSUMTTMOM14(6,IK)+CSUMTTMOM14(5,IK)-CSUMTTMOM14(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATOR 1
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE110(N4,ID)=(CSUMPLQ(1)+CSUMPLQ(2)+CSUMPLQ(3)+CSUMPLQ(4)&
        +CSUMPLQ(5)+CSUMPLQ(6)+CSUMPLQ(7)+CSUMPLQ(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM110(N4,ID,IK)=(CSUMPLQMOM(1,IK)+CSUMPLQMOM(2,IK)&
        +CSUMPLQMOM(3,IK)+CSUMPLQMOM(4,IK)+CSUMPLQMOM(5,IK)&
        +CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)+CSUMPLQMOM(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************      
        ALINE111(N4,ID)=(CSUMPLQ(1)+CSUMPLQ(2)+CSUMPLQ(3)+CSUMPLQ(4)&
        -CSUMPLQ(5)-CSUMPLQ(6)-CSUMPLQ(7)-CSUMPLQ(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM111(N4,ID,IK)=(CSUMPLQMOM(1,IK)+CSUMPLQMOM(2,IK)&
        +CSUMPLQMOM(3,IK)+CSUMPLQMOM(4,IK)-CSUMPLQMOM(5,IK)&
        -CSUMPLQMOM(6,IK)-CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************      
        ALINE112(N4,ID)=(CSUMPLQ(1)+GIOT*CSUMPLQ(2)-CSUMPLQ(3)&
        -GIOT*CSUMPLQ(4)&
        +(CSUMPLQ(6)+GIOT*CSUMPLQ(5)-CSUMPLQ(8)-GIOT*CSUMPLQ(7)))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************      
        do IK=1, 2
            ALINEMOM112(N4,ID,IK)=(CSUMPLQMOM(1,IK)+GIOT*CSUMPLQMOM(2,IK)&
        -CSUMPLQMOM(3,IK)-GIOT*CSUMPLQMOM(4,IK))*ADIV1
        enddo   
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************      
        ALINE113(N4,ID)=(CSUMPLQ(1)+GIOT*CSUMPLQ(2)-CSUMPLQ(3)&
        -GIOT*CSUMPLQ(4)&
        -(CSUMPLQ(6)+GIOT*CSUMPLQ(5)-CSUMPLQ(8)-GIOT*CSUMPLQ(7)))*ADIV2
    !**********************************************************************
    !     J=1, q=0
    !**********************************************************************      
        do IK=1, 2
            ALINEMOM113(N4,ID,IK)=(CSUMPLQMOM(6,IK)+GIOT*CSUMPLQMOM(5,IK)&
        -CSUMPLQMOM(8,IK)-GIOT*CSUMPLQMOM(7,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************      
        ALINE114(N4,ID)=(CSUMPLQ(1)-CSUMPLQ(2)+CSUMPLQ(3)-CSUMPLQ(4)&
        +CSUMPLQ(5)-CSUMPLQ(6)+CSUMPLQ(7)-CSUMPLQ(8))*ADIV2
    !**********************************************************************
    !     J=2, Pr=+, q=0
    !**********************************************************************      
        do IK=1, 2
            ALINEMOM114(N4,ID,IK)=(CSUMPLQMOM(1,IK)-CSUMPLQMOM(2,IK)&
        +CSUMPLQMOM(3,IK)-CSUMPLQMOM(4,IK)+CSUMPLQMOM(5,IK)&
        -CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************            
        ALINE115(N4,ID)=(CSUMPLQ(1)-CSUMPLQ(2)+CSUMPLQ(3)-CSUMPLQ(4)&
        -(CSUMPLQ(5)-CSUMPLQ(6)+CSUMPLQ(7)-CSUMPLQ(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pr=-, q=0
    !**********************************************************************            
        do IK=1, 2
            ALINEMOM115(N4,ID,IK)=(CSUMPLQMOM(1,IK)-CSUMPLQMOM(2,IK)&
        +CSUMPLQMOM(3,IK)-CSUMPLQMOM(4,IK)-(CSUMPLQMOM(5,IK)&
        -CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATOR 2
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !      
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************                  
        ALINE116(N4,ID)=(CSUMPLQ2(1)+CSUMPLQ2(2)+CSUMPLQ2(3)+CSUMPLQ2(4)&
        +CSUMPLQ2(5)+CSUMPLQ2(6)+CSUMPLQ2(7)+CSUMPLQ2(8))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM116(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+CSUMPLQMOM2(2,IK)&
        +CSUMPLQMOM2(3,IK)+CSUMPLQMOM2(4,IK)+CSUMPLQMOM2(5,IK)&
        +CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)+CSUMPLQMOM2(8,IK))*ADIV2
        enddo   
    !**********************************************************************
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************                  
        ALINE117(N4,ID)=(CSUMPLQ2(1)+CSUMPLQ2(2)+CSUMPLQ2(3)+CSUMPLQ2(4)&
        -(CSUMPLQ2(5)+CSUMPLQ2(6)+CSUMPLQ2(7)+CSUMPLQ2(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM117(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+CSUMPLQMOM2(2,IK)&
        +CSUMPLQMOM2(3,IK)+CSUMPLQMOM2(4,IK)-(CSUMPLQMOM2(5,IK)&
        +CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)+CSUMPLQMOM2(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************                        
        ALINE118(N4,ID)=(CSUMPLQ2(1)+GIOT*CSUMPLQ2(2)-CSUMPLQ2(3)&
        -GIOT*CSUMPLQ2(4)&
        +CSUMPLQ2(6)+GIOT*CSUMPLQ2(5)-CSUMPLQ2(8)-GIOT*CSUMPLQ2(7))*ADIV2
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM118(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+GIOT*CSUMPLQMOM2(2,IK)&
        -CSUMPLQMOM2(3,IK)-GIOT*CSUMPLQMOM2(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=-, q=0
    !**********************************************************************                            
        ALINE119(N4,ID)=(CSUMPLQ2(1)+GIOT*CSUMPLQ2(2)-CSUMPLQ2(3)&
        -GIOT*CSUMPLQ2(4)&
        -(CSUMPLQ2(6)+GIOT*CSUMPLQ2(5)-CSUMPLQ2(8)-GIOT*CSUMPLQ2(7)))&
        *ADIV2
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM119(N4,ID,IK)=(CSUMPLQMOM2(6,IK)+GIOT*CSUMPLQMOM2(5,IK)&
        -CSUMPLQMOM2(8,IK)-GIOT*CSUMPLQMOM2(7,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************                  
        ALINE120(N4,ID)=(CSUMPLQ2(1)-CSUMPLQ2(2)+CSUMPLQ2(3)-CSUMPLQ2(4)&
        +CSUMPLQ2(5)-CSUMPLQ2(6)+CSUMPLQ2(7)-CSUMPLQ2(8))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM120(N4,ID,IK)=(CSUMPLQMOM2(1,IK)-CSUMPLQMOM2(2,IK)&
        +CSUMPLQMOM2(3,IK)-CSUMPLQMOM2(4,IK)+CSUMPLQMOM2(5,IK)&
        -CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)-CSUMPLQMOM2(8,IK))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************                        
        ALINE121(N4,ID)=(CSUMPLQ2(1)-CSUMPLQ2(2)+CSUMPLQ2(3)-CSUMPLQ2(4)&
        -(CSUMPLQ2(5)-CSUMPLQ2(6)+CSUMPLQ2(7)-CSUMPLQ2(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM121(N4,ID,IK)=(CSUMPLQMOM2(1,IK)-CSUMPLQMOM2(2,IK)&
        +CSUMPLQMOM2(3,IK)-CSUMPLQMOM2(4,IK)-(CSUMPLQMOM2(5,IK)&
        -CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)-CSUMPLQMOM2(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! PLAQUETTE OPERATORS 3
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************                        
        ALINE122(N4,ID)=(CSUMPLQ3(1)+CSUMPLQ3(2)+CSUMPLQ3(3)&
        +CSUMPLQ3(4)+(CSUMPLQ3(5)+CSUMPLQ3(6)+CSUMPLQ3(7)&
        +CSUMPLQ3(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=+, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM122(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+CSUMPLQMOM3(2,IK)&
        +CSUMPLQMOM3(3,IK)+CSUMPLQMOM3(4,IK)&
        +(CSUMPLQMOM3(5,IK)+CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)&
        +CSUMPLQMOM3(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************                        
        ALINE123(N4,ID)=(CSUMPLQ3(1)+CSUMPLQ3(2)+CSUMPLQ3(3)&
        +CSUMPLQ3(4)-(CSUMPLQ3(5)+CSUMPLQ3(6)+CSUMPLQ3(7)&
        +CSUMPLQ3(8)))*ADIV2
    !**********************************************************************
    !     J=0, Pp=-, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM123(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+CSUMPLQMOM3(2,IK)&
        +CSUMPLQMOM3(3,IK)+CSUMPLQMOM3(4,IK)&
        -(CSUMPLQMOM3(5,IK)+CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)&
        +CSUMPLQMOM3(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************                            
        ALINE124(N4,ID)=(CSUMPLQ3(1)+GIOT*CSUMPLQ3(2)-CSUMPLQ3(3)&
        -GIOT*CSUMPLQ3(4))*ADIV1
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************                            
        do IK=1, 2
            ALINEMOM124(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+GIOT*CSUMPLQMOM3(2,IK)&
        -CSUMPLQMOM3(3,IK)-GIOT*CSUMPLQMOM3(4,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=1, Pr=+, q=0
    !**********************************************************************                            
        ALINE125(N4,ID)=(CSUMPLQ3(7)+GIOT*CSUMPLQ3(6)-CSUMPLQ3(5)&
        -GIOT*CSUMPLQ3(8))*ADIV1
    !**********************************************************************
    !     J=1, q=1,2
    !**********************************************************************                            
        do IK=1, 2
            ALINEMOM125(N4,ID,IK)=(CSUMPLQMOM3(7,IK)+GIOT*CSUMPLQMOM3(6,IK)&
        -CSUMPLQMOM3(5,IK)-GIOT*CSUMPLQMOM3(8,IK))*ADIV1
        enddo
    !**********************************************************************
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************                            
        ALINE126(N4,ID)=(CSUMPLQ3(1)-CSUMPLQ3(2)+CSUMPLQ3(3)&
        -CSUMPLQ3(4)+(CSUMPLQ3(5)-CSUMPLQ3(6)+CSUMPLQ3(7)&
        -CSUMPLQ3(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=+, q=1,2
    !**********************************************************************                            
        do IK=1, 2
            ALINEMOM126(N4,ID,IK)=(CSUMPLQMOM3(1,IK)-CSUMPLQMOM3(2,IK)&
        +CSUMPLQMOM3(3,IK)-CSUMPLQMOM3(4,IK)&
        +(CSUMPLQMOM3(5,IK)-CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)&
        -CSUMPLQMOM3(8,IK)))*ADIV2
        enddo
    !**********************************************************************
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************                            
        ALINE127(N4,ID)=(CSUMPLQ3(1)-CSUMPLQ3(2)+CSUMPLQ3(3)&
        -CSUMPLQ3(4)-(CSUMPLQ3(5)-CSUMPLQ3(6)+CSUMPLQ3(7)&
        -CSUMPLQ3(8)))*ADIV2
    !**********************************************************************
    !     J=2, Pp=-, q=1,2
    !**********************************************************************                            
        do IK=1, 2
            ALINEMOM127(N4,ID,IK)=(CSUMPLQMOM3(1,IK)-CSUMPLQMOM3(2,IK)&
        +CSUMPLQMOM3(3,IK)-CSUMPLQMOM3(4,IK)&
        -(CSUMPLQMOM3(5,IK)-CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)&
        -CSUMPLQMOM3(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 4
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************            
        ALINE128(N4,ID)=(CSUMPLQ4(1)+CSUMPLQ4(2)+CSUMPLQ4(3)&
        +CSUMPLQ4(4)+(CSUMPLQ4(5)+CSUMPLQ4(6)+CSUMPLQ4(7)&
        +CSUMPLQ4(8)))*ADIV2
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM128(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+CSUMPLQMOM4(2,IK)&
        +CSUMPLQMOM4(3,IK)+CSUMPLQMOM4(4,IK)&
        +(CSUMPLQMOM4(5,IK)+CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)&
        +CSUMPLQMOM4(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************                  
        ALINE129(N4,ID)=(CSUMPLQ4(1)+CSUMPLQ4(2)+CSUMPLQ4(3)&
        +CSUMPLQ4(4)-(CSUMPLQ4(5)+CSUMPLQ4(6)+CSUMPLQ4(7)&
        +CSUMPLQ4(8)))*ADIV2
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM129(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+CSUMPLQMOM4(2,IK)&
        +CSUMPLQMOM4(3,IK)+CSUMPLQMOM4(4,IK)&
        -(CSUMPLQMOM4(5,IK)+CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)&
        +CSUMPLQMOM4(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************                            
        ALINE130(N4,ID)=(CSUMPLQ4(1)+GIOT*CSUMPLQ4(2)-CSUMPLQ4(3)&
        -GIOT*CSUMPLQ4(4))*ADIV1
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                            
        do IK=1, 2
            ALINEMOM130(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+GIOT*CSUMPLQMOM4(2,IK)&
        -CSUMPLQMOM4(3,IK)-GIOT*CSUMPLQMOM4(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************                            
        ALINE131(N4,ID)=(CSUMPLQ4(7)+GIOT*CSUMPLQ4(6)-CSUMPLQ4(5)&
        -GIOT*CSUMPLQ4(8))*ADIV1
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                            
        do IK=1, 2
            ALINEMOM131(N4,ID,IK)=(CSUMPLQMOM4(7,IK)+GIOT*CSUMPLQMOM4(6,IK)&
        -CSUMPLQMOM4(5,IK)-GIOT*CSUMPLQMOM4(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************            
        ALINE132(N4,ID)=(CSUMPLQ4(1)-CSUMPLQ4(2)+CSUMPLQ4(3)&
        -CSUMPLQ4(4)+(CSUMPLQ4(5)-CSUMPLQ4(6)+CSUMPLQ4(7)&
        -CSUMPLQ4(8)))*ADIV2
    !**********************************************************************      
    !     J=2, Pp=+, q=1,2
    !**********************************************************************            
        do IK=1, 2
            ALINEMOM132(N4,ID,IK)=(CSUMPLQMOM4(1,IK)-CSUMPLQMOM4(2,IK)&
        +CSUMPLQMOM4(3,IK)-CSUMPLQMOM4(4,IK)&
        +(CSUMPLQMOM4(5,IK)-CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)&
        -CSUMPLQMOM4(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************                  
        ALINE133(N4,ID)=(CSUMPLQ4(1)-CSUMPLQ4(2)+CSUMPLQ4(3)&
        -CSUMPLQ4(4)-(CSUMPLQ4(5)-CSUMPLQ4(6)+CSUMPLQ4(7)&
        -CSUMPLQ4(8)))*ADIV2
    !**********************************************************************      
    !     J=2, Pp=-, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM133(N4,ID,IK)=(CSUMPLQMOM4(1,IK)-CSUMPLQMOM4(2,IK)&
        +CSUMPLQMOM4(3,IK)-CSUMPLQMOM4(4,IK)&
        -(CSUMPLQMOM4(5,IK)-CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)&
        -CSUMPLQMOM4(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! PLAQUETTE OPERATOR 5  
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************                        
        ALINE134(N4,ID)=(CSUMPLQ5(1)+CSUMPLQ5(2)+CSUMPLQ5(3)&
        +CSUMPLQ5(4)+(CSUMPLQ5(5)+CSUMPLQ5(6)+CSUMPLQ5(7)&
        +CSUMPLQ5(8)))*ADIV2
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************                            
        do IK=1, 2
            ALINEMOM134(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+CSUMPLQMOM5(2,IK)&
        +CSUMPLQMOM5(3,IK)+CSUMPLQMOM5(4,IK)&
        +(CSUMPLQMOM5(5,IK)+CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)&
        +CSUMPLQMOM5(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************                            
        ALINE135(N4,ID)=(CSUMPLQ5(1)+CSUMPLQ5(2)+CSUMPLQ5(3)&
        +CSUMPLQ5(4)-(CSUMPLQ5(5)+CSUMPLQ5(6)+CSUMPLQ5(7)&
        +CSUMPLQ5(8)))*ADIV2
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************      
        do IK=1, 2
            ALINEMOM135(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+CSUMPLQMOM5(2,IK)&
        +CSUMPLQMOM5(3,IK)+CSUMPLQMOM5(4,IK)&
        -(CSUMPLQMOM5(5,IK)+CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)&
        +CSUMPLQMOM5(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************            
        ALINE136(N4,ID)=(CSUMPLQ5(1)+GIOT*CSUMPLQ5(2)-CSUMPLQ5(3)&
        -GIOT*CSUMPLQ5(4))*ADIV1
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM136(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+GIOT*CSUMPLQMOM5(2,IK)&
        -CSUMPLQMOM5(3,IK)-GIOT*CSUMPLQMOM5(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE137(N4,ID)=(CSUMPLQ5(7)+GIOT*CSUMPLQ5(6)-CSUMPLQ5(5)&
        -GIOT*CSUMPLQ5(8))*ADIV1
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM137(N4,ID,IK)=(CSUMPLQMOM5(7,IK)+GIOT*CSUMPLQMOM5(6,IK)&
        -CSUMPLQMOM5(5,IK)-GIOT*CSUMPLQMOM5(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************                        
        ALINE138(N4,ID)=(CSUMPLQ5(1)-CSUMPLQ5(2)+CSUMPLQ5(3)&
        -CSUMPLQ5(4)+(CSUMPLQ5(5)-CSUMPLQ5(6)+CSUMPLQ5(7)&
        -CSUMPLQ5(8)))*ADIV2
    !**********************************************************************      
    !     J=2, Pp=+, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM138(N4,ID,IK)=(CSUMPLQMOM5(1,IK)-CSUMPLQMOM5(2,IK)&
        +CSUMPLQMOM5(3,IK)-CSUMPLQMOM5(4,IK)&
        +(CSUMPLQMOM5(5,IK)-CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)&
        -CSUMPLQMOM5(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************                            
        ALINE139(N4,ID)=(CSUMPLQ5(1)-CSUMPLQ5(2)+CSUMPLQ5(3)&
        -CSUMPLQ5(4)-(CSUMPLQ5(5)-CSUMPLQ5(6)+CSUMPLQ5(7)&
        -CSUMPLQ5(8)))*ADIV2
    !**********************************************************************      
    !     J=2, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM139(N4,ID,IK)=(CSUMPLQMOM5(1,IK)-CSUMPLQMOM5(2,IK)&
        +CSUMPLQMOM5(3,IK)-CSUMPLQMOM5(4,IK)&
        -(CSUMPLQMOM5(5,IK)-CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)&
        -CSUMPLQMOM5(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 6
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE140(N4,ID)=(CSUMPLQ6(1)+CSUMPLQ6(2)+CSUMPLQ6(3)&
        +CSUMPLQ6(4)+(CSUMPLQ6(5)+CSUMPLQ6(6)+CSUMPLQ6(7)&
        +CSUMPLQ6(8)))*ADIV2
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM140(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+CSUMPLQMOM6(2,IK)&
        +CSUMPLQMOM6(3,IK)+CSUMPLQMOM6(4,IK)&
        +(CSUMPLQMOM6(5,IK)+CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)&
        +CSUMPLQMOM6(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************      
        ALINE141(N4,ID)=(CSUMPLQ6(1)+CSUMPLQ6(2)+CSUMPLQ6(3)&
        +CSUMPLQ6(4)-(CSUMPLQ6(5)+CSUMPLQ6(6)+CSUMPLQ6(7)&
            +CSUMPLQ6(8)))*ADIV2
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************            
        do IK=1, 2
            ALINEMOM141(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+CSUMPLQMOM6(2,IK)&
        +CSUMPLQMOM6(3,IK)+CSUMPLQMOM6(4,IK)&
        -(CSUMPLQMOM6(5,IK)+CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)&
        +CSUMPLQMOM6(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE142(N4,ID)=(CSUMPLQ6(1)+GIOT*CSUMPLQ6(2)-CSUMPLQ6(3)&
        -GIOT*CSUMPLQ6(4)+(CSUMPLQ6(8)+GIOT*CSUMPLQ6(7)-CSUMPLQ6(6)&
        -GIOT*CSUMPLQ6(5)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM142(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+GIOT*CSUMPLQMOM6(2,IK)&
        -CSUMPLQMOM6(3,IK)-GIOT*CSUMPLQMOM6(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE143(N4,ID)=(CSUMPLQ6(1)+GIOT*CSUMPLQ6(2)-CSUMPLQ6(3)&
        -GIOT*CSUMPLQ6(4)-(CSUMPLQ6(8)+GIOT*CSUMPLQ6(7)-CSUMPLQ6(6)&
        -GIOT*CSUMPLQ6(5)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM143(N4,ID,IK)=(CSUMPLQMOM6(8,IK)+GIOT*CSUMPLQMOM6(7,IK)&
        -CSUMPLQMOM6(6,IK)-GIOT*CSUMPLQMOM6(5,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************                        
        ALINE144(N4,ID)=(CSUMPLQ6(1)-CSUMPLQ6(2)+CSUMPLQ6(3)&
        -CSUMPLQ6(4)+(CSUMPLQ6(5)-CSUMPLQ6(6)+CSUMPLQ6(7)&
        -CSUMPLQ6(8)))*ADIV2
    !**********************************************************************      
    !     J=2, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM144(N4,ID,IK)=(CSUMPLQMOM6(1,IK)-CSUMPLQMOM6(2,IK)&
        +CSUMPLQMOM6(3,IK)-CSUMPLQMOM6(4,IK)&
        +(CSUMPLQMOM6(5,IK)-CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)&
        -CSUMPLQMOM6(8,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE145(N4,ID)=(CSUMPLQ6(1)-CSUMPLQ6(2)+CSUMPLQ6(3)&
        -CSUMPLQ6(4)-(CSUMPLQ6(5)-CSUMPLQ6(6)+CSUMPLQ6(7)&
        -CSUMPLQ6(8)))*ADIV2
    !**********************************************************************      
    !     J=2, Pp=-, q=1,2
    !**********************************************************************      
        do IK=1, 2
            ALINEMOM145(N4,ID,IK)=(CSUMPLQMOM6(1,IK)-CSUMPLQMOM6(2,IK)&
        +CSUMPLQMOM6(3,IK)-CSUMPLQMOM6(4,IK)&
        -(CSUMPLQMOM6(5,IK)-CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)&
        -CSUMPLQMOM6(8,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 7
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE146(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)&
        +CSUMPLQ7(4)&
        +(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))&
        +(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))&
        +(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM146(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+CSUMPLQMOM7(2,IK)&
        +CSUMPLQMOM7(3,IK)+CSUMPLQMOM7(4,IK)&
        +(CSUMPLQMOM7(9,IK)+CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)&
        +CSUMPLQMOM7(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE147(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)&
        +CSUMPLQ7(4)&
        -(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))&
        +(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))&
        -(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM147(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+CSUMPLQMOM7(6,IK)&
        +CSUMPLQMOM7(7,IK)+CSUMPLQMOM7(8,IK)&
        +(CSUMPLQMOM7(13,IK)+CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)&
        +CSUMPLQMOM7(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE148(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)&
        +CSUMPLQ7(4)&
        +(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))&
        -(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))&
        -(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM148(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+CSUMPLQMOM7(2,IK)&
        +CSUMPLQMOM7(3,IK)+CSUMPLQMOM7(4,IK)&
        -(CSUMPLQMOM7(9,IK)+CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)&
        +CSUMPLQMOM7(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE149(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)&
        +CSUMPLQ7(4)&
        -(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))&
        -(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))&
        +(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM149(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+CSUMPLQMOM7(6,IK)&
        +CSUMPLQMOM7(7,IK)+CSUMPLQMOM7(8,IK)&
        -(CSUMPLQMOM7(13,IK)+CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)&
        +CSUMPLQMOM7(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE150(N4,ID)=(CSUMPLQ7(1)+GIOT*CSUMPLQ7(2)-CSUMPLQ7(3)&
        -GIOT*CSUMPLQ7(4)+(CSUMPLQ7(5)+GIOT*CSUMPLQ7(6)-CSUMPLQ7(7)&
        -GIOT*CSUMPLQ7(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM150(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+GIOT*CSUMPLQMOM7(2,IK)&
        -CSUMPLQMOM7(3,IK)-GIOT*CSUMPLQMOM7(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE151(N4,ID)=(CSUMPLQ7(1)+GIOT*CSUMPLQ7(2)-CSUMPLQ7(3)&
        -GIOT*CSUMPLQ7(4)-(CSUMPLQ7(5)+GIOT*CSUMPLQ7(6)-CSUMPLQ7(7)&
        -GIOT*CSUMPLQ7(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM151(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+GIOT*CSUMPLQMOM7(6,IK)&
        -CSUMPLQMOM7(7,IK)-GIOT*CSUMPLQMOM7(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE152(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)&
        -CSUMPLQ7(4)&
        +(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))&
        +(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))&
        +(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q
    !**********************************************************************
        do IK=1, 2
            ALINEMOM152(N4,ID,IK)=(CSUMPLQMOM7(1,IK)-CSUMPLQMOM7(2,IK)&
        +CSUMPLQMOM7(3,IK)-CSUMPLQMOM7(4,IK)&
        +(CSUMPLQMOM7(9,IK)-CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)&
        -CSUMPLQMOM7(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE153(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)&
        -CSUMPLQ7(4)&
        -(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))&
        +(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))&
        -(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM153(N4,ID,IK)=(CSUMPLQMOM7(5,IK)-CSUMPLQMOM7(6,IK)&
        +CSUMPLQMOM7(7,IK)-CSUMPLQMOM7(8,IK)&
        +(CSUMPLQMOM7(13,IK)-CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)&
        -CSUMPLQMOM7(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE154(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)&
        -CSUMPLQ7(4)&
        +(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))&
        -(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))&
        -(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM154(N4,ID,IK)=(CSUMPLQMOM7(1,IK)-CSUMPLQMOM7(2,IK)&
        +CSUMPLQMOM7(3,IK)-CSUMPLQMOM7(4,IK)&
        -(CSUMPLQMOM7(9,IK)-CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)&
        -CSUMPLQMOM7(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE155(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)&
        -CSUMPLQ7(4)&
        -(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))&
        -(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))&
        +(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM155(N4,ID,IK)=(CSUMPLQMOM7(5,IK)-CSUMPLQMOM7(6,IK)&
        +CSUMPLQMOM7(7,IK)-CSUMPLQMOM7(8,IK)&
        -(CSUMPLQMOM7(13,IK)-CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)&
        -CSUMPLQMOM7(16,IK)))*ADIV2
        enddo      
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 8
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE156(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)&
        +CSUMPLQ8(4)&
        +(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))&
        +(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))&
        +(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM156(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+CSUMPLQMOM8(2,IK)&
        +CSUMPLQMOM8(3,IK)+CSUMPLQMOM8(4,IK)&
        +(CSUMPLQMOM8(9,IK)+CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)&
        +CSUMPLQMOM8(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE157(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)&
        +CSUMPLQ8(4)&
        -(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))&
        +(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))&
        -(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM157(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+CSUMPLQMOM8(6,IK)&
        +CSUMPLQMOM8(7,IK)+CSUMPLQMOM8(8,IK)&
        +(CSUMPLQMOM8(13,IK)+CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)&
        +CSUMPLQMOM8(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE158(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)&
        +CSUMPLQ8(4)&
        +(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))&
        -(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))&
        -(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM158(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+CSUMPLQMOM8(2,IK)&
        +CSUMPLQMOM8(3,IK)+CSUMPLQMOM8(4,IK)&
        -(CSUMPLQMOM8(9,IK)+CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)&
        +CSUMPLQMOM8(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE159(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)&
        +CSUMPLQ8(4)&
        -(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))&
        -(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))&
        +(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM159(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+CSUMPLQMOM8(6,IK)&
        +CSUMPLQMOM8(7,IK)+CSUMPLQMOM8(8,IK)&
        -(CSUMPLQMOM8(13,IK)+CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)&
        +CSUMPLQMOM8(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE160(N4,ID)=(CSUMPLQ8(1)+GIOT*CSUMPLQ8(2)-CSUMPLQ8(3)&
        -GIOT*CSUMPLQ8(4)+(CSUMPLQ8(5)+GIOT*CSUMPLQ8(6)-CSUMPLQ8(7)&
        -GIOT*CSUMPLQ8(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM160(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+GIOT*CSUMPLQMOM8(2,IK)&
        -CSUMPLQMOM8(3,IK)-GIOT*CSUMPLQMOM8(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE161(N4,ID)=(CSUMPLQ8(1)+GIOT*CSUMPLQ8(2)-CSUMPLQ8(3)&
        -GIOT*CSUMPLQ8(4)-(CSUMPLQ8(5)+GIOT*CSUMPLQ8(6)-CSUMPLQ8(7)&
        -GIOT*CSUMPLQ8(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM161(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+GIOT*CSUMPLQMOM8(6,IK)&
        -CSUMPLQMOM8(7,IK)-GIOT*CSUMPLQMOM8(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE162(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)&
        -CSUMPLQ8(4)&
        +(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))&
        +(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))&
        +(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM162(N4,ID,IK)=(CSUMPLQMOM8(1,IK)-CSUMPLQMOM8(2,IK)&
        +CSUMPLQMOM8(3,IK)-CSUMPLQMOM8(4,IK)&
        +(CSUMPLQMOM8(9,IK)-CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)&
        -CSUMPLQMOM8(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE163(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)&
        -CSUMPLQ8(4)&
        -(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))&
        +(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))&
        -(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM163(N4,ID,IK)=(CSUMPLQMOM8(5,IK)-CSUMPLQMOM8(6,IK)&
        +CSUMPLQMOM8(7,IK)-CSUMPLQMOM8(8,IK)&
        +(CSUMPLQMOM8(13,IK)-CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)&
        -CSUMPLQMOM8(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE164(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)&
        -CSUMPLQ8(4)&
        +(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))&
        -(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))&
        -(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM164(N4,ID,IK)=(CSUMPLQMOM8(1,IK)-CSUMPLQMOM8(2,IK)&
        +CSUMPLQMOM8(3,IK)-CSUMPLQMOM8(4,IK)&
        -(CSUMPLQMOM8(9,IK)-CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)&
        -CSUMPLQMOM8(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE165(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)&
        -CSUMPLQ8(4)&
        -(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))&
        -(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))&
        +(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM165(N4,ID,IK)=(CSUMPLQMOM8(5,IK)-CSUMPLQMOM8(6,IK)&
        +CSUMPLQMOM8(7,IK)-CSUMPLQMOM8(8,IK)&
        -(CSUMPLQMOM8(13,IK)-CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)&
        -CSUMPLQMOM8(16,IK)))*ADIV2
        enddo      
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 9
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE166(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)&
        +CSUMPLQ9(4)&
        +(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))&
        +(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))&
        +(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM166(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+CSUMPLQMOM9(2,IK)&
        +CSUMPLQMOM9(3,IK)+CSUMPLQMOM9(4,IK)&
        +(CSUMPLQMOM9(9,IK)+CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)&
        +CSUMPLQMOM9(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE167(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)&
        +CSUMPLQ9(4)&
        -(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))&
        +(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))&
        -(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM167(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+CSUMPLQMOM9(6,IK)&
        +CSUMPLQMOM9(7,IK)+CSUMPLQMOM9(8,IK)&
        +(CSUMPLQMOM9(13,IK)+CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)&
        +CSUMPLQMOM9(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE168(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)&
        +CSUMPLQ9(4)&
        +(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))&
        -(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))&
        -(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM168(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+CSUMPLQMOM9(2,IK)&
        +CSUMPLQMOM9(3,IK)+CSUMPLQMOM9(4,IK)&
        -(CSUMPLQMOM9(9,IK)+CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)&
        +CSUMPLQMOM9(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE169(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)&
        +CSUMPLQ9(4)&
        -(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))&
        -(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))&
        +(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM169(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+CSUMPLQMOM9(6,IK)&
        +CSUMPLQMOM9(7,IK)+CSUMPLQMOM9(8,IK)&
        -(CSUMPLQMOM9(13,IK)+CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)&
        +CSUMPLQMOM9(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE170(N4,ID)=(CSUMPLQ9(1)+GIOT*CSUMPLQ9(2)-CSUMPLQ9(3)&
        -GIOT*CSUMPLQ9(4)+(CSUMPLQ9(5)+GIOT*CSUMPLQ9(6)-CSUMPLQ9(7)&
        -GIOT*CSUMPLQ9(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM170(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+GIOT*CSUMPLQMOM9(2,IK)&
        -CSUMPLQMOM9(3,IK)-GIOT*CSUMPLQMOM9(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE171(N4,ID)=(CSUMPLQ9(1)+GIOT*CSUMPLQ9(2)-CSUMPLQ9(3)&
        -GIOT*CSUMPLQ9(4)-(CSUMPLQ9(5)+GIOT*CSUMPLQ9(6)-CSUMPLQ9(7)&
        -GIOT*CSUMPLQ9(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM171(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+GIOT*CSUMPLQMOM9(6,IK)&
        -CSUMPLQMOM9(7,IK)-GIOT*CSUMPLQMOM9(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE172(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)&
        -CSUMPLQ9(4)&
        +(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))&
        +(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))&
        +(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM172(N4,ID,IK)=(CSUMPLQMOM9(1,IK)-CSUMPLQMOM9(2,IK)&
        +CSUMPLQMOM9(3,IK)-CSUMPLQMOM9(4,IK)&
        +(CSUMPLQMOM9(9,IK)-CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)&
        -CSUMPLQMOM9(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE173(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)&
        -CSUMPLQ9(4)&
        -(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))&
        +(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))&
        -(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM173(N4,ID,IK)=(CSUMPLQMOM9(5,IK)-CSUMPLQMOM9(6,IK)&
        +CSUMPLQMOM9(7,IK)-CSUMPLQMOM9(8,IK)&
        +(CSUMPLQMOM9(13,IK)-CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)&
        -CSUMPLQMOM9(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE174(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)&
        -CSUMPLQ9(4)&
        +(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))&
        -(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))&
        -(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM174(N4,ID,IK)=(CSUMPLQMOM9(1,IK)-CSUMPLQMOM9(2,IK)&
        +CSUMPLQMOM9(3,IK)-CSUMPLQMOM9(4,IK)&
        -(CSUMPLQMOM9(9,IK)-CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)&
        -CSUMPLQMOM9(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE175(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)&
        -CSUMPLQ9(4)&
        -(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))&
        -(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))&
        +(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM175(N4,ID,IK)=(CSUMPLQMOM9(5,IK)-CSUMPLQMOM9(6,IK)&
        +CSUMPLQMOM9(7,IK)-CSUMPLQMOM9(8,IK)&
        -(CSUMPLQMOM9(13,IK)-CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)&
        -CSUMPLQMOM9(16,IK)))*ADIV2
        enddo      
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 10
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE176(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)&
        +CSUMPLQ10(4)&
        +(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))&
        +(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))&
        +(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM176(N4,ID,IK)=(CSUMPLQMOM10(1,IK)+CSUMPLQMOM10(2,IK)&
        +CSUMPLQMOM10(3,IK)+CSUMPLQMOM10(4,IK)&
        +(CSUMPLQMOM10(9,IK)+CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)&
        +CSUMPLQMOM10(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE177(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)&
        +CSUMPLQ10(4)&
        -(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))&
        +(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))&
        -(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM177(N4,ID,IK)=(CSUMPLQMOM10(5,IK)+CSUMPLQMOM10(6,IK)&
        +CSUMPLQMOM10(7,IK)+CSUMPLQMOM10(8,IK)&
        +(CSUMPLQMOM10(13,IK)+CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)&
        +CSUMPLQMOM10(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE178(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)&
        +CSUMPLQ10(4)&
        +(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))&
        -(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))&
        -(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM178(N4,ID,IK)=(CSUMPLQMOM10(1,IK)+CSUMPLQMOM10(2,IK)&
        +CSUMPLQMOM10(3,IK)+CSUMPLQMOM10(4,IK)&
        -(CSUMPLQMOM10(9,IK)+CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)&
        +CSUMPLQMOM10(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE179(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)&
        +CSUMPLQ10(4)&
        -(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))&
        -(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))&
        +(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM179(N4,ID,IK)=(CSUMPLQMOM10(5,IK)+CSUMPLQMOM10(6,IK)&
        +CSUMPLQMOM10(7,IK)+CSUMPLQMOM10(8,IK)&
        -(CSUMPLQMOM10(13,IK)+CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)&
        +CSUMPLQMOM10(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE180(N4,ID)=(CSUMPLQ10(1)+GIOT*CSUMPLQ10(2)-CSUMPLQ10(3)&
        -GIOT*CSUMPLQ10(4)+(CSUMPLQ10(5)+GIOT*CSUMPLQ10(6)-CSUMPLQ10(7)&
        -GIOT*CSUMPLQ10(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM180(N4,ID,IK)=(CSUMPLQMOM10(1,IK)&
        +GIOT*CSUMPLQMOM10(2,IK)&
        -CSUMPLQMOM10(3,IK)-GIOT*CSUMPLQMOM10(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE181(N4,ID)=(CSUMPLQ10(1)+GIOT*CSUMPLQ10(2)-CSUMPLQ10(3)&
        -GIOT*CSUMPLQ10(4)-(CSUMPLQ10(5)+GIOT*CSUMPLQ10(6)-CSUMPLQ10(7)&
        -GIOT*CSUMPLQ10(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM181(N4,ID,IK)=(CSUMPLQMOM10(5,IK)&
        +GIOT*CSUMPLQMOM10(6,IK)&
        -CSUMPLQMOM10(7,IK)-GIOT*CSUMPLQMOM10(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE182(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)&
        -CSUMPLQ10(4)&
        +(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))&
        +(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))&
        +(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM182(N4,ID,IK)=(CSUMPLQMOM10(1,IK)-CSUMPLQMOM10(2,IK)&
        +CSUMPLQMOM10(3,IK)-CSUMPLQMOM10(4,IK)&
        +(CSUMPLQMOM10(9,IK)-CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)&
        -CSUMPLQMOM10(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE183(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)&
        -CSUMPLQ10(4)&
        -(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))&
        +(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))&
        -(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM183(N4,ID,IK)=(CSUMPLQMOM10(5,IK)-CSUMPLQMOM10(6,IK)&
        +CSUMPLQMOM10(7,IK)-CSUMPLQMOM10(8,IK)&
        +(CSUMPLQMOM10(13,IK)-CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)&
        -CSUMPLQMOM10(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE184(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)&
        -CSUMPLQ10(4)&
        +(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))&
        -(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))&
        -(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM184(N4,ID,IK)=(CSUMPLQMOM10(1,IK)-CSUMPLQMOM10(2,IK)&
        +CSUMPLQMOM10(3,IK)-CSUMPLQMOM10(4,IK)&
        -(CSUMPLQMOM10(9,IK)-CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)&
        -CSUMPLQMOM10(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE185(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)&
        -CSUMPLQ10(4)&
        -(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))&
        -(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))&
        +(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM185(N4,ID,IK)=(CSUMPLQMOM10(5,IK)-CSUMPLQMOM10(6,IK)&
        +CSUMPLQMOM10(7,IK)-CSUMPLQMOM10(8,IK)&
        -(CSUMPLQMOM10(13,IK)-CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)&
        -CSUMPLQMOM10(16,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 11
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE186(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)&
        +CSUMPLQ11(4)&
        +(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))&
        +(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))&
        +(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM186(N4,ID,IK)=(CSUMPLQMOM11(1,IK)+CSUMPLQMOM11(2,IK)&
        +CSUMPLQMOM11(3,IK)+CSUMPLQMOM11(4,IK)&
        +(CSUMPLQMOM11(9,IK)+CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)&
        +CSUMPLQMOM11(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE187(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)&
        +CSUMPLQ11(4)&
        -(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))&
        +(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))&
        -(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM187(N4,ID,IK)=(CSUMPLQMOM11(5,IK)+CSUMPLQMOM11(6,IK)&
        +CSUMPLQMOM11(7,IK)+CSUMPLQMOM11(8,IK)&
        +(CSUMPLQMOM11(13,IK)+CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)&
        +CSUMPLQMOM11(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE188(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)&
        +CSUMPLQ11(4)&
        +(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))&
        -(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))&
        -(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM188(N4,ID,IK)=(CSUMPLQMOM11(1,IK)+CSUMPLQMOM11(2,IK)&
        +CSUMPLQMOM11(3,IK)+CSUMPLQMOM11(4,IK)&
        -(CSUMPLQMOM11(9,IK)+CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)&
        +CSUMPLQMOM11(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE189(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)&
        +CSUMPLQ11(4)&
        -(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))&
        -(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))&
        +(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM189(N4,ID,IK)=(CSUMPLQMOM11(5,IK)+CSUMPLQMOM11(6,IK)&
        +CSUMPLQMOM11(7,IK)+CSUMPLQMOM11(8,IK)&
        -(CSUMPLQMOM11(13,IK)+CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)&
        +CSUMPLQMOM11(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE190(N4,ID)=(CSUMPLQ11(1)+GIOT*CSUMPLQ11(2)-CSUMPLQ11(3)&
        -GIOT*CSUMPLQ11(4)+(CSUMPLQ11(5)+GIOT*CSUMPLQ11(6)-CSUMPLQ11(7)&
        -GIOT*CSUMPLQ11(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM190(N4,ID,IK)=(CSUMPLQMOM11(1,IK)&
        +GIOT*CSUMPLQMOM11(2,IK)&
        -CSUMPLQMOM11(3,IK)-GIOT*CSUMPLQMOM11(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE191(N4,ID)=(CSUMPLQ11(1)+GIOT*CSUMPLQ11(2)-CSUMPLQ11(3)&
        -GIOT*CSUMPLQ11(4)-(CSUMPLQ11(5)+GIOT*CSUMPLQ11(6)-CSUMPLQ11(7)&
        -GIOT*CSUMPLQ11(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM191(N4,ID,IK)=(CSUMPLQMOM11(5,IK)&
        +GIOT*CSUMPLQMOM11(6,IK)&
        -CSUMPLQMOM11(7,IK)-GIOT*CSUMPLQMOM11(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE192(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)&
        -CSUMPLQ11(4)&
        +(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))&
        +(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))&
        +(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM192(N4,ID,IK)=(CSUMPLQMOM11(1,IK)-CSUMPLQMOM11(2,IK)&
        +CSUMPLQMOM11(3,IK)-CSUMPLQMOM11(4,IK)&
        +(CSUMPLQMOM11(9,IK)-CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)&
        -CSUMPLQMOM11(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE193(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)&
        -CSUMPLQ11(4)&
        -(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))&
        +(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))&
        -(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM193(N4,ID,IK)=(CSUMPLQMOM11(5,IK)-CSUMPLQMOM11(6,IK)&
        +CSUMPLQMOM11(7,IK)-CSUMPLQMOM11(8,IK)&
        +(CSUMPLQMOM11(13,IK)-CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)&
        -CSUMPLQMOM11(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE194(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)&
        -CSUMPLQ11(4)&
        +(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))&
        -(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))&
        -(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM194(N4,ID,IK)=(CSUMPLQMOM11(1,IK)-CSUMPLQMOM11(2,IK)&
        +CSUMPLQMOM11(3,IK)-CSUMPLQMOM11(4,IK)&
        -(CSUMPLQMOM11(9,IK)-CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)&
        -CSUMPLQMOM11(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE195(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)&
        -CSUMPLQ11(4)&
        -(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))&
        -(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))&
        +(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM195(N4,ID,IK)=(CSUMPLQMOM11(5,IK)-CSUMPLQMOM11(6,IK)&
        +CSUMPLQMOM11(7,IK)-CSUMPLQMOM11(8,IK)&
        -(CSUMPLQMOM11(13,IK)-CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)&
        -CSUMPLQMOM11(16,IK)))*ADIV2
        enddo      
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 12
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE196(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)&
        +CSUMPLQ12(4)&
        +(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))&
        +(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))&
        +(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM196(N4,ID,IK)=(CSUMPLQMOM12(1,IK)+CSUMPLQMOM12(2,IK)&
        +CSUMPLQMOM12(3,IK)+CSUMPLQMOM12(4,IK)&
        +(CSUMPLQMOM12(9,IK)+CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)&
        +CSUMPLQMOM12(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE197(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)&
        +CSUMPLQ12(4)&
        -(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))&
        +(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))&
        -(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM197(N4,ID,IK)=(CSUMPLQMOM12(5,IK)+CSUMPLQMOM12(6,IK)&
        +CSUMPLQMOM12(7,IK)+CSUMPLQMOM12(8,IK)&
        +(CSUMPLQMOM12(13,IK)+CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)&
        +CSUMPLQMOM12(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE198(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)&
        +CSUMPLQ12(4)&
        +(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))&
        -(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))&
        -(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM198(N4,ID,IK)=(CSUMPLQMOM12(1,IK)+CSUMPLQMOM12(2,IK)&
        +CSUMPLQMOM12(3,IK)+CSUMPLQMOM12(4,IK)&
        -(CSUMPLQMOM12(9,IK)+CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)&
        +CSUMPLQMOM12(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE199(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)&
        +CSUMPLQ12(4)&
        -(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))&
        -(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))&
        +(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM199(N4,ID,IK)=(CSUMPLQMOM12(5,IK)+CSUMPLQMOM12(6,IK)&
        +CSUMPLQMOM12(7,IK)+CSUMPLQMOM12(8,IK)&
        -(CSUMPLQMOM12(13,IK)+CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)&
        +CSUMPLQMOM12(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE200(N4,ID)=(CSUMPLQ12(1)+GIOT*CSUMPLQ12(2)-CSUMPLQ12(3)&
        -GIOT*CSUMPLQ12(4)+(CSUMPLQ12(5)+GIOT*CSUMPLQ12(6)-CSUMPLQ12(7)&
        -GIOT*CSUMPLQ12(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM200(N4,ID,IK)=(CSUMPLQMOM12(1,IK)&
        +GIOT*CSUMPLQMOM12(2,IK)&
        -CSUMPLQMOM12(3,IK)-GIOT*CSUMPLQMOM12(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE201(N4,ID)=(CSUMPLQ12(1)+GIOT*CSUMPLQ12(2)-CSUMPLQ12(3)&
        -GIOT*CSUMPLQ12(4)-(CSUMPLQ12(5)+GIOT*CSUMPLQ12(6)-CSUMPLQ12(7)&
        -GIOT*CSUMPLQ12(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM201(N4,ID,IK)=(CSUMPLQMOM12(5,IK)&
        +GIOT*CSUMPLQMOM12(6,IK)&
        -CSUMPLQMOM12(7,IK)-GIOT*CSUMPLQMOM12(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE202(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)&
        -CSUMPLQ12(4)&
        +(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))&
        +(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))&
        +(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM202(N4,ID,IK)=(CSUMPLQMOM12(1,IK)-CSUMPLQMOM12(2,IK)&
        +CSUMPLQMOM12(3,IK)-CSUMPLQMOM12(4,IK)&
        +(CSUMPLQMOM12(9,IK)-CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)&
        -CSUMPLQMOM12(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE203(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)&
        -CSUMPLQ12(4)&
        -(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))&
        +(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))&
        -(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM203(N4,ID,IK)=(CSUMPLQMOM12(5,IK)-CSUMPLQMOM12(6,IK)&
        +CSUMPLQMOM12(7,IK)-CSUMPLQMOM12(8,IK)&
        +(CSUMPLQMOM12(13,IK)-CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)&
        -CSUMPLQMOM12(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE204(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)&
        -CSUMPLQ12(4)&
        +(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))&
        -(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))&
        -(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM204(N4,ID,IK)=(CSUMPLQMOM12(1,IK)-CSUMPLQMOM12(2,IK)&
        +CSUMPLQMOM12(3,IK)-CSUMPLQMOM12(4,IK)&
        -(CSUMPLQMOM12(9,IK)-CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)&
        -CSUMPLQMOM12(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE205(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)&
        -CSUMPLQ12(4)&
        -(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))&
        -(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))&
        +(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM205(N4,ID,IK)=(CSUMPLQMOM12(5,IK)-CSUMPLQMOM12(6,IK)&
        +CSUMPLQMOM12(7,IK)-CSUMPLQMOM12(8,IK)&
        -(CSUMPLQMOM12(13,IK)-CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)&
        -CSUMPLQMOM12(16,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 13
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE206(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)&
        +CSUMPLQ13(4)&
        +(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))&
        +(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))&
        +(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM206(N4,ID,IK)=(CSUMPLQMOM13(1,IK)+CSUMPLQMOM13(2,IK)&
        +CSUMPLQMOM13(3,IK)+CSUMPLQMOM13(4,IK)&
        +(CSUMPLQMOM13(9,IK)+CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)&
        +CSUMPLQMOM13(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE207(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)&
        +CSUMPLQ13(4)&
        -(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))&
        +(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))&
        -(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM207(N4,ID,IK)=(CSUMPLQMOM13(5,IK)+CSUMPLQMOM13(6,IK)&
        +CSUMPLQMOM13(7,IK)+CSUMPLQMOM13(8,IK)&
        +(CSUMPLQMOM13(13,IK)+CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)&
        +CSUMPLQMOM13(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE208(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)&
        +CSUMPLQ13(4)&
        +(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))&
        -(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))&
        -(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM208(N4,ID,IK)=(CSUMPLQMOM13(1,IK)+CSUMPLQMOM13(2,IK)&
        +CSUMPLQMOM13(3,IK)+CSUMPLQMOM13(4,IK)&
        -(CSUMPLQMOM13(9,IK)+CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)&
        +CSUMPLQMOM13(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE209(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)&
        +CSUMPLQ13(4)&
        -(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))&
        -(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))&
        +(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM209(N4,ID,IK)=(CSUMPLQMOM13(5,IK)+CSUMPLQMOM13(6,IK)&
        +CSUMPLQMOM13(7,IK)+CSUMPLQMOM13(8,IK)&
        -(CSUMPLQMOM13(13,IK)+CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)&
        +CSUMPLQMOM13(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE210(N4,ID)=(CSUMPLQ13(1)+GIOT*CSUMPLQ13(2)-CSUMPLQ13(3)&
        -GIOT*CSUMPLQ13(4)+(CSUMPLQ13(5)+GIOT*CSUMPLQ13(6)-CSUMPLQ13(7)&
        -GIOT*CSUMPLQ13(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM210(N4,ID,IK)=(CSUMPLQMOM13(1,IK)&
        +GIOT*CSUMPLQMOM13(2,IK)&
        -CSUMPLQMOM13(3,IK)-GIOT*CSUMPLQMOM13(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE211(N4,ID)=(CSUMPLQ13(1)+GIOT*CSUMPLQ13(2)-CSUMPLQ13(3)&
        -GIOT*CSUMPLQ13(4)-(CSUMPLQ13(5)+GIOT*CSUMPLQ13(6)-CSUMPLQ13(7)&
        -GIOT*CSUMPLQ13(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM211(N4,ID,IK)=(CSUMPLQMOM13(5,IK)&
        +GIOT*CSUMPLQMOM13(6,IK)&
        -CSUMPLQMOM13(7,IK)-GIOT*CSUMPLQMOM13(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE212(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)&
        -CSUMPLQ13(4)&
        +(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))&
        +(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))&
        +(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM212(N4,ID,IK)=(CSUMPLQMOM13(1,IK)-CSUMPLQMOM13(2,IK)&
        +CSUMPLQMOM13(3,IK)-CSUMPLQMOM13(4,IK)&
        +(CSUMPLQMOM13(9,IK)-CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)&
        -CSUMPLQMOM13(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE213(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)&
        -CSUMPLQ13(4)&
        -(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))&
        +(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))&
        -(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM213(N4,ID,IK)=(CSUMPLQMOM13(5,IK)-CSUMPLQMOM13(6,IK)&
        +CSUMPLQMOM13(7,IK)-CSUMPLQMOM13(8,IK)&
        +(CSUMPLQMOM13(13,IK)-CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)&
        -CSUMPLQMOM13(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE214(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)&
        -CSUMPLQ13(4)&
        +(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))&
        -(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))&
        -(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM214(N4,ID,IK)=(CSUMPLQMOM13(1,IK)-CSUMPLQMOM13(2,IK)&
        +CSUMPLQMOM13(3,IK)-CSUMPLQMOM13(4,IK)&
        -(CSUMPLQMOM13(9,IK)-CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)&
        -CSUMPLQMOM13(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE215(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)&
        -CSUMPLQ13(4)&
        -(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))&
        -(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))&
        +(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM215(N4,ID,IK)=(CSUMPLQMOM13(5,IK)-CSUMPLQMOM13(6,IK)&
        +CSUMPLQMOM13(7,IK)-CSUMPLQMOM13(8,IK)&
        -(CSUMPLQMOM13(13,IK)-CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)&
        -CSUMPLQMOM13(16,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 14
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE216(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)&
        +CSUMPLQ14(4)&
        +(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))&
        +(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))&
        +(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM216(N4,ID,IK)=(CSUMPLQMOM14(1,IK)+CSUMPLQMOM14(2,IK)&
        +CSUMPLQMOM14(3,IK)+CSUMPLQMOM14(4,IK)&
        +(CSUMPLQMOM14(9,IK)+CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)&
        +CSUMPLQMOM14(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE217(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)&
        +CSUMPLQ14(4)&
        -(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))&
        +(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))&
        -(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM217(N4,ID,IK)=(CSUMPLQMOM14(5,IK)+CSUMPLQMOM14(6,IK)&
        +CSUMPLQMOM14(7,IK)+CSUMPLQMOM14(8,IK)&
        +(CSUMPLQMOM14(13,IK)+CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)&
        +CSUMPLQMOM14(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE218(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)&
        +CSUMPLQ14(4)&
        +(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))&
        -(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))&
        -(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM218(N4,ID,IK)=(CSUMPLQMOM14(1,IK)+CSUMPLQMOM14(2,IK)&
        +CSUMPLQMOM14(3,IK)+CSUMPLQMOM14(4,IK)&
        -(CSUMPLQMOM14(9,IK)+CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)&
        +CSUMPLQMOM14(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE219(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)&
        +CSUMPLQ14(4)&
        -(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))&
        -(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))&
        +(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM219(N4,ID,IK)=(CSUMPLQMOM14(5,IK)+CSUMPLQMOM14(6,IK)&
        +CSUMPLQMOM14(7,IK)+CSUMPLQMOM14(8,IK)&
        -(CSUMPLQMOM14(13,IK)+CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)&
        +CSUMPLQMOM14(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE220(N4,ID)=(CSUMPLQ14(1)+GIOT*CSUMPLQ14(2)-CSUMPLQ14(3)&
        -GIOT*CSUMPLQ14(4)+(CSUMPLQ14(5)+GIOT*CSUMPLQ14(6)-CSUMPLQ14(7)&
        -GIOT*CSUMPLQ14(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM220(N4,ID,IK)=(CSUMPLQMOM14(1,IK)&
        +GIOT*CSUMPLQMOM14(2,IK)&
        -CSUMPLQMOM14(3,IK)-GIOT*CSUMPLQMOM14(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE221(N4,ID)=(CSUMPLQ14(1)+GIOT*CSUMPLQ14(2)-CSUMPLQ14(3)&
        -GIOT*CSUMPLQ14(4)-(CSUMPLQ14(5)+GIOT*CSUMPLQ14(6)-CSUMPLQ14(7)&
        -GIOT*CSUMPLQ14(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM221(N4,ID,IK)=(CSUMPLQMOM14(5,IK)&
        +GIOT*CSUMPLQMOM14(6,IK)&
        -CSUMPLQMOM14(7,IK)-GIOT*CSUMPLQMOM14(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE222(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)&
        -CSUMPLQ14(4)&
        +(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))&
        +(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))&
        +(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM222(N4,ID,IK)=(CSUMPLQMOM14(1,IK)-CSUMPLQMOM14(2,IK)&
        +CSUMPLQMOM14(3,IK)-CSUMPLQMOM14(4,IK)&
        +(CSUMPLQMOM14(9,IK)-CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)&
        -CSUMPLQMOM14(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE223(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)&
        -CSUMPLQ14(4)&
        -(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))&
        +(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))&
        -(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM223(N4,ID,IK)=(CSUMPLQMOM14(5,IK)-CSUMPLQMOM14(6,IK)&
        +CSUMPLQMOM14(7,IK)-CSUMPLQMOM14(8,IK)&
        +(CSUMPLQMOM14(13,IK)-CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)&
        -CSUMPLQMOM14(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE224(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)&
        -CSUMPLQ14(4)&
        +(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))&
        -(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))&
        -(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM224(N4,ID,IK)=(CSUMPLQMOM14(1,IK)-CSUMPLQMOM14(2,IK)&
        +CSUMPLQMOM14(3,IK)-CSUMPLQMOM14(4,IK)&
        -(CSUMPLQMOM14(9,IK)-CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)&
        -CSUMPLQMOM14(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE225(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)&
        -CSUMPLQ14(4)&
        -(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))&
        -(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))&
        +(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM225(N4,ID,IK)=(CSUMPLQMOM14(5,IK)-CSUMPLQMOM14(6,IK)&
        +CSUMPLQMOM14(7,IK)-CSUMPLQMOM14(8,IK)&
        -(CSUMPLQMOM14(13,IK)-CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)&
        -CSUMPLQMOM14(16,IK)))*ADIV2
        enddo
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PLAQUETTE OPERATORS 15
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !**********************************************************************      
    !     J=0, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE226(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)&
        +CSUMPLQ15(4)&
        +(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))&
        +(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))&
        +(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM226(N4,ID,IK)=(CSUMPLQMOM15(1,IK)+CSUMPLQMOM15(2,IK)&
        +CSUMPLQMOM15(3,IK)+CSUMPLQMOM15(4,IK)&
        +(CSUMPLQMOM15(9,IK)+CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)&
        +CSUMPLQMOM15(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE227(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)&
        +CSUMPLQ15(4)&
        -(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))&
        +(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))&
        -(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=+, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM227(N4,ID,IK)=(CSUMPLQMOM15(5,IK)+CSUMPLQMOM15(6,IK)&
        +CSUMPLQMOM15(7,IK)+CSUMPLQMOM15(8,IK)&
        +(CSUMPLQMOM15(13,IK)+CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)&
        +CSUMPLQMOM15(16,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE228(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)&
        +CSUMPLQ15(4)&
        +(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))&
        -(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))&
        -(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM228(N4,ID,IK)=(CSUMPLQMOM15(1,IK)+CSUMPLQMOM15(2,IK)&
        +CSUMPLQMOM15(3,IK)+CSUMPLQMOM15(4,IK)&
        -(CSUMPLQMOM15(9,IK)+CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)&
        +CSUMPLQMOM15(12,IK)))*ADIV2
        enddo
    !**********************************************************************      
    !     J=0, Pp=-, Pr=-, q=0
    !**********************************************************************
        ALINE229(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)&
        +CSUMPLQ15(4)&
        -(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))&
        -(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))&
        +(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=0, Pp=-, q=1,2
    !**********************************************************************
        do IK=1, 2
            ALINEMOM229(N4,ID,IK)=(CSUMPLQMOM15(5,IK)+CSUMPLQMOM15(6,IK)&
        +CSUMPLQMOM15(7,IK)+CSUMPLQMOM15(8,IK)&
        -(CSUMPLQMOM15(13,IK)+CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)&
        +CSUMPLQMOM15(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=1, Pr=+, q=0
    !**********************************************************************            
        ALINE230(N4,ID)=(CSUMPLQ15(1)+GIOT*CSUMPLQ15(2)-CSUMPLQ15(3)&
        -GIOT*CSUMPLQ15(4)+(CSUMPLQ15(5)+GIOT*CSUMPLQ15(6)-CSUMPLQ15(7)&
        -GIOT*CSUMPLQ15(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2
    !**********************************************************************                  
        do IK=1, 2
            ALINEMOM230(N4,ID,IK)=(CSUMPLQMOM15(1,IK)&
        +GIOT*CSUMPLQMOM15(2,IK)&
        -CSUMPLQMOM15(3,IK)-GIOT*CSUMPLQMOM15(4,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=1, Pr=-, q=0
    !**********************************************************************                  
        ALINE231(N4,ID)=(CSUMPLQ15(1)+GIOT*CSUMPLQ15(2)-CSUMPLQ15(3)&
        -GIOT*CSUMPLQ15(4)-(CSUMPLQ15(5)+GIOT*CSUMPLQ15(6)-CSUMPLQ15(7)&
        -GIOT*CSUMPLQ15(8)))*ADIV2
    !**********************************************************************      
    !     J=1, q=1,2  Here!
    !**********************************************************************                        
        do IK=1, 2
            ALINEMOM231(N4,ID,IK)=(CSUMPLQMOM15(5,IK)&
        +GIOT*CSUMPLQMOM15(6,IK)&
        -CSUMPLQMOM15(7,IK)-GIOT*CSUMPLQMOM15(8,IK))*ADIV1
        enddo
    !**********************************************************************      
    !     J=2, Pp=+, Pr=+, q=0
    !**********************************************************************
        ALINE232(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)&
        -CSUMPLQ15(4)&
        +(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))&
        +(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))&
        +(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM232(N4,ID,IK)=(CSUMPLQMOM15(1,IK)-CSUMPLQMOM15(2,IK)&
        +CSUMPLQMOM15(3,IK)-CSUMPLQMOM15(4,IK)&
        +(CSUMPLQMOM15(9,IK)-CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)&
        -CSUMPLQMOM15(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=+, Pr=-, q=0
    !**********************************************************************
        ALINE233(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)&
        -CSUMPLQ15(4)&
        -(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))&
        +(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))&
        -(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=+, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM233(N4,ID,IK)=(CSUMPLQMOM15(5,IK)-CSUMPLQMOM15(6,IK)&
        +CSUMPLQMOM15(7,IK)-CSUMPLQMOM15(8,IK)&
        +(CSUMPLQMOM15(13,IK)-CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)&
        -CSUMPLQMOM15(16,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=+, q=0
    !**********************************************************************
        ALINE234(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)&
        -CSUMPLQ15(4)&
        +(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))&
        -(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))&
        -(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM234(N4,ID,IK)=(CSUMPLQMOM15(1,IK)-CSUMPLQMOM15(2,IK)&
        +CSUMPLQMOM15(3,IK)-CSUMPLQMOM15(4,IK)&
        -(CSUMPLQMOM15(9,IK)-CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)&
        -CSUMPLQMOM15(12,IK)))*ADIV2
        enddo      
    !**********************************************************************      
    !     J=2, Pp=-, Pr=-, q=0  test
    !**********************************************************************
        ALINE235(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)&
        -CSUMPLQ15(4)&
        -(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))&
        -(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))&
        +(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
    !**********************************************************************      
    !     J=2, Pp=-, q=0
    !**********************************************************************
        do IK=1, 2
            ALINEMOM235(N4,ID,IK)=(CSUMPLQMOM15(5,IK)-CSUMPLQMOM15(6,IK)&
        +CSUMPLQMOM15(7,IK)-CSUMPLQMOM15(8,IK)&
        -(CSUMPLQMOM15(13,IK)-CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)&
        -CSUMPLQMOM15(16,IK)))*ADIV2
        enddo      
    !**********************************************************************
    !**********************************************************************
        RETURN
        end
    !**********************************************************************
end module