C***************************************************************
C***************************************************************
C     ******  ******  ***    ******  **      ******  *    *
C       **    *    *  *  *   *       **      *    *  **   *
C       **    *    *  * *    ******  **      *    *  * *  *
C       **    *    *  *  *   *       **      *    *  *  * *
C       **    ******  *   *  ******  ******  ******  *   **
C***************************************************************
C***************************************************************
C     Code for calculation the correlation functions of torelons
C     So far the code includes operators for k=1 N-ality
C     In the future there will be k=2 N-ality operatorsx
C***************************************************************
C***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
      PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1)
C
      PARAMETER(ICALLG=1,NITER=2,NMEASL=NITER/ICALLG)
      PARAMETER(ICMIN=134764,ICMAX=NITER+ICMIN-1)     
      PARAMETER(IBLOK=5,IBING=1,NUMBIN=2)
      PARAMETER(PARBS=0.30,PARBDS=0.12)
C
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)      
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------
C******************************************************************************
C******************************************************************************
C      
      COMMON/ARRAYS/U11(NCOL2*LSIZE*4)
      COMPLEX U11
      REAL*8 rndnum
C
      CHARACTER :: homepath*43
      CHARACTER :: conf_directory1*24
      CHARACTER :: conf_directory2*35
      CHARACTER :: directory*102
      CHARACTER :: file_name*120
      CHARACTER :: copy_file*256
      CHARACTER :: trnsf_file*256
      CHARACTER :: confnum*12
      CHARACTER :: name_of_file*6
C                    
C      homepath=trim('/home/dp208/dp208/dc-athe1/AXIONS/NF2/b2.3/')
      homepath='/gpfs/scratch/ehpc598/torelons/'
      conf_directory1='cnfg/'
C      conf_directory1=trim('m-1.0/26x26x26x52/confs/')
      conf_directory2='run1_52x26x26x26nc2rADJnf2b2.300000'
      directory=trim(homepath) // trim(conf_directory1) // 
     &          trim(conf_directory2)
C      
      CALL SETUP
      ISEED=3591
      CALL RINI(ISEED)
C
      ITOT=NITER
C
      WRITE(*,*) "                                                "
      WRITE(6,90)
90    FORMAT(' *******************************************************')
      WRITE(*,*) "                                                "
      WRITE(*,*)
     &     " ******  ******  ***    ******  **      ******  *    *"
      WRITE(*,*)
     &     "   **    *    *  *  *   *       **      *    *  **   *"
      WRITE(*,*)
     &     "   **    *    *  * *    ******  **      *    *  * *  *"
      WRITE(*,*)
     &     "   **    *    *  *  *   *       **      *    *  *  * *"
      WRITE(*,*)
     &     "   **    ******  *   *  ******  ******  ******  *   **"
      WRITE(*,*) "                                                "
      WRITE(6,90)
      WRITE(6,91)
 91   FORMAT(' *')
      WRITE(6,92) LX1,LX2,LX3,LX4
 92   FORMAT("[Info][Lattice Size]  ",'             LX   = ',4I4)
      WRITE(6,193) BETAG
 193  FORMAT("[Info][Value of Beta] ",'           beta   =',F8.4)
      WRITE(6,194) NCOL
 194  FORMAT("[Info][Colors]        ",' Number of Colors = ',I6)
      WRITE(6,195) IHEAT
 195  FORMAT("[Info][Thermalisation]",'   Therml. sweeps = ',I6)
      WRITE(6,91)
      WRITE(6,90)
      WRITE(6,91)
      WRITE(6,196) ITOT
 196  FORMAT("[Info][Measurements]",'   Number of MC Iterations = ',I6)
      WRITE(6,197) ICALLG
 197  FORMAT("[Info][Measurements]",'   Sweeps per measurements = ',I6)
      WRITE(6,198) IBING
 198  FORMAT("[Info][Measurements]",'      Measurements per bin = ',I6)
      WRITE(6,91)
      WRITE(6,90)
 200  FORMAT("[Info][Number of Operators q=0]",'  Total Number = ',I4)
      WRITE(6,91)
      WRITE(6,205)
      WRITE(6,206)
      WRITE(6,210) NOPJ0PP
      WRITE(6,211) NOPJ0PM
      WRITE(6,212) NOPJ0MP
      WRITE(6,213) NOPJ0MM
      WRITE(6,214) NOPJ1P
      WRITE(6,216) NOPJ1M
      WRITE(6,218) NOPJ2PP
      WRITE(6,219) NOPJ2PM
      WRITE(6,220) NOPJ2MP
      WRITE(6,221) NOPJ2MM
      WRITE(6,207)
      WRITE(6,200) NTOTAL
 205  FORMAT("[Info][Number of Operators q=0]",'     J   P   R      #')
 206  FORMAT("[Info][Number of Operators q=0]",'     ----------------')
 210  FORMAT("[Info][Number of Operators q=0]",'     0   +   +   ',I4)
 211  FORMAT("[Info][Number of Operators q=0]",'     0   +   -   ',I4)
 212  FORMAT("[Info][Number of Operators q=0]",'     0   -   +   ',I4)
 213  FORMAT("[Info][Number of Operators q=0]",'     0   -   -   ',I4)
 214  FORMAT("[Info][Number of Operators q=0]",'     1   #   +   ',I4)
 216  FORMAT("[Info][Number of Operators q=0]",'     1   #   -   ',I4)
 218  FORMAT("[Info][Number of Operators q=0]",'     2   +   +   ',I4)
 219  FORMAT("[Info][Number of Operators q=0]",'     2   +   -   ',I4)
 220  FORMAT("[Info][Number of Operators q=0]",'     2   -   +   ',I4)
 221  FORMAT("[Info][Number of Operators q=0]",'     2   -   -   ',I4)
 207  FORMAT("[Info][Number of Operators q=0]",'     ----------------')
      WRITE(6,91)
      WRITE(6,90)
 300  FORMAT("[Info][Number of Operators q=1,2]",'  Total Number = ',I4)
      WRITE(6,91)
      WRITE(6,305)
      WRITE(6,306)
      WRITE(6,310) NOPJ0PMOM
      WRITE(6,312) NOPJ0MMOM
      WRITE(6,314) NOPJ1MOM
      WRITE(6,318) NOPJ2PMOM
      WRITE(6,320) NOPJ2MMOM
      WRITE(6,307)
      WRITE(6,300) NTOTALMOM
 305  FORMAT("[Info][Number of Operators q=1,2]",'     J   P   #')
 306  FORMAT("[Info][Number of Operators q=1,2]",'     ---------------')
 310  FORMAT("[Info][Number of Operators q=1,2]",'     0   + ',I4)
 312  FORMAT("[Info][Number of Operators q=1,2]",'     0   - ',I4)
 314  FORMAT("[Info][Number of Operators q=1,2]",'     1     ',I4)
 318  FORMAT("[Info][Number of Operators q=1,2]",'     2   + ',I4)
 320  FORMAT("[Info][Number of Operators q=1,2]",'     2   - ',I4)
 307  FORMAT("[Info][Number of Operators q=1,2]",'     ---------------')
      WRITE(6,91)
      WRITE(6,90)
C
      IFILE = ICMIN
      WRITE(6,*) ICMIN 
C - 1
      DO 20 ITER=1,NITER
         IFILE = IFILE + 16
C
         ist=IFILE
         write(confnum, '(i0)') ist
         file_name=trim(directory)//'m1.000000n'//trim(confnum)
C     
         write(6,*) "------------------------------"
         write(6,*) "Execution of external commands"
         write(6,*) "                              "
         write(6,*) "Commands to be executed"
         write(6,*) "------------------------------"
   
         copy_file=trim('cp '//file_name//' ./conf')
         write(6,*) copy_file
         write(6,*) "                              "
         write(6,*) "                              "
         write(6,*) "*************************************************"
         call execute_command_line(copy_file, WAIT=.true.)
         write(6,*) "*************************************************"
         write(6,*) "                              "
         write(6,*) "                              "
         name_of_file=trim('./conf')
         WRITE(6,*) 'Config to be read: ', file_name
C         
         CALL cpu_time(t1)
         CALL READ_GF(name_of_file)
         CALL cpu_time(t2)
         WRITE(6,902) sngl(t2-t1)
 902  FORMAT('[Info][Time]','    Configuration Read time:', F8.4)
C     
      CALL ACTION(0,ITER,TOTACT)
 511  format('[FM][0]Check plaq = ', f8.6)
      write (6,511) TOTACT 
      CALL POLYLOOP()
      CALL MEASURE(ITER,NITER)
      call execute_command_line("rm conf", WAIT=.true.)
C
 20   CONTINUE
C
      CALL TODISK1
      CALL TODISK2
C
      STOP
      END
C*********************************************************************
      function quaternion_to_matrix(quaternion)
      double precision, dimension(4) :: quaternion
      complex, dimension(2, 2) :: quaternion_to_matrix

      quaternion_to_matrix(1,1) = sngl(quaternion(1)) + 
     &sngl(quaternion(4))*(0, 1)
      quaternion_to_matrix(1,2) = -sngl(quaternion(3)) + 
     &sngl(quaternion(2))*(0, 1)
      quaternion_to_matrix(2,1) = sngl(quaternion(3)) + 
     &sngl(quaternion(2))*(0, 1)
      quaternion_to_matrix(2,2) = sngl(quaternion(1)) - 
     &sngl(quaternion(4))*(0, 1)

      return
      end function quaternion_to_matrix
C*********************************************************************
      SUBROUTINE READ_GF(filename)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(NCOL=2,ND=4)
C
      CHARACTER :: filename*6
      integer(kind=4) nc_read, nx_read, ny_read, nz_read, nt_read
      integer t, x, y, z, dir, dir_target
      double precision plaquette_read

      double precision :: quaternion(4)
      complex, dimension(2, 2) :: quaternion_to_matrix

      COMMON/ARRAYS/U11(NCOL,NCOL,LX1,LX2,LX3,LX4,4)

      COMPLEX U11

      open(unit=1, file=filename, access='stream', 
     &form='unformatted', convert="BIG_ENDIAN")
      read(1) nc_read
      read(1) nt_read
      read(1) nx_read
      read(1) ny_read
      read(1) nz_read
      read(1) plaquette_read

      WRITE(6,101) plaquette_read
 101  FORMAT('[I/O][Plaq]', 'Plaquette value:',F8.6) 
      WRITE(6,102) nc_read
 102  FORMAT('[I/O][Ncol]', 'Number of Colors:',I2.1) 
      WRITE(6,103) nt_read, nx_read, ny_read, nz_read 
 103  FORMAT('[I/O][Dim]', 'T x X x Y x Z=',I3.2, I3.2, I3.2, I3.2) 
C
      DO t = 1, LX4
         DO x = 1, LX1
            DO y = 1, LX2
               DO z = 1, LX3
                  DO dir = 1, ND
                     read(1) quaternion

                     if (dir == 1) then
                        dir_target = 4
                     else
                        dir_target = dir-1
                     endif

      U11(1, 1, x, y, z, t, dir_target) =  
     &complex(sngl(quaternion(1)),sngl(quaternion(4)))
      U11(1, 2, x, y, z, t, dir_target) = 
     &complex(-sngl(quaternion(3)),sngl(quaternion(2)))  
      U11(2, 1, x, y, z, t, dir_target) =  
     &complex(sngl(quaternion(3)),sngl(quaternion(2))) 
      U11(2, 2, x, y, z, t, dir_target) =  
     &complex(sngl(quaternion(1)),-sngl(quaternion(4)))
C
                  enddo
               enddo
            enddo
         enddo
      enddo

      rewind(1)
      end subroutine READ_GF
C*********************************************************************
C*********************************************************************
C*********************************************************************
      SUBROUTINE READ_GAUGE_ETMC(filename) ! to be fixed !
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
      
      COMMON/ARRAYS/U11(NCOL,NCOL,LX1,LX2,LX3,LX4,4)

      COMPLEX U11
      COMPLEX*16 UR11(LX4,LX3,LX2,LX1,4,NCOL,NCOL)
      CHARACTER :: filename*6
C
      INQUIRE(IOLENGTH=l_rec) UR11(1,1,1,1,1,1,1)
C
      IUN=60
C
      OPEN(IUN,FILE=trim(filename),access="direct",form='unformatted',
     &CONVERT='BIG_ENDIAN', RECL=l_rec)
C
c      REWIND(60)
C
      ICALL=1
      DO L4=1, LX4
         DO L3=1, LX3
            DO L2=1, LX2
               DO L1=1, LX1
                  DO MU=1, 4
                     DO IJ=1, NCOL
                        DO IK=1, NCOL
C     C
                      READ(IUN,REC=ICALL) UR11(L4,L3,L2,L1,MU,IJ,IK)
C     C
                           ICALL=ICALL+1
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO     
C
      DO IJ=1, NCOL
         DO IK=1, NCOL
            DO L1=1, LX1
               DO L2=1, LX2
                  DO L3=1, LX3
                     DO L4=1, LX4
                        DO MU=1, 4
C     
                           DREALCONF=DREAL(UR11(L4,L3,L2,L1,MU,IJ,IK))
                           DIMAGCONF=DIMAG(UR11(L4,L3,L2,L1,MU,IJ,IK))
                           SREALCONF=sngl(DREALCONF)
                           SIMAGCONF=sngl(DIMAGCONF)
        U11(IJ,IK,L1,L2,L3,L4,MU)=complex(SREALCONF,SIMAGCONF)
C
C        write(*,*) U11(IJ,IK,L1,L2,L3,L4,MU)
C                                                                                                                                                                        
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO      
C
      WRITE(6,*) 'Configuration read'
C
      RETURN
      END            
C***********************************************************************
C***********************************************************************
C***********************************************************************
C                          THERMAL LINES
C***********************************************************************
C***********************************************************************
      SUBROUTINE POLYLOOP()
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NITER=2,NUMBIN=2)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      DIMENSION DUM11(NCOL2)
     &,A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2)
      COMPLEX U11,A11,B11,C11,D11,DUM11,ACT
C
      DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN),AVACS(NUMBIN),AVACT(NUMBIN)
      DIMENSION VAL(NUMBIN),AV(NUMBIN)
      dimension icoordvect(4)
      complex cpol1
      complex*16 cpol(4)
C
      icoordvect(1) = LX1-1
      icoordvect(2) = LX2-1
      icoordvect(3) = LX3-1
      icoordvect(4) = LX4-1
C     
      cpol(:) = dcmplx(0.0d0,0.0d0) 
      dnorm = 1.0/dfloat(lsizeb*ncol)
c
      DO ipoint=1, lsizeb
         do idir=1, 4
            a11(:) = u11(:,ipoint,idir)
            ipoint1 = ipoint
            do icoord=2,icoordvect(idir) 
               ipoint2 = iup(ipoint1,idir) 
               b11 = u11(:,ipoint2,idir)
               call vmx(1,a11,b11,c11,1)
               a11(:) = c11(:)
               ipoint2 = ipoint1
            enddo
            ipoint2 = iup(ipoint1,idir) 
            b11 = u11(:,ipoint2,idir)
            call trvmx(1,a11,b11,cpol1,1)
            cpol(idir) = cpol(idir) + cpol1
c            write (*,*) cpol1, cpol(idir)
         enddo
      enddo
C
      cpol(:) = cpol(:)*dnorm
C
 888  format('Poly(' ,i1, ')  =  (',f16.12,','f16.12')')
      do idir=1,4
         write (92,888) idir,dreal(cpol(idir)),dimag(cpol(idir))
      enddo
C
      RETURN
      END
C*********************************************************************
C***************************SUBROUTINE TODISK*************************
C*********************************************************************
      SUBROUTINE TODISK1
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IBLOK=5,NUMBIN=2,NITER=2)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LMAX=LX4/2+1)
      PARAMETER(LMAXIR=3)
C******************************************************************************
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------
C******************************************************************************
C******************************************************************************
C******************************************************************************  
      COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLP(NUMBIN,NTOTAL)
      COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP)
     &,AVACLJ0PP(NUMBIN,NOPJ0PP)
      COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM)
     &,AVACLJ0PM(NUMBIN,NOPJ0PM)
      COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP)
     &,AVACLJ0MP(NUMBIN,NOPJ0MP)
      COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM)
     &,AVACLJ0MM(NUMBIN,NOPJ0MM)
      COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P)
     &,AVACLJ1P(NUMBIN,NOPJ1P)
      COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M)
     &,AVACLJ1M(NUMBIN,NOPJ1M)
      COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP)
     &,AVACLJ2PP(NUMBIN,NOPJ2PP)
      COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM)
     &,AVACLJ2PM(NUMBIN,NOPJ2PM)
      COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP)
     &,AVACLJ2MP(NUMBIN,NOPJ2MP)
      COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM)
     &,AVACLJ2MM(NUMBIN,NOPJ2MM)
C
      COMMON/PBLOCKJ0/ACORLJ0(NUMBIN,LMAX,NOPFULJ0,NOPFULJ0)
     &,AVACLJ0(NUMBIN,NOPFULJ0)
      COMMON/PBLOCKJ1/ACORLJ1(NUMBIN,LMAX,NOPFULJ1,NOPFULJ1)
     &,AVACLJ1(NUMBIN,NOPFULJ1)      
      COMMON/PBLOCKJ2/ACORLJ2(NUMBIN,LMAX,NOPFULJ2,NOPFULJ2)
     &,AVACLJ2(NUMBIN,NOPFULJ2)      
C
       COMMON/PBLOCKJALL/ACORLALL(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLALL(NUMBIN,NTOTAL)      
C
      COMMON/POLYT/TLINE(NITER),SQTLINE(NITER)
      COMMON/ASTORE/ACTN(NUMBIN,6),PLAQ(NITER,6)
C
      COMPLEX*16 ACORLP
      COMPLEX*16 AVACLP
      COMPLEX*16 ACORLJ0PP
      COMPLEX*16 AVACLJ0PP
      COMPLEX*16 ACORLJ0PM
      COMPLEX*16 AVACLJ0PM
      COMPLEX*16 ACORLJ0MP
      COMPLEX*16 AVACLJ0MP
      COMPLEX*16 ACORLJ0MM
      COMPLEX*16 AVACLJ0MM
      COMPLEX*16 ACORLJ1P
      COMPLEX*16 AVACLJ1P
      COMPLEX*16 ACORLJ1M
      COMPLEX*16 AVACLJ1M
      COMPLEX*16 ACORLJ2PP
      COMPLEX*16 AVACLJ2PP
      COMPLEX*16 ACORLJ2PM
      COMPLEX*16 AVACLJ2PM
      COMPLEX*16 ACORLJ2MP
      COMPLEX*16 AVACLJ2MP
      COMPLEX*16 ACORLJ2MM
      COMPLEX*16 AVACLJ2MM
C
      COMPLEX*16 ACORLJ0
      COMPLEX*16 AVACLJ0
      COMPLEX*16 ACORLJ1
      COMPLEX*16 AVACLJ1
      COMPLEX*16 ACORLJ2
      COMPLEX*16 AVACLJ2      
      COMPLEX*16 ACORLALL
      COMPLEX*16 AVACLALL
C     
      COMPLEX TLINE
C
      OPEN (11,FILE='ACORLQ0.DAT')
      OPEN (12,FILE='ACORLJ0PPQ0.DAT')
      OPEN (13,FILE='ACORLJ0PMQ0.DAT')
      OPEN (14,FILE='ACORLJ0MPQ0.DAT')
      OPEN (15,FILE='ACORLJ0MMQ0.DAT')
      OPEN (16,FILE='ACORLJ1PQ0.DAT')
      OPEN (17,FILE='ACORLJ1MQ0.DAT')
      OPEN (18,FILE='ACORLJ2PPQ0.DAT')
      OPEN (19,FILE='ACORLJ2PMQ0.DAT')
      OPEN (20,FILE='ACORLJ2MPQ0.DAT')
      OPEN (21,FILE='ACORLJ2MMQ0.DAT')
      OPEN (22,FILE='ACORLJ0Q0.DAT')
      OPEN (23,FILE='ACORLJ1Q0.DAT')
      OPEN (24,FILE='ACORLJ2Q0.DAT')
      OPEN (25,FILE='ACORLJALL.DAT')
C******************************************
      DO IN=1, NUMBIN

C
         DO 31 IL=1,LMAXIR
C     
            DO IX=1, NTOTAL
               DO IY=1, NTOTAL
                  WRITE(11,*) ACORLP(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NTOTAL
               DO IY=1, NTOTAL
                  WRITE(25,*) ACORLALL(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
 31      CONTINUE
C     
         DO IL=1,LMAX
C
C****************************************************
C*************  J0 CORRELATORS **********************
C****************************************************
            DO IX=1, NOPJ0PP
               DO IY=1, NOPJ0PP
                  WRITE(12,*) ACORLJ0PP(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPJ0PM
               DO IY=1, NOPJ0PM
                  WRITE(13,*) ACORLJ0PM(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPJ0MP
               DO IY=1, NOPJ0MP
                  WRITE(14,*) ACORLJ0MP(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPJ0MM
               DO IY=1, NOPJ0MM
                  WRITE(15,*) ACORLJ0MM(IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
C*************  J1 CORRELATORS **********************
C****************************************************
            DO IX=1, NOPJ1P
               DO IY=1, NOPJ1P
                  WRITE(16,*) ACORLJ1P(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPJ1M
               DO IY=1, NOPJ1M
                  WRITE(17,*) ACORLJ1M(IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
C*************  J2 CORRELATORS **********************
C****************************************************
            DO IX=1, NOPJ2PP
               DO IY=1, NOPJ2PP
                  WRITE(18,*) ACORLJ2PP(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPJ2PM
               DO IY=1, NOPJ2PM
                  WRITE(19,*) ACORLJ2PM(IN,IL,IY,IX)
               ENDDO
            ENDDO
C
            DO IX=1, NOPJ2MP
               DO IY=1, NOPJ2MP
                  WRITE(20,*) ACORLJ2MP(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPJ2MM
               DO IY=1, NOPJ2MM
                  WRITE(21,*) ACORLJ2MM(IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
C*************  FULL J CORRELATORS ******************
C****************************************************
            DO IX=1, NOPFULJ0
               DO IY=1, NOPFULJ0
                  WRITE(22,*) ACORLJ0(IN,IL,IY,IX)
               ENDDO
            ENDDO
C                  
            DO IX=1, NOPFULJ1
               DO IY=1, NOPFULJ1
                  WRITE(23,*) ACORLJ1(IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPFULJ2
               DO IY=1, NOPFULJ2
                  WRITE(24,*) ACORLJ2(IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
C****************************************************
C****************************************************               
         ENDDO
      ENDDO
C******************************************
      REWIND(52)
      WRITE(52) ACTN
      WRITE(52) TLINE,PLAQ,SQTLINE

      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)
      CLOSE(20)
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)
      CLOSE(25)
C
      RETURN
      END
C*********************************************************************
C***************************SUBROUTINE TODISK*************************
C*********************************************************************
      SUBROUTINE TODISK2
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IBLOK=5,NUMBIN=2,NITER=2)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LMAX=LX4/2+1)
      PARAMETER(LMAXIR=3)
C******************************************************************************
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)      
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------
C******************************************************************************
C******************************************************************************
C******************************************************************************  
      COMMON/PBLOCKMOMP/ACORLPMOM(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM)
     &,AVACLPMOM(2,NUMBIN,NTOTALMOM)
C
      COMMON/PBLOCKMOMJ0P/
     &ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM)
     &,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
      COMMON/PBLOCKMOMJ0M/
     &ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM)
     &,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)      
      COMMON/PBLOCKMOMJ1/
     &ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM)
     &,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)     
      COMMON/PBLOCKMOMJ2P/
     &ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM)
     &,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
      COMMON/PBLOCKMOMJ2M/
     &ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM)
     &,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
C
      COMMON/PBLOCKMOMJ0/
     &ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM)
     &,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
       COMMON/PBLOCKMOMJ2/
     &ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM)
     &,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)      
C
       COMMON/PBLOCKMOMJALL/
     &ACORLPMOMJALL(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM)
     &,AVACLMOMALL(2,NUMBIN,NTOTALMOM) 
C*****************************************************************
      COMPLEX*16 ACORLPMOM
      COMPLEX*16 AVACLPMOM
      COMPLEX*16 ACORLPMOMJ0P
      COMPLEX*16 AVACLMOMJ0P
      COMPLEX*16 ACORLPMOMJ0M
      COMPLEX*16 AVACLMOMJ0M
      COMPLEX*16 ACORLPMOMJ1
      COMPLEX*16 AVACLMOMJ1
      COMPLEX*16 ACORLPMOMJ2P
      COMPLEX*16 AVACLMOMJ2P
      COMPLEX*16 ACORLPMOMJ2M
      COMPLEX*16 AVACLMOMJ2M
      COMPLEX*16 ACORLPMOMJ0
      COMPLEX*16 AVACLMOMJ0
      COMPLEX*16 ACORLPMOMJ2
      COMPLEX*16 AVACLMOMJ2
      COMPLEX*16 ACORLPMOMJALL
      COMPLEX*16 AVACLMOMALL
C*****************************************************************
      OPEN (17, FILE='ACORLPMOM1.DAT')
      OPEN (18, FILE='ACORLPMOM2.DAT')
      OPEN (19, FILE='ACORLJ0PMOMQ1.DAT')
      OPEN (20, FILE='ACORLJ0PMOMQ2.DAT')
      OPEN (21, FILE='ACORLJ0MMOMQ1.DAT')
      OPEN (22, FILE='ACORLJ0MMOMQ2.DAT')
      OPEN (23, FILE='ACORLJ1MOMQ1.DAT')
      OPEN (24, FILE='ACORLJ1MOMQ2.DAT')
      OPEN (25, FILE='ACORLJ2PMOMQ1.DAT')
      OPEN (26, FILE='ACORLJ2PMOMQ2.DAT')
      OPEN (27, FILE='ACORLJ2MMOMQ1.DAT')
      OPEN (28, FILE='ACORLJ2MMOMQ2.DAT')
      OPEN (29, FILE='ACORLJ0MOMQ1.DAT')
      OPEN (30, FILE='ACORLJ0MOMQ2.DAT')
      OPEN (33, FILE='ACORLJ2MOMQ1.DAT')
      OPEN (34, FILE='ACORLJ2MOMQ2.DAT')
      OPEN (35, FILE='ACORLJALLMOMQ1.DAT')
      OPEN (36, FILE='ACORLJALLMOMQ2.DAT')
C*****************************************************************
      DO IN=1,NUMBIN

C
         DO 31 IL=1,LMAXIR
C            
            DO IX=1, NTOTALMOM
               DO IY=1, NTOTALMOM
                  WRITE(17,*) ACORLPMOM(1,IN,IL,IY,IX)
                  WRITE(18,*) ACORLPMOM(2,IN,IL,IY,IX)
                  WRITE(35,*) ACORLPMOMJALL(1,IN,IL,IY,IX)
                  WRITE(36,*) ACORLPMOMJALL(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
 31      CONTINUE
         
         DO IL=1,LMAX
C****************************************************
C*************  J0 CORRELATORS **********************
C****************************************************
            DO IX=1,NOPJ0PMOM
               DO IY=1, NOPJ0PMOM
                  WRITE(19,*) ACORLPMOMJ0P(1,IN,IL,IY,IX)
                  WRITE(20,*) ACORLPMOMJ0P(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1,NOPJ0MMOM
               DO IY=1, NOPJ0MMOM
                  WRITE(21,*) ACORLPMOMJ0M(1,IN,IL,IY,IX)
                  WRITE(22,*) ACORLPMOMJ0M(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
C*************  J1 CORRELATORS **********************
C****************************************************
            DO IX=1,NOPJ1MOM
               DO IY=1, NOPJ1MOM
                  WRITE(23,*) ACORLPMOMJ1(1,IN,IL,IY,IX)
                  WRITE(24,*) ACORLPMOMJ1(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
C*************  J2 CORRELATORS **********************
C****************************************************
            DO IX=1,NOPJ2PMOM
               DO IY=1, NOPJ2PMOM
                  WRITE(25,*) ACORLPMOMJ2P(1,IN,IL,IY,IX)
                  WRITE(26,*) ACORLPMOMJ2P(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1,NOPJ2MMOM
               DO IY=1, NOPJ2MMOM
                  WRITE(27,*) ACORLPMOMJ2M(1,IN,IL,IY,IX)
                  WRITE(28,*) ACORLPMOMJ2M(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
C*************  FULL J CORRELATORS ******************
C****************************************************
            DO IX=1, NOPFULJ0MOM
               DO IY=1, NOPFULJ0MOM
                  WRITE(29,*) ACORLPMOMJ0(1,IN,IL,IY,IX)
                  WRITE(30,*) ACORLPMOMJ0(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C     
            DO IX=1, NOPFULJ2MOM
               DO IY=1, NOPFULJ2MOM
                  WRITE(33,*) ACORLPMOMJ2(1,IN,IL,IY,IX)
                  WRITE(34,*) ACORLPMOMJ2(2,IN,IL,IY,IX)
               ENDDO
            ENDDO
C****************************************************
         ENDDO
      ENDDO
C *** CLOSE OUTPUT FILES ***
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)
      CLOSE(20)
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)
      CLOSE(25)
      CLOSE(26)
      CLOSE(27)
      CLOSE(28)
      CLOSE(29)
      CLOSE(30)
      CLOSE(31)
      CLOSE(32)
      CLOSE(33)
      CLOSE(34)
      CLOSE(35)
      CLOSE(36)
C*******************************************************************
      RETURN
      END
C*******************************************************************
C***************************SUBROUTINE MEASURE**********************
C*******************************************************************
      SUBROUTINE MEASURE(ITER,NTOT)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ICALLG=1,IBING=1,NUMBIN=2)
C
      COMMON/ITEM/ITERG,JBING,NTOTG
C
      ITERG=ITER/ICALLG
      JBING=IBING
      NTOTG=IBING*NTOT/(ICALLG*IBING)
C
      IF(ITER.EQ.(ITERG*ICALLG).AND.ITERG.LE.NTOTG) THEN
        PAR=1.0d0
C
        CALL SETUPB
C
c        WRITE(*,*) "CALLING BLOCK"
        CALL BLOCK
        WRITE(*,*) "Number of iteration:", ITER
      ENDIF
C
      IF(ITER.EQ.NTOTG*ICALLG) THEN
        CALL MLINE1
        CALL MLINE2
        CALL MLINE3
        CALL MLINE4
        CALL MLINE5
      ENDIF
C
      RETURN
      END
C*********************************************************************
C*************************SUBROUTINE MLINE1***************************
C*********************************************************************
      SUBROUTINE MLINE1
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
      PARAMETER(LMAXIR=3)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C******************************************************************************
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)      
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------      
C
      COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLP(NUMBIN,NTOTAL)
C
      COMMON/ITEM/ITR,IBIN,ITOT
      COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
      COMPLEX*16 ACORLP,AVACLP
C      
      DIMENSION CORP(NUMBIN,LMAX)
      DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
C
      OPEN (23, FILE='DIAGNALFULLQ0.DAT')
C
      NBIN=ITOT/IBIN
      LS(1)=LX1
      LS(2)=LX2
      LS(3)=LX3
      EPS=0.0000000001d0
C********************************************************************
C     OPERATORS
C********************************************************************
      WRITE(23,80)
      WRITE(23,80)
80    FORMAT('****************************************************')
      WRITE(23,181)
181   FORMAT('****************** OPERATORS WITH J = 0 ************')
      WRITE(23,80)
      WRITE(23,81)
      WRITE(23,80)
      WRITE(23,80)
      WRITE(23,81)
81    FORMAT(' *')
      WRITE(23,91)
91    FORMAT(' AVERAGE  X,Y  LINES ')
      WRITE(23,81)
82    FORMAT(' ******************* ')
C
      DO 10 ID=1, NTOTAL
            ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO 13 IB=1,NUMBIN
               AV(IB)=REALPART(AVACLP(IB,ID))/ACOUNT
13          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO 14 IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLP(IB,ID))/ACOUNT
14          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101         FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10    CONTINUE
C
      WRITE(23,80)
      WRITE(23,81)
      WRITE(23,92)
92    FORMAT(' AVERAGE   Y   LINES ')
      WRITE(23,81)

      DO 20 ID=1, NTOTAL
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
94       FORMAT(' OPERATOR = ',I4)
         WRITE(23,82)
         DO NT=1,LMAXIR
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLP(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO20
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAXIR,CORP,VAC)
c         DO 26 NT=1,MAXDTLS
         DO 26 NT=1,LMAXIR-1            
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102      FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
26       CONTINUE
20    CONTINUE
      CLOSE(23)
C
      RETURN
      END
C*********************************************************************
C*************************SUBROUTINE MLINE3***************************
C*********************************************************************
      SUBROUTINE MLINE2
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
      PARAMETER(LMAXIR=3)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C******************************************************************************
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)      
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------
C
      COMMON/PBLOCKMOMP/ACORLPMOM(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM)
     &,AVACLPMOM(2,NUMBIN,NTOTALMOM)
C
      COMMON/ITEM/ITR,IBIN,ITOT
      COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
      COMPLEX*16 ACORLPMOM,AVACLPMOM
      DIMENSION CORP(NUMBIN,LMAX)
      DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
C
      OPEN (24, FILE='DIAGNALFULLMOMQ0.DAT')
C
      NBIN=ITOT/IBIN
      LS(1)=LX1
      LS(2)=LX2
      LS(3)=LX3
      EPS=0.0000000001d0
C
      DO IJ=1,2
C********************************************************************
C     POSITIVE PARITY OPERATORS
C********************************************************************
         WRITE(24,80)
         WRITE(24,80)
 80      FORMAT('****************************************************')
 181     FORMAT('************ OPERATORS ************')
         WRITE(24,80)
         WRITE(24,80)
         WRITE(24,81)
         WRITE(24,80)
         WRITE(24,181)
         WRITE(24,80)
         WRITE(24,81)
81       FORMAT(' *')
         WRITE(24,91)
91       FORMAT(' AVERAGE  X,Y  LINES ')
         WRITE(24,81)
82       FORMAT(' ******************* ')
C
      DO 10 ID=1, NTOTALMOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
         DO 13 IB=1,NUMBIN
            AV(IB)=REALPART(AVACLPMOM(IJ,IB,ID))/ACOUNT
13       CONTINUE
         CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
         DO 14 IB=1,NUMBIN
            AV(IB)=IMAGPART(AVACLPMOM(IJ,IB,ID))/ACOUNT
14       CONTINUE
         CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
         WRITE(24,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101      FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10    CONTINUE
C     
      WRITE(24,80)
      WRITE(24,81)
      WRITE(24,92)
92    FORMAT(' AVERAGE   Y   LINES ')
      WRITE(24,81)
      
      DO 20 ID=1, NTOTALMOM
         WRITE(24,81)
         WRITE(24,82)
         WRITE(24,94) ID
 94      FORMAT(' OPERATOR = ',I4)
         WRITE(24,82)
         DO NT=1,LMAXIR
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOM(IJ,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO20
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAXIR,CORP,VAC)
c
c         DO 26 NT=1,MAXDTLS
         DO 26 NT=1, LMAXIR-1
         WRITE(24,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
 102     FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
26       CONTINUE
20    CONTINUE
C
      ENDDO

      CLOSE(24)
C
      RETURN
      END
C*********************************************************************
C*************************SUBROUTINE MLINE3***************************
C*********************************************************************
      SUBROUTINE MLINE3
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
      PARAMETER(LMAXIR=3)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C******************************************************************************
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)      
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------
C******************************************************************************
C******************************************************************************
C******************************************************************************   
      COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLP(NUMBIN,NTOTAL)
C
      COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP)
     &,AVACLJ0PP(NUMBIN,NOPJ0PP)
      COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM)
     &,AVACLJ0PM(NUMBIN,NOPJ0PM)
      COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP)
     &,AVACLJ0MP(NUMBIN,NOPJ0MP)
      COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM)
     &,AVACLJ0MM(NUMBIN,NOPJ0MM)
      COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P)
     &,AVACLJ1P(NUMBIN,NOPJ1P)
      COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M)
     &,AVACLJ1M(NUMBIN,NOPJ1M)
      COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP)
     &,AVACLJ2PP(NUMBIN,NOPJ2PP)
      COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM)
     &,AVACLJ2PM(NUMBIN,NOPJ2PM)
      COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP)
     &,AVACLJ2MP(NUMBIN,NOPJ2MP)
      COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM)
     &,AVACLJ2MM(NUMBIN,NOPJ2MM)
C      
      COMMON/ITEM/ITR,IBIN,ITOT
      COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
      COMPLEX*16 ACORLP,AVACLP
C***********************************************************
      COMPLEX*16 ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
      COMPLEX*16 AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
C***********************************************************
      COMPLEX*16 ACORLJ1P,ACORLJ1M
      COMPLEX*16 AVACLJ1P,AVACLJ1M
C***********************************************************
      COMPLEX*16 ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
      COMPLEX*16 AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
C***********************************************************
      DIMENSION CORP(NUMBIN,LMAX)
      DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
C
      OPEN (23, FILE='DIAGNALQ0IND.DAT')
C
      NBIN=ITOT/IBIN
      LS(1)=LX1
      LS(2)=LX2
      LS(3)=LX3
      EPS=0.0000000001d0
C********************************************************************
C     OPERATORS
C********************************************************************
      WRITE(23,80)
      WRITE(23,80)
80    FORMAT('****************************************************')
      WRITE(23,181)
181   FORMAT('****************** OPERATORS WITH Q = 0 ************')
      WRITE(23,80)
      WRITE(23,81)
      WRITE(23,80)
      WRITE(23,80)
      WRITE(23,81)
81    FORMAT(' *')
      WRITE(23,91)
91    FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+, Pr=+ ')
      WRITE(23,81)
82    FORMAT(' ******************* ')
C
      DO 10 ID=1, NOPJ0PP
            ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO 13 IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ0PP(IB,ID))/ACOUNT
13          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO 14 IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ0PP(IB,ID))/ACOUNT
14          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101         FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10    CONTINUE
      WRITE(23,81)
      WRITE(23,201)
 201  FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+, Pr=- ')
      WRITE(23,81)
C
      DO ID=1, NOPJ0PM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ0PM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ0PM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
         ENDDO
      WRITE(23,81)
      WRITE(23,202)
 202  FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=-, Pr=+ ')
       WRITE(23,81)
C
      DO ID=1, NOPJ0MP
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ0MP(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ0MP(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,203)
 203  FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=-, Pr=- ')
       WRITE(23,81)
C
      DO ID=1, NOPJ0MM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ0MM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ0MM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,204)
 204  FORMAT(' AVERAGE  X,Y  LINES, J=1, Pr=+ ')
      WRITE(23,81)
C
      DO ID=1, NOPJ1P
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ1P(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ1P(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,205)
 205  FORMAT(' AVERAGE  X,Y  LINES, J=1, Pr=- ')
      WRITE(23,81)
C
      DO ID=1, NOPJ1M
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ1M(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ1M(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,206)
 206  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+, Pr=+ ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2PP
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ2PP(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ2PP(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,208)
 208  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+, Pr=- ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2PM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ2PM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ2PM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,209)
 209  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=-, Pr=+ ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2MP
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ2MP(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ2MP(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,210)
 210  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=-, Pr=- ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2MM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLJ2MM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLJ2MM(IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO      
C**********************************************************************
C     DIAGONAL CORRELATORS
C**********************************************************************      
      WRITE(23,80)
      
      WRITE(23,81)
      WRITE(23,220)
 220  FORMAT(' DIAGONAL TORELONS J=0, Pp=+, Pr=+, q=0 ')
      WRITE(23,81)

      DO 20 ID=1, NOPJ0PP
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
94       FORMAT(' OPERATOR = ',I4)
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ0PP(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO20
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
 102        FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
         ENDDO
 20   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,211)
 211  FORMAT(' DIAGONAL TORELONS J=0, Pp=+, Pr=-, q=0 ')
      WRITE(23,81)
C
      DO 21 ID=1, NOPJ0PM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ0PM(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO21
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 21   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,212)
 212  FORMAT(' DIAGONAL TORELONS J=0, Pp=-, Pr=+, q=0 ')
      WRITE(23,81)
C
      DO 22 ID=1, NOPJ0MP
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ0MP(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO22
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 22   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,213)
 213  FORMAT(' DIAGONAL TORELONS J=0, Pp=-, Pr=-, q=0 ')
      WRITE(23,81)
C
      DO 23 ID=1, NOPJ0MM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ0MM(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO23
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 23   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,214)
 214  FORMAT(' DIAGONAL TORELONS J=1, Pr=+, q=0 ')
      WRITE(23,81)
C
      DO 24 ID=1, NOPJ1P
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ1P(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO24
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 24   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,215)
 215  FORMAT(' DIAGONAL TORELONS J=1, Pr=-, q=0 ')
      WRITE(23,81)
C
      DO 25 ID=1, NOPJ1M
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ1M(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO25
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 25   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,216)
 216  FORMAT(' DIAGONAL TORELONS J=2, Pp=+, Pr=+, q=0 ')
      WRITE(23,81)
C
      DO 26 ID=1, NOPJ2PP
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ2PP(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO26
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 26   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,217)
 217  FORMAT(' DIAGONAL TORELONS J=2, Pp=+, Pr=-, q=0 ')
      WRITE(23,81)
C
      DO 27 ID=1, NOPJ2PM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ2PM(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO27
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 27   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,218)
 218  FORMAT(' DIAGONAL TORELONS J=2, Pp=-, Pr=+, q=0 ')
      WRITE(23,81)
C
      DO 28 ID=1, NOPJ2MP
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ2MP(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO28
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 28   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,219)
 219  FORMAT(' DIAGONAL TORELONS J=2, Pp=-, Pr=-, q=0 ')
      WRITE(23,81)
C
      DO 29 ID=1, NOPJ2MM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLJ2MM(IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO29
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 29   CONTINUE
C
      CLOSE(23)
C*********************************************************************   
      RETURN
      END
C*********************************************************************
C*************************SUBROUTINE MLINE4***************************
C*********************************************************************
      SUBROUTINE MLINE4
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
      PARAMETER(LMAXIR=3)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C******************************************************************************
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)      
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------
C******************************************************************************
C******************************************************************************
C******************************************************************************   
      COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLP(NUMBIN,NTOTAL)
C
      COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP)
     &,AVACLJ0PP(NUMBIN,NOPJ0PP)
      COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM)
     &,AVACLJ0PM(NUMBIN,NOPJ0PM)
      COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP)
     &,AVACLJ0MP(NUMBIN,NOPJ0MP)
      COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM)
     &,AVACLJ0MM(NUMBIN,NOPJ0MM)
      COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P)
     &,AVACLJ1P(NUMBIN,NOPJ1P)
      COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M)
     &,AVACLJ1M(NUMBIN,NOPJ1M)
      COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP)
     &,AVACLJ2PP(NUMBIN,NOPJ2PP)
      COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM)
     &,AVACLJ2PM(NUMBIN,NOPJ2PM)
      COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP)
     &,AVACLJ2MP(NUMBIN,NOPJ2MP)
      COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM)
     &,AVACLJ2MM(NUMBIN,NOPJ2MM)
C
      COMMON/PBLOCKMOMJ0P/
     &ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM)
     &,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
      COMMON/PBLOCKMOMJ0M/
     &ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM)
     &,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)      
      COMMON/PBLOCKMOMJ1/
     &ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM)
     &,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)     
      COMMON/PBLOCKMOMJ2P/
     &ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM)
     &,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
      COMMON/PBLOCKMOMJ2M/
     &ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM)
     &,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
C
      COMMON/PBLOCKMOMJ0/
     &ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM)
     &,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
       COMMON/PBLOCKMOMJ2/
     &ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM)
     &,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)
C
      COMMON/ITEM/ITR,IBIN,ITOT
      COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
      COMPLEX*16 ACORLP,AVACLP
C***********************************************************
      COMPLEX*16 ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
      COMPLEX*16 AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
C***********************************************************
      COMPLEX*16 ACORLJ1P,ACORLJ1M
      COMPLEX*16 AVACLJ1P,AVACLJ1M
C***********************************************************
      COMPLEX*16 ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
      COMPLEX*16 AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
C***********************************************************
      COMPLEX*16 ACORLPMOMJ0P, ACORLPMOMJ0M
      COMPLEX*16 AVACLMOMJ0P, AVACLMOMJ0M
C***********************************************************
      COMPLEX*16 ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
      COMPLEX*16 AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
C***********************************************************
      COMPLEX*16 ACORLPMOMJ2P, ACORLPMOMJ2M
      COMPLEX*16 AVACLMOMJ2P, AVACLMOMJ2M
C***********************************************************
      DIMENSION CORP(NUMBIN,LMAX)
      DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
C
      OPEN (23, FILE='DIAGNALQ1IND.DAT')
C
      NBIN=ITOT/IBIN
      LS(1)=LX1
      LS(2)=LX2
      LS(3)=LX3
      EPS=0.0000000001d0
C********************************************************************
C     OPERATORS
C********************************************************************
      WRITE(23,80)
      WRITE(23,80)
80    FORMAT('****************************************************')
      WRITE(23,181)
181   FORMAT('****************** OPERATORS WITH Q = 1 ************')
      WRITE(23,80)
      WRITE(23,81)
      WRITE(23,80)
      WRITE(23,80)
      WRITE(23,81)
81    FORMAT(' *')
      WRITE(23,91)
91    FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+')
      WRITE(23,81)
82    FORMAT(' ******************* ')
C
      DO 10 ID=1, NOPJ0PMOM
            ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO 13 IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ0P(1,IB,ID))/ACOUNT
13          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO 14 IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ0P(1,IB,ID))/ACOUNT
14          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101         FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10    CONTINUE
      WRITE(23,81)
      WRITE(23,202)
 202  FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=- ')
       WRITE(23,81)
C
      DO ID=1, NOPJ0MMOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ0M(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ0M(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,204)
 204  FORMAT(' AVERAGE  X,Y  LINES, J=1 ')
      WRITE(23,81)
C
      DO ID=1, NOPJ1MOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ1(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ1(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,206)
 206  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+ ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2PMOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ2P(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ2P(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,207)
 207  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=- ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2MMOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ2M(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ2M(1,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
C**********************************************************************
C     DIAGONAL CORRELATORS smirloglou 00302102473443, 00306972428749
C**********************************************************************      
      WRITE(23,80)
      
      WRITE(23,81)
      WRITE(23,220)
 220  FORMAT(' DIAGONAL TORELONS J=0, Pp=+, q=1 ')
      WRITE(23,81)

      DO 20 ID=1, NOPJ0PMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
94       FORMAT(' OPERATOR = ',I4)
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ0P(1,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO20
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
 102        FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
         ENDDO
 20   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,211)
 211  FORMAT(' DIAGONAL TORELONS J=0, Pp=-, q=1 ')
      WRITE(23,81)
C
      DO 21 ID=1, NOPJ0MMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ0M(1,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO21
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 21   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,212)
 212  FORMAT(' DIAGONAL TORELONS J=1, q=1 ')
      WRITE(23,81)
C
      DO 22 ID=1, NOPJ1MOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ1(1,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO22
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 22   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,213)
 213  FORMAT(' DIAGONAL TORELONS J=2, Pp=+, q=1 ')
      WRITE(23,81)
C
      DO 23 ID=1, NOPJ2PMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ2P(1,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO23
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 23   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,214)
 214  FORMAT(' DIAGONAL TORELONS J=2, Pr=-, q=1 ')
      WRITE(23,81)
C
      DO 24 ID=1, NOPJ2MMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ2M(1,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO24
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 24   CONTINUE
C
      CLOSE(23)
C*********************************************************************   
      RETURN
      END
C*********************************************************************
C*************************SUBROUTINE MLINE5***************************
C*********************************************************************
      SUBROUTINE MLINE5
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
      PARAMETER(LMAXIR=3)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C******************************************************************************
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)      
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)            
C----------------------------------------------------------
C******************************************************************************
C******************************************************************************
C******************************************************************************   
      COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLP(NUMBIN,NTOTAL)
C
      COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP)
     &,AVACLJ0PP(NUMBIN,NOPJ0PP)
      COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM)
     &,AVACLJ0PM(NUMBIN,NOPJ0PM)
      COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP)
     &,AVACLJ0MP(NUMBIN,NOPJ0MP)
      COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM)
     &,AVACLJ0MM(NUMBIN,NOPJ0MM)
      COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P)
     &,AVACLJ1P(NUMBIN,NOPJ1P)
      COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M)
     &,AVACLJ1M(NUMBIN,NOPJ1M)
      COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP)
     &,AVACLJ2PP(NUMBIN,NOPJ2PP)
      COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM)
     &,AVACLJ2PM(NUMBIN,NOPJ2PM)
      COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP)
     &,AVACLJ2MP(NUMBIN,NOPJ2MP)
      COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM)
     &,AVACLJ2MM(NUMBIN,NOPJ2MM)
C
      COMMON/PBLOCKMOMJ0P/
     &ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM)
     &,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
      COMMON/PBLOCKMOMJ0M/
     &ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM)
     &,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)      
      COMMON/PBLOCKMOMJ1/
     &ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM)
     &,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)     
      COMMON/PBLOCKMOMJ2P/
     &ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM)
     &,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
      COMMON/PBLOCKMOMJ2M/
     &ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM)
     &,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
C
      COMMON/PBLOCKMOMJ0/
     &ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM)
     &,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
       COMMON/PBLOCKMOMJ2/
     &ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM)
     &,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)
C
      COMMON/ITEM/ITR,IBIN,ITOT
      COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
      COMPLEX*16 ACORLP,AVACLP
C***********************************************************
      COMPLEX*16 ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
      COMPLEX*16 AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
C***********************************************************
      COMPLEX*16 ACORLJ1P,ACORLJ1M
      COMPLEX*16 AVACLJ1P,AVACLJ1M
C***********************************************************
      COMPLEX*16 ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
      COMPLEX*16 AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
C***********************************************************
      COMPLEX*16 ACORLPMOMJ0P, ACORLPMOMJ0M
      COMPLEX*16 AVACLMOMJ0P, AVACLMOMJ0M
C***********************************************************
      COMPLEX*16 ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
      COMPLEX*16 AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
C***********************************************************
      COMPLEX*16 ACORLPMOMJ2P, ACORLPMOMJ2M
      COMPLEX*16 AVACLMOMJ2P, AVACLMOMJ2M
C***********************************************************
      DIMENSION CORP(NUMBIN,LMAX)
      DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
C
      OPEN (23, FILE='DIAGNALQ2IND.DAT')
C
      NBIN=ITOT/IBIN
      LS(1)=LX1
      LS(2)=LX2
      LS(3)=LX3
      EPS=0.0000000001d0
C********************************************************************
C     OPERATORS
C********************************************************************
      WRITE(23,80)
      WRITE(23,80)
80    FORMAT('****************************************************')
      WRITE(23,181)
181   FORMAT('****************** OPERATORS WITH Q = 2 ************')
      WRITE(23,80)
      WRITE(23,81)
      WRITE(23,80)
      WRITE(23,80)
      WRITE(23,81)
81    FORMAT(' *')
      WRITE(23,91)
91    FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+')
      WRITE(23,81)
82    FORMAT(' ******************* ')
C
      DO 10 ID=1, NOPJ0PMOM
            ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO 13 IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ0P(2,IB,ID))/ACOUNT
13          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO 14 IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ0P(2,IB,ID))/ACOUNT
14          CONTINUE
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101         FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10    CONTINUE
      WRITE(23,81)
      WRITE(23,202)
 202  FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=- ')
       WRITE(23,81)
C
      DO ID=1, NOPJ0MMOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ0M(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ0M(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,204)
 204  FORMAT(' AVERAGE  X,Y  LINES, J=1 ')
      WRITE(23,81)
C
      DO ID=1, NOPJ1MOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ1(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ1(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,206)
 206  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+ ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2PMOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ2P(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ2P(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
      WRITE(23,81)
      WRITE(23,207)
 207  FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=- ')
      WRITE(23,81)
C
      DO ID=1, NOPJ2MMOM
         ACOUNT=NCOL*IBIN*LSIZE/LS(1)
            DO IB=1,NUMBIN
               AV(IB)=REALPART(AVACLMOMJ2M(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
            DO IB=1,NUMBIN
               AV(IB)=IMAGPART(AVACLMOMJ2M(2,IB,ID))/ACOUNT
            ENDDO
            CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
            WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
      ENDDO
C**********************************************************************
C     DIAGONAL CORRELATORS smirloglou 00302102473443, 00306972428749
C**********************************************************************      
      WRITE(23,80)
      
      WRITE(23,81)
      WRITE(23,220)
 220  FORMAT(' DIAGONAL TORELONS J=0, Pp=+, q=2 ')
      WRITE(23,81)

      DO 20 ID=1, NOPJ0PMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
94       FORMAT(' OPERATOR = ',I4)
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ0P(2,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO20
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
 102        FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
         ENDDO
 20   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,211)
 211  FORMAT(' DIAGONAL TORELONS J=0, Pp=-, q=2 ')
      WRITE(23,81)
C
      DO 21 ID=1, NOPJ0MMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ0M(2,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO21
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 21   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,212)
 212  FORMAT(' DIAGONAL TORELONS J=1, q=2 ')
      WRITE(23,81)
C
      DO 22 ID=1, NOPJ1MOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ1(2,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO22
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 22   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,213)
 213  FORMAT(' DIAGONAL TORELONS J=2, Pp=+, q=2 ')
      WRITE(23,81)
C
      DO 23 ID=1, NOPJ2PMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ2P(2,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO23
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 23   CONTINUE
C*********************************************************************
      WRITE(23,81)
      WRITE(23,214)
 214  FORMAT(' DIAGONAL TORELONS J=2, Pr=-, q=2 ')
      WRITE(23,81)
C
      DO 24 ID=1, NOPJ2MMOM
         WRITE(23,81)
         WRITE(23,82)
         WRITE(23,94) ID
         WRITE(23,82)
         DO NT=1,LMAX
            DO IB=1,NBIN
               CORP(IB,NT)=REALPART(ACORLPMOMJ2M(2,IB,NT,ID,ID))
            ENDDO
         ENDDO
         IF(ABS(CORP(1,1)).LE.EPS)GOTO24
         CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
         DO NT=1,MAXDTLS
            WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
         ENDDO
 24   CONTINUE
C
      CLOSE(23)
C*********************************************************************   
      RETURN
      END      
C*********************************************************************
C************************ SUBROUTINE BLOCK *************************** 
C*********************************************************************
      SUBROUTINE BLOCK
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(IBLOK=5,NBLOK=2,NSMEAR=12)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL2,LSIZEB,LX4,4)
      COMMON/ARRAYB/UB11(NCOL2,LSIZEB,3,IBLOK)
      COMMON/ASMEAR1/UC11(NCOL2,LSIZEB,3)
      COMMON/ASMEAR2/IUP(LSIZEB,3),IDN(LSIZEB,3)
      COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
      DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2)
      COMPLEX U11,UB11,A11,B11,C11,UC11
C
      CALL cpu_time(t1)
C *** For every time slice
      DO 1 I4=1,LX4
C
         DO 3 IBL=1,IBLOK
C
C *** For the first blocking level

            IF(IBL.EQ.1)THEN
                DO MU=1,3
                    DO NN=1,LSIZEB
                        DO IJ=1,NCOL2
C
C *** First blocking level: a copy of the configuration ***
C                  
                            UB11(IJ,NN,MU,1)=U11(IJ,NN,I4,MU)
                        ENDDO
                    ENDDO
                ENDDO
C *** Calculate the thermal line and go to next blocking level                
            GOTO11
            ENDIF
C
         IBLM=IBL-1
         DO KK=1,3
            DO NN=1,LSIZEB
               IUP(NN,KK)=IUPB(NN,KK,IBLM)
               IDN(NN,KK)=IDNB(NN,KK,IBLM)
               DO IJ=1,NCOL2
                  UC11(IJ,NN,KK)=UB11(IJ,NN,KK,IBLM)
               ENDDO
            ENDDO
         ENDDO
         CALL SMEAR1
C
         DO MU=1,3
            DO NN=1,LSIZEB
               M1=NN
               M2=IUP(M1,MU)
C
               DO 22 IJ=1,NCOL2
                  A11(IJ)=UC11(IJ,M1,MU)
                  B11(IJ)=UC11(IJ,M2,MU)
 22            CONTINUE
               CALL VMX(1,A11,B11,C11,1)
               DO 24 IJ=1,NCOL2
                  UB11(IJ,NN,MU,IBL)=C11(IJ)
 24            CONTINUE
C
            ENDDO
         ENDDO
11       CONTINUE
C
         CALL THERML1(I4,IBL)
C
C************************************************
 3       CONTINUE
 1    CONTINUE
      CALL cpu_time(t2)
      WRITE(6,901) sngl(t2-t1)
 901  FORMAT('[Info][Time]','     Thermal Line time:', F8.4)      
C
      CALL cpu_time(t1)
      CALL POT
      CALL cpu_time(t2)
      WRITE(6,902) sngl(t2-t1)
 902  FORMAT('[Info][Time]','     Correlation Creation time:', F8.4)      
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE SMEAR1
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(PARBS=0.30,PARBDS=0.12)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
      PARAMETER(IDIAG=1)
C
      COMMON/ASMEAR1/UC11(NCOL2,LSIZEB,3)
      COMMON/ASMEAR2/IUP(LSIZEB,3),IDN(LSIZEB,3)
      COMMON/DIAGIN/UUC11(NCOL2,LSIZEB,3),
     &IIUP(LSIZEB,3),IIDN(LSIZEB,3)
      COMMON/DIAGOUT/UDD(NCOL2,LSIZEB,4),IDD(LSIZEB,4)
      DIMENSION UCC11(NCOL2,LSIZEB,3)
      DIMENSION UINT11(NCOL2),UREN11(NCOL2)
      DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2),DUM11(NCOL2)
      COMPLEX A11,B11,C11,UINT11,UREN11,UC11,DUM11,UCC11,UDD
      COMPLEX UUC11
C
      DO KK=1,3
         DO NN=1,LSIZEB
            IIUP(NN,KK)=IUP(NN,KK)
            IIDN(NN,KK)=IDN(NN,KK)
         ENDDO
      ENDDO
C
      DO MU=1,3
         DO NN=1,LSIZEB
            DO IJ=1,NCOL2
               UUC11(IJ,NN,MU)=UC11(IJ,NN,MU)
            ENDDO
         ENDDO
      ENDDO
C
      DO 1 MU=1,3
         IF(IDIAG.EQ.1) CALL DIAG(MU)
C
         DO 2 NN=1,LSIZEB
            M1=NN
C
            DO 10 IJ=1,NCOL2
               UINT11(IJ)=UC11(IJ,M1,MU)
10          CONTINUE
C
            DO 40 NU=1,3
            IF(MU.EQ.NU)GOTO40
C
            DO 26 IJ=1,NCOL2
               A11(IJ)=UC11(IJ,M1,NU)
 26         CONTINUE
            M2=IUP(M1,NU)
            DO 27 IJ=1,NCOL2
               B11(IJ)=UC11(IJ,M2,MU)
 27         CONTINUE
            CALL VMX(1,A11,B11,C11,1)
            M3=IUP(M1,MU)
            DO 28 IJ=1,NCOL2
               B11(IJ)=UC11(IJ,M3,NU)
 28         CONTINUE
            CALL HERM(1,B11,DUM11,1)
            CALL VMX(1,C11,B11,A11,1)
            DO 21 IJ=1,NCOL2
               UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARBS
 21         CONTINUE
C
            M2=IDN(M1,NU)
            DO 31 IJ=1,NCOL2
               A11(IJ)=UC11(IJ,M2,NU)
 31         CONTINUE
            DO 32 IJ=1,NCOL2
               B11(IJ)=UC11(IJ,M2,MU)
 32         CONTINUE
            CALL HERM(1,A11,DUM11,1)
            CALL VMX(1,A11,B11,C11,1)
            M3=IUP(M2,MU)
            DO 33 IJ=1,NCOL2
               B11(IJ)=UC11(IJ,M3,NU)
 33         CONTINUE
            CALL VMX(1,C11,B11,A11,1)
            DO 35 IJ=1,NCOL2
               UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARBS
 35         CONTINUE
C
 40      CONTINUE
C
         DO 50 MNU=1,4
            IF(IDIAG.NE.1)GOTO50
c
            DO 42 IJ=1,NCOL2
               A11(IJ)=UDD(IJ,M1,MNU)
 42         CONTINUE
            M2=IDD(M1,MNU)
            DO 43 IJ=1,NCOL2
               B11(IJ)=UC11(IJ,M2,MU)
 43         CONTINUE
            CALL VMX(1,A11,B11,C11,1)
            M3=IUP(M1,MU)
            DO 44 IJ=1,NCOL2
               B11(IJ)=UDD(IJ,M3,MNU)
 44         CONTINUE
            CALL HERM(1,B11,DUM11,1)
            CALL VMX(1,C11,B11,A11,1)
            DO 45 IJ=1,NCOL2
               UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARBDS
 45         CONTINUE
C
 50      CONTINUE
C
         CALL RENORMBS(UINT11,UREN11)
         DO 60 IJ=1,NCOL2
            UCC11(IJ,M1,MU)=UREN11(IJ)
 60      CONTINUE
C
 2    CONTINUE
 1    CONTINUE
C
      DO MU=1,3
         DO NN=1,LSIZEB
            DO IJ=1,NCOL2
               UC11(IJ,NN,MU)=UCC11(IJ,NN,MU)
            ENDDO
         ENDDO
      ENDDO
      CALL DODET
C
      RETURN
      END
C**********************************************************
C  THIS ROUTINE MAKES THE SMEARED LINKS HAVE DET=1
C**********************************************************
      SUBROUTINE DODET
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ASMEAR1/U11(NCOL,NCOL,LSIZEB,3)
      COMPLEX ADUM(NCOL,NCOL),U11,CSUM
C
      COMPLEX AA(NCOL,NCOL)
      complex DET,CNORM,CDET,DNCOL
C
      DO MU=1,3 
         DO NN=1,LSIZEB
C
            DO J=1,NCOL
               DO I=1,NCOL
                  ADUM(I,J)=U11(I,J,NN,MU)
               ENDDO
            ENDDO
C
            DO J=1, NCOL
               DO I=1,NCOL
                  AA(I,J)=U11(I,J,NN,MU)
               ENDDO
            ENDDO
            JMAT=NCOL
            CALL DETNANT(JMAT,DET,AA)
c
            DNCOL=CMPLX(1.0/NCOL)
            CDET=DET**DNCOL
            CNORM=1.0/CDET
            DO J=1,NCOL
               DO I=1, NCOL
                  AA(I,J)=ADUM(I,J)*CNORM
                  U11(I,J,NN,MU)=AA(I,J)
c
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
      RETURN 
      END
C*********************************************************************
C form 4 'diagonal' links in plane orth to MU, from each site 
C with end-pt in IDD using links and pointers in /DIAGIN/
C -- output matrices (not suN) and pointers in /DIAGOUT/
C*********************************************************************
      SUBROUTINE DIAG(MU)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZEB=LX1*LX2*LX3)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C     
      COMMON/DIAGOUT/UDD(NCOL2,LSIZEB,4),IDD(LSIZEB,4)
      COMMON/DIAGIN/UC11(NCOL2,LSIZEB,3),
     &     IUP(LSIZEB,3),IDN(LSIZEB,3)
      COMPLEX C11(NCOL2),F11(NCOL2),DUM11(NCOL2)
      COMPLEX A11(NCOL2),D11(NCOL2),E11(NCOL2)
      COMPLEX AA11(NCOL2),DD11(NCOL2),EE11(NCOL2)
      COMPLEX UDD,UC11
C     
      NU=MU+1
      IF(MU.EQ.3)NU=1
      KU=6-MU-NU
C     
      DO 1 NN=1,LSIZEB
         M1=NN
C     
         DO 2 IJ=1,NCOL2
            A11(IJ)=UC11(IJ,M1,NU)
 2       CONTINUE
         M2=IUP(M1,NU)
         DO 3 IJ=1,NCOL2
            AA11(IJ)=UC11(IJ,M2,KU)
 3       CONTINUE
         CALL VMX(1,A11,AA11,C11,1)
         DO 4 IJ=1,NCOL2
            D11(IJ)=UC11(IJ,M1,KU)
 4       CONTINUE
         M3=IUP(M1,KU)
         DO 5 IJ=1,NCOL2
            DD11(IJ)=UC11(IJ,M3,NU)
 5       CONTINUE
         CALL VMX(1,D11,DD11,F11,1)
         DO 6 IJ=1,NCOL2
            UDD(IJ,M1,1)=C11(IJ)+F11(IJ)
 6       CONTINUE
         M9=IUP(M3,NU)
         IDD(M1,1)=M9
C     
         DO 61 IJ=1,NCOL2
            AA11(IJ)=C11(IJ)+F11(IJ)
 61      CONTINUE
         CALL HERM(1,AA11,DUM11,1)
         DO 62 IJ=1,NCOL2
            UDD(IJ,M9,3)=AA11(IJ)
 62      CONTINUE
         IDD(M9,3)=M1
C     
         M4=IDN(M2,KU)
         DO 7 IJ=1,NCOL2
            AA11(IJ)=UC11(IJ,M4,KU)
 7       CONTINUE
         CALL HERM(1,AA11,DUM11,1)
         CALL VMX(1,A11,AA11,C11,1)
         M5=IDN(M1,KU)
         DO 8 IJ=1,NCOL2
            E11(IJ)=UC11(IJ,M5,KU)
 8       CONTINUE
         CALL HERM(1,E11,DUM11,1)
         DO 9 IJ=1,NCOL2
            EE11(IJ)=UC11(IJ,M5,NU)
 9       CONTINUE
         CALL VMX(1,E11,EE11,F11,1)
         DO 10 IJ=1,NCOL2
            UDD(IJ,M1,2)=C11(IJ)+F11(IJ)
 10      CONTINUE
         IDD(M1,2)=M4
C     
         DO 63 IJ=1,NCOL2
            AA11(IJ)=C11(IJ)+F11(IJ)
 63      CONTINUE
         CALL HERM(1,AA11,DUM11,1)
         DO 64 IJ=1,NCOL2
            UDD(IJ,M4,4)=AA11(IJ)
 64      CONTINUE
         IDD(M4,4)=M1
C     
 1    CONTINUE
C     
      RETURN
      END
C**********************************************************
C  THIS ROUTINE MAKES BLOCKED MATRICES SU(N) IN CRUDEST   *
C               POSSIBLE WAY - 1 COOL                     *
C**********************************************************
      SUBROUTINE RENORMBS(UU1,UREN11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C     
      DIMENSION UU1(NCOL,NCOL),UREN11(NCOL,NCOL)
      COMPLEX UB1,UU1,CSUM,ADUM(NCOL,NCOL),UREN11
      COMPLEX A11(NCOL2),B11(NCOL2),C11(NCOL2)
C
      DO N2=1,NCOL
         DO N1=1,NCOL
            ADUM(N1,N2)=UU1(N1,N2)
         ENDDO
      ENDDO
C     
      DO 20 N2=1,NCOL
C     
         DO 10 N3=1,N2-1
C     
            CSUM=(0.0,0.0)
            DO 5 N1=1,NCOL
               CSUM=CSUM+ADUM(N1,N2)*CONJG(ADUM(N1,N3))
 5          CONTINUE
            DO 6 N1=1,NCOL
               ADUM(N1,N2)=ADUM(N1,N2)-CSUM*ADUM(N1,N3)
 6          CONTINUE
C     
 10      CONTINUE
C     
         SUM=0.0
         DO 7 N1=1,NCOL
            SUM=SUM+ADUM(N1,N2)*CONJG(ADUM(N1,N2))
 7       CONTINUE
         ANORM=1.0/SQRT(SUM)
         DO 8 N1=1,NCOL
            ADUM(N1,N2)=ADUM(N1,N2)*ANORM
 8       CONTINUE
C     
 20   CONTINUE
C
      DO IK=1,NCOL
         DO IJ=1,NCOL
            IADD=NCOL*(IK-1)
            A11(IJ+IADD)=ADUM(IJ,IK)
         ENDDO
      ENDDO
C
      DO 22 J=1,NCOL2
         C11(J)=A11(J)
22    CONTINUE
      DO J=1,NCOL
         DO I=1,NCOL
            IJ=I+(J-1)*NCOL
            A11(IJ)=CONJG(UU1(J,I))
         ENDDO
      ENDDO
C
      CALL VMX(1,A11,C11,B11,1)
      CALL RUNGRB(B11,C11)
C
      DO J=1,NCOL
         DO I=1,NCOL
            IJ=I+(J-1)*NCOL
            UREN11(I,J)=C11(IJ)
         ENDDO
      ENDDO
C
      RETURN
      END
C *******************************************************************
C *************         SUBROUTINE RUNGRB           *****************
C *******************************************************************
      SUBROUTINE RUNGRB(B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C     
      COMPLEX B11(NCOL2),C11(NCOL2)
      COMMON/SUBB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
C     
      DO 10 LDU=1,NCOL-1
         DO 15 LDL=LDU+1,NCOL
            II1=LDU
            JJ1=LDU
            II2=LDL
            JJ2=LDL
            II3=LDL
            JJ3=LDU
            II4=LDU
            JJ4=LDL
            CALL SUBGRB(LDU,LDL,B11,C11)
 15      CONTINUE
 10   CONTINUE
C     
      RETURN
      END
C*********************************************************************
C*******************      SUBROUTINE SUBGRB      *********************
C*********************************************************************
      SUBROUTINE SUBGRB(LDU,LDL,B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C     
      COMMON/SUBB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
C     
      COMPLEX B11(NCOL,NCOL),C11(NCOL,NCOL)
      COMPLEX F11,F12,A11(NCOL,NCOL),S11(NCOL,NCOL)
      COMPLEX T11(NCOL,NCOL)
C     
      F11=(B11(II1,JJ1)+CONJG(B11(II2,JJ2)))*0.5
      F12=(B11(II3,JJ3)-CONJG(B11(II4,JJ4)))*0.5
      UMAG=SQRT(F11*CONJG(F11)+F12*CONJG(F12))
      UMAG=1./UMAG
      F11=F11*UMAG
      F12=F12*UMAG
C     
      A11(II1,JJ1)=CONJG(F11)
      A11(II4,JJ4)=CONJG(F12)
      A11(II3,JJ3)=-F12
      A11(II2,JJ2)=F11
C     
      DO 3 IJ1=1,NCOL
         S11(IJ1,LDU)=C11(IJ1,LDU)*(A11(LDU,LDU)-1.0)
     &        +C11(IJ1,LDL)*A11(LDL,LDU)
         S11(IJ1,LDL)=C11(IJ1,LDU)*A11(LDU,LDL)
     &        +C11(IJ1,LDL)*(A11(LDL,LDL)-1.0)
         T11(IJ1,LDU)=B11(IJ1,LDU)*(A11(LDU,LDU)-1.0)
     &        +B11(IJ1,LDL)*A11(LDL,LDU)
         T11(IJ1,LDL)=B11(IJ1,LDU)*A11(LDU,LDL)
     &        +B11(IJ1,LDL)*(A11(LDL,LDL)-1.0)
 3    CONTINUE
      DO 4 IJ1=1,NCOL
         C11(IJ1,LDU)=C11(IJ1,LDU)+S11(IJ1,LDU)
         C11(IJ1,LDL)=C11(IJ1,LDL)+S11(IJ1,LDL)
         B11(IJ1,LDU)=B11(IJ1,LDU)+T11(IJ1,LDU)
         B11(IJ1,LDL)=B11(IJ1,LDL)+T11(IJ1,LDL)
 4    CONTINUE
C     
      RETURN
      END
C*********************************************************************
C************* HERE I CREATE THE PULSE - LIKE OPERATORS **************
C*********************************************************************
      SUBROUTINE THERML1(N4,IBLL)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C********************************************************
      COMMON/LINES/ALINE1(LX4,IBLOK),ALINE2(LX4,IBLOK)
     &,ALINE3(LX4,IBLOK),ALINE4(LX4,IBLOK),ALINE5(LX4,IBLOK)
     &,ALINE6(LX4,IBLOK),ALINE7(LX4,IBLOK),ALINE8(LX4,IBLOK) 
     &,ALINE9(LX4,IBLOK),ALINE10(LX4,IBLOK),ALINE11(LX4,IBLOK)
     &,ALINE12(LX4,IBLOK),ALINE13(LX4,IBLOK),ALINE14(LX4,IBLOK)
     &,ALINE15(LX4,IBLOK),ALINE16(LX4,IBLOK),ALINE17(LX4,IBLOK)
     &,ALINE18(LX4,IBLOK),ALINE19(LX4,IBLOK),ALINE20(LX4,IBLOK)
     &,ALINE21(LX4,IBLOK),ALINE22(LX4,IBLOK),ALINE23(LX4,IBLOK)
     &,ALINE24(LX4,IBLOK),ALINE25(LX4,IBLOK),ALINE26(LX4,IBLOK)
     &,ALINE27(LX4,IBLOK),ALINE28(LX4,IBLOK),ALINE29(LX4,IBLOK)
     &,ALINE30(LX4,IBLOK),ALINE31(LX4,IBLOK),ALINE32(LX4,IBLOK)
     &,ALINE33(LX4,IBLOK),ALINE34(LX4,IBLOK),ALINE35(LX4,IBLOK)
     &,ALINE36(LX4,IBLOK),ALINE37(LX4,IBLOK),ALINE38(LX4,IBLOK)
     &,ALINE39(LX4,IBLOK),ALINE40(LX4,IBLOK),ALINE41(LX4,IBLOK)
     &,ALINE42(LX4,IBLOK),ALINE43(LX4,IBLOK),ALINE44(LX4,IBLOK)
     &,ALINE45(LX4,IBLOK),ALINE46(LX4,IBLOK),ALINE47(LX4,IBLOK)
     &,ALINE48(LX4,IBLOK),ALINE49(LX4,IBLOK),ALINE50(LX4,IBLOK)
     &,ALINE51(LX4,IBLOK),ALINE52(LX4,IBLOK),ALINE53(LX4,IBLOK)
     &,ALINE54(LX4,IBLOK),ALINE55(LX4,IBLOK),ALINE56(LX4,IBLOK)
     &,ALINE57(LX4,IBLOK),ALINE58(LX4,IBLOK),ALINE59(LX4,IBLOK)
     &,ALINE60(LX4,IBLOK),ALINE61(LX4,IBLOK),ALINE62(LX4,IBLOK)
     &,ALINE63(LX4,IBLOK),ALINE64(LX4,IBLOK),ALINE65(LX4,IBLOK)
     &,ALINE66(LX4,IBLOK),ALINE67(LX4,IBLOK),ALINE68(LX4,IBLOK)
     &,ALINE69(LX4,IBLOK),ALINE70(LX4,IBLOK),ALINE71(LX4,IBLOK)
     &,ALINE72(LX4,IBLOK),ALINE73(LX4,IBLOK),ALINE74(LX4,IBLOK)
     &,ALINE75(LX4,IBLOK),ALINE76(LX4,IBLOK),ALINE77(LX4,IBLOK)
     &,ALINE78(LX4,IBLOK),ALINE79(LX4,IBLOK),ALINE80(LX4,IBLOK)
     &,ALINE81(LX4,IBLOK),ALINE82(LX4,IBLOK),ALINE83(LX4,IBLOK)
     &,ALINE84(LX4,IBLOK),ALINE85(LX4,IBLOK),ALINE86(LX4,IBLOK)
     &,ALINE87(LX4,IBLOK),ALINE88(LX4,IBLOK),ALINE89(LX4,IBLOK)
     &,ALINE90(LX4,IBLOK),ALINE91(LX4,IBLOK),ALINE92(LX4,IBLOK)
     &,ALINE93(LX4,IBLOK),ALINE94(LX4,IBLOK),ALINE95(LX4,IBLOK)
     &,ALINE96(LX4,IBLOK),ALINE97(LX4,IBLOK),ALINE98(LX4,IBLOK)
     &,ALINE99(LX4,IBLOK),ALINE100(LX4,IBLOK),ALINE101(LX4,IBLOK)
     &,ALINE102(LX4,IBLOK),ALINE103(LX4,IBLOK),ALINE104(LX4,IBLOK)
     &,ALINE105(LX4,IBLOK),ALINE106(LX4,IBLOK),ALINE107(LX4,IBLOK)
     &,ALINE108(LX4,IBLOK),ALINE109(LX4,IBLOK),ALINE110(LX4,IBLOK)
     &,ALINE111(LX4,IBLOK),ALINE112(LX4,IBLOK),ALINE113(LX4,IBLOK)
     &,ALINE114(LX4,IBLOK),ALINE115(LX4,IBLOK),ALINE116(LX4,IBLOK)
     &,ALINE117(LX4,IBLOK),ALINE118(LX4,IBLOK),ALINE119(LX4,IBLOK)
     &,ALINE120(LX4,IBLOK),ALINE121(LX4,IBLOK),ALINE122(LX4,IBLOK)
     &,ALINE123(LX4,IBLOK),ALINE124(LX4,IBLOK),ALINE125(LX4,IBLOK)
     &,ALINE126(LX4,IBLOK),ALINE127(LX4,IBLOK),ALINE128(LX4,IBLOK)
     &,ALINE129(LX4,IBLOK),ALINE130(LX4,IBLOK),ALINE131(LX4,IBLOK)
     &,ALINE132(LX4,IBLOK),ALINE133(LX4,IBLOK),ALINE134(LX4,IBLOK)
     &,ALINE135(LX4,IBLOK),ALINE136(LX4,IBLOK),ALINE137(LX4,IBLOK)
     &,ALINE138(LX4,IBLOK),ALINE139(LX4,IBLOK),ALINE140(LX4,IBLOK)
     &,ALINE141(LX4,IBLOK),ALINE142(LX4,IBLOK),ALINE143(LX4,IBLOK)
     &,ALINE144(LX4,IBLOK),ALINE145(LX4,IBLOK),ALINE146(LX4,IBLOK)
     &,ALINE147(LX4,IBLOK),ALINE148(LX4,IBLOK),ALINE149(LX4,IBLOK)
     &,ALINE150(LX4,IBLOK),ALINE151(LX4,IBLOK),ALINE152(LX4,IBLOK)
     &,ALINE153(LX4,IBLOK),ALINE154(LX4,IBLOK),ALINE155(LX4,IBLOK)
     &,ALINE156(LX4,IBLOK),ALINE157(LX4,IBLOK),ALINE158(LX4,IBLOK)
     &,ALINE159(LX4,IBLOK),ALINE160(LX4,IBLOK),ALINE161(LX4,IBLOK)
     &,ALINE162(LX4,IBLOK),ALINE163(LX4,IBLOK),ALINE164(LX4,IBLOK)      
     &,ALINE165(LX4,IBLOK),ALINE166(LX4,IBLOK),ALINE167(LX4,IBLOK)
     &,ALINE168(LX4,IBLOK),ALINE169(LX4,IBLOK),ALINE170(LX4,IBLOK)      
     &,ALINE171(LX4,IBLOK),ALINE172(LX4,IBLOK),ALINE173(LX4,IBLOK)      
     &,ALINE174(LX4,IBLOK),ALINE175(LX4,IBLOK),ALINE176(LX4,IBLOK)      
     &,ALINE177(LX4,IBLOK),ALINE178(LX4,IBLOK),ALINE179(LX4,IBLOK)      
     &,ALINE180(LX4,IBLOK),ALINE181(LX4,IBLOK),ALINE182(LX4,IBLOK)      
     &,ALINE183(LX4,IBLOK),ALINE184(LX4,IBLOK),ALINE185(LX4,IBLOK)
     &,ALINE186(LX4,IBLOK),ALINE187(LX4,IBLOK),ALINE188(LX4,IBLOK)            
     &,ALINE189(LX4,IBLOK),ALINE190(LX4,IBLOK),ALINE191(LX4,IBLOK)      
     &,ALINE192(LX4,IBLOK),ALINE193(LX4,IBLOK),ALINE194(LX4,IBLOK)      
     &,ALINE195(LX4,IBLOK),ALINE196(LX4,IBLOK),ALINE197(LX4,IBLOK)      
     &,ALINE198(LX4,IBLOK),ALINE199(LX4,IBLOK),ALINE200(LX4,IBLOK)
     &,ALINE201(LX4,IBLOK),ALINE202(LX4,IBLOK),ALINE203(LX4,IBLOK)
     &,ALINE204(LX4,IBLOK),ALINE205(LX4,IBLOK),ALINE206(LX4,IBLOK)
     &,ALINE207(LX4,IBLOK),ALINE208(LX4,IBLOK),ALINE209(LX4,IBLOK)      
     &,ALINE210(LX4,IBLOK),ALINE211(LX4,IBLOK),ALINE212(LX4,IBLOK)
     &,ALINE213(LX4,IBLOK),ALINE214(LX4,IBLOK),ALINE215(LX4,IBLOK)
     &,ALINE216(LX4,IBLOK),ALINE217(LX4,IBLOK),ALINE218(LX4,IBLOK)      
     &,ALINE219(LX4,IBLOK),ALINE220(LX4,IBLOK),ALINE221(LX4,IBLOK)
     &,ALINE222(LX4,IBLOK),ALINE223(LX4,IBLOK),ALINE224(LX4,IBLOK)
     &,ALINE225(LX4,IBLOK),ALINE226(LX4,IBLOK),ALINE227(LX4,IBLOK)
     &,ALINE228(LX4,IBLOK),ALINE229(LX4,IBLOK),ALINE230(LX4,IBLOK)
     &,ALINE231(LX4,IBLOK),ALINE232(LX4,IBLOK),ALINE233(LX4,IBLOK)
     &,ALINE234(LX4,IBLOK),ALINE235(LX4,IBLOK)
C********************************************************
      COMMON/LINESMOM/ALINEMOM2(LX4,IBLOK,2),ALINEMOM3(LX4,IBLOK,2)
     &,ALINEMOM4(LX4,IBLOK,2),ALINEMOM5(LX4,IBLOK,2)
     &,ALINEMOM6(LX4,IBLOK,2),ALINEMOM7(LX4,IBLOK,2)
     &,ALINEMOM8(LX4,IBLOK,2),ALINEMOM9(LX4,IBLOK,2)
     &,ALINEMOM10(LX4,IBLOK,2),ALINEMOM11(LX4,IBLOK,2)
     &,ALINEMOM12(LX4,IBLOK,2),ALINEMOM13(LX4,IBLOK,2)
     &,ALINEMOM14(LX4,IBLOK,2),ALINEMOM15(LX4,IBLOK,2) 
     &,ALINEMOM16(LX4,IBLOK,2),ALINEMOM17(LX4,IBLOK,2)
     &,ALINEMOM18(LX4,IBLOK,2),ALINEMOM19(LX4,IBLOK,2)
     &,ALINEMOM20(LX4,IBLOK,2),ALINEMOM21(LX4,IBLOK,2)
     &,ALINEMOM22(LX4,IBLOK,2),ALINEMOM23(LX4,IBLOK,2)
     &,ALINEMOM24(LX4,IBLOK,2),ALINEMOM25(LX4,IBLOK,2)
     &,ALINEMOM26(LX4,IBLOK,2),ALINEMOM27(LX4,IBLOK,2)
     &,ALINEMOM28(LX4,IBLOK,2),ALINEMOM29(LX4,IBLOK,2)
     &,ALINEMOM30(LX4,IBLOK,2),ALINEMOM31(LX4,IBLOK,2)
     &,ALINEMOM32(LX4,IBLOK,2),ALINEMOM33(LX4,IBLOK,2)
     &,ALINEMOM34(LX4,IBLOK,2),ALINEMOM35(LX4,IBLOK,2)
     &,ALINEMOM36(LX4,IBLOK,2),ALINEMOM37(LX4,IBLOK,2)
     &,ALINEMOM38(LX4,IBLOK,2),ALINEMOM39(LX4,IBLOK,2)
     &,ALINEMOM40(LX4,IBLOK,2),ALINEMOM41(LX4,IBLOK,2)
     &,ALINEMOM42(LX4,IBLOK,2),ALINEMOM43(LX4,IBLOK,2)
     &,ALINEMOM44(LX4,IBLOK,2),ALINEMOM45(LX4,IBLOK,2)
     &,ALINEMOM46(LX4,IBLOK,2),ALINEMOM47(LX4,IBLOK,2)
     &,ALINEMOM48(LX4,IBLOK,2),ALINEMOM49(LX4,IBLOK,2)
     &,ALINEMOM50(LX4,IBLOK,2),ALINEMOM51(LX4,IBLOK,2)
     &,ALINEMOM52(LX4,IBLOK,2),ALINEMOM53(LX4,IBLOK,2)
     &,ALINEMOM54(LX4,IBLOK,2),ALINEMOM55(LX4,IBLOK,2)
     &,ALINEMOM56(LX4,IBLOK,2),ALINEMOM57(LX4,IBLOK,2)
     &,ALINEMOM58(LX4,IBLOK,2),ALINEMOM59(LX4,IBLOK,2)
     &,ALINEMOM60(LX4,IBLOK,2),ALINEMOM61(LX4,IBLOK,2)
     &,ALINEMOM62(LX4,IBLOK,2),ALINEMOM63(LX4,IBLOK,2)
     &,ALINEMOM64(LX4,IBLOK,2)
     &,ALINEMOM65(LX4,IBLOK,2),ALINEMOM66(LX4,IBLOK,2)
     &,ALINEMOM67(LX4,IBLOK,2),ALINEMOM68(LX4,IBLOK,2)
     &,ALINEMOM69(LX4,IBLOK,2),ALINEMOM70(LX4,IBLOK,2)
     &,ALINEMOM71(LX4,IBLOK,2),ALINEMOM72(LX4,IBLOK,2)
     &,ALINEMOM73(LX4,IBLOK,2),ALINEMOM74(LX4,IBLOK,2) 
     &,ALINEMOM75(LX4,IBLOK,2),ALINEMOM76(LX4,IBLOK,2)
     &,ALINEMOM77(LX4,IBLOK,2),ALINEMOM78(LX4,IBLOK,2)
     &,ALINEMOM79(LX4,IBLOK,2),ALINEMOM80(LX4,IBLOK,2)
     &,ALINEMOM81(LX4,IBLOK,2),ALINEMOM82(LX4,IBLOK,2)
     &,ALINEMOM83(LX4,IBLOK,2),ALINEMOM84(LX4,IBLOK,2) 
     &,ALINEMOM85(LX4,IBLOK,2),ALINEMOM86(LX4,IBLOK,2)
     &,ALINEMOM87(LX4,IBLOK,2),ALINEMOM88(LX4,IBLOK,2)
     &,ALINEMOM89(LX4,IBLOK,2),ALINEMOM90(LX4,IBLOK,2)
     &,ALINEMOM91(LX4,IBLOK,2),ALINEMOM92(LX4,IBLOK,2)
     &,ALINEMOM93(LX4,IBLOK,2),ALINEMOM94(LX4,IBLOK,2) 
     &,ALINEMOM95(LX4,IBLOK,2),ALINEMOM96(LX4,IBLOK,2)
     &,ALINEMOM97(LX4,IBLOK,2),ALINEMOM98(LX4,IBLOK,2)
     &,ALINEMOM99(LX4,IBLOK,2),ALINEMOM100(LX4,IBLOK,2)
     &,ALINEMOM101(LX4,IBLOK,2),ALINEMOM102(LX4,IBLOK,2)
     &,ALINEMOM103(LX4,IBLOK,2),ALINEMOM104(LX4,IBLOK,2) 
     &,ALINEMOM105(LX4,IBLOK,2),ALINEMOM106(LX4,IBLOK,2)
     &,ALINEMOM107(LX4,IBLOK,2),ALINEMOM108(LX4,IBLOK,2)
     &,ALINEMOM109(LX4,IBLOK,2),ALINEMOM110(LX4,IBLOK,2)
     &,ALINEMOM111(LX4,IBLOK,2),ALINEMOM112(LX4,IBLOK,2)
     &,ALINEMOM113(LX4,IBLOK,2),ALINEMOM114(LX4,IBLOK,2) 
     &,ALINEMOM115(LX4,IBLOK,2),ALINEMOM116(LX4,IBLOK,2)
     &,ALINEMOM117(LX4,IBLOK,2),ALINEMOM118(LX4,IBLOK,2)
     &,ALINEMOM119(LX4,IBLOK,2),ALINEMOM120(LX4,IBLOK,2)
     &,ALINEMOM121(LX4,IBLOK,2),ALINEMOM122(LX4,IBLOK,2)
     &,ALINEMOM123(LX4,IBLOK,2),ALINEMOM124(LX4,IBLOK,2) 
     &,ALINEMOM125(LX4,IBLOK,2),ALINEMOM126(LX4,IBLOK,2)
     &,ALINEMOM127(LX4,IBLOK,2),ALINEMOM128(LX4,IBLOK,2)
     &,ALINEMOM129(LX4,IBLOK,2),ALINEMOM130(LX4,IBLOK,2)
     &,ALINEMOM131(LX4,IBLOK,2),ALINEMOM132(LX4,IBLOK,2)
     &,ALINEMOM133(LX4,IBLOK,2),ALINEMOM134(LX4,IBLOK,2) 
     &,ALINEMOM135(LX4,IBLOK,2),ALINEMOM136(LX4,IBLOK,2)
     &,ALINEMOM137(LX4,IBLOK,2),ALINEMOM138(LX4,IBLOK,2)
     &,ALINEMOM139(LX4,IBLOK,2),ALINEMOM140(LX4,IBLOK,2)
     &,ALINEMOM141(LX4,IBLOK,2),ALINEMOM142(LX4,IBLOK,2)
     &,ALINEMOM143(LX4,IBLOK,2),ALINEMOM144(LX4,IBLOK,2) 
     &,ALINEMOM145(LX4,IBLOK,2),ALINEMOM146(LX4,IBLOK,2)
     &,ALINEMOM147(LX4,IBLOK,2),ALINEMOM148(LX4,IBLOK,2)
     &,ALINEMOM149(LX4,IBLOK,2),ALINEMOM150(LX4,IBLOK,2)      
     &,ALINEMOM151(LX4,IBLOK,2),ALINEMOM152(LX4,IBLOK,2)      
     &,ALINEMOM153(LX4,IBLOK,2),ALINEMOM154(LX4,IBLOK,2)      
     &,ALINEMOM155(LX4,IBLOK,2),ALINEMOM156(LX4,IBLOK,2)      
     &,ALINEMOM157(LX4,IBLOK,2),ALINEMOM158(LX4,IBLOK,2)      
     &,ALINEMOM159(LX4,IBLOK,2),ALINEMOM160(LX4,IBLOK,2)
     &,ALINEMOM161(LX4,IBLOK,2),ALINEMOM162(LX4,IBLOK,2)       
     &,ALINEMOM163(LX4,IBLOK,2),ALINEMOM164(LX4,IBLOK,2)       
     &,ALINEMOM165(LX4,IBLOK,2),ALINEMOM166(LX4,IBLOK,2)
     &,ALINEMOM167(LX4,IBLOK,2),ALINEMOM168(LX4,IBLOK,2)      
     &,ALINEMOM169(LX4,IBLOK,2),ALINEMOM170(LX4,IBLOK,2)
     &,ALINEMOM171(LX4,IBLOK,2),ALINEMOM172(LX4,IBLOK,2)
     &,ALINEMOM173(LX4,IBLOK,2),ALINEMOM174(LX4,IBLOK,2)
     &,ALINEMOM175(LX4,IBLOK,2),ALINEMOM176(LX4,IBLOK,2)
     &,ALINEMOM177(LX4,IBLOK,2),ALINEMOM178(LX4,IBLOK,2)
     &,ALINEMOM179(LX4,IBLOK,2),ALINEMOM180(LX4,IBLOK,2)
     &,ALINEMOM181(LX4,IBLOK,2),ALINEMOM182(LX4,IBLOK,2)
     &,ALINEMOM183(LX4,IBLOK,2),ALINEMOM184(LX4,IBLOK,2)
     &,ALINEMOM185(LX4,IBLOK,2),ALINEMOM186(LX4,IBLOK,2)
     &,ALINEMOM187(LX4,IBLOK,2),ALINEMOM188(LX4,IBLOK,2)
     &,ALINEMOM189(LX4,IBLOK,2),ALINEMOM190(LX4,IBLOK,2)
     &,ALINEMOM191(LX4,IBLOK,2),ALINEMOM192(LX4,IBLOK,2)
     &,ALINEMOM193(LX4,IBLOK,2),ALINEMOM194(LX4,IBLOK,2)
     &,ALINEMOM195(LX4,IBLOK,2),ALINEMOM196(LX4,IBLOK,2)
     &,ALINEMOM197(LX4,IBLOK,2),ALINEMOM198(LX4,IBLOK,2)
     &,ALINEMOM199(LX4,IBLOK,2),ALINEMOM200(LX4,IBLOK,2)
     &,ALINEMOM201(LX4,IBLOK,2),ALINEMOM202(LX4,IBLOK,2)
     &,ALINEMOM203(LX4,IBLOK,2),ALINEMOM204(LX4,IBLOK,2) 
     &,ALINEMOM205(LX4,IBLOK,2),ALINEMOM206(LX4,IBLOK,2)
     &,ALINEMOM207(LX4,IBLOK,2),ALINEMOM208(LX4,IBLOK,2)
     &,ALINEMOM209(LX4,IBLOK,2),ALINEMOM210(LX4,IBLOK,2)
     &,ALINEMOM211(LX4,IBLOK,2),ALINEMOM212(LX4,IBLOK,2)
     &,ALINEMOM213(LX4,IBLOK,2),ALINEMOM214(LX4,IBLOK,2) 
     &,ALINEMOM215(LX4,IBLOK,2),ALINEMOM216(LX4,IBLOK,2)
     &,ALINEMOM217(LX4,IBLOK,2),ALINEMOM218(LX4,IBLOK,2)
     &,ALINEMOM219(LX4,IBLOK,2),ALINEMOM220(LX4,IBLOK,2)
     &,ALINEMOM221(LX4,IBLOK,2),ALINEMOM222(LX4,IBLOK,2)
     &,ALINEMOM223(LX4,IBLOK,2),ALINEMOM224(LX4,IBLOK,2) 
     &,ALINEMOM225(LX4,IBLOK,2),ALINEMOM226(LX4,IBLOK,2)
     &,ALINEMOM227(LX4,IBLOK,2),ALINEMOM228(LX4,IBLOK,2)
     &,ALINEMOM229(LX4,IBLOK,2),ALINEMOM230(LX4,IBLOK,2)
     &,ALINEMOM231(LX4,IBLOK,2),ALINEMOM232(LX4,IBLOK,2)
     &,ALINEMOM233(LX4,IBLOK,2),ALINEMOM234(LX4,IBLOK,2) 
     &,ALINEMOM235(LX4,IBLOK,2)      
C********************************************************
      COMMON/ARRAYB/UB11(NCOL2,LSIZEB,3,IBLOK)
      COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
      COMMON/DUMMY/UC11(NCOL2,LSIZEB,3)
     &,IUP(LSIZEB,3),IDN(LSIZEB,3)
      DIMENSION LCNT(IBLOK),LB(IBLOK),IX(3),LS(3)
      DIMENSION ACT(3),AST(3)
      DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2),E11(NCOL2)
     &,F11(NCOL2),UINT11(NCOL2),DUM11(NCOL2)
      DIMENSION REM11(NCOL2)
C*********************************************************************
C     HERE I DEFINE THE SQUARE PULSES
C*********************************************************************
      DIMENSION SQUY1(NCOL2),SQDY1(NCOL2)
      DIMENSION SQUZ1(NCOL2),SQDZ1(NCOL2)
      DIMENSION SQUY2(NCOL2),SQDY2(NCOL2)
      DIMENSION SQUZ2(NCOL2),SQDZ2(NCOL2)
      DIMENSION SQUDY1(NCOL2),SQUDZ1(NCOL2)  !NEW!
      DIMENSION SQDUY1(NCOL2),SQDUZ1(NCOL2)  !NEW!
C*********************************************************************
C     HERE I DEFINE THE WAVE-LIKE SQUARE PULSES
C*********************************************************************
      DIMENSION WSQUY1(NCOL2),WSQUY2(NCOL2),WSQUY3(NCOL2),WSQUY4(NCOL2)
      DIMENSION WSQUZ1(NCOL2),WSQUZ2(NCOL2),WSQUZ3(NCOL2),WSQUZ4(NCOL2)
      DIMENSION WSQDY1(NCOL2),WSQDY2(NCOL2),WSQDY3(NCOL2),WSQDY4(NCOL2)
      DIMENSION WSQDZ1(NCOL2),WSQDZ2(NCOL2),WSQDZ3(NCOL2),WSQDZ4(NCOL2)
C*********************************************************************
C     HERE I DEFINE SOME TT-OPERATORS COMPONENTS
C*********************************************************************
      DIMENSION DUY1(NCOL2),DUZ1(NCOL2),DDY1(NCOL2),DDZ1(NCOL2)
      DIMENSION DUY2(NCOL2),DUZ2(NCOL2),DDY2(NCOL2),DDZ2(NCOL2)
C*********************************************************************
C     HERE I DEFINE SOME PLAQUETTE OPERATOR COMPONENTS
C*********************************************************************
      DIMENSION PQ1(NCOL2),PQ2(NCOL2),PQ3(NCOL2),PQ4(NCOL2)
      DIMENSION PQ5(NCOL2),PQ6(NCOL2),PQ7(NCOL2),PQ8(NCOL2)
C*********************************************************************
      DIMENSION G11(NCOL2)
C*********************************************************************
      DIMENSION WVUY1(NCOL2),WVUY2(NCOL2),WVDY1(NCOL2),WVDY2(NCOL2)
      DIMENSION WVUZ1(NCOL2),WVUZ2(NCOL2),WVDZ1(NCOL2),WVDZ2(NCOL2)
C*********************************************************************
      DIMENSION WANGL1(NCOL2),WANGL2(NCOL2),WANGL3(NCOL2),WANGL4(NCOL2) !new!
      DIMENSION WANGL5(NCOL2),WANGL6(NCOL2),WANGL7(NCOL2),WANGL8(NCOL2) !new!     
C*********************************************************************
      DIMENSION ZIG1W1(NCOL2),ZIG2W1(NCOL2),ZIG1W2(NCOL2),ZIG2W2(NCOL2)
      DIMENSION ZIG1W3(NCOL2),ZIG2W3(NCOL2),ZIG1W4(NCOL2),ZIG2W4(NCOL2)
      DIMENSION TIG1W1(NCOL2),TIG2W1(NCOL2),TIG1W2(NCOL2),TIG2W2(NCOL2)
      DIMENSION TIG1W3(NCOL2),TIG2W3(NCOL2),TIG1W4(NCOL2),TIG2W4(NCOL2)
C*********************************************************************
      DIMENSION LIN0(NCOL2),LIN1(NCOL2),LIN2(NCOL2),LIN4(NCOL2)
      DIMENSION PLQ1(NCOL2),PLQ2(NCOL2),PLQ3(NCOL2),PLQ4(NCOL2)
      DIMENSION PLQ5(NCOL2),PLQ6(NCOL2),PLQ7(NCOL2),PLQ8(NCOL2)
      DIMENSION DPLQ1(NCOL2),DPLQ2(NCOL2),DPLQ3(NCOL2),DPLQ4(NCOL2)
      DIMENSION DPLQ5(NCOL2),DPLQ6(NCOL2),DPLQ7(NCOL2),DPLQ8(NCOL2)
      DIMENSION PL(NCOL2)
C*********************************************************************
      COMPLEX REM11
      COMPLEX UB11,A11,B11,C11,D11,E11,F11,UC11,UINT11,DUM11
      COMPLEX SQUY1,SQDY1,SQUY2,SQDY2,G11
      COMPLEX SQUZ1,SQDZ1,SQUZ2,SQDZ2
      COMPLEX SQUDY1,SQUDZ1,SQDUY1,SQDUZ1  !new!
      COMPLEX WSQUY1,WSQUY2,WSQUY3,WSQUY4
      COMPLEX WSQUZ1,WSQUZ2,WSQUZ3,WSQUZ4
      COMPLEX WSQDY1,WSQDY2,WSQDY3,WSQDY4
      COMPLEX WSQDZ1,WSQDZ2,WSQDZ3,WSQDZ4
      COMPLEX WANGL1,WANGL2,WANGL3,WANGL4 !new!
      COMPLEX WANGL5,WANGL6,WANGL7,WANGL8 !new!     
      COMPLEX WVUY1,WVUY2,WVDY1,WVDY2
      COMPLEX WVUZ1,WVUZ2,WVDZ1,WVDZ2
      COMPLEX ZIG1W1,ZIG2W1,ZIG1W2,ZIG2W2
      COMPLEX ZIG1W3,ZIG2W3,ZIG1W4,ZIG2W4
      COMPLEX TIG1W1,TIG2W1,TIG1W2,TIG2W2
      COMPLEX TIG1W3,TIG2W3,TIG1W4,TIG2W4
      COMPLEX LIN0,LIN1,LIN2,LIN4
c     Plaquette operators on first cite
      COMPLEX PLQ1,PLQ2,PLQ3,PLQ4
      COMPLEX PLQ5,PLQ6,PLQ7,PLQ8      
c
      COMPLEX DPLQ1,DPLQ2,DPLQ3,DPLQ4
      COMPLEX DPLQ5,DPLQ6,DPLQ7,DPLQ8
c      
      COMPLEX PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,PQ8
      COMPLEX DUY1,DUZ1,DDY1,DDZ1,DUY2,DUZ2,DDY2,DDZ2
      COMPLEX PL
C**********************************************************************
      COMPLEX*16 CSUMN,CSUMS(4),CSUM2S(4),CSUM2WS(4)
      COMPLEX*16 CSUMW(4),CSUM2W(4),CSUM3W(4),CSUMUP(4),CSUMUD(4)
      COMPLEX*16 CSUMRC(4),CSUMRCW(4)
      COMPLEX*16 CSUMTT1(8),CSUMTT2(8),CSUMTT3(8),CSUMTT4(8),CSUMTT5(8)
      COMPLEX*16 CSUMTT6(8),CSUMTT7(8),CSUMTT8(8),CSUMTT9(8),CSUMTT10(8)
      COMPLEX*16 CSUMTT11(8),CSUMTT12(16),CSUMTT13(8),CSUMTT14(8)
C
      COMPLEX*16 CSUMPLQ(8),CSUMPLQ2(8),CSUMPLQ3(8),CSUMPLQ4(8)
      COMPLEX*16 CSUMPLQ5(8),CSUMPLQ6(8)
      COMPLEX*16 CSUMPLQ7(16),CSUMPLQ8(16),CSUMPLQ9(16),CSUMPLQ10(16)
      COMPLEX*16 CSUMPLQ11(16),CSUMPLQ12(16),CSUMPLQ13(16),CSUMPLQ14(16)
      COMPLEX*16 CSUMPLQ15(16)
C     
      COMPLEX*16 AKT1,ACT
C**********************************************************************
      COMPLEX*16 ALINE1,ALINE2,ALINE3,ALINE4,ALINE5,ALINE6,ALINE7
      COMPLEX*16 ALINE8,ALINE9,ALINE10,ALINE11,ALINE12,ALINE13
      COMPLEX*16 ALINE14,ALINE15,ALINE16,ALINE17,ALINE18,ALINE19
      COMPLEX*16 ALINE20,ALINE21,ALINE22,ALINE23,ALINE24,ALINE25
      COMPLEX*16 ALINE26,ALINE27,ALINE28,ALINE29,ALINE30
      COMPLEX*16 ALINE31,ALINE32,ALINE33,ALINE34,ALINE35
      COMPLEX*16 ALINE36,ALINE37,ALINE38,ALINE39,ALINE40
      COMPLEX*16 ALINE41,ALINE42,ALINE43,ALINE44,ALINE45
      COMPLEX*16 ALINE46,ALINE47,ALINE48,ALINE49,ALINE50
      COMPLEX*16 ALINE51,ALINE52,ALINE53,ALINE54,ALINE55
      COMPLEX*16 ALINE56,ALINE57,ALINE58,ALINE59,ALINE60
      COMPLEX*16 ALINE61,ALINE62,ALINE63,ALINE64,ALINE65
      COMPLEX*16 ALINE66,ALINE67,ALINE68,ALINE69,ALINE70
      COMPLEX*16 ALINE71,ALINE72,ALINE73,ALINE74,ALINE75
      COMPLEX*16 ALINE76,ALINE77,ALINE78,ALINE79,ALINE80
      COMPLEX*16 ALINE81,ALINE82,ALINE83,ALINE84,ALINE85
      COMPLEX*16 ALINE86,ALINE87,ALINE88,ALINE89,ALINE90
      COMPLEX*16 ALINE91,ALINE92,ALINE93,ALINE94,ALINE95
      COMPLEX*16 ALINE96,ALINE97,ALINE98,ALINE99,ALINE100
      COMPLEX*16 ALINE101,ALINE102,ALINE103,ALINE104,ALINE105
      COMPLEX*16 ALINE106,ALINE107,ALINE108,ALINE109,ALINE110
      COMPLEX*16 ALINE111,ALINE112,ALINE113,ALINE114,ALINE115
      COMPLEX*16 ALINE116,ALINE117,ALINE118,ALINE119,ALINE120
      COMPLEX*16 ALINE121,ALINE122,ALINE123,ALINE124,ALINE125
      COMPLEX*16 ALINE126,ALINE127,ALINE128,ALINE129,ALINE130
      COMPLEX*16 ALINE131,ALINE132,ALINE133,ALINE134,ALINE135
      COMPLEX*16 ALINE136,ALINE137,ALINE138,ALINE139,ALINE140
      COMPLEX*16 ALINE141,ALINE142,ALINE143,ALINE144,ALINE145
      COMPLEX*16 ALINE146,ALINE147,ALINE148,ALINE149,ALINE150
      COMPLEX*16 ALINE151,ALINE152,ALINE153,ALINE154,ALINE155
      COMPLEX*16 ALINE156,ALINE157,ALINE158,ALINE159,ALINE160
      COMPLEX*16 ALINE161,ALINE162,ALINE163,ALINE164,ALINE165
      COMPLEX*16 ALINE166,ALINE167,ALINE168,ALINE169,ALINE170
      COMPLEX*16 ALINE171,ALINE172,ALINE173,ALINE174,ALINE175
      COMPLEX*16 ALINE176,ALINE177,ALINE178,ALINE179,ALINE180
      COMPLEX*16 ALINE181,ALINE182,ALINE183,ALINE184,ALINE185
      COMPLEX*16 ALINE186,ALINE187,ALINE188,ALINE189,ALINE190
      COMPLEX*16 ALINE191,ALINE192,ALINE193,ALINE194,ALINE195
      COMPLEX*16 ALINE196,ALINE197,ALINE198,ALINE199,ALINE200
      COMPLEX*16 ALINE201,ALINE202,ALINE203,ALINE204,ALINE205
      COMPLEX*16 ALINE206,ALINE207,ALINE208,ALINE209,ALINE210
      COMPLEX*16 ALINE211,ALINE212,ALINE213,ALINE214,ALINE215
      COMPLEX*16 ALINE216,ALINE217,ALINE218,ALINE219,ALINE220
      COMPLEX*16 ALINE221,ALINE222,ALINE223,ALINE224,ALINE225
      COMPLEX*16 ALINE226,ALINE227,ALINE228,ALINE229,ALINE230
      COMPLEX*16 ALINE231,ALINE232,ALINE233,ALINE234,ALINE235
C**********************************************************************
      COMPLEX*16 ALINEMOM2,ALINEMOM3,ALINEMOM4
      COMPLEX*16 ALINEMOM5,ALINEMOM6,ALINEMOM7,ALINEMOM8
      COMPLEX*16 ALINEMOM9,ALINEMOM10,ALINEMOM11,ALINEMOM12
      COMPLEX*16 ALINEMOM13,ALINEMOM14,ALINEMOM15,ALINEMOM16
      COMPLEX*16 ALINEMOM17,ALINEMOM18,ALINEMOM19,ALINEMOM20
      COMPLEX*16 ALINEMOM21,ALINEMOM22,ALINEMOM23,ALINEMOM24
      COMPLEX*16 ALINEMOM25,ALINEMOM26,ALINEMOM27,ALINEMOM28
      COMPLEX*16 ALINEMOM29,ALINEMOM30,ALINEMOM31,ALINEMOM32
      COMPLEX*16 ALINEMOM33,ALINEMOM34,ALINEMOM35,ALINEMOM36
      COMPLEX*16 ALINEMOM37,ALINEMOM38,ALINEMOM39,ALINEMOM40
      COMPLEX*16 ALINEMOM41,ALINEMOM42,ALINEMOM43,ALINEMOM44
      COMPLEX*16 ALINEMOM45,ALINEMOM46,ALINEMOM47,ALINEMOM48
      COMPLEX*16 ALINEMOM49,ALINEMOM50,ALINEMOM51,ALINEMOM52
      COMPLEX*16 ALINEMOM53,ALINEMOM54,ALINEMOM55,ALINEMOM56
      COMPLEX*16 ALINEMOM57,ALINEMOM58,ALINEMOM59,ALINEMOM60
      COMPLEX*16 ALINEMOM61,ALINEMOM62,ALINEMOM63,ALINEMOM64
      COMPLEX*16 ALINEMOM65,ALINEMOM66,ALINEMOM67,ALINEMOM68
      COMPLEX*16 ALINEMOM69,ALINEMOM70,ALINEMOM71,ALINEMOM72
      COMPLEX*16 ALINEMOM73,ALINEMOM74,ALINEMOM75,ALINEMOM76
      COMPLEX*16 ALINEMOM77,ALINEMOM78,ALINEMOM79,ALINEMOM80
      COMPLEX*16 ALINEMOM81,ALINEMOM82,ALINEMOM83,ALINEMOM84
      COMPLEX*16 ALINEMOM85,ALINEMOM86,ALINEMOM87,ALINEMOM88
      COMPLEX*16 ALINEMOM89,ALINEMOM90,ALINEMOM91,ALINEMOM92
      COMPLEX*16 ALINEMOM93,ALINEMOM94,ALINEMOM95,ALINEMOM96
      COMPLEX*16 ALINEMOM97,ALINEMOM98,ALINEMOM99,ALINEMOM100 
      COMPLEX*16 ALINEMOM101,ALINEMOM102,ALINEMOM103,ALINEMOM104 
      COMPLEX*16 ALINEMOM105,ALINEMOM106,ALINEMOM107,ALINEMOM108
      COMPLEX*16 ALINEMOM109,ALINEMOM110
      COMPLEX*16 ALINEMOM111,ALINEMOM112,ALINEMOM113,ALINEMOM114
      COMPLEX*16 ALINEMOM115,ALINEMOM116,ALINEMOM117,ALINEMOM118
      COMPLEX*16 ALINEMOM119,ALINEMOM120,ALINEMOM121,ALINEMOM122
      COMPLEX*16 ALINEMOM123,ALINEMOM124,ALINEMOM125,ALINEMOM126
      COMPLEX*16 ALINEMOM127,ALINEMOM128,ALINEMOM129,ALINEMOM130
      COMPLEX*16 ALINEMOM131,ALINEMOM132,ALINEMOM133,ALINEMOM134
      COMPLEX*16 ALINEMOM135,ALINEMOM136,ALINEMOM137,ALINEMOM138
      COMPLEX*16 ALINEMOM139,ALINEMOM140,ALINEMOM141,ALINEMOM142
      COMPLEX*16 ALINEMOM143,ALINEMOM144,ALINEMOM145,ALINEMOM146
      COMPLEX*16 ALINEMOM147,ALINEMOM148,ALINEMOM149,ALINEMOM150
      COMPLEX*16 ALINEMOM151,ALINEMOM152,ALINEMOM153,ALINEMOM154
      COMPLEX*16 ALINEMOM155,ALINEMOM156,ALINEMOM157,ALINEMOM158
      COMPLEX*16 ALINEMOM159,ALINEMOM160,ALINEMOM161,ALINEMOM162
      COMPLEX*16 ALINEMOM163,ALINEMOM164,ALINEMOM165,ALINEMOM166
      COMPLEX*16 ALINEMOM167,ALINEMOM168,ALINEMOM169,ALINEMOM170
      COMPLEX*16 ALINEMOM171,ALINEMOM172,ALINEMOM173,ALINEMOM174
      COMPLEX*16 ALINEMOM175,ALINEMOM176,ALINEMOM177,ALINEMOM178
      COMPLEX*16 ALINEMOM179,ALINEMOM180,ALINEMOM181,ALINEMOM182
      COMPLEX*16 ALINEMOM183,ALINEMOM184,ALINEMOM185,ALINEMOM186
      COMPLEX*16 ALINEMOM187,ALINEMOM188,ALINEMOM189,ALINEMOM190
      COMPLEX*16 ALINEMOM191,ALINEMOM192,ALINEMOM193,ALINEMOM194
      COMPLEX*16 ALINEMOM195,ALINEMOM196,ALINEMOM197,ALINEMOM198
      COMPLEX*16 ALINEMOM199,ALINEMOM200,ALINEMOM201,ALINEMOM202      
      COMPLEX*16 ALINEMOM203,ALINEMOM204,ALINEMOM205,ALINEMOM206
      COMPLEX*16 ALINEMOM207,ALINEMOM208,ALINEMOM209,ALINEMOM210
      COMPLEX*16 ALINEMOM211,ALINEMOM212,ALINEMOM213,ALINEMOM214
      COMPLEX*16 ALINEMOM215,ALINEMOM216,ALINEMOM217,ALINEMOM218
      COMPLEX*16 ALINEMOM219,ALINEMOM220,ALINEMOM221,ALINEMOM222
      COMPLEX*16 ALINEMOM223,ALINEMOM224,ALINEMOM225,ALINEMOM226
      COMPLEX*16 ALINEMOM227,ALINEMOM228,ALINEMOM229,ALINEMOM230
      COMPLEX*16 ALINEMOM231,ALINEMOM232,ALINEMOM233,ALINEMOM234
      COMPLEX*16 ALINEMOM235
C**********************************************************************
      COMPLEX*16 CSUMSMOM(4,2),CSUM2SMOM(4,2),CSUM2WSMOM(4,2)
      COMPLEX*16 CSUMWMOM(4,2),CSUM2WMOM(4,2),CSUM3WMOM(4,2)
      COMPLEX*16 CSUMUPMOM(4,2),CSUMUDMOM(4,2)
C
      COMPLEX*16 CSUMPLQMOM(8,2),CSUMPLQMOM2(8,2),CSUMPLQMOM3(8,2)
      COMPLEX*16 CSUMPLQMOM4(8,2),CSUMPLQMOM5(8,2),CSUMPLQMOM6(8,2)
      COMPLEX*16 CSUMPLQMOM7(16,2), CSUMPLQMOM8(16,2)
      COMPLEX*16 CSUMPLQMOM9(16,2), CSUMPLQMOM10(16,2)
      COMPLEX*16 CSUMPLQMOM11(16,2), CSUMPLQMOM12(16,2)
      COMPLEX*16 CSUMPLQMOM13(16,2), CSUMPLQMOM14(16,2)
      COMPLEX*16 CSUMPLQMOM15(16,2)
C     
      COMPLEX*16 CSUMTTMOM1(8,2),CSUMTTMOM2(8,2),CSUMTTMOM3(8,2)
      COMPLEX*16 CSUMTTMOM4(8,2),CSUMTTMOM5(8,2),CSUMTTMOM6(8,2)
      COMPLEX*16 CSUMTTMOM7(8,2),CSUMTTMOM8(8,2),CSUMTTMOM9(8,2)
      COMPLEX*16 CSUMTTMOM10(8,2),CSUMTTMOM11(8,2)
      COMPLEX*16 CSUMTTMOM12(16,2)
      COMPLEX*16 CSUMTTMOM13(8,2),CSUMTTMOM14(8,2)
C
      COMPLEX*16 GIOT,PF(2)
      COMPLEX UREN11(NCOL2)
      COMPLEX AA(NCOL2)
      COMPLEX DET,CNORM,CDET,DNCOL
C**********************************************************************
      GIOT=(0.0,1.0)
      PI=4.0d0*DATAN(1.0d0)
      LS(1)=LX1
      LS(2)=LX2
      LS(3)=LX3
      LB(1)=1
C*********************************************************************     
      DO 1 IDD=2,IBLOK
         LB(IDD)=2*LB(IDD-1)
 1    CONTINUE
C********************************************************************* 
      DO 2 KU=1,1
         JU=KU+1
         IF(JU.GT.3)JU=JU-3
         IU=KU+2
         IF(IU.GT.3)IU=IU-3
C     
         ID=IBLL
C
         DO 3 IB=1,ID
            LCNT(IB)=0
 3       CONTINUE
C     
         LREST=LS(KU)
C     
         DO 4 IG=1,ID
            IDG=ID-IG+1
            LCNT(IDG)=LREST/LB(IDG)
            LREST=LREST-LCNT(IDG)*LB(IDG)
 4       CONTINUE
C     
         DO 5 IB=1,ID
            IF(LCNT(IB).GE.1)IDS=IB
 5       CONTINUE
C     
         IDSM1=IDS-1
         IF(IDS.EQ.1)IDSM1=IDS
C**********************************************************************
         DO KK=1,3
            DO NN=1,LSIZEB
               IUP(NN,KK)=IUPB(NN,KK,IDS)
               IDN(NN,KK)=IDNB(NN,KK,IDS)
            ENDDO
         ENDDO
C**********************************************************************    
         DO KK=1,3
            DO NN=1,LSIZEB
               DO IC=1,NCOL2
                  UC11(IC,NN,KK)=UB11(IC,NN,KK,IDS)
               ENDDO
            ENDDO
         ENDDO
C**********************************************************************
         CSUMN=(0.0,0.0)
C
         DO IJN=1,4
            CSUMS(IJN)=(0.0,0.0)
            CSUM2S(IJN)=(0.0,0.0)
            CSUM2WS(IJN)=(0.0,0.0)
            CSUMW(IJN)=(0.0,0.0)
            CSUM2W(IJN)=(0.0,0.0)
            CSUM3W(IJN)=(0.0,0.0)
            CSUMUP(IJN)=(0.0,0.0)
            CSUMUD(IJN)=(0.0,0.0)
         ENDDO
C
         DO IJN=1,8
            CSUMTT1(IJN)=(0.0,0.0)
            CSUMTT2(IJN)=(0.0,0.0)
            CSUMTT3(IJN)=(0.0,0.0)
            CSUMTT4(IJN)=(0.0,0.0)
            CSUMTT5(IJN)=(0.0,0.0)
            CSUMTT6(IJN)=(0.0,0.0)
            CSUMTT7(IJN)=(0.0,0.0)
            CSUMTT8(IJN)=(0.0,0.0)
            CSUMTT9(IJN)=(0.0,0.0)
            CSUMTT10(IJN)=(0.0,0.0)
            CSUMTT11(IJN)=(0.0,0.0)
            CSUMTT13(IJN)=(0.0,0.0)
            CSUMTT14(IJN)=(0.0,0.0)
            CSUMPLQ(IJN)=(0.0,0.0)
            CSUMPLQ2(IJN)=(0.0,0.0)
            CSUMPLQ3(IJN)=(0.0,0.0)
            CSUMPLQ4(IJN)=(0.0,0.0)
            CSUMPLQ5(IJN)=(0.0,0.0)
            CSUMPLQ6(IJN)=(0.0,0.0)
         ENDDO
C
         DO IJN=1,16
            CSUMTT12(IJN)=(0.0,0.0)
            CSUMPLQ7(IJN)=(0.0,0.0)
            CSUMPLQ8(IJN)=(0.0,0.0)
            CSUMPLQ9(IJN)=(0.0,0.0)
            CSUMPLQ10(IJN)=(0.0,0.0)
            CSUMPLQ11(IJN)=(0.0,0.0)
            CSUMPLQ12(IJN)=(0.0,0.0)
            CSUMPLQ13(IJN)=(0.0,0.0)
            CSUMPLQ14(IJN)=(0.0,0.0)
            CSUMPLQ15(IJN)=(0.0,0.0)
         ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO IJN=1, 4
            DO IJ=1, 2
               CSUMSMOM(IJN,IJ)=(0.0,0.0)
               CSUM2SMOM(IJN,IJ)=(0.0,0.0)
               CSUM2WSMOM(IJN,IJ)=(0.0,0.0)
               CSUMWMOM(IJN,IJ)=(0.0,0.0)
               CSUM2WMOM(IJN,IJ)=(0.0,0.0)
               CSUM3WMOM(IJN,IJ)=(0.0,0.0)
               CSUMUPMOM(IJN,IJ)=(0.0,0.0)
               CSUMUDMOM(IJN,IJ)=(0.0,0.0)
            ENDDO
         ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO IJN=1, 8
            DO IJ=1, 2
               CSUMTTMOM1(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM2(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM3(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM4(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM5(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM6(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM7(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM8(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM9(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM10(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM11(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM13(IJN,IJ)=(0.0,0.0)
               CSUMTTMOM14(IJN,IJ)=(0.0,0.0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               CSUMPLQMOM(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM2(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM3(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM4(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM5(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM6(IJN,IJ)=(0.0,0.0)
            ENDDO
         ENDDO
C
         DO IJN=1, 16
            DO IJ=1, 2
               CSUMTTMOM12(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM7(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM8(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM9(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM10(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM11(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM12(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM13(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM14(IJN,IJ)=(0.0,0.0)
               CSUMPLQMOM15(IJN,IJ)=(0.0,0.0)
            ENDDO
         ENDDO
C
C**********************************************************************
         NN=0
C**********************************************************************
C     OPENMP INTRINSICS
C     NK,LS,PF,KU
C     UC11 
C     A11, B11, C11, D11
C     SQUY1, SQUZ1, SQDY1, SQDZ1

         
         DO NK=1,LS(KU)
       PF(1)=DCOS((2.0*PI*NK)/LS(KU))+GIOT*DSIN((2.0*PI*NK)/LS(KU))
       PF(2)=DCOS((4.0*PI*NK)/LS(KU))+GIOT*DSIN((4.0*PI*NK)/LS(KU))
c       PF(3)=DCOS((6.0*PI*NK)/LS(KU))+GIOT*DSIN((6.0*PI*NK)/LS(KU))
c       PF(4)=DCOS((8.0*PI*NK)/LS(KU))+GIOT*DSIN((8.0*PI*NK)/LS(KU))
            DO NJ=1,LS(JU)
               DO NI=1,LS(IU)
               NN=NN+1
               IX(IU)=NI
               IX(JU)=NJ
               IX(KU)=NK
               MN=IX(1)+LS(1)*(IX(2)-1)+LS(1)*LS(2)*(IX(3)-1)
C**********************************************************************
               IF(IDS.EQ.ID) THEN
C**********************************************************************
               IREM=3
               IF(IDS.EQ.1) IREM=4
               DO 201 ILOOP=1,IREM
                  IF(ILOOP.EQ.1)LI=LCNT(IDS)
                  IF(ILOOP.GT.1)LI=LCNT(IDS)-2**(ILOOP-2)
                  IF(LI.LT.0) GOTO 201
C
                  M2=MN
                  IF(ILOOP.GT.1) THEN
                     DO I=1, 2**(ILOOP-2)
                        M3=IUP(M2,KU)
                        M2=M3
                     END DO
                  ENDIF
C
                  DO 202 IC=1,NCOL2                        
                     E11(IC)=(0.0,0.0)
 202              CONTINUE
C
                  DO 203 N1=1,NCOL                         
                     IJ=N1+NCOL*(N1-1)                     
                     E11(IJ)=(1.0,0.0)                     
 203              CONTINUE 
C
                  DO 204 NC=1,LI
                     DO 205 IC=1,NCOL2
                        B11(IC)=UC11(IC,M2,KU)
 205                 CONTINUE
                     CALL VMX(1,E11,B11,C11,1)
                     M3=IUP(M2,KU)
                     M2=M3
                     DO 206 IC=1,NCOL2
                        E11(IC)=C11(IC)
 206                 CONTINUE
 204              CONTINUE
C
                  ML=M2
                  IF(ILOOP.EQ.1) THEN
                     DO 207 IC=1,NCOL2
                        LIN0(IC)=E11(IC)
 207                 CONTINUE
                  ENDIF
C********************************************************************                  
                  IF(ILOOP.EQ.2) THEN
                     DO 208 IC=1,NCOL2
                        LIN1(IC)=E11(IC)   
 208                 CONTINUE
                  ENDIF
C********************************************************************                                 
                  IF(ILOOP.EQ.3) THEN
                     DO 209 IC=1,NCOL2
                        LIN2(IC)=E11(IC)   
 209                 CONTINUE
                  ENDIF
C********************************************************************
                  IF(ILOOP.EQ.4) THEN
                     DO 210 IC=1,NCOL2
                        LIN4(IC)=E11(IC)   
 210                 CONTINUE
                  ENDIF
C********************************************************************
 201           CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ELSE
               M2=MN
            ENDIF
C**********************************************************************
C*********************** REMAINING PIECE ******************************
C**********************************************************************            
            DO IC=1, NCOL2
               REM11(IC)=(0.0,0.0)
            ENDDO
C********************************************************************
            DO N1=1,NCOL
               IJ=N1+NCOL*(N1-1)                     
               REM11(IJ)=(1.0,0.0) 
            ENDDO
C********************************************************************
           DO 176 IG=1,IDS
               IDG=IDS-IG+1
               IF(IDG.EQ.ID)GOTO176
               DO 177 NC=1,LCNT(IDG)
                  DO 178 IC=1,NCOL2
                     UINT11(IC)=(0.0,0.0)
 178              CONTINUE
C     
                  DO 179 IJ=1,3
                     IF(IJ.EQ.KU)GOTO179
                     DO 180 IC=1,NCOL2
                        B11(IC)=UB11(IC,M2,IJ,IDSM1)
 180                 CONTINUE
                     M3=IUPB(M2,IJ,IDSM1)
                     DO 181 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDG)
 181                 CONTINUE
                     CALL VMX(1,B11,C11,D11,1)
                     M1=IUPB(M2,KU,IDG)
                     DO 182 IC=1,NCOL2
                        B11(IC)=UB11(IC,M1,IJ,IDSM1)
 182                 CONTINUE
                     CALL HERM(1,B11,DUM11,1)
                     CALL VMX(1,D11,B11,C11,1)
                     DO 183 IC=1,NCOL2
                        UINT11(IC)=UINT11(IC)+C11(IC)
 183                 CONTINUE
                     M3=IDNB(M2,IJ,IDSM1)
                     DO 184 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,IJ,IDSM1)
 184                 CONTINUE
                     DO 185 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDG)
 185                 CONTINUE
                     CALL HERM(1,B11,DUM11,1)
                     CALL VMX(1,B11,C11,D11,1)
                     M1=IUPB(M3,KU,IDG)
                     DO 186 IC=1,NCOL2
                        B11(IC)=UB11(IC,M1,IJ,IDSM1)
 186                 CONTINUE
                     CALL VMX(1,D11,B11,C11,1)
                     DO 187 IC=1,NCOL2
                        UINT11(IC)=UINT11(IC)+C11(IC)
 187                 CONTINUE
 179              CONTINUE
                  DO 188 IC=1,NCOL2
                     B11(IC)=UB11(IC,M2,KU,IDG)
 188              CONTINUE
                  DO 189 IC=1,NCOL2
                     UINT11(IC)=UINT11(IC)+B11(IC)
 189              CONTINUE
Cnew
                  CALL RENORMBS(UINT11,UREN11)
                  DO 441 IJ=1,NCOL2
                     AA(IJ)=UREN11(IJ)
 441              CONTINUE
                  JMAT=NCOL
                  CALL DETNANT(JMAT,DET,AA)
c
                  DNCOL=CMPLX(1.0/NCOL)
                  CDET=DET**DNCOL
                  CNORM=1.0/CDET
c
                  DO 644 IC=1,NCOL2
                     B11(IC)=UREN11(IC)*CNORM
 644               CONTINUE
c********************************************
                  CALL VMX(1,REM11,B11,C11,1)
                  M1=M2
                  DO 191 IC=1,NCOL2
                     REM11(IC)=C11(IC)
 191              CONTINUE
                  M2=IUPB(M1,KU,IDG)
 177           CONTINUE
 176        CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               DO 9 IDDD=1,337 !new!
C
                  M2=MN
C
                  DO 10 IJ=1,NCOL2
                     A11(IJ)=(0.0,0.0)
 10               CONTINUE
C
                  DO 11 N1=1,NCOL
                     IJ=N1+NCOL*(N1-1)
                     A11(IJ)=(1.0,0.0)
 11               CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  IF (IDS.EQ.ID) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               ICO=2
               IF(IDDD.LT.5) ICO=1
               IF(((IDDD.GT.12).AND.(IDDD.LT.17)).AND.(IDS.NE.1)) ICO=1
               IF(((IDDD.GT.16).AND.(IDDD.LT.21)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.20).AND.(IDDD.LT.25)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.24).AND.(IDDD.LT.29)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.28).AND.(IDDD.LT.33)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.40).AND.(IDDD.LT.49)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.48).AND.(IDDD.LT.57)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.56).AND.(IDDD.LT.65)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.64).AND.(IDDD.LT.69)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.68).AND.(IDDD.LT.73)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.72).AND.(IDDD.LT.81)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.80).AND.(IDDD.LT.89)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.88).AND.(IDDD.LT.97)).AND.(IDS.EQ.1)) ICO=4
               IF(((IDDD.GT.96).AND.(IDDD.LT.105)).AND.(IDS.EQ.1)) ICO=4
              IF(((IDDD.GT.104).AND.(IDDD.LT.113)).AND.(IDS.EQ.1)) ICO=4
              IF(((IDDD.GT.112).AND.(IDDD.LT.129)).AND.(IDS.EQ.1)) ICO=4
              IF(((IDDD.GT.128).AND.(IDDD.LT.137)).AND.(IDS.EQ.1)) ICO=4
              IF(((IDDD.GT.136).AND.(IDDD.LT.145)).AND.(IDS.NE.1)) ICO=1
               IF(IDDD.EQ.145) ICO=1
               IF((IDDD.GT.145).AND.(IDDD.LT.194)) ICO=1
               IF((IDDD.GT.209).AND.(IDDD.LT.338)) ICO=1 ! THIS NEEDS TO BE FIXED  !
C**********************************************************************C
                     IF(LCNT(IDS).GE.ICO) THEN
C**********************************************************************C
C                     UP SQUARE PULSE                                  C 
C**********************************************************************C
C                     UP Y
C**********************************************************************
                        IF(IDDD.EQ.1) THEN
C
                           DO 12 IC=1,NCOL2
                              B11(IC)=UC11(IC,M2,JU)
 12                        CONTINUE
                           M3=IUP(M2,JU)
                           DO 13 IC=1,NCOL2
                              C11(IC)=UC11(IC,M3,KU)
 13                        CONTINUE
                           CALL VMX(1,B11,C11,D11,1)
                           M3=IUP(M2,KU)
                           M2=M3
                           DO 14 IC=1,NCOL2
                              C11(IC)=UC11(IC,M2,JU)
 14                        CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,D11,C11,SQUY1,1)
                           CALL VMX(1,A11,SQUY1,C11,1)
C     
                           DO 15 IC=1,NCOL2
                              A11(IC)=C11(IC)
 15                        CONTINUE
C     
                           IEEE=1
C     
                        ENDIF
C**********************************************************************C
C     UP Z                                                             C 
C**********************************************************************C
                        IF(IDDD.EQ.2) THEN
C     
                           DO 16 IC=1,NCOL2
                              B11(IC)=UC11(IC,M2,IU)
 16                        CONTINUE
                           M3=IUP(M2,IU)
                           DO 17 IC=1,NCOL2
                              C11(IC)=UC11(IC,M3,KU)
 17                        CONTINUE
                           CALL VMX(1,B11,C11,D11,1)
                           M3=IUP(M2,KU)
                           M2=M3
                           DO 18 IC=1,NCOL2
                              C11(IC)=UC11(IC,M2,IU)
 18                        CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,D11,C11,SQUZ1,1)
                           CALL VMX(1,A11,SQUZ1,C11,1)
C     
                           DO 19 IC=1,NCOL2
                              A11(IC)=C11(IC)
 19                        CONTINUE
C     
                           IEEE=2
C     
                        ENDIF
C**********************************************************************C 
C                       DOWN Y                                         C
C**********************************************************************C
                        IF(IDDD.EQ.3) THEN
C     
                           M3=IDN(M2,JU)
                           DO 20 IC=1,NCOL2
                              D11(IC)=UC11(IC,M3,JU)
 20                        CONTINUE
                           CALL HERM(1,D11,DUM11,1)
                           DO 21 IC=1,NCOL2
                              C11(IC)=UC11(IC,M3,KU)
 21                        CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M4=IUP(M3,KU)
                           DO 22 IC=1,NCOL2
                              D11(IC)=UC11(IC,M4,JU)
 22                        CONTINUE
                           CALL VMX(1,B11,D11,SQDY1,1)
                           CALL VMX(1,A11,SQDY1,C11,1)
C     
                           DO 23 IC=1,NCOL2
                              A11(IC)=C11(IC)
 23                        CONTINUE
C     
                           IEEE=3
C     
                        ENDIF
C**********************************************************************C 
C                       DOWN Z                                 C
C**********************************************************************C
                        IF(IDDD.EQ.4) THEN
C     
                           M3=IDN(M2,IU)
                           DO 24 IC=1,NCOL2
                              D11(IC)=UC11(IC,M3,IU)
 24                        CONTINUE
                           CALL HERM(1,D11,DUM11,1)
                           DO 25 IC=1,NCOL2
                              C11(IC)=UC11(IC,M3,KU)
 25                        CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M4=IUP(M3,KU)
                           DO 26 IC=1,NCOL2
                              D11(IC)=UC11(IC,M4,IU)
 26                        CONTINUE
                           CALL VMX(1,B11,D11,SQDZ1,1)
                           CALL VMX(1,A11,SQDZ1,C11,1)
C     
                           DO 27 IC=1,NCOL2
                              A11(IC)=C11(IC)
 27                        CONTINUE
C     
                           IEEE=4
C     
                        ENDIF
C**********************************************************************C
C                       UP - UP SQUARE PULSES                         *C 
C**********************************************************************C
C                       UP Y
C**********************************************************************C
                        IF(IDDD.EQ.5) THEN
C     
                           M3=IUP(M2,KU)
                           DO 28 IC=1,NCOL2
                              D11(IC)=UC11(IC,M3,JU)
 28                        CONTINUE
                           M4=IUP(M3,JU)
                           DO 29 IC=1,NCOL2
                              C11(IC)=UC11(IC,M4,KU)
 29                        CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M2=IUP(M3,KU)
                           DO 30 IC=1,NCOL2
                              C11(IC)=UC11(IC,M2,JU)
 30                        CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,B11,C11,SQUY2,1)
                           CALL VMX(1,SQUY1,SQUY2,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C**********************************************************************C
C              UP Z
C**********************************************************************C
                        IF(IDDD.EQ.6) THEN
C     
                           M3=IUP(M2,KU)
                           DO 31 IC=1,NCOL2
                              D11(IC)=UC11(IC,M3,IU)
 31                        CONTINUE
                           M4=IUP(M3,IU)
                           DO 32 IC=1,NCOL2
                              C11(IC)=UC11(IC,M4,KU)
 32                        CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M2=IUP(M3,KU)
                           DO 33 IC=1,NCOL2
                              C11(IC)=UC11(IC,M2,IU)
 33                        CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,B11,C11,SQUZ2,1)
                           CALL VMX(1,SQUZ1,SQUZ2,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C**********************************************************************C 
C              DOWN Y                                              C
C**********************************************************************C
                        IF(IDDD.EQ.7) THEN
C     
                           M3=IUP(M2,KU)
                           M4=IDN(M3,JU)
                           DO 34 IC=1,NCOL2
                              D11(IC)=UC11(IC,M4,JU)
 34                        CONTINUE
                           CALL HERM(1,D11,DUM11,1)
                           DO 35 IC=1,NCOL2
                              C11(IC)=UC11(IC,M4,KU)
 35                        CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M1=IUP(M4,KU)
                           DO 36 IC=1,NCOL2
                              C11(IC)=UC11(IC,M1,JU)
 36                        CONTINUE
                           CALL VMX(1,B11,C11,SQDY2,1)
                           CALL VMX(1,SQDY1,SQDY2,A11,1)
                           M2=IUP(M3,KU)
C     
                           IEEE=3
C     
                        ENDIF
C**********************************************************************C 
C              DOWN Z                                              C
C**********************************************************************C
                        IF(IDDD.EQ.8) THEN
C     
                           M3=IUP(M2,KU)
                           M4=IDN(M3,IU)
                           DO 37 IC=1,NCOL2
                              D11(IC)=UC11(IC,M4,IU)
 37                        CONTINUE
                           CALL HERM(1,D11,DUM11,1)
                           DO 38 IC=1,NCOL2
                              C11(IC)=UC11(IC,M4,KU)
 38                        CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M1=IUP(M4,KU)
                           DO 39 IC=1,NCOL2
                              C11(IC)=UC11(IC,M1,IU)
 39                        CONTINUE
                           CALL VMX(1,B11,C11,SQDZ2,1)
                           CALL VMX(1,SQDZ1,SQDZ2,A11,1)
                           M2=IUP(M3,KU)
C     
                           IEEE=4
C     
                        ENDIF
C**********************************************************************C
C**********************************************************************C
C              UP - DOWN SQUARE PULSES                                 C 
C**********************************************************************C
C**********************************************************************C
C                 UP Y 
C**********************************************************************C
                        IF(IDDD.EQ.9) THEN
                           CALL VMX(1,SQUY1,SQDY2,SQUDY1,1)
                           DO IC=1, NCOL2
                              A11(IC)=SQUDY1(IC) !NEW!
                           ENDDO
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=1
C     
                        ENDIF
C**********************************************************************C
C                 UP Z
C**********************************************************************C
                        IF(IDDD.EQ.10) THEN    
                           CALL VMX(1,SQUZ1,SQDZ2,SQUDZ1,1)
                           DO IC=1, NCOL2
                              A11(IC)=SQUDZ1(IC) !NEW!
                           ENDDO
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=2
C     
                        ENDIF
C**********************************************************************C 
C                 DOWN Y
C**********************************************************************C
                        IF(IDDD.EQ.11) THEN
                           CALL VMX(1,SQDY1,SQUY2,SQDUY1,1)
                           DO IC=1, NCOL2
                              A11(IC)=SQDUY1(IC) !NEW!
                           ENDDO
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=3
C     
                        ENDIF
C**********************************************************************C 
C                 DOWN Z
C**********************************************************************C
                        IF(IDDD.EQ.12) THEN
                           CALL VMX(1,SQDZ1,SQUZ2,SQDUZ1,1)
                           DO IC=1, NCOL2
                              A11(IC)=SQDUZ1(IC) !NEW!
                           ENDDO
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C
                           IEEE=4
C     
                        ENDIF
C**********************************************************************C
C**********************************************************************C
C                 UP WAVE - LIKE PULSE                                *C
C**********************************************************************C
C**********************************************************************C
C                 UP Y
C**********************************************************************C
                  IF(IDDD.EQ.13) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C     
                     DO 40 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,JU,IDSW)
 40                  CONTINUE
                     M3=IUPB(M2,JU,IDSW)
                     DO 41 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 41                  CONTINUE
                     CALL VMX(1,D11,C11,B11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 42 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,JU,IDSW)
 42                  CONTINUE
                     CALL HERM(1,C11,DUM11,1)
                     CALL VMX(1,B11,C11,WSQUY1,1)
C                     
                     M3=IDNB(M2,JU,IDSW)
                     DO 43 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,JU,IDSW)
 43                  CONTINUE
                     CALL HERM(1,B11,DUM11,1)
                     DO 44 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 44                  CONTINUE
                     CALL VMX(1,B11,C11,D11,1)
                     M4=IUPB(M3,KU,IDSW)
                     DO 45 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,JU,IDSW)
 45                  CONTINUE
                     CALL VMX(1,D11,B11,WSQDY2,1)
                     CALL VMX(1,WSQUY1,WSQDY2,WVUY1,1)
                     DO 46 IC=1,NCOL2
                        A11(IC)=WVUY1(IC)
 46                  CONTINUE
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=1
C
                  ENDIF
C**********************************************************************C
C                 UP Z
C**********************************************************************C
                  IF(IDDD.EQ.14) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C     
                     DO 47 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,IU,IDSW)
 47                  CONTINUE
                     M3=IUPB(M2,IU,IDSW)
                     DO 48 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 48                  CONTINUE
                     CALL VMX(1,D11,C11,B11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 49 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,IU,IDSW)
 49                  CONTINUE
                     CALL HERM(1,C11,DUM11,1)
                     CALL VMX(1,B11,C11,WSQUZ1,1)
C                     
                     M3=IDNB(M2,IU,IDSW)
                     DO 50 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,IU,IDSW)
 50                  CONTINUE
                     CALL HERM(1,B11,DUM11,1)
                     DO 51 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 51                  CONTINUE
                     CALL VMX(1,B11,C11,D11,1)
                     M4=IUPB(M3,KU,IDSW)
                     DO 52 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,IU,IDSW)
 52                  CONTINUE
                     CALL VMX(1,D11,B11,WSQDZ2,1)
                     CALL VMX(1,WSQUZ1,WSQDZ2,WVUZ1,1)
                     DO 53 IC=1,NCOL2
                        A11(IC)=WVUZ1(IC)
 53                  CONTINUE
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=2
C
                  ENDIF
C**********************************************************************C
C                 DOWN Y                                              *C 
C**********************************************************************C
                  IF(IDDD.EQ.15) THEN
C     
                        IDSW=IDS-1
C
                        IF(IDSW.EQ.0) THEN
                           IDSW=1
                        ENDIF
C     
                        M3=IDNB(M2,JU,IDSW)
                        DO 54 IC=1,NCOL2
                           D11(IC)=UB11(IC,M3,JU,IDSW)
 54                     CONTINUE
                        CALL HERM(1,D11,DUM11,1)
                        DO 55 IC=1,NCOL2
                           C11(IC)=UB11(IC,M3,KU,IDSW)
 55                     CONTINUE
                        CALL VMX(1,D11,C11,B11,1)
                        M4=IUPB(M3,KU,IDSW)
                        DO 56 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,JU,IDSW)
 56                     CONTINUE
                        CALL VMX(1,B11,C11,WSQDY1,1)
                        M3=IUPB(M2,KU,IDSW)
                        M2=M3
C     
                        DO 57 IC=1,NCOL2
                           B11(IC)=UB11(IC,M2,JU,IDSW)
 57                     CONTINUE
                        M3=IUPB(M2,JU,IDSW)
                        DO 58 IC=1,NCOL2
                           C11(IC)=UB11(IC,M3,KU,IDSW)
 58                     CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IUPB(M2,KU,IDSW)
                        M2=M3
                        DO 59 IC=1,NCOL2
                           C11(IC)=UB11(IC,M2,JU,IDSW)
 59                     CONTINUE
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,D11,C11,WSQUY2,1)
                        CALL VMX(1,WSQDY1,WSQUY2,WVDY1,1)
C     
                        DO 60 IC=1, NCOL2
                           A11(IC)=WVDY1(IC)
 60                     CONTINUE
C
                        IEEE=3
C 
                  ENDIF
C**********************************************************************C
C                 DOWN Z                                              *C 
C**********************************************************************C
                  IF(IDDD.EQ.16) THEN
C     
                        IDSW=IDS-1
C
                        IF(IDSW.EQ.0) THEN
                           IDSW=1
                        ENDIF
C     
                        M3=IDNB(M2,IU,IDSW)
                        DO 61 IC=1,NCOL2
                           D11(IC)=UB11(IC,M3,IU,IDSW)
 61                     CONTINUE
                        CALL HERM(1,D11,DUM11,1)
                        DO 62 IC=1,NCOL2
                           C11(IC)=UB11(IC,M3,KU,IDSW)
 62                     CONTINUE
                        CALL VMX(1,D11,C11,B11,1)
                        M4=IUPB(M3,KU,IDSW)
                        DO 63 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,IU,IDSW)
 63                     CONTINUE
                        CALL VMX(1,B11,C11,WSQDZ1,1)
                        M3=IUPB(M2,KU,IDSW)
                        M2=M3
C     
                        DO 64 IC=1,NCOL2
                           B11(IC)=UB11(IC,M2,IU,IDSW)
 64                     CONTINUE
                        M3=IUPB(M2,IU,IDSW)
                        DO 65 IC=1,NCOL2
                           C11(IC)=UB11(IC,M3,KU,IDSW)
 65                     CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                        M3=IUPB(M2,KU,IDSW)
                        M2=M3
                        DO 66 IC=1,NCOL2
                           C11(IC)=UB11(IC,M2,IU,IDSW)
 66                     CONTINUE
                        CALL HERM(1,C11,DUM11,1)
                        CALL VMX(1,D11,C11,WSQUZ2,1)
                        CALL VMX(1,WSQDZ1,WSQUZ2,WVDZ1,1)
C     
                        DO 67 IC=1, NCOL2
                           A11(IC)=WVDZ1(IC)
 67                     CONTINUE
C
                        IEEE=4
C     
                  ENDIF
C**********************************************************************C
C**********************************************************************C
C              UP - UP WAVE-LIKE PULSE                                 C
C**********************************************************************C
C**********************************************************************C
C              UP Y
C**********************************************************************C
                  IF(IDDD.EQ.17) THEN
C
                     IDSW=IDS-1

                     IF(IDSW.EQ.0) THEN 
                        IDSW=1
                     ENDIF
C 
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C 
                     DO 68 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,JU,IDSW)
 68                  CONTINUE 
                     M3=IUPB(M2,JU,IDSW)
                     DO 69 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 69                  CONTINUE
                     CALL VMX(1,D11,C11,B11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 70 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,JU,IDSW)
 70                  CONTINUE
                     CALL HERM(1,C11,DUM11,1)
                     CALL VMX(1,B11,C11,WSQUY3,1)
C
                     M3=IDNB(M2,JU,IDSW)
                     DO 71 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,JU,IDSW)
 71                  CONTINUE
                     CALL HERM(1,B11,DUM11,1)
                     DO 72 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 72                  CONTINUE
                     CALL VMX(1,B11,C11,D11,1)
                     M4=IUPB(M3,KU,IDSW)
                     DO 73 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,JU,IDSW)
 73                  CONTINUE
                     CALL VMX(1,D11,B11,WSQDY4,1)
                     CALL VMX(1,WSQUY3,WSQDY4,WVUY2,1)
                     CALL VMX(1,WVUY1,WVUY2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=1
C     
                  ENDIF
C**********************************************************************C
C                 UP Z
C**********************************************************************C
                  IF(IDDD.EQ.18) THEN
C
                     IDSW=IDS-1

                     IF(IDSW.EQ.0) THEN 
                        IDSW=1
                     ENDIF
C 
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C 
                     DO 74 IC=1,NCOL2
                        D11(IC)=UB11(IC,M2,IU,IDSW)
 74                  CONTINUE 
                     M3=IUPB(M2,IU,IDSW)
                     DO 75 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 75                  CONTINUE
                     CALL VMX(1,D11,C11,B11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 76 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,IU,IDSW)
 76                  CONTINUE
                     CALL HERM(1,C11,DUM11,1)
                     CALL VMX(1,B11,C11,WSQUZ3,1)
C
                     M3=IDNB(M2,IU,IDSW)
                     DO 77 IC=1,NCOL2
                        B11(IC)=UB11(IC,M3,IU,IDSW)
 77                  CONTINUE
                     CALL HERM(1,B11,DUM11,1)
                     DO 78 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 78                  CONTINUE
                     CALL VMX(1,B11,C11,D11,1)
                     M4=IUPB(M3,KU,IDSW)
                     DO 79 IC=1,NCOL2
                        B11(IC)=UB11(IC,M4,IU,IDSW)
 79                  CONTINUE
                     CALL VMX(1,D11,B11,WSQDZ4,1)
                     CALL VMX(1,WSQUZ3,WSQDZ4,WVUZ2,1)
                     CALL VMX(1,WVUZ1,WVUZ2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=2
C     
                  ENDIF
C**********************************************************************C
C                 DOWN Y                                               C
C**********************************************************************C
                  IF(IDDD.EQ.19) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C     
                     M3=IDNB(M2,JU,IDSW)
                     DO 80 IC=1,NCOL2
                        D11(IC)=UB11(IC,M3,JU,IDSW)
 80                  CONTINUE
                     CALL HERM(1,D11,DUM11,1)
                     DO 81 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 81                  CONTINUE
                     CALL VMX(1,D11,C11,B11,1)
                     M4=IUPB(M3,KU,IDSW)
                     DO 82 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,JU,IDSW)
 82                  CONTINUE
                     CALL VMX(1,B11,C11,WSQDY3,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 83 IC=1,NCOL2
                        B11(IC)=UB11(IC,M2,JU,IDSW)
 83                  CONTINUE
                     M3=IUPB(M2,JU,IDSW)
                     DO 84 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 84                  CONTINUE
                     CALL VMX(1,B11,C11,D11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 85 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,JU,IDSW)
 85                  CONTINUE
                     CALL HERM(1,C11,DUM11,1)
                     CALL VMX(1,D11,C11,WSQUY4,1)
                     CALL VMX(1,WSQDY3,WSQUY4,WVDY2,1)
                     CALL VMX(1,WVDY1,WVDY2,A11,1)
C
                        IEEE=3
C     
                  ENDIF
C**********************************************************************C
C                 DOWN Y                                               C
C**********************************************************************C
                  IF(IDDD.EQ.20) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C     
                     M3=IDNB(M2,IU,IDSW)
                     DO 86 IC=1,NCOL2
                        D11(IC)=UB11(IC,M3,IU,IDSW)
 86                  CONTINUE
                     CALL HERM(1,D11,DUM11,1)
                     DO 87 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 87                  CONTINUE
                     CALL VMX(1,D11,C11,B11,1)
                     M4=IUPB(M3,KU,IDSW)
                     DO 88 IC=1,NCOL2
                        C11(IC)=UB11(IC,M4,IU,IDSW)
 88                  CONTINUE
                     CALL VMX(1,B11,C11,WSQDZ3,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 89 IC=1,NCOL2
                        B11(IC)=UB11(IC,M2,IU,IDSW)
 89                  CONTINUE
                     M3=IUPB(M2,IU,IDSW)
                     DO 90 IC=1,NCOL2
                        C11(IC)=UB11(IC,M3,KU,IDSW)
 90                  CONTINUE
                     CALL VMX(1,B11,C11,D11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     DO 91 IC=1,NCOL2
                        C11(IC)=UB11(IC,M2,IU,IDSW)
 91                  CONTINUE
                     CALL HERM(1,C11,DUM11,1)
                     CALL VMX(1,D11,C11,WSQUZ4,1)
                     CALL VMX(1,WSQDZ3,WSQUZ4,WVDZ2,1)
                     CALL VMX(1,WVDZ1,WVDZ2,A11,1)
C
                        IEEE=4
C     
                  ENDIF
C**********************************************************************C
C**********************************************************************C
C                 UP - DOWN WAVE-LIKE PULSE                           *C
C**********************************************************************C
C**********************************************************************C
C                 UP Y                                                 C
C**********************************************************************C
                  IF(IDDD.EQ.21) THEN
C     
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C     
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVUY1,WVDY2,A11,1)
C
                        IEEE=1
C
                  ENDIF
C**********************************************************************C
C                 UP Z                                                 C
C**********************************************************************C
                  IF(IDDD.EQ.22) THEN
C     
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C     
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVUZ1,WVDZ2,A11,1)
C
                        IEEE=2
C
                  ENDIF
C**********************************************************************C
C                 DOWN Y                                              *C
C**********************************************************************C
                  IF(IDDD.EQ.23) THEN
C     
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C     
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVDY1,WVUY2,A11,1)
C
                        IEEE=3
C     
                  ENDIF
C**********************************************************************C
C                 DOWN Z                                              *C
C**********************************************************************C
                  IF(IDDD.EQ.24) THEN
C     
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C     
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVDZ1,WVUZ2,A11,1)
C
                        IEEE=4
C     
                  ENDIF
C**********************************************************************C
C**********************************************************************C
C                   /\_____/\  UP PULSE                                C
C**********************************************************************C
C**********************************************************************C
C                  UP Y
C**********************************************************************C
                   IF(IDDD.EQ.25) THEN
C
                      IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                        M3=IUPB(M2,KU,IDSW)
                        DO 139 IC=1,NCOL2
                           B11(IC)=UB11(IC,M3,KU,IDSW)
 139                    CONTINUE
                        M4=IUPB(M3,KU,IDSW)
                        DO 140 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,KU,IDSW)
 140                    CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                     ELSE
                        M3=IUPB(M2,KU,IDSW)
                        DO 141 IC=1, NCOL2
                           D11(IC)=UB11(IC,M3,KU,IDSW+1)
 141                    CONTINUE
                     ENDIF
C
                     CALL VMX(1,WSQUY1,D11,C11,1)
                     CALL VMX(1,C11,WSQUY4,A11,1)
C
C
                        IEEE=1
C
                  ENDIF
C**********************************************************************C
C                  UP Z
C**********************************************************************C
                   IF(IDDD.EQ.26) THEN
C
                      IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                        M3=IUPB(M2,KU,IDSW)
                        DO 142 IC=1,NCOL2
                           B11(IC)=UB11(IC,M3,KU,IDSW)
 142                    CONTINUE
                        M4=IUPB(M3,KU,IDSW)
                        DO 143 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,KU,IDSW)
 143                    CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                     ELSE
                        M3=IUPB(M2,KU,IDSW)
                        DO 144 IC=1, NCOL2
                           D11(IC)=UB11(IC,M3,KU,IDSW+1)
 144                    CONTINUE
                     ENDIF
C
                     CALL VMX(1,WSQUZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQUZ4,A11,1)
C
                        IEEE=2
C
                  ENDIF
C**********************************************************************C
C                   DOWN Y                                             C
C**********************************************************************C
                   IF(IDDD.EQ.27) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 145 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 145                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 146 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 146                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 147 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 147                     CONTINUE
                      ENDIF
C     
                     CALL VMX(1,WSQDY1,D11,C11,1)
                     CALL VMX(1,C11,WSQDY4,A11,1)
C
                        IEEE=3
C
                  ENDIF
C**********************************************************************C
C                   DOWN Z                                             C
C**********************************************************************C
                   IF(IDDD.EQ.28) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 149 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 149                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 150 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 150                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 151 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 151                     CONTINUE
                      ENDIF
C     
                     CALL VMX(1,WSQDZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQDZ4,A11,1)
C
                        IEEE=4
C
                  ENDIF
C**********************************************************************C
C**********************************************************************C
C                   /\-----\/  PULSE                                   C
C**********************************************************************C
C**********************************************************************C
C                  UP Y
C**********************************************************************C
                   IF(IDDD.EQ.29) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 152 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 152                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 153 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 153                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 154 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 154                     CONTINUE
                      ENDIF
C
                     CALL VMX(1,WSQUY1,D11,C11,1)
                     CALL VMX(1,C11,WSQDY4,A11,1)
C
                        IEEE=1
C
                  ENDIF
C**********************************************************************C
C                  UP Z
C**********************************************************************C
                   IF(IDDD.EQ.30) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 155 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 155                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 156 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 156                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 157 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 157                     CONTINUE
                      ENDIF
C
                     CALL VMX(1,WSQUZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQDZ4,A11,1)
C
                        IEEE=2
C
                  ENDIF
C**********************************************************************C
C                 DOWN Y                                               C
C**********************************************************************C
                  IF(IDDD.EQ.31) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                        M3=IUPB(M2,KU,IDSW)
                        DO 158 IC=1,NCOL2
                           B11(IC)=UB11(IC,M3,KU,IDSW)
 158                    CONTINUE
                        M4=IUPB(M3,KU,IDSW)
                        DO 159 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,KU,IDSW)  
 159                    CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                     ELSE
                        M3=IUPB(M2,KU,IDSW)
                        DO 160 IC=1, NCOL2
                           D11(IC)=UB11(IC,M3,KU,IDSW+1)
 160                    CONTINUE
                     ENDIF
C     
                     CALL VMX(1,WSQDY1,D11,C11,1)
                     CALL VMX(1,C11,WSQUY4,A11,1)
C
                        IEEE=3
C     
                  ENDIF
C**********************************************************************C
C                 DOWN Z                                               C
C**********************************************************************C
                  IF(IDDD.EQ.32) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                        M3=IUPB(M2,KU,IDSW)
                        DO 161 IC=1,NCOL2
                           B11(IC)=UB11(IC,M3,KU,IDSW)
 161                    CONTINUE
                        M4=IUPB(M3,KU,IDSW)
                        DO 162 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,KU,IDSW)  
 162                    CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                     ELSE
                        M3=IUPB(M2,KU,IDSW)
                        DO 262 IC=1, NCOL2
                           D11(IC)=UB11(IC,M3,KU,IDSW+1)
 262                    CONTINUE
                     ENDIF
C     
                     CALL VMX(1,WSQDZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQUZ4,A11,1)
C
                        IEEE=4
C     
                  ENDIF
C**********************************************************************
C     TT-1 OPERATORS
C**********************************************************************
                        IF(IDDD.EQ.33)THEN
c     
                           CALL VMX(1,SQUY1,SQUZ2,A11,1)
C
                           DO IC=1, NCOL2
                              WANGL1(IC)=A11(IC) !NEW!
                           ENDDO
C                           
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=1
C     
                        ENDIF
C**********************************************************************
                        IF(IDDD.EQ.34) THEN
C     
                           CALL VMX(1,SQUZ1,SQDY2,A11,1)
C
                           DO IC=1, NCOL2
                              WANGL2(IC)=A11(IC) !NEW!
                           ENDDO
C                                                      
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=2
C     
                        ENDIF
C**********************************************************************
                        IF(IDDD.EQ.35) THEN
C     
                           CALL VMX(1,SQDY1,SQDZ2,A11,1)
C    
                           DO IC=1, NCOL2
                              WANGL3(IC)=A11(IC) !NEW!
                           ENDDO
C                           
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=3
C     
                        ENDIF
C**********************************************************************
                        IF(IDDD.EQ.36) THEN
C     
                           CALL VMX(1,SQDZ1,SQUY2,A11,1)
C                           
                           DO IC=1, NCOL2
                              WANGL4(IC)=A11(IC) !NEW!
                           ENDDO                           
C
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=4
C     
                        ENDIF
C**********************************************************************
                        IF(IDDD.EQ.37) THEN
C     
                           CALL VMX(1,SQUY1,SQDZ2,A11,1)
C                           
                           DO IC=1, NCOL2
                              WANGL5(IC)=A11(IC) !NEW!
                           ENDDO                           
C
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=5
C     
                        ENDIF
C**********************************************************************
                        IF(IDDD.EQ.38) THEN
C     
                           CALL VMX(1,SQUZ1,SQUY2,A11,1)
C                           
                           DO IC=1, NCOL2
                              WANGL6(IC)=A11(IC) !NEW!
                           ENDDO                           
C
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=6
C     
                        ENDIF
C**********************************************************************
                        IF(IDDD.EQ.39) THEN
C     
                           CALL VMX(1,SQDY1,SQUZ2,A11,1)
C                           
                           DO IC=1, NCOL2
                              WANGL7(IC)=A11(IC) !NEW!
                           ENDDO                           
C                           
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           IEEE=7
C     
                        ENDIF
C**********************************************************************
                        IF(IDDD.EQ.40) THEN
                           M3=IUP(M2,KU)
                           M2=IUP(M3,KU)
C     
                           CALL VMX(1,SQDZ1,SQDY2,A11,1)
C                           
                           DO IC=1, NCOL2
                              WANGL8(IC)=A11(IC) !NEW!
                           ENDDO                                                      
C     
                           IEEE=8
C     
                        ENDIF
C**********************************************************************
C     TT-2 OPERATORS
C**********************************************************************
                  IF(IDDD.EQ.41) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVUY1,WVUZ2,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.42) THEN
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVUZ1,WVDY2,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.43) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVDY1,WVDZ2,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.44) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVDZ1,WVUY2,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.45) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVUZ1,WVUY2,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.46) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVDY1,WVUZ2,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.47) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVDZ1,WVDY2,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.48) THEN
C
                     IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
C
                     M3=IUPB(M2,KU,IDSW)
                     M4=IUPB(M3,KU,IDSW)
                     M1=IUPB(M4,KU,IDSW)
                     M2=IUPB(M1,KU,IDSW)
C
                     CALL VMX(1,WVUY1,WVDZ2,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
C**********************************************************************C
C                   /\----/_/  UP PULSE-TT3                            C
C**********************************************************************C
C**********************************************************************C
C                  1
C**********************************************************************C
                   IF(IDDD.EQ.49) THEN
C
                      IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                        M3=IUPB(M2,KU,IDSW)
                        DO 301 IC=1,NCOL2
                           B11(IC)=UB11(IC,M3,KU,IDSW)
 301                    CONTINUE
                        M4=IUPB(M3,KU,IDSW)
                        DO 302 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,KU,IDSW)
 302                    CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                     ELSE
                        M3=IUPB(M2,KU,IDSW)
                        DO 303 IC=1, NCOL2
                           D11(IC)=UB11(IC,M3,KU,IDSW+1)
 303                    CONTINUE
                     ENDIF
C                    !   WARNING WE CAN SPEED UP THE COMPUTATION  !
                     CALL VMX(1,WSQUY1,D11,C11,1)
                     CALL VMX(1,C11,WSQUZ4,A11,1)
C
C
                        IEEE=1
C
                  ENDIF
C**********************************************************************C
C                  2
C**********************************************************************C
                   IF(IDDD.EQ.50) THEN
C
                      IDSW=IDS-1
C
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                        M3=IUPB(M2,KU,IDSW)
                        DO 304 IC=1,NCOL2
                           B11(IC)=UB11(IC,M3,KU,IDSW)
 304                    CONTINUE
                        M4=IUPB(M3,KU,IDSW)
                        DO 305 IC=1,NCOL2
                           C11(IC)=UB11(IC,M4,KU,IDSW)
 305                    CONTINUE
                        CALL VMX(1,B11,C11,D11,1)
                     ELSE
                        M3=IUPB(M2,KU,IDSW)
                        DO 306 IC=1, NCOL2
                           D11(IC)=UB11(IC,M3,KU,IDSW+1)
 306                    CONTINUE
                     ENDIF
C
                     CALL VMX(1,WSQUZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQDY4,A11,1)
C
                        IEEE=2
C
                  ENDIF
C**********************************************************************C
C                   3                                             C
C**********************************************************************C
                   IF(IDDD.EQ.51) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 307 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 307                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 308 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 308                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 309 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 309                     CONTINUE
                      ENDIF
C     
                     CALL VMX(1,WSQDY1,D11,C11,1)
                     CALL VMX(1,C11,WSQDZ4,A11,1)
C
                        IEEE=3
C
                  ENDIF
C**********************************************************************C
C                 4                                                    C
C**********************************************************************C
                   IF(IDDD.EQ.52) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 310 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 310                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 311 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 311                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 312 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 312                     CONTINUE
                      ENDIF
C     
                     CALL VMX(1,WSQDZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQUY4,A11,1)
C
                        IEEE=4
C
                  ENDIF
C**********************************************************************C
C                 5                                                    C
C**********************************************************************C
                   IF(IDDD.EQ.53) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 313 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 313                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 314 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 314                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 315 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 315                     CONTINUE
                      ENDIF
C
                     CALL VMX(1,WSQUY1,D11,C11,1)
                     CALL VMX(1,C11,WSQDZ4,A11,1)     
C
                        IEEE=5
C
                  ENDIF
C**********************************************************************C
C                 6                                                    C
C**********************************************************************C
                   IF(IDDD.EQ.54) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 316 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 316                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 317 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 317                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 318 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 318                     CONTINUE
                      ENDIF
C
                     CALL VMX(1,WSQUZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQUY4,A11,1)     
C
                        IEEE=6
C
                  ENDIF
C**********************************************************************C
C                 7                                                    C
C**********************************************************************C
                   IF(IDDD.EQ.55) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 319 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 319                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 320 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 320                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 321 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 321                     CONTINUE
                      ENDIF
C
                     CALL VMX(1,WSQDY1,D11,C11,1)
                     CALL VMX(1,C11,WSQUZ4,A11,1)     
C
                        IEEE=7
C
                  ENDIF
C**********************************************************************C
C                 8                                                    C
C**********************************************************************C
                   IF(IDDD.EQ.56) THEN
C
                      IDSW=IDS-1
C
                      IF(IDSW.EQ.0) THEN
                         IDSW=1
                         M3=IUPB(M2,KU,IDSW)
                         DO 322 IC=1,NCOL2
                            B11(IC)=UB11(IC,M3,KU,IDSW)
 322                     CONTINUE
                         M4=IUPB(M3,KU,IDSW)
                         DO 323 IC=1,NCOL2
                            C11(IC)=UB11(IC,M4,KU,IDSW)  
 323                     CONTINUE
                         CALL VMX(1,B11,C11,D11,1)
                      ELSE
                         M3=IUPB(M2,KU,IDSW)
                         DO 324 IC=1, NCOL2
                            D11(IC)=UB11(IC,M3,KU,IDSW+1)
 324                     CONTINUE
                      ENDIF
C
                     CALL VMX(1,WSQDZ1,D11,C11,1)
                     CALL VMX(1,C11,WSQDY4,A11,1)     
C
                        IEEE=8
C
                  ENDIF
**********************************************************************
C                 T-T4 OPERATORS
C**********************************************************************
                  IF(IDDD.EQ.57) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,ZIG1W1,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,ZIG2W1,1)
                     CALL VMX(1,ZIG1W1,ZIG2W1,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.58) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,ZIG1W2,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,ZIG2W2,1)
                     CALL VMX(1,ZIG1W2,ZIG2W2,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.59) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,ZIG1W3,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,ZIG2W3,1)
                     CALL VMX(1,ZIG1W3,ZIG2W3,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.60) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,ZIG1W4,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,ZIG2W4,1)
                     CALL VMX(1,ZIG1W4,ZIG2W4,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.61) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,TIG1W4,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,TIG2W4,1)
                     CALL VMX(1,TIG1W4,TIG2W4,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.62) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,TIG1W1,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,TIG2W1,1)
                     CALL VMX(1,TIG1W1,TIG2W1,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.63) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,TIG1W2,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,TIG2W2,1)
                     CALL VMX(1,TIG1W2,TIG2W2,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************
                  IF(IDDD.EQ.64) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,TIG1W3,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,TIG2W3,1)
                     CALL VMX(1,TIG1W3,TIG2W3,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************
C                 T-T5 OPERATORS
C**********************************************************************
                  IF(IDDD.EQ.65) THEN
C
                     CALL VMX(1,WSQUY1,WSQUY2,DUY1,1)
                     CALL VMX(1,WSQUY3,WSQUY4,DUY2,1)
                     CALL VMX(1,DUY1,DUY2,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.66) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUZ2,DUZ1,1)
                     CALL VMX(1,WSQUZ3,WSQUZ4,DUZ2,1)
                     CALL VMX(1,DUZ1,DUZ2,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.67) THEN
C
                     CALL VMX(1,WSQDY1,WSQDY2,DDY1,1)
                     CALL VMX(1,WSQDY3,WSQDY4,DDY2,1)
                     CALL VMX(1,DDY1,DDY2,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.68) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDZ2,DDZ1,1)
                     CALL VMX(1,WSQDZ3,WSQDZ4,DDZ2,1)
                     CALL VMX(1,DDZ1,DDZ2,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
C                 TT-6 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.69) THEN
C
                     CALL VMX(1,DUY1,DDY2,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.70) THEN
C
                     CALL VMX(1,DUZ1,DDZ2,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.71) THEN
C
                     CALL VMX(1,DDY1,DUY2,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.72) THEN
C
                     CALL VMX(1,DDZ1,DUZ2,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
C                 T-T 7 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.73) THEN
C
                     CALL VMX(1,DUY1,DUZ2,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.74) THEN
C
                     CALL VMX(1,DUZ1,DDY2,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.75) THEN
C
                     CALL VMX(1,DDY1,DDZ2,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.76) THEN
C
                     CALL VMX(1,DDZ1,DUY2,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.77) THEN
C
                     CALL VMX(1,DDZ1,DDY2,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.78) THEN
C
                     CALL VMX(1,DUY1,DDZ2,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.79) THEN
C
                     CALL VMX(1,DUZ1,DUY2,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.80) THEN
C
                     CALL VMX(1,DDY1,DUZ2,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
C                 T-T 8 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.81) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.82) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.83) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.84) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.85) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.86) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.87) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.88) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
C                 T-T 9 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.89) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.90) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.91) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.92) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.93) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.94) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.95) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.96) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
C                 T-T 10 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.97) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.98) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.99) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.100) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.101) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.102) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.103) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.104) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
C                 T-T 11 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.105) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.106) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.107) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.108) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.109) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.110) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.111) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.112) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
C                 T-T 12 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.113) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.114) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.115) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.116) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.117) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.118) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.119) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.120) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.121) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=9
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.122) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=10
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.123) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=11
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.124) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=12
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.125) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=13
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.126) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=14
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.127) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=15
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.128) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=16
C
                  ENDIF
C**********************************************************************C
C                 T-T 13 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.129) THEN
C
                     CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=1
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.130) THEN
C
                     CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=2
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.131) THEN
C
                     CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=3
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.132) THEN
C
                     CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.133) THEN
C
                     CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
                     CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=5
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.134) THEN
C
                     CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
                     CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=6
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.135) THEN
C
                     CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
                     CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=7
C
                  ENDIF
C**********************************************************************c
                  IF(IDDD.EQ.136) THEN
C
                     CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
                     CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
                     CALL VMX(1,B11,C11,A11,1)
C
                     IEEE=8
C
                  ENDIF
C**********************************************************************C
C**********************************************************************C
C                 T-T 14 OPERATORS
C**********************************************************************C
                  IF(IDDD.EQ.137) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQUY1,WSQUZ2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=1
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.138) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQUZ1,WSQDY2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=2
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.139) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQDY1,WSQDZ2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=3
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.140) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQDZ1,WSQUY2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=4
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.141) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQUY1,WSQDZ2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=5
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.142) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQUZ1,WSQUY2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=6
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.143) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQDY1,WSQUZ2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=7
C
                  ENDIF
C**********************************************************************C
                  IF(IDDD.EQ.144) THEN
C     
                     IDSW=IDS-1
C     
                     IF(IDSW.EQ.0) THEN
                        IDSW=1
                     ENDIF
                     CALL VMX(1,WSQDZ1,WSQDY2,A11,1)
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
                     M3=IUPB(M2,KU,IDSW)
                     M2=M3
C
                        IEEE=8
C
                  ENDIF
C**********************************************************************C
C                         NORMAL POLYAKOV LOOP                         C
C**********************************************************************C
                        IF (IDDD.EQ.145) THEN
                           DO 167 IC=1, NCOL2
                              A11(IC)=UC11(IC,M2,KU)
                              PL(IC)=UC11(IC,M2,KU)
 167                       CONTINUE
                        ENDIF
C**********************************************************************C
C***                      PLAQUETTE OPERATOR                        ***C
C**********************************************************************C
C***                      UP Y OPERATOR
C**********************************************************************C
                        IF (IDDD.EQ.146) THEN
                           DO 401 IC=1, NCOL2
                              B11(IC)=UC11(IC,M2,JU)
 401                       CONTINUE
                           M3=IUP(M2,JU)
                           DO 402 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,IU)
 402                       CONTINUE
                           CALL VMX(1,B11,C11,D11,1)
                           M3=IUP(M2,IU)
                           DO 403 IC=1, NCOL2
                              B11(IC)=UC11(IC,M3,JU)
 403                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           CALL VMX(1,D11,B11,C11,1)
                           DO 404 IC=1,NCOL2
                              D11(IC)=UC11(IC,M2,IU)
 404                       CONTINUE
                           CALL HERM(1,D11,DUM11,1)
                           CALL VMX(1,C11,D11,PLQ1,1)
C     
                           CALL VMX(1,PLQ1,PL,PQ1,1)
C
                           DO 505 IC=1,NCOL2
                              A11(IC)=PQ1(IC)
 505                       CONTINUE
C
                           IEEE=1
C     
                        ENDIF
C**********************************************************************C
C***                      UP Z OPERATOR
C**********************************************************************C
                        IF(IDDD.EQ.147) THEN
C     
                           DO 405 IC=1, NCOL2
                              B11(IC)=UC11(IC,M2,IU)
 405                       CONTINUE
                           M3=IUP(M2,IU)
                           M1=IDN(M3,JU)
                           DO 406 IC=1, NCOL2
                              C11(IC)=UC11(IC,M1,JU)
 406                       CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,B11,C11,D11,1)
                           M3=IDN(M2,JU)
                           DO 407 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,IU)
 407                       CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,D11,C11,B11,1)
                           DO 408 IC=1, NCOL2
                              D11(IC)=UC11(IC,M3,JU)
 408                       CONTINUE
                           CALL VMX(1,B11,D11,PLQ2,1)
C     
                           CALL VMX(1,PLQ2,PL,PQ2,1)
C     
                           DO 509 IC=1,NCOL2
                              A11(IC)=PQ2(IC)
 509                       CONTINUE
C
                           IEEE=2
C     
                        ENDIF
C**********************************************************************C
C***                      DOWN Y OPERATOR
C**********************************************************************C
                        IF(IDDD.EQ.148) THEN
C     
                           M3=IDN(M2,JU)
                           DO 409 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,JU)
 409                       CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           M1=IDN(M3,IU)
                           DO 410 IC=1, NCOL2
                              B11(IC)=UC11(IC,M1,IU)
 410                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           CALL VMX(1,C11,B11,D11,1)
                           DO 411 IC=1, NCOL2
                              C11(IC)=UC11(IC,M1,JU)
 411                       CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M3=IDN(M2,IU)
                           DO 412 IC=1, NCOL2
                              D11(IC)=UC11(IC,M3,IU)
 412                       CONTINUE
                           CALL VMX(1,B11,D11,PLQ3,1)
C     
                           CALL VMX(1,PLQ3,PL,PQ3,1)
C
                           DO 513 IC=1,NCOL2
                              A11(IC)=PQ3(IC)
 513                       CONTINUE
C     
                           IEEE=3
C     
                        ENDIF
C**********************************************************************C
C***                      DOWN Z OPERATOR
C**********************************************************************C
                        IF(IDDD.EQ.149) THEN
C     
                           M3=IDN(M2,IU)
                           DO 413 IC=1, NCOL2
                              B11(IC)=UC11(IC,M3,IU)
 413                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           DO 414 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,JU)
 414                       CONTINUE
                           CALL VMX(1,B11,C11,D11,1)
                           M1=IUP(M3,JU)
                           DO 415 IC=1, NCOL2
                              B11(IC)=UC11(IC,M1,IU)
 415                       CONTINUE
                           CALL VMX(1,D11,B11,C11,1)
                           DO 416 IC=1, NCOL2
                              B11(IC)=UC11(IC,M2,JU)
 416                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           CALL VMX(1,C11,B11,PLQ4,1)
C     
                           CALL VMX(1,PLQ4,PL,PQ4,1)
C
                           DO 517 IC=1,NCOL2
                              A11(IC)=PQ4(IC)
 517                       CONTINUE
C
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
C***                       UP Y OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.150)THEN
C
                           DO 801 IC=1, NCOL2
                              PLQ5(IC)=PLQ2(IC)
 801                       CONTINUE
C
                           CALL HERM(1,PLQ5,DUM11,1)
                           CALL VMX(1,PLQ5,PL,PQ5,1)
C
                           DO 518 IC=1,NCOL2
                              A11(IC)=PQ5(IC)
 518                       CONTINUE
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
C***                       UP Z OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.151)THEN
C     
                           DO 802 IC=1, NCOL2
                              PLQ6(IC)=PLQ1(IC)
 802                       CONTINUE
C
                           CALL HERM(1,PLQ6,DUM11,1)
                           CALL VMX(1,PLQ6,PL,PQ6,1)
C
                           DO 519 IC=1,NCOL2
                              A11(IC)=PQ6(IC)
 519                       CONTINUE
C
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
C***                       DOWN Y OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.152)THEN
C     
                           DO 803 IC=1, NCOL2
                              PLQ7(IC)=PLQ4(IC)
 803                       CONTINUE
C
                           CALL HERM(1,PLQ7,DUM11,1)
                           CALL VMX(1,PLQ7,PL,PQ7,1)
C     
                           DO 520 IC=1,NCOL2
                              A11(IC)=PQ7(IC)
 520                       CONTINUE
C
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
C***                       DOWN Z OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.153)THEN
C     
                           DO 804 IC=1, NCOL2
                              PLQ8(IC)=PLQ3(IC)
 804                       CONTINUE
C
                           CALL HERM(1,PLQ8,DUM11,1)
                           CALL VMX(1,PLQ8,PL,PQ8,1)
C
                           DO 521 IC=1,NCOL2
                              A11(IC)=PQ8(IC)
 521                       CONTINUE
C
                           IEEE=8
C
                        ENDIF
C**********************************************************************C
C***                      PLAQUETTE OPERATOR 2                      ***C
C**********************************************************************C
C***                      UP Y OPERATOR
C**********************************************************************C
                        IF (IDDD.EQ.154) THEN
                           M5=IUP(M2,KU)
                           DO 601 IC=1, NCOL2
                              B11(IC)=UC11(IC,M5,JU)
 601                       CONTINUE
                           M3=IUP(M5,JU)
                           DO 602 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,IU)
 602                       CONTINUE
                           CALL VMX(1,B11,C11,D11,1)
                           M3=IUP(M5,IU)
                           DO 603 IC=1, NCOL2
                              B11(IC)=UC11(IC,M3,JU)
 603                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           CALL VMX(1,D11,B11,C11,1)
                           DO 604 IC=1,NCOL2
                              D11(IC)=UC11(IC,M5,IU)
 604                       CONTINUE
                           CALL HERM(1,D11,DUM11,1)
                           CALL VMX(1,C11,D11,DPLQ1,1)
C     
                           CALL VMX(1,PQ1,DPLQ1,A11,1)
C
                           IEEE=1
C     
                        ENDIF
C**********************************************************************C
C***                      UP Z OPERATOR
C**********************************************************************C
                        IF(IDDD.EQ.155) THEN
C     
                           M5=IUP(M2,KU)
                           DO 605 IC=1, NCOL2
                              B11(IC)=UC11(IC,M5,IU)
 605                       CONTINUE
                           M3=IUP(M5,IU)
                           M1=IDN(M3,JU)
                           DO 606 IC=1, NCOL2
                              C11(IC)=UC11(IC,M1,JU)
 606                       CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,B11,C11,D11,1)
                           M3=IDN(M5,JU)
                           DO 607 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,IU)
 607                       CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,D11,C11,B11,1)
                           DO 608 IC=1, NCOL2
                              D11(IC)=UC11(IC,M3,JU)
 608                       CONTINUE
                           CALL VMX(1,B11,D11,DPLQ2,1)
C     
                           CALL VMX(1,PQ2,DPLQ2,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C**********************************************************************C
C***                      DOWN Y OPERATOR
C**********************************************************************C
                        IF(IDDD.EQ.156) THEN
C     
                           M5=IUP(M2,KU)
                           M3=IDN(M5,JU)
                           DO 609 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,JU)
 609                       CONTINUE
                           CALL HERM(1,C11,DUM11,1)
                           M1=IDN(M3,IU)
                           DO 610 IC=1, NCOL2
                              B11(IC)=UC11(IC,M1,IU)
 610                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           CALL VMX(1,C11,B11,D11,1)
                           DO 611 IC=1, NCOL2
                              C11(IC)=UC11(IC,M1,JU)
 611                       CONTINUE
                           CALL VMX(1,D11,C11,B11,1)
                           M3=IDN(M5,IU)
                           DO 612 IC=1, NCOL2
                              D11(IC)=UC11(IC,M3,IU)
 612                       CONTINUE
                           CALL VMX(1,B11,D11,DPLQ3,1)
C     
                           CALL VMX(1,PQ3,DPLQ3,A11,1)
C     
                           IEEE=3
C     
                        ENDIF
C**********************************************************************C
C***                      DOWN Z OPERATOR
C**********************************************************************C
                        IF(IDDD.EQ.157) THEN
C     
                           M5=IUP(M2,KU)
                           M3=IDN(M5,IU)
                           DO 613 IC=1, NCOL2
                              B11(IC)=UC11(IC,M3,IU)
 613                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           DO 614 IC=1, NCOL2
                              C11(IC)=UC11(IC,M3,JU)
 614                       CONTINUE
                           CALL VMX(1,B11,C11,D11,1)
                           M1=IUP(M3,JU)
                           DO 615 IC=1, NCOL2
                              B11(IC)=UC11(IC,M1,IU)
 615                       CONTINUE
                           CALL VMX(1,D11,B11,C11,1)
                           DO 616 IC=1, NCOL2
                              B11(IC)=UC11(IC,M5,JU)
 616                       CONTINUE
                           CALL HERM(1,B11,DUM11,1)
                           CALL VMX(1,C11,B11,DPLQ4,1)
C     
                           CALL VMX(1,PQ4,DPLQ4,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
C***                       UP Y OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.158)THEN
C
                           DO 701 IC=1, NCOL2
                              DPLQ5(IC)=DPLQ2(IC)
 701                       CONTINUE
C     
                           CALL HERM(1,DPLQ5,DUM11,1)
                           CALL VMX(1,PQ5,DPLQ5,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
C***                       UP Z OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.159)THEN
C     
                           DO 702 IC=1, NCOL2
                              DPLQ6(IC)=DPLQ1(IC)
 702                       CONTINUE
C
                           CALL HERM(1,DPLQ6,DUM11,1)
                           CALL VMX(1,PQ6,DPLQ6,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
C***  DOWN Y OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.160)THEN
C     
                           DO 703 IC=1, NCOL2
                             DPLQ7(IC)=DPLQ4(IC)
 703                       CONTINUE
C
                           CALL HERM(1,DPLQ7,DUM11,1)
                           CALL VMX(1,PQ7,DPLQ7,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
C***                       DOWN Z OPERATOR +
C******************************************************************C
                        IF(IDDD.EQ.161)THEN
C     
                           DO 704 IC=1, NCOL2
                              DPLQ8(IC)=DPLQ3(IC)
 704                       CONTINUE
C
                           CALL HERM(1,DPLQ8,DUM11,1)
                           CALL VMX(1,PQ8,DPLQ8,A11,1)
C
                           IEEE=8
C
                        ENDIF
C******************************************************************C
C***             PLAQUETTE OPERATPORS 3
C******************************************************************C
                        IF(IDDD.EQ.162)THEN
C
C                           DO 705 IC=1, NCOL2
C                              B11(IC)=DPLQ1(IC)
C 705                       CONTINUE
C
C                           CALL HERM(1,B11,DUM11,1)
C                           CALL VMX(1,PQ1,B11,A11,1)
                            CALL VMX(1,PQ1,DPLQ6,A11,1) ! New altered !                           
C     
                           IEEE=1
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.163)THEN
C
c                           DO 706 IC=1, NCOL2
c                              B11(IC)=DPLQ2(IC)
c 706                       CONTINUE
cC
c                          CALL HERM(1,B11,DUM11,1)
c                           CALL VMX(1,PQ2,B11,A11,1)                           
                           CALL VMX(1,PQ2,DPLQ5,A11,1) ! New altered !                           
C     
                           IEEE=2
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.164)THEN
C
c                           DO 707 IC=1, NCOL2
c                              B11(IC)=DPLQ3(IC)
c 707                       CONTINUE
cC
c                           CALL HERM(1,B11,DUM11,1)
c                           CALL VMX(1,PQ3,B11,A11,1)
                           
                            CALL VMX(1,PQ3,DPLQ8,A11,1) ! New altered !
C     
                           IEEE=3
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.165)THEN
C
c                           DO 708 IC=1, NCOL2
c                              B11(IC)=DPLQ4(IC)
c 708                       CONTINUE
cC
c                          CALL HERM(1,B11,DUM11,1)
c                           CALL VMX(1,PQ4,B11,A11,1)                           
                           CALL VMX(1,PQ4,DPLQ7,A11,1) ! New altered ! 
C     
                           IEEE=4
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.166)THEN
C
                           CALL VMX(1,PQ5,DPLQ2,A11,1)
C     
                           IEEE=5
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.167)THEN
C
                           CALL VMX(1,PQ6,DPLQ1,A11,1)
C     
                           IEEE=6
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.168)THEN
C
                           CALL VMX(1,PQ7,DPLQ4,A11,1)
C     
                           IEEE=7
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.169)THEN
C
                           CALL VMX(1,PQ8,DPLQ3,A11,1)
C     
                           IEEE=8
C
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 4 
C******************************************************************C
                        IF(IDDD.EQ.170)THEN
C
c                           DO 709 IC=1, NCOL2
c                              C11(IC)=DPLQ1(IC)
c 709                       CONTINUE
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ1,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ6,A11,1)
C
                           IEEE=1
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.171)THEN
C
c                           DO 710 IC=1, NCOL2
c                              C11(IC)=DPLQ2(IC)
c 710                       CONTINUE
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ2,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ5,A11,1)
C
                           IEEE=2
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.172)THEN
C
c                           DO 711 IC=1, NCOL2
c                              C11(IC)=DPLQ3(IC)
c 711                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ3,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ8,A11,1)
C
                           IEEE=3
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.173)THEN
C
c                           DO 712 IC=1, NCOL2
c                              C11(IC)=DPLQ4(IC)
c 712                       CONTINUE
C
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ4,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ7,A11,1)
C
                           IEEE=4
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.174)THEN
C
c                           DO 713 IC=1, NCOL2
c                              C11(IC)=PLQ2(IC)
c 713                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ5,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ2,A11,1)
C
                           IEEE=5
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.175)THEN
C
c                           DO 714 IC=1, NCOL2
c                              C11(IC)=PLQ1(IC)
c 714                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ6,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ1,A11,1)
C
                           IEEE=6
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.176)THEN
C
c                           DO 715 IC=1, NCOL2
c                              C11(IC)=PLQ4(IC)
c 715                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ7,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ4,A11,1)
C
                           IEEE=7
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.177)THEN
C
c                           DO 716 IC=1, NCOL2
c                              C11(IC)=PLQ3(IC)
c 716                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PLQ8,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ3,A11,1)
C
                           IEEE=8
C
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 5
C******************************************************************C
                        IF(IDDD.EQ.178)THEN
C
c                           DO 717 IC=1, NCOL2
c                              C11(IC)=DPLQ3(IC)
c 717                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ1,DPLQ8,A11,1)
C
                           IEEE=1
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.179)THEN
C
c                           DO 718 IC=1, NCOL2
c                              C11(IC)=DPLQ4(IC)
c 718                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ2,DPLQ7,A11,1)
C
                           IEEE=2
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.180)THEN
C
c                           DO 719 IC=1, NCOL2
c                              C11(IC)=DPLQ1(IC)
c 719                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ3,DPLQ6,A11,1)
C
                           IEEE=3
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.181)THEN
C
c                           DO 720 IC=1, NCOL2
c                              C11(IC)=DPLQ2(IC)
c 720                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ4,DPLQ5,A11,1)
C
                           IEEE=4
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.182)THEN
C
                           CALL VMX(1,PQ5,DPLQ4,A11,1)
C
                           IEEE=5
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.183)THEN
C
                           CALL VMX(1,PQ6,DPLQ3,A11,1)
C
                           IEEE=6
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.184)THEN
C
                           CALL VMX(1,PQ7,DPLQ2,A11,1)
C
                           IEEE=7
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.185)THEN
C
                           CALL VMX(1,PQ8,DPLQ1,A11,1)
C
                           IEEE=8
C
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 6
C******************************************************************C
                        IF(IDDD.EQ.186)THEN
C
                           CALL VMX(1,PQ1,DPLQ3,A11,1)
C     
                           IEEE=1
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.187)THEN
C
                           CALL VMX(1,PQ2,DPLQ4,A11,1)
C     
                           IEEE=2
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.188)THEN
C
                           CALL VMX(1,PQ3,DPLQ1,A11,1)
C     
                           IEEE=3
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.189)THEN
C
                           CALL VMX(1,PQ4,DPLQ2,A11,1)
C     
                           IEEE=4
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.190)THEN
C
c                           DO 721 IC=1, NCOL2
c                              C11(IC)=DPLQ4(IC)
c 721                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ5,DPLQ7,A11,1)
C     
                           IEEE=5
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.191)THEN
C
c                           DO 722 IC=1, NCOL2
c                              C11(IC)=DPLQ3(IC)
c 722                       CONTINUE
C
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ6,DPLQ8,A11,1)
C     
                           IEEE=6
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.192)THEN
C
c                           DO 723 IC=1, NCOL2
c                              C11(IC)=DPLQ2(IC)
c 723                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ7,DPLQ5,A11,1)
C     
                           IEEE=7
C
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.193)THEN
C
c                           DO 724 IC=1, NCOL2
c                              C11(IC)=DPLQ1(IC)
c 724                       CONTINUE
cC
c                           CALL HERM(1,C11,DUM11,1)
                           CALL VMX(1,PQ8,DPLQ6,A11,1)
C     
                           IEEE=8
C
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 7
C******************************************************************C
                        IF(IDDD.EQ.194)THEN
C
                           CALL VMX(1,SQUY1,DPLQ1,B11,1)
                           CALL VMX(1,B11,SQDZ2,A11,1)
C     
                           IEEE=1
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.195)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                           CALL VMX(1,B11,SQUY2,A11,1)
C     
                           IEEE=2
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.196)THEN
C
                           CALL VMX(1,SQDY1,DPLQ3,B11,1)
                           CALL VMX(1,B11,SQUZ2,A11,1)
C     
                           IEEE=3
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.197)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                           CALL VMX(1,B11,SQDY2,A11,1)
C     
                           IEEE=4
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.198)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ6,B11,1)
                           CALL VMX(1,B11,SQUY2,A11,1)
C     
                           IEEE=5
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.199)THEN
C
                           CALL VMX(1,SQUY1,DPLQ5,B11,1)
                           CALL VMX(1,B11,SQUZ2,A11,1)
C     
                           IEEE=6
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.200)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ8,B11,1)
                           CALL VMX(1,B11,SQDY2,A11,1)
C     
                           IEEE=7
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.201)THEN
C
                           CALL VMX(1,SQDY1,DPLQ7,B11,1)
                           CALL VMX(1,B11,SQDZ2,A11,1)
C     
                           IEEE=8
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.202)THEN
C
                           CALL VMX(1,SQDY1,DPLQ5,B11,1)
                           CALL VMX(1,B11,SQDZ2,A11,1)
C     
                           IEEE=9
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.203)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                           CALL VMX(1,B11,SQDY2,A11,1)
C     
                           IEEE=10
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.204)THEN
C
                           CALL VMX(1,SQUY1,DPLQ7,B11,1)
                           CALL VMX(1,B11,SQUZ2,A11,1)
C     
                           IEEE=11
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.205)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                           CALL VMX(1,B11,SQUY2,A11,1)
C     
                           IEEE=12
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.206)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ2,B11,1)
                           CALL VMX(1,B11,SQDY2,A11,1)
C     
                           IEEE=13
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.207)THEN
C
                           CALL VMX(1,SQDY1,DPLQ1,B11,1)
                           CALL VMX(1,B11,SQUZ2,A11,1)
C     
                           IEEE=14
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.208)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ4,B11,1)
                           CALL VMX(1,B11,SQUY2,A11,1)
C     
                           IEEE=15
C
                        ENDIF
C******************************************************************C                        
                        IF(IDDD.EQ.209)THEN
C
                           CALL VMX(1,SQUY1,DPLQ3,B11,1)
                           CALL VMX(1,B11,SQDZ2,A11,1)
C     
                           IEEE=16
C
                        ENDIF                            
C******************************************************************C
C     PLAQUETTE OPERATORS 8
C******************************************************************C
                        IF(IDDD.EQ.210)THEN
C
                           CALL VMX(1,SQUY1,DPLQ1,B11,1)
                           CALL VMX(1,B11,DPLQ2,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.211)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                           CALL VMX(1,B11,DPLQ3,A11,1)
C     
                           IEEE=2
C     
                        ENDIF 
C******************************************************************c
                        IF(IDDD.EQ.212)THEN
C
                           CALL VMX(1,SQDY1,DPLQ3,B11,1)
                           CALL VMX(1,B11,DPLQ4,A11,1)
C     
                           IEEE=3
C     
                        ENDIF 
C******************************************************************c
                        IF(IDDD.EQ.213)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                           CALL VMX(1,B11,DPLQ1,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.214)THEN
C
                           CALL VMX(1,PLQ5,PLQ6,B11,1)
                           CALL VMX(1,B11,SQUY1,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.215)THEN
C
                           CALL VMX(1,PLQ8,PLQ5,B11,1)
                           CALL VMX(1,B11,SQUZ1,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.216)THEN
C
                           CALL VMX(1,PLQ7,PLQ8,B11,1)
                           CALL VMX(1,B11,SQDY1,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.217)THEN
C
                           CALL VMX(1,PLQ6,PLQ7,B11,1)
                           CALL VMX(1,B11,SQDZ1,A11,1)
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.218)THEN
C
                           CALL VMX(1,SQDY1,DPLQ5,B11,1)
                           CALL VMX(1,B11,DPLQ6,A11,1)
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.219)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                           CALL VMX(1,B11,DPLQ7,A11,1)
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.220)THEN
C
                           CALL VMX(1,SQUY1,DPLQ7,B11,1)
                           CALL VMX(1,B11,DPLQ8,A11,1)
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.221)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                           CALL VMX(1,B11,DPLQ5,A11,1)
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.222)THEN
C
                           CALL VMX(1,PLQ1,PLQ2,B11,1)
                           CALL VMX(1,B11,SQDY1,A11,1)
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.223)THEN
C
                           CALL VMX(1,PLQ4,PLQ1,B11,1)
                           CALL VMX(1,B11,SQUZ1,A11,1)
C     
                           IEEE=14
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.224)THEN
C
                           CALL VMX(1,PLQ3,PLQ4,B11,1)
                           CALL VMX(1,B11,SQUY1,A11,1)
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************c
                        IF(IDDD.EQ.225)THEN
C
                           CALL VMX(1,PLQ2,PLQ3,B11,1)
                           CALL VMX(1,B11,SQDZ1,A11,1)
C     
                           IEEE=16
C     
                        ENDIF                        
C******************************************************************C
C     PLAQUETTE OPERATORS 9
C******************************************************************C
                        IF(IDDD.EQ.226)THEN
C     
                           CALL VMX(1,SQUY1,DPLQ1,B11,1)
                           CALL VMX(1,B11,DPLQ2,C11,1)
                           CALL VMX(1,C11,DPLQ3,A11,1)                           
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.227)THEN
C     
                           CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                           CALL VMX(1,B11,DPLQ3,C11,1)
                           CALL VMX(1,C11,DPLQ4,A11,1)                           
C     
                           IEEE=2
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.228)THEN
C     
                           CALL VMX(1,SQDY1,DPLQ3,B11,1)
                           CALL VMX(1,B11,DPLQ4,C11,1)
                           CALL VMX(1,C11,DPLQ1,A11,1)                           
C     
                           IEEE=3
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.229)THEN
C     
                           CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                           CALL VMX(1,B11,DPLQ1,C11,1)
                           CALL VMX(1,C11,DPLQ2,A11,1)                           
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.230)THEN
C     
                           CALL VMX(1,PLQ8,PLQ5,B11,1)
                           CALL VMX(1,B11,PLQ6,C11,1)
                           CALL VMX(1,C11,SQUY1,A11,1)                           
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.231)THEN
C     
                           CALL VMX(1,PLQ7,PLQ8,B11,1)
                           CALL VMX(1,B11,PLQ5,C11,1)
                           CALL VMX(1,C11,SQUZ1,A11,1)                           
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.232)THEN
C     
                           CALL VMX(1,PLQ6,PLQ7,B11,1)
                           CALL VMX(1,B11,PLQ8,C11,1)
                           CALL VMX(1,C11,SQDY1,A11,1)                           
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.233)THEN
C     
                           CALL VMX(1,PLQ5,PLQ6,B11,1)
                           CALL VMX(1,B11,PLQ7,C11,1)
                           CALL VMX(1,C11,SQDZ1,A11,1)                           
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.234)THEN
C
                           CALL VMX(1,SQDY1,DPLQ5,B11,1)
                           CALL VMX(1,B11,DPLQ6,C11,1)
                           CALL VMX(1,C11,DPLQ7,A11,1)    
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.235)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                           CALL VMX(1,B11,DPLQ7,C11,1)
                           CALL VMX(1,C11,DPLQ8,A11,1)    
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.236)THEN
C
                           CALL VMX(1,SQUY1,DPLQ7,B11,1)
                           CALL VMX(1,B11,DPLQ8,C11,1)
                           CALL VMX(1,C11,DPLQ5,A11,1)    
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.237)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                           CALL VMX(1,B11,DPLQ5,C11,1)
                           CALL VMX(1,C11,DPLQ6,A11,1)    
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.238)THEN
C
                           CALL VMX(1,PLQ4,PLQ1,B11,1)
                           CALL VMX(1,B11,PLQ2,C11,1)
                           CALL VMX(1,C11,SQDY1,A11,1)    
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.239)THEN
C
                           CALL VMX(1,PLQ3,PLQ4,B11,1)
                           CALL VMX(1,B11,PLQ1,C11,1)
                           CALL VMX(1,C11,SQUZ1,A11,1)    
C     
                           IEEE=14
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.240)THEN
C
                           CALL VMX(1,PLQ2,PLQ3,B11,1)
                           CALL VMX(1,B11,PLQ4,C11,1)
                           CALL VMX(1,C11,SQUY1,A11,1)    
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.241)THEN
C
                           CALL VMX(1,PLQ1,PLQ2,B11,1)
                           CALL VMX(1,B11,PLQ3,C11,1)
                           CALL VMX(1,C11,SQDZ1,A11,1)    
C     
                           IEEE=16
C     
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 10
C******************************************************************C
                        IF(IDDD.EQ.242)THEN
C     
                           CALL VMX(1,PLQ2,PLQ7,B11,1)
                           CALL VMX(1,B11,SQDZ1,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.243)THEN
C     
                           CALL VMX(1,PLQ3,PLQ6,B11,1)
                           CALL VMX(1,B11,SQUY1,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.244)THEN
C     
                           CALL VMX(1,PLQ4,PLQ5,B11,1)
                           CALL VMX(1,B11,SQUZ1,A11,1)
C     
                           IEEE=3
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.245)THEN
C     
                           CALL VMX(1,PLQ1,PLQ8,B11,1)
                           CALL VMX(1,B11,SQDY1,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.246)THEN
C     
                           CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                           CALL VMX(1,B11,DPLQ5,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.247)THEN
C     
                           CALL VMX(1,SQUY1,DPLQ1,B11,1)
                           CALL VMX(1,B11,DPLQ8,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.248)THEN
C     
                           CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                           CALL VMX(1,B11,DPLQ7,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.249)THEN
C     
                           CALL VMX(1,SQDY1,DPLQ3,B11,1)
                           CALL VMX(1,B11,DPLQ6,A11,1)
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.250)THEN
C     
                           CALL VMX(1,PLQ6,PLQ3,B11,1)
                           CALL VMX(1,B11,SQDZ1,A11,1)
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.251)THEN
C     
                           CALL VMX(1,PLQ7,PLQ2,B11,1)
                           CALL VMX(1,B11,SQDY1,A11,1)
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.252)THEN
C     
                           CALL VMX(1,PLQ8,PLQ1,B11,1)
                           CALL VMX(1,B11,SQUZ1,A11,1)
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.253)THEN
C     
                           CALL VMX(1,PLQ5,PLQ4,B11,1)
                           CALL VMX(1,B11,SQUY1,A11,1)
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.254)THEN
C     
                           CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                           CALL VMX(1,B11,DPLQ1,A11,1)
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.255)THEN
C     
                           CALL VMX(1,SQDY1,DPLQ5,B11,1)
                           CALL VMX(1,B11,DPLQ4,A11,1)
C     
                           IEEE=14
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.256)THEN
C     
                           CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                           CALL VMX(1,B11,DPLQ3,A11,1)
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.257)THEN
C     
                           CALL VMX(1,SQUY1,DPLQ7,B11,1)
                           CALL VMX(1,B11,DPLQ2,A11,1)
C     
                           IEEE=16
C     
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 11
C******************************************************************C
                        IF(IDDD.EQ.258)THEN
C     
                           CALL VMX(1,PLQ2,PLQ4,B11,1)
                           CALL VMX(1,B11,SQUY1,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.259)THEN
C     
                           CALL VMX(1,PLQ3,PLQ1,B11,1)
                           CALL VMX(1,B11,SQUZ1,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.260)THEN
C     
                           CALL VMX(1,PLQ4,PLQ2,B11,1)
                           CALL VMX(1,B11,SQDY1,A11,1)
C     
                           IEEE=3
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.261)THEN
C     
                           CALL VMX(1,PLQ1,PLQ3,B11,1)
                           CALL VMX(1,B11,SQDZ1,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.262)THEN
C
                           CALL VMX(1,SQUY1,DPLQ7,B11,1)
                           CALL VMX(1,B11,DPLQ5,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.263)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ6,B11,1)
                           CALL VMX(1,B11,DPLQ8,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.264)THEN
C
                           CALL VMX(1,SQDY1,DPLQ5,B11,1)
                           CALL VMX(1,B11,DPLQ7,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.265)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ8,B11,1)
                           CALL VMX(1,B11,DPLQ6,A11,1)
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.266)THEN
C     
                           CALL VMX(1,PLQ6,PLQ8,B11,1)
                           CALL VMX(1,B11,SQDY1,A11,1)
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.267)THEN
C     
                           CALL VMX(1,PLQ7,PLQ5,B11,1)
                           CALL VMX(1,B11,SQUZ1,A11,1)
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.268)THEN
C     
                           CALL VMX(1,PLQ8,PLQ6,B11,1)
                           CALL VMX(1,B11,SQUY1,A11,1)
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.269)THEN
C     
                           CALL VMX(1,PLQ5,PLQ7,B11,1)
                           CALL VMX(1,B11,SQDZ1,A11,1)
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.270)THEN
C
                           CALL VMX(1,SQDY1,DPLQ3,B11,1)
                           CALL VMX(1,B11,DPLQ1,A11,1)
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.271)THEN
C
                           CALL VMX(1,SQUZ1,DPLQ2,B11,1)
                           CALL VMX(1,B11,DPLQ4,A11,1)
C     
                           IEEE=14
C     
                        ENDIF                        
C******************************************************************C
                        IF(IDDD.EQ.272)THEN
C
                           CALL VMX(1,SQUY1,DPLQ1,B11,1)
                           CALL VMX(1,B11,DPLQ3,A11,1)
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.273)THEN
C
                           CALL VMX(1,SQDZ1,DPLQ4,B11,1)
                           CALL VMX(1,B11,DPLQ2,A11,1)
C     
                           IEEE=16
C     
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 12
C******************************************************************C
                        IF(IDDD.EQ.274)THEN
C     
                           CALL VMX(1,PLQ2,PLQ3,B11,1)
                           CALL VMX(1,B11,SQUY1,C11,1)
                           CALL VMX(1,C11,DPLQ7,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.275)THEN
C     
                           CALL VMX(1,PLQ3,PLQ4,B11,1)
                           CALL VMX(1,B11,SQUZ1,C11,1)
                           CALL VMX(1,C11,DPLQ6,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.276)THEN
C     
                           CALL VMX(1,PLQ4,PLQ1,B11,1)
                           CALL VMX(1,B11,SQDY1,C11,1)
                           CALL VMX(1,C11,DPLQ5,A11,1)
C     
                           IEEE=3
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.277)THEN
C     
                           CALL VMX(1,PLQ1,PLQ2,B11,1)
                           CALL VMX(1,B11,SQDZ1,C11,1)
                           CALL VMX(1,C11,DPLQ8,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.278)THEN
C     
                           CALL VMX(1,PLQ4,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ8,C11,1)
                           CALL VMX(1,C11,DPLQ5,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.279)THEN
C     
                           CALL VMX(1,PLQ1,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ7,C11,1)
                           CALL VMX(1,C11,DPLQ8,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.280)THEN
C     
                           CALL VMX(1,PLQ2,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ6,C11,1)
                           CALL VMX(1,C11,DPLQ7,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.281)THEN
C     
                           CALL VMX(1,PLQ3,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ5,C11,1)
                           CALL VMX(1,C11,DPLQ6,A11,1)
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.282)THEN
C     
                           CALL VMX(1,PLQ6,PLQ7,B11,1)
                           CALL VMX(1,B11,SQDY1,C11,1)
                           CALL VMX(1,C11,DPLQ3,A11,1)
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.283)THEN
C     
                           CALL VMX(1,PLQ7,PLQ8,B11,1)
                           CALL VMX(1,B11,SQUZ1,C11,1)
                           CALL VMX(1,C11,DPLQ2,A11,1)
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.284)THEN
C     
                           CALL VMX(1,PLQ8,PLQ5,B11,1)
                           CALL VMX(1,B11,SQUY1,C11,1)
                           CALL VMX(1,C11,DPLQ1,A11,1)
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.285)THEN
C     
                           CALL VMX(1,PLQ5,PLQ6,B11,1)
                           CALL VMX(1,B11,SQDZ1,C11,1)
                           CALL VMX(1,C11,DPLQ4,A11,1)
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.286)THEN
C     
                           CALL VMX(1,PLQ8,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ4,C11,1)
                           CALL VMX(1,C11,DPLQ1,A11,1)
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.287)THEN
C     
                           CALL VMX(1,PLQ5,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ3,C11,1)
                           CALL VMX(1,C11,DPLQ4,A11,1)
C     
                           IEEE=14
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.288)THEN
C     
                           CALL VMX(1,PLQ6,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ2,C11,1)
                           CALL VMX(1,C11,DPLQ3,A11,1)
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.289)THEN
C     
                           CALL VMX(1,PLQ7,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ1,C11,1)
                           CALL VMX(1,C11,DPLQ2,A11,1)
C     
                           IEEE=16
C     
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 13
C******************************************************************C
                        IF(IDDD.EQ.290)THEN
C     
                           CALL VMX(1,PLQ2,PLQ3,B11,1)
                           CALL VMX(1,B11,PLQ4,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ1,C11,1)
                           CALL VMX(1,C11,DPLQ2,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.291)THEN
C     
                           CALL VMX(1,PLQ3,PLQ4,B11,1)
                           CALL VMX(1,B11,PLQ1,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ2,C11,1)
                           CALL VMX(1,C11,DPLQ3,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.292)THEN
C     
                           CALL VMX(1,PLQ4,PLQ1,B11,1)
                           CALL VMX(1,B11,PLQ2,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ3,C11,1)
                           CALL VMX(1,C11,DPLQ4,A11,1)
C     
                           IEEE=3
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.293)THEN
C     
                           CALL VMX(1,PLQ1,PLQ2,B11,1)
                           CALL VMX(1,B11,PLQ3,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ4,C11,1)
                           CALL VMX(1,C11,DPLQ1,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.294)THEN
C     
                           CALL VMX(1,PLQ5,PLQ6,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ7,B11,1)
                           CALL VMX(1,B11,DPLQ8,C11,1)
                           CALL VMX(1,C11,DPLQ5,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.295)THEN
C     
                           CALL VMX(1,PLQ8,PLQ5,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ6,B11,1)
                           CALL VMX(1,B11,DPLQ7,C11,1)
                           CALL VMX(1,C11,DPLQ8,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.296)THEN
C     
                           CALL VMX(1,PLQ7,PLQ8,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ5,B11,1)
                           CALL VMX(1,B11,DPLQ6,C11,1)
                           CALL VMX(1,C11,DPLQ7,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.297)THEN
C     
                           CALL VMX(1,PLQ6,PLQ7,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ8,B11,1)
                           CALL VMX(1,B11,DPLQ5,C11,1)
                           CALL VMX(1,C11,DPLQ6,A11,1)
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.298)THEN
C     
                           CALL VMX(1,PLQ6,PLQ7,B11,1)
                           CALL VMX(1,B11,PLQ8,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ5,C11,1)
                           CALL VMX(1,C11,DPLQ6,A11,1)
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.299)THEN
C     
                           CALL VMX(1,PLQ7,PLQ8,B11,1)
                           CALL VMX(1,B11,PLQ5,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ6,C11,1)
                           CALL VMX(1,C11,DPLQ7,A11,1)
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.300)THEN
C     
                           CALL VMX(1,PLQ8,PLQ5,B11,1)
                           CALL VMX(1,B11,PLQ6,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ7,C11,1)
                           CALL VMX(1,C11,DPLQ8,A11,1)
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.301)THEN
C     
                           CALL VMX(1,PLQ5,PLQ6,B11,1)
                           CALL VMX(1,B11,PLQ7,C11,1)
                           CALL VMX(1,C11,PL,B11,1)
                           CALL VMX(1,B11,DPLQ8,C11,1)
                           CALL VMX(1,C11,DPLQ5,A11,1)
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.302)THEN
C     
                           CALL VMX(1,PLQ1,PLQ2,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ3,B11,1)
                           CALL VMX(1,B11,DPLQ4,C11,1)
                           CALL VMX(1,C11,DPLQ1,A11,1)
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.303)THEN
C     
                           CALL VMX(1,PLQ4,PLQ1,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ2,B11,1)
                           CALL VMX(1,B11,DPLQ3,C11,1)
                           CALL VMX(1,C11,DPLQ4,A11,1)
C     
                           IEEE=14
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.304)THEN
C     
                           CALL VMX(1,PLQ3,PLQ4,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ1,B11,1)
                           CALL VMX(1,B11,DPLQ2,C11,1)
                           CALL VMX(1,C11,DPLQ3,A11,1)
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.305)THEN
C     
                           CALL VMX(1,PLQ2,PLQ3,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ4,B11,1)
                           CALL VMX(1,B11,DPLQ1,C11,1)
                           CALL VMX(1,C11,DPLQ2,A11,1)
C     
                           IEEE=16
C     
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 14
C******************************************************************C
                        IF(IDDD.EQ.306)THEN
C     
                           CALL VMX(1,PLQ2,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ7,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.307)THEN
C     
                           CALL VMX(1,PLQ3,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ6,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.308)THEN
C     
                           CALL VMX(1,PLQ4,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ5,A11,1)
C     
                           IEEE=3
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.309)THEN
C     
                           CALL VMX(1,PLQ1,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ8,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.310)THEN
C     
                           CALL VMX(1,PLQ4,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ5,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.311)THEN
C     
                           CALL VMX(1,PLQ1,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ8,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.312)THEN
C     
                           CALL VMX(1,PLQ2,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ7,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.313)THEN
C     
                           CALL VMX(1,PLQ3,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ6,A11,1)
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.314)THEN
C     
                           CALL VMX(1,PLQ6,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ3,A11,1)
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.315)THEN
C     
                           CALL VMX(1,PLQ7,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ2,A11,1)
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.316)THEN
C     
                           CALL VMX(1,PLQ8,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ1,A11,1)
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.317)THEN
C     
                           CALL VMX(1,PLQ5,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ4,A11,1)
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.318)THEN
C     
                           CALL VMX(1,PLQ8,SQDY1,B11,1)
                           CALL VMX(1,B11,DPLQ1,A11,1)
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.319)THEN
C     
                           CALL VMX(1,PLQ5,SQUZ1,B11,1)
                           CALL VMX(1,B11,DPLQ4,A11,1)
C     
                           IEEE=14
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.320)THEN
C     
                           CALL VMX(1,PLQ6,SQUY1,B11,1)
                           CALL VMX(1,B11,DPLQ3,A11,1)
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.321)THEN
C     
                           CALL VMX(1,PLQ7,SQDZ1,B11,1)
                           CALL VMX(1,B11,DPLQ2,A11,1)
C     
                           IEEE=16
C     
                        ENDIF
C******************************************************************C
C     PLAQUETTE OPERATORS 15
C******************************************************************C
                        IF(IDDD.EQ.322)THEN
C     
                           CALL VMX(1,PLQ2,PLQ7,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ6,A11,1)
C     
                           IEEE=1
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.323)THEN
C     
                           CALL VMX(1,PLQ3,PLQ6,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ5,A11,1)
C     
                           IEEE=2
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.324)THEN
C     
                           CALL VMX(1,PLQ4,PLQ5,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ8,A11,1)
C     
                           IEEE=3
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.325)THEN
C     
                           CALL VMX(1,PLQ1,PLQ8,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ7,A11,1)
C     
                           IEEE=4
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.326)THEN
C     
                           CALL VMX(1,PLQ1,PL,B11,1)
                           CALL VMX(1,B11,DPLQ4,C11,1)
                           CALL VMX(1,C11,DPLQ5,A11,1)
C     
                           IEEE=5
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.327)THEN
C     
                           CALL VMX(1,PLQ2,PL,B11,1)
                           CALL VMX(1,B11,DPLQ1,C11,1)
                           CALL VMX(1,C11,DPLQ8,A11,1)
C     
                           IEEE=6
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.328)THEN
C     
                           CALL VMX(1,PLQ3,PL,B11,1)
                           CALL VMX(1,B11,DPLQ2,C11,1)
                           CALL VMX(1,C11,DPLQ7,A11,1)
C     
                           IEEE=7
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.329)THEN
C     
                           CALL VMX(1,PLQ4,PL,B11,1)
                           CALL VMX(1,B11,DPLQ3,C11,1)
                           CALL VMX(1,C11,DPLQ6,A11,1)
C     
                           IEEE=8
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.330)THEN
C     
                           CALL VMX(1,PLQ6,PLQ3,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ2,A11,1)
C     
                           IEEE=9
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.331)THEN
C     
                           CALL VMX(1,PLQ7,PLQ2,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ1,A11,1)
C     
                           IEEE=10
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.332)THEN
C     
                           CALL VMX(1,PLQ8,PLQ1,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ4,A11,1)
C     
                           IEEE=11
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.333)THEN
C     
                           CALL VMX(1,PLQ5,PLQ4,B11,1)
                           CALL VMX(1,B11,PL,C11,1)
                           CALL VMX(1,C11,DPLQ3,A11,1)
C     
                           IEEE=12
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.334)THEN
C     
                           CALL VMX(1,PLQ5,PL,B11,1)
                           CALL VMX(1,B11,DPLQ8,C11,1)
                           CALL VMX(1,C11,DPLQ1,A11,1)
C     
                           IEEE=13
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.335)THEN
C     
                           CALL VMX(1,PLQ6,PL,B11,1)
                           CALL VMX(1,B11,DPLQ5,C11,1)
                           CALL VMX(1,C11,DPLQ4,A11,1)
C     
                           IEEE=14
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.336)THEN
C     
                           CALL VMX(1,PLQ7,PL,B11,1)
                           CALL VMX(1,B11,DPLQ6,C11,1)
                           CALL VMX(1,C11,DPLQ3,A11,1)
C     
                           IEEE=15
C     
                        ENDIF
C******************************************************************C
                        IF(IDDD.EQ.337)THEN
C     
                           CALL VMX(1,PLQ8,PL,B11,1)
                           CALL VMX(1,B11,DPLQ7,C11,1)
                           CALL VMX(1,C11,DPLQ2,A11,1)
C     
                           IEEE=16
C     
                        ENDIF                        
C******************************************************************c
         ENDIF
C******************************************************************C
               IF(LCNT(IDS).LT.ICO) THEN
                  DO 168 IC=1, NCOL2
                     A11(IC)=LIN0(IC)
 168              CONTINUE
                  M2=ML
               ELSE
                  IF(ICO.EQ.1) THEN
                     CALL VMX(1,A11,LIN1,C11,1)
                     DO 169 IC=1,NCOL2
                        A11(IC)=C11(IC)
 169                 CONTINUE
                  ENDIF
C
                  IF(ICO.EQ.2) THEN
                     CALL VMX(1,A11,LIN2,C11,1)
                     DO 170 IC=1,NCOL2
                        A11(IC)=C11(IC)
 170                 CONTINUE
                  ENDIF
C
                  IF(ICO.EQ.4) THEN
                     CALL VMX(1,A11,LIN4,C11,1)
                     DO 171 IC=1,NCOL2
                        A11(IC)=C11(IC)
 171                 CONTINUE
                  ENDIF
                  M2=ML
               ENDIF
C***********************************************************************
            ENDIF
            CALL VMX(1,A11,REM11,B11,1)
C***********************************************************************
c            CALL RENORMBS(B11,UREN11)
c            DO IJ=1,NCOL2
c               AA(IJ)=UREN11(IJ)
c            ENDDO
c            JMAT=NCOL
c            CALL DETNANT(JMAT,DET,AA)
cc     
c            DNCOL=CMPLX(1.0/NCOL)
c            CDET=DET**DNCOL
c            CNORM=1.0/CDET
cc     
c            DO IC=1,NCOL2
c               B11(IC)=UREN11(IC)*CNORM
c            ENDDO
C***********************************************************************
            AKT1=(0.0,0.0)
C***********************************************************************
            DO 192 N1=1,NCOL
               IJ=N1+NCOL*(N1-1)
               AKT1=AKT1+B11(IJ)
 192        CONTINUE
C************ ***********************************************************     
            IF(IDDD.LE.4) THEN
               CSUMS(IEEE)=CSUMS(IEEE)+AKT1
               CSUMSMOM(IEEE,1)=CSUMSMOM(IEEE,1)+AKT1*PF(1)
               CSUMSMOM(IEEE,2)=CSUMSMOM(IEEE,2)+AKT1*PF(2)
c               CSUMSMOM(IEEE,3)=CSUMSMOM(IEEE,3)+AKT1*PF(3)
c               CSUMSMOM(IEEE,4)=CSUMSMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.4).AND.(IDDD.LT.9)) THEN
               CSUM2S(IEEE)=CSUM2S(IEEE)+AKT1
               CSUM2SMOM(IEEE,1)=CSUM2SMOM(IEEE,1)+AKT1*PF(1)
               CSUM2SMOM(IEEE,2)=CSUM2SMOM(IEEE,2)+AKT1*PF(2)
c               CSUM2SMOM(IEEE,3)=CSUM2SMOM(IEEE,3)+AKT1*PF(3)
c               CSUM2SMOM(IEEE,4)=CSUM2SMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.8).AND.(IDDD.LT.13)) THEN
               CSUM2WS(IEEE)=CSUM2WS(IEEE)+AKT1
               CSUM2WSMOM(IEEE,1)=CSUM2WSMOM(IEEE,1)+AKT1*PF(1)
               CSUM2WSMOM(IEEE,2)=CSUM2WSMOM(IEEE,2)+AKT1*PF(2)
c               CSUM2WSMOM(IEEE,3)=CSUM2WSMOM(IEEE,3)+AKT1*PF(3)
c               CSUM2WSMOM(IEEE,4)=CSUM2WSMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.12).AND.(IDDD.LT.17)) THEN
               CSUMW(IEEE)=CSUMW(IEEE)+AKT1
               CSUMWMOM(IEEE,1)=CSUMWMOM(IEEE,1)+AKT1*PF(1)
               CSUMWMOM(IEEE,2)=CSUMWMOM(IEEE,2)+AKT1*PF(2)
c               CSUMWMOM(IEEE,3)=CSUMWMOM(IEEE,3)+AKT1*PF(3)
c               CSUMWMOM(IEEE,4)=CSUMWMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.16).AND.(IDDD.LT.21)) THEN
               CSUM2W(IEEE)=CSUM2W(IEEE)+AKT1
               CSUM2WMOM(IEEE,1)=CSUM2WMOM(IEEE,1)+AKT1*PF(1)
               CSUM2WMOM(IEEE,2)=CSUM2WMOM(IEEE,2)+AKT1*PF(2)
c               CSUM2WMOM(IEEE,3)=CSUM2WMOM(IEEE,3)+AKT1*PF(3)
c               CSUM2WMOM(IEEE,4)=CSUM2WMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.20).AND.(IDDD.LT.25)) THEN
               CSUM3W(IEEE)=CSUM3W(IEEE)+AKT1
               CSUM3WMOM(IEEE,1)=CSUM3WMOM(IEEE,1)+AKT1*PF(1)
               CSUM3WMOM(IEEE,2)=CSUM3WMOM(IEEE,2)+AKT1*PF(2)
c               CSUM3WMOM(IEEE,3)=CSUM3WMOM(IEEE,3)+AKT1*PF(3)
c               CSUM3WMOM(IEEE,4)=CSUM3WMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.24).AND.(IDDD.LT.29)) THEN
               CSUMUP(IEEE)=CSUMUP(IEEE)+AKT1
               CSUMUPMOM(IEEE,1)=CSUMUPMOM(IEEE,1)+AKT1*PF(1)
               CSUMUPMOM(IEEE,2)=CSUMUPMOM(IEEE,2)+AKT1*PF(2)
c               CSUMUPMOM(IEEE,3)=CSUMUPMOM(IEEE,3)+AKT1*PF(3)
c               CSUMUPMOM(IEEE,4)=CSUMUPMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.28).AND.(IDDD.LT.33)) THEN
               CSUMUD(IEEE)=CSUMUD(IEEE)+AKT1
               CSUMUDMOM(IEEE,1)=CSUMUDMOM(IEEE,1)+AKT1*PF(1)
               CSUMUDMOM(IEEE,2)=CSUMUDMOM(IEEE,2)+AKT1*PF(2)
c               CSUMUDMOM(IEEE,3)=CSUMUDMOM(IEEE,3)+AKT1*PF(3)
c               CSUMUDMOM(IEEE,4)=CSUMUDMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.32).AND.(IDDD.LT.41)) THEN
               CSUMTT1(IEEE)=CSUMTT1(IEEE)+AKT1
               CSUMTTMOM1(IEEE,1)=CSUMTTMOM1(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM1(IEEE,2)=CSUMTTMOM1(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM1(IEEE,3)=CSUMTTMOM1(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM1(IEEE,4)=CSUMTTMOM1(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.40).AND.(IDDD.LT.49)) THEN
               CSUMTT2(IEEE)=CSUMTT2(IEEE)+AKT1
               CSUMTTMOM2(IEEE,1)=CSUMTTMOM2(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM2(IEEE,2)=CSUMTTMOM2(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM2(IEEE,3)=CSUMTTMOM2(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM2(IEEE,4)=CSUMTTMOM2(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.48).AND.(IDDD.LT.57)) THEN
               CSUMTT3(IEEE)=CSUMTT3(IEEE)+AKT1
               CSUMTTMOM3(IEEE,1)=CSUMTTMOM3(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM3(IEEE,2)=CSUMTTMOM3(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM3(IEEE,3)=CSUMTTMOM3(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM3(IEEE,4)=CSUMTTMOM3(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.56).AND.(IDDD.LT.65)) THEN
               CSUMTT4(IEEE)=CSUMTT4(IEEE)+AKT1
               CSUMTTMOM4(IEEE,1)=CSUMTTMOM4(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM4(IEEE,2)=CSUMTTMOM4(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM4(IEEE,3)=CSUMTTMOM4(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM4(IEEE,4)=CSUMTTMOM4(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.64).AND.(IDDD.LT.69)) THEN
               CSUMTT5(IEEE)=CSUMTT5(IEEE)+AKT1
               CSUMTTMOM5(IEEE,1)=CSUMTTMOM5(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM5(IEEE,2)=CSUMTTMOM5(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM5(IEEE,3)=CSUMTTMOM5(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM5(IEEE,4)=CSUMTTMOM5(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.68).AND.(IDDD.LT.73)) THEN
               CSUMTT6(IEEE)=CSUMTT6(IEEE)+AKT1
               CSUMTTMOM6(IEEE,1)=CSUMTTMOM6(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM6(IEEE,2)=CSUMTTMOM6(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM6(IEEE,3)=CSUMTTMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM6(IEEE,4)=CSUMTTMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF
C            
            IF((IDDD.GT.72).AND.(IDDD.LT.81)) THEN
               CSUMTT7(IEEE)=CSUMTT7(IEEE)+AKT1
               CSUMTTMOM7(IEEE,1)=CSUMTTMOM7(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM7(IEEE,2)=CSUMTTMOM7(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM7(IEEE,3)=CSUMTTMOM7(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM7(IEEE,4)=CSUMTTMOM7(IEEE,4)+AKT1*PF(4)
            ENDIF
C            
            IF((IDDD.GT.80).AND.(IDDD.LT.89)) THEN
               CSUMTT8(IEEE)=CSUMTT8(IEEE)+AKT1
               CSUMTTMOM8(IEEE,1)=CSUMTTMOM8(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM8(IEEE,2)=CSUMTTMOM8(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM8(IEEE,3)=CSUMTTMOM8(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM8(IEEE,4)=CSUMTTMOM8(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.88).AND.(IDDD.LT.97)) THEN
               CSUMTT9(IEEE)=CSUMTT9(IEEE)+AKT1
               CSUMTTMOM9(IEEE,1)=CSUMTTMOM9(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM9(IEEE,2)=CSUMTTMOM9(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM9(IEEE,3)=CSUMTTMOM9(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM9(IEEE,4)=CSUMTTMOM9(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.96).AND.(IDDD.LT.105)) THEN
               CSUMTT10(IEEE)=CSUMTT10(IEEE)+AKT1
               CSUMTTMOM10(IEEE,1)=CSUMTTMOM10(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM10(IEEE,2)=CSUMTTMOM10(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM10(IEEE,3)=CSUMTTMOM10(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM10(IEEE,4)=CSUMTTMOM10(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.104).AND.(IDDD.LT.113)) THEN
               CSUMTT11(IEEE)=CSUMTT11(IEEE)+AKT1
               CSUMTTMOM11(IEEE,1)=CSUMTTMOM11(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM11(IEEE,2)=CSUMTTMOM11(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM11(IEEE,3)=CSUMTTMOM11(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM11(IEEE,4)=CSUMTTMOM11(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.112).AND.(IDDD.LT.129)) THEN
               CSUMTT12(IEEE)=CSUMTT12(IEEE)+AKT1
               CSUMTTMOM12(IEEE,1)=CSUMTTMOM12(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM12(IEEE,2)=CSUMTTMOM12(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM12(IEEE,3)=CSUMTTMOM12(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM12(IEEE,4)=CSUMTTMOM12(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.128).AND.(IDDD.LT.137)) THEN
               CSUMTT13(IEEE)=CSUMTT13(IEEE)+AKT1
               CSUMTTMOM13(IEEE,1)=CSUMTTMOM13(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM13(IEEE,2)=CSUMTTMOM13(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM13(IEEE,3)=CSUMTTMOM13(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM13(IEEE,4)=CSUMTTMOM13(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.136).AND.(IDDD.LT.145)) THEN
               CSUMTT14(IEEE)=CSUMTT14(IEEE)+AKT1
               CSUMTTMOM14(IEEE,1)=CSUMTTMOM14(IEEE,1)+AKT1*PF(1)
               CSUMTTMOM14(IEEE,2)=CSUMTTMOM14(IEEE,2)+AKT1*PF(2)
c               CSUMTTMOM14(IEEE,3)=CSUMTTMOM14(IEEE,3)+AKT1*PF(3)
c               CSUMTTMOM14(IEEE,4)=CSUMTTMOM14(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF(IDDD.EQ.145) THEN
            CSUMN=CSUMN+AKT1
            ENDIF
C
            IF((IDDD.GT.145).AND.(IDDD.LT.154)) THEN
               CSUMPLQ(IEEE)=CSUMPLQ(IEEE)+AKT1
               CSUMPLQMOM(IEEE,1)=CSUMPLQMOM(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM(IEEE,2)=CSUMPLQMOM(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM(IEEE,3)=CSUMPLQMOM(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM(IEEE,4)=CSUMPLQMOM(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.153).AND.(IDDD.LT.162)) THEN
               CSUMPLQ2(IEEE)=CSUMPLQ2(IEEE)+AKT1
               CSUMPLQMOM2(IEEE,1)=CSUMPLQMOM2(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM2(IEEE,2)=CSUMPLQMOM2(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM2(IEEE,3)=CSUMPLQMOM2(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM2(IEEE,4)=CSUMPLQMOM2(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.161).AND.(IDDD.LT.170)) THEN
               CSUMPLQ3(IEEE)=CSUMPLQ3(IEEE)+AKT1
               CSUMPLQMOM3(IEEE,1)=CSUMPLQMOM3(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM3(IEEE,2)=CSUMPLQMOM3(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM3(IEEE,3)=CSUMPLQMOM3(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM3(IEEE,4)=CSUMPLQMOM3(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.169).AND.(IDDD.LT.178)) THEN
               CSUMPLQ4(IEEE)=CSUMPLQ4(IEEE)+AKT1
               CSUMPLQMOM4(IEEE,1)=CSUMPLQMOM4(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM4(IEEE,2)=CSUMPLQMOM4(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM4(IEEE,3)=CSUMPLQMOM4(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM4(IEEE,4)=CSUMPLQMOM4(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.177).AND.(IDDD.LT.186)) THEN
               CSUMPLQ5(IEEE)=CSUMPLQ5(IEEE)+AKT1
               CSUMPLQMOM5(IEEE,1)=CSUMPLQMOM5(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM5(IEEE,2)=CSUMPLQMOM5(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM5(IEEE,3)=CSUMPLQMOM5(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM5(IEEE,4)=CSUMPLQMOM5(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.185).AND.(IDDD.LT.194)) THEN
               CSUMPLQ6(IEEE)=CSUMPLQ6(IEEE)+AKT1
               CSUMPLQMOM6(IEEE,1)=CSUMPLQMOM6(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM6(IEEE,2)=CSUMPLQMOM6(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF
C
            IF((IDDD.GT.193).AND.(IDDD.LT.210)) THEN
               CSUMPLQ7(IEEE)=CSUMPLQ7(IEEE)+AKT1
               CSUMPLQMOM7(IEEE,1)=CSUMPLQMOM7(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM7(IEEE,2)=CSUMPLQMOM7(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF
c
            IF((IDDD.GT.209).AND.(IDDD.LT.226)) THEN 
               CSUMPLQ8(IEEE)=CSUMPLQ8(IEEE)+AKT1
               CSUMPLQMOM8(IEEE,1)=CSUMPLQMOM8(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM8(IEEE,2)=CSUMPLQMOM8(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF
c
            IF((IDDD.GT.225).AND.(IDDD.LT.242)) THEN 
               CSUMPLQ9(IEEE)=CSUMPLQ9(IEEE)+AKT1
               CSUMPLQMOM9(IEEE,1)=CSUMPLQMOM9(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM9(IEEE,2)=CSUMPLQMOM9(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF
c
            IF((IDDD.GT.241).AND.(IDDD.LT.258)) THEN 
               CSUMPLQ10(IEEE)=CSUMPLQ10(IEEE)+AKT1
               CSUMPLQMOM10(IEEE,1)=CSUMPLQMOM10(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM10(IEEE,2)=CSUMPLQMOM10(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF   
c
            IF((IDDD.GT.257).AND.(IDDD.LT.274)) THEN 
               CSUMPLQ11(IEEE)=CSUMPLQ11(IEEE)+AKT1
               CSUMPLQMOM11(IEEE,1)=CSUMPLQMOM11(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM11(IEEE,2)=CSUMPLQMOM11(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF
c
            IF((IDDD.GT.273).AND.(IDDD.LT.290)) THEN 
               CSUMPLQ12(IEEE)=CSUMPLQ12(IEEE)+AKT1
               CSUMPLQMOM12(IEEE,1)=CSUMPLQMOM12(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM12(IEEE,2)=CSUMPLQMOM12(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF  
c
            IF((IDDD.GT.289).AND.(IDDD.LT.306)) THEN 
               CSUMPLQ13(IEEE)=CSUMPLQ13(IEEE)+AKT1
               CSUMPLQMOM13(IEEE,1)=CSUMPLQMOM13(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM13(IEEE,2)=CSUMPLQMOM13(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF  
c
            IF((IDDD.GT.305).AND.(IDDD.LT.322)) THEN 
               CSUMPLQ14(IEEE)=CSUMPLQ14(IEEE)+AKT1
               CSUMPLQMOM14(IEEE,1)=CSUMPLQMOM14(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM14(IEEE,2)=CSUMPLQMOM14(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF  
c
            IF((IDDD.GT.321).AND.(IDDD.LT.338)) THEN 
               CSUMPLQ15(IEEE)=CSUMPLQ15(IEEE)+AKT1
               CSUMPLQMOM15(IEEE,1)=CSUMPLQMOM15(IEEE,1)+AKT1*PF(1)
               CSUMPLQMOM15(IEEE,2)=CSUMPLQMOM15(IEEE,2)+AKT1*PF(2)
c               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
c               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
            ENDIF  
C     
 9       CONTINUE 
c
      ENDDO
      ENDDO
      ENDDO
C**********************************************************************     
      ADIV1=1.0/(4.0*LS(KU))
      ADIV2=1.0/(8.0*LS(KU))
      ADIV3=1.0/(16.0*LS(KU))
      ADIVN=1.0/LS(KU)
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE1(N4,ID)=CSUMN*ADIVN
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE2(N4,ID)=(CSUMS(1)+CSUMS(2)+CSUMS(3)+CSUMS(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM2(N4,ID,1)=(CSUMSMOM(1,1)+CSUMSMOM(2,1)
     &+CSUMSMOM(3,1)+CSUMSMOM(4,1))*ADIV1
      ALINEMOM2(N4,ID,2)=(CSUMSMOM(1,2)+CSUMSMOM(2,2)
     &+CSUMSMOM(3,2)+CSUMSMOM(4,2))*ADIV1
c      ALINEMOM2(N4,ID,3)=(CSUMSMOM(1,3)+CSUMSMOM(2,3)
c     &+CSUMSMOM(3,3)+CSUMSMOM(4,3))*ADIV1
c      ALINEMOM2(N4,ID,4)=(CSUMSMOM(1,4)+CSUMSMOM(2,4)
c     &+CSUMSMOM(3,4)+CSUMSMOM(4,4))*ADIV1
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE3(N4,ID)=(CSUMS(1)+GIOT*CSUMS(2)-CSUMS(3)-GIOT*CSUMS(4))
     &*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM3(N4,ID,1)=(CSUMSMOM(1,1)+GIOT*CSUMSMOM(2,1)
     &-CSUMSMOM(3,1)-GIOT*CSUMSMOM(4,1))*ADIV1
      ALINEMOM3(N4,ID,2)=(CSUMSMOM(1,2)+GIOT*CSUMSMOM(2,2)
     &-CSUMSMOM(3,2)-GIOT*CSUMSMOM(4,2))*ADIV1
c      ALINEMOM3(N4,ID,3)=(CSUMSMOM(1,3)+GIOT*CSUMSMOM(2,3)
c     &-CSUMSMOM(3,3)-GIOT*CSUMSMOM(4,3))*ADIV1
c      ALINEMOM3(N4,ID,4)=(CSUMSMOM(1,4)+GIOT*CSUMSMOM(2,4)
c     &-CSUMSMOM(3,4)-GIOT*CSUMSMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE4(N4,ID)=(CSUMS(1)-CSUMS(2)+CSUMS(3)-CSUMS(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM4(N4,ID,1)=(CSUMSMOM(1,1)-CSUMSMOM(2,1)
     &+CSUMSMOM(3,1)-CSUMSMOM(4,1))*ADIV1
      ALINEMOM4(N4,ID,2)=(CSUMSMOM(1,2)-CSUMSMOM(2,2)
     &+CSUMSMOM(3,2)-CSUMSMOM(4,2))*ADIV1
c      ALINEMOM4(N4,ID,3)=(CSUMSMOM(1,3)-CSUMSMOM(2,3)
c     &+CSUMSMOM(3,3)-CSUMSMOM(4,3))*ADIV1
c      ALINEMOM4(N4,ID,4)=(CSUMSMOM(1,4)-CSUMSMOM(2,4)
c     &+CSUMSMOM(3,4)-CSUMSMOM(4,4))*ADIV1
c**********************************************************************
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE5(N4,ID)=(CSUM2S(1)+CSUM2S(2)+CSUM2S(3)+CSUM2S(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM5(N4,ID,1)=(CSUM2SMOM(1,1)+CSUM2SMOM(2,1)
     &+CSUM2SMOM(3,1)+CSUM2SMOM(4,1))*ADIV1
      ALINEMOM5(N4,ID,2)=(CSUM2SMOM(1,2)+CSUM2SMOM(2,2)
     &+CSUM2SMOM(3,2)+CSUM2SMOM(4,2))*ADIV1
c      ALINEMOM5(N4,ID,3)=(CSUM2SMOM(1,3)+CSUM2SMOM(2,3)
c     &+CSUM2SMOM(3,3)+CSUM2SMOM(4,3))*ADIV1
c      ALINEMOM5(N4,ID,4)=(CSUM2SMOM(1,4)+CSUM2SMOM(2,4)
c     &+CSUM2SMOM(3,4)+CSUM2SMOM(4,4))*ADIV1
C***********************************************************************
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE6(N4,ID)=(CSUM2S(1)+GIOT*CSUM2S(2)-CSUM2S(3)
     &-GIOT*CSUM2S(4))*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM6(N4,ID,1)=(CSUM2SMOM(1,1)+GIOT*CSUM2SMOM(2,1)
     &-CSUM2SMOM(3,1)-GIOT*CSUM2SMOM(4,1))*ADIV1
      ALINEMOM6(N4,ID,2)=(CSUM2SMOM(1,2)+GIOT*CSUM2SMOM(2,2)
     &-CSUM2SMOM(3,2)-GIOT*CSUM2SMOM(4,2))*ADIV1
c      ALINEMOM6(N4,ID,3)=(CSUM2SMOM(1,3)+GIOT*CSUM2SMOM(2,3)
c     &-CSUM2SMOM(3,3)-GIOT*CSUM2SMOM(4,3))*ADIV1
c      ALINEMOM6(N4,ID,4)=(CSUM2SMOM(1,4)+GIOT*CSUM2SMOM(2,4)
c     &-CSUM2SMOM(3,4)-GIOT*CSUM2SMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE7(N4,ID)=(CSUM2S(1)-CSUM2S(2)+CSUM2S(3)-CSUM2S(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM7(N4,ID,1)=(CSUM2SMOM(1,1)-CSUM2SMOM(2,1)
     &+CSUM2SMOM(3,1)-CSUM2SMOM(4,1))*ADIV1
      ALINEMOM7(N4,ID,2)=(CSUM2SMOM(1,2)-CSUM2SMOM(2,2)
     &+CSUM2SMOM(3,2)-CSUM2SMOM(4,2))*ADIV1
c      ALINEMOM7(N4,ID,3)=(CSUM2SMOM(1,3)-CSUM2SMOM(2,3)
c     &+CSUM2SMOM(3,3)-CSUM2SMOM(4,3))*ADIV1
c      ALINEMOM7(N4,ID,4)=(CSUM2SMOM(1,4)-CSUM2SMOM(2,4)
c     &+CSUM2SMOM(3,4)-CSUM2SMOM(4,4))*ADIV1
C**********************************************************************
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE8(N4,ID)=(CSUM2WS(1)+CSUM2WS(2)+CSUM2WS(3)+CSUM2WS(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM8(N4,ID,1)=(CSUM2WSMOM(1,1)+CSUM2WSMOM(2,1)
     &+CSUM2WSMOM(3,1)+CSUM2WSMOM(4,1))*ADIV1
      ALINEMOM8(N4,ID,2)=(CSUM2WSMOM(1,2)+CSUM2WSMOM(2,2)
     &+CSUM2WSMOM(3,2)+CSUM2WSMOM(4,2))*ADIV1
c      ALINEMOM8(N4,ID,3)=(CSUM2WSMOM(1,3)+CSUM2WSMOM(2,3)
c     &+CSUM2WSMOM(3,3)+CSUM2WSMOM(4,3))*ADIV1
c      ALINEMOM8(N4,ID,4)=(CSUM2WSMOM(1,4)+CSUM2WSMOM(2,4)
c     &+CSUM2WSMOM(3,4)+CSUM2WSMOM(4,4))*ADIV1
C************************************************************************
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE9(N4,ID)=(CSUM2WS(1)+GIOT*CSUM2WS(2)-CSUM2WS(3)
     &-GIOT*CSUM2WS(4))*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM9(N4,ID,1)=(CSUM2WSMOM(1,1)+GIOT*CSUM2WSMOM(2,1)
     &-CSUM2WSMOM(3,1)-GIOT*CSUM2WSMOM(4,1))*ADIV1
      ALINEMOM9(N4,ID,2)=(CSUM2WSMOM(1,2)+GIOT*CSUM2WSMOM(2,2)
     &-CSUM2WSMOM(3,2)-GIOT*CSUM2WSMOM(4,2))*ADIV1
c      ALINEMOM9(N4,ID,3)=(CSUM2WSMOM(1,3)+GIOT*CSUM2WSMOM(2,3)
c     &-CSUM2WSMOM(3,3)-GIOT*CSUM2WSMOM(4,3))*ADIV1
c      ALINEMOM9(N4,ID,4)=(CSUM2WSMOM(1,4)+GIOT*CSUM2WSMOM(2,4)
c     &-CSUM2WSMOM(3,4)-GIOT*CSUM2WSMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE10(N4,ID)=(CSUM2WS(1)-CSUM2WS(2)+CSUM2WS(3)-CSUM2WS(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM10(N4,ID,1)=(CSUM2WSMOM(1,1)-CSUM2WSMOM(2,1)
     &+CSUM2WSMOM(3,1)-CSUM2WSMOM(4,1))*ADIV1
      ALINEMOM10(N4,ID,2)=(CSUM2WSMOM(1,2)-CSUM2WSMOM(2,2)
     &+CSUM2WSMOM(3,2)-CSUM2WSMOM(4,2))*ADIV1
c      ALINEMOM10(N4,ID,3)=(CSUM2WSMOM(1,3)-CSUM2WSMOM(2,3)
c     &+CSUM2WSMOM(3,3)-CSUM2WSMOM(4,3))*ADIV1
c      ALINEMOM10(N4,ID,4)=(CSUM2WSMOM(1,4)-CSUM2WSMOM(2,4)
c     &+CSUM2WSMOM(3,4)-CSUM2WSMOM(4,4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE11(N4,ID)=(CSUMW(1)+CSUMW(2)+CSUMW(3)+CSUMW(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM11(N4,ID,1)=(CSUMWMOM(1,1)+CSUMWMOM(2,1)
     &+CSUMWMOM(3,1)+CSUMWMOM(4,1))*ADIV1
      ALINEMOM11(N4,ID,2)=(CSUMWMOM(1,2)+CSUMWMOM(2,2)
     &+CSUMWMOM(3,2)+CSUMWMOM(4,2))*ADIV1
c      ALINEMOM11(N4,ID,3)=(CSUMWMOM(1,3)+CSUMWMOM(2,3)
c     &+CSUMWMOM(3,3)+CSUMWMOM(4,3))*ADIV1
c      ALINEMOM11(N4,ID,4)=(CSUMWMOM(1,4)+CSUMWMOM(2,4)
c     &+CSUMWMOM(3,4)+CSUMWMOM(4,4))*ADIV1
C***********************************************************************
c     J=1, Pr=-, q=0
c***********************************************************************
      ALINE12(N4,ID)=(CSUMW(1)+GIOT*CSUMW(2)-CSUMW(3)
     &-GIOT*CSUMW(4))*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM12(N4,ID,1)=(CSUMWMOM(1,1)+GIOT*CSUMWMOM(2,1)
     &-CSUMWMOM(3,1)-GIOT*CSUMWMOM(4,1))*ADIV1
      ALINEMOM12(N4,ID,2)=(CSUMWMOM(1,2)+GIOT*CSUMWMOM(2,2)
     &-CSUMWMOM(3,2)-GIOT*CSUMWMOM(4,2))*ADIV1
c      ALINEMOM12(N4,ID,3)=(CSUMWMOM(1,3)+GIOT*CSUMWMOM(2,3)
c     &-CSUMWMOM(3,3)-GIOT*CSUMWMOM(4,3))*ADIV1
c      ALINEMOM12(N4,ID,4)=(CSUMWMOM(1,4)+GIOT*CSUMWMOM(2,4)
c     &-CSUMWMOM(3,4)-GIOT*CSUMWMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE13(N4,ID)=(CSUMW(1)-CSUMW(2)+CSUMW(3)-CSUMW(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM13(N4,ID,1)=(CSUMWMOM(1,1)-CSUMWMOM(2,1)
     &+CSUMWMOM(3,1)-CSUMWMOM(4,1))*ADIV1
      ALINEMOM13(N4,ID,2)=(CSUMWMOM(1,2)-CSUMWMOM(2,2)
     &+CSUMWMOM(3,2)-CSUMWMOM(4,2))*ADIV1
c      ALINEMOM13(N4,ID,3)=(CSUMWMOM(1,3)-CSUMWMOM(2,3)
c     &+CSUMWMOM(3,3)-CSUMWMOM(4,3))*ADIV1
c      ALINEMOM13(N4,ID,4)=(CSUMWMOM(1,4)-CSUMWMOM(2,4)
c     &+CSUMWMOM(3,4)-CSUMWMOM(4,4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE14(N4,ID)=(CSUM2W(1)+CSUM2W(2)+CSUM2W(3)+CSUM2W(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM14(N4,ID,1)=(CSUM2WMOM(1,1)+CSUM2WMOM(2,1)
     &+CSUM2WMOM(3,1)+CSUM2WMOM(4,1))*ADIV1
      ALINEMOM14(N4,ID,2)=(CSUM2WMOM(1,2)+CSUM2WMOM(2,2)
     &+CSUM2WMOM(3,2)+CSUM2WMOM(4,2))*ADIV1
c      ALINEMOM14(N4,ID,3)=(CSUM2WMOM(1,3)+CSUM2WMOM(2,3)
c     &+CSUM2WMOM(3,3)+CSUM2WMOM(4,3))*ADIV1
c      ALINEMOM14(N4,ID,4)=(CSUM2WMOM(1,4)+CSUM2WMOM(2,4)
c     &+CSUM2WMOM(3,4)+CSUM2WMOM(4,4))*ADIV1
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE15(N4,ID)=(CSUM2W(1)+GIOT*CSUM2W(2)-CSUM2W(3)
     &-GIOT*CSUM2W(4))*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM15(N4,ID,1)=(CSUM2WMOM(1,1)+GIOT*CSUM2WMOM(2,1)
     &-CSUM2WMOM(3,1)-GIOT*CSUM2WMOM(4,1))*ADIV1
      ALINEMOM15(N4,ID,2)=(CSUM2WMOM(1,2)+GIOT*CSUM2WMOM(2,2)
     &-CSUM2WMOM(3,2)-GIOT*CSUM2WMOM(4,2))*ADIV1
c      ALINEMOM15(N4,ID,3)=(CSUM2WMOM(1,3)+GIOT*CSUM2WMOM(2,3)
c     &-CSUM2WMOM(3,3)-GIOT*CSUM2WMOM(4,3))*ADIV1
c      ALINEMOM15(N4,ID,4)=(CSUM2WMOM(1,4)+GIOT*CSUM2WMOM(2,4)
c     &-CSUM2WMOM(3,4)-GIOT*CSUM2WMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE16(N4,ID)=(CSUM2W(1)-CSUM2W(2)+CSUM2W(3)-CSUM2W(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM16(N4,ID,1)=(CSUM2WMOM(1,1)-CSUM2WMOM(2,1)
     &+CSUM2WMOM(3,1)-CSUM2WMOM(4,1))*ADIV1
      ALINEMOM16(N4,ID,2)=(CSUM2WMOM(1,2)-CSUM2WMOM(2,2)
     &+CSUM2WMOM(3,2)-CSUM2WMOM(4,2))*ADIV1
c      ALINEMOM16(N4,ID,3)=(CSUM2WMOM(1,3)-CSUM2WMOM(2,3)
c     &+CSUM2WMOM(3,3)-CSUM2WMOM(4,3))*ADIV1
c      ALINEMOM16(N4,ID,4)=(CSUM2WMOM(1,4)-CSUM2WMOM(2,4)
c     &+CSUM2WMOM(3,4)-CSUM2WMOM(4,4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE17(N4,ID)=(CSUM3W(1)+CSUM3W(2)+CSUM3W(3)+CSUM3W(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM17(N4,ID,1)=(CSUM3WMOM(1,1)+CSUM3WMOM(2,1)
     &+CSUM3WMOM(3,1)+CSUM3WMOM(4,1))*ADIV1
      ALINEMOM17(N4,ID,2)=(CSUM3WMOM(1,2)+CSUM3WMOM(2,2)
     &+CSUM3WMOM(3,2)+CSUM3WMOM(4,2))*ADIV1
c      ALINEMOM17(N4,ID,3)=(CSUM3WMOM(1,3)+CSUM3WMOM(2,3)
c     &+CSUM3WMOM(3,3)+CSUM3WMOM(4,3))*ADIV1
c      ALINEMOM17(N4,ID,4)=(CSUM3WMOM(1,4)+CSUM3WMOM(2,4)
c     &+CSUM3WMOM(3,4)+CSUM3WMOM(4,4))*ADIV1
C***********************************************************************
c     J=1, Pr=+ q=0
c**********************************************************************
      ALINE18(N4,ID)=(CSUM3W(1)+GIOT*CSUM3W(2)-CSUM3W(3)-GIOT*CSUM3W(4))
     &*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM18(N4,ID,1)=(CSUM3WMOM(1,1)+GIOT*CSUM3WMOM(2,1)
     &-CSUM3WMOM(3,1)-GIOT*CSUM3WMOM(4,1))*ADIV1
      ALINEMOM18(N4,ID,2)=(CSUM3WMOM(1,2)+GIOT*CSUM3WMOM(2,2)
     &-CSUM3WMOM(3,2)-GIOT*CSUM3WMOM(4,2))*ADIV1
c      ALINEMOM18(N4,ID,3)=(CSUM3WMOM(1,3)+GIOT*CSUM3WMOM(2,3)
c     &-CSUM3WMOM(3,3)-GIOT*CSUM3WMOM(4,3))*ADIV1
c      ALINEMOM18(N4,ID,4)=(CSUM3WMOM(1,4)+GIOT*CSUM3WMOM(2,4)
c     &-CSUM3WMOM(3,4)-GIOT*CSUM3WMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE19(N4,ID)=(CSUM3W(1)-CSUM3W(2)+CSUM3W(3)-CSUM3W(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM19(N4,ID,1)=(CSUM3WMOM(1,1)-CSUM3WMOM(2,1)
     &+CSUM3WMOM(3,1)-CSUM3WMOM(4,1))*ADIV1
      ALINEMOM19(N4,ID,2)=(CSUM3WMOM(1,2)-CSUM3WMOM(2,2)
     &+CSUM3WMOM(3,2)-CSUM3WMOM(4,2))*ADIV1
c      ALINEMOM19(N4,ID,3)=(CSUM3WMOM(1,3)-CSUM3WMOM(2,3)
c     &+CSUM3WMOM(3,3)-CSUM3WMOM(4,3))*ADIV1
c      ALINEMOM19(N4,ID,4)=(CSUM3WMOM(1,4)-CSUM3WMOM(2,4)
c     &+CSUM3WMOM(3,4)-CSUM3WMOM(4,4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE20(N4,ID)=(CSUMUP(1)+CSUMUP(2)+CSUMUP(3)+CSUMUP(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM20(N4,ID,1)=(CSUMUPMOM(1,1)+CSUMUPMOM(2,1)
     &+CSUMUPMOM(3,1)+CSUMUPMOM(4,1))*ADIV1
      ALINEMOM20(N4,ID,2)=(CSUMUPMOM(1,2)+CSUMUPMOM(2,2)
     &+CSUMUPMOM(3,2)+CSUMUPMOM(4,2))*ADIV1
c      ALINEMOM20(N4,ID,3)=(CSUMUPMOM(1,3)+CSUMUPMOM(2,3)
c     &+CSUMUPMOM(3,3)+CSUMUPMOM(4,3))*ADIV1
c      ALINEMOM20(N4,ID,4)=(CSUMUPMOM(1,4)+CSUMUPMOM(2,4)
c     &+CSUMUPMOM(3,4)+CSUMUPMOM(4,4))*ADIV1
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE21(N4,ID)=(CSUMUP(1)+GIOT*CSUMUP(2)-CSUMUP(3)-GIOT*CSUMUP(4))
     &*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM21(N4,ID,1)=(CSUMUPMOM(1,1)+GIOT*CSUMUPMOM(2,1)
     &-CSUMUPMOM(3,1)-GIOT*CSUMUPMOM(4,1))*ADIV1
      ALINEMOM21(N4,ID,2)=(CSUMUPMOM(1,2)+GIOT*CSUMUPMOM(2,2)
     &-CSUMUPMOM(3,2)-GIOT*CSUMUPMOM(4,2))*ADIV1
c      ALINEMOM21(N4,ID,3)=(CSUMUPMOM(1,3)+GIOT*CSUMUPMOM(2,3)
c     &-CSUMUPMOM(3,3)-GIOT*CSUMUPMOM(4,3))*ADIV1
c      ALINEMOM21(N4,ID,4)=(CSUMUPMOM(1,4)+GIOT*CSUMUPMOM(2,4)
c     &-CSUMUPMOM(3,4)-GIOT*CSUMUPMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE22(N4,ID)=(CSUMUP(1)-CSUMUP(2)+CSUMUP(3)-CSUMUP(4))*ADIV1
C***********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c***********************************************************************
      ALINEMOM22(N4,ID,1)=(CSUMUPMOM(1,1)-CSUMUPMOM(2,1)
     &+CSUMUPMOM(3,1)-CSUMUPMOM(4,1))*ADIV1
      ALINEMOM22(N4,ID,2)=(CSUMUPMOM(1,2)-CSUMUPMOM(2,2)
     &+CSUMUPMOM(3,2)-CSUMUPMOM(4,2))*ADIV1
c      ALINEMOM22(N4,ID,3)=(CSUMUPMOM(1,3)-CSUMUPMOM(2,3)
c     &+CSUMUPMOM(3,3)-CSUMUPMOM(4,3))*ADIV1
c      ALINEMOM22(N4,ID,4)=(CSUMUPMOM(1,4)-CSUMUPMOM(2,4)
c     &+CSUMUPMOM(3,4)-CSUMUPMOM(4,4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE23(N4,ID)=(CSUMUD(1)+CSUMUD(2)+CSUMUD(3)+CSUMUD(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM23(N4,ID,1)=(CSUMUDMOM(1,1)+CSUMUDMOM(2,1)
     &+CSUMUDMOM(3,1)+CSUMUDMOM(4,1))*ADIV1
      ALINEMOM23(N4,ID,2)=(CSUMUDMOM(1,2)+CSUMUDMOM(2,2)
     &+CSUMUDMOM(3,2)+CSUMUDMOM(4,2))*ADIV1
c      ALINEMOM23(N4,ID,3)=(CSUMUDMOM(1,3)+CSUMUDMOM(2,3)
c     &+CSUMUDMOM(3,3)+CSUMUDMOM(4,3))*ADIV1
c      ALINEMOM23(N4,ID,4)=(CSUMUDMOM(1,4)+CSUMUDMOM(2,4)
c     &+CSUMUDMOM(3,4)+CSUMUDMOM(4,4))*ADIV1
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE24(N4,ID)=(CSUMUD(1)+GIOT*CSUMUD(2)-CSUMUD(3)
     &-GIOT*CSUMUD(4))*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c***********************************************************************
      ALINEMOM24(N4,ID,1)=(CSUMUDMOM(1,1)+GIOT*CSUMUDMOM(2,1)
     &-CSUMUDMOM(3,1)-GIOT*CSUMUDMOM(4,1))*ADIV1
      ALINEMOM24(N4,ID,2)=(CSUMUDMOM(1,2)+GIOT*CSUMUDMOM(2,2)
     &-CSUMUDMOM(3,2)-GIOT*CSUMUDMOM(4,2))*ADIV1
c      ALINEMOM24(N4,ID,3)=(CSUMUDMOM(1,3)+GIOT*CSUMUDMOM(2,3)
c     &-CSUMUDMOM(3,3)-GIOT*CSUMUDMOM(4,3))*ADIV1
c      ALINEMOM24(N4,ID,4)=(CSUMUDMOM(1,4)+GIOT*CSUMUDMOM(2,4)
c     &-CSUMUDMOM(3,4)-GIOT*CSUMUDMOM(4,4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE25(N4,ID)=(CSUMUD(1)-CSUMUD(2)+CSUMUD(3)-CSUMUD(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      ALINEMOM25(N4,ID,1)=(CSUMUDMOM(1,1)-CSUMUDMOM(2,1)
     &+CSUMUDMOM(3,1)-CSUMUDMOM(4,1))*ADIV1
      ALINEMOM25(N4,ID,2)=(CSUMUDMOM(1,2)-CSUMUDMOM(2,2)
     &+CSUMUDMOM(3,2)-CSUMUDMOM(4,2))*ADIV1
c      ALINEMOM25(N4,ID,3)=(CSUMUDMOM(1,3)-CSUMUDMOM(2,3)
c     &+CSUMUDMOM(3,3)-CSUMUDMOM(4,3))*ADIV1
c      ALINEMOM25(N4,ID,4)=(CSUMUDMOM(1,4)-CSUMUDMOM(2,4)
c     &+CSUMUDMOM(3,4)-CSUMUDMOM(4,4))*ADIV1
c**********************************************************************
C**********************************************************************
C**********************************************************************
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE26(N4,ID)=(CSUMTT1(1)+CSUMTT1(2)+CSUMTT1(3)+CSUMTT1(4)+
     &(CSUMTT1(7)+CSUMTT1(6)+CSUMTT1(5)+CSUMTT1(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM26(N4,ID,IK)=(CSUMTTMOM1(1,IK)+CSUMTTMOM1(2,IK)
     &+CSUMTTMOM1(3,IK)+CSUMTTMOM1(4,IK)+CSUMTTMOM1(7,IK)
     &+CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)
     &+CSUMTTMOM1(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE27(N4,ID)=(CSUMTT1(1)+CSUMTT1(2)+CSUMTT1(3)+CSUMTT1(4)-
     &(CSUMTT1(7)+CSUMTT1(6)+CSUMTT1(5)+CSUMTT1(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM27(N4,ID,IK)=(CSUMTTMOM1(1,IK)+CSUMTTMOM1(2,IK)
     &+CSUMTTMOM1(3,IK)+CSUMTTMOM1(4,IK)-(CSUMTTMOM1(7,IK)
     &+CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)
     &+CSUMTTMOM1(8,IK)))*ADIV2
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE28(N4,ID)=(CSUMTT1(1)+GIOT*CSUMTT1(2)-CSUMTT1(3)
     &-GIOT*CSUMTT1(4)+(CSUMTT1(6)+GIOT*CSUMTT1(7)-CSUMTT1(8)
     &-GIOT*CSUMTT1(5)))*ADIV2
C***********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM28(N4,ID,IK)=(CSUMTTMOM1(1,IK)+GIOT*CSUMTTMOM1(2,IK)
     &-CSUMTTMOM1(3,IK)-GIOT*CSUMTTMOM1(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE29(N4,ID)=(CSUMTT1(1)+GIOT*CSUMTT1(2)-CSUMTT1(3)
     &-GIOT*CSUMTT1(4)-(CSUMTT1(6)+GIOT*CSUMTT1(7)-CSUMTT1(8)
     &-GIOT*CSUMTT1(5)))*ADIV2
C***********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM29(N4,ID,IK)=(CSUMTTMOM1(6,IK)+GIOT*CSUMTTMOM1(7,IK)
     &-CSUMTTMOM1(8,IK)-GIOT*CSUMTTMOM1(5,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE30(N4,ID)=(CSUMTT1(1)-CSUMTT1(2)+CSUMTT1(3)-CSUMTT1(4)+
     &CSUMTT1(7)-CSUMTT1(6)+CSUMTT1(5)-CSUMTT1(8))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM30(N4,ID,IK)=(CSUMTTMOM1(1,IK)-CSUMTTMOM1(2,IK)
     &+CSUMTTMOM1(3,IK)-CSUMTTMOM1(4,IK)+CSUMTTMOM1(7,IK)
     &-CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)-CSUMTTMOM1(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE31(N4,ID)=(CSUMTT1(1)-CSUMTT1(2)+CSUMTT1(3)-CSUMTT1(4)-
     &(CSUMTT1(7)-CSUMTT1(6)+CSUMTT1(5)-CSUMTT1(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM31(N4,ID,IK)=(CSUMTTMOM1(1,IK)-CSUMTTMOM1(2,IK)
     &+CSUMTTMOM1(3,IK)-CSUMTTMOM1(4,IK)-(CSUMTTMOM1(7,IK)
     &-CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)-CSUMTTMOM1(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE32(N4,ID)=(CSUMTT2(1)+CSUMTT2(2)+CSUMTT2(3)+CSUMTT2(4)+
     &CSUMTT2(5)+CSUMTT2(6)+CSUMTT2(7)+CSUMTT2(8))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM32(N4,ID,IK)=(CSUMTTMOM2(1,IK)+CSUMTTMOM2(2,IK)
     &+CSUMTTMOM2(3,IK)+CSUMTTMOM2(4,IK)+CSUMTTMOM2(5,IK)
     &+CSUMTTMOM2(6,IK)+CSUMTTMOM2(7,IK)
     &+CSUMTTMOM2(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE33(N4,ID)=(CSUMTT2(1)+CSUMTT2(2)+CSUMTT2(3)+CSUMTT2(4)-
     &(CSUMTT2(5)+CSUMTT2(6)+CSUMTT2(7)+CSUMTT2(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q
c**********************************************************************
      DO IK=1,2
      ALINEMOM33(N4,ID,IK)=(CSUMTTMOM2(1,IK)+CSUMTTMOM2(2,IK)
     &+CSUMTTMOM2(3,IK)+CSUMTTMOM2(4,IK)-(CSUMTTMOM2(5,IK)
     &+CSUMTTMOM2(6,IK)+CSUMTTMOM2(7,IK)
     &+CSUMTTMOM2(8,IK)))*ADIV2
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE34(N4,ID)=(CSUMTT2(1)+GIOT*CSUMTT2(2)-CSUMTT2(3)
     &-GIOT*CSUMTT2(4)+(CSUMTT2(7)+GIOT*CSUMTT2(8)-CSUMTT2(5)
     &-GIOT*CSUMTT2(6)))*ADIV2
C***********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM34(N4,ID,IK)=(CSUMTTMOM2(1,IK)+GIOT*CSUMTTMOM2(2,IK)
     &-CSUMTTMOM2(3,IK)-GIOT*CSUMTTMOM2(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE35(N4,ID)=(CSUMTT2(1)+GIOT*CSUMTT2(2)-CSUMTT2(3)
     &-GIOT*CSUMTT2(4)-(CSUMTT2(7)+GIOT*CSUMTT2(8)-CSUMTT2(5)
     &-GIOT*CSUMTT2(6)))*ADIV2
C***********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM35(N4,ID,IK)=(CSUMTTMOM2(7,IK)+GIOT*CSUMTTMOM2(8,IK)
     &-CSUMTTMOM2(5,IK)-GIOT*CSUMTTMOM2(6,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE36(N4,ID)=(CSUMTT2(1)-CSUMTT2(2)+CSUMTT2(3)-CSUMTT2(4)+
     &CSUMTT2(6)-CSUMTT2(5)+CSUMTT2(8)-CSUMTT2(7))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM36(N4,ID,IK)=(CSUMTTMOM2(1,IK)-CSUMTTMOM2(2,IK)
     &+CSUMTTMOM2(3,IK)-CSUMTTMOM2(4,IK)+CSUMTTMOM2(6,IK)
     &-CSUMTTMOM2(5,IK)+CSUMTTMOM2(8,IK)
     &-CSUMTTMOM2(7,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE37(N4,ID)=(CSUMTT2(1)-CSUMTT2(2)+CSUMTT2(3)-CSUMTT2(4)-
     &(CSUMTT2(6)-CSUMTT2(5)+CSUMTT2(8)-CSUMTT2(7)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM37(N4,ID,IK)=(CSUMTTMOM2(1,IK)-CSUMTTMOM2(2,IK)
     &+CSUMTTMOM2(3,IK)-CSUMTTMOM2(4,IK)-(CSUMTTMOM2(6,IK)
     &-CSUMTTMOM2(5,IK)+CSUMTTMOM2(8,IK)
     &-CSUMTTMOM2(7,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE38(N4,ID)=(CSUMTT3(1)+CSUMTT3(2)+CSUMTT3(3)+CSUMTT3(4)+
     &CSUMTT3(5)+CSUMTT3(6)+CSUMTT3(7)+CSUMTT3(8))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM38(N4,ID,IK)=(CSUMTTMOM3(1,IK)+CSUMTTMOM3(2,IK)
     &+CSUMTTMOM3(3,IK)+CSUMTTMOM3(4,IK)+CSUMTTMOM3(5,IK)
     &+CSUMTTMOM3(6,IK)+CSUMTTMOM3(7,IK)
     &+CSUMTTMOM3(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE39(N4,ID)=(CSUMTT3(1)+CSUMTT3(2)+CSUMTT3(3)+CSUMTT3(4)-
     &(CSUMTT3(5)+CSUMTT3(6)+CSUMTT3(7)+CSUMTT3(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM39(N4,ID,IK)=(CSUMTTMOM3(1,IK)+CSUMTTMOM3(2,IK)
     &+CSUMTTMOM3(3,IK)+CSUMTTMOM3(4,IK)-(CSUMTTMOM3(5,IK)
     &+CSUMTTMOM3(6,IK)+CSUMTTMOM3(7,IK)
     &+CSUMTTMOM3(8,IK)))*ADIV2
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE40(N4,ID)=(CSUMTT3(1)+GIOT*CSUMTT3(2)-CSUMTT3(3)
     &-GIOT*CSUMTT3(4)+(CSUMTT3(6)+GIOT*CSUMTT3(7)-CSUMTT3(8)
     &-GIOT*CSUMTT3(5)))*ADIV2
C***********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM40(N4,ID,IK)=(CSUMTTMOM3(1,IK)+GIOT*CSUMTTMOM3(2,IK)
     &-CSUMTTMOM3(3,IK)-GIOT*CSUMTTMOM3(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE41(N4,ID)=(CSUMTT3(1)+GIOT*CSUMTT3(2)-CSUMTT3(3)
     &-GIOT*CSUMTT3(4)-(CSUMTT3(6)+GIOT*CSUMTT3(7)-CSUMTT3(8)
     &-GIOT*CSUMTT3(5)))*ADIV2
C***********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM41(N4,ID,IK)=(CSUMTTMOM3(6,IK)+GIOT*CSUMTTMOM3(7,IK)
     &-CSUMTTMOM3(8,IK)-GIOT*CSUMTTMOM3(5,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE42(N4,ID)=(CSUMTT3(1)-CSUMTT3(2)+CSUMTT3(3)-CSUMTT3(4)+
     &(CSUMTT3(7)-CSUMTT3(6)+CSUMTT3(5)-CSUMTT3(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM42(N4,ID,IK)=(CSUMTTMOM3(1,IK)-CSUMTTMOM3(2,IK)
     &+CSUMTTMOM3(3,IK)-CSUMTTMOM3(4,IK)+(CSUMTTMOM3(7,IK)
     &-CSUMTTMOM3(6,IK)+CSUMTTMOM3(5,IK)-CSUMTTMOM3(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE43(N4,ID)=(CSUMTT3(1)-CSUMTT3(2)+CSUMTT3(3)-CSUMTT3(4)-
     &(CSUMTT3(7)-CSUMTT3(6)+CSUMTT3(5)-CSUMTT3(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM43(N4,ID,IK)=(CSUMTTMOM3(1,IK)-CSUMTTMOM3(2,IK)
     &+CSUMTTMOM3(3,IK)-CSUMTTMOM3(4,IK)-(CSUMTTMOM3(7,IK)
     &-CSUMTTMOM3(6,IK)+CSUMTTMOM3(5,IK)-CSUMTTMOM3(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE44(N4,ID)=(CSUMTT4(1)+CSUMTT4(2)+CSUMTT4(3)+CSUMTT4(4)+
     &CSUMTT4(5)+CSUMTT4(6)+CSUMTT4(7)+CSUMTT4(8))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM44(N4,ID,IK)=(CSUMTTMOM4(1,IK)+CSUMTTMOM4(2,IK)
     &+CSUMTTMOM4(3,IK)+CSUMTTMOM4(4,IK)+CSUMTTMOM4(5,IK)
     &+CSUMTTMOM4(6,IK)+CSUMTTMOM4(7,IK)
     &+CSUMTTMOM4(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE45(N4,ID)=(CSUMTT4(1)+CSUMTT4(2)+CSUMTT4(3)+CSUMTT4(4)-
     &(CSUMTT4(5)+CSUMTT4(6)+CSUMTT4(7)+CSUMTT4(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM45(N4,ID,IK)=(CSUMTTMOM4(1,IK)+CSUMTTMOM4(2,IK)
     &+CSUMTTMOM4(3,IK)+CSUMTTMOM4(4,IK)-(CSUMTTMOM4(5,IK)
     &+CSUMTTMOM4(6,IK)+CSUMTTMOM4(7,IK)
     &+CSUMTTMOM4(8,IK)))*ADIV2
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE46(N4,ID)=(CSUMTT4(1)+GIOT*CSUMTT4(2)-CSUMTT4(3)
     &-GIOT*CSUMTT4(4)+(CSUMTT4(7)+GIOT*CSUMTT4(8)-CSUMTT4(5)
     &-GIOT*CSUMTT4(6)))*ADIV2
C***********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM46(N4,ID,IK)=(CSUMTTMOM4(1,IK)+GIOT*CSUMTTMOM4(2,IK)
     &-CSUMTTMOM4(3,IK)-GIOT*CSUMTTMOM4(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE47(N4,ID)=(CSUMTT4(1)+GIOT*CSUMTT4(2)-CSUMTT4(3)
     &-GIOT*CSUMTT4(4)-(CSUMTT4(7)+GIOT*CSUMTT4(8)-CSUMTT4(5)
     &-GIOT*CSUMTT4(6)))*ADIV2
C***********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM47(N4,ID,IK)=(CSUMTTMOM4(7,IK)+GIOT*CSUMTTMOM4(8,IK)
     &-CSUMTTMOM4(5,IK)-GIOT*CSUMTTMOM4(6,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE48(N4,ID)=(CSUMTT4(1)-CSUMTT4(2)+CSUMTT4(3)-CSUMTT4(4)+
     &CSUMTT4(8)-CSUMTT4(7)+CSUMTT4(6)-CSUMTT4(5))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM48(N4,ID,IK)=(CSUMTTMOM4(1,IK)-CSUMTTMOM4(2,IK)
     &+CSUMTTMOM4(3,IK)-CSUMTTMOM4(4,IK)+CSUMTTMOM4(8,IK)
     &-CSUMTTMOM4(7,IK)+CSUMTTMOM4(6,IK)-CSUMTTMOM4(5,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE49(N4,ID)=(CSUMTT4(1)-CSUMTT4(2)+CSUMTT4(3)-CSUMTT4(4)-
     &(CSUMTT4(8)-CSUMTT4(7)+CSUMTT4(6)-CSUMTT4(5)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM49(N4,ID,IK)=(CSUMTTMOM4(1,IK)-CSUMTTMOM4(2,IK)
     &+CSUMTTMOM4(3,IK)-CSUMTTMOM4(4,IK)-(CSUMTTMOM4(8,IK)
     &-CSUMTTMOM4(7,IK)+CSUMTTMOM4(6,IK)-CSUMTTMOM4(5,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-5 OPERATORS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE50(N4,ID)=(CSUMTT5(1)+CSUMTT5(2)+CSUMTT5(3)+CSUMTT5(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM50(N4,ID,IK)=(CSUMTTMOM5(1,IK)+CSUMTTMOM5(2,IK)
     &+CSUMTTMOM5(3,IK)+CSUMTTMOM5(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE51(N4,ID)=(CSUMTT5(1)+GIOT*CSUMTT5(2)-CSUMTT5(3)
     &-GIOT*CSUMTT5(4))*ADIV1
C***********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM51(N4,ID,IK)=(CSUMTTMOM5(1,IK)+GIOT*CSUMTTMOM5(2,IK)
     &-CSUMTTMOM5(3,IK)-GIOT*CSUMTTMOM5(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE52(N4,ID)=(CSUMTT5(1)-CSUMTT5(2)+CSUMTT5(3)-CSUMTT5(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM52(N4,ID,IK)=(CSUMTTMOM5(1,IK)-CSUMTTMOM5(2,IK)
     &+CSUMTTMOM5(3,IK)-CSUMTTMOM5(4,IK))*ADIV1
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 6
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-6 OPERATORS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE53(N4,ID)=(CSUMTT6(1)+CSUMTT6(2)+CSUMTT6(3)+CSUMTT6(4))*ADIV1
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM53(N4,ID,IK)=(CSUMTTMOM6(1,IK)+CSUMTTMOM6(2,IK)
     &+CSUMTTMOM6(3,IK)+CSUMTTMOM6(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE54(N4,ID)=(CSUMTT6(1)+GIOT*CSUMTT6(2)-CSUMTT6(3)
     &-GIOT*CSUMTT6(4))*ADIV1
C***********************************************************************
c     J=1, q
c**********************************************************************
      DO IK=1,2
      ALINEMOM54(N4,ID,IK)=(CSUMTTMOM6(1,IK)+GIOT*CSUMTTMOM6(2,IK)
     &-CSUMTTMOM6(3,IK)-GIOT*CSUMTTMOM6(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE55(N4,ID)=(CSUMTT6(1)-CSUMTT6(2)+CSUMTT6(3)-CSUMTT6(4))*ADIV1
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM55(N4,ID,IK)=(CSUMTTMOM6(1,IK)-CSUMTTMOM6(2,IK)
     &+CSUMTTMOM6(3,IK)-CSUMTTMOM6(4,IK))*ADIV1
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 7
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     TT-7 OPERATORS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE56(N4,ID)=(CSUMTT7(1)+CSUMTT7(2)+CSUMTT7(3)+CSUMTT7(4)+
     &(CSUMTT7(5)+CSUMTT7(6)+CSUMTT7(7)+CSUMTT7(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM56(N4,ID,IK)=(CSUMTTMOM7(1,IK)+CSUMTTMOM7(2,IK)
     &+CSUMTTMOM7(3,IK)+CSUMTTMOM7(4,IK)+(CSUMTTMOM7(5,IK)
     &+CSUMTTMOM7(6,IK)+CSUMTTMOM7(7,IK)
     &+CSUMTTMOM7(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE57(N4,ID)=(CSUMTT7(1)+CSUMTT7(2)+CSUMTT7(3)+CSUMTT7(4)-
     &(CSUMTT7(5)+CSUMTT7(6)+CSUMTT7(7)+CSUMTT7(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM57(N4,ID,IK)=(CSUMTTMOM7(1,IK)+CSUMTTMOM7(2,IK)
     &+CSUMTTMOM7(3,IK)+CSUMTTMOM7(4,IK)-(CSUMTTMOM7(5,IK)
     &+CSUMTTMOM7(6,IK)+CSUMTTMOM7(7,IK)
     &+CSUMTTMOM7(8,IK)))*ADIV2
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE58(N4,ID)=(CSUMTT7(1)+GIOT*CSUMTT7(2)-CSUMTT7(3)
     &-GIOT*CSUMTT7(4)+
     &(CSUMTT7(7)+GIOT*CSUMTT7(8)-CSUMTT7(5)-GIOT*CSUMTT7(6)))*ADIV2
C***********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM58(N4,ID,IK)=(CSUMTTMOM7(1,IK)+GIOT*CSUMTTMOM7(2,IK)
     &-CSUMTTMOM7(3,IK)-GIOT*CSUMTTMOM7(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE59(N4,ID)=(CSUMTT7(1)+GIOT*CSUMTT7(2)-CSUMTT7(3)
     &-GIOT*CSUMTT7(4)-
     &(CSUMTT7(7)+GIOT*CSUMTT7(8)-CSUMTT7(5)-GIOT*CSUMTT7(6)))*ADIV2
C***********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM59(N4,ID,IK)=(CSUMTTMOM7(7,IK)+GIOT*CSUMTTMOM7(8,IK)
     &-CSUMTTMOM7(5,IK)-GIOT*CSUMTTMOM7(6,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE60(N4,ID)=(CSUMTT7(1)-CSUMTT7(2)+CSUMTT7(3)-CSUMTT7(4)+
     &(CSUMTT7(8)-CSUMTT7(7)+CSUMTT7(6)-CSUMTT7(5)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM60(N4,ID,IK)=(CSUMTTMOM7(1,IK)-CSUMTTMOM7(2,IK)
     &+CSUMTTMOM7(3,IK)-CSUMTTMOM7(4,IK)+(CSUMTTMOM7(8,IK)
     &-CSUMTTMOM7(7,IK)+CSUMTTMOM7(6,IK)
     &-CSUMTTMOM7(5,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE61(N4,ID)=(CSUMTT7(1)-CSUMTT7(2)+CSUMTT7(3)-CSUMTT7(4)-
     &(CSUMTT7(8)-CSUMTT7(7)+CSUMTT7(6)-CSUMTT7(5)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM61(N4,ID,IK)=(CSUMTTMOM7(1,IK)-CSUMTTMOM7(2,IK)
     &+CSUMTTMOM7(3,IK)-CSUMTTMOM7(4,IK)-(CSUMTTMOM7(8,IK)
     &-CSUMTTMOM7(7,IK)+CSUMTTMOM7(6,IK)
     &-CSUMTTMOM7(5,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 8
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-8 OPERATORS CP=+,Pr=+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE62(N4,ID)=(CSUMTT8(1)+CSUMTT8(2)+CSUMTT8(3)+CSUMTT8(4)
     &+(CSUMTT8(5)+CSUMTT8(6)+CSUMTT8(7)+CSUMTT8(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM62(N4,ID,IK)=(CSUMTTMOM8(1,IK)+CSUMTTMOM8(2,IK)
     &+CSUMTTMOM8(3,IK)+CSUMTTMOM8(4,IK)+(CSUMTTMOM8(5,IK)
     &+CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)
     &+CSUMTTMOM8(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE63(N4,ID)=(CSUMTT8(1)+CSUMTT8(2)+CSUMTT8(3)+CSUMTT8(4)
     &-(CSUMTT8(5)+CSUMTT8(6)+CSUMTT8(7)+CSUMTT8(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM63(N4,ID,IK)=(CSUMTTMOM8(1,IK)+CSUMTTMOM8(2,IK)
     &+CSUMTTMOM8(3,IK)+CSUMTTMOM8(4,IK)-(CSUMTTMOM8(5,IK)
     &+CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)
     &+CSUMTTMOM8(8,IK)))*ADIV2
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE64(N4,ID)=(CSUMTT8(1)+GIOT*CSUMTT8(2)-CSUMTT8(3)
     &-GIOT*CSUMTT8(4))*ADIV1
C***********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
         ALINEMOM64(N4,ID,IK)=(CSUMTTMOM8(1,IK)+GIOT*CSUMTTMOM8(2,IK)
     &-CSUMTTMOM8(3,IK)-GIOT*CSUMTTMOM8(4,IK))*ADIV1
      ENDDO
C***********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE65(N4,ID)=(CSUMTT8(7)+GIOT*CSUMTT8(6)-CSUMTT8(5)
     &-GIOT*CSUMTT8(8))*ADIV1
C***********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
         ALINEMOM65(N4,ID,IK)=(CSUMTTMOM8(7,IK)+GIOT*CSUMTTMOM8(6,IK)
     &-CSUMTTMOM8(5,IK)-GIOT*CSUMTTMOM8(8,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE66(N4,ID)=(CSUMTT8(1)-CSUMTT8(2)+CSUMTT8(3)-CSUMTT8(4)
     &+(CSUMTT8(5)-CSUMTT8(6)+CSUMTT8(7)-CSUMTT8(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM66(N4,ID,IK)=(CSUMTTMOM8(1,IK)-CSUMTTMOM8(2,IK)
     &+CSUMTTMOM8(3,IK)-CSUMTTMOM8(4,IK)+(CSUMTTMOM8(5,IK)
     &-CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)
     &-CSUMTTMOM8(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE67(N4,ID)=(CSUMTT8(1)-CSUMTT8(2)+CSUMTT8(3)-CSUMTT8(4)
     &-(CSUMTT8(5)-CSUMTT8(6)+CSUMTT8(7)-CSUMTT8(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM67(N4,ID,IK)=(CSUMTTMOM8(1,IK)-CSUMTTMOM8(2,IK)
     &+CSUMTTMOM8(3,IK)-CSUMTTMOM8(4,IK)-(CSUMTTMOM8(5,IK)
     &-CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK)
     &-CSUMTTMOM8(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 9
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-9 OPERATORS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE68(N4,ID)=(CSUMTT9(1)+CSUMTT9(2)+CSUMTT9(3)+CSUMTT9(4)+
     &(CSUMTT9(5)+CSUMTT9(6)+CSUMTT9(7)+CSUMTT9(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1,2
         ALINEMOM68(N4,ID,IK)=(CSUMTTMOM9(1,IK)+CSUMTTMOM9(2,IK)+
     &CSUMTTMOM9(3,IK)+CSUMTTMOM9(4,IK)
     &+(CSUMTTMOM9(5,IK)+CSUMTTMOM9(6,IK)+CSUMTTMOM9(7,IK)
     &+CSUMTTMOM9(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-9 OPERATORS CP=-,Pz=-,J=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE69(N4,ID)=(CSUMTT9(1)+CSUMTT9(2)+CSUMTT9(3)+CSUMTT9(4)
     &-(CSUMTT9(5)+CSUMTT9(6)+CSUMTT9(7)+CSUMTT9(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1,2
         ALINEMOM69(N4,ID,IK)=(CSUMTTMOM9(1,IK)+CSUMTTMOM9(2,IK)+
     &CSUMTTMOM9(3,IK)+CSUMTTMOM9(4,IK)
     &-(CSUMTTMOM9(5,IK)+CSUMTTMOM9(6,IK)+CSUMTTMOM9(7,IK)
     &+CSUMTTMOM9(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************      
      ALINE70(N4,ID)=(CSUMTT9(1)+GIOT*CSUMTT9(2)-CSUMTT9(3)
     &-GIOT*CSUMTT9(4)+(CSUMTT9(7)+GIOT*CSUMTT9(8)-CSUMTT9(5)
     &-GIOT*CSUMTT9(6)))*ADIV2
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************      
      DO IK=1,2
         ALINEMOM70(N4,ID,IK)=(CSUMTTMOM9(1,IK)+GIOT*CSUMTTMOM9(2,IK)
     &-CSUMTTMOM9(3,IK)-GIOT*CSUMTTMOM9(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************            
      ALINE71(N4,ID)=(CSUMTT9(1)+GIOT*CSUMTT9(2)-CSUMTT9(3)
     &-GIOT*CSUMTT9(4)-(CSUMTT9(7)+GIOT*CSUMTT9(8)-CSUMTT9(5)
     &-GIOT*CSUMTT9(6)))*ADIV2
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************            
      DO IK=1,2
         ALINEMOM71(N4,ID,IK)=(CSUMTTMOM9(7,IK)+GIOT*CSUMTTMOM9(8,IK)
     &-CSUMTTMOM9(5,IK)-GIOT*CSUMTTMOM9(6,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************            
      ALINE72(N4,ID)=(CSUMTT9(1)-CSUMTT9(2)+CSUMTT9(3)-CSUMTT9(4)
     &+(CSUMTT9(7)-CSUMTT9(6)+CSUMTT9(5)-CSUMTT9(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2
c**********************************************************************            
      DO IK=1,2
         ALINEMOM72(N4,ID,IK)=(CSUMTTMOM9(1,IK)-CSUMTTMOM9(2,IK)
     &+CSUMTTMOM9(3,IK)-CSUMTTMOM9(4,IK)
     &+(CSUMTTMOM9(7,IK)-CSUMTTMOM9(6,IK)+CSUMTTMOM9(5,IK)
     &-CSUMTTMOM9(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=-, q=0 !Here!
c**********************************************************************                  
      ALINE73(N4,ID)=(CSUMTT9(1)-CSUMTT9(2)+CSUMTT9(3)-CSUMTT9(4)
     &-(CSUMTT9(7)-CSUMTT9(6)+CSUMTT9(5)-CSUMTT9(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2
c**********************************************************************                  
      DO IK=1,2
         ALINEMOM73(N4,ID,IK)=(CSUMTTMOM9(1,IK)-CSUMTTMOM9(2,IK)
     &+CSUMTTMOM9(3,IK)-CSUMTTMOM9(4,IK)
     &-(CSUMTTMOM9(7,IK)-CSUMTTMOM9(6,IK)+CSUMTTMOM9(5,IK)
     &-CSUMTTMOM9(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-10 OPERATORS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************                  
      ALINE74(N4,ID)=(CSUMTT10(1)+CSUMTT10(2)+CSUMTT10(3)+CSUMTT10(4)
     &+(CSUMTT10(5)+CSUMTT10(6)+CSUMTT10(7)+CSUMTT10(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2
c**********************************************************************                  
      DO IK=1,2
         ALINEMOM74(N4,ID,IK)=(CSUMTTMOM10(1,IK)+CSUMTTMOM10(2,IK)
     &+CSUMTTMOM10(3,IK)+CSUMTTMOM10(4,IK)
     &+(CSUMTTMOM10(5,IK)+CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)
     &+CSUMTTMOM10(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************                        
      ALINE75(N4,ID)=(CSUMTT10(1)+CSUMTT10(2)+CSUMTT10(3)+CSUMTT10(4)
     &-(CSUMTT10(5)+CSUMTT10(6)+CSUMTT10(7)+CSUMTT10(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2
c**********************************************************************                        
      DO IK=1,2
         ALINEMOM75(N4,ID,IK)=(CSUMTTMOM10(1,IK)+CSUMTTMOM10(2,IK)
     &+CSUMTTMOM10(3,IK)+CSUMTTMOM10(4,IK)
     &-(CSUMTTMOM10(5,IK)+CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)
     &+CSUMTTMOM10(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************                  
      ALINE76(N4,ID)=(CSUMTT10(1)+GIOT*CSUMTT10(2)-CSUMTT10(3)
     &-GIOT*CSUMTT10(4)+(CSUMTT10(7)+GIOT*CSUMTT10(6)-CSUMTT10(5)
     &-GIOT*CSUMTT10(8)))*ADIV2
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1,2
         ALINEMOM76(N4,ID,IK)=(CSUMTTMOM10(1,IK)+GIOT*CSUMTTMOM10(2,IK)
     &-CSUMTTMOM10(3,IK)-GIOT*CSUMTTMOM10(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE77(N4,ID)=(CSUMTT10(1)+GIOT*CSUMTT10(2)-CSUMTT10(3)
     &-GIOT*CSUMTT10(4)-(CSUMTT10(7)+GIOT*CSUMTT10(6)-CSUMTT10(5)
     &-GIOT*CSUMTT10(8)))*ADIV2
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************
      DO IK=1,2
         ALINEMOM77(N4,ID,IK)=(CSUMTTMOM10(7,IK)+GIOT*CSUMTTMOM10(6,IK)
     &-CSUMTTMOM10(5,IK)-GIOT*CSUMTTMOM10(8,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE78(N4,ID)=(CSUMTT10(1)-CSUMTT10(2)+CSUMTT10(3)-CSUMTT10(4)
     &+(CSUMTT10(5)-CSUMTT10(6)+CSUMTT10(7)-CSUMTT10(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2
c**********************************************************************
      DO IK=1,2
         ALINEMOM78(N4,ID,IK)=(CSUMTTMOM10(1,IK)-CSUMTTMOM10(2,IK)
     &+CSUMTTMOM10(3,IK)-CSUMTTMOM10(4,IK)
     &+(CSUMTTMOM10(5,IK)-CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)
     &-CSUMTTMOM10(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************      
      ALINE79(N4,ID)=(CSUMTT10(1)-CSUMTT10(2)+CSUMTT10(3)-CSUMTT10(4)
     &-(CSUMTT10(5)-CSUMTT10(6)+CSUMTT10(7)-CSUMTT10(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2
c**********************************************************************      
      DO IK=1,2
         ALINEMOM79(N4,ID,IK)=(CSUMTTMOM10(1,IK)-CSUMTTMOM10(2,IK)
     &+CSUMTTMOM10(3,IK)-CSUMTTMOM10(4,IK)
     &-(CSUMTTMOM10(5,IK)-CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK)
     &-CSUMTTMOM10(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-11 OPERATORS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************            
      ALINE80(N4,ID)=(CSUMTT11(1)+CSUMTT11(2)+CSUMTT11(3)+CSUMTT11(4)
     &+(CSUMTT11(5)+CSUMTT11(6)+CSUMTT11(7)+CSUMTT11(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2
c**********************************************************************            
      DO IK=1,2
         ALINEMOM80(N4,ID,IK)=(CSUMTTMOM11(1,IK)+CSUMTTMOM11(2,IK)
     &+CSUMTTMOM11(3,IK)+CSUMTTMOM11(4,IK)
     &+(CSUMTTMOM11(5,IK)+CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)
     &+CSUMTTMOM11(8,IK)))*ADIV2
      ENDDO      
c**********************************************************************
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************            
      ALINE81(N4,ID)=(CSUMTT11(1)+CSUMTT11(2)+CSUMTT11(3)+CSUMTT11(4)
     &-(CSUMTT11(5)+CSUMTT11(6)+CSUMTT11(7)+CSUMTT11(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2
c**********************************************************************            
      DO IK=1,2
         ALINEMOM81(N4,ID,IK)=(CSUMTTMOM11(1,IK)+CSUMTTMOM11(2,IK)
     &+CSUMTTMOM11(3,IK)+CSUMTTMOM11(4,IK)
     &-(CSUMTTMOM11(5,IK)+CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)
     &+CSUMTTMOM11(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE82(N4,ID)=(CSUMTT11(1)+GIOT*CSUMTT11(2)-CSUMTT11(3)
     &-GIOT*CSUMTT11(4))*ADIV1
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1,2
         ALINEMOM82(N4,ID,IK)=(CSUMTTMOM11(1,IK)+GIOT*CSUMTTMOM11(2,IK)
     &-CSUMTTMOM11(3,IK)-GIOT*CSUMTTMOM11(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************                        
      ALINE83(N4,ID)=(CSUMTT11(7)+GIOT*CSUMTT11(6)-CSUMTT11(5)
     &-GIOT*CSUMTT11(8))*ADIV1
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1,2
         ALINEMOM83(N4,ID,IK)=(CSUMTTMOM11(7,IK)+GIOT*CSUMTTMOM11(6,IK)
     &-CSUMTTMOM11(5,IK)-GIOT*CSUMTTMOM11(8,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE84(N4,ID)=(CSUMTT11(1)-CSUMTT11(2)+CSUMTT11(3)-CSUMTT11(4)
     &+(CSUMTT11(5)-CSUMTT11(6)+CSUMTT11(7)-CSUMTT11(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2
c**********************************************************************                  
      DO IK=1,2
         ALINEMOM84(N4,ID,IK)=(CSUMTTMOM11(1,IK)-CSUMTTMOM11(2,IK)
     &+CSUMTTMOM11(3,IK)-CSUMTTMOM11(4,IK)
     &+(CSUMTTMOM11(5,IK)-CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)
     &-CSUMTTMOM11(8,IK)))*ADIV2
      ENDDO      
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************                        
      ALINE85(N4,ID)=(CSUMTT11(1)-CSUMTT11(2)+CSUMTT11(3)-CSUMTT11(4)
     &-(CSUMTT11(5)-CSUMTT11(6)+CSUMTT11(7)-CSUMTT11(8)))*ADIV2
Cc**********************************************************************
c     J=2, Pp=-, q=1,2
c**********************************************************************                        
      DO IK=1,2
         ALINEMOM85(N4,ID,IK)=(CSUMTTMOM11(1,IK)-CSUMTTMOM11(2,IK)
     &+CSUMTTMOM11(3,IK)-CSUMTTMOM11(4,IK)
     &-(CSUMTTMOM11(5,IK)-CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK)
     &-CSUMTTMOM11(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     12 THIS!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-12 OPERATORS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE86(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)
     &+(CSUMTT12(5)+CSUMTT12(6)+CSUMTT12(7)+CSUMTT12(8))
     &+(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))
     &+(CSUMTT12(13)+CSUMTT12(14)+CSUMTT12(15)+CSUMTT12(16)))*ADIV3
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM86(N4,ID,IK)=(CSUMTTMOM12(1,IK)+CSUMTTMOM12(2,IK)
     &+CSUMTTMOM12(3,IK)+CSUMTTMOM12(4,IK)
     &+(CSUMTTMOM12(9,IK)+CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)
     &+CSUMTTMOM12(12,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE87(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)
     &+(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))
     &-(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6))
     &-(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM87(N4,ID,IK)=(CSUMTTMOM12(7,IK)+CSUMTTMOM12(8,IK)
     &+CSUMTTMOM12(5,IK)+CSUMTTMOM12(6,IK)+(CSUMTTMOM12(15,IK)
     &+CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)+CSUMTTMOM12(14,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE88(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)
     &-(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))
     &+(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6))
     &-(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM88(N4,ID,IK)=(CSUMTTMOM12(1,IK)+CSUMTTMOM12(2,IK)
     &+CSUMTTMOM12(3,IK)+CSUMTTMOM12(4,IK)-(CSUMTTMOM12(9,IK)
     &+CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)+CSUMTTMOM12(12,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE89(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4)
     &-(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12))
     &-(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6))
     &+(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM89(N4,ID,IK)=(CSUMTTMOM12(7,IK)+CSUMTTMOM12(8,IK)
     &+CSUMTTMOM12(5,IK)+CSUMTTMOM12(6,IK)-(CSUMTTMOM12(15,IK)
     &+CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)+CSUMTTMOM12(14,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE90(N4,ID)=(CSUMTT12(1)+GIOT*CSUMTT12(2)-CSUMTT12(3)
     &-GIOT*CSUMTT12(4)
     &+(CSUMTT12(7)+GIOT*CSUMTT12(8)-CSUMTT12(5)-GIOT*CSUMTT12(6))
     &)*ADIV2
c**********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM90(N4,ID,IK)=(CSUMTTMOM12(1,IK)+GIOT*CSUMTTMOM12(2,IK)
     &-CSUMTTMOM12(3,IK)-GIOT*CSUMTTMOM12(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE91(N4,ID)=(CSUMTT12(9)+GIOT*CSUMTT12(10)-CSUMTT12(11)
     &-GIOT*CSUMTT12(12)
     &+(CSUMTT12(15)+GIOT*CSUMTT12(16)-CSUMTT12(13)-GIOT*CSUMTT12(14))
     &)*ADIV2
c**********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM91(N4,ID,IK)=(CSUMTTMOM12(9,IK)+GIOT*CSUMTTMOM12(10,IK)
     &-CSUMTTMOM12(11,IK)-GIOT*CSUMTTMOM12(12,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE92(N4,ID)=(CSUMTT12(1)+GIOT*CSUMTT12(2)-CSUMTT12(3)
     &-GIOT*CSUMTT12(4)
     &-(CSUMTT12(7)+GIOT*CSUMTT12(8)-CSUMTT12(5)-GIOT*CSUMTT12(6))
     &)*ADIV2
c**********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM92(N4,ID,IK)=(CSUMTTMOM12(7,IK)+GIOT*CSUMTTMOM12(8,IK)
     &-CSUMTTMOM12(5,IK)-GIOT*CSUMTTMOM12(6,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE93(N4,ID)=(CSUMTT12(9)+GIOT*CSUMTT12(10)-CSUMTT12(11)
     &-GIOT*CSUMTT12(12)
     &-(CSUMTT12(15)+GIOT*CSUMTT12(16)-CSUMTT12(13)-GIOT*CSUMTT12(14))
     &)*ADIV2
c**********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM93(N4,ID,IK)=(CSUMTTMOM12(15,IK)
     &+GIOT*CSUMTTMOM12(16,IK)
     &-CSUMTTMOM12(13,IK)-GIOT*CSUMTTMOM12(14,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE94(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)
     &+(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))
     &+(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))
     &+(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM94(N4,ID,IK)=(CSUMTTMOM12(1,IK)-CSUMTTMOM12(2,IK)
     &+CSUMTTMOM12(3,IK)-CSUMTTMOM12(4,IK)
     &+(CSUMTTMOM12(9,IK)-CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)
     &-CSUMTTMOM12(12,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE95(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)
     &+(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))
     &-(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))
     &-(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM95(N4,ID,IK)=(CSUMTTMOM12(7,IK)-CSUMTTMOM12(8,IK)
     &+CSUMTTMOM12(5,IK)-CSUMTTMOM12(6,IK)+(CSUMTTMOM12(15,IK)
     &-CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)-CSUMTTMOM12(14,IK))
     &)*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE96(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)
     &-(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))
     &+(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))
     &-(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM96(N4,ID,IK)=(CSUMTTMOM12(1,IK)-CSUMTTMOM12(2,IK)
     &+CSUMTTMOM12(3,IK)-CSUMTTMOM12(4,IK)-(CSUMTTMOM12(9,IK)
     &-CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)-CSUMTTMOM12(12,IK))
     &)*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE97(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4)
     &-(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12))
     &-(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6))
     &+(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM97(N4,ID,IK)=(CSUMTTMOM12(7,IK)-CSUMTTMOM12(8,IK)
     &+CSUMTTMOM12(5,IK)-CSUMTTMOM12(6,IK)-(CSUMTTMOM12(15,IK)
     &-CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)-CSUMTTMOM12(14,IK))
     &)*ADIV2
      ENDDO  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     13
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     TT-13 OPERATORS CP=+,Pz=+,J=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE98(N4,ID)=(CSUMTT13(1)+CSUMTT13(2)+CSUMTT13(3)+CSUMTT13(4)+
     &CSUMTT13(5)+CSUMTT13(6)+CSUMTT13(7)+CSUMTT13(8))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM98(N4,ID,IK)=(CSUMTTMOM13(1,IK)+CSUMTTMOM13(2,IK)
     &+CSUMTTMOM13(3,IK)+CSUMTTMOM13(4,IK)+CSUMTTMOM13(5,IK)
     &+CSUMTTMOM13(6,IK)+CSUMTTMOM13(7,IK)
     &+CSUMTTMOM13(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE99(N4,ID)=(CSUMTT13(1)+CSUMTT13(2)+CSUMTT13(3)+CSUMTT13(4)
     &-(CSUMTT13(5)+CSUMTT13(6)+CSUMTT13(7)+CSUMTT13(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM99(N4,ID,IK)=(CSUMTTMOM13(1,IK)+CSUMTTMOM13(2,IK)
     &+CSUMTTMOM13(3,IK)+CSUMTTMOM13(4,IK)-(CSUMTTMOM13(5,IK)
     &+CSUMTTMOM13(6,IK)+CSUMTTMOM13(7,IK)
     &+CSUMTTMOM13(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE100(N4,ID)=(CSUMTT13(1)+GIOT*CSUMTT13(2)-CSUMTT13(3)
     &-GIOT*CSUMTT13(4)+(CSUMTT13(7)+GIOT*CSUMTT13(8)-CSUMTT13(5)
     &-GIOT*CSUMTT13(6)))*ADIV2
c**********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM100(N4,ID,IK)=(CSUMTTMOM13(1,IK)+GIOT*CSUMTTMOM13(2,IK)
     &-CSUMTTMOM13(3,IK)-GIOT*CSUMTTMOM13(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE101(N4,ID)=(CSUMTT13(1)+GIOT*CSUMTT13(2)-CSUMTT13(3)
     &-GIOT*CSUMTT13(4)
     &-(CSUMTT13(7)+GIOT*CSUMTT13(8)-CSUMTT13(5)
     &-GIOT*CSUMTT13(6)))*ADIV2
c**********************************************************************
c     J=1, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM101(N4,ID,IK)=(CSUMTTMOM13(7,IK)+GIOT*CSUMTTMOM13(8,IK)
     &-CSUMTTMOM13(5,IK)-GIOT*CSUMTTMOM13(6,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE102(N4,ID)=(CSUMTT13(1)-CSUMTT13(2)+CSUMTT13(3)-CSUMTT13(4)+
     &CSUMTT13(6)-CSUMTT13(5)+CSUMTT13(8)-CSUMTT13(7))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM102(N4,ID,IK)=(CSUMTTMOM13(1,IK)-CSUMTTMOM13(2,IK)
     &+CSUMTTMOM13(3,IK)-CSUMTTMOM13(4,IK)+CSUMTTMOM13(6,IK)
     &-CSUMTTMOM13(5,IK)+CSUMTTMOM13(8,IK)
     &-CSUMTTMOM13(7,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE103(N4,ID)=(CSUMTT13(1)-CSUMTT13(2)+CSUMTT13(3)-CSUMTT13(4)
     &-(CSUMTT13(6)-CSUMTT13(5)+CSUMTT13(8)-CSUMTT13(7)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM103(N4,ID,IK)=(CSUMTTMOM13(1,IK)-CSUMTTMOM13(2,IK)
     &+CSUMTTMOM13(3,IK)-CSUMTTMOM13(4,IK)-(CSUMTTMOM13(6,IK)
     &-CSUMTTMOM13(5,IK)+CSUMTTMOM13(8,IK)
     &-CSUMTTMOM13(7,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     14 THISSSSSSSS!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE104(N4,ID)=(CSUMTT14(1)+CSUMTT14(2)+CSUMTT14(3)+CSUMTT14(4)+
     &CSUMTT14(5)+CSUMTT14(6)+CSUMTT14(7)+CSUMTT14(8))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=0
c**********************************************************************
      DO IK=1,2
      ALINEMOM104(N4,ID,IK)=(CSUMTTMOM14(1,IK)+CSUMTTMOM14(2,IK)
     &+CSUMTTMOM14(3,IK)+CSUMTTMOM14(4,IK)+CSUMTTMOM14(5,IK)
     &+CSUMTTMOM14(6,IK)+CSUMTTMOM14(7,IK)
     &+CSUMTTMOM14(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE105(N4,ID)=(CSUMTT14(1)+CSUMTT14(2)+CSUMTT14(3)+CSUMTT14(4)-
     &(CSUMTT14(5)+CSUMTT14(6)+CSUMTT14(7)+CSUMTT14(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM105(N4,ID,IK)=(CSUMTTMOM14(1,IK)+CSUMTTMOM14(2,IK)
     &+CSUMTTMOM14(3,IK)+CSUMTTMOM14(4,IK)-(CSUMTTMOM14(5,IK)
     &+CSUMTTMOM14(6,IK)+CSUMTTMOM14(7,IK)
     &+CSUMTTMOM14(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************
      ALINE106(N4,ID)=(CSUMTT14(1)+GIOT*CSUMTT14(2)-CSUMTT14(3)
     &-GIOT*CSUMTT14(4)+(CSUMTT14(6)+GIOT*CSUMTT14(7)-CSUMTT14(8)
     &-GIOT*CSUMTT14(5)))*ADIV2
c**********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
      ALINEMOM106(N4,ID,IK)=(CSUMTTMOM14(1,IK)+GIOT*CSUMTTMOM14(2,IK)
     &-CSUMTTMOM14(3,IK)-GIOT*CSUMTTMOM14(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************
      ALINE107(N4,ID)=(CSUMTT14(1)+GIOT*CSUMTT14(2)-CSUMTT14(3)
     &-GIOT*CSUMTT14(4)-(CSUMTT14(6)+GIOT*CSUMTT14(7)-CSUMTT14(8)
     &-GIOT*CSUMTT14(5)))*ADIV2
c**********************************************************************
c     J=1, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM107(N4,ID,IK)=(CSUMTTMOM14(6,IK)+GIOT*CSUMTTMOM14(7,IK)
     &-CSUMTTMOM14(8,IK)-GIOT*CSUMTTMOM14(5,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE108(N4,ID)=(CSUMTT14(1)-CSUMTT14(2)+CSUMTT14(3)-CSUMTT14(4)+
     &(CSUMTT14(7)-CSUMTT14(6)+CSUMTT14(5)-CSUMTT14(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2,3,4
c**********************************************************************
      DO IK=1,2
         ALINEMOM108(N4,ID,IK)=(CSUMTTMOM14(1,IK)-CSUMTTMOM14(2,IK)
     &+CSUMTTMOM14(3,IK)-CSUMTTMOM14(4,IK)+(CSUMTTMOM14(7,IK)
     &-CSUMTTMOM14(6,IK)+CSUMTTMOM14(5,IK)-CSUMTTMOM14(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE109(N4,ID)=(CSUMTT14(1)-CSUMTT14(2)+CSUMTT14(3)-CSUMTT14(4)-
     &(CSUMTT14(7)-CSUMTT14(6)+CSUMTT14(5)-CSUMTT14(8)))*ADIV2
C**********************************************************************
c     J=2, Pp=-, q=1,2,3,4
c**********************************************************************
      DO IK=1,2 
      ALINEMOM109(N4,ID,IK)=(CSUMTTMOM14(1,IK)-CSUMTTMOM14(2,IK)
     &+CSUMTTMOM14(3,IK)-CSUMTTMOM14(4,IK)-(CSUMTTMOM14(7,IK)
     &-CSUMTTMOM14(6,IK)+CSUMTTMOM14(5,IK)-CSUMTTMOM14(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATOR 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE110(N4,ID)=(CSUMPLQ(1)+CSUMPLQ(2)+CSUMPLQ(3)+CSUMPLQ(4)
     &+CSUMPLQ(5)+CSUMPLQ(6)+CSUMPLQ(7)+CSUMPLQ(8))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM110(N4,ID,IK)=(CSUMPLQMOM(1,IK)+CSUMPLQMOM(2,IK)
     &+CSUMPLQMOM(3,IK)+CSUMPLQMOM(4,IK)+CSUMPLQMOM(5,IK)
     &+CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)+CSUMPLQMOM(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************      
      ALINE111(N4,ID)=(CSUMPLQ(1)+CSUMPLQ(2)+CSUMPLQ(3)+CSUMPLQ(4)
     &-CSUMPLQ(5)-CSUMPLQ(6)-CSUMPLQ(7)-CSUMPLQ(8))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM111(N4,ID,IK)=(CSUMPLQMOM(1,IK)+CSUMPLQMOM(2,IK)
     &+CSUMPLQMOM(3,IK)+CSUMPLQMOM(4,IK)-CSUMPLQMOM(5,IK)
     &-CSUMPLQMOM(6,IK)-CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************      
      ALINE112(N4,ID)=(CSUMPLQ(1)+GIOT*CSUMPLQ(2)-CSUMPLQ(3)
     &-GIOT*CSUMPLQ(4)
     &+(CSUMPLQ(6)+GIOT*CSUMPLQ(5)-CSUMPLQ(8)-GIOT*CSUMPLQ(7)))*ADIV2
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************      
      DO IK=1, 2
         ALINEMOM112(N4,ID,IK)=(CSUMPLQMOM(1,IK)+GIOT*CSUMPLQMOM(2,IK)
     &-CSUMPLQMOM(3,IK)-GIOT*CSUMPLQMOM(4,IK))*ADIV1
      ENDDO   
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************      
      ALINE113(N4,ID)=(CSUMPLQ(1)+GIOT*CSUMPLQ(2)-CSUMPLQ(3)
     &-GIOT*CSUMPLQ(4)
     &-(CSUMPLQ(6)+GIOT*CSUMPLQ(5)-CSUMPLQ(8)-GIOT*CSUMPLQ(7)))*ADIV2
c**********************************************************************
c     J=1, q=0
c**********************************************************************      
      DO IK=1, 2
         ALINEMOM113(N4,ID,IK)=(CSUMPLQMOM(6,IK)+GIOT*CSUMPLQMOM(5,IK)
     &-CSUMPLQMOM(8,IK)-GIOT*CSUMPLQMOM(7,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************      
      ALINE114(N4,ID)=(CSUMPLQ(1)-CSUMPLQ(2)+CSUMPLQ(3)-CSUMPLQ(4)
     &+CSUMPLQ(5)-CSUMPLQ(6)+CSUMPLQ(7)-CSUMPLQ(8))*ADIV2
c**********************************************************************
c     J=2, Pr=+, q=0
c**********************************************************************      
      DO IK=1, 2
         ALINEMOM114(N4,ID,IK)=(CSUMPLQMOM(1,IK)-CSUMPLQMOM(2,IK)
     &+CSUMPLQMOM(3,IK)-CSUMPLQMOM(4,IK)+CSUMPLQMOM(5,IK)
     &-CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************            
      ALINE115(N4,ID)=(CSUMPLQ(1)-CSUMPLQ(2)+CSUMPLQ(3)-CSUMPLQ(4)
     &-(CSUMPLQ(5)-CSUMPLQ(6)+CSUMPLQ(7)-CSUMPLQ(8)))*ADIV2
c**********************************************************************
c     J=2, Pr=-, q=0
c**********************************************************************            
      DO IK=1, 2
         ALINEMOM115(N4,ID,IK)=(CSUMPLQMOM(1,IK)-CSUMPLQMOM(2,IK)
     &+CSUMPLQMOM(3,IK)-CSUMPLQMOM(4,IK)-(CSUMPLQMOM(5,IK)
     &-CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATOR 2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************                  
      ALINE116(N4,ID)=(CSUMPLQ2(1)+CSUMPLQ2(2)+CSUMPLQ2(3)+CSUMPLQ2(4)
     &+CSUMPLQ2(5)+CSUMPLQ2(6)+CSUMPLQ2(7)+CSUMPLQ2(8))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM116(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+CSUMPLQMOM2(2,IK)
     &+CSUMPLQMOM2(3,IK)+CSUMPLQMOM2(4,IK)+CSUMPLQMOM2(5,IK)
     &+CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)+CSUMPLQMOM2(8,IK))*ADIV2
      ENDDO   
c**********************************************************************
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************                  
      ALINE117(N4,ID)=(CSUMPLQ2(1)+CSUMPLQ2(2)+CSUMPLQ2(3)+CSUMPLQ2(4)
     &-(CSUMPLQ2(5)+CSUMPLQ2(6)+CSUMPLQ2(7)+CSUMPLQ2(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM117(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+CSUMPLQMOM2(2,IK)
     &+CSUMPLQMOM2(3,IK)+CSUMPLQMOM2(4,IK)-(CSUMPLQMOM2(5,IK)
     &+CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)+CSUMPLQMOM2(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************                        
      ALINE118(N4,ID)=(CSUMPLQ2(1)+GIOT*CSUMPLQ2(2)-CSUMPLQ2(3)
     &-GIOT*CSUMPLQ2(4)
     &+CSUMPLQ2(6)+GIOT*CSUMPLQ2(5)-CSUMPLQ2(8)-GIOT*CSUMPLQ2(7))*ADIV2
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM118(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+GIOT*CSUMPLQMOM2(2,IK)
     &-CSUMPLQMOM2(3,IK)-GIOT*CSUMPLQMOM2(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=-, q=0
c**********************************************************************                            
      ALINE119(N4,ID)=(CSUMPLQ2(1)+GIOT*CSUMPLQ2(2)-CSUMPLQ2(3)
     &-GIOT*CSUMPLQ2(4)
     &-(CSUMPLQ2(6)+GIOT*CSUMPLQ2(5)-CSUMPLQ2(8)-GIOT*CSUMPLQ2(7)))
     &*ADIV2
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM119(N4,ID,IK)=(CSUMPLQMOM2(6,IK)+GIOT*CSUMPLQMOM2(5,IK)
     &-CSUMPLQMOM2(8,IK)-GIOT*CSUMPLQMOM2(7,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************                  
      ALINE120(N4,ID)=(CSUMPLQ2(1)-CSUMPLQ2(2)+CSUMPLQ2(3)-CSUMPLQ2(4)
     &+CSUMPLQ2(5)-CSUMPLQ2(6)+CSUMPLQ2(7)-CSUMPLQ2(8))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM120(N4,ID,IK)=(CSUMPLQMOM2(1,IK)-CSUMPLQMOM2(2,IK)
     &+CSUMPLQMOM2(3,IK)-CSUMPLQMOM2(4,IK)+CSUMPLQMOM2(5,IK)
     &-CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)-CSUMPLQMOM2(8,IK))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************                        
      ALINE121(N4,ID)=(CSUMPLQ2(1)-CSUMPLQ2(2)+CSUMPLQ2(3)-CSUMPLQ2(4)
     &-(CSUMPLQ2(5)-CSUMPLQ2(6)+CSUMPLQ2(7)-CSUMPLQ2(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM121(N4,ID,IK)=(CSUMPLQMOM2(1,IK)-CSUMPLQMOM2(2,IK)
     &+CSUMPLQMOM2(3,IK)-CSUMPLQMOM2(4,IK)-(CSUMPLQMOM2(5,IK)
     &-CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)-CSUMPLQMOM2(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLAQUETTE OPERATORS 3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************                        
      ALINE122(N4,ID)=(CSUMPLQ3(1)+CSUMPLQ3(2)+CSUMPLQ3(3)
     &+CSUMPLQ3(4)+(CSUMPLQ3(5)+CSUMPLQ3(6)+CSUMPLQ3(7)
     &+CSUMPLQ3(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=+, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM122(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+CSUMPLQMOM3(2,IK)
     &+CSUMPLQMOM3(3,IK)+CSUMPLQMOM3(4,IK)
     &+(CSUMPLQMOM3(5,IK)+CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)
     &+CSUMPLQMOM3(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************                        
      ALINE123(N4,ID)=(CSUMPLQ3(1)+CSUMPLQ3(2)+CSUMPLQ3(3)
     &+CSUMPLQ3(4)-(CSUMPLQ3(5)+CSUMPLQ3(6)+CSUMPLQ3(7)
     &+CSUMPLQ3(8)))*ADIV2
c**********************************************************************
c     J=0, Pp=-, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM123(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+CSUMPLQMOM3(2,IK)
     &+CSUMPLQMOM3(3,IK)+CSUMPLQMOM3(4,IK)
     &-(CSUMPLQMOM3(5,IK)+CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)
     &+CSUMPLQMOM3(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************                            
      ALINE124(N4,ID)=(CSUMPLQ3(1)+GIOT*CSUMPLQ3(2)-CSUMPLQ3(3)
     &-GIOT*CSUMPLQ3(4))*ADIV1
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************                            
      DO IK=1, 2
         ALINEMOM124(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+GIOT*CSUMPLQMOM3(2,IK)
     &-CSUMPLQMOM3(3,IK)-GIOT*CSUMPLQMOM3(4,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=1, Pr=+, q=0
c**********************************************************************                            
      ALINE125(N4,ID)=(CSUMPLQ3(7)+GIOT*CSUMPLQ3(6)-CSUMPLQ3(5)
     &-GIOT*CSUMPLQ3(8))*ADIV1
c**********************************************************************
c     J=1, q=1,2
c**********************************************************************                            
      DO IK=1, 2
         ALINEMOM125(N4,ID,IK)=(CSUMPLQMOM3(7,IK)+GIOT*CSUMPLQMOM3(6,IK)
     &-CSUMPLQMOM3(5,IK)-GIOT*CSUMPLQMOM3(8,IK))*ADIV1
      ENDDO
c**********************************************************************
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************                            
      ALINE126(N4,ID)=(CSUMPLQ3(1)-CSUMPLQ3(2)+CSUMPLQ3(3)
     &-CSUMPLQ3(4)+(CSUMPLQ3(5)-CSUMPLQ3(6)+CSUMPLQ3(7)
     &-CSUMPLQ3(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=+, q=1,2
c**********************************************************************                            
      DO IK=1, 2
         ALINEMOM126(N4,ID,IK)=(CSUMPLQMOM3(1,IK)-CSUMPLQMOM3(2,IK)
     &+CSUMPLQMOM3(3,IK)-CSUMPLQMOM3(4,IK)
     &+(CSUMPLQMOM3(5,IK)-CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)
     &-CSUMPLQMOM3(8,IK)))*ADIV2
      ENDDO
c**********************************************************************
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************                            
      ALINE127(N4,ID)=(CSUMPLQ3(1)-CSUMPLQ3(2)+CSUMPLQ3(3)
     &-CSUMPLQ3(4)-(CSUMPLQ3(5)-CSUMPLQ3(6)+CSUMPLQ3(7)
     &-CSUMPLQ3(8)))*ADIV2
c**********************************************************************
c     J=2, Pp=-, q=1,2
c**********************************************************************                            
      DO IK=1, 2
         ALINEMOM127(N4,ID,IK)=(CSUMPLQMOM3(1,IK)-CSUMPLQMOM3(2,IK)
     &+CSUMPLQMOM3(3,IK)-CSUMPLQMOM3(4,IK)
     &-(CSUMPLQMOM3(5,IK)-CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK)
     &-CSUMPLQMOM3(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c**********************************************************************      
C     J=0, Pp=+, Pr=+, q=0
c**********************************************************************            
      ALINE128(N4,ID)=(CSUMPLQ4(1)+CSUMPLQ4(2)+CSUMPLQ4(3)
     &+CSUMPLQ4(4)+(CSUMPLQ4(5)+CSUMPLQ4(6)+CSUMPLQ4(7)
     &+CSUMPLQ4(8)))*ADIV2
c**********************************************************************      
C     J=0, Pp=+, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM128(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+CSUMPLQMOM4(2,IK)
     &+CSUMPLQMOM4(3,IK)+CSUMPLQMOM4(4,IK)
     &+(CSUMPLQMOM4(5,IK)+CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)
     &+CSUMPLQMOM4(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
C     J=0, Pp=-, Pr=+, q=0
c**********************************************************************                  
      ALINE129(N4,ID)=(CSUMPLQ4(1)+CSUMPLQ4(2)+CSUMPLQ4(3)
     &+CSUMPLQ4(4)-(CSUMPLQ4(5)+CSUMPLQ4(6)+CSUMPLQ4(7)
     &+CSUMPLQ4(8)))*ADIV2
c**********************************************************************      
C     J=0, Pp=-, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM129(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+CSUMPLQMOM4(2,IK)
     &+CSUMPLQMOM4(3,IK)+CSUMPLQMOM4(4,IK)
     &-(CSUMPLQMOM4(5,IK)+CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)
     &+CSUMPLQMOM4(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
C     J=1, Pr=+, q=0
c**********************************************************************                            
      ALINE130(N4,ID)=(CSUMPLQ4(1)+GIOT*CSUMPLQ4(2)-CSUMPLQ4(3)
     &-GIOT*CSUMPLQ4(4))*ADIV1
c**********************************************************************      
C     J=1, q=1,2
c**********************************************************************                            
      DO IK=1, 2
         ALINEMOM130(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+GIOT*CSUMPLQMOM4(2,IK)
     &-CSUMPLQMOM4(3,IK)-GIOT*CSUMPLQMOM4(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
C     J=1, Pr=+, q=0
c**********************************************************************                            
      ALINE131(N4,ID)=(CSUMPLQ4(7)+GIOT*CSUMPLQ4(6)-CSUMPLQ4(5)
     &-GIOT*CSUMPLQ4(8))*ADIV1
c**********************************************************************      
C     J=1, q=1,2
c**********************************************************************                            
      DO IK=1, 2
         ALINEMOM131(N4,ID,IK)=(CSUMPLQMOM4(7,IK)+GIOT*CSUMPLQMOM4(6,IK)
     &-CSUMPLQMOM4(5,IK)-GIOT*CSUMPLQMOM4(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************            
      ALINE132(N4,ID)=(CSUMPLQ4(1)-CSUMPLQ4(2)+CSUMPLQ4(3)
     &-CSUMPLQ4(4)+(CSUMPLQ4(5)-CSUMPLQ4(6)+CSUMPLQ4(7)
     &-CSUMPLQ4(8)))*ADIV2
c**********************************************************************      
c     J=2, Pp=+, q=1,2
c**********************************************************************            
      DO IK=1, 2
         ALINEMOM132(N4,ID,IK)=(CSUMPLQMOM4(1,IK)-CSUMPLQMOM4(2,IK)
     &+CSUMPLQMOM4(3,IK)-CSUMPLQMOM4(4,IK)
     &+(CSUMPLQMOM4(5,IK)-CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)
     &-CSUMPLQMOM4(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************                  
      ALINE133(N4,ID)=(CSUMPLQ4(1)-CSUMPLQ4(2)+CSUMPLQ4(3)
     &-CSUMPLQ4(4)-(CSUMPLQ4(5)-CSUMPLQ4(6)+CSUMPLQ4(7)
     &-CSUMPLQ4(8)))*ADIV2
c**********************************************************************      
c     J=2, Pp=-, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM133(N4,ID,IK)=(CSUMPLQMOM4(1,IK)-CSUMPLQMOM4(2,IK)
     &+CSUMPLQMOM4(3,IK)-CSUMPLQMOM4(4,IK)
     &-(CSUMPLQMOM4(5,IK)-CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK)
     &-CSUMPLQMOM4(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLAQUETTE OPERATOR 5  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************                        
      ALINE134(N4,ID)=(CSUMPLQ5(1)+CSUMPLQ5(2)+CSUMPLQ5(3)
     &+CSUMPLQ5(4)+(CSUMPLQ5(5)+CSUMPLQ5(6)+CSUMPLQ5(7)
     &+CSUMPLQ5(8)))*ADIV2
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************                            
      DO IK=1, 2
         ALINEMOM134(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+CSUMPLQMOM5(2,IK)
     &+CSUMPLQMOM5(3,IK)+CSUMPLQMOM5(4,IK)
     &+(CSUMPLQMOM5(5,IK)+CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)
     &+CSUMPLQMOM5(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************                            
      ALINE135(N4,ID)=(CSUMPLQ5(1)+CSUMPLQ5(2)+CSUMPLQ5(3)
     &+CSUMPLQ5(4)-(CSUMPLQ5(5)+CSUMPLQ5(6)+CSUMPLQ5(7)
     &+CSUMPLQ5(8)))*ADIV2
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************      
      DO IK=1, 2
         ALINEMOM135(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+CSUMPLQMOM5(2,IK)
     &+CSUMPLQMOM5(3,IK)+CSUMPLQMOM5(4,IK)
     &-(CSUMPLQMOM5(5,IK)+CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)
     &+CSUMPLQMOM5(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************            
      ALINE136(N4,ID)=(CSUMPLQ5(1)+GIOT*CSUMPLQ5(2)-CSUMPLQ5(3)
     &-GIOT*CSUMPLQ5(4))*ADIV1
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM136(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+GIOT*CSUMPLQMOM5(2,IK)
     &-CSUMPLQMOM5(3,IK)-GIOT*CSUMPLQMOM5(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE137(N4,ID)=(CSUMPLQ5(7)+GIOT*CSUMPLQ5(6)-CSUMPLQ5(5)
     &-GIOT*CSUMPLQ5(8))*ADIV1
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM137(N4,ID,IK)=(CSUMPLQMOM5(7,IK)+GIOT*CSUMPLQMOM5(6,IK)
     &-CSUMPLQMOM5(5,IK)-GIOT*CSUMPLQMOM5(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************                        
      ALINE138(N4,ID)=(CSUMPLQ5(1)-CSUMPLQ5(2)+CSUMPLQ5(3)
     &-CSUMPLQ5(4)+(CSUMPLQ5(5)-CSUMPLQ5(6)+CSUMPLQ5(7)
     &-CSUMPLQ5(8)))*ADIV2
c**********************************************************************      
c     J=2, Pp=+, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM138(N4,ID,IK)=(CSUMPLQMOM5(1,IK)-CSUMPLQMOM5(2,IK)
     &+CSUMPLQMOM5(3,IK)-CSUMPLQMOM5(4,IK)
     &+(CSUMPLQMOM5(5,IK)-CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)
     &-CSUMPLQMOM5(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************                            
      ALINE139(N4,ID)=(CSUMPLQ5(1)-CSUMPLQ5(2)+CSUMPLQ5(3)
     &-CSUMPLQ5(4)-(CSUMPLQ5(5)-CSUMPLQ5(6)+CSUMPLQ5(7)
     &-CSUMPLQ5(8)))*ADIV2
c**********************************************************************      
c     J=2, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM139(N4,ID,IK)=(CSUMPLQMOM5(1,IK)-CSUMPLQMOM5(2,IK)
     &+CSUMPLQMOM5(3,IK)-CSUMPLQMOM5(4,IK)
     &-(CSUMPLQMOM5(5,IK)-CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK)
     &-CSUMPLQMOM5(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 6
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE140(N4,ID)=(CSUMPLQ6(1)+CSUMPLQ6(2)+CSUMPLQ6(3)
     &+CSUMPLQ6(4)+(CSUMPLQ6(5)+CSUMPLQ6(6)+CSUMPLQ6(7)
     &+CSUMPLQ6(8)))*ADIV2
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM140(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+CSUMPLQMOM6(2,IK)
     &+CSUMPLQMOM6(3,IK)+CSUMPLQMOM6(4,IK)
     &+(CSUMPLQMOM6(5,IK)+CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)
     &+CSUMPLQMOM6(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************      
      ALINE141(N4,ID)=(CSUMPLQ6(1)+CSUMPLQ6(2)+CSUMPLQ6(3)
     &+CSUMPLQ6(4)-(CSUMPLQ6(5)+CSUMPLQ6(6)+CSUMPLQ6(7)
     &     +CSUMPLQ6(8)))*ADIV2
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************            
      DO IK=1, 2
         ALINEMOM141(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+CSUMPLQMOM6(2,IK)
     &+CSUMPLQMOM6(3,IK)+CSUMPLQMOM6(4,IK)
     &-(CSUMPLQMOM6(5,IK)+CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)
     &+CSUMPLQMOM6(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE142(N4,ID)=(CSUMPLQ6(1)+GIOT*CSUMPLQ6(2)-CSUMPLQ6(3)
     &-GIOT*CSUMPLQ6(4)+(CSUMPLQ6(8)+GIOT*CSUMPLQ6(7)-CSUMPLQ6(6)
     &-GIOT*CSUMPLQ6(5)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM142(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+GIOT*CSUMPLQMOM6(2,IK)
     &-CSUMPLQMOM6(3,IK)-GIOT*CSUMPLQMOM6(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE143(N4,ID)=(CSUMPLQ6(1)+GIOT*CSUMPLQ6(2)-CSUMPLQ6(3)
     &-GIOT*CSUMPLQ6(4)-(CSUMPLQ6(8)+GIOT*CSUMPLQ6(7)-CSUMPLQ6(6)
     &-GIOT*CSUMPLQ6(5)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM143(N4,ID,IK)=(CSUMPLQMOM6(8,IK)+GIOT*CSUMPLQMOM6(7,IK)
     &-CSUMPLQMOM6(6,IK)-GIOT*CSUMPLQMOM6(5,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************                        
      ALINE144(N4,ID)=(CSUMPLQ6(1)-CSUMPLQ6(2)+CSUMPLQ6(3)
     &-CSUMPLQ6(4)+(CSUMPLQ6(5)-CSUMPLQ6(6)+CSUMPLQ6(7)
     &-CSUMPLQ6(8)))*ADIV2
c**********************************************************************      
c     J=2, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM144(N4,ID,IK)=(CSUMPLQMOM6(1,IK)-CSUMPLQMOM6(2,IK)
     &+CSUMPLQMOM6(3,IK)-CSUMPLQMOM6(4,IK)
     &+(CSUMPLQMOM6(5,IK)-CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)
     &-CSUMPLQMOM6(8,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE145(N4,ID)=(CSUMPLQ6(1)-CSUMPLQ6(2)+CSUMPLQ6(3)
     &-CSUMPLQ6(4)-(CSUMPLQ6(5)-CSUMPLQ6(6)+CSUMPLQ6(7)
     &-CSUMPLQ6(8)))*ADIV2
c**********************************************************************      
c     J=2, Pp=-, q=1,2
c**********************************************************************      
      DO IK=1, 2
         ALINEMOM145(N4,ID,IK)=(CSUMPLQMOM6(1,IK)-CSUMPLQMOM6(2,IK)
     &+CSUMPLQMOM6(3,IK)-CSUMPLQMOM6(4,IK)
     &-(CSUMPLQMOM6(5,IK)-CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK)
     &-CSUMPLQMOM6(8,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 7
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE146(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)
     &+CSUMPLQ7(4)
     &+(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))
     &+(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))
     &+(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM146(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+CSUMPLQMOM7(2,IK)
     &+CSUMPLQMOM7(3,IK)+CSUMPLQMOM7(4,IK)
     &+(CSUMPLQMOM7(9,IK)+CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)
     &+CSUMPLQMOM7(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE147(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)
     &+CSUMPLQ7(4)
     &-(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))
     &+(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))
     &-(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM147(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+CSUMPLQMOM7(6,IK)
     &+CSUMPLQMOM7(7,IK)+CSUMPLQMOM7(8,IK)
     &+(CSUMPLQMOM7(13,IK)+CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)
     &+CSUMPLQMOM7(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE148(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)
     &+CSUMPLQ7(4)
     &+(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))
     &-(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))
     &-(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM148(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+CSUMPLQMOM7(2,IK)
     &+CSUMPLQMOM7(3,IK)+CSUMPLQMOM7(4,IK)
     &-(CSUMPLQMOM7(9,IK)+CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)
     &+CSUMPLQMOM7(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE149(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3)
     &+CSUMPLQ7(4)
     &-(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8))
     &-(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12))
     &+(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM149(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+CSUMPLQMOM7(6,IK)
     &+CSUMPLQMOM7(7,IK)+CSUMPLQMOM7(8,IK)
     &-(CSUMPLQMOM7(13,IK)+CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)
     &+CSUMPLQMOM7(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE150(N4,ID)=(CSUMPLQ7(1)+GIOT*CSUMPLQ7(2)-CSUMPLQ7(3)
     &-GIOT*CSUMPLQ7(4)+(CSUMPLQ7(5)+GIOT*CSUMPLQ7(6)-CSUMPLQ7(7)
     &-GIOT*CSUMPLQ7(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM150(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+GIOT*CSUMPLQMOM7(2,IK)
     &-CSUMPLQMOM7(3,IK)-GIOT*CSUMPLQMOM7(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE151(N4,ID)=(CSUMPLQ7(1)+GIOT*CSUMPLQ7(2)-CSUMPLQ7(3)
     &-GIOT*CSUMPLQ7(4)-(CSUMPLQ7(5)+GIOT*CSUMPLQ7(6)-CSUMPLQ7(7)
     &-GIOT*CSUMPLQ7(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM151(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+GIOT*CSUMPLQMOM7(6,IK)
     &-CSUMPLQMOM7(7,IK)-GIOT*CSUMPLQMOM7(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE152(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)
     &-CSUMPLQ7(4)
     &+(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))
     &+(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))
     &+(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q
c**********************************************************************
      DO IK=1, 2
         ALINEMOM152(N4,ID,IK)=(CSUMPLQMOM7(1,IK)-CSUMPLQMOM7(2,IK)
     &+CSUMPLQMOM7(3,IK)-CSUMPLQMOM7(4,IK)
     &+(CSUMPLQMOM7(9,IK)-CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)
     &-CSUMPLQMOM7(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE153(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)
     &-CSUMPLQ7(4)
     &-(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))
     &+(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))
     &-(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM153(N4,ID,IK)=(CSUMPLQMOM7(5,IK)-CSUMPLQMOM7(6,IK)
     &+CSUMPLQMOM7(7,IK)-CSUMPLQMOM7(8,IK)
     &+(CSUMPLQMOM7(13,IK)-CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)
     &-CSUMPLQMOM7(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE154(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)
     &-CSUMPLQ7(4)
     &+(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))
     &-(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))
     &-(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM154(N4,ID,IK)=(CSUMPLQMOM7(1,IK)-CSUMPLQMOM7(2,IK)
     &+CSUMPLQMOM7(3,IK)-CSUMPLQMOM7(4,IK)
     &-(CSUMPLQMOM7(9,IK)-CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK)
     &-CSUMPLQMOM7(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE155(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3)
     &-CSUMPLQ7(4)
     &-(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8))
     &-(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12))
     &+(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM155(N4,ID,IK)=(CSUMPLQMOM7(5,IK)-CSUMPLQMOM7(6,IK)
     &+CSUMPLQMOM7(7,IK)-CSUMPLQMOM7(8,IK)
     &-(CSUMPLQMOM7(13,IK)-CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK)
     &-CSUMPLQMOM7(16,IK)))*ADIV2
      ENDDO      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 8
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE156(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)
     &+CSUMPLQ8(4)
     &+(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))
     &+(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))
     &+(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM156(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+CSUMPLQMOM8(2,IK)
     &+CSUMPLQMOM8(3,IK)+CSUMPLQMOM8(4,IK)
     &+(CSUMPLQMOM8(9,IK)+CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)
     &+CSUMPLQMOM8(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE157(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)
     &+CSUMPLQ8(4)
     &-(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))
     &+(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))
     &-(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM157(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+CSUMPLQMOM8(6,IK)
     &+CSUMPLQMOM8(7,IK)+CSUMPLQMOM8(8,IK)
     &+(CSUMPLQMOM8(13,IK)+CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)
     &+CSUMPLQMOM8(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE158(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)
     &+CSUMPLQ8(4)
     &+(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))
     &-(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))
     &-(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM158(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+CSUMPLQMOM8(2,IK)
     &+CSUMPLQMOM8(3,IK)+CSUMPLQMOM8(4,IK)
     &-(CSUMPLQMOM8(9,IK)+CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)
     &+CSUMPLQMOM8(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE159(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3)
     &+CSUMPLQ8(4)
     &-(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8))
     &-(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12))
     &+(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM159(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+CSUMPLQMOM8(6,IK)
     &+CSUMPLQMOM8(7,IK)+CSUMPLQMOM8(8,IK)
     &-(CSUMPLQMOM8(13,IK)+CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)
     &+CSUMPLQMOM8(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE160(N4,ID)=(CSUMPLQ8(1)+GIOT*CSUMPLQ8(2)-CSUMPLQ8(3)
     &-GIOT*CSUMPLQ8(4)+(CSUMPLQ8(5)+GIOT*CSUMPLQ8(6)-CSUMPLQ8(7)
     &-GIOT*CSUMPLQ8(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM160(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+GIOT*CSUMPLQMOM8(2,IK)
     &-CSUMPLQMOM8(3,IK)-GIOT*CSUMPLQMOM8(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE161(N4,ID)=(CSUMPLQ8(1)+GIOT*CSUMPLQ8(2)-CSUMPLQ8(3)
     &-GIOT*CSUMPLQ8(4)-(CSUMPLQ8(5)+GIOT*CSUMPLQ8(6)-CSUMPLQ8(7)
     &-GIOT*CSUMPLQ8(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM161(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+GIOT*CSUMPLQMOM8(6,IK)
     &-CSUMPLQMOM8(7,IK)-GIOT*CSUMPLQMOM8(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE162(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)
     &-CSUMPLQ8(4)
     &+(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))
     &+(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))
     &+(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM162(N4,ID,IK)=(CSUMPLQMOM8(1,IK)-CSUMPLQMOM8(2,IK)
     &+CSUMPLQMOM8(3,IK)-CSUMPLQMOM8(4,IK)
     &+(CSUMPLQMOM8(9,IK)-CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)
     &-CSUMPLQMOM8(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE163(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)
     &-CSUMPLQ8(4)
     &-(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))
     &+(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))
     &-(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM163(N4,ID,IK)=(CSUMPLQMOM8(5,IK)-CSUMPLQMOM8(6,IK)
     &+CSUMPLQMOM8(7,IK)-CSUMPLQMOM8(8,IK)
     &+(CSUMPLQMOM8(13,IK)-CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)
     &-CSUMPLQMOM8(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE164(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)
     &-CSUMPLQ8(4)
     &+(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))
     &-(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))
     &-(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM164(N4,ID,IK)=(CSUMPLQMOM8(1,IK)-CSUMPLQMOM8(2,IK)
     &+CSUMPLQMOM8(3,IK)-CSUMPLQMOM8(4,IK)
     &-(CSUMPLQMOM8(9,IK)-CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK)
     &-CSUMPLQMOM8(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE165(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3)
     &-CSUMPLQ8(4)
     &-(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8))
     &-(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12))
     &+(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM165(N4,ID,IK)=(CSUMPLQMOM8(5,IK)-CSUMPLQMOM8(6,IK)
     &+CSUMPLQMOM8(7,IK)-CSUMPLQMOM8(8,IK)
     &-(CSUMPLQMOM8(13,IK)-CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK)
     &-CSUMPLQMOM8(16,IK)))*ADIV2
      ENDDO      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 9
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE166(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)
     &+CSUMPLQ9(4)
     &+(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))
     &+(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))
     &+(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM166(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+CSUMPLQMOM9(2,IK)
     &+CSUMPLQMOM9(3,IK)+CSUMPLQMOM9(4,IK)
     &+(CSUMPLQMOM9(9,IK)+CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)
     &+CSUMPLQMOM9(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE167(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)
     &+CSUMPLQ9(4)
     &-(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))
     &+(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))
     &-(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM167(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+CSUMPLQMOM9(6,IK)
     &+CSUMPLQMOM9(7,IK)+CSUMPLQMOM9(8,IK)
     &+(CSUMPLQMOM9(13,IK)+CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)
     &+CSUMPLQMOM9(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE168(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)
     &+CSUMPLQ9(4)
     &+(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))
     &-(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))
     &-(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM168(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+CSUMPLQMOM9(2,IK)
     &+CSUMPLQMOM9(3,IK)+CSUMPLQMOM9(4,IK)
     &-(CSUMPLQMOM9(9,IK)+CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)
     &+CSUMPLQMOM9(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE169(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3)
     &+CSUMPLQ9(4)
     &-(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8))
     &-(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12))
     &+(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM169(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+CSUMPLQMOM9(6,IK)
     &+CSUMPLQMOM9(7,IK)+CSUMPLQMOM9(8,IK)
     &-(CSUMPLQMOM9(13,IK)+CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)
     &+CSUMPLQMOM9(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE170(N4,ID)=(CSUMPLQ9(1)+GIOT*CSUMPLQ9(2)-CSUMPLQ9(3)
     &-GIOT*CSUMPLQ9(4)+(CSUMPLQ9(5)+GIOT*CSUMPLQ9(6)-CSUMPLQ9(7)
     &-GIOT*CSUMPLQ9(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM170(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+GIOT*CSUMPLQMOM9(2,IK)
     &-CSUMPLQMOM9(3,IK)-GIOT*CSUMPLQMOM9(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE171(N4,ID)=(CSUMPLQ9(1)+GIOT*CSUMPLQ9(2)-CSUMPLQ9(3)
     &-GIOT*CSUMPLQ9(4)-(CSUMPLQ9(5)+GIOT*CSUMPLQ9(6)-CSUMPLQ9(7)
     &-GIOT*CSUMPLQ9(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM171(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+GIOT*CSUMPLQMOM9(6,IK)
     &-CSUMPLQMOM9(7,IK)-GIOT*CSUMPLQMOM9(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE172(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)
     &-CSUMPLQ9(4)
     &+(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))
     &+(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))
     &+(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM172(N4,ID,IK)=(CSUMPLQMOM9(1,IK)-CSUMPLQMOM9(2,IK)
     &+CSUMPLQMOM9(3,IK)-CSUMPLQMOM9(4,IK)
     &+(CSUMPLQMOM9(9,IK)-CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)
     &-CSUMPLQMOM9(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE173(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)
     &-CSUMPLQ9(4)
     &-(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))
     &+(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))
     &-(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM173(N4,ID,IK)=(CSUMPLQMOM9(5,IK)-CSUMPLQMOM9(6,IK)
     &+CSUMPLQMOM9(7,IK)-CSUMPLQMOM9(8,IK)
     &+(CSUMPLQMOM9(13,IK)-CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)
     &-CSUMPLQMOM9(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE174(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)
     &-CSUMPLQ9(4)
     &+(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))
     &-(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))
     &-(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM174(N4,ID,IK)=(CSUMPLQMOM9(1,IK)-CSUMPLQMOM9(2,IK)
     &+CSUMPLQMOM9(3,IK)-CSUMPLQMOM9(4,IK)
     &-(CSUMPLQMOM9(9,IK)-CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK)
     &-CSUMPLQMOM9(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE175(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3)
     &-CSUMPLQ9(4)
     &-(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8))
     &-(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12))
     &+(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM175(N4,ID,IK)=(CSUMPLQMOM9(5,IK)-CSUMPLQMOM9(6,IK)
     &+CSUMPLQMOM9(7,IK)-CSUMPLQMOM9(8,IK)
     &-(CSUMPLQMOM9(13,IK)-CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK)
     &-CSUMPLQMOM9(16,IK)))*ADIV2
      ENDDO      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 10
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE176(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)
     &+CSUMPLQ10(4)
     &+(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))
     &+(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))
     &+(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM176(N4,ID,IK)=(CSUMPLQMOM10(1,IK)+CSUMPLQMOM10(2,IK)
     &+CSUMPLQMOM10(3,IK)+CSUMPLQMOM10(4,IK)
     &+(CSUMPLQMOM10(9,IK)+CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)
     &+CSUMPLQMOM10(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE177(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)
     &+CSUMPLQ10(4)
     &-(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))
     &+(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))
     &-(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM177(N4,ID,IK)=(CSUMPLQMOM10(5,IK)+CSUMPLQMOM10(6,IK)
     &+CSUMPLQMOM10(7,IK)+CSUMPLQMOM10(8,IK)
     &+(CSUMPLQMOM10(13,IK)+CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)
     &+CSUMPLQMOM10(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE178(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)
     &+CSUMPLQ10(4)
     &+(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))
     &-(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))
     &-(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM178(N4,ID,IK)=(CSUMPLQMOM10(1,IK)+CSUMPLQMOM10(2,IK)
     &+CSUMPLQMOM10(3,IK)+CSUMPLQMOM10(4,IK)
     &-(CSUMPLQMOM10(9,IK)+CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)
     &+CSUMPLQMOM10(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE179(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3)
     &+CSUMPLQ10(4)
     &-(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8))
     &-(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12))
     &+(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM179(N4,ID,IK)=(CSUMPLQMOM10(5,IK)+CSUMPLQMOM10(6,IK)
     &+CSUMPLQMOM10(7,IK)+CSUMPLQMOM10(8,IK)
     &-(CSUMPLQMOM10(13,IK)+CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)
     &+CSUMPLQMOM10(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE180(N4,ID)=(CSUMPLQ10(1)+GIOT*CSUMPLQ10(2)-CSUMPLQ10(3)
     &-GIOT*CSUMPLQ10(4)+(CSUMPLQ10(5)+GIOT*CSUMPLQ10(6)-CSUMPLQ10(7)
     &-GIOT*CSUMPLQ10(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM180(N4,ID,IK)=(CSUMPLQMOM10(1,IK)
     &+GIOT*CSUMPLQMOM10(2,IK)
     &-CSUMPLQMOM10(3,IK)-GIOT*CSUMPLQMOM10(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE181(N4,ID)=(CSUMPLQ10(1)+GIOT*CSUMPLQ10(2)-CSUMPLQ10(3)
     &-GIOT*CSUMPLQ10(4)-(CSUMPLQ10(5)+GIOT*CSUMPLQ10(6)-CSUMPLQ10(7)
     &-GIOT*CSUMPLQ10(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM181(N4,ID,IK)=(CSUMPLQMOM10(5,IK)
     &+GIOT*CSUMPLQMOM10(6,IK)
     &-CSUMPLQMOM10(7,IK)-GIOT*CSUMPLQMOM10(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE182(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)
     &-CSUMPLQ10(4)
     &+(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))
     &+(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))
     &+(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM182(N4,ID,IK)=(CSUMPLQMOM10(1,IK)-CSUMPLQMOM10(2,IK)
     &+CSUMPLQMOM10(3,IK)-CSUMPLQMOM10(4,IK)
     &+(CSUMPLQMOM10(9,IK)-CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)
     &-CSUMPLQMOM10(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE183(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)
     &-CSUMPLQ10(4)
     &-(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))
     &+(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))
     &-(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM183(N4,ID,IK)=(CSUMPLQMOM10(5,IK)-CSUMPLQMOM10(6,IK)
     &+CSUMPLQMOM10(7,IK)-CSUMPLQMOM10(8,IK)
     &+(CSUMPLQMOM10(13,IK)-CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)
     &-CSUMPLQMOM10(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE184(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)
     &-CSUMPLQ10(4)
     &+(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))
     &-(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))
     &-(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM184(N4,ID,IK)=(CSUMPLQMOM10(1,IK)-CSUMPLQMOM10(2,IK)
     &+CSUMPLQMOM10(3,IK)-CSUMPLQMOM10(4,IK)
     &-(CSUMPLQMOM10(9,IK)-CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK)
     &-CSUMPLQMOM10(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE185(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3)
     &-CSUMPLQ10(4)
     &-(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8))
     &-(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12))
     &+(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM185(N4,ID,IK)=(CSUMPLQMOM10(5,IK)-CSUMPLQMOM10(6,IK)
     &+CSUMPLQMOM10(7,IK)-CSUMPLQMOM10(8,IK)
     &-(CSUMPLQMOM10(13,IK)-CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK)
     &-CSUMPLQMOM10(16,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 11
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE186(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)
     &+CSUMPLQ11(4)
     &+(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))
     &+(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))
     &+(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM186(N4,ID,IK)=(CSUMPLQMOM11(1,IK)+CSUMPLQMOM11(2,IK)
     &+CSUMPLQMOM11(3,IK)+CSUMPLQMOM11(4,IK)
     &+(CSUMPLQMOM11(9,IK)+CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)
     &+CSUMPLQMOM11(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE187(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)
     &+CSUMPLQ11(4)
     &-(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))
     &+(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))
     &-(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM187(N4,ID,IK)=(CSUMPLQMOM11(5,IK)+CSUMPLQMOM11(6,IK)
     &+CSUMPLQMOM11(7,IK)+CSUMPLQMOM11(8,IK)
     &+(CSUMPLQMOM11(13,IK)+CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)
     &+CSUMPLQMOM11(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE188(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)
     &+CSUMPLQ11(4)
     &+(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))
     &-(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))
     &-(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM188(N4,ID,IK)=(CSUMPLQMOM11(1,IK)+CSUMPLQMOM11(2,IK)
     &+CSUMPLQMOM11(3,IK)+CSUMPLQMOM11(4,IK)
     &-(CSUMPLQMOM11(9,IK)+CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)
     &+CSUMPLQMOM11(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE189(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3)
     &+CSUMPLQ11(4)
     &-(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8))
     &-(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12))
     &+(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM189(N4,ID,IK)=(CSUMPLQMOM11(5,IK)+CSUMPLQMOM11(6,IK)
     &+CSUMPLQMOM11(7,IK)+CSUMPLQMOM11(8,IK)
     &-(CSUMPLQMOM11(13,IK)+CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)
     &+CSUMPLQMOM11(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE190(N4,ID)=(CSUMPLQ11(1)+GIOT*CSUMPLQ11(2)-CSUMPLQ11(3)
     &-GIOT*CSUMPLQ11(4)+(CSUMPLQ11(5)+GIOT*CSUMPLQ11(6)-CSUMPLQ11(7)
     &-GIOT*CSUMPLQ11(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM190(N4,ID,IK)=(CSUMPLQMOM11(1,IK)
     &+GIOT*CSUMPLQMOM11(2,IK)
     &-CSUMPLQMOM11(3,IK)-GIOT*CSUMPLQMOM11(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE191(N4,ID)=(CSUMPLQ11(1)+GIOT*CSUMPLQ11(2)-CSUMPLQ11(3)
     &-GIOT*CSUMPLQ11(4)-(CSUMPLQ11(5)+GIOT*CSUMPLQ11(6)-CSUMPLQ11(7)
     &-GIOT*CSUMPLQ11(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM191(N4,ID,IK)=(CSUMPLQMOM11(5,IK)
     &+GIOT*CSUMPLQMOM11(6,IK)
     &-CSUMPLQMOM11(7,IK)-GIOT*CSUMPLQMOM11(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE192(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)
     &-CSUMPLQ11(4)
     &+(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))
     &+(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))
     &+(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM192(N4,ID,IK)=(CSUMPLQMOM11(1,IK)-CSUMPLQMOM11(2,IK)
     &+CSUMPLQMOM11(3,IK)-CSUMPLQMOM11(4,IK)
     &+(CSUMPLQMOM11(9,IK)-CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)
     &-CSUMPLQMOM11(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE193(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)
     &-CSUMPLQ11(4)
     &-(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))
     &+(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))
     &-(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM193(N4,ID,IK)=(CSUMPLQMOM11(5,IK)-CSUMPLQMOM11(6,IK)
     &+CSUMPLQMOM11(7,IK)-CSUMPLQMOM11(8,IK)
     &+(CSUMPLQMOM11(13,IK)-CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)
     &-CSUMPLQMOM11(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE194(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)
     &-CSUMPLQ11(4)
     &+(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))
     &-(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))
     &-(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM194(N4,ID,IK)=(CSUMPLQMOM11(1,IK)-CSUMPLQMOM11(2,IK)
     &+CSUMPLQMOM11(3,IK)-CSUMPLQMOM11(4,IK)
     &-(CSUMPLQMOM11(9,IK)-CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK)
     &-CSUMPLQMOM11(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE195(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3)
     &-CSUMPLQ11(4)
     &-(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8))
     &-(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12))
     &+(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM195(N4,ID,IK)=(CSUMPLQMOM11(5,IK)-CSUMPLQMOM11(6,IK)
     &+CSUMPLQMOM11(7,IK)-CSUMPLQMOM11(8,IK)
     &-(CSUMPLQMOM11(13,IK)-CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK)
     &-CSUMPLQMOM11(16,IK)))*ADIV2
      ENDDO      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 12
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE196(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)
     &+CSUMPLQ12(4)
     &+(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))
     &+(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))
     &+(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM196(N4,ID,IK)=(CSUMPLQMOM12(1,IK)+CSUMPLQMOM12(2,IK)
     &+CSUMPLQMOM12(3,IK)+CSUMPLQMOM12(4,IK)
     &+(CSUMPLQMOM12(9,IK)+CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)
     &+CSUMPLQMOM12(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE197(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)
     &+CSUMPLQ12(4)
     &-(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))
     &+(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))
     &-(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM197(N4,ID,IK)=(CSUMPLQMOM12(5,IK)+CSUMPLQMOM12(6,IK)
     &+CSUMPLQMOM12(7,IK)+CSUMPLQMOM12(8,IK)
     &+(CSUMPLQMOM12(13,IK)+CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)
     &+CSUMPLQMOM12(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE198(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)
     &+CSUMPLQ12(4)
     &+(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))
     &-(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))
     &-(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM198(N4,ID,IK)=(CSUMPLQMOM12(1,IK)+CSUMPLQMOM12(2,IK)
     &+CSUMPLQMOM12(3,IK)+CSUMPLQMOM12(4,IK)
     &-(CSUMPLQMOM12(9,IK)+CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)
     &+CSUMPLQMOM12(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE199(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3)
     &+CSUMPLQ12(4)
     &-(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8))
     &-(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12))
     &+(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM199(N4,ID,IK)=(CSUMPLQMOM12(5,IK)+CSUMPLQMOM12(6,IK)
     &+CSUMPLQMOM12(7,IK)+CSUMPLQMOM12(8,IK)
     &-(CSUMPLQMOM12(13,IK)+CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)
     &+CSUMPLQMOM12(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE200(N4,ID)=(CSUMPLQ12(1)+GIOT*CSUMPLQ12(2)-CSUMPLQ12(3)
     &-GIOT*CSUMPLQ12(4)+(CSUMPLQ12(5)+GIOT*CSUMPLQ12(6)-CSUMPLQ12(7)
     &-GIOT*CSUMPLQ12(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM200(N4,ID,IK)=(CSUMPLQMOM12(1,IK)
     &+GIOT*CSUMPLQMOM12(2,IK)
     &-CSUMPLQMOM12(3,IK)-GIOT*CSUMPLQMOM12(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE201(N4,ID)=(CSUMPLQ12(1)+GIOT*CSUMPLQ12(2)-CSUMPLQ12(3)
     &-GIOT*CSUMPLQ12(4)-(CSUMPLQ12(5)+GIOT*CSUMPLQ12(6)-CSUMPLQ12(7)
     &-GIOT*CSUMPLQ12(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM201(N4,ID,IK)=(CSUMPLQMOM12(5,IK)
     &+GIOT*CSUMPLQMOM12(6,IK)
     &-CSUMPLQMOM12(7,IK)-GIOT*CSUMPLQMOM12(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE202(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)
     &-CSUMPLQ12(4)
     &+(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))
     &+(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))
     &+(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM202(N4,ID,IK)=(CSUMPLQMOM12(1,IK)-CSUMPLQMOM12(2,IK)
     &+CSUMPLQMOM12(3,IK)-CSUMPLQMOM12(4,IK)
     &+(CSUMPLQMOM12(9,IK)-CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)
     &-CSUMPLQMOM12(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE203(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)
     &-CSUMPLQ12(4)
     &-(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))
     &+(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))
     &-(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM203(N4,ID,IK)=(CSUMPLQMOM12(5,IK)-CSUMPLQMOM12(6,IK)
     &+CSUMPLQMOM12(7,IK)-CSUMPLQMOM12(8,IK)
     &+(CSUMPLQMOM12(13,IK)-CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)
     &-CSUMPLQMOM12(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE204(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)
     &-CSUMPLQ12(4)
     &+(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))
     &-(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))
     &-(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM204(N4,ID,IK)=(CSUMPLQMOM12(1,IK)-CSUMPLQMOM12(2,IK)
     &+CSUMPLQMOM12(3,IK)-CSUMPLQMOM12(4,IK)
     &-(CSUMPLQMOM12(9,IK)-CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK)
     &-CSUMPLQMOM12(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE205(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3)
     &-CSUMPLQ12(4)
     &-(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8))
     &-(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12))
     &+(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM205(N4,ID,IK)=(CSUMPLQMOM12(5,IK)-CSUMPLQMOM12(6,IK)
     &+CSUMPLQMOM12(7,IK)-CSUMPLQMOM12(8,IK)
     &-(CSUMPLQMOM12(13,IK)-CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK)
     &-CSUMPLQMOM12(16,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 13
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE206(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)
     &+CSUMPLQ13(4)
     &+(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))
     &+(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))
     &+(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM206(N4,ID,IK)=(CSUMPLQMOM13(1,IK)+CSUMPLQMOM13(2,IK)
     &+CSUMPLQMOM13(3,IK)+CSUMPLQMOM13(4,IK)
     &+(CSUMPLQMOM13(9,IK)+CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)
     &+CSUMPLQMOM13(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE207(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)
     &+CSUMPLQ13(4)
     &-(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))
     &+(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))
     &-(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM207(N4,ID,IK)=(CSUMPLQMOM13(5,IK)+CSUMPLQMOM13(6,IK)
     &+CSUMPLQMOM13(7,IK)+CSUMPLQMOM13(8,IK)
     &+(CSUMPLQMOM13(13,IK)+CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)
     &+CSUMPLQMOM13(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE208(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)
     &+CSUMPLQ13(4)
     &+(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))
     &-(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))
     &-(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM208(N4,ID,IK)=(CSUMPLQMOM13(1,IK)+CSUMPLQMOM13(2,IK)
     &+CSUMPLQMOM13(3,IK)+CSUMPLQMOM13(4,IK)
     &-(CSUMPLQMOM13(9,IK)+CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)
     &+CSUMPLQMOM13(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE209(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3)
     &+CSUMPLQ13(4)
     &-(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8))
     &-(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12))
     &+(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM209(N4,ID,IK)=(CSUMPLQMOM13(5,IK)+CSUMPLQMOM13(6,IK)
     &+CSUMPLQMOM13(7,IK)+CSUMPLQMOM13(8,IK)
     &-(CSUMPLQMOM13(13,IK)+CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)
     &+CSUMPLQMOM13(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE210(N4,ID)=(CSUMPLQ13(1)+GIOT*CSUMPLQ13(2)-CSUMPLQ13(3)
     &-GIOT*CSUMPLQ13(4)+(CSUMPLQ13(5)+GIOT*CSUMPLQ13(6)-CSUMPLQ13(7)
     &-GIOT*CSUMPLQ13(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM210(N4,ID,IK)=(CSUMPLQMOM13(1,IK)
     &+GIOT*CSUMPLQMOM13(2,IK)
     &-CSUMPLQMOM13(3,IK)-GIOT*CSUMPLQMOM13(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE211(N4,ID)=(CSUMPLQ13(1)+GIOT*CSUMPLQ13(2)-CSUMPLQ13(3)
     &-GIOT*CSUMPLQ13(4)-(CSUMPLQ13(5)+GIOT*CSUMPLQ13(6)-CSUMPLQ13(7)
     &-GIOT*CSUMPLQ13(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM211(N4,ID,IK)=(CSUMPLQMOM13(5,IK)
     &+GIOT*CSUMPLQMOM13(6,IK)
     &-CSUMPLQMOM13(7,IK)-GIOT*CSUMPLQMOM13(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE212(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)
     &-CSUMPLQ13(4)
     &+(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))
     &+(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))
     &+(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM212(N4,ID,IK)=(CSUMPLQMOM13(1,IK)-CSUMPLQMOM13(2,IK)
     &+CSUMPLQMOM13(3,IK)-CSUMPLQMOM13(4,IK)
     &+(CSUMPLQMOM13(9,IK)-CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)
     &-CSUMPLQMOM13(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE213(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)
     &-CSUMPLQ13(4)
     &-(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))
     &+(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))
     &-(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM213(N4,ID,IK)=(CSUMPLQMOM13(5,IK)-CSUMPLQMOM13(6,IK)
     &+CSUMPLQMOM13(7,IK)-CSUMPLQMOM13(8,IK)
     &+(CSUMPLQMOM13(13,IK)-CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)
     &-CSUMPLQMOM13(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE214(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)
     &-CSUMPLQ13(4)
     &+(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))
     &-(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))
     &-(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM214(N4,ID,IK)=(CSUMPLQMOM13(1,IK)-CSUMPLQMOM13(2,IK)
     &+CSUMPLQMOM13(3,IK)-CSUMPLQMOM13(4,IK)
     &-(CSUMPLQMOM13(9,IK)-CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK)
     &-CSUMPLQMOM13(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE215(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3)
     &-CSUMPLQ13(4)
     &-(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8))
     &-(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12))
     &+(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM215(N4,ID,IK)=(CSUMPLQMOM13(5,IK)-CSUMPLQMOM13(6,IK)
     &+CSUMPLQMOM13(7,IK)-CSUMPLQMOM13(8,IK)
     &-(CSUMPLQMOM13(13,IK)-CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK)
     &-CSUMPLQMOM13(16,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 14
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE216(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)
     &+CSUMPLQ14(4)
     &+(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))
     &+(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))
     &+(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM216(N4,ID,IK)=(CSUMPLQMOM14(1,IK)+CSUMPLQMOM14(2,IK)
     &+CSUMPLQMOM14(3,IK)+CSUMPLQMOM14(4,IK)
     &+(CSUMPLQMOM14(9,IK)+CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)
     &+CSUMPLQMOM14(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE217(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)
     &+CSUMPLQ14(4)
     &-(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))
     &+(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))
     &-(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM217(N4,ID,IK)=(CSUMPLQMOM14(5,IK)+CSUMPLQMOM14(6,IK)
     &+CSUMPLQMOM14(7,IK)+CSUMPLQMOM14(8,IK)
     &+(CSUMPLQMOM14(13,IK)+CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)
     &+CSUMPLQMOM14(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE218(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)
     &+CSUMPLQ14(4)
     &+(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))
     &-(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))
     &-(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM218(N4,ID,IK)=(CSUMPLQMOM14(1,IK)+CSUMPLQMOM14(2,IK)
     &+CSUMPLQMOM14(3,IK)+CSUMPLQMOM14(4,IK)
     &-(CSUMPLQMOM14(9,IK)+CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)
     &+CSUMPLQMOM14(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE219(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3)
     &+CSUMPLQ14(4)
     &-(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8))
     &-(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12))
     &+(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM219(N4,ID,IK)=(CSUMPLQMOM14(5,IK)+CSUMPLQMOM14(6,IK)
     &+CSUMPLQMOM14(7,IK)+CSUMPLQMOM14(8,IK)
     &-(CSUMPLQMOM14(13,IK)+CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)
     &+CSUMPLQMOM14(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE220(N4,ID)=(CSUMPLQ14(1)+GIOT*CSUMPLQ14(2)-CSUMPLQ14(3)
     &-GIOT*CSUMPLQ14(4)+(CSUMPLQ14(5)+GIOT*CSUMPLQ14(6)-CSUMPLQ14(7)
     &-GIOT*CSUMPLQ14(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM220(N4,ID,IK)=(CSUMPLQMOM14(1,IK)
     &+GIOT*CSUMPLQMOM14(2,IK)
     &-CSUMPLQMOM14(3,IK)-GIOT*CSUMPLQMOM14(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE221(N4,ID)=(CSUMPLQ14(1)+GIOT*CSUMPLQ14(2)-CSUMPLQ14(3)
     &-GIOT*CSUMPLQ14(4)-(CSUMPLQ14(5)+GIOT*CSUMPLQ14(6)-CSUMPLQ14(7)
     &-GIOT*CSUMPLQ14(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM221(N4,ID,IK)=(CSUMPLQMOM14(5,IK)
     &+GIOT*CSUMPLQMOM14(6,IK)
     &-CSUMPLQMOM14(7,IK)-GIOT*CSUMPLQMOM14(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE222(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)
     &-CSUMPLQ14(4)
     &+(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))
     &+(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))
     &+(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM222(N4,ID,IK)=(CSUMPLQMOM14(1,IK)-CSUMPLQMOM14(2,IK)
     &+CSUMPLQMOM14(3,IK)-CSUMPLQMOM14(4,IK)
     &+(CSUMPLQMOM14(9,IK)-CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)
     &-CSUMPLQMOM14(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE223(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)
     &-CSUMPLQ14(4)
     &-(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))
     &+(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))
     &-(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM223(N4,ID,IK)=(CSUMPLQMOM14(5,IK)-CSUMPLQMOM14(6,IK)
     &+CSUMPLQMOM14(7,IK)-CSUMPLQMOM14(8,IK)
     &+(CSUMPLQMOM14(13,IK)-CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)
     &-CSUMPLQMOM14(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE224(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)
     &-CSUMPLQ14(4)
     &+(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))
     &-(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))
     &-(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM224(N4,ID,IK)=(CSUMPLQMOM14(1,IK)-CSUMPLQMOM14(2,IK)
     &+CSUMPLQMOM14(3,IK)-CSUMPLQMOM14(4,IK)
     &-(CSUMPLQMOM14(9,IK)-CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK)
     &-CSUMPLQMOM14(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE225(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3)
     &-CSUMPLQ14(4)
     &-(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8))
     &-(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12))
     &+(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM225(N4,ID,IK)=(CSUMPLQMOM14(5,IK)-CSUMPLQMOM14(6,IK)
     &+CSUMPLQMOM14(7,IK)-CSUMPLQMOM14(8,IK)
     &-(CSUMPLQMOM14(13,IK)-CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK)
     &-CSUMPLQMOM14(16,IK)))*ADIV2
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PLAQUETTE OPERATORS 15
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c**********************************************************************      
c     J=0, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE226(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)
     &+CSUMPLQ15(4)
     &+(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))
     &+(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))
     &+(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM226(N4,ID,IK)=(CSUMPLQMOM15(1,IK)+CSUMPLQMOM15(2,IK)
     &+CSUMPLQMOM15(3,IK)+CSUMPLQMOM15(4,IK)
     &+(CSUMPLQMOM15(9,IK)+CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)
     &+CSUMPLQMOM15(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE227(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)
     &+CSUMPLQ15(4)
     &-(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))
     &+(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))
     &-(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=+, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM227(N4,ID,IK)=(CSUMPLQMOM15(5,IK)+CSUMPLQMOM15(6,IK)
     &+CSUMPLQMOM15(7,IK)+CSUMPLQMOM15(8,IK)
     &+(CSUMPLQMOM15(13,IK)+CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)
     &+CSUMPLQMOM15(16,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE228(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)
     &+CSUMPLQ15(4)
     &+(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))
     &-(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))
     &-(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM228(N4,ID,IK)=(CSUMPLQMOM15(1,IK)+CSUMPLQMOM15(2,IK)
     &+CSUMPLQMOM15(3,IK)+CSUMPLQMOM15(4,IK)
     &-(CSUMPLQMOM15(9,IK)+CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)
     &+CSUMPLQMOM15(12,IK)))*ADIV2
      ENDDO
c**********************************************************************      
c     J=0, Pp=-, Pr=-, q=0
c**********************************************************************
      ALINE229(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3)
     &+CSUMPLQ15(4)
     &-(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8))
     &-(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12))
     &+(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=0, Pp=-, q=1,2
c**********************************************************************
      DO IK=1, 2
         ALINEMOM229(N4,ID,IK)=(CSUMPLQMOM15(5,IK)+CSUMPLQMOM15(6,IK)
     &+CSUMPLQMOM15(7,IK)+CSUMPLQMOM15(8,IK)
     &-(CSUMPLQMOM15(13,IK)+CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)
     &+CSUMPLQMOM15(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=1, Pr=+, q=0
c**********************************************************************            
      ALINE230(N4,ID)=(CSUMPLQ15(1)+GIOT*CSUMPLQ15(2)-CSUMPLQ15(3)
     &-GIOT*CSUMPLQ15(4)+(CSUMPLQ15(5)+GIOT*CSUMPLQ15(6)-CSUMPLQ15(7)
     &-GIOT*CSUMPLQ15(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2
c**********************************************************************                  
      DO IK=1, 2
         ALINEMOM230(N4,ID,IK)=(CSUMPLQMOM15(1,IK)
     &+GIOT*CSUMPLQMOM15(2,IK)
     &-CSUMPLQMOM15(3,IK)-GIOT*CSUMPLQMOM15(4,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=1, Pr=-, q=0
c**********************************************************************                  
      ALINE231(N4,ID)=(CSUMPLQ15(1)+GIOT*CSUMPLQ15(2)-CSUMPLQ15(3)
     &-GIOT*CSUMPLQ15(4)-(CSUMPLQ15(5)+GIOT*CSUMPLQ15(6)-CSUMPLQ15(7)
     &-GIOT*CSUMPLQ15(8)))*ADIV2
c**********************************************************************      
c     J=1, q=1,2  Here!
c**********************************************************************                        
      DO IK=1, 2
         ALINEMOM231(N4,ID,IK)=(CSUMPLQMOM15(5,IK)
     &+GIOT*CSUMPLQMOM15(6,IK)
     &-CSUMPLQMOM15(7,IK)-GIOT*CSUMPLQMOM15(8,IK))*ADIV1
      ENDDO
c**********************************************************************      
c     J=2, Pp=+, Pr=+, q=0
c**********************************************************************
      ALINE232(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)
     &-CSUMPLQ15(4)
     &+(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))
     &+(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))
     &+(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM232(N4,ID,IK)=(CSUMPLQMOM15(1,IK)-CSUMPLQMOM15(2,IK)
     &+CSUMPLQMOM15(3,IK)-CSUMPLQMOM15(4,IK)
     &+(CSUMPLQMOM15(9,IK)-CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)
     &-CSUMPLQMOM15(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=+, Pr=-, q=0
c**********************************************************************
      ALINE233(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)
     &-CSUMPLQ15(4)
     &-(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))
     &+(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))
     &-(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=+, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM233(N4,ID,IK)=(CSUMPLQMOM15(5,IK)-CSUMPLQMOM15(6,IK)
     &+CSUMPLQMOM15(7,IK)-CSUMPLQMOM15(8,IK)
     &+(CSUMPLQMOM15(13,IK)-CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)
     &-CSUMPLQMOM15(16,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=+, q=0
c**********************************************************************
      ALINE234(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)
     &-CSUMPLQ15(4)
     &+(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))
     &-(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))
     &-(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM234(N4,ID,IK)=(CSUMPLQMOM15(1,IK)-CSUMPLQMOM15(2,IK)
     &+CSUMPLQMOM15(3,IK)-CSUMPLQMOM15(4,IK)
     &-(CSUMPLQMOM15(9,IK)-CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK)
     &-CSUMPLQMOM15(12,IK)))*ADIV2
      ENDDO      
c**********************************************************************      
c     J=2, Pp=-, Pr=-, q=0  test
c**********************************************************************
      ALINE235(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3)
     &-CSUMPLQ15(4)
     &-(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8))
     &-(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12))
     &+(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
c**********************************************************************      
c     J=2, Pp=-, q=0
c**********************************************************************
      DO IK=1, 2
         ALINEMOM235(N4,ID,IK)=(CSUMPLQMOM15(5,IK)-CSUMPLQMOM15(6,IK)
     &+CSUMPLQMOM15(7,IK)-CSUMPLQMOM15(8,IK)
     &-(CSUMPLQMOM15(13,IK)-CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK)
     &-CSUMPLQMOM15(16,IK)))*ADIV2
      ENDDO      
C**********************************************************************
 2    CONTINUE
C**********************************************************************
      RETURN
      END
C**********************************************************************
C**********************************************************************
C     HERE I CREATE THE DIAGONAL OPERATORS
C**********************************************************************
C********************************************************************
C   HERE, WE CREATE HERMITIAN CORRELATORS!!!!                  
C   COR(I,J)=F(I,T+DT)*F^+(J,T)+F^+(J,T+DT)*F(I,T)             
C   COR(J,I)=F(J,T+DT)*F^+(I,T)+F^+(I,T+DT)*F(J,T)             
C   COR^+(J,I)=F^+(J,T+DT)*F(I,T)+F(I,T+DT)*F^+(J,T)           
C   => COR(I,J)=COR^+(J,I)!!                                   
C*********************************************************************
      SUBROUTINE POT
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NUMBIN=2,LMAX=LX4/2+1)
      PARAMETER(LMAXIR=3)
C******************************************************************************
C******************************************************************************     
      PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
      PARAMETER(I1P=34,I1M=32) !TBC!
      PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!      
C----------------------------------------------
      PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
      PARAMETER(I1MOM=66) !TBC!
      PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
C-----------------------------------
      PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
      PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
      PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
      PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
      PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
C-----------------------------------------------------
      PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
      PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
      PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
C-------------------------------------------------------------
      PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
      PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
      PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
      PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
C-------------------------------------------------
      PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
      PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
      PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
      PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
C----------------------------------------------------------
C******************************************************************************
C******************************************************************************
C******************************************************************************  
C     
      COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLP(NUMBIN,NTOTAL)
C
      COMMON/PBLOCKMOMP/ACORLPMOM(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM)
     &,AVACLPMOM(2,NUMBIN,NTOTALMOM)
C
      COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP)
     &,AVACLJ0PP(NUMBIN,NOPJ0PP)
      COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM)
     &,AVACLJ0PM(NUMBIN,NOPJ0PM)
      COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP)
     &,AVACLJ0MP(NUMBIN,NOPJ0MP)
      COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM)
     &,AVACLJ0MM(NUMBIN,NOPJ0MM)
      COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P)
     &,AVACLJ1P(NUMBIN,NOPJ1P)
      COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M)
     &,AVACLJ1M(NUMBIN,NOPJ1M)
      COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP)
     &,AVACLJ2PP(NUMBIN,NOPJ2PP)
      COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM)
     &,AVACLJ2PM(NUMBIN,NOPJ2PM)
      COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP)
     &,AVACLJ2MP(NUMBIN,NOPJ2MP)
      COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM)
     &,AVACLJ2MM(NUMBIN,NOPJ2MM)
C
      COMMON/PBLOCKJ0/ACORLJ0(NUMBIN,LMAX,NOPFULJ0,NOPFULJ0)
     &,AVACLJ0(NUMBIN,NOPFULJ0)
      COMMON/PBLOCKJ1/ACORLJ1(NUMBIN,LMAX,NOPFULJ1,NOPFULJ1)
     &,AVACLJ1(NUMBIN,NOPFULJ1)      
      COMMON/PBLOCKJ2/ACORLJ2(NUMBIN,LMAX,NOPFULJ2,NOPFULJ2)
     &,AVACLJ2(NUMBIN,NOPFULJ2)
C
      COMMON/PBLOCKJALL/ACORLALL(NUMBIN,LMAXIR,NTOTAL,NTOTAL)
     &,AVACLALL(NUMBIN,NTOTAL)      
C
      COMMON/PBLOCKMOMJ0P/
     &ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM)
     &,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
      COMMON/PBLOCKMOMJ0M/
     &ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM)
     &,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)      
      COMMON/PBLOCKMOMJ1/
     &ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM)
     &,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)     
      COMMON/PBLOCKMOMJ2P/
     &ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM)
     &,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
      COMMON/PBLOCKMOMJ2M/
     &ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM)
     &,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
C
      COMMON/PBLOCKMOMJ0/
     &ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM)
     &,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
       COMMON/PBLOCKMOMJ2/
     &ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM)
     &,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)
C
       COMMON/PBLOCKMOMJALL/
     &ACORLPMOMJALL(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM)
     &,AVACLMOMALL(2,NUMBIN,NTOTALMOM) 
C      
C***************************************************************
      COMMON/LINES/ALINE1(LX4,IBLOK),ALINE2(LX4,IBLOK)
     &,ALINE3(LX4,IBLOK),ALINE4(LX4,IBLOK),ALINE5(LX4,IBLOK)
     &,ALINE6(LX4,IBLOK),ALINE7(LX4,IBLOK),ALINE8(LX4,IBLOK) 
     &,ALINE9(LX4,IBLOK),ALINE10(LX4,IBLOK),ALINE11(LX4,IBLOK)
     &,ALINE12(LX4,IBLOK),ALINE13(LX4,IBLOK),ALINE14(LX4,IBLOK)
     &,ALINE15(LX4,IBLOK),ALINE16(LX4,IBLOK),ALINE17(LX4,IBLOK)
     &,ALINE18(LX4,IBLOK),ALINE19(LX4,IBLOK),ALINE20(LX4,IBLOK)
     &,ALINE21(LX4,IBLOK),ALINE22(LX4,IBLOK),ALINE23(LX4,IBLOK)
     &,ALINE24(LX4,IBLOK),ALINE25(LX4,IBLOK),ALINE26(LX4,IBLOK)
     &,ALINE27(LX4,IBLOK),ALINE28(LX4,IBLOK),ALINE29(LX4,IBLOK)
     &,ALINE30(LX4,IBLOK),ALINE31(LX4,IBLOK),ALINE32(LX4,IBLOK)
     &,ALINE33(LX4,IBLOK),ALINE34(LX4,IBLOK),ALINE35(LX4,IBLOK)
     &,ALINE36(LX4,IBLOK),ALINE37(LX4,IBLOK),ALINE38(LX4,IBLOK)
     &,ALINE39(LX4,IBLOK),ALINE40(LX4,IBLOK),ALINE41(LX4,IBLOK)
     &,ALINE42(LX4,IBLOK),ALINE43(LX4,IBLOK),ALINE44(LX4,IBLOK)
     &,ALINE45(LX4,IBLOK),ALINE46(LX4,IBLOK),ALINE47(LX4,IBLOK)
     &,ALINE48(LX4,IBLOK),ALINE49(LX4,IBLOK),ALINE50(LX4,IBLOK)
     &,ALINE51(LX4,IBLOK),ALINE52(LX4,IBLOK),ALINE53(LX4,IBLOK)
     &,ALINE54(LX4,IBLOK),ALINE55(LX4,IBLOK),ALINE56(LX4,IBLOK)
     &,ALINE57(LX4,IBLOK),ALINE58(LX4,IBLOK),ALINE59(LX4,IBLOK)
     &,ALINE60(LX4,IBLOK),ALINE61(LX4,IBLOK),ALINE62(LX4,IBLOK)
     &,ALINE63(LX4,IBLOK),ALINE64(LX4,IBLOK),ALINE65(LX4,IBLOK)
     &,ALINE66(LX4,IBLOK),ALINE67(LX4,IBLOK),ALINE68(LX4,IBLOK)
     &,ALINE69(LX4,IBLOK),ALINE70(LX4,IBLOK),ALINE71(LX4,IBLOK)
     &,ALINE72(LX4,IBLOK),ALINE73(LX4,IBLOK),ALINE74(LX4,IBLOK)
     &,ALINE75(LX4,IBLOK),ALINE76(LX4,IBLOK),ALINE77(LX4,IBLOK)
     &,ALINE78(LX4,IBLOK),ALINE79(LX4,IBLOK),ALINE80(LX4,IBLOK)
     &,ALINE81(LX4,IBLOK),ALINE82(LX4,IBLOK),ALINE83(LX4,IBLOK)
     &,ALINE84(LX4,IBLOK),ALINE85(LX4,IBLOK),ALINE86(LX4,IBLOK)
     &,ALINE87(LX4,IBLOK),ALINE88(LX4,IBLOK),ALINE89(LX4,IBLOK)
     &,ALINE90(LX4,IBLOK),ALINE91(LX4,IBLOK),ALINE92(LX4,IBLOK)
     &,ALINE93(LX4,IBLOK),ALINE94(LX4,IBLOK),ALINE95(LX4,IBLOK)
     &,ALINE96(LX4,IBLOK),ALINE97(LX4,IBLOK),ALINE98(LX4,IBLOK)
     &,ALINE99(LX4,IBLOK),ALINE100(LX4,IBLOK),ALINE101(LX4,IBLOK)
     &,ALINE102(LX4,IBLOK),ALINE103(LX4,IBLOK),ALINE104(LX4,IBLOK)
     &,ALINE105(LX4,IBLOK),ALINE106(LX4,IBLOK),ALINE107(LX4,IBLOK)
     &,ALINE108(LX4,IBLOK),ALINE109(LX4,IBLOK),ALINE110(LX4,IBLOK)
     &,ALINE111(LX4,IBLOK),ALINE112(LX4,IBLOK),ALINE113(LX4,IBLOK)
     &,ALINE114(LX4,IBLOK),ALINE115(LX4,IBLOK),ALINE116(LX4,IBLOK)
     &,ALINE117(LX4,IBLOK),ALINE118(LX4,IBLOK),ALINE119(LX4,IBLOK)
     &,ALINE120(LX4,IBLOK),ALINE121(LX4,IBLOK),ALINE122(LX4,IBLOK)
     &,ALINE123(LX4,IBLOK),ALINE124(LX4,IBLOK),ALINE125(LX4,IBLOK)
     &,ALINE126(LX4,IBLOK),ALINE127(LX4,IBLOK),ALINE128(LX4,IBLOK)
     &,ALINE129(LX4,IBLOK),ALINE130(LX4,IBLOK),ALINE131(LX4,IBLOK)
     &,ALINE132(LX4,IBLOK),ALINE133(LX4,IBLOK),ALINE134(LX4,IBLOK)
     &,ALINE135(LX4,IBLOK),ALINE136(LX4,IBLOK),ALINE137(LX4,IBLOK)
     &,ALINE138(LX4,IBLOK),ALINE139(LX4,IBLOK),ALINE140(LX4,IBLOK)
     &,ALINE141(LX4,IBLOK),ALINE142(LX4,IBLOK),ALINE143(LX4,IBLOK)
     &,ALINE144(LX4,IBLOK),ALINE145(LX4,IBLOK),ALINE146(LX4,IBLOK)
     &,ALINE147(LX4,IBLOK),ALINE148(LX4,IBLOK),ALINE149(LX4,IBLOK)
     &,ALINE150(LX4,IBLOK),ALINE151(LX4,IBLOK),ALINE152(LX4,IBLOK)
     &,ALINE153(LX4,IBLOK),ALINE154(LX4,IBLOK),ALINE155(LX4,IBLOK)
     &,ALINE156(LX4,IBLOK),ALINE157(LX4,IBLOK),ALINE158(LX4,IBLOK)
     &,ALINE159(LX4,IBLOK),ALINE160(LX4,IBLOK),ALINE161(LX4,IBLOK)
     &,ALINE162(LX4,IBLOK),ALINE163(LX4,IBLOK),ALINE164(LX4,IBLOK)      
     &,ALINE165(LX4,IBLOK),ALINE166(LX4,IBLOK),ALINE167(LX4,IBLOK)
     &,ALINE168(LX4,IBLOK),ALINE169(LX4,IBLOK),ALINE170(LX4,IBLOK)      
     &,ALINE171(LX4,IBLOK),ALINE172(LX4,IBLOK),ALINE173(LX4,IBLOK)      
     &,ALINE174(LX4,IBLOK),ALINE175(LX4,IBLOK),ALINE176(LX4,IBLOK)      
     &,ALINE177(LX4,IBLOK),ALINE178(LX4,IBLOK),ALINE179(LX4,IBLOK)      
     &,ALINE180(LX4,IBLOK),ALINE181(LX4,IBLOK),ALINE182(LX4,IBLOK)      
     &,ALINE183(LX4,IBLOK),ALINE184(LX4,IBLOK),ALINE185(LX4,IBLOK)
     &,ALINE186(LX4,IBLOK),ALINE187(LX4,IBLOK),ALINE188(LX4,IBLOK)            
     &,ALINE189(LX4,IBLOK),ALINE190(LX4,IBLOK),ALINE191(LX4,IBLOK)      
     &,ALINE192(LX4,IBLOK),ALINE193(LX4,IBLOK),ALINE194(LX4,IBLOK)      
     &,ALINE195(LX4,IBLOK),ALINE196(LX4,IBLOK),ALINE197(LX4,IBLOK)      
     &,ALINE198(LX4,IBLOK),ALINE199(LX4,IBLOK),ALINE200(LX4,IBLOK)
     &,ALINE201(LX4,IBLOK),ALINE202(LX4,IBLOK),ALINE203(LX4,IBLOK)
     &,ALINE204(LX4,IBLOK),ALINE205(LX4,IBLOK),ALINE206(LX4,IBLOK)
     &,ALINE207(LX4,IBLOK),ALINE208(LX4,IBLOK),ALINE209(LX4,IBLOK)      
     &,ALINE210(LX4,IBLOK),ALINE211(LX4,IBLOK),ALINE212(LX4,IBLOK)
     &,ALINE213(LX4,IBLOK),ALINE214(LX4,IBLOK),ALINE215(LX4,IBLOK)
     &,ALINE216(LX4,IBLOK),ALINE217(LX4,IBLOK),ALINE218(LX4,IBLOK)      
     &,ALINE219(LX4,IBLOK),ALINE220(LX4,IBLOK),ALINE221(LX4,IBLOK)
     &,ALINE222(LX4,IBLOK),ALINE223(LX4,IBLOK),ALINE224(LX4,IBLOK)
     &,ALINE225(LX4,IBLOK),ALINE226(LX4,IBLOK),ALINE227(LX4,IBLOK)
     &,ALINE228(LX4,IBLOK),ALINE229(LX4,IBLOK),ALINE230(LX4,IBLOK)
     &,ALINE231(LX4,IBLOK),ALINE232(LX4,IBLOK),ALINE233(LX4,IBLOK)
     &,ALINE234(LX4,IBLOK),ALINE235(LX4,IBLOK) 
C********************************************************
      COMMON/LINESMOM/ALINEMOM2(LX4,IBLOK,2),ALINEMOM3(LX4,IBLOK,2)
     &,ALINEMOM4(LX4,IBLOK,2),ALINEMOM5(LX4,IBLOK,2)
     &,ALINEMOM6(LX4,IBLOK,2),ALINEMOM7(LX4,IBLOK,2)
     &,ALINEMOM8(LX4,IBLOK,2),ALINEMOM9(LX4,IBLOK,2)
     &,ALINEMOM10(LX4,IBLOK,2),ALINEMOM11(LX4,IBLOK,2)
     &,ALINEMOM12(LX4,IBLOK,2),ALINEMOM13(LX4,IBLOK,2)
     &,ALINEMOM14(LX4,IBLOK,2),ALINEMOM15(LX4,IBLOK,2) 
     &,ALINEMOM16(LX4,IBLOK,2),ALINEMOM17(LX4,IBLOK,2)
     &,ALINEMOM18(LX4,IBLOK,2),ALINEMOM19(LX4,IBLOK,2)
     &,ALINEMOM20(LX4,IBLOK,2),ALINEMOM21(LX4,IBLOK,2)
     &,ALINEMOM22(LX4,IBLOK,2),ALINEMOM23(LX4,IBLOK,2)
     &,ALINEMOM24(LX4,IBLOK,2),ALINEMOM25(LX4,IBLOK,2)
     &,ALINEMOM26(LX4,IBLOK,2),ALINEMOM27(LX4,IBLOK,2)
     &,ALINEMOM28(LX4,IBLOK,2),ALINEMOM29(LX4,IBLOK,2)
     &,ALINEMOM30(LX4,IBLOK,2),ALINEMOM31(LX4,IBLOK,2)
     &,ALINEMOM32(LX4,IBLOK,2),ALINEMOM33(LX4,IBLOK,2)
     &,ALINEMOM34(LX4,IBLOK,2),ALINEMOM35(LX4,IBLOK,2)
     &,ALINEMOM36(LX4,IBLOK,2),ALINEMOM37(LX4,IBLOK,2)
     &,ALINEMOM38(LX4,IBLOK,2),ALINEMOM39(LX4,IBLOK,2)
     &,ALINEMOM40(LX4,IBLOK,2),ALINEMOM41(LX4,IBLOK,2)
     &,ALINEMOM42(LX4,IBLOK,2),ALINEMOM43(LX4,IBLOK,2)
     &,ALINEMOM44(LX4,IBLOK,2),ALINEMOM45(LX4,IBLOK,2)
     &,ALINEMOM46(LX4,IBLOK,2),ALINEMOM47(LX4,IBLOK,2)
     &,ALINEMOM48(LX4,IBLOK,2),ALINEMOM49(LX4,IBLOK,2)
     &,ALINEMOM50(LX4,IBLOK,2),ALINEMOM51(LX4,IBLOK,2)
     &,ALINEMOM52(LX4,IBLOK,2),ALINEMOM53(LX4,IBLOK,2)
     &,ALINEMOM54(LX4,IBLOK,2),ALINEMOM55(LX4,IBLOK,2)
     &,ALINEMOM56(LX4,IBLOK,2),ALINEMOM57(LX4,IBLOK,2)
     &,ALINEMOM58(LX4,IBLOK,2),ALINEMOM59(LX4,IBLOK,2)
     &,ALINEMOM60(LX4,IBLOK,2),ALINEMOM61(LX4,IBLOK,2)
     &,ALINEMOM62(LX4,IBLOK,2),ALINEMOM63(LX4,IBLOK,2)
     &,ALINEMOM64(LX4,IBLOK,2)
     &,ALINEMOM65(LX4,IBLOK,2),ALINEMOM66(LX4,IBLOK,2)
     &,ALINEMOM67(LX4,IBLOK,2),ALINEMOM68(LX4,IBLOK,2)
     &,ALINEMOM69(LX4,IBLOK,2),ALINEMOM70(LX4,IBLOK,2)
     &,ALINEMOM71(LX4,IBLOK,2),ALINEMOM72(LX4,IBLOK,2)
     &,ALINEMOM73(LX4,IBLOK,2),ALINEMOM74(LX4,IBLOK,2) 
     &,ALINEMOM75(LX4,IBLOK,2),ALINEMOM76(LX4,IBLOK,2)
     &,ALINEMOM77(LX4,IBLOK,2),ALINEMOM78(LX4,IBLOK,2)
     &,ALINEMOM79(LX4,IBLOK,2),ALINEMOM80(LX4,IBLOK,2)
     &,ALINEMOM81(LX4,IBLOK,2),ALINEMOM82(LX4,IBLOK,2)
     &,ALINEMOM83(LX4,IBLOK,2),ALINEMOM84(LX4,IBLOK,2) 
     &,ALINEMOM85(LX4,IBLOK,2),ALINEMOM86(LX4,IBLOK,2)
     &,ALINEMOM87(LX4,IBLOK,2),ALINEMOM88(LX4,IBLOK,2)
     &,ALINEMOM89(LX4,IBLOK,2),ALINEMOM90(LX4,IBLOK,2)
     &,ALINEMOM91(LX4,IBLOK,2),ALINEMOM92(LX4,IBLOK,2)
     &,ALINEMOM93(LX4,IBLOK,2),ALINEMOM94(LX4,IBLOK,2) 
     &,ALINEMOM95(LX4,IBLOK,2),ALINEMOM96(LX4,IBLOK,2)
     &,ALINEMOM97(LX4,IBLOK,2),ALINEMOM98(LX4,IBLOK,2)
     &,ALINEMOM99(LX4,IBLOK,2),ALINEMOM100(LX4,IBLOK,2)
     &,ALINEMOM101(LX4,IBLOK,2),ALINEMOM102(LX4,IBLOK,2)
     &,ALINEMOM103(LX4,IBLOK,2),ALINEMOM104(LX4,IBLOK,2) 
     &,ALINEMOM105(LX4,IBLOK,2),ALINEMOM106(LX4,IBLOK,2)
     &,ALINEMOM107(LX4,IBLOK,2),ALINEMOM108(LX4,IBLOK,2)
     &,ALINEMOM109(LX4,IBLOK,2),ALINEMOM110(LX4,IBLOK,2)
     &,ALINEMOM111(LX4,IBLOK,2),ALINEMOM112(LX4,IBLOK,2)
     &,ALINEMOM113(LX4,IBLOK,2),ALINEMOM114(LX4,IBLOK,2) 
     &,ALINEMOM115(LX4,IBLOK,2),ALINEMOM116(LX4,IBLOK,2)
     &,ALINEMOM117(LX4,IBLOK,2),ALINEMOM118(LX4,IBLOK,2)
     &,ALINEMOM119(LX4,IBLOK,2),ALINEMOM120(LX4,IBLOK,2)
     &,ALINEMOM121(LX4,IBLOK,2),ALINEMOM122(LX4,IBLOK,2)
     &,ALINEMOM123(LX4,IBLOK,2),ALINEMOM124(LX4,IBLOK,2) 
     &,ALINEMOM125(LX4,IBLOK,2),ALINEMOM126(LX4,IBLOK,2)
     &,ALINEMOM127(LX4,IBLOK,2),ALINEMOM128(LX4,IBLOK,2)
     &,ALINEMOM129(LX4,IBLOK,2),ALINEMOM130(LX4,IBLOK,2)
     &,ALINEMOM131(LX4,IBLOK,2),ALINEMOM132(LX4,IBLOK,2)
     &,ALINEMOM133(LX4,IBLOK,2),ALINEMOM134(LX4,IBLOK,2) 
     &,ALINEMOM135(LX4,IBLOK,2),ALINEMOM136(LX4,IBLOK,2)
     &,ALINEMOM137(LX4,IBLOK,2),ALINEMOM138(LX4,IBLOK,2)
     &,ALINEMOM139(LX4,IBLOK,2),ALINEMOM140(LX4,IBLOK,2)
     &,ALINEMOM141(LX4,IBLOK,2),ALINEMOM142(LX4,IBLOK,2)
     &,ALINEMOM143(LX4,IBLOK,2),ALINEMOM144(LX4,IBLOK,2) 
     &,ALINEMOM145(LX4,IBLOK,2),ALINEMOM146(LX4,IBLOK,2)
     &,ALINEMOM147(LX4,IBLOK,2),ALINEMOM148(LX4,IBLOK,2)
     &,ALINEMOM149(LX4,IBLOK,2),ALINEMOM150(LX4,IBLOK,2)      
     &,ALINEMOM151(LX4,IBLOK,2),ALINEMOM152(LX4,IBLOK,2)      
     &,ALINEMOM153(LX4,IBLOK,2),ALINEMOM154(LX4,IBLOK,2)      
     &,ALINEMOM155(LX4,IBLOK,2),ALINEMOM156(LX4,IBLOK,2)      
     &,ALINEMOM157(LX4,IBLOK,2),ALINEMOM158(LX4,IBLOK,2)      
     &,ALINEMOM159(LX4,IBLOK,2),ALINEMOM160(LX4,IBLOK,2)
     &,ALINEMOM161(LX4,IBLOK,2),ALINEMOM162(LX4,IBLOK,2)       
     &,ALINEMOM163(LX4,IBLOK,2),ALINEMOM164(LX4,IBLOK,2)       
     &,ALINEMOM165(LX4,IBLOK,2),ALINEMOM166(LX4,IBLOK,2)
     &,ALINEMOM167(LX4,IBLOK,2),ALINEMOM168(LX4,IBLOK,2)      
     &,ALINEMOM169(LX4,IBLOK,2),ALINEMOM170(LX4,IBLOK,2)
     &,ALINEMOM171(LX4,IBLOK,2),ALINEMOM172(LX4,IBLOK,2)
     &,ALINEMOM173(LX4,IBLOK,2),ALINEMOM174(LX4,IBLOK,2)
     &,ALINEMOM175(LX4,IBLOK,2),ALINEMOM176(LX4,IBLOK,2)
     &,ALINEMOM177(LX4,IBLOK,2),ALINEMOM178(LX4,IBLOK,2)
     &,ALINEMOM179(LX4,IBLOK,2),ALINEMOM180(LX4,IBLOK,2)
     &,ALINEMOM181(LX4,IBLOK,2),ALINEMOM182(LX4,IBLOK,2)
     &,ALINEMOM183(LX4,IBLOK,2),ALINEMOM184(LX4,IBLOK,2)
     &,ALINEMOM185(LX4,IBLOK,2),ALINEMOM186(LX4,IBLOK,2)
     &,ALINEMOM187(LX4,IBLOK,2),ALINEMOM188(LX4,IBLOK,2)
     &,ALINEMOM189(LX4,IBLOK,2),ALINEMOM190(LX4,IBLOK,2)
     &,ALINEMOM191(LX4,IBLOK,2),ALINEMOM192(LX4,IBLOK,2)
     &,ALINEMOM193(LX4,IBLOK,2),ALINEMOM194(LX4,IBLOK,2)
     &,ALINEMOM195(LX4,IBLOK,2),ALINEMOM196(LX4,IBLOK,2)
     &,ALINEMOM197(LX4,IBLOK,2),ALINEMOM198(LX4,IBLOK,2)
     &,ALINEMOM199(LX4,IBLOK,2),ALINEMOM200(LX4,IBLOK,2)
     &,ALINEMOM201(LX4,IBLOK,2),ALINEMOM202(LX4,IBLOK,2)
     &,ALINEMOM203(LX4,IBLOK,2),ALINEMOM204(LX4,IBLOK,2) 
     &,ALINEMOM205(LX4,IBLOK,2),ALINEMOM206(LX4,IBLOK,2)
     &,ALINEMOM207(LX4,IBLOK,2),ALINEMOM208(LX4,IBLOK,2)
     &,ALINEMOM209(LX4,IBLOK,2),ALINEMOM210(LX4,IBLOK,2)
     &,ALINEMOM211(LX4,IBLOK,2),ALINEMOM212(LX4,IBLOK,2)
     &,ALINEMOM213(LX4,IBLOK,2),ALINEMOM214(LX4,IBLOK,2) 
     &,ALINEMOM215(LX4,IBLOK,2),ALINEMOM216(LX4,IBLOK,2)
     &,ALINEMOM217(LX4,IBLOK,2),ALINEMOM218(LX4,IBLOK,2)
     &,ALINEMOM219(LX4,IBLOK,2),ALINEMOM220(LX4,IBLOK,2)
     &,ALINEMOM221(LX4,IBLOK,2),ALINEMOM222(LX4,IBLOK,2)
     &,ALINEMOM223(LX4,IBLOK,2),ALINEMOM224(LX4,IBLOK,2) 
     &,ALINEMOM225(LX4,IBLOK,2),ALINEMOM226(LX4,IBLOK,2)
     &,ALINEMOM227(LX4,IBLOK,2),ALINEMOM228(LX4,IBLOK,2)
     &,ALINEMOM229(LX4,IBLOK,2),ALINEMOM230(LX4,IBLOK,2)
     &,ALINEMOM231(LX4,IBLOK,2),ALINEMOM232(LX4,IBLOK,2)
     &,ALINEMOM233(LX4,IBLOK,2),ALINEMOM234(LX4,IBLOK,2) 
     &,ALINEMOM235(LX4,IBLOK,2)
C***************************************************************
      COMMON/ITEM/ITR,IBIN,ITOT
C***********************************************************
      COMPLEX*16 ALLP
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMPLEX*16 ALLPJ0PP, ALLPJ0PM, ALLPJ0MP, ALLPJ0MM
      COMPLEX*16 ALLPJ2PP, ALLPJ2PM, ALLPJ2MP, ALLPJ2MM
      COMPLEX*16 ALLPJ1P, ALLPJ1M
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMPLEX*16 ALLPJ0, ALLPJ1, ALLPJ2, ALLPJALL
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      COMPLEX*16 ALINE1,ALINE2,ALINE3,ALINE4,ALINE5,ALINE6,ALINE7
      COMPLEX*16 ALINE8,ALINE9,ALINE10,ALINE11,ALINE12,ALINE13
      COMPLEX*16 ALINE14,ALINE15,ALINE16,ALINE17,ALINE18,ALINE19
      COMPLEX*16 ALINE20,ALINE21,ALINE22,ALINE23,ALINE24,ALINE25
      COMPLEX*16 ALINE26,ALINE27,ALINE28,ALINE29,ALINE30
      COMPLEX*16 ALINE31,ALINE32,ALINE33,ALINE34,ALINE35
      COMPLEX*16 ALINE36,ALINE37,ALINE38,ALINE39,ALINE40
      COMPLEX*16 ALINE41,ALINE42,ALINE43,ALINE44,ALINE45
      COMPLEX*16 ALINE46,ALINE47,ALINE48,ALINE49,ALINE50
      COMPLEX*16 ALINE51,ALINE52,ALINE53,ALINE54,ALINE55
      COMPLEX*16 ALINE56,ALINE57,ALINE58,ALINE59,ALINE60
      COMPLEX*16 ALINE61,ALINE62,ALINE63,ALINE64,ALINE65
      COMPLEX*16 ALINE66,ALINE67,ALINE68,ALINE69,ALINE70
      COMPLEX*16 ALINE71,ALINE72,ALINE73,ALINE74,ALINE75
      COMPLEX*16 ALINE76,ALINE77,ALINE78,ALINE79,ALINE80
      COMPLEX*16 ALINE81,ALINE82,ALINE83,ALINE84,ALINE85
      COMPLEX*16 ALINE86,ALINE87,ALINE88,ALINE89,ALINE90
      COMPLEX*16 ALINE91,ALINE92,ALINE93,ALINE94,ALINE95
      COMPLEX*16 ALINE96,ALINE97,ALINE98,ALINE99,ALINE100
      COMPLEX*16 ALINE101,ALINE102,ALINE103,ALINE104,ALINE105
      COMPLEX*16 ALINE106,ALINE107,ALINE108,ALINE109,ALINE110
      COMPLEX*16 ALINE111,ALINE112,ALINE113,ALINE114,ALINE115
      COMPLEX*16 ALINE116,ALINE117,ALINE118,ALINE119,ALINE120
      COMPLEX*16 ALINE121,ALINE122,ALINE123,ALINE124,ALINE125
      COMPLEX*16 ALINE126,ALINE127,ALINE128,ALINE129,ALINE130
      COMPLEX*16 ALINE131,ALINE132,ALINE133,ALINE134,ALINE135
      COMPLEX*16 ALINE136,ALINE137,ALINE138,ALINE139,ALINE140
      COMPLEX*16 ALINE141,ALINE142,ALINE143,ALINE144,ALINE145
      COMPLEX*16 ALINE146,ALINE147,ALINE148,ALINE149,ALINE150
      COMPLEX*16 ALINE151,ALINE152,ALINE153,ALINE154,ALINE155
      COMPLEX*16 ALINE156,ALINE157,ALINE158,ALINE159,ALINE160
      COMPLEX*16 ALINE161,ALINE162,ALINE163,ALINE164,ALINE165
      COMPLEX*16 ALINE166,ALINE167,ALINE168,ALINE169,ALINE170
      COMPLEX*16 ALINE171,ALINE172,ALINE173,ALINE174,ALINE175
      COMPLEX*16 ALINE176,ALINE177,ALINE178,ALINE179,ALINE180
      COMPLEX*16 ALINE181,ALINE182,ALINE183,ALINE184,ALINE185
      COMPLEX*16 ALINE186,ALINE187,ALINE188,ALINE189,ALINE190
      COMPLEX*16 ALINE191,ALINE192,ALINE193,ALINE194,ALINE195
      COMPLEX*16 ALINE196,ALINE197,ALINE198,ALINE199,ALINE200
      COMPLEX*16 ALINE201,ALINE202,ALINE203,ALINE204,ALINE205
      COMPLEX*16 ALINE206,ALINE207,ALINE208,ALINE209,ALINE210
      COMPLEX*16 ALINE211,ALINE212,ALINE213,ALINE214,ALINE215
      COMPLEX*16 ALINE216,ALINE217,ALINE218,ALINE219,ALINE220
      COMPLEX*16 ALINE221,ALINE222,ALINE223,ALINE224,ALINE225
      COMPLEX*16 ALINE226,ALINE227,ALINE228,ALINE229,ALINE230
      COMPLEX*16 ALINE231,ALINE232,ALINE233,ALINE234,ALINE235
C***********************************************************
      COMPLEX*16 ALLPMOM
      COMPLEX*16 ALLPMOMJ0P, ALLPMOMJ0M
      COMPLEX*16 ALLPMOMJ1
      COMPLEX*16 ALLPMOMJ2P, ALLPMOMJ2M
      COMPLEX*16 ALLPMOMJ0, ALLPMOMJ2
      COMPLEX*16 ALLPMOMJALL
C***********************************************************
      COMPLEX*16 ALINEMOM2,ALINEMOM3,ALINEMOM4
      COMPLEX*16 ALINEMOM5,ALINEMOM6,ALINEMOM7,ALINEMOM8
      COMPLEX*16 ALINEMOM9,ALINEMOM10,ALINEMOM11,ALINEMOM12
      COMPLEX*16 ALINEMOM13,ALINEMOM14,ALINEMOM15,ALINEMOM16
      COMPLEX*16 ALINEMOM17,ALINEMOM18,ALINEMOM19,ALINEMOM20
      COMPLEX*16 ALINEMOM21,ALINEMOM22,ALINEMOM23,ALINEMOM24
      COMPLEX*16 ALINEMOM25,ALINEMOM26,ALINEMOM27,ALINEMOM28
      COMPLEX*16 ALINEMOM29,ALINEMOM30,ALINEMOM31,ALINEMOM32
      COMPLEX*16 ALINEMOM33,ALINEMOM34,ALINEMOM35,ALINEMOM36
      COMPLEX*16 ALINEMOM37,ALINEMOM38,ALINEMOM39,ALINEMOM40
      COMPLEX*16 ALINEMOM41,ALINEMOM42,ALINEMOM43,ALINEMOM44
      COMPLEX*16 ALINEMOM45,ALINEMOM46,ALINEMOM47,ALINEMOM48
      COMPLEX*16 ALINEMOM49,ALINEMOM50,ALINEMOM51,ALINEMOM52
      COMPLEX*16 ALINEMOM53,ALINEMOM54,ALINEMOM55,ALINEMOM56
      COMPLEX*16 ALINEMOM57,ALINEMOM58,ALINEMOM59,ALINEMOM60
      COMPLEX*16 ALINEMOM61,ALINEMOM62,ALINEMOM63,ALINEMOM64
      COMPLEX*16 ALINEMOM65,ALINEMOM66,ALINEMOM67,ALINEMOM68
      COMPLEX*16 ALINEMOM69,ALINEMOM70,ALINEMOM71,ALINEMOM72
      COMPLEX*16 ALINEMOM73,ALINEMOM74,ALINEMOM75,ALINEMOM76
      COMPLEX*16 ALINEMOM77,ALINEMOM78,ALINEMOM79,ALINEMOM80
      COMPLEX*16 ALINEMOM81,ALINEMOM82,ALINEMOM83,ALINEMOM84
      COMPLEX*16 ALINEMOM85,ALINEMOM86,ALINEMOM87,ALINEMOM88
      COMPLEX*16 ALINEMOM89,ALINEMOM90,ALINEMOM91,ALINEMOM92
      COMPLEX*16 ALINEMOM93,ALINEMOM94,ALINEMOM95,ALINEMOM96
      COMPLEX*16 ALINEMOM97,ALINEMOM98,ALINEMOM99,ALINEMOM100
      COMPLEX*16 ALINEMOM101,ALINEMOM102,ALINEMOM103,ALINEMOM104 
      COMPLEX*16 ALINEMOM105,ALINEMOM106,ALINEMOM107,ALINEMOM108
      COMPLEX*16 ALINEMOM109,ALINEMOM110
      COMPLEX*16 ALINEMOM111,ALINEMOM112,ALINEMOM113,ALINEMOM114
      COMPLEX*16 ALINEMOM115,ALINEMOM116,ALINEMOM117,ALINEMOM118
      COMPLEX*16 ALINEMOM119,ALINEMOM120,ALINEMOM121,ALINEMOM122
      COMPLEX*16 ALINEMOM123,ALINEMOM124,ALINEMOM125,ALINEMOM126
      COMPLEX*16 ALINEMOM127,ALINEMOM128,ALINEMOM129,ALINEMOM130
      COMPLEX*16 ALINEMOM131,ALINEMOM132,ALINEMOM133,ALINEMOM134
      COMPLEX*16 ALINEMOM135,ALINEMOM136,ALINEMOM137,ALINEMOM138
      COMPLEX*16 ALINEMOM139,ALINEMOM140,ALINEMOM141,ALINEMOM142
      COMPLEX*16 ALINEMOM143,ALINEMOM144,ALINEMOM145,ALINEMOM146
      COMPLEX*16 ALINEMOM147,ALINEMOM148,ALINEMOM149,ALINEMOM150
      COMPLEX*16 ALINEMOM151,ALINEMOM152,ALINEMOM153,ALINEMOM154
      COMPLEX*16 ALINEMOM155,ALINEMOM156,ALINEMOM157,ALINEMOM158
      COMPLEX*16 ALINEMOM159,ALINEMOM160,ALINEMOM161,ALINEMOM162
      COMPLEX*16 ALINEMOM163,ALINEMOM164,ALINEMOM165,ALINEMOM166
      COMPLEX*16 ALINEMOM167,ALINEMOM168,ALINEMOM169,ALINEMOM170
      COMPLEX*16 ALINEMOM171,ALINEMOM172,ALINEMOM173,ALINEMOM174
      COMPLEX*16 ALINEMOM175,ALINEMOM176,ALINEMOM177,ALINEMOM178
      COMPLEX*16 ALINEMOM179,ALINEMOM180,ALINEMOM181,ALINEMOM182
      COMPLEX*16 ALINEMOM183,ALINEMOM184,ALINEMOM185,ALINEMOM186
      COMPLEX*16 ALINEMOM187,ALINEMOM188,ALINEMOM189,ALINEMOM190
      COMPLEX*16 ALINEMOM191,ALINEMOM192,ALINEMOM193,ALINEMOM194
      COMPLEX*16 ALINEMOM195,ALINEMOM196,ALINEMOM197,ALINEMOM198
      COMPLEX*16 ALINEMOM199,ALINEMOM200,ALINEMOM201,ALINEMOM202      
      COMPLEX*16 ALINEMOM203,ALINEMOM204,ALINEMOM205,ALINEMOM206
      COMPLEX*16 ALINEMOM207,ALINEMOM208,ALINEMOM209,ALINEMOM210
      COMPLEX*16 ALINEMOM211,ALINEMOM212,ALINEMOM213,ALINEMOM214
      COMPLEX*16 ALINEMOM215,ALINEMOM216,ALINEMOM217,ALINEMOM218
      COMPLEX*16 ALINEMOM219,ALINEMOM220,ALINEMOM221,ALINEMOM222
      COMPLEX*16 ALINEMOM223,ALINEMOM224,ALINEMOM225,ALINEMOM226
      COMPLEX*16 ALINEMOM227,ALINEMOM228,ALINEMOM229,ALINEMOM230
      COMPLEX*16 ALINEMOM231,ALINEMOM232,ALINEMOM233,ALINEMOM234
      COMPLEX*16 ALINEMOM235
C***********************************************************
      COMPLEX*16 ACORLP
      COMPLEX*16 AVACLP
C***********************************************************
      COMPLEX*16 ACORLPMOM
      COMPLEX*16 AVACLPMOM
C***********************************************************
      COMPLEX*16 ACORLPMOMJ0P, ACORLPMOMJ0M
      COMPLEX*16 AVACLMOMJ0P, AVACLMOMJ0M
C***********************************************************
      COMPLEX*16 ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
      COMPLEX*16 AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
      COMPLEX*16 ACORLPMOMJALL
      COMPLEX*16 AVACLMOMALL
C***********************************************************
      COMPLEX*16 ACORLPMOMJ2P, ACORLPMOMJ2M
      COMPLEX*16 AVACLMOMJ2P, AVACLMOMJ2M
C***********************************************************
      COMPLEX*16 ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
      COMPLEX*16 AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
C***********************************************************
      COMPLEX*16 ACORLJ1P,ACORLJ1M
      COMPLEX*16 AVACLJ1P,AVACLJ1M
C***********************************************************
      COMPLEX*16 ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
      COMPLEX*16 AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
C***********************************************************
      COMPLEX*16 ACORLJ0, ACORLJ1, ACORLJ2, ACORLALL
      COMPLEX*16 AVACLJ0, AVACLJ1, AVACLJ2, AVACLALL
C***********************************************************      
      COMPLEX*16 ALINEP(LX4,NTOTAL)
      COMPLEX*16 ALINEPMOM(2,LX4,NTOTALMOM)
C***********************************************************
C     NEW CORRELATORS
C***********************************************************
      COMPLEX*16 ALINEJ0PP(LX4,NOPJ0PP)
      COMPLEX*16 ALINEJ0PM(LX4,NOPJ0PM)
      COMPLEX*16 ALINEJ0MP(LX4,NOPJ0MP)
      COMPLEX*16 ALINEJ0MM(LX4,NOPJ0MM)
      COMPLEX*16 ALINEJ2PP(LX4,NOPJ2PP)
      COMPLEX*16 ALINEJ2PM(LX4,NOPJ2PM)
      COMPLEX*16 ALINEJ2MP(LX4,NOPJ2MP)
      COMPLEX*16 ALINEJ2MM(LX4,NOPJ2MM)
      COMPLEX*16 ALINEJ1P(LX4,NOPJ1P)
      COMPLEX*16 ALINEJ1M(LX4,NOPJ1M)
C***********************************************************
      COMPLEX*16 ALINEJ0(LX4,NOPFULJ0)
      COMPLEX*16 ALINEJ2(LX4,NOPFULJ2)
      COMPLEX*16 ALINEJ1(LX4,NOPFULJ1)
      COMPLEX*16 ALINEJALL(LX4,NTOTAL)
C***********************************************************
      COMPLEX*16 ALINEMOMJ0P(2,LX4,NOPJ0PMOM)
      COMPLEX*16 ALINEMOMJ0M(2,LX4,NOPJ0MMOM)
      COMPLEX*16 ALINEMOMJ2P(2,LX4,NOPJ2PMOM)
      COMPLEX*16 ALINEMOMJ2M(2,LX4,NOPJ2MMOM)
      COMPLEX*16 ALINEMOMJ1(2,LX4,NOPJ1MOM)
C***********************************************************
      COMPLEX*16 ALINEMOMJ0(2,LX4,NOPFULJ0MOM)
      COMPLEX*16 ALINEMOMJ2(2,LX4,NOPFULJ2MOM)
C***********************************************************
      COMPLEX*16 ALINEMOMJALL(2,LX4,NTOTALMOM)      
C************************************************************
      JBIN=(ITR-1)/IBIN+1
C************************************************************
      DO 2 N4=1,LX4
C
         DO 3 IDL=1,IBLOK
            ALINEP(N4,IDL)=ALINE1(N4,IDL)
            ALINEJ0PP(N4,IDL)=ALINE1(N4,IDL)
            ALINEJ0PM(N4,IDL)=ALINE87(N4,IDL)
            ALINEJ0MP(N4,IDL)=ALINE63(N4,IDL)            
            ALINEJ0MM(N4,IDL)=ALINE27(N4,IDL)            
c
            ALINEJ1P(N4,IDL)=ALINE3(N4,IDL)
            ALINEJ1M(N4,IDL)=ALINE9(N4,IDL)
c
            ALINEJ2PP(N4,IDL)=ALINE4(N4,IDL)
            ALINEJ2PM(N4,IDL)=ALINE30(N4,IDL)
            ALINEJ2MP(N4,IDL)=ALINE31(N4,IDL)            
            ALINEJ2MM(N4,IDL)=ALINE97(N4,IDL)            
c            
 3       CONTINUE
C
         ID=IBLOK
         DO 4 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE2(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE2(N4,IDL)
            ALINEJ0PM(N4,ID)=ALINE147(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE81(N4,IDL)                        
            ALINEJ0MM(N4,ID)=ALINE33(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE6(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE12(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE7(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE36(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE37(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE155(N4,IDL)            
c            
 4       CONTINUE
C
         ID=2*IBLOK
         DO 5 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE3(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE5(N4,IDL)
            ALINEJ0PM(N4,ID)=ALINE157(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE88(N4,IDL)                        
            ALINEJ0MM(N4,ID)=ALINE39(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE18(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE15(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE10(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE42(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE43(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE165(N4,IDL)            
c
 5       CONTINUE
C
         ID=3*IBLOK
         DO 6 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE4(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE8(N4,IDL)            
            ALINEJ0PM(N4,ID)=ALINE167(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE123(N4,IDL)            
            ALINEJ0MM(N4,ID)=ALINE45(N4,IDL)            
c
            ALINEJ1P(N4,ID)=ALINE21(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE24(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE13(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE48(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE49(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE175(N4,IDL)            
c            
 6       CONTINUE
C
         ID=4*IBLOK
         DO 7 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE5(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE11(N4,IDL)
            ALINEJ0PM(N4,ID)=ALINE177(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE129(N4,IDL)            
            ALINEJ0MM(N4,ID)=ALINE57(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE28(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE29(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE16(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE60(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE61(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE185(N4,IDL)            
c
 7       CONTINUE
C
         ID=5*IBLOK
         DO 8 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE6(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE14(N4,IDL)
            ALINEJ0PM(N4,ID)=ALINE187(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE135(N4,IDL)            
            ALINEJ0MM(N4,ID)=ALINE69(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE34(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE35(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE19(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE95(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE67(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE195(N4,IDL)            
c
 8       CONTINUE
C
         ID=6*IBLOK
         DO 9 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE7(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE17(N4,IDL)
            ALINEJ0PM(N4,ID)=ALINE197(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE148(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE75(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE40(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE41(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE22(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE102(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE85(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE205(N4,IDL)            
c
 9       CONTINUE
C
         ID=7*IBLOK
         DO 10 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE8(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE20(N4,IDL)
            ALINEJ0PM(N4,ID)=ALINE207(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE158(N4,IDL)            
            ALINEJ0MM(N4,ID)=ALINE89(N4,IDL)                        
c
            ALINEJ1P(N4,ID)=ALINE46(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE47(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE25(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE108(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE96(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE215(N4,IDL)            
c
 10      CONTINUE
C     
         ID=8*IBLOK
         DO 11 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE9(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE23(N4,IDL)
            ALINEJ0PM(N4,ID)=ALINE217(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE168(N4,IDL)            
            ALINEJ0MM(N4,ID)=ALINE99(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE51(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE54(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE52(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE114(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE103(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE225(N4,IDL)            
c
 11      CONTINUE
C     
         ID=9*IBLOK
         DO 12 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE10(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE26(N4,IDL)            
            ALINEJ0PM(N4,ID)=ALINE227(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE178(N4,IDL)            
            ALINEJ0MM(N4,ID)=ALINE105(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE58(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE59(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE55(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE120(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE109(N4,IDL)            
            ALINEJ2MM(N4,ID)=ALINE235(N4,IDL)            
c
 12      CONTINUE
C     
         ID=10*IBLOK
         DO 13 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE11(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE32(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE188(N4,IDL)            
            ALINEJ0MM(N4,ID)=ALINE111(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE64(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE71(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE66(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE144(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE115(N4,IDL)
            ALINEJ2MM(N4,ID)=ALINE73(N4,IDL)
c
 13      CONTINUE
C
         ID=11*IBLOK
         DO 14 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE12(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE38(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE198(N4,IDL) 
            ALINEJ0MM(N4,ID)=ALINE117(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE65(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE77(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE72(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE153(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE121(N4,IDL)
            ALINEJ2MM(N4,ID)=ALINE79(N4,IDL)
c
 14      CONTINUE
C
         ID=12*IBLOK
         DO 15 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE13(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE44(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE208(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE141(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE70(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE82(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE78(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE163(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE127(N4,IDL)            
c
 15      CONTINUE
C
         ID=13*IBLOK
         DO 16 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE14(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE50(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE218(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE149(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE76(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE83(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE84(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE173(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE133(N4,IDL)            
c
 16      CONTINUE
C
         ID=14*IBLOK
         DO 17 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE15(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE53(N4,IDL)
            ALINEJ0MP(N4,ID)=ALINE228(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE159(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE90(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE92(N4,IDL)            
c
            ALINEJ2PP(N4,ID)=ALINE94(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE183(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE139(N4,IDL)            
c
 17      CONTINUE
C
         ID=15*IBLOK
         DO 18 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE16(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE56(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE169(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE91(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE93(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE126(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE193(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE145(N4,IDL)            
c
 18      CONTINUE
C
         ID=16*IBLOK
         DO 19 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE17(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE62(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE179(N4,IDL)            
c
            ALINEJ1P(N4,ID)=ALINE100(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE101(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE132(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE203(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE154(N4,IDL)            
c
 19      CONTINUE         
C
         ID=17*IBLOK
         DO 20 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE18(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE68(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE189(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE106(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE107(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE138(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE213(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE164(N4,IDL)            
c
 20      CONTINUE
C
         ID=18*IBLOK
         DO 21 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE19(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE74(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE199(N4,IDL)            
c
            ALINEJ1P(N4,ID)=ALINE112(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE113(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE152(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE223(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE174(N4,IDL)            
c  
 21      CONTINUE
C
         ID=19*IBLOK
         DO 22 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE20(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE80(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE209(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE118(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE136(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE162(N4,IDL)
            ALINEJ2PM(N4,ID)=ALINE233(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE184(N4,IDL)            
c
 22      CONTINUE
C
         ID=20*IBLOK
         DO 23 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE21(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE86(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE219(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE124(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE137(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE172(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE194(N4,IDL)            
c            
 23      CONTINUE
C
         ID=21*IBLOK
         DO 24 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE22(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE98(N4,IDL)
            ALINEJ0MM(N4,ID)=ALINE229(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE125(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE143(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE182(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE204(N4,IDL)            
c
 24      CONTINUE
C
         ID=22*IBLOK
         DO 25 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE23(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE104(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE130(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE151(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE192(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE214(N4,IDL)            
c
 25      CONTINUE
C
         ID=23*IBLOK
         DO 26 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE24(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE110(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE131(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE161(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE202(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE224(N4,IDL)            
c
 26      CONTINUE
C
         ID=24*IBLOK
         DO 27 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE25(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE116(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE142(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE171(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE212(N4,IDL)
            ALINEJ2MP(N4,ID)=ALINE234(N4,IDL)            
c
 27      CONTINUE
C
         ID=25*IBLOK
         DO 28 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE26(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE122(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE150(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE181(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE222(N4,IDL)
c
 28      CONTINUE
C
         ID=26*IBLOK
         DO 29 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE27(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE128(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE160(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE191(N4,IDL)
c
            ALINEJ2PP(N4,ID)=ALINE232(N4,IDL)
c
 29      CONTINUE
C
         ID=27*IBLOK
         DO 31 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE28(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE134(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE170(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE201(N4,IDL)
c
 31      CONTINUE
C     
         ID=28*IBLOK
         DO 32 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE29(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE140(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE180(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE211(N4,IDL)
c
 32      CONTINUE
C
         ID=29*IBLOK
         DO 33 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE30(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE146(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE190(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE221(N4,IDL)
c
 33      CONTINUE
C
         ID=30*IBLOK
         DO 34 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE31(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE156(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE200(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE231(N4,IDL)            
c
 34      CONTINUE
C
         ID=31*IBLOK
         DO 35 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE32(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE166(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE210(N4,IDL)
            ALINEJ1M(N4,ID)=ALINE119(N4,IDL)
c     
 35      CONTINUE
C
         ID=32*IBLOK
         DO 36 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE33(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE176(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE220(N4,IDL)
c
 36      CONTINUE
C
         ID=33*IBLOK
         DO 37 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE34(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE186(N4,IDL)
c
            ALINEJ1P(N4,ID)=ALINE230(N4,IDL)
c
 37      CONTINUE
C     
         ID=34*IBLOK
         DO 38 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE35(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE196(N4,IDL)
c
c
 38      CONTINUE
C
         ID=35*IBLOK
         DO 39 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE36(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE206(N4,IDL)            
 39      CONTINUE
C     
         ID=36*IBLOK
         DO 40 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE37(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE216(N4,IDL)
 40      CONTINUE
C
         ID=37*IBLOK
         DO 41 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE38(N4,IDL)
            ALINEJ0PP(N4,ID)=ALINE226(N4,IDL)                        
 41      CONTINUE
C
         ID=38*IBLOK
         DO 42 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE39(N4,IDL)
 42      CONTINUE
C
         ID=39*IBLOK
         DO 43 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE40(N4,IDL)
 43      CONTINUE
C
         ID=40*IBLOK
         DO 44 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE41(N4,IDL)
 44      CONTINUE
C
         ID=41*IBLOK
         DO 45 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE42(N4,IDL)
 45      CONTINUE
C
         ID=42*IBLOK
         DO 46 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE43(N4,IDL)
 46      CONTINUE
C
         ID=43*IBLOK
         DO 47 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE44(N4,IDL)
 47      CONTINUE
C
         ID=44*IBLOK
         DO 48 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE45(N4,IDL)
 48      CONTINUE
C     
         ID=45*IBLOK
         DO 49 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE46(N4,IDL)
 49      CONTINUE
C
         ID=46*IBLOK
         DO 50 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE47(N4,IDL)
 50      CONTINUE
C
         ID=47*IBLOK
         DO 51 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE48(N4,IDL)
 51      CONTINUE
C
         ID=48*IBLOK
         DO 52 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE49(N4,IDL)
 52      CONTINUE
C
         ID=49*IBLOK
         DO 53 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE50(N4,IDL)
 53      CONTINUE
C
         ID=50*IBLOK
         DO 54 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE51(N4,IDL)
 54      CONTINUE
C
         ID=51*IBLOK
         DO 55 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE52(N4,IDL)
 55      CONTINUE
C
         ID=52*IBLOK
         DO 56 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE53(N4,IDL)
 56      CONTINUE
C
         ID=53*IBLOK
         DO 57 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE54(N4,IDL)
 57      CONTINUE
C
         ID=54*IBLOK
         DO 58 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE55(N4,IDL)
 58      CONTINUE
C
         ID=55*IBLOK
         DO 59 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE56(N4,IDL)
 59      CONTINUE
C
         ID=56*IBLOK
         DO 60 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE57(N4,IDL)
 60      CONTINUE
C
         ID=57*IBLOK
         DO 61 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE58(N4,IDL)
 61      CONTINUE
C
         ID=58*IBLOK
         DO 62 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE59(N4,IDL)
 62      CONTINUE
C
         ID=59*IBLOK
         DO 63 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE60(N4,IDL)
 63      CONTINUE
C
         ID=60*IBLOK
         DO 64 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE61(N4,IDL)
 64      CONTINUE
C
         ID=61*IBLOK
         DO 65 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE62(N4,IDL)
 65      CONTINUE
C
         ID=62*IBLOK
         DO 66 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE63(N4,IDL)
 66      CONTINUE
C
         ID=63*IBLOK
         DO 67 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE64(N4,IDL)
 67      CONTINUE
C
         ID=64*IBLOK
         DO 68 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE65(N4,IDL)
 68      CONTINUE
C
         ID=65*IBLOK
         DO 69 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE66(N4,IDL)
 69      CONTINUE
C
         ID=66*IBLOK
         DO 70 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE67(N4,IDL)
 70      CONTINUE
C
         ID=67*IBLOK
         DO 71 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE68(N4,IDL)
 71      CONTINUE
C
         ID=68*IBLOK
         DO 72 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE69(N4,IDL)
 72      CONTINUE
C
         ID=69*IBLOK
         DO 73 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE70(N4,IDL)
 73      CONTINUE
C
         ID=70*IBLOK
         DO 74 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE71(N4,IDL)
 74      CONTINUE
C
         ID=71*IBLOK
         DO 75 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE72(N4,IDL)
 75      CONTINUE
C
         ID=72*IBLOK
         DO 76 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE73(N4,IDL)
 76      CONTINUE
C
         ID=73*IBLOK
         DO 77 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE74(N4,IDL)
 77      CONTINUE
C
         ID=74*IBLOK
         DO 78 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE75(N4,IDL)
 78      CONTINUE
C
         ID=75*IBLOK
         DO 79 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE76(N4,IDL)
 79      CONTINUE
C
         ID=76*IBLOK
         DO 80 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE77(N4,IDL)
 80      CONTINUE
C
         ID=77*IBLOK
         DO 81 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE78(N4,IDL)
 81      CONTINUE
C     
         ID=78*IBLOK
         DO 82 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE79(N4,IDL)
 82      CONTINUE
C
         ID=79*IBLOK
         DO 83 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE80(N4,IDL)
 83      CONTINUE
C
         ID=80*IBLOK
         DO 84 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE81(N4,IDL)
 84      CONTINUE
C
         ID=81*IBLOK
         DO 85 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE82(N4,IDL)
 85      CONTINUE
C
         ID=82*IBLOK
         DO 86 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE83(N4,IDL)
 86      CONTINUE
C
         ID=83*IBLOK
         DO 87 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE84(N4,IDL)
 87      CONTINUE
C
         ID=84*IBLOK
         DO 88 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE85(N4,IDL)
 88      CONTINUE
C
         ID=85*IBLOK
         DO 89 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE86(N4,IDL)
 89      CONTINUE
C
         ID=86*IBLOK
         DO 90 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE87(N4,IDL)
 90      CONTINUE
C
         ID=87*IBLOK
         DO 91 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE88(N4,IDL)
 91      CONTINUE
C
         ID=88*IBLOK
         DO 92 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE89(N4,IDL)
 92      CONTINUE
C
         ID=89*IBLOK
         DO 93 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE90(N4,IDL)
 93      CONTINUE
C
         ID=90*IBLOK
         DO 94 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE91(N4,IDL)
 94      CONTINUE
C
         ID=91*IBLOK
         DO 95 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE92(N4,IDL)
 95      CONTINUE
C
         ID=92*IBLOK
         DO 96 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE93(N4,IDL)
 96      CONTINUE
C
         ID=93*IBLOK
         DO  97 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE94(N4,IDL)
 97      CONTINUE
C
         ID=94*IBLOK
         DO 98 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE95(N4,IDL)
 98      CONTINUE
C
         ID=95*IBLOK
         DO 99 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE96(N4,IDL)
 99      CONTINUE
C
         ID=96*IBLOK
         DO 300 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE97(N4,IDL)
 300     CONTINUE
C
         ID=97*IBLOK
         DO 301 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE98(N4,IDL)
 301     CONTINUE
C
         ID=98*IBLOK
         DO 302 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE99(N4,IDL)
 302     CONTINUE
C
         ID=99*IBLOK
         DO 303 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE100(N4,IDL)
 303     CONTINUE
C
         ID=100*IBLOK
         DO 304 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE101(N4,IDL)
 304     CONTINUE
C
         ID=101*IBLOK
         DO 305 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE102(N4,IDL)
 305     CONTINUE
C
         ID=102*IBLOK
         DO 306 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE103(N4,IDL)
 306     CONTINUE
C
         ID=103*IBLOK
         DO 307 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE104(N4,IDL)
 307     CONTINUE
C
         ID=104*IBLOK
         DO 308 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE105(N4,IDL)
 308     CONTINUE
C
         ID=105*IBLOK
         DO 309 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE106(N4,IDL)
 309     CONTINUE
C
         ID=106*IBLOK
         DO 310 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE107(N4,IDL)
 310     CONTINUE         
C
         ID=107*IBLOK
         DO 311 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE108(N4,IDL)
 311     CONTINUE
C
         ID=108*IBLOK
         DO 312 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE109(N4,IDL)
 312     CONTINUE
C
         ID=109*IBLOK
         DO 313 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE110(N4,IDL)
 313     CONTINUE
C
         ID=110*IBLOK
         DO 314 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE111(N4,IDL)
 314     CONTINUE
C
         ID=111*IBLOK
         DO 315 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE112(N4,IDL)
 315     CONTINUE
C
         ID=112*IBLOK
         DO 316 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE113(N4,IDL)
 316     CONTINUE
C
         ID=113*IBLOK
         DO 317 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE114(N4,IDL)
 317     CONTINUE
C
         ID=114*IBLOK
         DO 318 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE115(N4,IDL)
 318     CONTINUE
C
         ID=115*IBLOK
         DO 319 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE116(N4,IDL)
 319     CONTINUE
C
         ID=116*IBLOK
         DO 320 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE117(N4,IDL)
 320     CONTINUE
C
         ID=117*IBLOK
         DO 321 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE118(N4,IDL)
 321     CONTINUE
C     
         ID=118*IBLOK
         DO 322 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE119(N4,IDL)
 322     CONTINUE
C
         ID=119*IBLOK
         DO 323 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE120(N4,IDL)
 323      CONTINUE
C
         ID=120*IBLOK
         DO 324 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE121(N4,IDL)
 324      CONTINUE
C
         ID=121*IBLOK
         DO 325 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE122(N4,IDL)
 325      CONTINUE
C
         ID=122*IBLOK
         DO 326 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE123(N4,IDL)
 326      CONTINUE
C
         ID=123*IBLOK
         DO 327 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE124(N4,IDL)
 327      CONTINUE
C     
         ID=124*IBLOK
         DO 328 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE125(N4,IDL)
 328      CONTINUE
C
         ID=125*IBLOK
         DO 329 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE126(N4,IDL)
 329      CONTINUE
C     
         ID=126*IBLOK
         DO 330 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE127(N4,IDL)
 330      CONTINUE
C
         ID=127*IBLOK
         DO 331 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE128(N4,IDL)
 331      CONTINUE
C
         ID=128*IBLOK
         DO 332 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE129(N4,IDL)
 332      CONTINUE
C
         ID=129*IBLOK
         DO 333 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE130(N4,IDL)
 333      CONTINUE
C
         ID=130*IBLOK
         DO 334 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE131(N4,IDL)
 334      CONTINUE
C
         ID=131*IBLOK
         DO 335 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE132(N4,IDL)
 335      CONTINUE
C
         ID=132*IBLOK
         DO 336 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE133(N4,IDL)
 336      CONTINUE
C
         ID=133*IBLOK
         DO 337 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE134(N4,IDL)
 337      CONTINUE
C
         ID=134*IBLOK
         DO 338 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE135(N4,IDL)
 338      CONTINUE
C     
         ID=135*IBLOK
         DO 339 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE136(N4,IDL)
 339      CONTINUE
C
          ID=136*IBLOK
         DO 340 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE137(N4,IDL)
 340      CONTINUE
C
         ID=137*IBLOK
         DO 341 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE138(N4,IDL)
 341      CONTINUE
C
         ID=138*IBLOK
         DO 342 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE139(N4,IDL)
 342      CONTINUE
C
         ID=139*IBLOK
         DO 343 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE140(N4,IDL)
 343      CONTINUE
C
         ID=140*IBLOK
         DO 344 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE141(N4,IDL)
 344      CONTINUE
C
         ID=141*IBLOK
         DO 345 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE142(N4,IDL)
 345      CONTINUE
C
         ID=142*IBLOK
         DO 346 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE143(N4,IDL)
 346      CONTINUE
C
         ID=143*IBLOK
         DO 347 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE144(N4,IDL)
 347      CONTINUE
C
         ID=144*IBLOK
         DO 348 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE145(N4,IDL)
 348      CONTINUE
C
         ID=145*IBLOK
         DO 349 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE146(N4,IDL)
 349     CONTINUE
C
         ID=146*IBLOK
         DO 350 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE147(N4,IDL)
 350     CONTINUE
C
         ID=147*IBLOK
         DO 351 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE148(N4,IDL)
 351     CONTINUE
C
         ID=148*IBLOK
         DO 352 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE149(N4,IDL)
 352     CONTINUE
C
         ID=149*IBLOK
         DO 353 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE150(N4,IDL)
 353     CONTINUE
C
         ID=150*IBLOK
         DO 354 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE151(N4,IDL)
 354     CONTINUE
C
         ID=151*IBLOK
         DO 355 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE152(N4,IDL)
 355     CONTINUE
C
         ID=152*IBLOK
         DO 356 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE153(N4,IDL)
 356      CONTINUE
C
         ID=153*IBLOK
         DO 357 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE154(N4,IDL)
 357      CONTINUE
C
         ID=154*IBLOK
         DO 358 IDL=1, IBLOK
            ID=ID+1
           ALINEP(N4,ID)=ALINE155(N4,IDL)
 358      CONTINUE
C
         ID=155*IBLOK
         DO 359 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE156(N4,IDL)
 359      CONTINUE
C
         ID=156*IBLOK
         DO 360 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE157(N4,IDL)
 360      CONTINUE
C
         ID=157*IBLOK
         DO 361 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE158(N4,IDL)
 361      CONTINUE
C
         ID=158*IBLOK
         DO 362 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE159(N4,IDL)
 362      CONTINUE
C
         ID=159*IBLOK
         DO 363 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE160(N4,IDL)
 363      CONTINUE
C
         ID=160*IBLOK
         DO 364 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE161(N4,IDL)
 364      CONTINUE
C
         ID=161*IBLOK
         DO 365 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE162(N4,IDL)
 365      CONTINUE
C
         ID=162*IBLOK
         DO 366 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE163(N4,IDL)
 366      CONTINUE
C
         ID=163*IBLOK
         DO 367 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE164(N4,IDL)
 367      CONTINUE
C
         ID=164*IBLOK
         DO 368 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE165(N4,IDL)
 368      CONTINUE
C
         ID=165*IBLOK
         DO 369 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE166(N4,IDL)
 369      CONTINUE
C
         ID=166*IBLOK
         DO 370 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE167(N4,IDL)
 370      CONTINUE
C
         ID=167*IBLOK
         DO 371 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE168(N4,IDL)
 371      CONTINUE
C
         ID=168*IBLOK
         DO 372 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE169(N4,IDL)
 372      CONTINUE
C
         ID=169*IBLOK
         DO 373 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE170(N4,IDL)
 373      CONTINUE
C
         ID=170*IBLOK
         DO 374 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE171(N4,IDL)
 374      CONTINUE
C
         ID=171*IBLOK
         DO 375 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE172(N4,IDL)
 375      CONTINUE
C
         ID=172*IBLOK
         DO 376 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE173(N4,IDL)
 376      CONTINUE
C
         ID=173*IBLOK
         DO 377 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE174(N4,IDL)
 377      CONTINUE
C
         ID=174*IBLOK
         DO 378 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE175(N4,IDL)
 378      CONTINUE
C
         ID=175*IBLOK
         DO 379 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE176(N4,IDL)
 379      CONTINUE
C
         ID=176*IBLOK
         DO 380 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE177(N4,IDL)
 380      CONTINUE
C
         ID=177*IBLOK
         DO 381 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE178(N4,IDL)
 381      CONTINUE
C     
         ID=178*IBLOK
         DO 382 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE179(N4,IDL)
 382      CONTINUE
C
         ID=179*IBLOK
         DO 383 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE180(N4,IDL)
 383      CONTINUE
C
         ID=180*IBLOK
         DO 384 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE181(N4,IDL)
 384      CONTINUE
C
         ID=181*IBLOK
         DO 385 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE182(N4,IDL)
 385      CONTINUE
C
         ID=182*IBLOK
         DO 386 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE183(N4,IDL)
 386      CONTINUE
C
         ID=183*IBLOK
         DO 387 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE184(N4,IDL)
 387      CONTINUE
C
         ID=184*IBLOK
         DO 388 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE185(N4,IDL)
 388      CONTINUE
C
         ID=185*IBLOK
         DO 389 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE186(N4,IDL)
 389      CONTINUE
C
         ID=186*IBLOK
         DO 390 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE187(N4,IDL)
 390      CONTINUE
C
         ID=187*IBLOK
         DO 391 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE188(N4,IDL)
 391      CONTINUE
C
         ID=188*IBLOK
         DO 392 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE189(N4,IDL)
 392      CONTINUE
C
         ID=189*IBLOK
         DO 393 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE190(N4,IDL)
 393      CONTINUE
C
         ID=190*IBLOK
         DO 394 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE191(N4,IDL)
 394      CONTINUE
C
         ID=191*IBLOK
         DO 395 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE192(N4,IDL)
 395      CONTINUE
C
         ID=192*IBLOK
         DO 396 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE193(N4,IDL)
 396      CONTINUE
C
         ID=193*IBLOK
         DO 397 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE194(N4,IDL)
 397      CONTINUE
C
         ID=194*IBLOK
         DO 398 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE195(N4,IDL)
 398      CONTINUE
C
         ID=195*IBLOK
         DO 399 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE196(N4,IDL)
 399      CONTINUE
C
         ID=196*IBLOK
         DO 400 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE197(N4,IDL)
 400     CONTINUE
C
         ID=197*IBLOK
         DO 401 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE198(N4,IDL)
 401     CONTINUE
C
         ID=198*IBLOK
         DO 402 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE199(N4,IDL)
 402     CONTINUE
C
         ID=199*IBLOK
         DO 403 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE200(N4,IDL)
 403     CONTINUE
C
         ID=200*IBLOK
         DO 404 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE201(N4,IDL)
 404     CONTINUE
C
         ID=201*IBLOK
         DO 405 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE202(N4,IDL)
 405     CONTINUE
C
         ID=202*IBLOK
         DO 406 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE203(N4,IDL)
 406     CONTINUE
C
         ID=203*IBLOK
         DO 407 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE204(N4,IDL)
 407     CONTINUE
C
         ID=204*IBLOK
         DO 408 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE205(N4,IDL)
 408     CONTINUE
C
         ID=205*IBLOK
         DO 409 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE206(N4,IDL)
 409     CONTINUE
C
         ID=206*IBLOK
         DO 410 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE207(N4,IDL)
 410     CONTINUE         
C
         ID=207*IBLOK
         DO 411 IDL=1,IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE208(N4,IDL)
 411     CONTINUE
C
         ID=208*IBLOK
         DO 412 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE209(N4,IDL)
 412     CONTINUE
C
         ID=209*IBLOK
         DO 413 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE210(N4,IDL)
 413     CONTINUE
C
         ID=210*IBLOK
         DO 414 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE211(N4,IDL)
 414     CONTINUE
C
         ID=211*IBLOK
         DO 415 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE212(N4,IDL)
 415     CONTINUE
C
         ID=212*IBLOK
         DO 416 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE213(N4,IDL)
 416     CONTINUE
C
         ID=213*IBLOK
         DO 417 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE214(N4,IDL)
 417     CONTINUE
C
         ID=214*IBLOK
         DO 418 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE215(N4,IDL)
 418     CONTINUE
C
         ID=215*IBLOK
         DO 419 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE216(N4,IDL)
 419     CONTINUE
C
         ID=216*IBLOK
         DO 420 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE217(N4,IDL)
 420     CONTINUE
C
         ID=217*IBLOK
         DO 421 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE218(N4,IDL)
 421     CONTINUE
C     
         ID=218*IBLOK
         DO 422 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE219(N4,IDL)
 422     CONTINUE
C
         ID=219*IBLOK
         DO 423 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE220(N4,IDL)
 423      CONTINUE
C
         ID=220*IBLOK
         DO 424 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE221(N4,IDL)
 424      CONTINUE
C
         ID=221*IBLOK
         DO 425 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE222(N4,IDL)
 425      CONTINUE
C
         ID=222*IBLOK
         DO 426 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE223(N4,IDL)
 426      CONTINUE
C
         ID=223*IBLOK
         DO 427 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE224(N4,IDL)
 427      CONTINUE
C     
         ID=224*IBLOK
         DO 428 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE225(N4,IDL)
 428      CONTINUE
C
         ID=225*IBLOK
         DO 429 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE226(N4,IDL)
 429      CONTINUE
C     
         ID=226*IBLOK
         DO 430 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE227(N4,IDL)
 430      CONTINUE
C
         ID=227*IBLOK
         DO 431 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE228(N4,IDL)
 431      CONTINUE
C
         ID=228*IBLOK
         DO 432 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE229(N4,IDL)
 432      CONTINUE
C
         ID=229*IBLOK
         DO 433 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE230(N4,IDL)
 433      CONTINUE
C
         ID=230*IBLOK
         DO 434 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE231(N4,IDL)
 434      CONTINUE
C
         ID=231*IBLOK
         DO 435 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE232(N4,IDL)
 435      CONTINUE
C
         ID=232*IBLOK
         DO 436 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE233(N4,IDL)
 436      CONTINUE
C
         ID=233*IBLOK
         DO 437 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE234(N4,IDL)
 437      CONTINUE
C
         ID=234*IBLOK
         DO 438 IDL=1, IBLOK
            ID=ID+1
            ALINEP(N4,ID)=ALINE235(N4,IDL)
 438      CONTINUE 
C          
 2    CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC FULL J LINES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO N4=1, LX4
C         
         DO IDL=1, NOPJ0PP
            ALINEJ0(N4,IDL)=ALINEJ0PP(N4,IDL)
         ENDDO
C
         ID=NOPJ0PP
         DO IDL=1, NOPJ0PM
            ID=ID+1
            ALINEJ0(N4,ID)=ALINEJ0PM(N4,IDL)
         ENDDO
C
         ID=NOPJ0PP+NOPJ0PM
         DO IDL=1, NOPJ0MP
            ID=ID+1
            ALINEJ0(N4,ID)=ALINEJ0MP(N4,IDL)
         ENDDO         
C
         ID=NOPJ0PP+NOPJ0PM+NOPJ0MP
         DO IDL=1, NOPJ0MM
            ID=ID+1
            ALINEJ0(N4,ID)=ALINEJ0MM(N4,IDL)
         ENDDO         
C
C     *************************************************************
C
         DO IDL=1, NOPJ1P
            ALINEJ1(N4,IDL)=ALINEJ1P(N4,IDL)
         ENDDO
C
         ID=NOPJ1P
         DO IDL=1, NOPJ1M
            ID=ID+1
            ALINEJ1(N4,ID)=ALINEJ1M(N4,IDL)
         ENDDO
C
C     *************************************************************
C
C         
         DO IDL=1, NOPJ2PP
            ALINEJ2(N4,IDL)=ALINEJ2PP(N4,IDL)
         ENDDO
C
         ID=NOPJ2PP
         DO IDL=1, NOPJ2PM
            ID=ID+1
            ALINEJ2(N4,ID)=ALINEJ2PM(N4,IDL)
         ENDDO
C
         ID=NOPJ2PP+NOPJ2PM
         DO IDL=1, NOPJ2MP
            ID=ID+1
            ALINEJ2(N4,ID)=ALINEJ2MP(N4,IDL)
         ENDDO         
C
         ID=NOPJ2PP+NOPJ2PM+NOPJ2MP
         DO IDL=1, NOPJ2MM
            ID=ID+1
            ALINEJ2(N4,ID)=ALINEJ2MM(N4,IDL)
         ENDDO                  
C     
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC FULL J LINES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO N4=1, LX4
C
         DO IDL=1, NOPFULJ0
            ALINEJALL(N4,IDL)=ALINEJ0(N4,IDL)
         ENDDO
C
         ID=NOPFULJ0
         DO IDL=1, NOPFULJ1
            ID=ID+1
            ALINEJALL(N4,ID)=ALINEJ1(N4,IDL)
         ENDDO
C
         ID=NOPFULJ0+NOPFULJ1
         DO IDL=1, NOPFULJ2
            ID=ID+1
            ALINEJALL(N4,ID)=ALINEJ2(N4,IDL)
         ENDDO
C         
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 102 N4=1,LX4
C     
         DO 103 IDL=1,IBLOK
            DO IJ=1,2
               ALINEPMOM(IJ,N4,IDL)=ALINEMOM2(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,IDL)=ALINEMOM2(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,IDL)=ALINEMOM27(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,IDL)=ALINEMOM3(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,IDL)=ALINEMOM4(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,IDL)=ALINEMOM31(N4,IDL,IJ)               
c     
            ENDDO
 103     CONTINUE
C
         ID=IBLOK
         DO 104 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM3(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM5(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM33(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM6(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM7(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM37(N4,IDL,IJ)
c 
            ENDDO
 104     CONTINUE
C
         ID=2*IBLOK
         DO 105 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM4(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM8(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM39(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM9(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM10(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM43(N4,IDL,IJ)
c 
            ENDDO
 105     CONTINUE
C     
         ID=3*IBLOK
         DO 106 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM5(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM11(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM45(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM12(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM13(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM49(N4,IDL,IJ)
c 
            ENDDO
 106     CONTINUE
C
         ID=4*IBLOK
         DO 107 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM6(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM14(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM57(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM15(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM16(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM61(N4,IDL,IJ)
c 
            ENDDO
 107     CONTINUE
C
         ID=5*IBLOK
         DO 108 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM7(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM17(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM63(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM18(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM19(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM67(N4,IDL,IJ)
c 
            ENDDO
 108     CONTINUE
C
         ID=6*IBLOK
         DO 109 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM8(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM20(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM69(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM21(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM22(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM73(N4,IDL,IJ)
c 
            ENDDO
 109     CONTINUE
C
         ID=7*IBLOK
         DO 110 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM9(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM23(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM75(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM24(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM25(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM79(N4,IDL,IJ)
c 
            ENDDO
 110     CONTINUE
C
         ID=8*IBLOK
         DO 111 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM10(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM26(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM81(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM28(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM30(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM85(N4,IDL,IJ)
c               
            ENDDO
 111     CONTINUE
C
         ID=9*IBLOK
         DO 112 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM11(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM32(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM88(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM29(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM36(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM96(N4,IDL,IJ)
c 
            ENDDO
 112     CONTINUE
C
         ID=10*IBLOK
         DO 113 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM12(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM38(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM89(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM34(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM42(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM97(N4,IDL,IJ)
c
            ENDDO
 113     CONTINUE
C
         ID=11*IBLOK
         DO 114 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM13(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM44(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM105(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM35(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM48(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM103(N4,IDL,IJ)
c
            ENDDO
 114     CONTINUE
C
         ID=12*IBLOK
         DO 115 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM14(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM50(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM111(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM40(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM52(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM109(N4,IDL,IJ)
c
            ENDDO
 115     CONTINUE
C
         ID=13*IBLOK
         DO 116 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM15(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM53(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM117(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM41(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM55(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM115(N4,IDL,IJ)
c
            ENDDO
 116     CONTINUE
c            
         ID=14*IBLOK
         DO 117 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM16(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM56(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM123(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM46(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM60(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM121(N4,IDL,IJ)
c
            ENDDO
 117     CONTINUE
c
         ID=15*IBLOK
         DO 118 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM17(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM62(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM129(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM47(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM66(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM127(N4,IDL,IJ)
c
            ENDDO
 118     CONTINUE
c
         ID=16*IBLOK
         DO 119 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM18(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM68(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM135(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM51(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM72(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM133(N4,IDL,IJ)
c
            ENDDO
 119     CONTINUE
c
         ID=17*IBLOK
         DO 120 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM19(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM74(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM141(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM54(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM78(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM139(N4,IDL,IJ)
c
            ENDDO
 120     CONTINUE
c
         ID=18*IBLOK
         DO 121 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM20(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM80(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM148(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM58(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM84(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM145(N4,IDL,IJ)
c               
            ENDDO
 121     CONTINUE
c
         ID=19*IBLOK
         DO 122 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM21(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM86(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM149(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM59(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM94(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM154(N4,IDL,IJ)
c
            ENDDO
 122     CONTINUE
c
         ID=20*IBLOK
         DO 123 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM22(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM87(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM158(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM64(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM95(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM155(N4,IDL,IJ)
c
            ENDDO
 123     CONTINUE
C
         ID=21*IBLOK
         DO 724 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM23(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM98(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM159(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM65(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM102(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM164(N4,IDL,IJ)
c               
            ENDDO
 724     CONTINUE
C
         ID=22*IBLOK
         DO 124 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM24(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM227(N4,IDL,IJ)               
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM168(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM70(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM108(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM165(N4,IDL,IJ)
c
            ENDDO
 124     CONTINUE
C
         ID=23*IBLOK
         DO 125 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM25(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM104(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM169(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM71(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM114(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM174(N4,IDL,IJ)
c
            ENDDO
 125     CONTINUE
C
         ID=24*IBLOK
         DO 126 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM26(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM110(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM178(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM76(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM120(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM175(N4,IDL,IJ)
c
            ENDDO
 126     CONTINUE
C
         ID=25*IBLOK
         DO 127 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM27(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM116(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM179(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM77(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM126(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM184(N4,IDL,IJ)
c
            ENDDO
 127     CONTINUE
C
         ID=26*IBLOK
         DO 128 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM28(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM122(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM188(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM82(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM132(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM185(N4,IDL,IJ)
c
            ENDDO
 128     CONTINUE
C
         ID=27*IBLOK
         DO 129 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM29(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM128(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM189(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM83(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM138(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM194(N4,IDL,IJ)
c
            ENDDO
 129     CONTINUE
C
         ID=28*IBLOK
         DO 130 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM30(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM134(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM198(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM90(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM144(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM195(N4,IDL,IJ)
c
            ENDDO
 130     CONTINUE
C
         ID=29*IBLOK
         DO 131 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM31(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM140(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM199(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM91(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM152(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM204(N4,IDL,IJ)
c              
            ENDDO
 131     CONTINUE
C
         ID=30*IBLOK
         DO 132 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM32(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM146(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM208(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM92(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM153(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM205(N4,IDL,IJ)
c
            ENDDO
 132     CONTINUE
C
         ID=31*IBLOK
         DO 133 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM33(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM147(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM209(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM93(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM162(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM214(N4,IDL,IJ)
c
            ENDDO
 133     CONTINUE
C
         ID=32*IBLOK
         DO 134 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM34(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM156(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM218(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM100(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM163(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM215(N4,IDL,IJ)
c
            ENDDO
 134     CONTINUE
C
         ID=33*IBLOK
         DO 135 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM35(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM157(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM219(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM101(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM172(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM224(N4,IDL,IJ)
c
            ENDDO
 135     CONTINUE
C
         ID=34*IBLOK
         DO 136 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM36(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM166(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM228(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM106(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM173(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM225(N4,IDL,IJ)
c
            ENDDO
 136     CONTINUE
C
         ID=35*IBLOK
         DO 137 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM37(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM167(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM229(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM107(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM182(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM234(N4,IDL,IJ)
c
            ENDDO
 137     CONTINUE
C
         ID=36*IBLOK
         DO 138 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM38(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM176(N4,IDL,IJ)
               ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM99(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM112(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM183(N4,IDL,IJ)
               ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM235(N4,IDL,IJ)
c
            ENDDO
 138     CONTINUE
C
         ID=37*IBLOK
         DO 139 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM39(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM177(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM113(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM192(N4,IDL,IJ)
c
            ENDDO
 139     CONTINUE
C
         ID=38*IBLOK
         DO 140 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM40(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM186(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM118(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM193(N4,IDL,IJ)
c
            ENDDO
 140     CONTINUE
C
         ID=39*IBLOK
         DO 141 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM41(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM187(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM119(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM202(N4,IDL,IJ)
c
            ENDDO
 141     CONTINUE
C
         ID=40*IBLOK
         DO 142 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM42(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM196(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM124(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM203(N4,IDL,IJ)
c
            ENDDO
 142     CONTINUE
C
         ID=41*IBLOK
         DO 143 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM43(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM197(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM125(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM212(N4,IDL,IJ)
c
            ENDDO
 143     CONTINUE
C
         ID=42*IBLOK
         DO 144 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM44(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM206(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM130(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM213(N4,IDL,IJ)
c
            ENDDO
 144     CONTINUE
C
         ID=43*IBLOK
         DO 145 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM45(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM207(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM131(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM222(N4,IDL,IJ)
c
            ENDDO
 145     CONTINUE
C
         ID=44*IBLOK
         DO 146 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM46(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM216(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM136(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM223(N4,IDL,IJ)
c
            ENDDO
 146     CONTINUE
C
         ID=45*IBLOK
         DO 147 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM47(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM217(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM137(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM232(N4,IDL,IJ)
c
            ENDDO
 147     CONTINUE
C
         ID=46*IBLOK
         DO 148 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM48(N4,IDL,IJ)
c
               ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM226(N4,IDL,IJ)
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM142(N4,IDL,IJ)
               ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM233(N4,IDL,IJ)
c
            ENDDO
 148     CONTINUE
C
         ID=47*IBLOK
         DO 149 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM49(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM143(N4,IDL,IJ)
c
            ENDDO
 149     CONTINUE
C
         ID=48*IBLOK
         DO 150 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM50(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM150(N4,IDL,IJ)
c
            ENDDO
 150     CONTINUE
C
         ID=49*IBLOK
         DO 151 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM51(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM151(N4,IDL,IJ)
c
            ENDDO
 151     CONTINUE
C
         ID=50*IBLOK
         DO 152 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM52(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM160(N4,IDL,IJ)
c
            ENDDO
 152     CONTINUE
C
         ID=51*IBLOK
         DO 153 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM53(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM161(N4,IDL,IJ)
c
            ENDDO
 153     CONTINUE
C
         ID=52*IBLOK
         DO 154 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM54(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM170(N4,IDL,IJ)
c
            ENDDO
 154     CONTINUE
C
         ID=53*IBLOK
         DO 155 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM55(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM171(N4,IDL,IJ)
c
            ENDDO
 155     CONTINUE
C
         ID=54*IBLOK
         DO 156 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM56(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM180(N4,IDL,IJ)
c
            ENDDO
 156     CONTINUE
C
         ID=55*IBLOK
         DO 157 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM57(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM181(N4,IDL,IJ)
c
            ENDDO
 157     CONTINUE
C
         ID=56*IBLOK
         DO 158 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM58(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM190(N4,IDL,IJ)
c
            ENDDO
 158     CONTINUE
C
         ID=57*IBLOK
         DO 159 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM59(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM191(N4,IDL,IJ)
c
            ENDDO
 159     CONTINUE         
C     
         ID=58*IBLOK
         DO 160 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM60(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM200(N4,IDL,IJ)
c
            ENDDO
 160     CONTINUE
C
         ID=59*IBLOK
         DO 161 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM61(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM201(N4,IDL,IJ)
c
            ENDDO
 161     CONTINUE
C
         ID=60*IBLOK
         DO 162 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM62(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM210(N4,IDL,IJ)
c
            ENDDO
 162     CONTINUE
C
         ID=61*IBLOK
         DO 163 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM63(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM211(N4,IDL,IJ)
c
            ENDDO
 163     CONTINUE
C
         ID=62*IBLOK
         DO 164 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM64(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM220(N4,IDL,IJ)
c
            ENDDO
 164     CONTINUE
C
         ID=63*IBLOK
         DO 165 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM65(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM221(N4,IDL,IJ)
c
            ENDDO
 165     CONTINUE
C
         ID=64*IBLOK
         DO 166 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM66(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM230(N4,IDL,IJ)
c
            ENDDO
 166     CONTINUE
C
         ID=65*IBLOK
         DO 167 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM67(N4,IDL,IJ)
c
               ALINEMOMJ1(IJ,N4,ID)=ALINEMOM231(N4,IDL,IJ)
c
            ENDDO
 167     CONTINUE
C
         ID=66*IBLOK
         DO 168 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM68(N4,IDL,IJ)
            ENDDO
 168     CONTINUE
C
         ID=67*IBLOK
         DO 169 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM69(N4,IDL,IJ)
            ENDDO
 169     CONTINUE
C
         ID=68*IBLOK
         DO 170 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM70(N4,IDL,IJ)
            ENDDO
 170     CONTINUE
C
         ID=69*IBLOK
         DO 171 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM71(N4,IDL,IJ)
            ENDDO
 171     CONTINUE
C
         ID=70*IBLOK
         DO 172 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM72(N4,IDL,IJ)
            ENDDO
 172     CONTINUE
C
         ID=71*IBLOK
         DO 173 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM73(N4,IDL,IJ)
            ENDDO
 173     CONTINUE
c
         ID=72*IBLOK
         DO 174 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM74(N4,IDL,IJ)
            ENDDO
 174     CONTINUE
c
         ID=73*IBLOK
         DO 175 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM75(N4,IDL,IJ)
            ENDDO
 175     CONTINUE
C
         ID=74*IBLOK
         DO 176 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM76(N4,IDL,IJ)
            ENDDO
 176     CONTINUE         
C
         ID=75*IBLOK
         DO 177 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM77(N4,IDL,IJ)
            ENDDO
 177     CONTINUE
C
         ID=76*IBLOK
         DO 178 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM78(N4,IDL,IJ)
            ENDDO
 178     CONTINUE
C
         ID=77*IBLOK
         DO 179 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM79(N4,IDL,IJ)
            ENDDO
 179     CONTINUE
C
         ID=78*IBLOK
         DO 180 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM80(N4,IDL,IJ)
            ENDDO
 180     CONTINUE
C
         ID=79*IBLOK
         DO 181 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM81(N4,IDL,IJ)
            ENDDO
 181     CONTINUE
C
         ID=80*IBLOK
         DO 182 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM82(N4,IDL,IJ)
            ENDDO
 182     CONTINUE
C
         ID=81*IBLOK
         DO 183 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM83(N4,IDL,IJ)
            ENDDO
 183     CONTINUE
c
         ID=81*IBLOK
         DO 184 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM84(N4,IDL,IJ)
            ENDDO
 184     CONTINUE
c
         ID=83*IBLOK
         DO 185 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM85(N4,IDL,IJ)
            ENDDO
 185     CONTINUE
C
         ID=84*IBLOK
         DO 186 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM86(N4,IDL,IJ)
            ENDDO
 186     CONTINUE         
C
         ID=85*IBLOK
         DO 187 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM87(N4,IDL,IJ)
            ENDDO
 187     CONTINUE
C
         ID=86*IBLOK
         DO 188 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM88(N4,IDL,IJ)
            ENDDO
 188     CONTINUE
C
         ID=87*IBLOK
         DO 189 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM89(N4,IDL,IJ)
            ENDDO
 189     CONTINUE
C
         ID=88*IBLOK
         DO 190 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM90(N4,IDL,IJ)
            ENDDO
 190     CONTINUE
C
         ID=89*IBLOK
         DO 191 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM91(N4,IDL,IJ)
            ENDDO
 191     CONTINUE
C
         ID=90*IBLOK
         DO 192 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM92(N4,IDL,IJ)
            ENDDO
 192     CONTINUE
C
         ID=91*IBLOK
         DO 193 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM93(N4,IDL,IJ)
            ENDDO
 193     CONTINUE
c
         ID=92*IBLOK
         DO 194 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM94(N4,IDL,IJ)
            ENDDO
 194     CONTINUE
c
         ID=93*IBLOK
         DO 195 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM95(N4,IDL,IJ)
            ENDDO
 195     CONTINUE
C
         ID=94*IBLOK
         DO 196 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM96(N4,IDL,IJ)
            ENDDO
 196     CONTINUE         
C
         ID=95*IBLOK
         DO 197 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM97(N4,IDL,IJ)
            ENDDO
 197     CONTINUE
C
         ID=96*IBLOK
         DO 198 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM98(N4,IDL,IJ)
            ENDDO
 198     CONTINUE
C
         ID=97*IBLOK
         DO 199 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM99(N4,IDL,IJ)
            ENDDO
 199     CONTINUE
C
         ID=98*IBLOK
         DO 200 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM100(N4,IDL,IJ)
            ENDDO
 200       CONTINUE
C
         ID=99*IBLOK
         DO 201 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM101(N4,IDL,IJ)
            ENDDO
 201     CONTINUE
C
         ID=100*IBLOK
         DO 202 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM102(N4,IDL,IJ)
            ENDDO
 202     CONTINUE
C
         ID=101*IBLOK
         DO 203 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM103(N4,IDL,IJ)
            ENDDO
 203     CONTINUE
C
         ID=102*IBLOK
         DO 204 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM104(N4,IDL,IJ)
            ENDDO
 204     CONTINUE
C
         ID=103*IBLOK
         DO 205 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM105(N4,IDL,IJ)
            ENDDO
 205     CONTINUE
C
         ID=104*IBLOK
         DO 206 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM106(N4,IDL,IJ)
            ENDDO
 206     CONTINUE
C
         ID=105*IBLOK
         DO 207 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM107(N4,IDL,IJ)
            ENDDO
 207     CONTINUE
C
         ID=106*IBLOK
         DO 208 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM108(N4,IDL,IJ)
            ENDDO
 208     CONTINUE
C
         ID=107*IBLOK
         DO 209 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM109(N4,IDL,IJ)
            ENDDO
 209     CONTINUE
C
         ID=108*IBLOK
         DO 210 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM110(N4,IDL,IJ)
            ENDDO
 210     CONTINUE
C
         ID=109*IBLOK
         DO 211 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM111(N4,IDL,IJ)
            ENDDO
 211     CONTINUE
C
         ID=110*IBLOK
         DO 212 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM112(N4,IDL,IJ)
            ENDDO
 212     CONTINUE
C
         ID=111*IBLOK
         DO 213 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM113(N4,IDL,IJ)
            ENDDO
 213     CONTINUE
C
         ID=112*IBLOK
         DO 214 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM114(N4,IDL,IJ)
            ENDDO
 214     CONTINUE
C
         ID=113*IBLOK
         DO 215 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM115(N4,IDL,IJ)
            ENDDO
 215     CONTINUE
C
         ID=114*IBLOK
         DO 216 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM116(N4,IDL,IJ)
            ENDDO
 216     CONTINUE
C
         ID=115*IBLOK
         DO 217 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM117(N4,IDL,IJ)
            ENDDO
 217     CONTINUE
C
         ID=116*IBLOK
         DO 218 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM118(N4,IDL,IJ)
            ENDDO
 218     CONTINUE
C
         ID=117*IBLOK
         DO 219 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM119(N4,IDL,IJ)
            ENDDO
 219     CONTINUE
C
         ID=118*IBLOK
         DO 220 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM120(N4,IDL,IJ)
            ENDDO
 220     CONTINUE
C
         ID=119*IBLOK
         DO 221 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM121(N4,IDL,IJ)
            ENDDO
 221     CONTINUE
C
         ID=120*IBLOK
         DO 222 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM122(N4,IDL,IJ)
            ENDDO
 222     CONTINUE
C
         ID=121*IBLOK
         DO 223 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM123(N4,IDL,IJ)
            ENDDO
 223     CONTINUE
C
         ID=122*IBLOK
         DO 224 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM124(N4,IDL,IJ)
            ENDDO
 224     CONTINUE
C
         ID=123*IBLOK
         DO 225 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM125(N4,IDL,IJ)
            ENDDO
 225     CONTINUE
C
         ID=124*IBLOK
         DO 226 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM126(N4,IDL,IJ)
            ENDDO
 226     CONTINUE
C
         ID=125*IBLOK
         DO 227 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM127(N4,IDL,IJ)
            ENDDO
 227     CONTINUE
C
         ID=126*IBLOK
         DO 228 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM128(N4,IDL,IJ)
            ENDDO
 228     CONTINUE
C
         ID=127*IBLOK
         DO 229 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM129(N4,IDL,IJ)
            ENDDO
 229     CONTINUE
C
         ID=128*IBLOK
         DO 230 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM130(N4,IDL,IJ)
            ENDDO
 230     CONTINUE
C
         ID=129*IBLOK
         DO 231 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM131(N4,IDL,IJ)
            ENDDO
 231     CONTINUE
C
         ID=130*IBLOK
         DO 232 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM132(N4,IDL,IJ)
            ENDDO
 232     CONTINUE
C
         ID=131*IBLOK
         DO 233 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM133(N4,IDL,IJ)
            ENDDO
 233     CONTINUE
C
         ID=132*IBLOK
         DO 234 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM134(N4,IDL,IJ)
            ENDDO
 234     CONTINUE
C
         ID=133*IBLOK
         DO 235 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM135(N4,IDL,IJ)
            ENDDO
 235     CONTINUE
C
         ID=134*IBLOK
         DO 236 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM136(N4,IDL,IJ)
            ENDDO
 236     CONTINUE
C
         ID=135*IBLOK
         DO 237 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM137(N4,IDL,IJ)
            ENDDO
 237     CONTINUE
C
         ID=136*IBLOK
         DO 238 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM138(N4,IDL,IJ)
            ENDDO
 238     CONTINUE
C
         ID=137*IBLOK
         DO 239 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM139(N4,IDL,IJ)
            ENDDO
 239     CONTINUE
C
         ID=138*IBLOK
         DO 240 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM140(N4,IDL,IJ)
            ENDDO
 240     CONTINUE
C
         ID=139*IBLOK
         DO 241 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM141(N4,IDL,IJ)
            ENDDO
 241     CONTINUE
C
         ID=140*IBLOK
         DO 242 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM142(N4,IDL,IJ)
            ENDDO
 242     CONTINUE         
C
         ID=141*IBLOK
         DO 243 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM143(N4,IDL,IJ)
            ENDDO
 243     CONTINUE
C
         ID=142*IBLOK
         DO 244 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM144(N4,IDL,IJ)
            ENDDO
 244     CONTINUE
C
         ID=143*IBLOK
         DO 245 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM145(N4,IDL,IJ)
            ENDDO
 245     CONTINUE
CCCCCCCCCCCCCCCCCCC
         ID=144*IBLOK
         DO 246 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM146(N4,IDL,IJ)
            ENDDO
 246     CONTINUE
C
         ID=145*IBLOK
         DO 247 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM147(N4,IDL,IJ)
            ENDDO
 247     CONTINUE
C
         ID=146*IBLOK
         DO 248 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM148(N4,IDL,IJ)
            ENDDO
 248     CONTINUE
C
         ID=147*IBLOK
         DO 249 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM149(N4,IDL,IJ)
            ENDDO
 249     CONTINUE
C
         ID=148*IBLOK
         DO 250 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM150(N4,IDL,IJ)
            ENDDO
 250     CONTINUE
C
         ID=149*IBLOK
         DO 251 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM151(N4,IDL,IJ)
            ENDDO
 251     CONTINUE
C
         ID=150*IBLOK
         DO 252 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM152(N4,IDL,IJ)
            ENDDO
 252     CONTINUE
C
         ID=151*IBLOK
         DO 253 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM153(N4,IDL,IJ)
            ENDDO
 253     CONTINUE
C
         ID=152*IBLOK
         DO 254 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM154(N4,IDL,IJ)
            ENDDO
 254     CONTINUE
C
         ID=153*IBLOK
         DO 255 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM155(N4,IDL,IJ)
            ENDDO
 255     CONTINUE
C
         ID=154*IBLOK
         DO 256 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM156(N4,IDL,IJ)
            ENDDO
 256     CONTINUE
C
         ID=155*IBLOK
         DO 257 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM157(N4,IDL,IJ)
            ENDDO
 257     CONTINUE
C
         ID=156*IBLOK
         DO 258 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM158(N4,IDL,IJ)
            ENDDO
 258     CONTINUE
C
         ID=157*IBLOK
         DO 259 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM159(N4,IDL,IJ)
            ENDDO
 259     CONTINUE         
C     
         ID=158*IBLOK
         DO 260 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM160(N4,IDL,IJ)
            ENDDO
 260     CONTINUE
C
         ID=159*IBLOK
         DO 261 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM161(N4,IDL,IJ)
            ENDDO
 261     CONTINUE
C
         ID=160*IBLOK
         DO 262 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM162(N4,IDL,IJ)
            ENDDO
 262     CONTINUE
C
         ID=161*IBLOK
         DO 263 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM163(N4,IDL,IJ)
            ENDDO
 263     CONTINUE
C
         ID=162*IBLOK
         DO 264 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM164(N4,IDL,IJ)
            ENDDO
 264     CONTINUE
C
         ID=163*IBLOK
         DO 265 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM165(N4,IDL,IJ)
            ENDDO
 265     CONTINUE
C
         ID=164*IBLOK
         DO 266 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM166(N4,IDL,IJ)
            ENDDO
 266     CONTINUE
C
         ID=165*IBLOK
         DO 267 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM167(N4,IDL,IJ)
            ENDDO
 267     CONTINUE
C
         ID=166*IBLOK
         DO 268 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM168(N4,IDL,IJ)
            ENDDO
 268     CONTINUE
C
         ID=167*IBLOK
         DO 269 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM169(N4,IDL,IJ)
            ENDDO
 269     CONTINUE
C
         ID=168*IBLOK
         DO 270 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM170(N4,IDL,IJ)
            ENDDO
 270     CONTINUE
C
         ID=169*IBLOK
         DO 271 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM171(N4,IDL,IJ)
            ENDDO
 271     CONTINUE
C
         ID=170*IBLOK
         DO 272 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM172(N4,IDL,IJ)
            ENDDO
 272     CONTINUE
C
         ID=171*IBLOK
         DO 273 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM173(N4,IDL,IJ)
            ENDDO
 273     CONTINUE
c
         ID=172*IBLOK
         DO 274 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM174(N4,IDL,IJ)
            ENDDO
 274     CONTINUE
c
         ID=173*IBLOK
         DO 275 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM175(N4,IDL,IJ)
            ENDDO
 275     CONTINUE
C
         ID=174*IBLOK
         DO 276 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM176(N4,IDL,IJ)
            ENDDO
 276     CONTINUE         
C
         ID=175*IBLOK
         DO 277 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM177(N4,IDL,IJ)
            ENDDO
 277     CONTINUE
C
         ID=176*IBLOK
         DO 278 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM178(N4,IDL,IJ)
            ENDDO
 278     CONTINUE
C
         ID=177*IBLOK
         DO 279 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM179(N4,IDL,IJ)
            ENDDO
 279     CONTINUE
C
         ID=178*IBLOK
         DO 280 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM180(N4,IDL,IJ)
            ENDDO
 280     CONTINUE
C
         ID=179*IBLOK
         DO 281 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM181(N4,IDL,IJ)
            ENDDO
 281     CONTINUE
C
         ID=180*IBLOK
         DO 282 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM182(N4,IDL,IJ)
            ENDDO
 282     CONTINUE
C
         ID=181*IBLOK
         DO 283 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM183(N4,IDL,IJ)
            ENDDO
 283     CONTINUE
c
         ID=181*IBLOK
         DO 284 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM184(N4,IDL,IJ)
            ENDDO
 284     CONTINUE
c
         ID=183*IBLOK
         DO 285 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM185(N4,IDL,IJ)
            ENDDO
 285     CONTINUE
C
         ID=184*IBLOK
         DO 286 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM186(N4,IDL,IJ)
            ENDDO
 286     CONTINUE         
C
         ID=185*IBLOK
         DO 287 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM187(N4,IDL,IJ)
            ENDDO
 287     CONTINUE
C
         ID=186*IBLOK
         DO 288 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM188(N4,IDL,IJ)
            ENDDO
 288     CONTINUE
C
         ID=187*IBLOK
         DO 289 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM189(N4,IDL,IJ)
            ENDDO
 289     CONTINUE
C
         ID=188*IBLOK
         DO 290 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM190(N4,IDL,IJ)
            ENDDO
 290     CONTINUE
C
         ID=189*IBLOK
         DO 291 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM191(N4,IDL,IJ)
            ENDDO
 291     CONTINUE
C
         ID=190*IBLOK
         DO 292 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM192(N4,IDL,IJ)
            ENDDO
 292     CONTINUE
C
         ID=191*IBLOK
         DO 293 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM193(N4,IDL,IJ)
            ENDDO
 293     CONTINUE
c
         ID=192*IBLOK
         DO 294 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM194(N4,IDL,IJ)
            ENDDO
 294     CONTINUE
c
         ID=193*IBLOK
         DO 295 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM195(N4,IDL,IJ)
            ENDDO
 295     CONTINUE
C
         ID=194*IBLOK
         DO 296 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM196(N4,IDL,IJ)
            ENDDO
 296     CONTINUE         
C
         ID=195*IBLOK
         DO 297 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM197(N4,IDL,IJ)
            ENDDO
 297     CONTINUE
C
         ID=196*IBLOK
         DO 298 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM198(N4,IDL,IJ)
            ENDDO
 298     CONTINUE
C
         ID=197*IBLOK
         DO 299 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM199(N4,IDL,IJ)
            ENDDO
 299     CONTINUE
C
         ID=198*IBLOK
         DO 500 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM200(N4,IDL,IJ)
            ENDDO
 500       CONTINUE
C
         ID=199*IBLOK
         DO 501 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM201(N4,IDL,IJ)
            ENDDO
 501     CONTINUE
C
         ID=200*IBLOK
         DO 502 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM202(N4,IDL,IJ)
            ENDDO
 502     CONTINUE
C
         ID=201*IBLOK
         DO 503 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM203(N4,IDL,IJ)
            ENDDO
 503     CONTINUE
C
         ID=202*IBLOK
         DO 504 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM204(N4,IDL,IJ)
            ENDDO
 504     CONTINUE
C
         ID=203*IBLOK
         DO 505 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM205(N4,IDL,IJ)
            ENDDO
 505     CONTINUE
C
         ID=204*IBLOK
         DO 506 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM206(N4,IDL,IJ)
            ENDDO
 506     CONTINUE
C
         ID=205*IBLOK
         DO 507 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM207(N4,IDL,IJ)
            ENDDO
 507     CONTINUE
C
         ID=206*IBLOK
         DO 508 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM208(N4,IDL,IJ)
            ENDDO
 508     CONTINUE
C
         ID=207*IBLOK
         DO 509 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM209(N4,IDL,IJ)
            ENDDO
 509     CONTINUE
C
         ID=208*IBLOK
         DO 510 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM210(N4,IDL,IJ)
            ENDDO
 510     CONTINUE
C
         ID=209*IBLOK
         DO 511 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM211(N4,IDL,IJ)
            ENDDO
 511     CONTINUE
C
         ID=210*IBLOK
         DO 512 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM212(N4,IDL,IJ)
            ENDDO
 512     CONTINUE
C
         ID=211*IBLOK
         DO 513 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM213(N4,IDL,IJ)
            ENDDO
 513     CONTINUE
C
         ID=212*IBLOK
         DO 514 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM214(N4,IDL,IJ)
            ENDDO
 514     CONTINUE
C
         ID=213*IBLOK
         DO 515 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM215(N4,IDL,IJ)
            ENDDO
 515     CONTINUE
C
         ID=214*IBLOK
         DO 516 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM216(N4,IDL,IJ)
            ENDDO
 516     CONTINUE
C
         ID=215*IBLOK
         DO 517 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM217(N4,IDL,IJ)
            ENDDO
 517     CONTINUE
C
         ID=216*IBLOK
         DO 518 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM218(N4,IDL,IJ)
            ENDDO
 518     CONTINUE
C
         ID=217*IBLOK
         DO 519 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM219(N4,IDL,IJ)
            ENDDO
 519     CONTINUE
C
         ID=218*IBLOK
         DO 520 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM220(N4,IDL,IJ)
            ENDDO
 520     CONTINUE
C
         ID=219*IBLOK
         DO 521 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM221(N4,IDL,IJ)
            ENDDO
 521     CONTINUE
C
         ID=220*IBLOK
         DO 522 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM222(N4,IDL,IJ)
            ENDDO
 522     CONTINUE
C
         ID=221*IBLOK
         DO 523 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM223(N4,IDL,IJ)
            ENDDO
 523     CONTINUE
C
         ID=222*IBLOK
         DO 524 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM224(N4,IDL,IJ)
            ENDDO
 524     CONTINUE
C
         ID=223*IBLOK
         DO 525 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM225(N4,IDL,IJ)
            ENDDO
 525     CONTINUE
C
         ID=224*IBLOK
         DO 526 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM226(N4,IDL,IJ)
            ENDDO
 526     CONTINUE
C
         ID=225*IBLOK
         DO 527 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM227(N4,IDL,IJ)
            ENDDO
 527     CONTINUE
C
         ID=226*IBLOK
         DO 528 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM228(N4,IDL,IJ)
            ENDDO
 528     CONTINUE
C
         ID=227*IBLOK
         DO 529 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM229(N4,IDL,IJ)
            ENDDO
 529     CONTINUE
C
         ID=228*IBLOK
         DO 530 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM230(N4,IDL,IJ)
            ENDDO
 530     CONTINUE
C
         ID=229*IBLOK
         DO 531 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM231(N4,IDL,IJ)
            ENDDO
 531     CONTINUE
C
         ID=230*IBLOK
         DO 532 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM232(N4,IDL,IJ)
            ENDDO
 532     CONTINUE
C
         ID=231*IBLOK
         DO 533 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM233(N4,IDL,IJ)
            ENDDO
 533     CONTINUE
C
         ID=232*IBLOK
         DO 534 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM234(N4,IDL,IJ)
            ENDDO
 534     CONTINUE
C
         ID=233*IBLOK
         DO 535 IDL=1,IBLOK
            ID=ID+1
            DO IJ=1,2
               ALINEPMOM(IJ,N4,ID)=ALINEMOM235(N4,IDL,IJ)
            ENDDO
 535     CONTINUE
C
C
 102  CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC FULL J MOM-LINES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO N4=1, LX4
C         
         DO IDL=1, NOPJ0PMOM
            DO IJ=1, 2
               ALINEMOMJ0(IJ,N4,IDL)=ALINEMOMJ0P(IJ,N4,IDL)
            ENDDO 
         ENDDO
C
         ID=NOPJ0PMOM
         DO IDL=1, NOPJ0MMOM
            ID=ID+1
            DO IJ=1, 2
               ALINEMOMJ0(IJ,N4,ID)=ALINEMOMJ0M(IJ,N4,IDL)
            ENDDO 
         ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO IDL=1, NOPJ2PMOM
            DO IJ=1, 2
               ALINEMOMJ2(IJ,N4,IDL)=ALINEMOMJ2P(IJ,N4,IDL)
            ENDDO 
         ENDDO
C
         ID=NOPJ2PMOM
         DO IDL=1, NOPJ2MMOM
            ID=ID+1
            DO IJ=1, 2
               ALINEMOMJ2(IJ,N4,ID)=ALINEMOMJ2M(IJ,N4,IDL)
            ENDDO 
         ENDDO
C         
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC FULL J ALL MOM-LINES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO N4=1, LX4
C
         DO IDL=1, NOPFULJ0MOM
            DO IJ=1, 2
               ALINEMOMJALL(IJ,N4,IDL)=ALINEMOMJ0(IJ,N4,IDL)
            ENDDO 
         ENDDO
C
         ID=NOPFULJ0MOM
         DO IDL=1, NOPFULJ1MOM
            ID=ID+1
            DO IJ=1, 2
               ALINEMOMJALL(IJ,N4,ID)=ALINEMOMJ1(IJ,N4,IDL)
            ENDDO 
         ENDDO         
C
         ID=NOPFULJ0MOM+NOPFULJ1MOM
         DO IDL=1, NOPFULJ2MOM
            ID=ID+1
            DO IJ=1, 2
               ALINEMOMJALL(IJ,N4,ID)=ALINEMOMJ2(IJ,N4,IDL)
            ENDDO 
         ENDDO   
C         
      ENDDO
C***************************************************************
      DO 605 N4=1,LX4
C***************************************************************     
         DO 606 ID=1, NTOTAL
            ALLP=ALINEP(N4,ID)
            AVACLP(JBIN,ID)=AVACLP(JBIN,ID)+ALLP
 606     CONTINUE
C
         DO 607 ID=1, NOPJ0PP
            ALLPJ0PP=ALINEJ0PP(N4,ID)
            AVACLJ0PP(JBIN,ID)=AVACLJ0PP(JBIN,ID)+ALLPJ0PP
 607     CONTINUE         
C
         DO 608 ID=1, NOPJ0PM
            ALLPJ0PM=ALINEJ0PM(N4,ID)
            AVACLJ0PM(JBIN,ID)=AVACLJ0PM(JBIN,ID)+ALLPJ0PM
 608     CONTINUE                  
C
         DO 609 ID=1, NOPJ0MP
            ALLPJ0MP=ALINEJ0MP(N4,ID)
            AVACLJ0MP(JBIN,ID)=AVACLJ0MP(JBIN,ID)+ALLPJ0MP
 609     CONTINUE         
C
         DO 610 ID=1, NOPJ0MM
            ALLPJ0MM=ALINEJ0MM(N4,ID)
            AVACLJ0MM(JBIN,ID)=AVACLJ0MM(JBIN,ID)+ALLPJ0MM
 610     CONTINUE         
C
         DO 611 ID=1, NOPJ1P
            ALLPJ1P=ALINEJ1P(N4,ID)
            AVACLJ1P(JBIN,ID)=AVACLJ1P(JBIN,ID)+ALLPJ1P
 611     CONTINUE 
C
         DO 612 ID=1, NOPJ1M
            ALLPJ1M=ALINEJ1M(N4,ID)
            AVACLJ1M(JBIN,ID)=AVACLJ1M(JBIN,ID)+ALLPJ1M
 612     CONTINUE          
C
         DO 613 ID=1, NOPJ2PP
            ALLPJ2PP=ALINEJ2PP(N4,ID)
            AVACLJ2PP(JBIN,ID)=AVACLJ2PP(JBIN,ID)+ALLPJ2PP
 613     CONTINUE         
C
         DO 614 ID=1, NOPJ2PM
            ALLPJ2PM=ALINEJ2PM(N4,ID)
            AVACLJ2PM(JBIN,ID)=AVACLJ2PM(JBIN,ID)+ALLPJ2PM
 614     CONTINUE                  
C
         DO 615 ID=1, NOPJ2MP
            ALLPJ2MP=ALINEJ2MP(N4,ID)
            AVACLJ2MP(JBIN,ID)=AVACLJ2MP(JBIN,ID)+ALLPJ2MP
 615     CONTINUE         
C
         DO 616 ID=1, NOPJ2MM
            ALLPJ2MM=ALINEJ2MM(N4,ID)
            AVACLJ2MM(JBIN,ID)=AVACLJ2MM(JBIN,ID)+ALLPJ2MM
 616     CONTINUE          
C
         DO 617 ID=1, NOPFULJ0
            ALLPJ0=ALINEJ0(N4,ID)
            AVACLJ0(JBIN,ID)=AVACLJ0(JBIN,ID)+ALLPJ0
 617     CONTINUE
C
         DO 618 ID=1, NOPFULJ1
            ALLPJ1=ALINEJ1(N4,ID)
            AVACLJ1(JBIN,ID)=AVACLJ1(JBIN,ID)+ALLPJ1
 618     CONTINUE         
C
         DO 619 ID=1, NOPFULJ2
            ALLPJ2=ALINEJ2(N4,ID)
            AVACLJ2(JBIN,ID)=AVACLJ2(JBIN,ID)+ALLPJ2
 619     CONTINUE 
C
         DO 620 ID=1, NTOTAL
            ALLPJALL=ALINEJALL(N4,ID)
            AVACLALL(JBIN,ID)=AVACLALL(JBIN,ID)+ALLPJALL
 620     CONTINUE          
C         
         DO 807 ID=1, NTOTALMOM
            DO IJ=1,2
               ALLPMOM=ALINEPMOM(IJ,N4,ID)
               AVACLPMOM(IJ,JBIN,ID)=AVACLPMOM(IJ,JBIN,ID)+ALLPMOM
            ENDDO
 807     CONTINUE
C
         DO 808 ID=1, NOPJ0PMOM
            DO IJ=1,2
               ALLPMOMJ0P=ALINEMOMJ0P(IJ,N4,ID)
               AVACLMOMJ0P(IJ,JBIN,ID)=AVACLMOMJ0P(IJ,JBIN,ID)
     &+ALLPMOMJ0P
            ENDDO
 808     CONTINUE
C
         DO 809 ID=1, NOPJ0MMOM
            DO IJ=1,2
               ALLPMOMJ0M=ALINEMOMJ0M(IJ,N4,ID)
               AVACLMOMJ0M(IJ,JBIN,ID)=AVACLMOMJ0M(IJ,JBIN,ID)
     &+ALLPMOMJ0M
            ENDDO
 809     CONTINUE         
C
         DO 810 ID=1, NOPJ1MOM
            DO IJ=1,2
               ALLPMOMJ1=ALINEMOMJ1(IJ,N4,ID)
               AVACLMOMJ1(IJ,JBIN,ID)=AVACLMOMJ1(IJ,JBIN,ID)
     &+ALLPMOMJ1
            ENDDO
 810     CONTINUE
C
         DO 811 ID=1, NOPJ2PMOM
            DO IJ=1,2
               ALLPMOMJ2P=ALINEMOMJ2P(IJ,N4,ID)
               AVACLMOMJ2P(IJ,JBIN,ID)=AVACLMOMJ2P(IJ,JBIN,ID)
     &+ALLPMOMJ2P
            ENDDO
 811     CONTINUE
C
         DO 812 ID=1, NOPJ2MMOM
            DO IJ=1,2
               ALLPMOMJ2M=ALINEMOMJ2M(IJ,N4,ID)
               AVACLMOMJ2M(IJ,JBIN,ID)=AVACLMOMJ2M(IJ,JBIN,ID)
     &+ALLPMOMJ2M
            ENDDO
 812     CONTINUE
C
         DO 813 ID=1, NOPFULJ0MOM
            DO IJ=1,2
               ALLPMOMJ0=ALINEMOMJ0(IJ,N4,ID)
               AVACLMOMJ0(IJ,JBIN,ID)=AVACLMOMJ0(IJ,JBIN,ID)
     &+ALLPMOMJ0
            ENDDO
 813     CONTINUE
C
         DO 814 ID=1, NOPFULJ2MOM
            DO IJ=1,2
               ALLPMOMJ2=ALINEMOMJ2(IJ,N4,ID)
               AVACLMOMJ2(IJ,JBIN,ID)=AVACLMOMJ2(IJ,JBIN,ID)
     &+ALLPMOMJ2
            ENDDO
 814     CONTINUE
C
         DO 815 ID=1, NTOTALMOM
            DO IJ=1,2
               ALLPMOMJALL=ALINEMOMJALL(IJ,N4,ID)
               AVACLMOMALL(IJ,JBIN,ID)=AVACLMOMALL(IJ,JBIN,ID)
     &+ALLPMOMJALL
            ENDDO
 815     CONTINUE
C         
 605  CONTINUE
C***************************************************************
C     FULL CORRELATOR
C***************************************************************      
      INO=0
      DO ID2=1, NTOTAL
         DO ID1=1, NTOTAL
C
            DO N4=1,LX4
               DO NT=1,LMAXIR
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLP(JBIN,NT,ID1,ID2)=ACORLP(JBIN,NT,ID1,ID2)
     &            +(ALINEP(N4,ID1)*CONJG(ALINEP(N4X,ID2))
     &            +CONJG(ALINEP(N4,ID2))*ALINEP(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=0, Pp=+, Pr=+, q=0 correlator
C***************************************************************
      INO=0
      DO ID2=1,NOPJ0PP
         DO ID1=1,NOPJ0PP
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ0PP(JBIN,NT,ID1,ID2)=ACORLJ0PP(JBIN,NT,ID1,ID2)
     &            +(ALINEJ0PP(N4,ID1)*CONJG(ALINEJ0PP(N4X,ID2))
     &            +CONJG(ALINEJ0PP(N4,ID2))*ALINEJ0PP(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=0, Pp=+, Pr=-, q=0 correlator
C***************************************************************
      INO=0   
      DO ID2=1,NOPJ0PM
         DO ID1=1,NOPJ0PM
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ0PM(JBIN,NT,ID1,ID2)=ACORLJ0PM(JBIN,NT,ID1,ID2)
     &            +(ALINEJ0PM(N4,ID1)*CONJG(ALINEJ0PM(N4X,ID2))
     &            +CONJG(ALINEJ0PM(N4,ID2))*ALINEJ0PM(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=0, Pp=-, Pr=+, q=0 correlator
C***************************************************************
      INO=0 
      DO ID2=1,NOPJ0MP
         DO ID1=1,NOPJ0MP
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ0MP(JBIN,NT,ID1,ID2)=ACORLJ0MP(JBIN,NT,ID1,ID2)
     &            +(ALINEJ0MP(N4,ID1)*CONJG(ALINEJ0MP(N4X,ID2))
     &            +CONJG(ALINEJ0MP(N4,ID2))*ALINEJ0MP(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=0, Pp=-, Pr=-, q=0 correlator
C***************************************************************
      INO=0 
      DO ID2=1,NOPJ0MM
         DO ID1=1,NOPJ0MM
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ0MM(JBIN,NT,ID1,ID2)=ACORLJ0MM(JBIN,NT,ID1,ID2)
     &            +(ALINEJ0MM(N4,ID1)*CONJG(ALINEJ0MM(N4X,ID2))
     &            +CONJG(ALINEJ0MM(N4,ID2))*ALINEJ0MM(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=1, Pr=+, q=0 correlator
C***************************************************************
      INO=0       
      DO ID2=1,NOPJ1P
         DO ID1=1,NOPJ1P
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ1P(JBIN,NT,ID1,ID2)=ACORLJ1P(JBIN,NT,ID1,ID2)
     &            +(ALINEJ1P(N4,ID1)*CONJG(ALINEJ1P(N4X,ID2))
     &            +CONJG(ALINEJ1P(N4,ID2))*ALINEJ1P(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=1, Pr=-, q=0 correlator
C***************************************************************
      INO=0              
      DO ID2=1,NOPJ1M
         DO ID1=1,NOPJ1M
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ1M(JBIN,NT,ID1,ID2)=ACORLJ1M(JBIN,NT,ID1,ID2)
     &            +(ALINEJ1M(N4,ID1)*CONJG(ALINEJ1M(N4X,ID2))
     &            +CONJG(ALINEJ1M(N4,ID2))*ALINEJ1M(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=2, Pp=+, Pr=+, q=0 correlator
C***************************************************************
      INO=0   
      DO ID2=1,NOPJ2PP
         DO ID1=1,NOPJ2PP
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ2PP(JBIN,NT,ID1,ID2)=ACORLJ2PP(JBIN,NT,ID1,ID2)
     &            +(ALINEJ2PP(N4,ID1)*CONJG(ALINEJ2PP(N4X,ID2))
     &            +CONJG(ALINEJ2PP(N4,ID2))*ALINEJ2PP(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=2, Pp=+, Pr=-, q=0 correlator
C***************************************************************
      INO=0       
      DO ID2=1,NOPJ2PM
         DO ID1=1,NOPJ2PM
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ2PM(JBIN,NT,ID1,ID2)=ACORLJ2PM(JBIN,NT,ID1,ID2)
     &            +(ALINEJ2PM(N4,ID1)*CONJG(ALINEJ2PM(N4X,ID2))
     &            +CONJG(ALINEJ2PM(N4,ID2))*ALINEJ2PM(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C     
         ENDDO
      ENDDO
C***************************************************************
C     J=2, Pp=-, Pr=+, q=0 correlator
C***************************************************************
      INO=0   
      DO ID2=1,NOPJ2MP
         DO ID1=1,NOPJ2MP
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ2MP(JBIN,NT,ID1,ID2)=ACORLJ2MP(JBIN,NT,ID1,ID2)
     &            +(ALINEJ2MP(N4,ID1)*CONJG(ALINEJ2MP(N4X,ID2))
     &            +CONJG(ALINEJ2MP(N4,ID2))*ALINEJ2MP(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=2, Pp=-, Pr=-, q=0 correlator
C***************************************************************
      INO=0  
      DO ID2=1,NOPJ2MM
         DO ID1=1,NOPJ2MM
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ2MM(JBIN,NT,ID1,ID2)=ACORLJ2MM(JBIN,NT,ID1,ID2)
     &            +(ALINEJ2MM(N4,ID1)*CONJG(ALINEJ2MM(N4X,ID2))
     &            +CONJG(ALINEJ2MM(N4,ID2))*ALINEJ2MM(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=0 FULL, q=0 correlator
C***************************************************************
      INO=0        
      DO ID2=1,NOPFULJ0
         DO ID1=1,NOPFULJ0
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ0(JBIN,NT,ID1,ID2)=ACORLJ0(JBIN,NT,ID1,ID2)
     &            +(ALINEJ0(N4,ID1)*CONJG(ALINEJ0(N4X,ID2))
     &            +CONJG(ALINEJ0(N4,ID2))*ALINEJ0(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=1 FULL, q=0 correlator
C***************************************************************      
      INO=0   
      DO ID2=1,NOPFULJ1
         DO ID1=1,NOPFULJ1
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ1(JBIN,NT,ID1,ID2)=ACORLJ1(JBIN,NT,ID1,ID2)
     &            +(ALINEJ1(N4,ID1)*CONJG(ALINEJ1(N4X,ID2))
     &            +CONJG(ALINEJ1(N4,ID2))*ALINEJ1(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=2 FULL, q=0 correlator
C***************************************************************
      INO=0
      DO ID2=1,NOPFULJ2
         DO ID1=1,NOPFULJ2
C
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLJ2(JBIN,NT,ID1,ID2)=ACORLJ2(JBIN,NT,ID1,ID2)
     &            +(ALINEJ2(N4,ID1)*CONJG(ALINEJ2(N4X,ID2))
     &            +CONJG(ALINEJ2(N4,ID2))*ALINEJ2(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C     J=ALL FULL, q=0 correlator
C***************************************************************      
      INO=0
      DO ID2=1,NTOTAL
         DO ID1=1,NTOTAL
C
            DO N4=1,LX4
               DO NT=1,LMAXIR
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  ACORLALL(JBIN,NT,ID1,ID2)=ACORLALL(JBIN,NT,ID1,ID2)
     &            +(ALINEJALL(N4,ID1)*CONJG(ALINEJALL(N4X,ID2))
     &            +CONJG(ALINEJALL(N4,ID2))*ALINEJALL(N4X,ID1))*0.5
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C***************************************************************
C***************************************************************      
C***************************************************************
CCCCCC   Momentum
C***************************************************************
C***************************************************************      
C***************************************************************      
C     ID12=0
      INO=0
      DO ID2=1, NTOTALMOM
         DO ID1=1, NTOTALMOM
C     ID12=ID12+1
C     
            DO N4=1,LX4
               DO NT=1,LMAXIR
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOM(IJ,JBIN,NT,ID1,ID2)=ACORLPMOM(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEPMOM(IJ,N4,ID1)*CONJG(ALINEPMOM(IJ,N4X,ID2))
     &+CONJG(ALINEPMOM(IJ,N4,ID2))*ALINEPMOM(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
c
         ENDDO
      ENDDO
C****************************************************************
C     J=0, Pp=+, q=1,2
C****************************************************************
      INO=0             
      DO ID2=1,NOPJ0PMOM
         DO ID1=1,NOPJ0PMOM
C     
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOMJ0P(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ0P(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJ0P(IJ,N4,ID1)*CONJG(ALINEMOMJ0P(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJ0P(IJ,N4,ID2))*ALINEMOMJ0P(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C               
         ENDDO
      ENDDO
C****************************************************************
C     J=0, Pp=-, q=1,2
C****************************************************************
      INO=0    
      DO ID2=1,NOPJ0MMOM
         DO ID1=1,NOPJ0MMOM
C     
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOMJ0M(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ0M(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJ0M(IJ,N4,ID1)*CONJG(ALINEMOMJ0M(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJ0M(IJ,N4,ID2))*ALINEMOMJ0M(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C****************************************************************
C     J=1, q=1,2
C****************************************************************
      INO=0          
      DO ID2=1,NOPJ1MOM
         DO ID1=1,NOPJ1MOM
C     
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOMJ1(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ1(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJ1(IJ,N4,ID1)*CONJG(ALINEMOMJ1(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJ1(IJ,N4,ID2))*ALINEMOMJ1(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C****************************************************************
C     J=2, Pp=+, q=1,2
C****************************************************************
      INO=0           
      DO ID2=1,NOPJ2PMOM
         DO ID1=1,NOPJ2PMOM
C     
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOMJ2P(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ2P(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJ2P(IJ,N4,ID1)*CONJG(ALINEMOMJ2P(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJ2P(IJ,N4,ID2))*ALINEMOMJ2P(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C                
         ENDDO
      ENDDO
C****************************************************************
C     J=2, Pp=-, q=1,2
C****************************************************************
      INO=0         
      DO ID2=1,NOPJ2MMOM
         DO ID1=1,NOPJ2MMOM
C     
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOMJ2M(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ2M(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJ2M(IJ,N4,ID1)*CONJG(ALINEMOMJ2M(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJ2M(IJ,N4,ID2))*ALINEMOMJ2M(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C****************************************************************
C     J=0, q=1,2
C****************************************************************
      INO=0         
      DO ID2=1,NOPFULJ0MOM
         DO ID1=1,NOPFULJ0MOM
C     
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOMJ0(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ0(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJ0(IJ,N4,ID1)*CONJG(ALINEMOMJ0(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJ0(IJ,N4,ID2))*ALINEMOMJ0(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C****************************************************************
C     J=2, q=1,2
C****************************************************************
      INO=0         
      DO ID2=1,NOPFULJ2MOM
         DO ID1=1,NOPFULJ2MOM
C     
            DO N4=1,LX4
               DO NT=1,LMAX
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
      ACORLPMOMJ2(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ2(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJ2(IJ,N4,ID1)*CONJG(ALINEMOMJ2(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJ2(IJ,N4,ID2))*ALINEMOMJ2(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C  
         ENDDO
      ENDDO
C****************************************************************
C     J=ALL, q=1,2
C****************************************************************
      INO=0         
      DO ID2=1, NTOTALMOM
         DO ID1=1, NTOTALMOM
C     
            DO N4=1,LX4
               DO NT=1,LMAXIR
                  N4X=N4+NT-1
                  IF(N4X.GT.LX4) N4X=N4X-LX4
C*******************************************
                  DO IJ=1,2
C*******************************************
                     ACORLPMOMJALL(IJ,JBIN,NT,ID1,ID2)=
     &ACORLPMOMJALL(IJ,JBIN,NT,ID1,ID2)
     &+(ALINEMOMJALL(IJ,N4,ID1)*CONJG(ALINEMOMJALL(IJ,N4X,ID2))
     &+CONJG(ALINEMOMJALL(IJ,N4,ID2))*ALINEMOMJALL(IJ,N4X,ID1))*0.5
C*******************************************
                  ENDDO
C*******************************************
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C********************************************************************
C********************************************************************
C********************************************************************
C********************************************************************
      RETURN
      END
C*********************************************************************
C*********************************************************************
C  Here we set up blocked link pointers...
C*********************************************************************
C*********************************************************************
      SUBROUTINE SETUPB
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
C
      COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
      DIMENSION LB(IBLOK+1),LB1(IBLOK+1),LB12(IBLOK+1)
C
      LX12=LX1*LX2
      LX123=LX12*LX3
      LB(1)=1
      LB1(1)=LX1
      LB12(1)=LX12
      DO 1 IBL=2,IBLOK+1
         LB(IBL)=2*LB(IBL-1)
         LB1(IBL)=2*LB1(IBL-1)
         LB12(IBL)=2*LB12(IBL-1)
 1    CONTINUE
C
      NN=0
      DO L3=1,LX3
         DO L2=1,LX2
            DO L1=1,LX1
               NN=NN+1
               DO 3 ID=1,IBLOK+1
C     
                  NU1=NN+LB(ID)
                  IF((L1+LB(ID)).GT.(LX1*1)) NU1=NU1-LX1
                  IF((L1+LB(ID)).GT.(LX1*2)) NU1=NU1-LX1
                  IF((L1+LB(ID)).GT.(LX1*3)) NU1=NU1-LX1
                  IF((L1+LB(ID)).GT.(LX1*4)) NU1=NU1-LX1
                  IUPB(NN,1,ID)=NU1
                  ND1=NN-LB(ID)
                  IF((L1-LB(ID)).LT.(1-LX1*0)) ND1=ND1+LX1
                  IF((L1-LB(ID)).LT.(1-LX1*1)) ND1=ND1+LX1
                  IF((L1-LB(ID)).LT.(1-LX1*2)) ND1=ND1+LX1
                  IF((L1-LB(ID)).LT.(1-LX1*3)) ND1=ND1+LX1
                  IDNB(NN,1,ID)=ND1
C     
                  NU2=NN+LB1(ID)
                  IF((L2+LB(ID)).GT.(LX2*1)) NU2=NU2-LX12
                  IF((L2+LB(ID)).GT.(LX2*2)) NU2=NU2-LX12
                  IF((L2+LB(ID)).GT.(LX2*3)) NU2=NU2-LX12
                  IF((L2+LB(ID)).GT.(LX2*4)) NU2=NU2-LX12
                  IUPB(NN,2,ID)=NU2
                  ND2=NN-LB1(ID)
                  IF((L2-LB(ID)).LT.(1-LX2*0)) ND2=ND2+LX12
                  IF((L2-LB(ID)).LT.(1-LX2*1)) ND2=ND2+LX12
                  IF((L2-LB(ID)).LT.(1-LX2*2)) ND2=ND2+LX12
                  IF((L2-LB(ID)).LT.(1-LX2*3)) ND2=ND2+LX12
                  IDNB(NN,2,ID)=ND2
C     
                  NU3=NN+LB12(ID)
                  IF((L3+LB(ID)).GT.(LX3*1)) NU3=NU3-LX123
                  IF((L3+LB(ID)).GT.(LX3*2)) NU3=NU3-LX123
                  IF((L3+LB(ID)).GT.(LX3*3)) NU3=NU3-LX123
                  IF((L3+LB(ID)).GT.(LX3*4)) NU3=NU3-LX123
                  IUPB(NN,3,ID)=NU3
                  ND3=NN-LB12(ID)
                  IF((L3-LB(ID)).LT.(1-LX3*0)) ND3=ND3+LX123
                  IF((L3-LB(ID)).LT.(1-LX3*1)) ND3=ND3+LX123
                  IF((L3-LB(ID)).LT.(1-LX3*2)) ND3=ND3+LX123
                  IF((L3-LB(ID)).LT.(1-LX3*3)) ND3=ND3+LX123
                  IDNB(NN,3,ID)=ND3
C     
 3             CONTINUE
            ENDDO
         ENDDO
      ENDDO
C     
      RETURN
      END
C*********************************************************************
C*********************************************************************
C                        SUBROUTINE ACTION                           *
C*********************************************************************
C*********************************************************************
      SUBROUTINE ACTION(IPR,ITER,TOTACT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(NITER=2,NUMBIN=2)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ASTORE/ACTN(NUMBIN,6),PLAQ(NITER,6)
      COMMON/COMMUNICATE/AVACS(NUMBIN),AVACT(NUMBIN)
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      DIMENSION DUM11(NCOL2)
     &,A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2)
      COMPLEX U11,A11,B11,C11,D11,DUM11,ACT
C
      DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN)
      DIMENSION VAL(NUMBIN),AV(NUMBIN)
C
      IBIN=NITER/NUMBIN
C
      IF(IPR.EQ.0)THEN
      JBIN=(ITER-1)/IBIN+1
C
      IF(ITER.EQ.1)THEN
      DO 3 NB=1,NUMBIN
      AVAC(NB)=0.0d0
      AVACSQ(NB)=0.0d0
      AVACS(NB)=0.0d0
      AVACT(NB)=0.0d0
3     CONTINUE
      ENDIF
C
      VACTS=0.0
      VACTT=0.0
C
      DO NN=1,LSIZE
         DO MU=1,3
            M1=NN
C     
            DO 6  IJ=1,NCOL2
               A11(IJ)=U11(IJ,M1,MU)
 6          CONTINUE
            M2=IUP(M1,MU)
C     
            DO 11 NU=MU+1,4
               IPLAQ=6-NU-MU+5*(NU/4)
C     
               DO 7  IJ=1,NCOL2
                  B11(IJ)=U11(IJ,M2,NU)
 7             CONTINUE
               CALL VMX(1,A11,B11,C11,1)
               M3=IUP(M1,NU)
               DO 8  IJ=1,NCOL2
                  D11(IJ)=U11(IJ,M3,MU)
 8             CONTINUE
               CALL HERM(1,D11,DUM11,1)
               CALL VMX(1,C11,D11,B11,1)
               DO 9  IJ=1,NCOL2
                  D11(IJ)=U11(IJ,M1,NU)
 9             CONTINUE
               CALL HERM(1,D11,DUM11,1)
               CALL TRVMX(1,B11,D11,ACT,1)
               ANN=1.0/NCOL
               ACT=ANN*REAL(ACT)
               ACTN(JBIN,IPLAQ)=ACTN(JBIN,IPLAQ)+ACT
               IF(NU.NE.4)THEN
                  VACTS=VACTS+ACT
               ELSE
                  VACTT=VACTT+ACT
               ENDIF
C     
 11         CONTINUE
         ENDDO
      ENDDO
C
      VACTS=VACTS/(3.0*LSIZE)
      VACTT=VACTT/(3.0*LSIZE)
      TOTACT = 0.5*(VACTS+VACTT)
      write (91,*) VACTS, VACTT, TOTACT 
C
      AVACS(JBIN)=AVACS(JBIN)+VACTS
      AVACT(JBIN)=AVACT(JBIN)+VACTT
C
      IF(ITER.EQ.NITER)THEN
      DO IP=1,6
         DO JB=1,NUMBIN
            ACTN(JB,IP)=ACTN(JB,IP)/(IBIN*LSIZE)
         ENDDO
      ENDDO
      ENDIF
C
      ENDIF
C
      IF(IPR.EQ.1)THEN
C
      WRITE(6,80)
80    FORMAT(' *****************************************************')
      WRITE(6,100)
100   FORMAT('                 ACTION            ')
      WRITE(6,80)
      WRITE(6,81)
81    FORMAT(' *  ')
      DO 20 NB=1,NUMBIN
      AV(NB)=AVACS(NB)/IBIN
20    CONTINUE
      CALL JACKM(NUMBIN,AV,AVRS,ERS)
      DO 22 NB=1,NUMBIN
      AV(NB)=AVACT(NB)/IBIN
22    CONTINUE
      CALL JACKM(NUMBIN,AV,AVRT,ERT)
      WRITE(6,110) AVRS,ERS
110   FORMAT('   AVER,ERR   SPACE ACTION =',2F12.8)
      WRITE(6,112) AVRT,ERT
112   FORMAT('   AVER,ERR   TIME  ACTION =',2F12.8)
      WRITE(6,81)
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
C***********************************************************************
C                          THERMAL LINES
C***********************************************************************
C***********************************************************************
      SUBROUTINE POLY(IPR,ITER)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(NITER=2,NUMBIN=2)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/POLYT/TLINE(NITER),SQTLINE(NITER)
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMMON/COMMUNICATE2/AVAC(NUMBIN)
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2)
      DIMENSION AV1(NUMBIN),AV2(NUMBIN)
C
      COMPLEX U11,A11,B11,C11,AKT1,TLINE,AVAC,CSUM
C
      IBIN=NITER/NUMBIN
      IF(IPR.EQ.1)GOTO100
C
      JBIN=(ITER-1)/IBIN+1
C
      IF(ITER.EQ.1)THEN
      DO 3 NB=1,NUMBIN
         AVAC(NB)=(0.0,0.0)
3     CONTINUE
      ENDIF
C
      LX12=LX1*LX2
      LX123=LX1*LX2*LX3
C
      AKT1=(0.0,0.0)
      AKT2=0.0
      DO N4=1,LX4
         DO N3=1,LX3
            DO N2=1,LX2
               M2=1+LX1*(N2-1)+LX12*(N3-1)+LX123*(N4-1)
C
            DO 2 IC=1,NCOL2
               A11(IC)=(0.0,0.0)
2           CONTINUE
            DO 4 NC=1,NCOL
               IC=NC+NCOL*(NC-1)
               A11(IC)=(1.0,0.0)
4           CONTINUE
C
        DO 10 N1=1,LX1
           DO 6 IC=1,NCOL2
              B11(IC)=U11(IC,M2,1)
6          CONTINUE
           CALL VMX(1,A11,B11,C11,1)
           M3=IUP(M2,1)
           M2=M3
           DO 8 IC=1,NCOL2
              A11(IC)=C11(IC)
8          CONTINUE
10      CONTINUE
C
        CSUM=(0.0,0.0)
        DO 12 NC=1,NCOL
           IC=NC+NCOL*(NC-1)
           CSUM=CSUM+A11(IC)
12      CONTINUE
        CSUM=CSUM/NCOL
        AKT1=AKT1+CSUM
        AKT2=AKT2+REALPART(CSUM*CONJG(CSUM))
C
      ENDDO
      ENDDO
      ENDDO
C
      AKT1=AKT1/(LX2*LX3)
      AKT2=AKT2/(LX2*LX3)
      TLINE(ITER)=AKT1
      SQTLINE(ITER)=AKT2
      AVAC(JBIN)=AVAC(JBIN)+AKT1
C
      GOTO200
100   CONTINUE
C
      WRITE(6,80)
80    FORMAT(' *****************************************************')
      WRITE(6,108)
108   FORMAT('     THERMAL LINES           ')
      WRITE(6,80)
      WRITE(6,81)
81    FORMAT(' *  ')
      DO 20 NB=1,NUMBIN
         AV1(NB)=REALPART(AVAC(NB))/IBIN
         AV2(NB)=IMAGPART(AVAC(NB))/IBIN
         WRITE(6,109) NB,AV1(NB),AV2(NB)
109      FORMAT(' BIN =',I4,'   REAL,IMAG  POLY =',2F9.4)
20    CONTINUE
      CALL JACKM(NUMBIN,AV1,AVR,ERR)
      CALL JACKM(NUMBIN,AV2,AVI,ERI)
      WRITE(6,81)
      WRITE(6,110) AVR,ERR
110   FORMAT('   AV,ER  REAL POLY =',2F12.8)
      WRITE(6,111) AVI,ERI
111   FORMAT('   AV,ER  IMAG POLY =',2F12.8)
      WRITE(6,81)
C
200   RETURN
      END
C**********************************************************
C  THIS ROUTINE REIMPOSES THE UNITARITY CONSTRAINTS ON
C  OUR SU(3) MATRICES
C**********************************************************
      SUBROUTINE RENORM
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL,NCOL,LSIZE,4)
      COMPLEX ADUM(NCOL,NCOL),U11,CSUM
C
      COMPLEX aa(ncol,ncol)
      complex det,cnorm,cdet,dncol
C
      DO MU=1,4
         DO NN=1,LSIZE
C
            DO J=1,NCOL
               DO I=1,NCOL
                  ADUM(I,J)=U11(I,J,NN,MU)
               ENDDO
            ENDDO
C
            DO 20 N2=1,NCOL
C
               DO 10 N3=1,N2-1
C
                  CSUM=(0.0,0.0)
                  DO 5 N1=1,NCOL
                     CSUM=CSUM+ADUM(N1,N2)*CONJG(ADUM(N1,N3))
 5                CONTINUE
                  DO 6 N1=1,NCOL
                     ADUM(N1,N2)=ADUM(N1,N2)-CSUM*ADUM(N1,N3)
 6                CONTINUE
C     
 10            CONTINUE
C
               SUM=0.0
               DO 7 N1=1,NCOL
                  SUM=SUM+ADUM(N1,N2)*CONJG(ADUM(N1,N2))
 7             CONTINUE
                  ANORM=1.0/SQRT(SUM)
                  DO 8 N1=1,NCOL
                     ADUM(N1,N2)=ADUM(N1,N2)*ANORM
8                 CONTINUE
C
20          CONTINUE
C
            DO J=1,NCOL
               DO I=1,NCOL
                  U11(I,J,NN,MU)=ADUM(I,J)
               ENDDO
            ENDDO
C
            do j=1,ncol
               do i=1,ncol
                  aa(i,j)=adum(i,j)
               enddo
            enddo
            jmat=ncol
            CALL DETNANT(jmat,det,aa)
c
            dncol=cmplx(1.0/ncol)
            cdet=det**dncol
            cnorm=1.0/cdet
            do j=1,ncol
               do i=1,ncol
                  aa(i,j)=adum(i,j)*cnorm
                  u11(i,j,nn,mu)=aa(i,j)
               enddo
            enddo
C
         ENDDO
      ENDDO
      RETURN
      END
c********************************************************
c********************************************************
      SUBROUTINE DETNANT(NUMOP,DET,AA)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      complex AA(NUMOP,NUMOP),B(NUMOP,NUMOP),det
C
      NP=NUMOP
C
      DO I=1,NP
         DO J=1,NP
            B(I,J)=AA(I,J)
         ENDDO
      ENDDO
      CALL DETMAT(B,NP,DET)
C     
      RETURN
      END
C***************************************************************
      SUBROUTINE DETMAT(A,NP,DET)
      IMPLICIT REAL*8 (A-H,O-Z)
C     
      DIMENSION INDX(NP)
      COMPLEX A(NP,NP),DET,DD
C     
      N=NP
      CALL LUDCMP(A,N,NP,INDX,D)
      DD=D
      DO 13 J=1,NP
         DD=DD*A(J,J)
 13   CONTINUE
      DET=DD
C     
      RETURN
      END
C***************************************************************
C***************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2)
      PARAMETER(NMAX=NCOL,TINY=1.0E-20)
C     
      complex A(NP,NP),SUM,CDUM
      DIMENSION INDX(N),VV(NMAX)
C     
      D=1.
      DO 12 I=1,N
         AAMAX=0.
         DO 11 J=1,N
            IF(CABS(A(I,J)).GT.AAMAX) AAMAX=CABS(A(I,J))
 11      CONTINUE
C     IF(AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX IN LUDCMP'
         IF(AAMAX.EQ.0.) return
         VV(I)=1./AAMAX
 12   CONTINUE
C     
      DO 19 J=1,N
         DO 14 I=1,J-1
            SUM=A(I,J)
            DO 13 K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
 13         CONTINUE
            A(I,J)=SUM
 14      CONTINUE
C     
         AAMAX=0.
         DO 16 I=J,N
            SUM=A(I,J)
            DO 15 K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
 15         CONTINUE
            A(I,J)=SUM
            DUM=VV(I)*CABS(SUM)
            IF(DUM.GE.AAMAX)THEN
               IMAX=I
               AAMAX=DUM        
            ENDIF
 16      CONTINUE
C     
         IF(J.NE.IMAX)THEN
            DO 17 K=1,N
               CDUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=CDUM
 17         CONTINUE
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
C     
         INDX(J)=IMAX
         IF(CABS(A(J,J)).EQ.0.) A(J,J)=TINY
         IF(J.NE.N)THEN
            CDUM=1./A(J,J)
            DO 18 I=J+1,N
               A(I,J)=A(I,J)*CDUM
 18         CONTINUE
         ENDIF
C     
 19   CONTINUE
C     
      RETURN
      END
C****************************************************************
C****************************************************************
C     HERE WE SET UP  LINK POINTERS FOR FULL LATTICE            *
C****************************************************************
C****************************************************************
      SUBROUTINE SETUP
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
C
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
C
      LX12=LX1*LX2
      LX123=LX12*LX3
      LX1234=LX123*LX4
C
      NN=0
      DO L4=1,LX4
         DO L3=1,LX3
            DO L2=1,LX2
               DO L1=1,LX1
                  NN=NN+1
C     
                  NU1=NN+1
                  IF((L1+1).GT.LX1) NU1=NU1-LX1
                  IUP(NN,1)=NU1
                  ND1=NN-1
                  IF((L1-1).LT.1) ND1=ND1+LX1
                  IDN(NN,1)=ND1
C     
                  NU2=NN+LX1
                  IF((L2+1).GT.LX2) NU2=NU2-LX12
                  IUP(NN,2)=NU2
                  ND2=NN-LX1
                  IF((L2-1).LT.1) ND2=ND2+LX12
                  IDN(NN,2)=ND2
C     
                  NU3=NN+LX12
                  IF((L3+1).GT.LX3) NU3=NU3-LX123
                  IUP(NN,3)=NU3
                  ND3=NN-LX12
                  IF((L3-1).LT.1) ND3=ND3+LX123
                  IDN(NN,3)=ND3
C     
                  NU4=NN+LX123
                  IF((L4+1).GT.LX4) NU4=NU4-LX1234
                  IUP(NN,4)=NU4
                  ND4=NN-LX123
                  IF((L4-1).LT.1) ND4=ND4+LX1234
                  IDN(NN,4)=ND4
C     
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C***********************************************************************
C                  jack-knife errors for av(i=1,..,num)                *
C***********************************************************************
      SUBROUTINE JACKM(NUM,AV,AVER,ERR)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AV(*)
C
      SUM=0.0d0
      DO 10 N=1,NUM
         SUM=SUM+AV(N)
10    CONTINUE
      SUM=SUM/NUM
C
      ESUM=0.0d0
      DO 12 N=1,NUM
         DIF=SUM-AV(N)
         ESUM=ESUM+DIF*DIF
12    CONTINUE
      ESUM=ESUM/NUM
C
      AVER=SUM
      ERR=SQRT(ESUM*NUM)/(NUM-1)
C
      RETURN
      END
C***********************************************************************
C                 jack-knife errors for ratio avu to avd               *
C***********************************************************************
      SUBROUTINE JACKMR(NUM,NMAX,AVU,AVD,AVER,ERR)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AVU(NMAX),AVD(NMAX)
      DIMENSION DIFU(NMAX),DIFD(NMAX)
C
      AVER=0.0d0
      ERR=0.0d0
C
      SUMU=0.0d0
      SUMD=0.0d0
      DO 10 N=1,NUM
         SUMU=SUMU+AVU(N)
         SUMD=SUMD+AVD(N)
10    CONTINUE
      DO 12 N=1,NUM
         DIFU(N)=SUMU-AVU(N)
         DIFD(N)=SUMD-AVD(N)
12    CONTINUE
      DO 14 N=1,NUM
         IF(DIFD(N).EQ.0.0d0)THEN
         ERR=99.0d0
         GOTO99
         ENDIF
14    CONTINUE
C
      ASUM=0.0d0
      ESUM=0.0d0
      DO 16 N=1,NUM
         DIF=DIFU(N)/DIFD(N)
         ASUM=ASUM+DIF
         ESUM=ESUM+DIF*DIF
16    CONTINUE
      ASUM=ASUM/NUM
      ESUM=ESUM/NUM
      ESUM=ESUM-ASUM*ASUM
C
      IF(SUMD.NE.0.0) AVER=SUMU/SUMD
      IF(ESUM.GT.0.0) ERR=SQRT(ESUM*NUM)
C
99    CONTINUE
      RETURN
      END
C***********************************************************************
C                           fitting masses                             *
C***********************************************************************
      SUBROUTINE JACK(NBIN,ISUB,LXI,NUMBIN,LMAX,COR,VAC)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
      DIMENSION COR(NUMBIN,LMAX),VAC(NUMBIN)
      DIMENSION CORR(200),AW(200)
C
      DO 1 I4=1,LMAX
         ACOR(I4)=0.0d0
         SCOR(I4)=0.0d0
         AMM(I4)=0.0d0
         SMM(I4)=0.0d0
 1    CONTINUE
C
      DO 2 IB=1,NBIN
C
         IF(ISUB.EQ.0)THEN
            VACC=0.0d0
            DO 4 JB=1,NBIN
               IF(JB.EQ.IB)GOTO4
               VACC=VACC+VAC(JB)
 4          CONTINUE
         ENDIF
c
         DO 6 NT=1,LMAX
            CORR(NT)=0.0d0
            DO 7 JB=1,NBIN
               IF(JB.EQ.IB)GOTO7
               CORR(NT)=CORR(NT)+COR(JB,NT)
 7          CONTINUE
            IF(ISUB.EQ.0) CORR(NT)=CORR(NT)-VACC*VACC/(NBIN-1)
 6       CONTINUE
         ANORM=CORR(1)
         IF(ANORM.EQ.0.0)GOTO99
         DO 8 NT=1,LMAX
            CORR(NT)=CORR(NT)/ANORM
 8       CONTINUE
         CALL FITT(LXI,CORR,AW,NTMAX)
         DO 9 NT=1,LMAX
            ACOR(NT)=ACOR(NT)+CORR(NT)
            SCOR(NT)=SCOR(NT)+CORR(NT)**2
            AMM(NT)=AMM(NT)+AW(NT)
            SMM(NT)=SMM(NT)+AW(NT)**2
 9       CONTINUE
 2    CONTINUE
C  C
      DO 10 NT=1,LMAX
         ACOR(NT)=ACOR(NT)/NBIN
         SCOR(NT)=SCOR(NT)/NBIN
         SCOR(NT)=(SCOR(NT)-ACOR(NT)*ACOR(NT))*NBIN
         IF(SCOR(NT).GT.0.0) SCOR(NT)=SQRT(SCOR(NT))
         AMM(NT)=AMM(NT)/NBIN
         SMM(NT)=SMM(NT)/NBIN
         SMM(NT)=(SMM(NT)-AMM(NT)*AMM(NT))*NBIN
         IF(SMM(NT).GT.0.0) SMM(NT)=SQRT(SMM(NT))
 10   CONTINUE
C     
      IF(ISUB.EQ.0)THEN
         VACC=0.0d0
         DO 22 JB=1,NBIN
            VACC=VACC+VAC(JB)
 22      CONTINUE
      ENDIF
      DO 24 NT=1,LMAX
         CORR(NT)=0.0d0
         DO 26 JB=1,NBIN
            CORR(NT)=CORR(NT)+COR(JB,NT)
 26      CONTINUE
         IF(ISUB.EQ.0) CORR(NT)=CORR(NT)-VACC*VACC/NBIN
 24   CONTINUE
      ANORM=CORR(1)
      IF(ANORM.EQ.0.0)GOTO99
      DO 28 NT=1,LMAX
         CORR(NT)=CORR(NT)/ANORM
 28   CONTINUE
      CALL FITT(LXI,CORR,AW,NTMAX)
      DO 29 NT=1,LMAX
         ACOR(NT)=CORR(NT)
         AMM(NT)=AW(NT)
 29   CONTINUE
C     C
 99   RETURN
      END
C*********************************************************************
C effective energies using a local cosh fit - lattice length LT.
C AW(timedif+1) is input corrln function;
C BW(nt) is eff energy from  time differences nt-1 to nt.
C*********************************************************************
      SUBROUTINE FITT(LT,AW,BW,NTMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AW(200),BW(200)
C
      DO 1 I=1,LT/2+1
         BW(I)=0.0d0
1     CONTINUE
C
      NTMAX=LT/2+1
      DO 2 I4=1,LT/2
C  C
         IF(AW(I4).LE.0.000001.OR.AW(I4+1).LE.0.000001)THEN
         NTMAX=I4-1
         GOTO99
         ENDIF
         AMSM=0.0d0
         FTM=(AW(I4)/AW(I4+1))
         IF(FTM.GT.1.0)THEN
         AML=DLOG(FTM)
         AMU=DLOG(2.0d0*FTM)
         DO 13 NS=1,20
            AMS=(AML+AMU)/2
            FTS=(EXP(-AMS*(I4-1))+EXP(-(LT-I4+1)*AMS))
     &           /(EXP(-AMS*(I4))+EXP(-(LT-I4)*AMS))
            IF(FTS.LT.FTM)THEN
               AML=AMS
            ELSE
               AMU=AMS
            ENDIF
 13      CONTINUE
         AMSM=(AML+AMU)/2
      ELSE
         NTMAX=I4-1
         GOTO99
      ENDIF
      BW(I4)=AMSM
C     
 2    CONTINUE
C
 99   CONTINUE
      RETURN
      END
C***********************************************************************
C                 VECTOR MATRIX MULTIPLY ... 5*5 COMPLEX               *
C***********************************************************************
      SUBROUTINE VMX(NNN1,A,B,C,NNN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMPLEX A(NCOL,NCOL),B(NCOL,NCOL),C(NCOL,NCOL),CSUM
C
      DO J=1,NCOL
         DO I=1,NCOL
            CSUM=(0.0,0.0)
            DO 2 K=1,NCOL
               CSUM=CSUM+A(I,K)*B(K,J)
2           CONTINUE
            C(I,J)=CSUM
         ENDDO
      ENDDO
C
      RETURN
      END
C***********************************************************************
C***********************************************************************
C                   TRACE PRODUCT .. Ncol*Ncol COMPLEX                 *
C***********************************************************************
C***********************************************************************
      SUBROUTINE TRVMX(NNN1,A,B,CC,NNN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMPLEX A(NCOL,NCOL),B(NCOL,NCOL),CC
C
      CC=(0.0,0.0)
      DO I=1,NCOL
         DO K=1,NCOL
            CC=CC+A(I,K)*B(K,I)
         ENDDO
      ENDDO
C
      RETURN
      END
C***********************************************************************
C***********************************************************************
C                           HERMITIAN CONJUGATE                        *
C***********************************************************************
C***********************************************************************
      SUBROUTINE HERM(NNN1,A11,DUM11,NNN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMPLEX A11(NCOL,NCOL),DUM11(NCOL,NCOL)
C
      DO I=1,NCOL
         DO J=1,NCOL
            DUM11(I,J)=CONJG(A11(J,I))
         ENDDO
      ENDDO
      DO I=1,NCOL
         DO J=1,NCOL
            A11(I,J)=DUM11(I,J)
         ENDDO
      ENDDO
C
      RETURN
      END
