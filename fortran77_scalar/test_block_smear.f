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
      PARAMETER(IBLOK=5,IBING=109,NUMBIN=2)
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
      homepath='/gpfs/scratch/ehpc598/torelons/'
      conf_directory1='cnfg/'
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
      CALL MEASURE(ITER,NITER)
      call execute_command_line("rm conf", WAIT=.true.)
C
 20   CONTINUE
C
      STOP
      END
C*********************************************************************
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
C***************************SUBROUTINE MEASURE**********************
C*******************************************************************
      SUBROUTINE MEASURE(ITER,NTOT)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ICALLG=1,IBING=109,NUMBIN=2)
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

C **** SAVE THE SMEARED FUNCTION
      OPEN (11,FILE='SMEARED_SAVE.DAT',FORM='UNFORMATTED',
     &status='REPLACE',ACCESS='STREAM')

      WRITE(11) UC11

      CLOSE(11)
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
C *** SAVE BLOCKED 
      OPEN (11,FILE='BLOCKED_SAVE.DAT',FORM='UNFORMATTED',
     &status='REPLACE',ACCESS='STREAM')

      WRITE(11) UB11

      CLOSE(11)
11       CONTINUE
C
C
C************************************************
 3       CONTINUE
 1    CONTINUE
      CALL cpu_time(t2)
      WRITE(6,901) sngl(t2-t1)
 901  FORMAT('[Info][Time]','     Thermal Line time:', F8.4)      
C
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
C***********************************************************************
C                          THERMAL LINES
C***********************************************************************
C***********************************************************************
      SUBROUTINE POLY(IPR,ITER)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
      PARAMETER(NITER=218,NUMBIN=2)
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
