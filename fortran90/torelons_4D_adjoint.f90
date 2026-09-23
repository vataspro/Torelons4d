!***************************************************************
!***************************************************************
!     ******  ******  ***    ******  **      ******  *    *
!       **    *    *  *  *   *       **      *    *  **   *
!       **    *    *  * *    ******  **      *    *  * *  *
!       **    *    *  *  *   *       **      *    *  *  * *
!       **    ******  *   *  ******  ******  ******  *   **
!***************************************************************
!***************************************************************
!     Code for calculation the correlation functions of torelons
!     So far the code includes operators for k=1 N-ality
!     In the future there will be k=2 N-ality operatorsx
!***************************************************************
!***************************************************************
module torelons_binary_io
  use iso_fortran_env, only : int8, int32, int64, real64
  implicit none
  private
  public :: read_be_int32, read_be_real64

contains

  function read_be_int32(unit) result(value)
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

  function read_be_real64(unit) result(value)
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

end module torelons_binary_io

program torelons_4d_adjoint
use iso_fortran_env, only : real32, real64, int32
use luxury, only : rluxgo
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1)
!
PARAMETER(ICALLG=1,NITER=2,NMEASL=NITER/ICALLG)
PARAMETER(ICMIN=1,ICMAX=NITER+ICMIN-1)
PARAMETER(IBLOK=5,IBING=1,NUMBIN=2)
PARAMETER(PARBS=0.30,PARBDS=0.12)
!
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!******************************************************************************
!******************************************************************************
!      
COMMON/ARRAYS/U11(NCOL2*LSIZE*4)
complex(real32) U11
real(real64) rndnum
!
character(len=43) :: homepath
character(len=24) :: conf_directory1
character(len=35) :: conf_directory2
character(len=102) :: directory
character(len=120) :: file_name
character(len=256) :: copy_file
character(len=256) :: trnsf_file
character(len=5) :: confnum
character(len=6) :: name_of_file
!                    
homepath=trim('/home/dp208/dp208/dc-athe1/AXIONS/NF2/b2.3/')
conf_directory1=trim('m-1.0/26x26x26x52/confs/')
conf_directory2=trim('run1_52x26x26x26nc2rADJnf2b2.300000')
directory=trim(homepath//conf_directory1//conf_directory2)
!      
CALL SETUP
ISEED=3591
call rluxgo(3, iseed, 0, 0)
!
ITOT=NITER
!
WRITE(*,*) "                                                "
WRITE(6,90)
90 FORMAT(' *******************************************************')
WRITE(*,*) "                                                "
WRITE(*,*) &
  & " ******  ******  ***    ******  **      ******  *    *"
WRITE(*,*) &
  & "   **    *    *  *  *   *       **      *    *  **   *"
WRITE(*,*) &
  & "   **    *    *  * *    ******  **      *    *  * *  *"
WRITE(*,*) &
  & "   **    *    *  *  *   *       **      *    *  *  * *"
WRITE(*,*) &
  & "   **    ******  *   *  ******  ******  ******  *   **"
WRITE(*,*) "                                                "
WRITE(6,90)
WRITE(6,91)
91 FORMAT(' *')
WRITE(6,92) LX1,LX2,LX3,LX4
92 FORMAT("[Info][Lattice Size]  ",'             LX   = ',4I4)
WRITE(6,193) BETAG
193 FORMAT("[Info][Value of Beta] ",'           beta   =',F8.4)
WRITE(6,194) NCOL
194 FORMAT("[Info][Colors]        ",' Number of Colors = ',I6)
WRITE(6,195) IHEAT
195 FORMAT("[Info][Thermalisation]",'   Therml. sweeps = ',I6)
WRITE(6,91)
WRITE(6,90)
WRITE(6,91)
WRITE(6,196) ITOT
196 FORMAT("[Info][Measurements]",'   Number of MC Iterations = ',I6)
WRITE(6,197) ICALLG
197 FORMAT("[Info][Measurements]",'   Sweeps per measurements = ',I6)
WRITE(6,198) IBING
198 FORMAT("[Info][Measurements]",'      Measurements per bin = ',I6)
WRITE(6,91)
WRITE(6,90)
200 FORMAT("[Info][Number of Operators q=0]",'  Total Number = ',I4)
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
205 FORMAT("[Info][Number of Operators q=0]",'     J   P   R      #')
206 FORMAT("[Info][Number of Operators q=0]",'     ----------------')
210 FORMAT("[Info][Number of Operators q=0]",'     0   +   +   ',I4)
211 FORMAT("[Info][Number of Operators q=0]",'     0   +   -   ',I4)
212 FORMAT("[Info][Number of Operators q=0]",'     0   -   +   ',I4)
213 FORMAT("[Info][Number of Operators q=0]",'     0   -   -   ',I4)
214 FORMAT("[Info][Number of Operators q=0]",'     1   #   +   ',I4)
216 FORMAT("[Info][Number of Operators q=0]",'     1   #   -   ',I4)
218 FORMAT("[Info][Number of Operators q=0]",'     2   +   +   ',I4)
219 FORMAT("[Info][Number of Operators q=0]",'     2   +   -   ',I4)
220 FORMAT("[Info][Number of Operators q=0]",'     2   -   +   ',I4)
221 FORMAT("[Info][Number of Operators q=0]",'     2   -   -   ',I4)
207 FORMAT("[Info][Number of Operators q=0]",'     ----------------')
WRITE(6,91)
WRITE(6,90)
300 FORMAT("[Info][Number of Operators q=1,2]",'  Total Number = ',I4)
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
305 FORMAT("[Info][Number of Operators q=1,2]",'     J   P   #')
306 FORMAT("[Info][Number of Operators q=1,2]",'     ---------------')
310 FORMAT("[Info][Number of Operators q=1,2]",'     0   + ',I4)
312 FORMAT("[Info][Number of Operators q=1,2]",'     0   - ',I4)
314 FORMAT("[Info][Number of Operators q=1,2]",'     1     ',I4)
318 FORMAT("[Info][Number of Operators q=1,2]",'     2   + ',I4)
320 FORMAT("[Info][Number of Operators q=1,2]",'     2   - ',I4)
307 FORMAT("[Info][Number of Operators q=1,2]",'     ---------------')
WRITE(6,91)
WRITE(6,90)
!
IFILE = ICMIN - 1
do ITER=1,NITER
IFILE = IFILE + 1
!
ist=IFILE
if (ist < 10) then
write(confnum,'(i1)') ist
file_name=trim(directory//'m1.000000n'//trim(confnum))
else if (ist < 100) then
write(confnum,'(i2)') ist
file_name=trim(directory//'m1.000000n'//trim(confnum))
else if (ist < 1000) then
write(confnum,'(i3)') ist
file_name=trim(directory//'m1.000000n'//trim(confnum))
else if (ist < 10000) then
write(confnum,'(i4)') ist
file_name=trim(directory//'m1.000000n'//trim(confnum))
else if (ist < 100000) then
write(confnum,'(i5)') ist
file_name=trim(directory//'m1.000000n'//trim(confnum))
end if
!     
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
!         
CALL cpu_time(t1)
CALL READ_GF(name_of_file)
CALL cpu_time(t2)
WRITE(6,902) real(t2-t1, kind=real32)
902 FORMAT('[Info][Time]','    Configuration Read time:', F8.4)
!     
CALL ACTION(0,ITER,TOTACT)
511 format('[FM][0]Check plaq = ', f8.6)
write (6,511) TOTACT
CALL POLYLOOP()
CALL MEASURE(ITER,NITER)
call execute_command_line("rm conf", WAIT=.true.)
!
20 CONTINUE
end do
!
CALL TODISK1
CALL TODISK2
!
STOP
end program torelons_4d_adjoint
!*********************************************************************
function quaternion_to_matrix(quaternion)
use iso_fortran_env, only : real32, real64, int32
real(real64), dimension(4) :: quaternion
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
!*********************************************************************
SUBROUTINE READ_GF(filename)
use iso_fortran_env, only : real32, real64, int32
use torelons_binary_io, only : read_be_int32, read_be_real64
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(NCOL=2,ND=4)
!
character(len=6) :: filename
integer(int32) nc_read, nx_read, ny_read, nz_read, nt_read
integer :: t, x, y, z, dir, dir_target, iun, iq
real(real64) :: plaquette_read

real(real64) :: quaternion(4)
complex(real32), dimension(2, 2) :: quaternion_to_matrix

COMMON/ARRAYS/U11(NCOL,NCOL,LX1,LX2,LX3,LX4,4)

complex(real32) U11

open(newunit=iun, file=filename, access='stream', form='unformatted', &
     status='old', action='read')
nc_read      = read_be_int32(iun)
nt_read      = read_be_int32(iun)
nx_read      = read_be_int32(iun)
ny_read      = read_be_int32(iun)
nz_read      = read_be_int32(iun)
plaquette_read = read_be_real64(iun)

WRITE(6,101) plaquette_read
101 FORMAT('[I/O][Plaq]', 'Plaquette value:',F8.6)
WRITE(6,102) nc_read
102 FORMAT('[I/O][Ncol]', 'Number of Colors:',I2.1)
WRITE(6,103) nt_read, nx_read, ny_read, nz_read
103 FORMAT('[I/O][Dim]', 'T x X x Y x Z=',I3.2, I3.2, I3.2, I3.2)
!
DO t = 1, LX4
DO x = 1, LX1
DO y = 1, LX2
DO z = 1, LX3
DO dir = 1, ND
do iq = 1, size(quaternion)
  quaternion(iq) = read_be_real64(iun)
end do

if (dir == 1) then
dir_target = 4
else
dir_target = dir-1
end if

U11(1, 1, x, y, z, t, dir_target) = &
  & cmplx(real(quaternion(1), kind=real32),real(quaternion(4), kind=real32))
U11(1, 2, x, y, z, t, dir_target) = &
  & cmplx(-real(quaternion(3), kind=real32),real(quaternion(2), kind=real32))
U11(2, 1, x, y, z, t, dir_target) = &
  & cmplx(real(quaternion(3), kind=real32),real(quaternion(2), kind=real32))
U11(2, 2, x, y, z, t, dir_target) = &
  & cmplx(real(quaternion(1), kind=real32),-real(quaternion(4), kind=real32))
!
end do
end do
end do
end do
end do

close(iun)
end subroutine READ_GF
!*********************************************************************
!*********************************************************************
!*********************************************************************
SUBROUTINE READ_GAUGE_ETMC(filename) ! to be fixed !
use torelons_binary_io, only : read_be_int32, read_be_real64
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)

COMMON/ARRAYS/U11(NCOL,NCOL,LX1,LX2,LX3,LX4,4)

complex(real32) U11
complex(real64) UR11(LX4,LX3,LX2,LX1,4,NCOL,NCOL)
character(len=6) :: filename
integer :: iun
real(real64) :: rpart, ipart
!
open(newunit=iun, file=trim(filename), access='stream', &
     form='unformatted', status='old', action='read')
!
!      REWIND(60)
!
ICALL=1
DO L4=1, LX4
DO L3=1, LX3
DO L2=1, LX2
DO L1=1, LX1
DO MU=1, 4
DO IJ=1, NCOL
DO IK=1, NCOL
!     C
rpart = read_be_real64(iun)
ipart = read_be_real64(iun)
UR11(L4,L3,L2,L1,MU,IJ,IK) = cmplx(rpart, ipart, kind=real64)
!     C
ICALL=ICALL+1
end do
end do
end do
end do
end do
end do
end do
!
DO IJ=1, NCOL
DO IK=1, NCOL
DO L1=1, LX1
DO L2=1, LX2
DO L3=1, LX3
DO L4=1, LX4
DO MU=1, 4
!     
DREALCONF=real(UR11(L4,L3,L2,L1,MU,IJ,IK))
DIMAGCONF=aimag(UR11(L4,L3,L2,L1,MU,IJ,IK))
SREALCONF=real(DREALCONF, kind=real32)
SIMAGCONF=real(DIMAGCONF, kind=real32)
U11(IJ,IK,L1,L2,L3,L4,MU)=cmplx(SREALCONF,SIMAGCONF)
!
!        write(*,*) U11(IJ,IK,L1,L2,L3,L4,MU)
!                                                                                                                                                                        
end do
end do
end do
end do
end do
end do
end do
!
WRITE(6,*) 'Configuration read'
close(iun)
!
RETURN
end subroutine READ_GAUGE_ETMC
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                          THERMAL LINES
!***********************************************************************
!***********************************************************************
SUBROUTINE POLYLOOP()
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(NITER=2,NUMBIN=2)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
DIMENSION DUM11(NCOL2) &
  & ,A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2)
complex(real32) U11,A11,B11,C11,D11,DUM11,ACT
!
DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN),AVACS(NUMBIN),AVACT(NUMBIN)
DIMENSION VAL(NUMBIN),AV(NUMBIN)
dimension icoordvect(4)
complex(real32) cpol1
complex(real64) cpol(4)
!
icoordvect(1) = LX1-1
icoordvect(2) = LX2-1
icoordvect(3) = LX3-1
icoordvect(4) = LX4-1
!     
cpol(:) = cmplx(0.0_real64, 0.0_real64, kind=real64)
dnorm = 1.0/real(lsizeb*ncol, kind=real64)
!
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
end do
ipoint2 = iup(ipoint1,idir)
b11 = u11(:,ipoint2,idir)
call trvmx(1,a11,b11,cpol1,1)
cpol(idir) = cpol(idir) + cpol1
!            write (*,*) cpol1, cpol(idir)
end do
end do
!
cpol(:) = cpol(:)*dnorm
!
888 format('Poly(' ,i1, ')  =  (',f16.12,',',f16.12,')')
do idir=1,4
write (92,888) idir,real(cpol(idir)),aimag(cpol(idir))
end do
!
RETURN
end subroutine POLYLOOP
!*********************************************************************
!***************************SUBROUTINE TODISK*************************
!*********************************************************************
SUBROUTINE TODISK1
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(IBLOK=5,NUMBIN=2,NITER=2)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LMAX=LX4/2+1)
PARAMETER(LMAXIR=3)
!******************************************************************************
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************  
COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLP(NUMBIN,NTOTAL)
COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP) &
  & ,AVACLJ0PP(NUMBIN,NOPJ0PP)
COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM) &
  & ,AVACLJ0PM(NUMBIN,NOPJ0PM)
COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP) &
  & ,AVACLJ0MP(NUMBIN,NOPJ0MP)
COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM) &
  & ,AVACLJ0MM(NUMBIN,NOPJ0MM)
COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P) &
  & ,AVACLJ1P(NUMBIN,NOPJ1P)
COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M) &
  & ,AVACLJ1M(NUMBIN,NOPJ1M)
COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP) &
  & ,AVACLJ2PP(NUMBIN,NOPJ2PP)
COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM) &
  & ,AVACLJ2PM(NUMBIN,NOPJ2PM)
COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP) &
  & ,AVACLJ2MP(NUMBIN,NOPJ2MP)
COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM) &
  & ,AVACLJ2MM(NUMBIN,NOPJ2MM)
!
COMMON/PBLOCKJ0/ACORLJ0(NUMBIN,LMAX,NOPFULJ0,NOPFULJ0) &
  & ,AVACLJ0(NUMBIN,NOPFULJ0)
COMMON/PBLOCKJ1/ACORLJ1(NUMBIN,LMAX,NOPFULJ1,NOPFULJ1) &
  & ,AVACLJ1(NUMBIN,NOPFULJ1)
COMMON/PBLOCKJ2/ACORLJ2(NUMBIN,LMAX,NOPFULJ2,NOPFULJ2) &
  & ,AVACLJ2(NUMBIN,NOPFULJ2)
!
COMMON/PBLOCKJALL/ACORLALL(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLALL(NUMBIN,NTOTAL)
!
COMMON/POLYT/TLINE(NITER),SQTLINE(NITER)
COMMON/ASTORE/ACTN(NUMBIN,6),PLAQ(NITER,6)
!
complex(real64) ACORLP
complex(real64) AVACLP
complex(real64) ACORLJ0PP
complex(real64) AVACLJ0PP
complex(real64) ACORLJ0PM
complex(real64) AVACLJ0PM
complex(real64) ACORLJ0MP
complex(real64) AVACLJ0MP
complex(real64) ACORLJ0MM
complex(real64) AVACLJ0MM
complex(real64) ACORLJ1P
complex(real64) AVACLJ1P
complex(real64) ACORLJ1M
complex(real64) AVACLJ1M
complex(real64) ACORLJ2PP
complex(real64) AVACLJ2PP
complex(real64) ACORLJ2PM
complex(real64) AVACLJ2PM
complex(real64) ACORLJ2MP
complex(real64) AVACLJ2MP
complex(real64) ACORLJ2MM
complex(real64) AVACLJ2MM
!
complex(real64) ACORLJ0
complex(real64) AVACLJ0
complex(real64) ACORLJ1
complex(real64) AVACLJ1
complex(real64) ACORLJ2
complex(real64) AVACLJ2
complex(real64) ACORLALL
complex(real64) AVACLALL
!     
complex(real32) TLINE
!
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
!******************************************
DO IN=1, NUMBIN

!
do IL=1,LMAXIR
!     
DO IX=1, NTOTAL
DO IY=1, NTOTAL
WRITE(11,*) ACORLP(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NTOTAL
DO IY=1, NTOTAL
WRITE(25,*) ACORLALL(IN,IL,IY,IX)
end do
end do
!     
31 CONTINUE
end do
!     
DO IL=1,LMAX
!
!****************************************************
!*************  J0 CORRELATORS **********************
!****************************************************
DO IX=1, NOPJ0PP
DO IY=1, NOPJ0PP
WRITE(12,*) ACORLJ0PP(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPJ0PM
DO IY=1, NOPJ0PM
WRITE(13,*) ACORLJ0PM(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPJ0MP
DO IY=1, NOPJ0MP
WRITE(14,*) ACORLJ0MP(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPJ0MM
DO IY=1, NOPJ0MM
WRITE(15,*) ACORLJ0MM(IN,IL,IY,IX)
end do
end do
!****************************************************
!*************  J1 CORRELATORS **********************
!****************************************************
DO IX=1, NOPJ1P
DO IY=1, NOPJ1P
WRITE(16,*) ACORLJ1P(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPJ1M
DO IY=1, NOPJ1M
WRITE(17,*) ACORLJ1M(IN,IL,IY,IX)
end do
end do
!****************************************************
!*************  J2 CORRELATORS **********************
!****************************************************
DO IX=1, NOPJ2PP
DO IY=1, NOPJ2PP
WRITE(18,*) ACORLJ2PP(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPJ2PM
DO IY=1, NOPJ2PM
WRITE(19,*) ACORLJ2PM(IN,IL,IY,IX)
end do
end do
!
DO IX=1, NOPJ2MP
DO IY=1, NOPJ2MP
WRITE(20,*) ACORLJ2MP(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPJ2MM
DO IY=1, NOPJ2MM
WRITE(21,*) ACORLJ2MM(IN,IL,IY,IX)
end do
end do
!****************************************************
!*************  FULL J CORRELATORS ******************
!****************************************************
DO IX=1, NOPFULJ0
DO IY=1, NOPFULJ0
WRITE(22,*) ACORLJ0(IN,IL,IY,IX)
end do
end do
!                  
DO IX=1, NOPFULJ1
DO IY=1, NOPFULJ1
WRITE(23,*) ACORLJ1(IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPFULJ2
DO IY=1, NOPFULJ2
WRITE(24,*) ACORLJ2(IN,IL,IY,IX)
end do
end do
!****************************************************
!****************************************************
!****************************************************               
end do
end do
!******************************************
REWIND(52)
WRITE(52) ACTN
WRITE(52) TLINE,PLAQ,SQTLINE
! Close all output files
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
RETURN
end subroutine TODISK1
!*********************************************************************
!***************************SUBROUTINE TODISK*************************
!*********************************************************************
SUBROUTINE TODISK2
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(IBLOK=5,NUMBIN=2,NITER=2)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LMAX=LX4/2+1)
PARAMETER(LMAXIR=3)
!******************************************************************************
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************  
COMMON/PBLOCKMOMP/ACORLPMOM(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM) &
  & ,AVACLPMOM(2,NUMBIN,NTOTALMOM)
!
COMMON/PBLOCKMOMJ0P/ &
  & ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM) &
  & ,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
COMMON/PBLOCKMOMJ0M/ &
  & ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM) &
  & ,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)
COMMON/PBLOCKMOMJ1/ &
  & ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM) &
  & ,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)
COMMON/PBLOCKMOMJ2P/ &
  & ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM) &
  & ,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
COMMON/PBLOCKMOMJ2M/ &
  & ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM) &
  & ,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
!
COMMON/PBLOCKMOMJ0/ &
  & ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM) &
  & ,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
COMMON/PBLOCKMOMJ2/ &
  & ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM) &
  & ,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)
!
COMMON/PBLOCKMOMJALL/ &
  & ACORLPMOMJALL(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM) &
  & ,AVACLMOMALL(2,NUMBIN,NTOTALMOM)
!*****************************************************************
complex(real64) ACORLPMOM
complex(real64) AVACLPMOM
complex(real64) ACORLPMOMJ0P
complex(real64) AVACLMOMJ0P
complex(real64) ACORLPMOMJ0M
complex(real64) AVACLMOMJ0M
complex(real64) ACORLPMOMJ1
complex(real64) AVACLMOMJ1
complex(real64) ACORLPMOMJ2P
complex(real64) AVACLMOMJ2P
complex(real64) ACORLPMOMJ2M
complex(real64) AVACLMOMJ2M
complex(real64) ACORLPMOMJ0
complex(real64) AVACLMOMJ0
complex(real64) ACORLPMOMJ2
complex(real64) AVACLMOMJ2
complex(real64) ACORLPMOMJALL
complex(real64) AVACLMOMALL
!*****************************************************************
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
!*****************************************************************
DO IN=1,NUMBIN

!
do IL=1,LMAXIR
!            
DO IX=1, NTOTALMOM
DO IY=1, NTOTALMOM
WRITE(17,*) ACORLPMOM(1,IN,IL,IY,IX)
WRITE(18,*) ACORLPMOM(2,IN,IL,IY,IX)
WRITE(35,*) ACORLPMOMJALL(1,IN,IL,IY,IX)
WRITE(36,*) ACORLPMOMJALL(2,IN,IL,IY,IX)
end do
end do
!     
31 CONTINUE
end do

DO IL=1,LMAX
!****************************************************
!*************  J0 CORRELATORS **********************
!****************************************************
DO IX=1,NOPJ0PMOM
DO IY=1, NOPJ0PMOM
WRITE(19,*) ACORLPMOMJ0P(1,IN,IL,IY,IX)
WRITE(20,*) ACORLPMOMJ0P(2,IN,IL,IY,IX)
end do
end do
!     
DO IX=1,NOPJ0MMOM
DO IY=1, NOPJ0MMOM
WRITE(21,*) ACORLPMOMJ0M(1,IN,IL,IY,IX)
WRITE(22,*) ACORLPMOMJ0M(2,IN,IL,IY,IX)
end do
end do
!****************************************************
!*************  J1 CORRELATORS **********************
!****************************************************
DO IX=1,NOPJ1MOM
DO IY=1, NOPJ1MOM
WRITE(23,*) ACORLPMOMJ1(1,IN,IL,IY,IX)
WRITE(24,*) ACORLPMOMJ1(2,IN,IL,IY,IX)
end do
end do
!****************************************************
!*************  J2 CORRELATORS **********************
!****************************************************
DO IX=1,NOPJ2PMOM
DO IY=1, NOPJ2PMOM
WRITE(25,*) ACORLPMOMJ2P(1,IN,IL,IY,IX)
WRITE(26,*) ACORLPMOMJ2P(2,IN,IL,IY,IX)
end do
end do
!     
DO IX=1,NOPJ2MMOM
DO IY=1, NOPJ2MMOM
WRITE(27,*) ACORLPMOMJ2M(1,IN,IL,IY,IX)
WRITE(28,*) ACORLPMOMJ2M(2,IN,IL,IY,IX)
end do
end do
!****************************************************
!*************  FULL J CORRELATORS ******************
!****************************************************
DO IX=1, NOPFULJ0MOM
DO IY=1, NOPFULJ0MOM
WRITE(29,*) ACORLPMOMJ0(1,IN,IL,IY,IX)
WRITE(30,*) ACORLPMOMJ0(2,IN,IL,IY,IX)
end do
end do
!     
DO IX=1, NOPFULJ2MOM
DO IY=1, NOPFULJ2MOM
WRITE(33,*) ACORLPMOMJ2(1,IN,IL,IY,IX)
WRITE(34,*) ACORLPMOMJ2(2,IN,IL,IY,IX)
end do
end do
!****************************************************
end do
end do
!*******************************************************************
RETURN
end subroutine TODISK2
!*******************************************************************
!***************************SUBROUTINE MEASURE**********************
!*******************************************************************
SUBROUTINE MEASURE(ITER,NTOT)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
!
PARAMETER(ICALLG=1,IBING=1,NUMBIN=2)
!
COMMON/ITEM/ITERG,JBING,NTOTG
!
ITERG=ITER/ICALLG
JBING=IBING
NTOTG=IBING*NTOT/(ICALLG*IBING)
!
IF(ITER.EQ.(ITERG*ICALLG).AND.ITERG.LE.NTOTG) THEN
PAR=1.0d0
!
CALL SETUPB
!
!        WRITE(*,*) "CALLING BLOCK"
CALL BLOCK
WRITE(*,*) "Number of iteration:", ITER
end if
!
IF(ITER.EQ.NTOTG*ICALLG) THEN
CALL MLINE1
CALL MLINE2
CALL MLINE3
CALL MLINE4
CALL MLINE5
end if
! Close output files
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
RETURN
end subroutine MEASURE
!*********************************************************************
!*************************SUBROUTINE MLINE1***************************
!*********************************************************************
SUBROUTINE MLINE1
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
PARAMETER(LMAXIR=3)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!******************************************************************************
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------      
!
COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLP(NUMBIN,NTOTAL)
!
COMMON/ITEM/ITR,IBIN,ITOT
COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
!
complex(real64) ACORLP,AVACLP
!      
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE='DIAGNALFULLQ0.DAT')
!
NBIN=ITOT/IBIN
LS(1)=LX1
LS(2)=LX2
LS(3)=LX3
EPS=0.0000000001d0
!********************************************************************
!     OPERATORS
!********************************************************************
WRITE(23,80)
WRITE(23,80)
80 FORMAT('****************************************************')
WRITE(23,181)
181 FORMAT('****************** OPERATORS WITH J = 0 ************')
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(' *')
WRITE(23,91)
91 FORMAT(' AVERAGE  X,Y  LINES ')
WRITE(23,81)
82 FORMAT(' ******************* ')
!
do ID=1, NTOTAL
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
do IB=1,NUMBIN
AV(IB)=real(AVACLP(IB,ID))/ACOUNT
13 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
do IB=1,NUMBIN
AV(IB)=aimag(AVACLP(IB,ID))/ACOUNT
14 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101 FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10 CONTINUE
end do
!
WRITE(23,80)
WRITE(23,81)
WRITE(23,92)
92 FORMAT(' AVERAGE   Y   LINES ')
WRITE(23,81)

do ID=1, NTOTAL
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(' OPERATOR = ',I4)
WRITE(23,82)
DO NT=1,LMAXIR
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLP(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAXIR,CORP,VAC)
!         DO 26 NT=1,MAXDTLS
do NT=1,LMAXIR-1
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
26 CONTINUE
end do
20 CONTINUE
end do
!
CLOSE(23)
RETURN
end subroutine MLINE1
!*********************************************************************
!*************************SUBROUTINE MLINE3***************************
!*********************************************************************
SUBROUTINE MLINE2
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
PARAMETER(LMAXIR=3)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!******************************************************************************
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!
COMMON/PBLOCKMOMP/ACORLPMOM(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM) &
  & ,AVACLPMOM(2,NUMBIN,NTOTALMOM)
!
COMMON/ITEM/ITR,IBIN,ITOT
COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
!
complex(real64) ACORLPMOM,AVACLPMOM
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (24, FILE='DIAGNALFULLMOMQ0.DAT')
!
NBIN=ITOT/IBIN
LS(1)=LX1
LS(2)=LX2
LS(3)=LX3
EPS=0.0000000001d0
!
DO IJ=1,2
!********************************************************************
!     POSITIVE PARITY OPERATORS
!********************************************************************
WRITE(24,80)
WRITE(24,80)
80 FORMAT('****************************************************')
181 FORMAT('************ OPERATORS ************')
WRITE(24,80)
WRITE(24,80)
WRITE(24,81)
WRITE(24,80)
WRITE(24,181)
WRITE(24,80)
WRITE(24,81)
81 FORMAT(' *')
WRITE(24,91)
91 FORMAT(' AVERAGE  X,Y  LINES ')
WRITE(24,81)
82 FORMAT(' ******************* ')
!
do ID=1, NTOTALMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
do IB=1,NUMBIN
AV(IB)=real(AVACLPMOM(IJ,IB,ID))/ACOUNT
13 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
do IB=1,NUMBIN
AV(IB)=aimag(AVACLPMOM(IJ,IB,ID))/ACOUNT
14 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(24,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101 FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10 CONTINUE
end do
!     
WRITE(24,80)
WRITE(24,81)
WRITE(24,92)
92 FORMAT(' AVERAGE   Y   LINES ')
WRITE(24,81)

do ID=1, NTOTALMOM
WRITE(24,81)
WRITE(24,82)
WRITE(24,94) ID
94 FORMAT(' OPERATOR = ',I4)
WRITE(24,82)
DO NT=1,LMAXIR
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOM(IJ,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAXIR,CORP,VAC)
!
!         DO 26 NT=1,MAXDTLS
do NT=1, LMAXIR-1
WRITE(24,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
26 CONTINUE
end do
20 CONTINUE
end do
!
end do
!
CLOSE(24)
RETURN
end subroutine MLINE2
!*********************************************************************
!*************************SUBROUTINE MLINE3***************************
!*********************************************************************
SUBROUTINE MLINE3
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
PARAMETER(LMAXIR=3)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!******************************************************************************
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************   
COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLP(NUMBIN,NTOTAL)
!
COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP) &
  & ,AVACLJ0PP(NUMBIN,NOPJ0PP)
COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM) &
  & ,AVACLJ0PM(NUMBIN,NOPJ0PM)
COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP) &
  & ,AVACLJ0MP(NUMBIN,NOPJ0MP)
COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM) &
  & ,AVACLJ0MM(NUMBIN,NOPJ0MM)
COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P) &
  & ,AVACLJ1P(NUMBIN,NOPJ1P)
COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M) &
  & ,AVACLJ1M(NUMBIN,NOPJ1M)
COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP) &
  & ,AVACLJ2PP(NUMBIN,NOPJ2PP)
COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM) &
  & ,AVACLJ2PM(NUMBIN,NOPJ2PM)
COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP) &
  & ,AVACLJ2MP(NUMBIN,NOPJ2MP)
COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM) &
  & ,AVACLJ2MM(NUMBIN,NOPJ2MM)
!      
COMMON/ITEM/ITR,IBIN,ITOT
COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
!
complex(real64) ACORLP,AVACLP
!***********************************************************
complex(real64) ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
complex(real64) AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
!***********************************************************
complex(real64) ACORLJ1P,ACORLJ1M
complex(real64) AVACLJ1P,AVACLJ1M
!***********************************************************
complex(real64) ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
complex(real64) AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
!***********************************************************
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE='DIAGNALQ0IND.DAT')
!
NBIN=ITOT/IBIN
LS(1)=LX1
LS(2)=LX2
LS(3)=LX3
EPS=0.0000000001d0
!********************************************************************
!     OPERATORS
!********************************************************************
WRITE(23,80)
WRITE(23,80)
80 FORMAT('****************************************************')
WRITE(23,181)
181 FORMAT('****************** OPERATORS WITH Q = 0 ************')
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(' *')
WRITE(23,91)
91 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+, Pr=+ ')
WRITE(23,81)
82 FORMAT(' ******************* ')
!
do ID=1, NOPJ0PP
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
do IB=1,NUMBIN
AV(IB)=real(AVACLJ0PP(IB,ID))/ACOUNT
13 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
do IB=1,NUMBIN
AV(IB)=aimag(AVACLJ0PP(IB,ID))/ACOUNT
14 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101 FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10 CONTINUE
end do
WRITE(23,81)
WRITE(23,201)
201 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+, Pr=- ')
WRITE(23,81)
!
DO ID=1, NOPJ0PM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ0PM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ0PM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,202)
202 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=-, Pr=+ ')
WRITE(23,81)
!
DO ID=1, NOPJ0MP
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ0MP(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ0MP(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,203)
203 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=-, Pr=- ')
WRITE(23,81)
!
DO ID=1, NOPJ0MM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ0MM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ0MM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,204)
204 FORMAT(' AVERAGE  X,Y  LINES, J=1, Pr=+ ')
WRITE(23,81)
!
DO ID=1, NOPJ1P
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ1P(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ1P(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,205)
205 FORMAT(' AVERAGE  X,Y  LINES, J=1, Pr=- ')
WRITE(23,81)
!
DO ID=1, NOPJ1M
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ1M(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ1M(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,206)
206 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+, Pr=+ ')
WRITE(23,81)
!
DO ID=1, NOPJ2PP
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ2PP(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ2PP(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,208)
208 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+, Pr=- ')
WRITE(23,81)
!
DO ID=1, NOPJ2PM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ2PM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ2PM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,209)
209 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=-, Pr=+ ')
WRITE(23,81)
!
DO ID=1, NOPJ2MP
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ2MP(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ2MP(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,210)
210 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=-, Pr=- ')
WRITE(23,81)
!
DO ID=1, NOPJ2MM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLJ2MM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLJ2MM(IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
!**********************************************************************
!     DIAGONAL CORRELATORS
!**********************************************************************      
WRITE(23,80)

WRITE(23,81)
WRITE(23,220)
220 FORMAT(' DIAGONAL TORELONS J=0, Pp=+, Pr=+, q=0 ')
WRITE(23,81)

do ID=1, NOPJ0PP
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(' OPERATOR = ',I4)
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ0PP(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
end do
20 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,211)
211 FORMAT(' DIAGONAL TORELONS J=0, Pp=+, Pr=-, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ0PM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ0PM(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 21
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
21 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,212)
212 FORMAT(' DIAGONAL TORELONS J=0, Pp=-, Pr=+, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ0MP
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ0MP(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 22
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
22 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,213)
213 FORMAT(' DIAGONAL TORELONS J=0, Pp=-, Pr=-, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ0MM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ0MM(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 23
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
23 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,214)
214 FORMAT(' DIAGONAL TORELONS J=1, Pr=+, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ1P
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ1P(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 24
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
24 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,215)
215 FORMAT(' DIAGONAL TORELONS J=1, Pr=-, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ1M
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ1M(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 25
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
25 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,216)
216 FORMAT(' DIAGONAL TORELONS J=2, Pp=+, Pr=+, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ2PP
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ2PP(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 26
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
26 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,217)
217 FORMAT(' DIAGONAL TORELONS J=2, Pp=+, Pr=-, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ2PM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ2PM(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 27
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
27 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,218)
218 FORMAT(' DIAGONAL TORELONS J=2, Pp=-, Pr=+, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ2MP
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ2MP(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 28
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
28 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,219)
219 FORMAT(' DIAGONAL TORELONS J=2, Pp=-, Pr=-, q=0 ')
WRITE(23,81)
!
do ID=1, NOPJ2MM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ2MM(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 29
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
29 CONTINUE
end do
CLOSE(23)
!*********************************************************************   
RETURN
end subroutine MLINE3
!*********************************************************************
!*************************SUBROUTINE MLINE4***************************
!*********************************************************************
SUBROUTINE MLINE4
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
PARAMETER(LMAXIR=3)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!******************************************************************************
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************   
COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLP(NUMBIN,NTOTAL)
!
COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP) &
  & ,AVACLJ0PP(NUMBIN,NOPJ0PP)
COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM) &
  & ,AVACLJ0PM(NUMBIN,NOPJ0PM)
COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP) &
  & ,AVACLJ0MP(NUMBIN,NOPJ0MP)
COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM) &
  & ,AVACLJ0MM(NUMBIN,NOPJ0MM)
COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P) &
  & ,AVACLJ1P(NUMBIN,NOPJ1P)
COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M) &
  & ,AVACLJ1M(NUMBIN,NOPJ1M)
COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP) &
  & ,AVACLJ2PP(NUMBIN,NOPJ2PP)
COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM) &
  & ,AVACLJ2PM(NUMBIN,NOPJ2PM)
COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP) &
  & ,AVACLJ2MP(NUMBIN,NOPJ2MP)
COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM) &
  & ,AVACLJ2MM(NUMBIN,NOPJ2MM)
!
COMMON/PBLOCKMOMJ0P/ &
  & ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM) &
  & ,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
COMMON/PBLOCKMOMJ0M/ &
  & ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM) &
  & ,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)
COMMON/PBLOCKMOMJ1/ &
  & ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM) &
  & ,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)
COMMON/PBLOCKMOMJ2P/ &
  & ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM) &
  & ,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
COMMON/PBLOCKMOMJ2M/ &
  & ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM) &
  & ,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
!
COMMON/PBLOCKMOMJ0/ &
  & ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM) &
  & ,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
COMMON/PBLOCKMOMJ2/ &
  & ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM) &
  & ,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)
!
COMMON/ITEM/ITR,IBIN,ITOT
COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
!
complex(real64) ACORLP,AVACLP
!***********************************************************
complex(real64) ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
complex(real64) AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
!***********************************************************
complex(real64) ACORLJ1P,ACORLJ1M
complex(real64) AVACLJ1P,AVACLJ1M
!***********************************************************
complex(real64) ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
complex(real64) AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
!***********************************************************
complex(real64) ACORLPMOMJ0P, ACORLPMOMJ0M
complex(real64) AVACLMOMJ0P, AVACLMOMJ0M
!***********************************************************
complex(real64) ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
complex(real64) AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
!***********************************************************
complex(real64) ACORLPMOMJ2P, ACORLPMOMJ2M
complex(real64) AVACLMOMJ2P, AVACLMOMJ2M
!***********************************************************
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE='DIAGNALQ1IND.DAT')
!
NBIN=ITOT/IBIN
LS(1)=LX1
LS(2)=LX2
LS(3)=LX3
EPS=0.0000000001d0
!********************************************************************
!     OPERATORS
!********************************************************************
WRITE(23,80)
WRITE(23,80)
80 FORMAT('****************************************************')
WRITE(23,181)
181 FORMAT('****************** OPERATORS WITH Q = 1 ************')
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(' *')
WRITE(23,91)
91 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+')
WRITE(23,81)
82 FORMAT(' ******************* ')
!
do ID=1, NOPJ0PMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
do IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ0P(1,IB,ID))/ACOUNT
13 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
do IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ0P(1,IB,ID))/ACOUNT
14 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101 FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10 CONTINUE
end do
WRITE(23,81)
WRITE(23,202)
202 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=- ')
WRITE(23,81)
!
DO ID=1, NOPJ0MMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ0M(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ0M(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,204)
204 FORMAT(' AVERAGE  X,Y  LINES, J=1 ')
WRITE(23,81)
!
DO ID=1, NOPJ1MOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ1(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ1(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,206)
206 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+ ')
WRITE(23,81)
!
DO ID=1, NOPJ2PMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ2P(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ2P(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,207)
207 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=- ')
WRITE(23,81)
!
DO ID=1, NOPJ2MMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ2M(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ2M(1,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
!**********************************************************************
!     DIAGONAL CORRELATORS smirloglou 00302102473443, 00306972428749
!**********************************************************************      
WRITE(23,80)

WRITE(23,81)
WRITE(23,220)
220 FORMAT(' DIAGONAL TORELONS J=0, Pp=+, q=1 ')
WRITE(23,81)

do ID=1, NOPJ0PMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(' OPERATOR = ',I4)
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ0P(1,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
end do
20 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,211)
211 FORMAT(' DIAGONAL TORELONS J=0, Pp=-, q=1 ')
WRITE(23,81)
!
do ID=1, NOPJ0MMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ0M(1,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 21
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
21 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,212)
212 FORMAT(' DIAGONAL TORELONS J=1, q=1 ')
WRITE(23,81)
!
do ID=1, NOPJ1MOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ1(1,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 22
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
22 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,213)
213 FORMAT(' DIAGONAL TORELONS J=2, Pp=+, q=1 ')
WRITE(23,81)
!
do ID=1, NOPJ2PMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ2P(1,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 23
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
23 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,214)
214 FORMAT(' DIAGONAL TORELONS J=2, Pr=-, q=1 ')
WRITE(23,81)
!
do ID=1, NOPJ2MMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ2M(1,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 24
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
24 CONTINUE
end do
CLOSE(23)
!*********************************************************************   
RETURN
end subroutine MLINE4
!*********************************************************************
!*************************SUBROUTINE MLINE5***************************
!*********************************************************************
SUBROUTINE MLINE5
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1,NUMBIN=2)
PARAMETER(LMAXIR=3)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!******************************************************************************
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************   
COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLP(NUMBIN,NTOTAL)
!
COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP) &
  & ,AVACLJ0PP(NUMBIN,NOPJ0PP)
COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM) &
  & ,AVACLJ0PM(NUMBIN,NOPJ0PM)
COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP) &
  & ,AVACLJ0MP(NUMBIN,NOPJ0MP)
COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM) &
  & ,AVACLJ0MM(NUMBIN,NOPJ0MM)
COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P) &
  & ,AVACLJ1P(NUMBIN,NOPJ1P)
COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M) &
  & ,AVACLJ1M(NUMBIN,NOPJ1M)
COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP) &
  & ,AVACLJ2PP(NUMBIN,NOPJ2PP)
COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM) &
  & ,AVACLJ2PM(NUMBIN,NOPJ2PM)
COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP) &
  & ,AVACLJ2MP(NUMBIN,NOPJ2MP)
COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM) &
  & ,AVACLJ2MM(NUMBIN,NOPJ2MM)
!
COMMON/PBLOCKMOMJ0P/ &
  & ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM) &
  & ,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
COMMON/PBLOCKMOMJ0M/ &
  & ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM) &
  & ,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)
COMMON/PBLOCKMOMJ1/ &
  & ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM) &
  & ,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)
COMMON/PBLOCKMOMJ2P/ &
  & ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM) &
  & ,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
COMMON/PBLOCKMOMJ2M/ &
  & ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM) &
  & ,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
!
COMMON/PBLOCKMOMJ0/ &
  & ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM) &
  & ,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
COMMON/PBLOCKMOMJ2/ &
  & ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM) &
  & ,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)
!
COMMON/ITEM/ITR,IBIN,ITOT
COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
!
complex(real64) ACORLP,AVACLP
!***********************************************************
complex(real64) ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
complex(real64) AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
!***********************************************************
complex(real64) ACORLJ1P,ACORLJ1M
complex(real64) AVACLJ1P,AVACLJ1M
!***********************************************************
complex(real64) ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
complex(real64) AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
!***********************************************************
complex(real64) ACORLPMOMJ0P, ACORLPMOMJ0M
complex(real64) AVACLMOMJ0P, AVACLMOMJ0M
!***********************************************************
complex(real64) ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
complex(real64) AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
!***********************************************************
complex(real64) ACORLPMOMJ2P, ACORLPMOMJ2M
complex(real64) AVACLMOMJ2P, AVACLMOMJ2M
!***********************************************************
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE='DIAGNALQ2IND.DAT')
!
NBIN=ITOT/IBIN
LS(1)=LX1
LS(2)=LX2
LS(3)=LX3
EPS=0.0000000001d0
!********************************************************************
!     OPERATORS
!********************************************************************
WRITE(23,80)
WRITE(23,80)
80 FORMAT('****************************************************')
WRITE(23,181)
181 FORMAT('****************** OPERATORS WITH Q = 2 ************')
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(' *')
WRITE(23,91)
91 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=+')
WRITE(23,81)
82 FORMAT(' ******************* ')
!
do ID=1, NOPJ0PMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
do IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ0P(2,IB,ID))/ACOUNT
13 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
do IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ0P(2,IB,ID))/ACOUNT
14 CONTINUE
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
101 FORMAT('OP =',I4,'  AV,ER R,I LINES =',4F9.4)
10 CONTINUE
end do
WRITE(23,81)
WRITE(23,202)
202 FORMAT(' AVERAGE  X,Y  LINES, J=0, Pp=- ')
WRITE(23,81)
!
DO ID=1, NOPJ0MMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ0M(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ0M(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,204)
204 FORMAT(' AVERAGE  X,Y  LINES, J=1 ')
WRITE(23,81)
!
DO ID=1, NOPJ1MOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ1(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ1(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,206)
206 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=+ ')
WRITE(23,81)
!
DO ID=1, NOPJ2PMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ2P(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ2P(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
WRITE(23,81)
WRITE(23,207)
207 FORMAT(' AVERAGE  X,Y  LINES, J=2, Pp=- ')
WRITE(23,81)
!
DO ID=1, NOPJ2MMOM
ACOUNT=NCOL*IBIN*LSIZE/LS(1)
DO IB=1,NUMBIN
AV(IB)=real(AVACLMOMJ2M(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQRP,SPLAQRP)
DO IB=1,NUMBIN
AV(IB)=aimag(AVACLMOMJ2M(2,IB,ID))/ACOUNT
end do
CALL JACKM(NUMBIN,AV,APLAQIP,SPLAQIP)
WRITE(23,101) ID,APLAQRP,APLAQIP,SPLAQRP,SPLAQIP
end do
!**********************************************************************
!     DIAGONAL CORRELATORS smirloglou 00302102473443, 00306972428749
!**********************************************************************      
WRITE(23,80)

WRITE(23,81)
WRITE(23,220)
220 FORMAT(' DIAGONAL TORELONS J=0, Pp=+, q=2 ')
WRITE(23,81)

do ID=1, NOPJ0PMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(' OPERATOR = ',I4)
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ0P(2,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
end do
20 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,211)
211 FORMAT(' DIAGONAL TORELONS J=0, Pp=-, q=2 ')
WRITE(23,81)
!
do ID=1, NOPJ0MMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ0M(2,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 21
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
21 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,212)
212 FORMAT(' DIAGONAL TORELONS J=1, q=2 ')
WRITE(23,81)
!
do ID=1, NOPJ1MOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ1(2,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 22
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
22 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,213)
213 FORMAT(' DIAGONAL TORELONS J=2, Pp=+, q=2 ')
WRITE(23,81)
!
do ID=1, NOPJ2PMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ2P(2,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 23
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
23 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,214)
214 FORMAT(' DIAGONAL TORELONS J=2, Pr=-, q=2 ')
WRITE(23,81)
!
do ID=1, NOPJ2MMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ2M(2,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1)).LE.EPS)go to 24
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
24 CONTINUE
end do
CLOSE(23)
!*********************************************************************   
RETURN
end subroutine MLINE5
!*********************************************************************
!************************ SUBROUTINE BLOCK *************************** 
!*********************************************************************
SUBROUTINE BLOCK
use iso_fortran_env, only : real32, real64, int32
use correlator_construction_mod, only : pot
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(IBLOK=5,NBLOK=2,NSMEAR=12)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
COMMON/ARRAYS/U11(NCOL2,LSIZEB,LX4,4)
COMMON/ARRAYB/UB11(NCOL2,LSIZEB,3,IBLOK)
COMMON/ASMEAR1/UC11(NCOL2,LSIZEB,3)
COMMON/ASMEAR2/IUP(LSIZEB,3),IDN(LSIZEB,3)
COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2)
complex(real32) U11,UB11,A11,B11,C11,UC11
!
CALL cpu_time(t1)
do I4=1,LX4
!
do IBL=1,IBLOK
!
IF(IBL.EQ.1)THEN
DO MU=1,3
DO NN=1,LSIZEB
DO IJ=1,NCOL2
UB11(IJ,NN,MU,1)=U11(IJ,NN,I4,MU)
end do
end do
end do
go to 11
end if
!
IBLM=IBL-1
DO KK=1,3
DO NN=1,LSIZEB
IUP(NN,KK)=IUPB(NN,KK,IBLM)
IDN(NN,KK)=IDNB(NN,KK,IBLM)
DO IJ=1,NCOL2
UC11(IJ,NN,KK)=UB11(IJ,NN,KK,IBLM)
end do
end do
end do
CALL SMEAR1
!
DO MU=1,3
DO NN=1,LSIZEB
M1=NN
M2=IUP(M1,MU)
!
do IJ=1,NCOL2
A11(IJ)=UC11(IJ,M1,MU)
B11(IJ)=UC11(IJ,M2,MU)
22 CONTINUE
end do
CALL VMX(1,A11,B11,C11,1)
do IJ=1,NCOL2
UB11(IJ,NN,MU,IBL)=C11(IJ)
24 CONTINUE
end do
!
end do
end do
11 CONTINUE
!
CALL THERML1(I4,IBL)
!
!************************************************
3 CONTINUE
end do
1 CONTINUE
end do
CALL cpu_time(t2)
WRITE(6,901) real(t2-t1, kind=real32)
901 FORMAT('[Info][Time]','     Thermal Line time:', F8.4)
!
CALL cpu_time(t1)
CALL POT
CALL cpu_time(t2)
WRITE(6,902) real(t2-t1, kind=real32)
902 FORMAT('[Info][Time]','     Correlation Creation time:', F8.4)
!
RETURN
end subroutine BLOCK
!*********************************************************************
SUBROUTINE SMEAR1
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(PARBS=0.30,PARBDS=0.12)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
PARAMETER(IDIAG=1)
!
COMMON/ASMEAR1/UC11(NCOL2,LSIZEB,3)
COMMON/ASMEAR2/IUP(LSIZEB,3),IDN(LSIZEB,3)
COMMON/DIAGIN/UUC11(NCOL2,LSIZEB,3), &
  & IIUP(LSIZEB,3),IIDN(LSIZEB,3)
COMMON/DIAGOUT/UDD(NCOL2,LSIZEB,4),IDD(LSIZEB,4)
DIMENSION UCC11(NCOL2,LSIZEB,3)
DIMENSION UINT11(NCOL2),UREN11(NCOL2)
DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2),DUM11(NCOL2)
complex(real32) A11,B11,C11,UINT11,UREN11,UC11,DUM11,UCC11,UDD
complex(real32) UUC11
!
DO KK=1,3
DO NN=1,LSIZEB
IIUP(NN,KK)=IUP(NN,KK)
IIDN(NN,KK)=IDN(NN,KK)
end do
end do
!
DO MU=1,3
DO NN=1,LSIZEB
DO IJ=1,NCOL2
UUC11(IJ,NN,MU)=UC11(IJ,NN,MU)
end do
end do
end do
!
do MU=1,3
IF(IDIAG.EQ.1) CALL DIAG(MU)
!
do NN=1,LSIZEB
M1=NN
!
do IJ=1,NCOL2
UINT11(IJ)=UC11(IJ,M1,MU)
10 CONTINUE
end do
!
do NU=1,3
IF(MU.EQ.NU)go to 40
!
do IJ=1,NCOL2
A11(IJ)=UC11(IJ,M1,NU)
26 CONTINUE
end do
M2=IUP(M1,NU)
do IJ=1,NCOL2
B11(IJ)=UC11(IJ,M2,MU)
27 CONTINUE
end do
CALL VMX(1,A11,B11,C11,1)
M3=IUP(M1,MU)
do IJ=1,NCOL2
B11(IJ)=UC11(IJ,M3,NU)
28 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,C11,B11,A11,1)
do IJ=1,NCOL2
UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARBS
21 CONTINUE
end do
!
M2=IDN(M1,NU)
do IJ=1,NCOL2
A11(IJ)=UC11(IJ,M2,NU)
31 CONTINUE
end do
do IJ=1,NCOL2
B11(IJ)=UC11(IJ,M2,MU)
32 CONTINUE
end do
CALL HERM(1,A11,DUM11,1)
CALL VMX(1,A11,B11,C11,1)
M3=IUP(M2,MU)
do IJ=1,NCOL2
B11(IJ)=UC11(IJ,M3,NU)
33 CONTINUE
end do
CALL VMX(1,C11,B11,A11,1)
do IJ=1,NCOL2
UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARBS
35 CONTINUE
end do
!
40 CONTINUE
end do
!
do MNU=1,4
IF(IDIAG.NE.1)go to 50
!
do IJ=1,NCOL2
A11(IJ)=UDD(IJ,M1,MNU)
42 CONTINUE
end do
M2=IDD(M1,MNU)
do IJ=1,NCOL2
B11(IJ)=UC11(IJ,M2,MU)
43 CONTINUE
end do
CALL VMX(1,A11,B11,C11,1)
M3=IUP(M1,MU)
do IJ=1,NCOL2
B11(IJ)=UDD(IJ,M3,MNU)
44 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,C11,B11,A11,1)
do IJ=1,NCOL2
UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARBDS
45 CONTINUE
end do
!
50 CONTINUE
end do
!
CALL RENORMBS(UINT11,UREN11)
do IJ=1,NCOL2
UCC11(IJ,M1,MU)=UREN11(IJ)
60 CONTINUE
end do
!
2 CONTINUE
end do
1 CONTINUE
end do
!
DO MU=1,3
DO NN=1,LSIZEB
DO IJ=1,NCOL2
UC11(IJ,NN,MU)=UCC11(IJ,NN,MU)
end do
end do
end do
CALL DODET
!
RETURN
end subroutine SMEAR1
!**********************************************************
!  THIS ROUTINE MAKES THE SMEARED LINKS HAVE DET=1
!**********************************************************
SUBROUTINE DODET
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
COMMON/ASMEAR1/U11(NCOL,NCOL,LSIZEB,3)
complex(real32) ADUM(NCOL,NCOL),U11,CSUM
!
complex(real32) AA(NCOL,NCOL)
complex(real32) DET,CNORM,CDET,DNCOL
!
DO MU=1,3
DO NN=1,LSIZEB
!
DO J=1,NCOL
DO I=1,NCOL
ADUM(I,J)=U11(I,J,NN,MU)
end do
end do
!
DO J=1, NCOL
DO I=1,NCOL
AA(I,J)=U11(I,J,NN,MU)
end do
end do
JMAT=NCOL
CALL DETNANT(JMAT,DET,AA)
!
DNCOL=CMPLX(1.0/NCOL)
CDET=DET**DNCOL
CNORM=1.0/CDET
DO J=1,NCOL
DO I=1, NCOL
AA(I,J)=ADUM(I,J)*CNORM
U11(I,J,NN,MU)=AA(I,J)
!
end do
end do
!
end do
end do
RETURN
end subroutine DODET
!*********************************************************************
! form 4 'diagonal' links in plane orth to MU, from each site 
! with end-pt in IDD using links and pointers in /DIAGIN/
! -- output matrices (not suN) and pointers in /DIAGOUT/
!*********************************************************************
SUBROUTINE DIAG(MU)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZEB=LX1*LX2*LX3)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!     
COMMON/DIAGOUT/UDD(NCOL2,LSIZEB,4),IDD(LSIZEB,4)
COMMON/DIAGIN/UC11(NCOL2,LSIZEB,3), &
  & IUP(LSIZEB,3),IDN(LSIZEB,3)
complex(real32) C11(NCOL2),F11(NCOL2),DUM11(NCOL2)
complex(real32) A11(NCOL2),D11(NCOL2),E11(NCOL2)
complex(real32) AA11(NCOL2),DD11(NCOL2),EE11(NCOL2)
complex(real32) UDD,UC11
!     
NU=MU+1
IF(MU.EQ.3)NU=1
KU=6-MU-NU
!     
do NN=1,LSIZEB
M1=NN
!     
do IJ=1,NCOL2
A11(IJ)=UC11(IJ,M1,NU)
2 CONTINUE
end do
M2=IUP(M1,NU)
do IJ=1,NCOL2
AA11(IJ)=UC11(IJ,M2,KU)
3 CONTINUE
end do
CALL VMX(1,A11,AA11,C11,1)
do IJ=1,NCOL2
D11(IJ)=UC11(IJ,M1,KU)
4 CONTINUE
end do
M3=IUP(M1,KU)
do IJ=1,NCOL2
DD11(IJ)=UC11(IJ,M3,NU)
5 CONTINUE
end do
CALL VMX(1,D11,DD11,F11,1)
do IJ=1,NCOL2
UDD(IJ,M1,1)=C11(IJ)+F11(IJ)
6 CONTINUE
end do
M9=IUP(M3,NU)
IDD(M1,1)=M9
!     
do IJ=1,NCOL2
AA11(IJ)=C11(IJ)+F11(IJ)
61 CONTINUE
end do
CALL HERM(1,AA11,DUM11,1)
do IJ=1,NCOL2
UDD(IJ,M9,3)=AA11(IJ)
62 CONTINUE
end do
IDD(M9,3)=M1
!     
M4=IDN(M2,KU)
do IJ=1,NCOL2
AA11(IJ)=UC11(IJ,M4,KU)
7 CONTINUE
end do
CALL HERM(1,AA11,DUM11,1)
CALL VMX(1,A11,AA11,C11,1)
M5=IDN(M1,KU)
do IJ=1,NCOL2
E11(IJ)=UC11(IJ,M5,KU)
8 CONTINUE
end do
CALL HERM(1,E11,DUM11,1)
do IJ=1,NCOL2
EE11(IJ)=UC11(IJ,M5,NU)
9 CONTINUE
end do
CALL VMX(1,E11,EE11,F11,1)
do IJ=1,NCOL2
UDD(IJ,M1,2)=C11(IJ)+F11(IJ)
10 CONTINUE
end do
IDD(M1,2)=M4
!     
do IJ=1,NCOL2
AA11(IJ)=C11(IJ)+F11(IJ)
63 CONTINUE
end do
CALL HERM(1,AA11,DUM11,1)
do IJ=1,NCOL2
UDD(IJ,M4,4)=AA11(IJ)
64 CONTINUE
end do
IDD(M4,4)=M1
!     
1 CONTINUE
end do
!     
RETURN
end subroutine DIAG
!**********************************************************
!  THIS ROUTINE MAKES BLOCKED MATRICES SU(N) IN CRUDEST   *
!               POSSIBLE WAY - 1 COOL                     *
!**********************************************************
SUBROUTINE RENORMBS(UU1,UREN11)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!     
DIMENSION UU1(NCOL,NCOL),UREN11(NCOL,NCOL)
complex(real32) UB1,UU1,CSUM,ADUM(NCOL,NCOL),UREN11
complex(real32) A11(NCOL2),B11(NCOL2),C11(NCOL2)
!
DO N2=1,NCOL
DO N1=1,NCOL
ADUM(N1,N2)=UU1(N1,N2)
end do
end do
!     
do N2=1,NCOL
!     
do N3=1,N2-1
!     
CSUM=(0.0,0.0)
do N1=1,NCOL
CSUM=CSUM+ADUM(N1,N2)*CONJG(ADUM(N1,N3))
5 CONTINUE
end do
do N1=1,NCOL
ADUM(N1,N2)=ADUM(N1,N2)-CSUM*ADUM(N1,N3)
6 CONTINUE
end do
!     
10 CONTINUE
end do
!     
SUM=0.0
do N1=1,NCOL
SUM=SUM+ADUM(N1,N2)*CONJG(ADUM(N1,N2))
7 CONTINUE
end do
ANORM=1.0/SQRT(SUM)
do N1=1,NCOL
ADUM(N1,N2)=ADUM(N1,N2)*ANORM
8 CONTINUE
end do
!     
20 CONTINUE
end do
!
DO IK=1,NCOL
DO IJ=1,NCOL
IADD=NCOL*(IK-1)
A11(IJ+IADD)=ADUM(IJ,IK)
end do
end do
!
do J=1,NCOL2
C11(J)=A11(J)
22 CONTINUE
end do
DO J=1,NCOL
DO I=1,NCOL
IJ=I+(J-1)*NCOL
A11(IJ)=CONJG(UU1(J,I))
end do
end do
!
CALL VMX(1,A11,C11,B11,1)
CALL RUNGRB(B11,C11)
!
DO J=1,NCOL
DO I=1,NCOL
IJ=I+(J-1)*NCOL
UREN11(I,J)=C11(IJ)
end do
end do
!
RETURN
end subroutine RENORMBS
! *******************************************************************
! *************         SUBROUTINE RUNGRB           *****************
! *******************************************************************
SUBROUTINE RUNGRB(B11,C11)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!     
complex(real32) B11(NCOL2),C11(NCOL2)
COMMON/SUBB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
!     
do LDU=1,NCOL-1
do LDL=LDU+1,NCOL
II1=LDU
JJ1=LDU
II2=LDL
JJ2=LDL
II3=LDL
JJ3=LDU
II4=LDU
JJ4=LDL
CALL SUBGRB(LDU,LDL,B11,C11)
15 CONTINUE
end do
10 CONTINUE
end do
!     
RETURN
end subroutine RUNGRB
!*********************************************************************
!*******************      SUBROUTINE SUBGRB      *********************
!*********************************************************************
SUBROUTINE SUBGRB(LDU,LDL,B11,C11)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!     
COMMON/SUBB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
!     
complex(real32) B11(NCOL,NCOL),C11(NCOL,NCOL)
complex(real32) F11,F12,A11(NCOL,NCOL),S11(NCOL,NCOL)
complex(real32) T11(NCOL,NCOL)
!     
F11=(B11(II1,JJ1)+CONJG(B11(II2,JJ2)))*0.5
F12=(B11(II3,JJ3)-CONJG(B11(II4,JJ4)))*0.5
UMAG=SQRT(F11*CONJG(F11)+F12*CONJG(F12))
UMAG=1./UMAG
F11=F11*UMAG
F12=F12*UMAG
!     
A11(II1,JJ1)=CONJG(F11)
A11(II4,JJ4)=CONJG(F12)
A11(II3,JJ3)=-F12
A11(II2,JJ2)=F11
!     
do IJ1=1,NCOL
S11(IJ1,LDU)=C11(IJ1,LDU)*(A11(LDU,LDU)-1.0) &
  & +C11(IJ1,LDL)*A11(LDL,LDU)
S11(IJ1,LDL)=C11(IJ1,LDU)*A11(LDU,LDL) &
  & +C11(IJ1,LDL)*(A11(LDL,LDL)-1.0)
T11(IJ1,LDU)=B11(IJ1,LDU)*(A11(LDU,LDU)-1.0) &
  & +B11(IJ1,LDL)*A11(LDL,LDU)
T11(IJ1,LDL)=B11(IJ1,LDU)*A11(LDU,LDL) &
  & +B11(IJ1,LDL)*(A11(LDL,LDL)-1.0)
3 CONTINUE
end do
do IJ1=1,NCOL
C11(IJ1,LDU)=C11(IJ1,LDU)+S11(IJ1,LDU)
C11(IJ1,LDL)=C11(IJ1,LDL)+S11(IJ1,LDL)
B11(IJ1,LDU)=B11(IJ1,LDU)+T11(IJ1,LDU)
B11(IJ1,LDL)=B11(IJ1,LDL)+T11(IJ1,LDL)
4 CONTINUE
end do
!     
RETURN
end subroutine SUBGRB
!*********************************************************************
!************* HERE I CREATE THE PULSE - LIKE OPERATORS **************
!*********************************************************************
SUBROUTINE THERML1(N4,IBLL)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!********************************************************
COMMON/LINES/ALINE1(LX4,IBLOK),ALINE2(LX4,IBLOK) &
  & ,ALINE3(LX4,IBLOK),ALINE4(LX4,IBLOK),ALINE5(LX4,IBLOK) &
  & ,ALINE6(LX4,IBLOK),ALINE7(LX4,IBLOK),ALINE8(LX4,IBLOK) &
  & ,ALINE9(LX4,IBLOK),ALINE10(LX4,IBLOK),ALINE11(LX4,IBLOK) &
  & ,ALINE12(LX4,IBLOK),ALINE13(LX4,IBLOK),ALINE14(LX4,IBLOK) &
  & ,ALINE15(LX4,IBLOK),ALINE16(LX4,IBLOK),ALINE17(LX4,IBLOK) &
  & ,ALINE18(LX4,IBLOK),ALINE19(LX4,IBLOK),ALINE20(LX4,IBLOK) &
  & ,ALINE21(LX4,IBLOK),ALINE22(LX4,IBLOK),ALINE23(LX4,IBLOK) &
  & ,ALINE24(LX4,IBLOK),ALINE25(LX4,IBLOK),ALINE26(LX4,IBLOK) &
  & ,ALINE27(LX4,IBLOK),ALINE28(LX4,IBLOK),ALINE29(LX4,IBLOK) &
  & ,ALINE30(LX4,IBLOK),ALINE31(LX4,IBLOK),ALINE32(LX4,IBLOK) &
  & ,ALINE33(LX4,IBLOK),ALINE34(LX4,IBLOK),ALINE35(LX4,IBLOK) &
  & ,ALINE36(LX4,IBLOK),ALINE37(LX4,IBLOK),ALINE38(LX4,IBLOK) &
  & ,ALINE39(LX4,IBLOK),ALINE40(LX4,IBLOK),ALINE41(LX4,IBLOK) &
  & ,ALINE42(LX4,IBLOK),ALINE43(LX4,IBLOK),ALINE44(LX4,IBLOK) &
  & ,ALINE45(LX4,IBLOK),ALINE46(LX4,IBLOK),ALINE47(LX4,IBLOK) &
  & ,ALINE48(LX4,IBLOK),ALINE49(LX4,IBLOK),ALINE50(LX4,IBLOK) &
  & ,ALINE51(LX4,IBLOK),ALINE52(LX4,IBLOK),ALINE53(LX4,IBLOK) &
  & ,ALINE54(LX4,IBLOK),ALINE55(LX4,IBLOK),ALINE56(LX4,IBLOK) &
  & ,ALINE57(LX4,IBLOK),ALINE58(LX4,IBLOK),ALINE59(LX4,IBLOK) &
  & ,ALINE60(LX4,IBLOK),ALINE61(LX4,IBLOK),ALINE62(LX4,IBLOK) &
  & ,ALINE63(LX4,IBLOK),ALINE64(LX4,IBLOK),ALINE65(LX4,IBLOK) &
  & ,ALINE66(LX4,IBLOK),ALINE67(LX4,IBLOK),ALINE68(LX4,IBLOK) &
  & ,ALINE69(LX4,IBLOK),ALINE70(LX4,IBLOK),ALINE71(LX4,IBLOK) &
  & ,ALINE72(LX4,IBLOK),ALINE73(LX4,IBLOK),ALINE74(LX4,IBLOK) &
  & ,ALINE75(LX4,IBLOK),ALINE76(LX4,IBLOK),ALINE77(LX4,IBLOK) &
  & ,ALINE78(LX4,IBLOK),ALINE79(LX4,IBLOK),ALINE80(LX4,IBLOK) &
  & ,ALINE81(LX4,IBLOK),ALINE82(LX4,IBLOK),ALINE83(LX4,IBLOK) &
  & ,ALINE84(LX4,IBLOK),ALINE85(LX4,IBLOK),ALINE86(LX4,IBLOK) &
  & ,ALINE87(LX4,IBLOK),ALINE88(LX4,IBLOK),ALINE89(LX4,IBLOK) &
  & ,ALINE90(LX4,IBLOK),ALINE91(LX4,IBLOK),ALINE92(LX4,IBLOK) &
  & ,ALINE93(LX4,IBLOK),ALINE94(LX4,IBLOK),ALINE95(LX4,IBLOK) &
  & ,ALINE96(LX4,IBLOK),ALINE97(LX4,IBLOK),ALINE98(LX4,IBLOK) &
  & ,ALINE99(LX4,IBLOK),ALINE100(LX4,IBLOK),ALINE101(LX4,IBLOK) &
  & ,ALINE102(LX4,IBLOK),ALINE103(LX4,IBLOK),ALINE104(LX4,IBLOK) &
  & ,ALINE105(LX4,IBLOK),ALINE106(LX4,IBLOK),ALINE107(LX4,IBLOK) &
  & ,ALINE108(LX4,IBLOK),ALINE109(LX4,IBLOK),ALINE110(LX4,IBLOK) &
  & ,ALINE111(LX4,IBLOK),ALINE112(LX4,IBLOK),ALINE113(LX4,IBLOK) &
  & ,ALINE114(LX4,IBLOK),ALINE115(LX4,IBLOK),ALINE116(LX4,IBLOK) &
  & ,ALINE117(LX4,IBLOK),ALINE118(LX4,IBLOK),ALINE119(LX4,IBLOK) &
  & ,ALINE120(LX4,IBLOK),ALINE121(LX4,IBLOK),ALINE122(LX4,IBLOK) &
  & ,ALINE123(LX4,IBLOK),ALINE124(LX4,IBLOK),ALINE125(LX4,IBLOK) &
  & ,ALINE126(LX4,IBLOK),ALINE127(LX4,IBLOK),ALINE128(LX4,IBLOK) &
  & ,ALINE129(LX4,IBLOK),ALINE130(LX4,IBLOK),ALINE131(LX4,IBLOK) &
  & ,ALINE132(LX4,IBLOK),ALINE133(LX4,IBLOK),ALINE134(LX4,IBLOK) &
  & ,ALINE135(LX4,IBLOK),ALINE136(LX4,IBLOK),ALINE137(LX4,IBLOK) &
  & ,ALINE138(LX4,IBLOK),ALINE139(LX4,IBLOK),ALINE140(LX4,IBLOK) &
  & ,ALINE141(LX4,IBLOK),ALINE142(LX4,IBLOK),ALINE143(LX4,IBLOK) &
  & ,ALINE144(LX4,IBLOK),ALINE145(LX4,IBLOK),ALINE146(LX4,IBLOK) &
  & ,ALINE147(LX4,IBLOK),ALINE148(LX4,IBLOK),ALINE149(LX4,IBLOK) &
  & ,ALINE150(LX4,IBLOK),ALINE151(LX4,IBLOK),ALINE152(LX4,IBLOK) &
  & ,ALINE153(LX4,IBLOK),ALINE154(LX4,IBLOK),ALINE155(LX4,IBLOK) &
  & ,ALINE156(LX4,IBLOK),ALINE157(LX4,IBLOK),ALINE158(LX4,IBLOK) &
  & ,ALINE159(LX4,IBLOK),ALINE160(LX4,IBLOK),ALINE161(LX4,IBLOK) &
  & ,ALINE162(LX4,IBLOK),ALINE163(LX4,IBLOK),ALINE164(LX4,IBLOK) &
  & ,ALINE165(LX4,IBLOK),ALINE166(LX4,IBLOK),ALINE167(LX4,IBLOK) &
  & ,ALINE168(LX4,IBLOK),ALINE169(LX4,IBLOK),ALINE170(LX4,IBLOK) &
  & ,ALINE171(LX4,IBLOK),ALINE172(LX4,IBLOK),ALINE173(LX4,IBLOK) &
  & ,ALINE174(LX4,IBLOK),ALINE175(LX4,IBLOK),ALINE176(LX4,IBLOK) &
  & ,ALINE177(LX4,IBLOK),ALINE178(LX4,IBLOK),ALINE179(LX4,IBLOK) &
  & ,ALINE180(LX4,IBLOK),ALINE181(LX4,IBLOK),ALINE182(LX4,IBLOK) &
  & ,ALINE183(LX4,IBLOK),ALINE184(LX4,IBLOK),ALINE185(LX4,IBLOK) &
  & ,ALINE186(LX4,IBLOK),ALINE187(LX4,IBLOK),ALINE188(LX4,IBLOK) &
  & ,ALINE189(LX4,IBLOK),ALINE190(LX4,IBLOK),ALINE191(LX4,IBLOK) &
  & ,ALINE192(LX4,IBLOK),ALINE193(LX4,IBLOK),ALINE194(LX4,IBLOK) &
  & ,ALINE195(LX4,IBLOK),ALINE196(LX4,IBLOK),ALINE197(LX4,IBLOK) &
  & ,ALINE198(LX4,IBLOK),ALINE199(LX4,IBLOK),ALINE200(LX4,IBLOK) &
  & ,ALINE201(LX4,IBLOK),ALINE202(LX4,IBLOK),ALINE203(LX4,IBLOK) &
  & ,ALINE204(LX4,IBLOK),ALINE205(LX4,IBLOK),ALINE206(LX4,IBLOK) &
  & ,ALINE207(LX4,IBLOK),ALINE208(LX4,IBLOK),ALINE209(LX4,IBLOK) &
  & ,ALINE210(LX4,IBLOK),ALINE211(LX4,IBLOK),ALINE212(LX4,IBLOK) &
  & ,ALINE213(LX4,IBLOK),ALINE214(LX4,IBLOK),ALINE215(LX4,IBLOK) &
  & ,ALINE216(LX4,IBLOK),ALINE217(LX4,IBLOK),ALINE218(LX4,IBLOK) &
  & ,ALINE219(LX4,IBLOK),ALINE220(LX4,IBLOK),ALINE221(LX4,IBLOK) &
  & ,ALINE222(LX4,IBLOK),ALINE223(LX4,IBLOK),ALINE224(LX4,IBLOK) &
  & ,ALINE225(LX4,IBLOK),ALINE226(LX4,IBLOK),ALINE227(LX4,IBLOK) &
  & ,ALINE228(LX4,IBLOK),ALINE229(LX4,IBLOK),ALINE230(LX4,IBLOK) &
  & ,ALINE231(LX4,IBLOK),ALINE232(LX4,IBLOK),ALINE233(LX4,IBLOK) &
  & ,ALINE234(LX4,IBLOK),ALINE235(LX4,IBLOK)
!********************************************************
COMMON/LINESMOM/ALINEMOM2(LX4,IBLOK,2),ALINEMOM3(LX4,IBLOK,2) &
  & ,ALINEMOM4(LX4,IBLOK,2),ALINEMOM5(LX4,IBLOK,2) &
  & ,ALINEMOM6(LX4,IBLOK,2),ALINEMOM7(LX4,IBLOK,2) &
  & ,ALINEMOM8(LX4,IBLOK,2),ALINEMOM9(LX4,IBLOK,2) &
  & ,ALINEMOM10(LX4,IBLOK,2),ALINEMOM11(LX4,IBLOK,2) &
  & ,ALINEMOM12(LX4,IBLOK,2),ALINEMOM13(LX4,IBLOK,2) &
  & ,ALINEMOM14(LX4,IBLOK,2),ALINEMOM15(LX4,IBLOK,2) &
  & ,ALINEMOM16(LX4,IBLOK,2),ALINEMOM17(LX4,IBLOK,2) &
  & ,ALINEMOM18(LX4,IBLOK,2),ALINEMOM19(LX4,IBLOK,2) &
  & ,ALINEMOM20(LX4,IBLOK,2),ALINEMOM21(LX4,IBLOK,2) &
  & ,ALINEMOM22(LX4,IBLOK,2),ALINEMOM23(LX4,IBLOK,2) &
  & ,ALINEMOM24(LX4,IBLOK,2),ALINEMOM25(LX4,IBLOK,2) &
  & ,ALINEMOM26(LX4,IBLOK,2),ALINEMOM27(LX4,IBLOK,2) &
  & ,ALINEMOM28(LX4,IBLOK,2),ALINEMOM29(LX4,IBLOK,2) &
  & ,ALINEMOM30(LX4,IBLOK,2),ALINEMOM31(LX4,IBLOK,2) &
  & ,ALINEMOM32(LX4,IBLOK,2),ALINEMOM33(LX4,IBLOK,2) &
  & ,ALINEMOM34(LX4,IBLOK,2),ALINEMOM35(LX4,IBLOK,2) &
  & ,ALINEMOM36(LX4,IBLOK,2),ALINEMOM37(LX4,IBLOK,2) &
  & ,ALINEMOM38(LX4,IBLOK,2),ALINEMOM39(LX4,IBLOK,2) &
  & ,ALINEMOM40(LX4,IBLOK,2),ALINEMOM41(LX4,IBLOK,2) &
  & ,ALINEMOM42(LX4,IBLOK,2),ALINEMOM43(LX4,IBLOK,2) &
  & ,ALINEMOM44(LX4,IBLOK,2),ALINEMOM45(LX4,IBLOK,2) &
  & ,ALINEMOM46(LX4,IBLOK,2),ALINEMOM47(LX4,IBLOK,2) &
  & ,ALINEMOM48(LX4,IBLOK,2),ALINEMOM49(LX4,IBLOK,2) &
  & ,ALINEMOM50(LX4,IBLOK,2),ALINEMOM51(LX4,IBLOK,2) &
  & ,ALINEMOM52(LX4,IBLOK,2),ALINEMOM53(LX4,IBLOK,2) &
  & ,ALINEMOM54(LX4,IBLOK,2),ALINEMOM55(LX4,IBLOK,2) &
  & ,ALINEMOM56(LX4,IBLOK,2),ALINEMOM57(LX4,IBLOK,2) &
  & ,ALINEMOM58(LX4,IBLOK,2),ALINEMOM59(LX4,IBLOK,2) &
  & ,ALINEMOM60(LX4,IBLOK,2),ALINEMOM61(LX4,IBLOK,2) &
  & ,ALINEMOM62(LX4,IBLOK,2),ALINEMOM63(LX4,IBLOK,2) &
  & ,ALINEMOM64(LX4,IBLOK,2) &
  & ,ALINEMOM65(LX4,IBLOK,2),ALINEMOM66(LX4,IBLOK,2) &
  & ,ALINEMOM67(LX4,IBLOK,2),ALINEMOM68(LX4,IBLOK,2) &
  & ,ALINEMOM69(LX4,IBLOK,2),ALINEMOM70(LX4,IBLOK,2) &
  & ,ALINEMOM71(LX4,IBLOK,2),ALINEMOM72(LX4,IBLOK,2) &
  & ,ALINEMOM73(LX4,IBLOK,2),ALINEMOM74(LX4,IBLOK,2) &
  & ,ALINEMOM75(LX4,IBLOK,2),ALINEMOM76(LX4,IBLOK,2) &
  & ,ALINEMOM77(LX4,IBLOK,2),ALINEMOM78(LX4,IBLOK,2) &
  & ,ALINEMOM79(LX4,IBLOK,2),ALINEMOM80(LX4,IBLOK,2) &
  & ,ALINEMOM81(LX4,IBLOK,2),ALINEMOM82(LX4,IBLOK,2) &
  & ,ALINEMOM83(LX4,IBLOK,2),ALINEMOM84(LX4,IBLOK,2) &
  & ,ALINEMOM85(LX4,IBLOK,2),ALINEMOM86(LX4,IBLOK,2) &
  & ,ALINEMOM87(LX4,IBLOK,2),ALINEMOM88(LX4,IBLOK,2) &
  & ,ALINEMOM89(LX4,IBLOK,2),ALINEMOM90(LX4,IBLOK,2) &
  & ,ALINEMOM91(LX4,IBLOK,2),ALINEMOM92(LX4,IBLOK,2) &
  & ,ALINEMOM93(LX4,IBLOK,2),ALINEMOM94(LX4,IBLOK,2) &
  & ,ALINEMOM95(LX4,IBLOK,2),ALINEMOM96(LX4,IBLOK,2) &
  & ,ALINEMOM97(LX4,IBLOK,2),ALINEMOM98(LX4,IBLOK,2) &
  & ,ALINEMOM99(LX4,IBLOK,2),ALINEMOM100(LX4,IBLOK,2) &
  & ,ALINEMOM101(LX4,IBLOK,2),ALINEMOM102(LX4,IBLOK,2) &
  & ,ALINEMOM103(LX4,IBLOK,2),ALINEMOM104(LX4,IBLOK,2) &
  & ,ALINEMOM105(LX4,IBLOK,2),ALINEMOM106(LX4,IBLOK,2) &
  & ,ALINEMOM107(LX4,IBLOK,2),ALINEMOM108(LX4,IBLOK,2) &
  & ,ALINEMOM109(LX4,IBLOK,2),ALINEMOM110(LX4,IBLOK,2) &
  & ,ALINEMOM111(LX4,IBLOK,2),ALINEMOM112(LX4,IBLOK,2) &
  & ,ALINEMOM113(LX4,IBLOK,2),ALINEMOM114(LX4,IBLOK,2) &
  & ,ALINEMOM115(LX4,IBLOK,2),ALINEMOM116(LX4,IBLOK,2) &
  & ,ALINEMOM117(LX4,IBLOK,2),ALINEMOM118(LX4,IBLOK,2) &
  & ,ALINEMOM119(LX4,IBLOK,2),ALINEMOM120(LX4,IBLOK,2) &
  & ,ALINEMOM121(LX4,IBLOK,2),ALINEMOM122(LX4,IBLOK,2) &
  & ,ALINEMOM123(LX4,IBLOK,2),ALINEMOM124(LX4,IBLOK,2) &
  & ,ALINEMOM125(LX4,IBLOK,2),ALINEMOM126(LX4,IBLOK,2) &
  & ,ALINEMOM127(LX4,IBLOK,2),ALINEMOM128(LX4,IBLOK,2) &
  & ,ALINEMOM129(LX4,IBLOK,2),ALINEMOM130(LX4,IBLOK,2) &
  & ,ALINEMOM131(LX4,IBLOK,2),ALINEMOM132(LX4,IBLOK,2) &
  & ,ALINEMOM133(LX4,IBLOK,2),ALINEMOM134(LX4,IBLOK,2) &
  & ,ALINEMOM135(LX4,IBLOK,2),ALINEMOM136(LX4,IBLOK,2) &
  & ,ALINEMOM137(LX4,IBLOK,2),ALINEMOM138(LX4,IBLOK,2) &
  & ,ALINEMOM139(LX4,IBLOK,2),ALINEMOM140(LX4,IBLOK,2) &
  & ,ALINEMOM141(LX4,IBLOK,2),ALINEMOM142(LX4,IBLOK,2) &
  & ,ALINEMOM143(LX4,IBLOK,2),ALINEMOM144(LX4,IBLOK,2) &
  & ,ALINEMOM145(LX4,IBLOK,2),ALINEMOM146(LX4,IBLOK,2) &
  & ,ALINEMOM147(LX4,IBLOK,2),ALINEMOM148(LX4,IBLOK,2) &
  & ,ALINEMOM149(LX4,IBLOK,2),ALINEMOM150(LX4,IBLOK,2) &
  & ,ALINEMOM151(LX4,IBLOK,2),ALINEMOM152(LX4,IBLOK,2) &
  & ,ALINEMOM153(LX4,IBLOK,2),ALINEMOM154(LX4,IBLOK,2) &
  & ,ALINEMOM155(LX4,IBLOK,2),ALINEMOM156(LX4,IBLOK,2) &
  & ,ALINEMOM157(LX4,IBLOK,2),ALINEMOM158(LX4,IBLOK,2) &
  & ,ALINEMOM159(LX4,IBLOK,2),ALINEMOM160(LX4,IBLOK,2) &
  & ,ALINEMOM161(LX4,IBLOK,2),ALINEMOM162(LX4,IBLOK,2) &
  & ,ALINEMOM163(LX4,IBLOK,2),ALINEMOM164(LX4,IBLOK,2) &
  & ,ALINEMOM165(LX4,IBLOK,2),ALINEMOM166(LX4,IBLOK,2) &
  & ,ALINEMOM167(LX4,IBLOK,2),ALINEMOM168(LX4,IBLOK,2) &
  & ,ALINEMOM169(LX4,IBLOK,2),ALINEMOM170(LX4,IBLOK,2) &
  & ,ALINEMOM171(LX4,IBLOK,2),ALINEMOM172(LX4,IBLOK,2) &
  & ,ALINEMOM173(LX4,IBLOK,2),ALINEMOM174(LX4,IBLOK,2) &
  & ,ALINEMOM175(LX4,IBLOK,2),ALINEMOM176(LX4,IBLOK,2) &
  & ,ALINEMOM177(LX4,IBLOK,2),ALINEMOM178(LX4,IBLOK,2) &
  & ,ALINEMOM179(LX4,IBLOK,2),ALINEMOM180(LX4,IBLOK,2) &
  & ,ALINEMOM181(LX4,IBLOK,2),ALINEMOM182(LX4,IBLOK,2) &
  & ,ALINEMOM183(LX4,IBLOK,2),ALINEMOM184(LX4,IBLOK,2) &
  & ,ALINEMOM185(LX4,IBLOK,2),ALINEMOM186(LX4,IBLOK,2) &
  & ,ALINEMOM187(LX4,IBLOK,2),ALINEMOM188(LX4,IBLOK,2) &
  & ,ALINEMOM189(LX4,IBLOK,2),ALINEMOM190(LX4,IBLOK,2) &
  & ,ALINEMOM191(LX4,IBLOK,2),ALINEMOM192(LX4,IBLOK,2) &
  & ,ALINEMOM193(LX4,IBLOK,2),ALINEMOM194(LX4,IBLOK,2) &
  & ,ALINEMOM195(LX4,IBLOK,2),ALINEMOM196(LX4,IBLOK,2) &
  & ,ALINEMOM197(LX4,IBLOK,2),ALINEMOM198(LX4,IBLOK,2) &
  & ,ALINEMOM199(LX4,IBLOK,2),ALINEMOM200(LX4,IBLOK,2) &
  & ,ALINEMOM201(LX4,IBLOK,2),ALINEMOM202(LX4,IBLOK,2) &
  & ,ALINEMOM203(LX4,IBLOK,2),ALINEMOM204(LX4,IBLOK,2) &
  & ,ALINEMOM205(LX4,IBLOK,2),ALINEMOM206(LX4,IBLOK,2) &
  & ,ALINEMOM207(LX4,IBLOK,2),ALINEMOM208(LX4,IBLOK,2) &
  & ,ALINEMOM209(LX4,IBLOK,2),ALINEMOM210(LX4,IBLOK,2) &
  & ,ALINEMOM211(LX4,IBLOK,2),ALINEMOM212(LX4,IBLOK,2) &
  & ,ALINEMOM213(LX4,IBLOK,2),ALINEMOM214(LX4,IBLOK,2) &
  & ,ALINEMOM215(LX4,IBLOK,2),ALINEMOM216(LX4,IBLOK,2) &
  & ,ALINEMOM217(LX4,IBLOK,2),ALINEMOM218(LX4,IBLOK,2) &
  & ,ALINEMOM219(LX4,IBLOK,2),ALINEMOM220(LX4,IBLOK,2) &
  & ,ALINEMOM221(LX4,IBLOK,2),ALINEMOM222(LX4,IBLOK,2) &
  & ,ALINEMOM223(LX4,IBLOK,2),ALINEMOM224(LX4,IBLOK,2) &
  & ,ALINEMOM225(LX4,IBLOK,2),ALINEMOM226(LX4,IBLOK,2) &
  & ,ALINEMOM227(LX4,IBLOK,2),ALINEMOM228(LX4,IBLOK,2) &
  & ,ALINEMOM229(LX4,IBLOK,2),ALINEMOM230(LX4,IBLOK,2) &
  & ,ALINEMOM231(LX4,IBLOK,2),ALINEMOM232(LX4,IBLOK,2) &
  & ,ALINEMOM233(LX4,IBLOK,2),ALINEMOM234(LX4,IBLOK,2) &
  & ,ALINEMOM235(LX4,IBLOK,2)
!********************************************************
COMMON/ARRAYB/UB11(NCOL2,LSIZEB,3,IBLOK)
COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
COMMON/DUMMY/UC11(NCOL2,LSIZEB,3) &
  & ,IUP(LSIZEB,3),IDN(LSIZEB,3)
DIMENSION LCNT(IBLOK),LB(IBLOK),IX(3),LS(3)
DIMENSION ACT(3),AST(3)
DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2),E11(NCOL2) &
  & ,F11(NCOL2),UINT11(NCOL2),DUM11(NCOL2)
DIMENSION REM11(NCOL2)
!*********************************************************************
!     HERE I DEFINE THE SQUARE PULSES
!*********************************************************************
DIMENSION SQUY1(NCOL2),SQDY1(NCOL2)
DIMENSION SQUZ1(NCOL2),SQDZ1(NCOL2)
DIMENSION SQUY2(NCOL2),SQDY2(NCOL2)
DIMENSION SQUZ2(NCOL2),SQDZ2(NCOL2)
DIMENSION SQUDY1(NCOL2),SQUDZ1(NCOL2)  !NEW!
DIMENSION SQDUY1(NCOL2),SQDUZ1(NCOL2)  !NEW!
!*********************************************************************
!     HERE I DEFINE THE WAVE-LIKE SQUARE PULSES
!*********************************************************************
DIMENSION WSQUY1(NCOL2),WSQUY2(NCOL2),WSQUY3(NCOL2),WSQUY4(NCOL2)
DIMENSION WSQUZ1(NCOL2),WSQUZ2(NCOL2),WSQUZ3(NCOL2),WSQUZ4(NCOL2)
DIMENSION WSQDY1(NCOL2),WSQDY2(NCOL2),WSQDY3(NCOL2),WSQDY4(NCOL2)
DIMENSION WSQDZ1(NCOL2),WSQDZ2(NCOL2),WSQDZ3(NCOL2),WSQDZ4(NCOL2)
!*********************************************************************
!     HERE I DEFINE SOME TT-OPERATORS COMPONENTS
!*********************************************************************
DIMENSION DUY1(NCOL2),DUZ1(NCOL2),DDY1(NCOL2),DDZ1(NCOL2)
DIMENSION DUY2(NCOL2),DUZ2(NCOL2),DDY2(NCOL2),DDZ2(NCOL2)
!*********************************************************************
!     HERE I DEFINE SOME PLAQUETTE OPERATOR COMPONENTS
!*********************************************************************
DIMENSION PQ1(NCOL2),PQ2(NCOL2),PQ3(NCOL2),PQ4(NCOL2)
DIMENSION PQ5(NCOL2),PQ6(NCOL2),PQ7(NCOL2),PQ8(NCOL2)
!*********************************************************************
DIMENSION G11(NCOL2)
!*********************************************************************
DIMENSION WVUY1(NCOL2),WVUY2(NCOL2),WVDY1(NCOL2),WVDY2(NCOL2)
DIMENSION WVUZ1(NCOL2),WVUZ2(NCOL2),WVDZ1(NCOL2),WVDZ2(NCOL2)
!*********************************************************************
DIMENSION WANGL1(NCOL2),WANGL2(NCOL2),WANGL3(NCOL2),WANGL4(NCOL2) !new!
DIMENSION WANGL5(NCOL2),WANGL6(NCOL2),WANGL7(NCOL2),WANGL8(NCOL2) !new!
!*********************************************************************
DIMENSION ZIG1W1(NCOL2),ZIG2W1(NCOL2),ZIG1W2(NCOL2),ZIG2W2(NCOL2)
DIMENSION ZIG1W3(NCOL2),ZIG2W3(NCOL2),ZIG1W4(NCOL2),ZIG2W4(NCOL2)
DIMENSION TIG1W1(NCOL2),TIG2W1(NCOL2),TIG1W2(NCOL2),TIG2W2(NCOL2)
DIMENSION TIG1W3(NCOL2),TIG2W3(NCOL2),TIG1W4(NCOL2),TIG2W4(NCOL2)
!*********************************************************************
DIMENSION LIN0(NCOL2),LIN1(NCOL2),LIN2(NCOL2),LIN4(NCOL2)
DIMENSION PLQ1(NCOL2),PLQ2(NCOL2),PLQ3(NCOL2),PLQ4(NCOL2)
DIMENSION PLQ5(NCOL2),PLQ6(NCOL2),PLQ7(NCOL2),PLQ8(NCOL2)
DIMENSION DPLQ1(NCOL2),DPLQ2(NCOL2),DPLQ3(NCOL2),DPLQ4(NCOL2)
DIMENSION DPLQ5(NCOL2),DPLQ6(NCOL2),DPLQ7(NCOL2),DPLQ8(NCOL2)
DIMENSION PL(NCOL2)
!*********************************************************************
complex(real32) REM11
complex(real32) UB11,A11,B11,C11,D11,E11,F11,UC11,UINT11,DUM11
complex(real32) SQUY1,SQDY1,SQUY2,SQDY2,G11
complex(real32) SQUZ1,SQDZ1,SQUZ2,SQDZ2
complex(real32) SQUDY1,SQUDZ1,SQDUY1,SQDUZ1  !new!
complex(real32) WSQUY1,WSQUY2,WSQUY3,WSQUY4
complex(real32) WSQUZ1,WSQUZ2,WSQUZ3,WSQUZ4
complex(real32) WSQDY1,WSQDY2,WSQDY3,WSQDY4
complex(real32) WSQDZ1,WSQDZ2,WSQDZ3,WSQDZ4
complex(real32) WANGL1,WANGL2,WANGL3,WANGL4 !new!
complex(real32) WANGL5,WANGL6,WANGL7,WANGL8 !new!
complex(real32) WVUY1,WVUY2,WVDY1,WVDY2
complex(real32) WVUZ1,WVUZ2,WVDZ1,WVDZ2
complex(real32) ZIG1W1,ZIG2W1,ZIG1W2,ZIG2W2
complex(real32) ZIG1W3,ZIG2W3,ZIG1W4,ZIG2W4
complex(real32) TIG1W1,TIG2W1,TIG1W2,TIG2W2
complex(real32) TIG1W3,TIG2W3,TIG1W4,TIG2W4
complex(real32) LIN0,LIN1,LIN2,LIN4
!     Plaquette operators on first cite
complex(real32) PLQ1,PLQ2,PLQ3,PLQ4
complex(real32) PLQ5,PLQ6,PLQ7,PLQ8
!
complex(real32) DPLQ1,DPLQ2,DPLQ3,DPLQ4
complex(real32) DPLQ5,DPLQ6,DPLQ7,DPLQ8
!      
complex(real32) PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,PQ8
complex(real32) DUY1,DUZ1,DDY1,DDZ1,DUY2,DUZ2,DDY2,DDZ2
complex(real32) PL
!**********************************************************************
complex(real64) CSUMN,CSUMS(4),CSUM2S(4),CSUM2WS(4)
complex(real64) CSUMW(4),CSUM2W(4),CSUM3W(4),CSUMUP(4),CSUMUD(4)
complex(real64) CSUMRC(4),CSUMRCW(4)
complex(real64) CSUMTT1(8),CSUMTT2(8),CSUMTT3(8),CSUMTT4(8),CSUMTT5(8)
complex(real64) CSUMTT6(8),CSUMTT7(8),CSUMTT8(8),CSUMTT9(8),CSUMTT10(8)
complex(real64) CSUMTT11(8),CSUMTT12(16),CSUMTT13(8),CSUMTT14(8)
!
complex(real64) CSUMPLQ(8),CSUMPLQ2(8),CSUMPLQ3(8),CSUMPLQ4(8)
complex(real64) CSUMPLQ5(8),CSUMPLQ6(8)
complex(real64) CSUMPLQ7(16),CSUMPLQ8(16),CSUMPLQ9(16),CSUMPLQ10(16)
complex(real64) CSUMPLQ11(16),CSUMPLQ12(16),CSUMPLQ13(16),CSUMPLQ14(16)
complex(real64) CSUMPLQ15(16)
!     
complex(real64) AKT1,ACT
!**********************************************************************
complex(real64) ALINE1,ALINE2,ALINE3,ALINE4,ALINE5,ALINE6,ALINE7
complex(real64) ALINE8,ALINE9,ALINE10,ALINE11,ALINE12,ALINE13
complex(real64) ALINE14,ALINE15,ALINE16,ALINE17,ALINE18,ALINE19
complex(real64) ALINE20,ALINE21,ALINE22,ALINE23,ALINE24,ALINE25
complex(real64) ALINE26,ALINE27,ALINE28,ALINE29,ALINE30
complex(real64) ALINE31,ALINE32,ALINE33,ALINE34,ALINE35
complex(real64) ALINE36,ALINE37,ALINE38,ALINE39,ALINE40
complex(real64) ALINE41,ALINE42,ALINE43,ALINE44,ALINE45
complex(real64) ALINE46,ALINE47,ALINE48,ALINE49,ALINE50
complex(real64) ALINE51,ALINE52,ALINE53,ALINE54,ALINE55
complex(real64) ALINE56,ALINE57,ALINE58,ALINE59,ALINE60
complex(real64) ALINE61,ALINE62,ALINE63,ALINE64,ALINE65
complex(real64) ALINE66,ALINE67,ALINE68,ALINE69,ALINE70
complex(real64) ALINE71,ALINE72,ALINE73,ALINE74,ALINE75
complex(real64) ALINE76,ALINE77,ALINE78,ALINE79,ALINE80
complex(real64) ALINE81,ALINE82,ALINE83,ALINE84,ALINE85
complex(real64) ALINE86,ALINE87,ALINE88,ALINE89,ALINE90
complex(real64) ALINE91,ALINE92,ALINE93,ALINE94,ALINE95
complex(real64) ALINE96,ALINE97,ALINE98,ALINE99,ALINE100
complex(real64) ALINE101,ALINE102,ALINE103,ALINE104,ALINE105
complex(real64) ALINE106,ALINE107,ALINE108,ALINE109,ALINE110
complex(real64) ALINE111,ALINE112,ALINE113,ALINE114,ALINE115
complex(real64) ALINE116,ALINE117,ALINE118,ALINE119,ALINE120
complex(real64) ALINE121,ALINE122,ALINE123,ALINE124,ALINE125
complex(real64) ALINE126,ALINE127,ALINE128,ALINE129,ALINE130
complex(real64) ALINE131,ALINE132,ALINE133,ALINE134,ALINE135
complex(real64) ALINE136,ALINE137,ALINE138,ALINE139,ALINE140
complex(real64) ALINE141,ALINE142,ALINE143,ALINE144,ALINE145
complex(real64) ALINE146,ALINE147,ALINE148,ALINE149,ALINE150
complex(real64) ALINE151,ALINE152,ALINE153,ALINE154,ALINE155
complex(real64) ALINE156,ALINE157,ALINE158,ALINE159,ALINE160
complex(real64) ALINE161,ALINE162,ALINE163,ALINE164,ALINE165
complex(real64) ALINE166,ALINE167,ALINE168,ALINE169,ALINE170
complex(real64) ALINE171,ALINE172,ALINE173,ALINE174,ALINE175
complex(real64) ALINE176,ALINE177,ALINE178,ALINE179,ALINE180
complex(real64) ALINE181,ALINE182,ALINE183,ALINE184,ALINE185
complex(real64) ALINE186,ALINE187,ALINE188,ALINE189,ALINE190
complex(real64) ALINE191,ALINE192,ALINE193,ALINE194,ALINE195
complex(real64) ALINE196,ALINE197,ALINE198,ALINE199,ALINE200
complex(real64) ALINE201,ALINE202,ALINE203,ALINE204,ALINE205
complex(real64) ALINE206,ALINE207,ALINE208,ALINE209,ALINE210
complex(real64) ALINE211,ALINE212,ALINE213,ALINE214,ALINE215
complex(real64) ALINE216,ALINE217,ALINE218,ALINE219,ALINE220
complex(real64) ALINE221,ALINE222,ALINE223,ALINE224,ALINE225
complex(real64) ALINE226,ALINE227,ALINE228,ALINE229,ALINE230
complex(real64) ALINE231,ALINE232,ALINE233,ALINE234,ALINE235
!**********************************************************************
complex(real64) ALINEMOM2,ALINEMOM3,ALINEMOM4
complex(real64) ALINEMOM5,ALINEMOM6,ALINEMOM7,ALINEMOM8
complex(real64) ALINEMOM9,ALINEMOM10,ALINEMOM11,ALINEMOM12
complex(real64) ALINEMOM13,ALINEMOM14,ALINEMOM15,ALINEMOM16
complex(real64) ALINEMOM17,ALINEMOM18,ALINEMOM19,ALINEMOM20
complex(real64) ALINEMOM21,ALINEMOM22,ALINEMOM23,ALINEMOM24
complex(real64) ALINEMOM25,ALINEMOM26,ALINEMOM27,ALINEMOM28
complex(real64) ALINEMOM29,ALINEMOM30,ALINEMOM31,ALINEMOM32
complex(real64) ALINEMOM33,ALINEMOM34,ALINEMOM35,ALINEMOM36
complex(real64) ALINEMOM37,ALINEMOM38,ALINEMOM39,ALINEMOM40
complex(real64) ALINEMOM41,ALINEMOM42,ALINEMOM43,ALINEMOM44
complex(real64) ALINEMOM45,ALINEMOM46,ALINEMOM47,ALINEMOM48
complex(real64) ALINEMOM49,ALINEMOM50,ALINEMOM51,ALINEMOM52
complex(real64) ALINEMOM53,ALINEMOM54,ALINEMOM55,ALINEMOM56
complex(real64) ALINEMOM57,ALINEMOM58,ALINEMOM59,ALINEMOM60
complex(real64) ALINEMOM61,ALINEMOM62,ALINEMOM63,ALINEMOM64
complex(real64) ALINEMOM65,ALINEMOM66,ALINEMOM67,ALINEMOM68
complex(real64) ALINEMOM69,ALINEMOM70,ALINEMOM71,ALINEMOM72
complex(real64) ALINEMOM73,ALINEMOM74,ALINEMOM75,ALINEMOM76
complex(real64) ALINEMOM77,ALINEMOM78,ALINEMOM79,ALINEMOM80
complex(real64) ALINEMOM81,ALINEMOM82,ALINEMOM83,ALINEMOM84
complex(real64) ALINEMOM85,ALINEMOM86,ALINEMOM87,ALINEMOM88
complex(real64) ALINEMOM89,ALINEMOM90,ALINEMOM91,ALINEMOM92
complex(real64) ALINEMOM93,ALINEMOM94,ALINEMOM95,ALINEMOM96
complex(real64) ALINEMOM97,ALINEMOM98,ALINEMOM99,ALINEMOM100
complex(real64) ALINEMOM101,ALINEMOM102,ALINEMOM103,ALINEMOM104
complex(real64) ALINEMOM105,ALINEMOM106,ALINEMOM107,ALINEMOM108
complex(real64) ALINEMOM109,ALINEMOM110
complex(real64) ALINEMOM111,ALINEMOM112,ALINEMOM113,ALINEMOM114
complex(real64) ALINEMOM115,ALINEMOM116,ALINEMOM117,ALINEMOM118
complex(real64) ALINEMOM119,ALINEMOM120,ALINEMOM121,ALINEMOM122
complex(real64) ALINEMOM123,ALINEMOM124,ALINEMOM125,ALINEMOM126
complex(real64) ALINEMOM127,ALINEMOM128,ALINEMOM129,ALINEMOM130
complex(real64) ALINEMOM131,ALINEMOM132,ALINEMOM133,ALINEMOM134
complex(real64) ALINEMOM135,ALINEMOM136,ALINEMOM137,ALINEMOM138
complex(real64) ALINEMOM139,ALINEMOM140,ALINEMOM141,ALINEMOM142
complex(real64) ALINEMOM143,ALINEMOM144,ALINEMOM145,ALINEMOM146
complex(real64) ALINEMOM147,ALINEMOM148,ALINEMOM149,ALINEMOM150
complex(real64) ALINEMOM151,ALINEMOM152,ALINEMOM153,ALINEMOM154
complex(real64) ALINEMOM155,ALINEMOM156,ALINEMOM157,ALINEMOM158
complex(real64) ALINEMOM159,ALINEMOM160,ALINEMOM161,ALINEMOM162
complex(real64) ALINEMOM163,ALINEMOM164,ALINEMOM165,ALINEMOM166
complex(real64) ALINEMOM167,ALINEMOM168,ALINEMOM169,ALINEMOM170
complex(real64) ALINEMOM171,ALINEMOM172,ALINEMOM173,ALINEMOM174
complex(real64) ALINEMOM175,ALINEMOM176,ALINEMOM177,ALINEMOM178
complex(real64) ALINEMOM179,ALINEMOM180,ALINEMOM181,ALINEMOM182
complex(real64) ALINEMOM183,ALINEMOM184,ALINEMOM185,ALINEMOM186
complex(real64) ALINEMOM187,ALINEMOM188,ALINEMOM189,ALINEMOM190
complex(real64) ALINEMOM191,ALINEMOM192,ALINEMOM193,ALINEMOM194
complex(real64) ALINEMOM195,ALINEMOM196,ALINEMOM197,ALINEMOM198
complex(real64) ALINEMOM199,ALINEMOM200,ALINEMOM201,ALINEMOM202
complex(real64) ALINEMOM203,ALINEMOM204,ALINEMOM205,ALINEMOM206
complex(real64) ALINEMOM207,ALINEMOM208,ALINEMOM209,ALINEMOM210
complex(real64) ALINEMOM211,ALINEMOM212,ALINEMOM213,ALINEMOM214
complex(real64) ALINEMOM215,ALINEMOM216,ALINEMOM217,ALINEMOM218
complex(real64) ALINEMOM219,ALINEMOM220,ALINEMOM221,ALINEMOM222
complex(real64) ALINEMOM223,ALINEMOM224,ALINEMOM225,ALINEMOM226
complex(real64) ALINEMOM227,ALINEMOM228,ALINEMOM229,ALINEMOM230
complex(real64) ALINEMOM231,ALINEMOM232,ALINEMOM233,ALINEMOM234
complex(real64) ALINEMOM235
!**********************************************************************
complex(real64) CSUMSMOM(4,2),CSUM2SMOM(4,2),CSUM2WSMOM(4,2)
complex(real64) CSUMWMOM(4,2),CSUM2WMOM(4,2),CSUM3WMOM(4,2)
complex(real64) CSUMUPMOM(4,2),CSUMUDMOM(4,2)
!
complex(real64) CSUMPLQMOM(8,2),CSUMPLQMOM2(8,2),CSUMPLQMOM3(8,2)
complex(real64) CSUMPLQMOM4(8,2),CSUMPLQMOM5(8,2),CSUMPLQMOM6(8,2)
complex(real64) CSUMPLQMOM7(16,2), CSUMPLQMOM8(16,2)
complex(real64) CSUMPLQMOM9(16,2), CSUMPLQMOM10(16,2)
complex(real64) CSUMPLQMOM11(16,2), CSUMPLQMOM12(16,2)
complex(real64) CSUMPLQMOM13(16,2), CSUMPLQMOM14(16,2)
complex(real64) CSUMPLQMOM15(16,2)
!     
complex(real64) CSUMTTMOM1(8,2),CSUMTTMOM2(8,2),CSUMTTMOM3(8,2)
complex(real64) CSUMTTMOM4(8,2),CSUMTTMOM5(8,2),CSUMTTMOM6(8,2)
complex(real64) CSUMTTMOM7(8,2),CSUMTTMOM8(8,2),CSUMTTMOM9(8,2)
complex(real64) CSUMTTMOM10(8,2),CSUMTTMOM11(8,2)
complex(real64) CSUMTTMOM12(16,2)
complex(real64) CSUMTTMOM13(8,2),CSUMTTMOM14(8,2)
!
complex(real64) GIOT,PF(2)
complex(real32) UREN11(NCOL2)
complex(real32) AA(NCOL2)
complex(real32) DET,CNORM,CDET,DNCOL
!**********************************************************************
GIOT=(0.0,1.0)
PI=4.0d0*DATAN(1.0d0)
LS(1)=LX1
LS(2)=LX2
LS(3)=LX3
LB(1)=1
!*********************************************************************     
do IDD=2,IBLOK
LB(IDD)=2*LB(IDD-1)
1 CONTINUE
end do
!********************************************************************* 
do KU=1,1
JU=KU+1
IF(JU.GT.3)JU=JU-3
IU=KU+2
IF(IU.GT.3)IU=IU-3
!     
ID=IBLL
!
do IB=1,ID
LCNT(IB)=0
3 CONTINUE
end do
!     
LREST=LS(KU)
!     
do IG=1,ID
IDG=ID-IG+1
LCNT(IDG)=LREST/LB(IDG)
LREST=LREST-LCNT(IDG)*LB(IDG)
4 CONTINUE
end do
!     
do IB=1,ID
IF(LCNT(IB).GE.1)IDS=IB
5 CONTINUE
end do
!     
IDSM1=IDS-1
IF(IDS.EQ.1)IDSM1=IDS
!**********************************************************************
DO KK=1,3
DO NN=1,LSIZEB
IUP(NN,KK)=IUPB(NN,KK,IDS)
IDN(NN,KK)=IDNB(NN,KK,IDS)
end do
end do
!**********************************************************************    
DO KK=1,3
DO NN=1,LSIZEB
DO IC=1,NCOL2
UC11(IC,NN,KK)=UB11(IC,NN,KK,IDS)
end do
end do
end do
!**********************************************************************
CSUMN=(0.0,0.0)
!
DO IJN=1,4
CSUMS(IJN)=(0.0,0.0)
CSUM2S(IJN)=(0.0,0.0)
CSUM2WS(IJN)=(0.0,0.0)
CSUMW(IJN)=(0.0,0.0)
CSUM2W(IJN)=(0.0,0.0)
CSUM3W(IJN)=(0.0,0.0)
CSUMUP(IJN)=(0.0,0.0)
CSUMUD(IJN)=(0.0,0.0)
end do
!
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
end do
!
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
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
end do
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CSUMPLQMOM(IJN,IJ)=(0.0,0.0)
CSUMPLQMOM2(IJN,IJ)=(0.0,0.0)
CSUMPLQMOM3(IJN,IJ)=(0.0,0.0)
CSUMPLQMOM4(IJN,IJ)=(0.0,0.0)
CSUMPLQMOM5(IJN,IJ)=(0.0,0.0)
CSUMPLQMOM6(IJN,IJ)=(0.0,0.0)
end do
end do
!
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
end do
end do
!
!**********************************************************************
NN=0
!**********************************************************************
!     OPENMP INTRINSICS
!     NK,LS,PF,KU
!     UC11 
!     A11, B11, C11, D11
!     SQUY1, SQUZ1, SQDY1, SQDZ1


DO NK=1,LS(KU)
PF(1)=DCOS((2.0*PI*NK)/LS(KU))+GIOT*DSIN((2.0*PI*NK)/LS(KU))
PF(2)=DCOS((4.0*PI*NK)/LS(KU))+GIOT*DSIN((4.0*PI*NK)/LS(KU))
!       PF(3)=DCOS((6.0*PI*NK)/LS(KU))+GIOT*DSIN((6.0*PI*NK)/LS(KU))
!       PF(4)=DCOS((8.0*PI*NK)/LS(KU))+GIOT*DSIN((8.0*PI*NK)/LS(KU))
DO NJ=1,LS(JU)
DO NI=1,LS(IU)
NN=NN+1
IX(IU)=NI
IX(JU)=NJ
IX(KU)=NK
MN=IX(1)+LS(1)*(IX(2)-1)+LS(1)*LS(2)*(IX(3)-1)
!**********************************************************************
IF(IDS.EQ.ID) THEN
!**********************************************************************
IREM=3
IF(IDS.EQ.1) IREM=4
do ILOOP=1,IREM
IF(ILOOP.EQ.1)LI=LCNT(IDS)
IF(ILOOP.GT.1)LI=LCNT(IDS)-2**(ILOOP-2)
IF(LI.LT.0) go to 201
!
M2=MN
IF(ILOOP.GT.1) THEN
DO I=1, 2**(ILOOP-2)
M3=IUP(M2,KU)
M2=M3
END DO
end if
!
do IC=1,NCOL2
E11(IC)=(0.0,0.0)
202 CONTINUE
end do
!
do N1=1,NCOL
IJ=N1+NCOL*(N1-1)
E11(IJ)=(1.0,0.0)
203 CONTINUE
end do
!
do NC=1,LI
do IC=1,NCOL2
B11(IC)=UC11(IC,M2,KU)
205 CONTINUE
end do
CALL VMX(1,E11,B11,C11,1)
M3=IUP(M2,KU)
M2=M3
do IC=1,NCOL2
E11(IC)=C11(IC)
206 CONTINUE
end do
204 CONTINUE
end do
!
ML=M2
IF(ILOOP.EQ.1) THEN
do IC=1,NCOL2
LIN0(IC)=E11(IC)
207 CONTINUE
end do
end if
!********************************************************************                  
IF(ILOOP.EQ.2) THEN
do IC=1,NCOL2
LIN1(IC)=E11(IC)
208 CONTINUE
end do
end if
!********************************************************************                                 
IF(ILOOP.EQ.3) THEN
do IC=1,NCOL2
LIN2(IC)=E11(IC)
209 CONTINUE
end do
end if
!********************************************************************
IF(ILOOP.EQ.4) THEN
do IC=1,NCOL2
LIN4(IC)=E11(IC)
210 CONTINUE
end do
end if
!********************************************************************
201 CONTINUE
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ELSE
M2=MN
end if
!**********************************************************************
!*********************** REMAINING PIECE ******************************
!**********************************************************************            
DO IC=1, NCOL2
REM11(IC)=(0.0,0.0)
end do
!********************************************************************
DO N1=1,NCOL
IJ=N1+NCOL*(N1-1)
REM11(IJ)=(1.0,0.0)
end do
!********************************************************************
do IG=1,IDS
IDG=IDS-IG+1
IF(IDG.EQ.ID)go to 176
do NC=1,LCNT(IDG)
do IC=1,NCOL2
UINT11(IC)=(0.0,0.0)
178 CONTINUE
end do
!     
do IJ=1,3
IF(IJ.EQ.KU)go to 179
do IC=1,NCOL2
B11(IC)=UB11(IC,M2,IJ,IDSM1)
180 CONTINUE
end do
M3=IUPB(M2,IJ,IDSM1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDG)
181 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M1=IUPB(M2,KU,IDG)
do IC=1,NCOL2
B11(IC)=UB11(IC,M1,IJ,IDSM1)
182 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,D11,B11,C11,1)
do IC=1,NCOL2
UINT11(IC)=UINT11(IC)+C11(IC)
183 CONTINUE
end do
M3=IDNB(M2,IJ,IDSM1)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,IJ,IDSM1)
184 CONTINUE
end do
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDG)
185 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,B11,C11,D11,1)
M1=IUPB(M3,KU,IDG)
do IC=1,NCOL2
B11(IC)=UB11(IC,M1,IJ,IDSM1)
186 CONTINUE
end do
CALL VMX(1,D11,B11,C11,1)
do IC=1,NCOL2
UINT11(IC)=UINT11(IC)+C11(IC)
187 CONTINUE
end do
179 CONTINUE
end do
do IC=1,NCOL2
B11(IC)=UB11(IC,M2,KU,IDG)
188 CONTINUE
end do
do IC=1,NCOL2
UINT11(IC)=UINT11(IC)+B11(IC)
189 CONTINUE
end do
!new
CALL RENORMBS(UINT11,UREN11)
do IJ=1,NCOL2
AA(IJ)=UREN11(IJ)
441 CONTINUE
end do
JMAT=NCOL
CALL DETNANT(JMAT,DET,AA)
!
DNCOL=CMPLX(1.0/NCOL)
CDET=DET**DNCOL
CNORM=1.0/CDET
!
do IC=1,NCOL2
B11(IC)=UREN11(IC)*CNORM
644 CONTINUE
end do
!********************************************
CALL VMX(1,REM11,B11,C11,1)
M1=M2
do IC=1,NCOL2
REM11(IC)=C11(IC)
191 CONTINUE
end do
M2=IUPB(M1,KU,IDG)
177 CONTINUE
end do
176 CONTINUE
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
do IDDD=1,337 !new!
!
M2=MN
!
do IJ=1,NCOL2
A11(IJ)=(0.0,0.0)
10 CONTINUE
end do
!
do N1=1,NCOL
IJ=N1+NCOL*(N1-1)
A11(IJ)=(1.0,0.0)
11 CONTINUE
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
IF (IDS.EQ.ID) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
!**********************************************************************C
IF(LCNT(IDS).GE.ICO) THEN
!**********************************************************************C
!                     UP SQUARE PULSE                                  C 
!**********************************************************************C
!                     UP Y
!**********************************************************************
IF(IDDD.EQ.1) THEN
!
do IC=1,NCOL2
B11(IC)=UC11(IC,M2,JU)
12 CONTINUE
end do
M3=IUP(M2,JU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M3,KU)
13 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUP(M2,KU)
M2=M3
do IC=1,NCOL2
C11(IC)=UC11(IC,M2,JU)
14 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,SQUY1,1)
CALL VMX(1,A11,SQUY1,C11,1)
!     
do IC=1,NCOL2
A11(IC)=C11(IC)
15 CONTINUE
end do
!     
IEEE=1
!     
end if
!**********************************************************************C
!     UP Z                                                             C 
!**********************************************************************C
IF(IDDD.EQ.2) THEN
!     
do IC=1,NCOL2
B11(IC)=UC11(IC,M2,IU)
16 CONTINUE
end do
M3=IUP(M2,IU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M3,KU)
17 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUP(M2,KU)
M2=M3
do IC=1,NCOL2
C11(IC)=UC11(IC,M2,IU)
18 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,SQUZ1,1)
CALL VMX(1,A11,SQUZ1,C11,1)
!     
do IC=1,NCOL2
A11(IC)=C11(IC)
19 CONTINUE
end do
!     
IEEE=2
!     
end if
!**********************************************************************C 
!                       DOWN Y                                         C
!**********************************************************************C
IF(IDDD.EQ.3) THEN
!     
M3=IDN(M2,JU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M3,JU)
20 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UC11(IC,M3,KU)
21 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M4=IUP(M3,KU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M4,JU)
22 CONTINUE
end do
CALL VMX(1,B11,D11,SQDY1,1)
CALL VMX(1,A11,SQDY1,C11,1)
!     
do IC=1,NCOL2
A11(IC)=C11(IC)
23 CONTINUE
end do
!     
IEEE=3
!     
end if
!**********************************************************************C 
!                       DOWN Z                                 C
!**********************************************************************C
IF(IDDD.EQ.4) THEN
!     
M3=IDN(M2,IU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M3,IU)
24 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UC11(IC,M3,KU)
25 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M4=IUP(M3,KU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M4,IU)
26 CONTINUE
end do
CALL VMX(1,B11,D11,SQDZ1,1)
CALL VMX(1,A11,SQDZ1,C11,1)
!     
do IC=1,NCOL2
A11(IC)=C11(IC)
27 CONTINUE
end do
!     
IEEE=4
!     
end if
!**********************************************************************C
!                       UP - UP SQUARE PULSES                         *C 
!**********************************************************************C
!                       UP Y
!**********************************************************************C
IF(IDDD.EQ.5) THEN
!     
M3=IUP(M2,KU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M3,JU)
28 CONTINUE
end do
M4=IUP(M3,JU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M4,KU)
29 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M2=IUP(M3,KU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M2,JU)
30 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,SQUY2,1)
CALL VMX(1,SQUY1,SQUY2,A11,1)
!     
IEEE=1
!     
end if
!**********************************************************************C
!              UP Z
!**********************************************************************C
IF(IDDD.EQ.6) THEN
!     
M3=IUP(M2,KU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M3,IU)
31 CONTINUE
end do
M4=IUP(M3,IU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M4,KU)
32 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M2=IUP(M3,KU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M2,IU)
33 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,SQUZ2,1)
CALL VMX(1,SQUZ1,SQUZ2,A11,1)
!     
IEEE=2
!     
end if
!**********************************************************************C 
!              DOWN Y                                              C
!**********************************************************************C
IF(IDDD.EQ.7) THEN
!     
M3=IUP(M2,KU)
M4=IDN(M3,JU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M4,JU)
34 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UC11(IC,M4,KU)
35 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M1=IUP(M4,KU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M1,JU)
36 CONTINUE
end do
CALL VMX(1,B11,C11,SQDY2,1)
CALL VMX(1,SQDY1,SQDY2,A11,1)
M2=IUP(M3,KU)
!     
IEEE=3
!     
end if
!**********************************************************************C 
!              DOWN Z                                              C
!**********************************************************************C
IF(IDDD.EQ.8) THEN
!     
M3=IUP(M2,KU)
M4=IDN(M3,IU)
do IC=1,NCOL2
D11(IC)=UC11(IC,M4,IU)
37 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UC11(IC,M4,KU)
38 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M1=IUP(M4,KU)
do IC=1,NCOL2
C11(IC)=UC11(IC,M1,IU)
39 CONTINUE
end do
CALL VMX(1,B11,C11,SQDZ2,1)
CALL VMX(1,SQDZ1,SQDZ2,A11,1)
M2=IUP(M3,KU)
!     
IEEE=4
!     
end if
!**********************************************************************C
!**********************************************************************C
!              UP - DOWN SQUARE PULSES                                 C 
!**********************************************************************C
!**********************************************************************C
!                 UP Y 
!**********************************************************************C
IF(IDDD.EQ.9) THEN
CALL VMX(1,SQUY1,SQDY2,SQUDY1,1)
DO IC=1, NCOL2
A11(IC)=SQUDY1(IC) !NEW!
end do
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=1
!     
end if
!**********************************************************************C
!                 UP Z
!**********************************************************************C
IF(IDDD.EQ.10) THEN
CALL VMX(1,SQUZ1,SQDZ2,SQUDZ1,1)
DO IC=1, NCOL2
A11(IC)=SQUDZ1(IC) !NEW!
end do
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=2
!     
end if
!**********************************************************************C 
!                 DOWN Y
!**********************************************************************C
IF(IDDD.EQ.11) THEN
CALL VMX(1,SQDY1,SQUY2,SQDUY1,1)
DO IC=1, NCOL2
A11(IC)=SQDUY1(IC) !NEW!
end do
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=3
!     
end if
!**********************************************************************C 
!                 DOWN Z
!**********************************************************************C
IF(IDDD.EQ.12) THEN
CALL VMX(1,SQDZ1,SQUZ2,SQDUZ1,1)
DO IC=1, NCOL2
A11(IC)=SQDUZ1(IC) !NEW!
end do
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!
IEEE=4
!     
end if
!**********************************************************************C
!**********************************************************************C
!                 UP WAVE - LIKE PULSE                                *C
!**********************************************************************C
!**********************************************************************C
!                 UP Y
!**********************************************************************C
IF(IDDD.EQ.13) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
do IC=1,NCOL2
D11(IC)=UB11(IC,M2,JU,IDSW)
40 CONTINUE
end do
M3=IUPB(M2,JU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
41 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,JU,IDSW)
42 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,WSQUY1,1)
!                     
M3=IDNB(M2,JU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,JU,IDSW)
43 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
44 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M4,JU,IDSW)
45 CONTINUE
end do
CALL VMX(1,D11,B11,WSQDY2,1)
CALL VMX(1,WSQUY1,WSQDY2,WVUY1,1)
do IC=1,NCOL2
A11(IC)=WVUY1(IC)
46 CONTINUE
end do
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=1
!
end if
!**********************************************************************C
!                 UP Z
!**********************************************************************C
IF(IDDD.EQ.14) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
do IC=1,NCOL2
D11(IC)=UB11(IC,M2,IU,IDSW)
47 CONTINUE
end do
M3=IUPB(M2,IU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
48 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,IU,IDSW)
49 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,WSQUZ1,1)
!                     
M3=IDNB(M2,IU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,IU,IDSW)
50 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
51 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M4,IU,IDSW)
52 CONTINUE
end do
CALL VMX(1,D11,B11,WSQDZ2,1)
CALL VMX(1,WSQUZ1,WSQDZ2,WVUZ1,1)
do IC=1,NCOL2
A11(IC)=WVUZ1(IC)
53 CONTINUE
end do
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=2
!
end if
!**********************************************************************C
!                 DOWN Y                                              *C 
!**********************************************************************C
IF(IDDD.EQ.15) THEN
!     
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
M3=IDNB(M2,JU,IDSW)
do IC=1,NCOL2
D11(IC)=UB11(IC,M3,JU,IDSW)
54 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
55 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,JU,IDSW)
56 CONTINUE
end do
CALL VMX(1,B11,C11,WSQDY1,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
!     
do IC=1,NCOL2
B11(IC)=UB11(IC,M2,JU,IDSW)
57 CONTINUE
end do
M3=IUPB(M2,JU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
58 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,JU,IDSW)
59 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,WSQUY2,1)
CALL VMX(1,WSQDY1,WSQUY2,WVDY1,1)
!     
do IC=1, NCOL2
A11(IC)=WVDY1(IC)
60 CONTINUE
end do
!
IEEE=3
! 
end if
!**********************************************************************C
!                 DOWN Z                                              *C 
!**********************************************************************C
IF(IDDD.EQ.16) THEN
!     
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
M3=IDNB(M2,IU,IDSW)
do IC=1,NCOL2
D11(IC)=UB11(IC,M3,IU,IDSW)
61 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
62 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,IU,IDSW)
63 CONTINUE
end do
CALL VMX(1,B11,C11,WSQDZ1,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
!     
do IC=1,NCOL2
B11(IC)=UB11(IC,M2,IU,IDSW)
64 CONTINUE
end do
M3=IUPB(M2,IU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
65 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,IU,IDSW)
66 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,WSQUZ2,1)
CALL VMX(1,WSQDZ1,WSQUZ2,WVDZ1,1)
!     
do IC=1, NCOL2
A11(IC)=WVDZ1(IC)
67 CONTINUE
end do
!
IEEE=4
!     
end if
!**********************************************************************C
!**********************************************************************C
!              UP - UP WAVE-LIKE PULSE                                 C
!**********************************************************************C
!**********************************************************************C
!              UP Y
!**********************************************************************C
IF(IDDD.EQ.17) THEN
!
IDSW=IDS-1

IF(IDSW.EQ.0) THEN
IDSW=1
end if
! 
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
! 
do IC=1,NCOL2
D11(IC)=UB11(IC,M2,JU,IDSW)
68 CONTINUE
end do
M3=IUPB(M2,JU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
69 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,JU,IDSW)
70 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,WSQUY3,1)
!
M3=IDNB(M2,JU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,JU,IDSW)
71 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
72 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M4,JU,IDSW)
73 CONTINUE
end do
CALL VMX(1,D11,B11,WSQDY4,1)
CALL VMX(1,WSQUY3,WSQDY4,WVUY2,1)
CALL VMX(1,WVUY1,WVUY2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=1
!     
end if
!**********************************************************************C
!                 UP Z
!**********************************************************************C
IF(IDDD.EQ.18) THEN
!
IDSW=IDS-1

IF(IDSW.EQ.0) THEN
IDSW=1
end if
! 
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
! 
do IC=1,NCOL2
D11(IC)=UB11(IC,M2,IU,IDSW)
74 CONTINUE
end do
M3=IUPB(M2,IU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
75 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,IU,IDSW)
76 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,WSQUZ3,1)
!
M3=IDNB(M2,IU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,IU,IDSW)
77 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
78 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M4,IU,IDSW)
79 CONTINUE
end do
CALL VMX(1,D11,B11,WSQDZ4,1)
CALL VMX(1,WSQUZ3,WSQDZ4,WVUZ2,1)
CALL VMX(1,WVUZ1,WVUZ2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=2
!     
end if
!**********************************************************************C
!                 DOWN Y                                               C
!**********************************************************************C
IF(IDDD.EQ.19) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!     
M3=IDNB(M2,JU,IDSW)
do IC=1,NCOL2
D11(IC)=UB11(IC,M3,JU,IDSW)
80 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
81 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,JU,IDSW)
82 CONTINUE
end do
CALL VMX(1,B11,C11,WSQDY3,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
B11(IC)=UB11(IC,M2,JU,IDSW)
83 CONTINUE
end do
M3=IUPB(M2,JU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
84 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,JU,IDSW)
85 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,WSQUY4,1)
CALL VMX(1,WSQDY3,WSQUY4,WVDY2,1)
CALL VMX(1,WVDY1,WVDY2,A11,1)
!
IEEE=3
!     
end if
!**********************************************************************C
!                 DOWN Y                                               C
!**********************************************************************C
IF(IDDD.EQ.20) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!     
M3=IDNB(M2,IU,IDSW)
do IC=1,NCOL2
D11(IC)=UB11(IC,M3,IU,IDSW)
86 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
87 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,IU,IDSW)
88 CONTINUE
end do
CALL VMX(1,B11,C11,WSQDZ3,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
B11(IC)=UB11(IC,M2,IU,IDSW)
89 CONTINUE
end do
M3=IUPB(M2,IU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M3,KU,IDSW)
90 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
do IC=1,NCOL2
C11(IC)=UB11(IC,M2,IU,IDSW)
91 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,WSQUZ4,1)
CALL VMX(1,WSQDZ3,WSQUZ4,WVDZ2,1)
CALL VMX(1,WVDZ1,WVDZ2,A11,1)
!
IEEE=4
!     
end if
!**********************************************************************C
!**********************************************************************C
!                 UP - DOWN WAVE-LIKE PULSE                           *C
!**********************************************************************C
!**********************************************************************C
!                 UP Y                                                 C
!**********************************************************************C
IF(IDDD.EQ.21) THEN
!     
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVUY1,WVDY2,A11,1)
!
IEEE=1
!
end if
!**********************************************************************C
!                 UP Z                                                 C
!**********************************************************************C
IF(IDDD.EQ.22) THEN
!     
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVUZ1,WVDZ2,A11,1)
!
IEEE=2
!
end if
!**********************************************************************C
!                 DOWN Y                                              *C
!**********************************************************************C
IF(IDDD.EQ.23) THEN
!     
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVDY1,WVUY2,A11,1)
!
IEEE=3
!     
end if
!**********************************************************************C
!                 DOWN Z                                              *C
!**********************************************************************C
IF(IDDD.EQ.24) THEN
!     
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!     
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVDZ1,WVUZ2,A11,1)
!
IEEE=4
!     
end if
!**********************************************************************C
!**********************************************************************C
!                   /\_____/\  UP PULSE                                C
!**********************************************************************C
!**********************************************************************C
!                  UP Y
!**********************************************************************C
IF(IDDD.EQ.25) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
139 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
140 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
141 CONTINUE
end do
end if
!
CALL VMX(1,WSQUY1,D11,C11,1)
CALL VMX(1,C11,WSQUY4,A11,1)
!
!
IEEE=1
!
end if
!**********************************************************************C
!                  UP Z
!**********************************************************************C
IF(IDDD.EQ.26) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
142 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
143 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
144 CONTINUE
end do
end if
!
CALL VMX(1,WSQUZ1,D11,C11,1)
CALL VMX(1,C11,WSQUZ4,A11,1)
!
IEEE=2
!
end if
!**********************************************************************C
!                   DOWN Y                                             C
!**********************************************************************C
IF(IDDD.EQ.27) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
145 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
146 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
147 CONTINUE
end do
end if
!     
CALL VMX(1,WSQDY1,D11,C11,1)
CALL VMX(1,C11,WSQDY4,A11,1)
!
IEEE=3
!
end if
!**********************************************************************C
!                   DOWN Z                                             C
!**********************************************************************C
IF(IDDD.EQ.28) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
149 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
150 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
151 CONTINUE
end do
end if
!     
CALL VMX(1,WSQDZ1,D11,C11,1)
CALL VMX(1,C11,WSQDZ4,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
!**********************************************************************C
!                   /\-----\/  PULSE                                   C
!**********************************************************************C
!**********************************************************************C
!                  UP Y
!**********************************************************************C
IF(IDDD.EQ.29) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
152 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
153 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
154 CONTINUE
end do
end if
!
CALL VMX(1,WSQUY1,D11,C11,1)
CALL VMX(1,C11,WSQDY4,A11,1)
!
IEEE=1
!
end if
!**********************************************************************C
!                  UP Z
!**********************************************************************C
IF(IDDD.EQ.30) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
155 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
156 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
157 CONTINUE
end do
end if
!
CALL VMX(1,WSQUZ1,D11,C11,1)
CALL VMX(1,C11,WSQDZ4,A11,1)
!
IEEE=2
!
end if
!**********************************************************************C
!                 DOWN Y                                               C
!**********************************************************************C
IF(IDDD.EQ.31) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
158 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
159 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
160 CONTINUE
end do
end if
!     
CALL VMX(1,WSQDY1,D11,C11,1)
CALL VMX(1,C11,WSQUY4,A11,1)
!
IEEE=3
!     
end if
!**********************************************************************C
!                 DOWN Z                                               C
!**********************************************************************C
IF(IDDD.EQ.32) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
161 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
162 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
262 CONTINUE
end do
end if
!     
CALL VMX(1,WSQDZ1,D11,C11,1)
CALL VMX(1,C11,WSQUZ4,A11,1)
!
IEEE=4
!     
end if
!**********************************************************************
!     TT-1 OPERATORS
!**********************************************************************
IF(IDDD.EQ.33)THEN
!     
CALL VMX(1,SQUY1,SQUZ2,A11,1)
!
DO IC=1, NCOL2
WANGL1(IC)=A11(IC) !NEW!
end do
!                           
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=1
!     
end if
!**********************************************************************
IF(IDDD.EQ.34) THEN
!     
CALL VMX(1,SQUZ1,SQDY2,A11,1)
!
DO IC=1, NCOL2
WANGL2(IC)=A11(IC) !NEW!
end do
!                                                      
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=2
!     
end if
!**********************************************************************
IF(IDDD.EQ.35) THEN
!     
CALL VMX(1,SQDY1,SQDZ2,A11,1)
!    
DO IC=1, NCOL2
WANGL3(IC)=A11(IC) !NEW!
end do
!                           
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=3
!     
end if
!**********************************************************************
IF(IDDD.EQ.36) THEN
!     
CALL VMX(1,SQDZ1,SQUY2,A11,1)
!                           
DO IC=1, NCOL2
WANGL4(IC)=A11(IC) !NEW!
end do
!
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=4
!     
end if
!**********************************************************************
IF(IDDD.EQ.37) THEN
!     
CALL VMX(1,SQUY1,SQDZ2,A11,1)
!                           
DO IC=1, NCOL2
WANGL5(IC)=A11(IC) !NEW!
end do
!
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=5
!     
end if
!**********************************************************************
IF(IDDD.EQ.38) THEN
!     
CALL VMX(1,SQUZ1,SQUY2,A11,1)
!                           
DO IC=1, NCOL2
WANGL6(IC)=A11(IC) !NEW!
end do
!
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=6
!     
end if
!**********************************************************************
IF(IDDD.EQ.39) THEN
!     
CALL VMX(1,SQDY1,SQUZ2,A11,1)
!                           
DO IC=1, NCOL2
WANGL7(IC)=A11(IC) !NEW!
end do
!                           
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
IEEE=7
!     
end if
!**********************************************************************
IF(IDDD.EQ.40) THEN
M3=IUP(M2,KU)
M2=IUP(M3,KU)
!     
CALL VMX(1,SQDZ1,SQDY2,A11,1)
!                           
DO IC=1, NCOL2
WANGL8(IC)=A11(IC) !NEW!
end do
!     
IEEE=8
!     
end if
!**********************************************************************
!     TT-2 OPERATORS
!**********************************************************************
IF(IDDD.EQ.41) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVUY1,WVUZ2,A11,1)
!
IEEE=1
!
end if
!**********************************************************************
IF(IDDD.EQ.42) THEN
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVUZ1,WVDY2,A11,1)
!
IEEE=2
!
end if
!**********************************************************************
IF(IDDD.EQ.43) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVDY1,WVDZ2,A11,1)
!
IEEE=3
!
end if
!**********************************************************************
IF(IDDD.EQ.44) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVDZ1,WVUY2,A11,1)
!
IEEE=4
!
end if
!**********************************************************************
IF(IDDD.EQ.45) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVUZ1,WVUY2,A11,1)
!
IEEE=5
!
end if
!**********************************************************************
IF(IDDD.EQ.46) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVDY1,WVUZ2,A11,1)
!
IEEE=6
!
end if
!**********************************************************************
IF(IDDD.EQ.47) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVDZ1,WVDY2,A11,1)
!
IEEE=7
!
end if
!**********************************************************************
IF(IDDD.EQ.48) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
end if
!
M3=IUPB(M2,KU,IDSW)
M4=IUPB(M3,KU,IDSW)
M1=IUPB(M4,KU,IDSW)
M2=IUPB(M1,KU,IDSW)
!
CALL VMX(1,WVUY1,WVDZ2,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
!**********************************************************************C
!                   /\----/_/  UP PULSE-TT3                            C
!**********************************************************************C
!**********************************************************************C
!                  1
!**********************************************************************C
IF(IDDD.EQ.49) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
301 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
302 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
303 CONTINUE
end do
end if
!                    !   WARNING WE CAN SPEED UP THE COMPUTATION  !
CALL VMX(1,WSQUY1,D11,C11,1)
CALL VMX(1,C11,WSQUZ4,A11,1)
!
!
IEEE=1
!
end if
!**********************************************************************C
!                  2
!**********************************************************************C
IF(IDDD.EQ.50) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
304 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
305 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
306 CONTINUE
end do
end if
!
CALL VMX(1,WSQUZ1,D11,C11,1)
CALL VMX(1,C11,WSQDY4,A11,1)
!
IEEE=2
!
end if
!**********************************************************************C
!                   3                                             C
!**********************************************************************C
IF(IDDD.EQ.51) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
307 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
308 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
309 CONTINUE
end do
end if
!     
CALL VMX(1,WSQDY1,D11,C11,1)
CALL VMX(1,C11,WSQDZ4,A11,1)
!
IEEE=3
!
end if
!**********************************************************************C
!                 4                                                    C
!**********************************************************************C
IF(IDDD.EQ.52) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
310 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
311 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
312 CONTINUE
end do
end if
!     
CALL VMX(1,WSQDZ1,D11,C11,1)
CALL VMX(1,C11,WSQUY4,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
!                 5                                                    C
!**********************************************************************C
IF(IDDD.EQ.53) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
313 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
314 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
315 CONTINUE
end do
end if
!
CALL VMX(1,WSQUY1,D11,C11,1)
CALL VMX(1,C11,WSQDZ4,A11,1)
!
IEEE=5
!
end if
!**********************************************************************C
!                 6                                                    C
!**********************************************************************C
IF(IDDD.EQ.54) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
316 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
317 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
318 CONTINUE
end do
end if
!
CALL VMX(1,WSQUZ1,D11,C11,1)
CALL VMX(1,C11,WSQUY4,A11,1)
!
IEEE=6
!
end if
!**********************************************************************C
!                 7                                                    C
!**********************************************************************C
IF(IDDD.EQ.55) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
319 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
320 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
321 CONTINUE
end do
end if
!
CALL VMX(1,WSQDY1,D11,C11,1)
CALL VMX(1,C11,WSQUZ4,A11,1)
!
IEEE=7
!
end if
!**********************************************************************C
!                 8                                                    C
!**********************************************************************C
IF(IDDD.EQ.56) THEN
!
IDSW=IDS-1
!
IF(IDSW.EQ.0) THEN
IDSW=1
M3=IUPB(M2,KU,IDSW)
do IC=1,NCOL2
B11(IC)=UB11(IC,M3,KU,IDSW)
322 CONTINUE
end do
M4=IUPB(M3,KU,IDSW)
do IC=1,NCOL2
C11(IC)=UB11(IC,M4,KU,IDSW)
323 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
ELSE
M3=IUPB(M2,KU,IDSW)
do IC=1, NCOL2
D11(IC)=UB11(IC,M3,KU,IDSW+1)
324 CONTINUE
end do
end if
!
CALL VMX(1,WSQDZ1,D11,C11,1)
CALL VMX(1,C11,WSQDY4,A11,1)
!
IEEE=8
!
end if
!*********************************************************************
!                 T-T4 OPERATORS
!**********************************************************************
IF(IDDD.EQ.57) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,ZIG1W1,1)
CALL VMX(1,WSQUY3,WSQUZ4,ZIG2W1,1)
CALL VMX(1,ZIG1W1,ZIG2W1,A11,1)
!
IEEE=1
!
end if
!**********************************************************************
IF(IDDD.EQ.58) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,ZIG1W2,1)
CALL VMX(1,WSQUZ3,WSQDY4,ZIG2W2,1)
CALL VMX(1,ZIG1W2,ZIG2W2,A11,1)
!
IEEE=2
!
end if
!**********************************************************************
IF(IDDD.EQ.59) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,ZIG1W3,1)
CALL VMX(1,WSQDY3,WSQDZ4,ZIG2W3,1)
CALL VMX(1,ZIG1W3,ZIG2W3,A11,1)
!
IEEE=3
!
end if
!**********************************************************************
IF(IDDD.EQ.60) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,ZIG1W4,1)
CALL VMX(1,WSQDZ3,WSQUY4,ZIG2W4,1)
CALL VMX(1,ZIG1W4,ZIG2W4,A11,1)
!
IEEE=4
!
end if
!**********************************************************************
IF(IDDD.EQ.61) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,TIG1W4,1)
CALL VMX(1,WSQDZ3,WSQDY4,TIG2W4,1)
CALL VMX(1,TIG1W4,TIG2W4,A11,1)
!
IEEE=5
!
end if
!**********************************************************************
IF(IDDD.EQ.62) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,TIG1W1,1)
CALL VMX(1,WSQUY3,WSQDZ4,TIG2W1,1)
CALL VMX(1,TIG1W1,TIG2W1,A11,1)
!
IEEE=6
!
end if
!**********************************************************************
IF(IDDD.EQ.63) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,TIG1W2,1)
CALL VMX(1,WSQUZ3,WSQUY4,TIG2W2,1)
CALL VMX(1,TIG1W2,TIG2W2,A11,1)
!
IEEE=7
!
end if
!**********************************************************************
IF(IDDD.EQ.64) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,TIG1W3,1)
CALL VMX(1,WSQDY3,WSQUZ4,TIG2W3,1)
CALL VMX(1,TIG1W3,TIG2W3,A11,1)
!
IEEE=8
!
end if
!**********************************************************************
!                 T-T5 OPERATORS
!**********************************************************************
IF(IDDD.EQ.65) THEN
!
CALL VMX(1,WSQUY1,WSQUY2,DUY1,1)
CALL VMX(1,WSQUY3,WSQUY4,DUY2,1)
CALL VMX(1,DUY1,DUY2,A11,1)
!
IEEE=1
!
end if
!**********************************************************************C
IF(IDDD.EQ.66) THEN
!
CALL VMX(1,WSQUZ1,WSQUZ2,DUZ1,1)
CALL VMX(1,WSQUZ3,WSQUZ4,DUZ2,1)
CALL VMX(1,DUZ1,DUZ2,A11,1)
!
IEEE=2
!
end if
!**********************************************************************C
IF(IDDD.EQ.67) THEN
!
CALL VMX(1,WSQDY1,WSQDY2,DDY1,1)
CALL VMX(1,WSQDY3,WSQDY4,DDY2,1)
CALL VMX(1,DDY1,DDY2,A11,1)
!
IEEE=3
!
end if
!**********************************************************************C
IF(IDDD.EQ.68) THEN
!
CALL VMX(1,WSQDZ1,WSQDZ2,DDZ1,1)
CALL VMX(1,WSQDZ3,WSQDZ4,DDZ2,1)
CALL VMX(1,DDZ1,DDZ2,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
!                 TT-6 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.69) THEN
!
CALL VMX(1,DUY1,DDY2,A11,1)
!
IEEE=1
!
end if
!**********************************************************************C
IF(IDDD.EQ.70) THEN
!
CALL VMX(1,DUZ1,DDZ2,A11,1)
!
IEEE=2
!
end if
!**********************************************************************C
IF(IDDD.EQ.71) THEN
!
CALL VMX(1,DDY1,DUY2,A11,1)
!
IEEE=3
!
end if
!**********************************************************************C
IF(IDDD.EQ.72) THEN
!
CALL VMX(1,DDZ1,DUZ2,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
!                 T-T 7 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.73) THEN
!
CALL VMX(1,DUY1,DUZ2,A11,1)
!
IEEE=1
!
end if
!**********************************************************************C
IF(IDDD.EQ.74) THEN
!
CALL VMX(1,DUZ1,DDY2,A11,1)
!
IEEE=2
!
end if
!**********************************************************************C
IF(IDDD.EQ.75) THEN
!
CALL VMX(1,DDY1,DDZ2,A11,1)
!
IEEE=3
!
end if
!**********************************************************************C
IF(IDDD.EQ.76) THEN
!
CALL VMX(1,DDZ1,DUY2,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.77) THEN
!
CALL VMX(1,DDZ1,DDY2,A11,1)
!
IEEE=5
!
end if
!**********************************************************************C
IF(IDDD.EQ.78) THEN
!
CALL VMX(1,DUY1,DDZ2,A11,1)
!
IEEE=6
!
end if
!**********************************************************************C
IF(IDDD.EQ.79) THEN
!
CALL VMX(1,DUZ1,DUY2,A11,1)
!
IEEE=7
!
end if
!**********************************************************************C
IF(IDDD.EQ.80) THEN
!
CALL VMX(1,DDY1,DUZ2,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
!                 T-T 8 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.81) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=1
!
end if
!**********************************************************************c
IF(IDDD.EQ.82) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=2
!
end if
!**********************************************************************c
IF(IDDD.EQ.83) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=3
!
end if
!**********************************************************************c
IF(IDDD.EQ.84) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.85) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=5
!
end if
!**********************************************************************c
IF(IDDD.EQ.86) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=6
!
end if
!**********************************************************************c
IF(IDDD.EQ.87) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=7
!
end if
!**********************************************************************c
IF(IDDD.EQ.88) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
!                 T-T 9 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.89) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=1
!
end if
!**********************************************************************c
IF(IDDD.EQ.90) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=2
!
end if
!**********************************************************************c
IF(IDDD.EQ.91) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=3
!
end if
!**********************************************************************c
IF(IDDD.EQ.92) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.93) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=5
!
end if
!**********************************************************************c
IF(IDDD.EQ.94) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=6
!
end if
!**********************************************************************c
IF(IDDD.EQ.95) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=7
!
end if
!**********************************************************************c
IF(IDDD.EQ.96) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
!                 T-T 10 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.97) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=1
!
end if
!**********************************************************************c
IF(IDDD.EQ.98) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=2
!
end if
!**********************************************************************c
IF(IDDD.EQ.99) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=3
!
end if
!**********************************************************************c
IF(IDDD.EQ.100) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.101) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=5
!
end if
!**********************************************************************c
IF(IDDD.EQ.102) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=6
!
end if
!**********************************************************************c
IF(IDDD.EQ.103) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=7
!
end if
!**********************************************************************c
IF(IDDD.EQ.104) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
!                 T-T 11 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.105) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=1
!
end if
!**********************************************************************c
IF(IDDD.EQ.106) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=2
!
end if
!**********************************************************************c
IF(IDDD.EQ.107) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=3
!
end if
!**********************************************************************c
IF(IDDD.EQ.108) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.109) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=5
!
end if
!**********************************************************************c
IF(IDDD.EQ.110) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=6
!
end if
!**********************************************************************c
IF(IDDD.EQ.111) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=7
!
end if
!**********************************************************************c
IF(IDDD.EQ.112) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
!                 T-T 12 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.113) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=1
!
end if
!**********************************************************************c
IF(IDDD.EQ.114) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=2
!
end if
!**********************************************************************c
IF(IDDD.EQ.115) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=3
!
end if
!**********************************************************************c
IF(IDDD.EQ.116) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.117) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=5
!
end if
!**********************************************************************c
IF(IDDD.EQ.118) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=6
!
end if
!**********************************************************************c
IF(IDDD.EQ.119) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=7
!
end if
!**********************************************************************c
IF(IDDD.EQ.120) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
IF(IDDD.EQ.121) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=9
!
end if
!**********************************************************************c
IF(IDDD.EQ.122) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=10
!
end if
!**********************************************************************c
IF(IDDD.EQ.123) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=11
!
end if
!**********************************************************************c
IF(IDDD.EQ.124) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=12
!
end if
!**********************************************************************C
IF(IDDD.EQ.125) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=13
!
end if
!**********************************************************************c
IF(IDDD.EQ.126) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=14
!
end if
!**********************************************************************c
IF(IDDD.EQ.127) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=15
!
end if
!**********************************************************************c
IF(IDDD.EQ.128) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=16
!
end if
!**********************************************************************C
!                 T-T 13 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.129) THEN
!
CALL VMX(1,WSQUY1,WSQUZ2,B11,1)
CALL VMX(1,WSQDY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=1
!
end if
!**********************************************************************c
IF(IDDD.EQ.130) THEN
!
CALL VMX(1,WSQUZ1,WSQDY2,B11,1)
CALL VMX(1,WSQDZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=2
!
end if
!**********************************************************************c
IF(IDDD.EQ.131) THEN
!
CALL VMX(1,WSQDY1,WSQDZ2,B11,1)
CALL VMX(1,WSQUY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=3
!
end if
!**********************************************************************c
IF(IDDD.EQ.132) THEN
!
CALL VMX(1,WSQDZ1,WSQUY2,B11,1)
CALL VMX(1,WSQUZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.133) THEN
!
CALL VMX(1,WSQUZ1,WSQUY2,B11,1)
CALL VMX(1,WSQDZ3,WSQDY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=5
!
end if
!**********************************************************************c
IF(IDDD.EQ.134) THEN
!
CALL VMX(1,WSQDY1,WSQUZ2,B11,1)
CALL VMX(1,WSQUY3,WSQDZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=6
!
end if
!**********************************************************************c
IF(IDDD.EQ.135) THEN
!
CALL VMX(1,WSQDZ1,WSQDY2,B11,1)
CALL VMX(1,WSQUZ3,WSQUY4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=7
!
end if
!**********************************************************************c
IF(IDDD.EQ.136) THEN
!
CALL VMX(1,WSQUY1,WSQDZ2,B11,1)
CALL VMX(1,WSQDY3,WSQUZ4,C11,1)
CALL VMX(1,B11,C11,A11,1)
!
IEEE=8
!
end if
!**********************************************************************C
!**********************************************************************C
!                 T-T 14 OPERATORS
!**********************************************************************C
IF(IDDD.EQ.137) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQUY1,WSQUZ2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=1
!
end if
!**********************************************************************C
IF(IDDD.EQ.138) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQUZ1,WSQDY2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=2
!
end if
!**********************************************************************C
IF(IDDD.EQ.139) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQDY1,WSQDZ2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=3
!
end if
!**********************************************************************C
IF(IDDD.EQ.140) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQDZ1,WSQUY2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=4
!
end if
!**********************************************************************C
IF(IDDD.EQ.141) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQUY1,WSQDZ2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=5
!
end if
!**********************************************************************C
IF(IDDD.EQ.142) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQUZ1,WSQUY2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=6
!
end if
!**********************************************************************C
IF(IDDD.EQ.143) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQDY1,WSQUZ2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=7
!
end if
!**********************************************************************C
IF(IDDD.EQ.144) THEN
!     
IDSW=IDS-1
!     
IF(IDSW.EQ.0) THEN
IDSW=1
end if
CALL VMX(1,WSQDZ1,WSQDY2,A11,1)
M3=IUPB(M2,KU,IDSW)
M2=M3
M3=IUPB(M2,KU,IDSW)
M2=M3
!
IEEE=8
!
end if
!**********************************************************************C
!                         NORMAL POLYAKOV LOOP                         C
!**********************************************************************C
IF (IDDD.EQ.145) THEN
do IC=1, NCOL2
A11(IC)=UC11(IC,M2,KU)
PL(IC)=UC11(IC,M2,KU)
167 CONTINUE
end do
end if
!**********************************************************************C
!***                      PLAQUETTE OPERATOR                        ***C
!**********************************************************************C
!***                      UP Y OPERATOR
!**********************************************************************C
IF (IDDD.EQ.146) THEN
do IC=1, NCOL2
B11(IC)=UC11(IC,M2,JU)
401 CONTINUE
end do
M3=IUP(M2,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,IU)
402 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUP(M2,IU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M3,JU)
403 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,D11,B11,C11,1)
do IC=1,NCOL2
D11(IC)=UC11(IC,M2,IU)
404 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
CALL VMX(1,C11,D11,PLQ1,1)
!     
CALL VMX(1,PLQ1,PL,PQ1,1)
!
do IC=1,NCOL2
A11(IC)=PQ1(IC)
505 CONTINUE
end do
!
IEEE=1
!     
end if
!**********************************************************************C
!***                      UP Z OPERATOR
!**********************************************************************C
IF(IDDD.EQ.147) THEN
!     
do IC=1, NCOL2
B11(IC)=UC11(IC,M2,IU)
405 CONTINUE
end do
M3=IUP(M2,IU)
M1=IDN(M3,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M1,JU)
406 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,D11,1)
M3=IDN(M2,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,IU)
407 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,B11,1)
do IC=1, NCOL2
D11(IC)=UC11(IC,M3,JU)
408 CONTINUE
end do
CALL VMX(1,B11,D11,PLQ2,1)
!     
CALL VMX(1,PLQ2,PL,PQ2,1)
!     
do IC=1,NCOL2
A11(IC)=PQ2(IC)
509 CONTINUE
end do
!
IEEE=2
!     
end if
!**********************************************************************C
!***                      DOWN Y OPERATOR
!**********************************************************************C
IF(IDDD.EQ.148) THEN
!     
M3=IDN(M2,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,JU)
409 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
M1=IDN(M3,IU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M1,IU)
410 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,C11,B11,D11,1)
do IC=1, NCOL2
C11(IC)=UC11(IC,M1,JU)
411 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M3=IDN(M2,IU)
do IC=1, NCOL2
D11(IC)=UC11(IC,M3,IU)
412 CONTINUE
end do
CALL VMX(1,B11,D11,PLQ3,1)
!     
CALL VMX(1,PLQ3,PL,PQ3,1)
!
do IC=1,NCOL2
A11(IC)=PQ3(IC)
513 CONTINUE
end do
!     
IEEE=3
!     
end if
!**********************************************************************C
!***                      DOWN Z OPERATOR
!**********************************************************************C
IF(IDDD.EQ.149) THEN
!     
M3=IDN(M2,IU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M3,IU)
413 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,JU)
414 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M1=IUP(M3,JU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M1,IU)
415 CONTINUE
end do
CALL VMX(1,D11,B11,C11,1)
do IC=1, NCOL2
B11(IC)=UC11(IC,M2,JU)
416 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,C11,B11,PLQ4,1)
!     
CALL VMX(1,PLQ4,PL,PQ4,1)
!
do IC=1,NCOL2
A11(IC)=PQ4(IC)
517 CONTINUE
end do
!
IEEE=4
!     
end if
!******************************************************************C
!***                       UP Y OPERATOR +
!******************************************************************C
IF(IDDD.EQ.150)THEN
!
do IC=1, NCOL2
PLQ5(IC)=PLQ2(IC)
801 CONTINUE
end do
!
CALL HERM(1,PLQ5,DUM11,1)
CALL VMX(1,PLQ5,PL,PQ5,1)
!
do IC=1,NCOL2
A11(IC)=PQ5(IC)
518 CONTINUE
end do
!     
IEEE=5
!     
end if
!******************************************************************C
!***                       UP Z OPERATOR +
!******************************************************************C
IF(IDDD.EQ.151)THEN
!     
do IC=1, NCOL2
PLQ6(IC)=PLQ1(IC)
802 CONTINUE
end do
!
CALL HERM(1,PLQ6,DUM11,1)
CALL VMX(1,PLQ6,PL,PQ6,1)
!
do IC=1,NCOL2
A11(IC)=PQ6(IC)
519 CONTINUE
end do
!
IEEE=6
!     
end if
!******************************************************************C
!***                       DOWN Y OPERATOR +
!******************************************************************C
IF(IDDD.EQ.152)THEN
!     
do IC=1, NCOL2
PLQ7(IC)=PLQ4(IC)
803 CONTINUE
end do
!
CALL HERM(1,PLQ7,DUM11,1)
CALL VMX(1,PLQ7,PL,PQ7,1)
!     
do IC=1,NCOL2
A11(IC)=PQ7(IC)
520 CONTINUE
end do
!
IEEE=7
!     
end if
!******************************************************************C
!***                       DOWN Z OPERATOR +
!******************************************************************C
IF(IDDD.EQ.153)THEN
!     
do IC=1, NCOL2
PLQ8(IC)=PLQ3(IC)
804 CONTINUE
end do
!
CALL HERM(1,PLQ8,DUM11,1)
CALL VMX(1,PLQ8,PL,PQ8,1)
!
do IC=1,NCOL2
A11(IC)=PQ8(IC)
521 CONTINUE
end do
!
IEEE=8
!
end if
!**********************************************************************C
!***                      PLAQUETTE OPERATOR 2                      ***C
!**********************************************************************C
!***                      UP Y OPERATOR
!**********************************************************************C
IF (IDDD.EQ.154) THEN
M5=IUP(M2,KU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M5,JU)
601 CONTINUE
end do
M3=IUP(M5,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,IU)
602 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M3=IUP(M5,IU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M3,JU)
603 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,D11,B11,C11,1)
do IC=1,NCOL2
D11(IC)=UC11(IC,M5,IU)
604 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
CALL VMX(1,C11,D11,DPLQ1,1)
!     
CALL VMX(1,PQ1,DPLQ1,A11,1)
!
IEEE=1
!     
end if
!**********************************************************************C
!***                      UP Z OPERATOR
!**********************************************************************C
IF(IDDD.EQ.155) THEN
!     
M5=IUP(M2,KU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M5,IU)
605 CONTINUE
end do
M3=IUP(M5,IU)
M1=IDN(M3,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M1,JU)
606 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,B11,C11,D11,1)
M3=IDN(M5,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,IU)
607 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
CALL VMX(1,D11,C11,B11,1)
do IC=1, NCOL2
D11(IC)=UC11(IC,M3,JU)
608 CONTINUE
end do
CALL VMX(1,B11,D11,DPLQ2,1)
!     
CALL VMX(1,PQ2,DPLQ2,A11,1)
!     
IEEE=2
!     
end if
!**********************************************************************C
!***                      DOWN Y OPERATOR
!**********************************************************************C
IF(IDDD.EQ.156) THEN
!     
M5=IUP(M2,KU)
M3=IDN(M5,JU)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,JU)
609 CONTINUE
end do
CALL HERM(1,C11,DUM11,1)
M1=IDN(M3,IU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M1,IU)
610 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,C11,B11,D11,1)
do IC=1, NCOL2
C11(IC)=UC11(IC,M1,JU)
611 CONTINUE
end do
CALL VMX(1,D11,C11,B11,1)
M3=IDN(M5,IU)
do IC=1, NCOL2
D11(IC)=UC11(IC,M3,IU)
612 CONTINUE
end do
CALL VMX(1,B11,D11,DPLQ3,1)
!     
CALL VMX(1,PQ3,DPLQ3,A11,1)
!     
IEEE=3
!     
end if
!**********************************************************************C
!***                      DOWN Z OPERATOR
!**********************************************************************C
IF(IDDD.EQ.157) THEN
!     
M5=IUP(M2,KU)
M3=IDN(M5,IU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M3,IU)
613 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
do IC=1, NCOL2
C11(IC)=UC11(IC,M3,JU)
614 CONTINUE
end do
CALL VMX(1,B11,C11,D11,1)
M1=IUP(M3,JU)
do IC=1, NCOL2
B11(IC)=UC11(IC,M1,IU)
615 CONTINUE
end do
CALL VMX(1,D11,B11,C11,1)
do IC=1, NCOL2
B11(IC)=UC11(IC,M5,JU)
616 CONTINUE
end do
CALL HERM(1,B11,DUM11,1)
CALL VMX(1,C11,B11,DPLQ4,1)
!     
CALL VMX(1,PQ4,DPLQ4,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
!***                       UP Y OPERATOR +
!******************************************************************C
IF(IDDD.EQ.158)THEN
!
do IC=1, NCOL2
DPLQ5(IC)=DPLQ2(IC)
701 CONTINUE
end do
!     
CALL HERM(1,DPLQ5,DUM11,1)
CALL VMX(1,PQ5,DPLQ5,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
!***                       UP Z OPERATOR +
!******************************************************************C
IF(IDDD.EQ.159)THEN
!     
do IC=1, NCOL2
DPLQ6(IC)=DPLQ1(IC)
702 CONTINUE
end do
!
CALL HERM(1,DPLQ6,DUM11,1)
CALL VMX(1,PQ6,DPLQ6,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
!***  DOWN Y OPERATOR +
!******************************************************************C
IF(IDDD.EQ.160)THEN
!     
do IC=1, NCOL2
DPLQ7(IC)=DPLQ4(IC)
703 CONTINUE
end do
!
CALL HERM(1,DPLQ7,DUM11,1)
CALL VMX(1,PQ7,DPLQ7,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
!***                       DOWN Z OPERATOR +
!******************************************************************C
IF(IDDD.EQ.161)THEN
!     
do IC=1, NCOL2
DPLQ8(IC)=DPLQ3(IC)
704 CONTINUE
end do
!
CALL HERM(1,DPLQ8,DUM11,1)
CALL VMX(1,PQ8,DPLQ8,A11,1)
!
IEEE=8
!
end if
!******************************************************************C
!***             PLAQUETTE OPERATPORS 3
!******************************************************************C
IF(IDDD.EQ.162)THEN
!
!                           DO 705 IC=1, NCOL2
!                              B11(IC)=DPLQ1(IC)
! 705                       CONTINUE
!
!                           CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ1,B11,A11,1)
CALL VMX(1,PQ1,DPLQ6,A11,1) ! New altered !
!     
IEEE=1
!
end if
!******************************************************************C
IF(IDDD.EQ.163)THEN
!
!                           DO 706 IC=1, NCOL2
!                              B11(IC)=DPLQ2(IC)
! 706                       CONTINUE
!C
!                          CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ2,B11,A11,1)                           
CALL VMX(1,PQ2,DPLQ5,A11,1) ! New altered !
!     
IEEE=2
!
end if
!******************************************************************C
IF(IDDD.EQ.164)THEN
!
!                           DO 707 IC=1, NCOL2
!                              B11(IC)=DPLQ3(IC)
! 707                       CONTINUE
!C
!                           CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ3,B11,A11,1)

CALL VMX(1,PQ3,DPLQ8,A11,1) ! New altered !
!     
IEEE=3
!
end if
!******************************************************************C
IF(IDDD.EQ.165)THEN
!
!                           DO 708 IC=1, NCOL2
!                              B11(IC)=DPLQ4(IC)
! 708                       CONTINUE
!C
!                          CALL HERM(1,B11,DUM11,1)
!                           CALL VMX(1,PQ4,B11,A11,1)                           
CALL VMX(1,PQ4,DPLQ7,A11,1) ! New altered !
!     
IEEE=4
!
end if
!******************************************************************C
IF(IDDD.EQ.166)THEN
!
CALL VMX(1,PQ5,DPLQ2,A11,1)
!     
IEEE=5
!
end if
!******************************************************************C
IF(IDDD.EQ.167)THEN
!
CALL VMX(1,PQ6,DPLQ1,A11,1)
!     
IEEE=6
!
end if
!******************************************************************C
IF(IDDD.EQ.168)THEN
!
CALL VMX(1,PQ7,DPLQ4,A11,1)
!     
IEEE=7
!
end if
!******************************************************************C
IF(IDDD.EQ.169)THEN
!
CALL VMX(1,PQ8,DPLQ3,A11,1)
!     
IEEE=8
!
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 4 
!******************************************************************C
IF(IDDD.EQ.170)THEN
!
!                           DO 709 IC=1, NCOL2
!                              C11(IC)=DPLQ1(IC)
! 709                       CONTINUE
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ1,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ6,A11,1)
!
IEEE=1
!
end if
!******************************************************************C
IF(IDDD.EQ.171)THEN
!
!                           DO 710 IC=1, NCOL2
!                              C11(IC)=DPLQ2(IC)
! 710                       CONTINUE
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ2,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ5,A11,1)
!
IEEE=2
!
end if
!******************************************************************C
IF(IDDD.EQ.172)THEN
!
!                           DO 711 IC=1, NCOL2
!                              C11(IC)=DPLQ3(IC)
! 711                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ3,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ8,A11,1)
!
IEEE=3
!
end if
!******************************************************************C
IF(IDDD.EQ.173)THEN
!
!                           DO 712 IC=1, NCOL2
!                              C11(IC)=DPLQ4(IC)
! 712                       CONTINUE
!
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ4,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ7,A11,1)
!
IEEE=4
!
end if
!******************************************************************C
IF(IDDD.EQ.174)THEN
!
!                           DO 713 IC=1, NCOL2
!                              C11(IC)=PLQ2(IC)
! 713                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ5,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ2,A11,1)
!
IEEE=5
!
end if
!******************************************************************C
IF(IDDD.EQ.175)THEN
!
!                           DO 714 IC=1, NCOL2
!                              C11(IC)=PLQ1(IC)
! 714                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ6,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ1,A11,1)
!
IEEE=6
!
end if
!******************************************************************C
IF(IDDD.EQ.176)THEN
!
!                           DO 715 IC=1, NCOL2
!                              C11(IC)=PLQ4(IC)
! 715                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ7,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ4,A11,1)
!
IEEE=7
!
end if
!******************************************************************C
IF(IDDD.EQ.177)THEN
!
!                           DO 716 IC=1, NCOL2
!                              C11(IC)=PLQ3(IC)
! 716                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PLQ8,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ3,A11,1)
!
IEEE=8
!
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 5
!******************************************************************C
IF(IDDD.EQ.178)THEN
!
!                           DO 717 IC=1, NCOL2
!                              C11(IC)=DPLQ3(IC)
! 717                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ1,DPLQ8,A11,1)
!
IEEE=1
!
end if
!******************************************************************C
IF(IDDD.EQ.179)THEN
!
!                           DO 718 IC=1, NCOL2
!                              C11(IC)=DPLQ4(IC)
! 718                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ2,DPLQ7,A11,1)
!
IEEE=2
!
end if
!******************************************************************C
IF(IDDD.EQ.180)THEN
!
!                           DO 719 IC=1, NCOL2
!                              C11(IC)=DPLQ1(IC)
! 719                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ3,DPLQ6,A11,1)
!
IEEE=3
!
end if
!******************************************************************C
IF(IDDD.EQ.181)THEN
!
!                           DO 720 IC=1, NCOL2
!                              C11(IC)=DPLQ2(IC)
! 720                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ4,DPLQ5,A11,1)
!
IEEE=4
!
end if
!******************************************************************C
IF(IDDD.EQ.182)THEN
!
CALL VMX(1,PQ5,DPLQ4,A11,1)
!
IEEE=5
!
end if
!******************************************************************C
IF(IDDD.EQ.183)THEN
!
CALL VMX(1,PQ6,DPLQ3,A11,1)
!
IEEE=6
!
end if
!******************************************************************C
IF(IDDD.EQ.184)THEN
!
CALL VMX(1,PQ7,DPLQ2,A11,1)
!
IEEE=7
!
end if
!******************************************************************C
IF(IDDD.EQ.185)THEN
!
CALL VMX(1,PQ8,DPLQ1,A11,1)
!
IEEE=8
!
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 6
!******************************************************************C
IF(IDDD.EQ.186)THEN
!
CALL VMX(1,PQ1,DPLQ3,A11,1)
!     
IEEE=1
!
end if
!******************************************************************C
IF(IDDD.EQ.187)THEN
!
CALL VMX(1,PQ2,DPLQ4,A11,1)
!     
IEEE=2
!
end if
!******************************************************************C
IF(IDDD.EQ.188)THEN
!
CALL VMX(1,PQ3,DPLQ1,A11,1)
!     
IEEE=3
!
end if
!******************************************************************C
IF(IDDD.EQ.189)THEN
!
CALL VMX(1,PQ4,DPLQ2,A11,1)
!     
IEEE=4
!
end if
!******************************************************************C
IF(IDDD.EQ.190)THEN
!
!                           DO 721 IC=1, NCOL2
!                              C11(IC)=DPLQ4(IC)
! 721                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ5,DPLQ7,A11,1)
!     
IEEE=5
!
end if
!******************************************************************C
IF(IDDD.EQ.191)THEN
!
!                           DO 722 IC=1, NCOL2
!                              C11(IC)=DPLQ3(IC)
! 722                       CONTINUE
!
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ6,DPLQ8,A11,1)
!     
IEEE=6
!
end if
!******************************************************************C
IF(IDDD.EQ.192)THEN
!
!                           DO 723 IC=1, NCOL2
!                              C11(IC)=DPLQ2(IC)
! 723                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ7,DPLQ5,A11,1)
!     
IEEE=7
!
end if
!******************************************************************C
IF(IDDD.EQ.193)THEN
!
!                           DO 724 IC=1, NCOL2
!                              C11(IC)=DPLQ1(IC)
! 724                       CONTINUE
!C
!                           CALL HERM(1,C11,DUM11,1)
CALL VMX(1,PQ8,DPLQ6,A11,1)
!     
IEEE=8
!
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 7
!******************************************************************C
IF(IDDD.EQ.194)THEN
!
CALL VMX(1,SQUY1,DPLQ1,B11,1)
CALL VMX(1,B11,SQDZ2,A11,1)
!     
IEEE=1
!
end if
!******************************************************************C                        
IF(IDDD.EQ.195)THEN
!
CALL VMX(1,SQUZ1,DPLQ2,B11,1)
CALL VMX(1,B11,SQUY2,A11,1)
!     
IEEE=2
!
end if
!******************************************************************C                        
IF(IDDD.EQ.196)THEN
!
CALL VMX(1,SQDY1,DPLQ3,B11,1)
CALL VMX(1,B11,SQUZ2,A11,1)
!     
IEEE=3
!
end if
!******************************************************************C                        
IF(IDDD.EQ.197)THEN
!
CALL VMX(1,SQDZ1,DPLQ4,B11,1)
CALL VMX(1,B11,SQDY2,A11,1)
!     
IEEE=4
!
end if
!******************************************************************C                        
IF(IDDD.EQ.198)THEN
!
CALL VMX(1,SQDZ1,DPLQ6,B11,1)
CALL VMX(1,B11,SQUY2,A11,1)
!     
IEEE=5
!
end if
!******************************************************************C                        
IF(IDDD.EQ.199)THEN
!
CALL VMX(1,SQUY1,DPLQ5,B11,1)
CALL VMX(1,B11,SQUZ2,A11,1)
!     
IEEE=6
!
end if
!******************************************************************C                        
IF(IDDD.EQ.200)THEN
!
CALL VMX(1,SQUZ1,DPLQ8,B11,1)
CALL VMX(1,B11,SQDY2,A11,1)
!     
IEEE=7
!
end if
!******************************************************************C                        
IF(IDDD.EQ.201)THEN
!
CALL VMX(1,SQDY1,DPLQ7,B11,1)
CALL VMX(1,B11,SQDZ2,A11,1)
!     
IEEE=8
!
end if
!******************************************************************C                        
IF(IDDD.EQ.202)THEN
!
CALL VMX(1,SQDY1,DPLQ5,B11,1)
CALL VMX(1,B11,SQDZ2,A11,1)
!     
IEEE=9
!
end if
!******************************************************************C                        
IF(IDDD.EQ.203)THEN
!
CALL VMX(1,SQUZ1,DPLQ6,B11,1)
CALL VMX(1,B11,SQDY2,A11,1)
!     
IEEE=10
!
end if
!******************************************************************C                        
IF(IDDD.EQ.204)THEN
!
CALL VMX(1,SQUY1,DPLQ7,B11,1)
CALL VMX(1,B11,SQUZ2,A11,1)
!     
IEEE=11
!
end if
!******************************************************************C                        
IF(IDDD.EQ.205)THEN
!
CALL VMX(1,SQDZ1,DPLQ8,B11,1)
CALL VMX(1,B11,SQUY2,A11,1)
!     
IEEE=12
!
end if
!******************************************************************C                        
IF(IDDD.EQ.206)THEN
!
CALL VMX(1,SQDZ1,DPLQ2,B11,1)
CALL VMX(1,B11,SQDY2,A11,1)
!     
IEEE=13
!
end if
!******************************************************************C                        
IF(IDDD.EQ.207)THEN
!
CALL VMX(1,SQDY1,DPLQ1,B11,1)
CALL VMX(1,B11,SQUZ2,A11,1)
!     
IEEE=14
!
end if
!******************************************************************C                        
IF(IDDD.EQ.208)THEN
!
CALL VMX(1,SQUZ1,DPLQ4,B11,1)
CALL VMX(1,B11,SQUY2,A11,1)
!     
IEEE=15
!
end if
!******************************************************************C                        
IF(IDDD.EQ.209)THEN
!
CALL VMX(1,SQUY1,DPLQ3,B11,1)
CALL VMX(1,B11,SQDZ2,A11,1)
!     
IEEE=16
!
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 8
!******************************************************************C
IF(IDDD.EQ.210)THEN
!
CALL VMX(1,SQUY1,DPLQ1,B11,1)
CALL VMX(1,B11,DPLQ2,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************c
IF(IDDD.EQ.211)THEN
!
CALL VMX(1,SQUZ1,DPLQ2,B11,1)
CALL VMX(1,B11,DPLQ3,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************c
IF(IDDD.EQ.212)THEN
!
CALL VMX(1,SQDY1,DPLQ3,B11,1)
CALL VMX(1,B11,DPLQ4,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************c
IF(IDDD.EQ.213)THEN
!
CALL VMX(1,SQDZ1,DPLQ4,B11,1)
CALL VMX(1,B11,DPLQ1,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************c
IF(IDDD.EQ.214)THEN
!
CALL VMX(1,PLQ5,PLQ6,B11,1)
CALL VMX(1,B11,SQUY1,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************c
IF(IDDD.EQ.215)THEN
!
CALL VMX(1,PLQ8,PLQ5,B11,1)
CALL VMX(1,B11,SQUZ1,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************c
IF(IDDD.EQ.216)THEN
!
CALL VMX(1,PLQ7,PLQ8,B11,1)
CALL VMX(1,B11,SQDY1,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************c
IF(IDDD.EQ.217)THEN
!
CALL VMX(1,PLQ6,PLQ7,B11,1)
CALL VMX(1,B11,SQDZ1,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************c
IF(IDDD.EQ.218)THEN
!
CALL VMX(1,SQDY1,DPLQ5,B11,1)
CALL VMX(1,B11,DPLQ6,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************c
IF(IDDD.EQ.219)THEN
!
CALL VMX(1,SQUZ1,DPLQ6,B11,1)
CALL VMX(1,B11,DPLQ7,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************c
IF(IDDD.EQ.220)THEN
!
CALL VMX(1,SQUY1,DPLQ7,B11,1)
CALL VMX(1,B11,DPLQ8,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************c
IF(IDDD.EQ.221)THEN
!
CALL VMX(1,SQDZ1,DPLQ8,B11,1)
CALL VMX(1,B11,DPLQ5,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************c
IF(IDDD.EQ.222)THEN
!
CALL VMX(1,PLQ1,PLQ2,B11,1)
CALL VMX(1,B11,SQDY1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************c
IF(IDDD.EQ.223)THEN
!
CALL VMX(1,PLQ4,PLQ1,B11,1)
CALL VMX(1,B11,SQUZ1,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************c
IF(IDDD.EQ.224)THEN
!
CALL VMX(1,PLQ3,PLQ4,B11,1)
CALL VMX(1,B11,SQUY1,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************c
IF(IDDD.EQ.225)THEN
!
CALL VMX(1,PLQ2,PLQ3,B11,1)
CALL VMX(1,B11,SQDZ1,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 9
!******************************************************************C
IF(IDDD.EQ.226)THEN
!     
CALL VMX(1,SQUY1,DPLQ1,B11,1)
CALL VMX(1,B11,DPLQ2,C11,1)
CALL VMX(1,C11,DPLQ3,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************C
IF(IDDD.EQ.227)THEN
!     
CALL VMX(1,SQUZ1,DPLQ2,B11,1)
CALL VMX(1,B11,DPLQ3,C11,1)
CALL VMX(1,C11,DPLQ4,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************C
IF(IDDD.EQ.228)THEN
!     
CALL VMX(1,SQDY1,DPLQ3,B11,1)
CALL VMX(1,B11,DPLQ4,C11,1)
CALL VMX(1,C11,DPLQ1,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************C
IF(IDDD.EQ.229)THEN
!     
CALL VMX(1,SQDZ1,DPLQ4,B11,1)
CALL VMX(1,B11,DPLQ1,C11,1)
CALL VMX(1,C11,DPLQ2,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
IF(IDDD.EQ.230)THEN
!     
CALL VMX(1,PLQ8,PLQ5,B11,1)
CALL VMX(1,B11,PLQ6,C11,1)
CALL VMX(1,C11,SQUY1,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
IF(IDDD.EQ.231)THEN
!     
CALL VMX(1,PLQ7,PLQ8,B11,1)
CALL VMX(1,B11,PLQ5,C11,1)
CALL VMX(1,C11,SQUZ1,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
IF(IDDD.EQ.232)THEN
!     
CALL VMX(1,PLQ6,PLQ7,B11,1)
CALL VMX(1,B11,PLQ8,C11,1)
CALL VMX(1,C11,SQDY1,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
IF(IDDD.EQ.233)THEN
!     
CALL VMX(1,PLQ5,PLQ6,B11,1)
CALL VMX(1,B11,PLQ7,C11,1)
CALL VMX(1,C11,SQDZ1,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************C
IF(IDDD.EQ.234)THEN
!
CALL VMX(1,SQDY1,DPLQ5,B11,1)
CALL VMX(1,B11,DPLQ6,C11,1)
CALL VMX(1,C11,DPLQ7,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************C
IF(IDDD.EQ.235)THEN
!
CALL VMX(1,SQUZ1,DPLQ6,B11,1)
CALL VMX(1,B11,DPLQ7,C11,1)
CALL VMX(1,C11,DPLQ8,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************C
IF(IDDD.EQ.236)THEN
!
CALL VMX(1,SQUY1,DPLQ7,B11,1)
CALL VMX(1,B11,DPLQ8,C11,1)
CALL VMX(1,C11,DPLQ5,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************C
IF(IDDD.EQ.237)THEN
!
CALL VMX(1,SQDZ1,DPLQ8,B11,1)
CALL VMX(1,B11,DPLQ5,C11,1)
CALL VMX(1,C11,DPLQ6,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************C
IF(IDDD.EQ.238)THEN
!
CALL VMX(1,PLQ4,PLQ1,B11,1)
CALL VMX(1,B11,PLQ2,C11,1)
CALL VMX(1,C11,SQDY1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************C
IF(IDDD.EQ.239)THEN
!
CALL VMX(1,PLQ3,PLQ4,B11,1)
CALL VMX(1,B11,PLQ1,C11,1)
CALL VMX(1,C11,SQUZ1,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************C
IF(IDDD.EQ.240)THEN
!
CALL VMX(1,PLQ2,PLQ3,B11,1)
CALL VMX(1,B11,PLQ4,C11,1)
CALL VMX(1,C11,SQUY1,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************C
IF(IDDD.EQ.241)THEN
!
CALL VMX(1,PLQ1,PLQ2,B11,1)
CALL VMX(1,B11,PLQ3,C11,1)
CALL VMX(1,C11,SQDZ1,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 10
!******************************************************************C
IF(IDDD.EQ.242)THEN
!     
CALL VMX(1,PLQ2,PLQ7,B11,1)
CALL VMX(1,B11,SQDZ1,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************C
IF(IDDD.EQ.243)THEN
!     
CALL VMX(1,PLQ3,PLQ6,B11,1)
CALL VMX(1,B11,SQUY1,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************C
IF(IDDD.EQ.244)THEN
!     
CALL VMX(1,PLQ4,PLQ5,B11,1)
CALL VMX(1,B11,SQUZ1,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************C
IF(IDDD.EQ.245)THEN
!     
CALL VMX(1,PLQ1,PLQ8,B11,1)
CALL VMX(1,B11,SQDY1,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
IF(IDDD.EQ.246)THEN
!     
CALL VMX(1,SQDZ1,DPLQ4,B11,1)
CALL VMX(1,B11,DPLQ5,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
IF(IDDD.EQ.247)THEN
!     
CALL VMX(1,SQUY1,DPLQ1,B11,1)
CALL VMX(1,B11,DPLQ8,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
IF(IDDD.EQ.248)THEN
!     
CALL VMX(1,SQUZ1,DPLQ2,B11,1)
CALL VMX(1,B11,DPLQ7,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
IF(IDDD.EQ.249)THEN
!     
CALL VMX(1,SQDY1,DPLQ3,B11,1)
CALL VMX(1,B11,DPLQ6,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************C
IF(IDDD.EQ.250)THEN
!     
CALL VMX(1,PLQ6,PLQ3,B11,1)
CALL VMX(1,B11,SQDZ1,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************C
IF(IDDD.EQ.251)THEN
!     
CALL VMX(1,PLQ7,PLQ2,B11,1)
CALL VMX(1,B11,SQDY1,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************C
IF(IDDD.EQ.252)THEN
!     
CALL VMX(1,PLQ8,PLQ1,B11,1)
CALL VMX(1,B11,SQUZ1,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************C
IF(IDDD.EQ.253)THEN
!     
CALL VMX(1,PLQ5,PLQ4,B11,1)
CALL VMX(1,B11,SQUY1,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************C
IF(IDDD.EQ.254)THEN
!     
CALL VMX(1,SQDZ1,DPLQ8,B11,1)
CALL VMX(1,B11,DPLQ1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************C
IF(IDDD.EQ.255)THEN
!     
CALL VMX(1,SQDY1,DPLQ5,B11,1)
CALL VMX(1,B11,DPLQ4,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************C
IF(IDDD.EQ.256)THEN
!     
CALL VMX(1,SQUZ1,DPLQ6,B11,1)
CALL VMX(1,B11,DPLQ3,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************C
IF(IDDD.EQ.257)THEN
!     
CALL VMX(1,SQUY1,DPLQ7,B11,1)
CALL VMX(1,B11,DPLQ2,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 11
!******************************************************************C
IF(IDDD.EQ.258)THEN
!     
CALL VMX(1,PLQ2,PLQ4,B11,1)
CALL VMX(1,B11,SQUY1,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************C
IF(IDDD.EQ.259)THEN
!     
CALL VMX(1,PLQ3,PLQ1,B11,1)
CALL VMX(1,B11,SQUZ1,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************C
IF(IDDD.EQ.260)THEN
!     
CALL VMX(1,PLQ4,PLQ2,B11,1)
CALL VMX(1,B11,SQDY1,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************C
IF(IDDD.EQ.261)THEN
!     
CALL VMX(1,PLQ1,PLQ3,B11,1)
CALL VMX(1,B11,SQDZ1,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
IF(IDDD.EQ.262)THEN
!
CALL VMX(1,SQUY1,DPLQ7,B11,1)
CALL VMX(1,B11,DPLQ5,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
IF(IDDD.EQ.263)THEN
!
CALL VMX(1,SQUZ1,DPLQ6,B11,1)
CALL VMX(1,B11,DPLQ8,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
IF(IDDD.EQ.264)THEN
!
CALL VMX(1,SQDY1,DPLQ5,B11,1)
CALL VMX(1,B11,DPLQ7,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
IF(IDDD.EQ.265)THEN
!
CALL VMX(1,SQDZ1,DPLQ8,B11,1)
CALL VMX(1,B11,DPLQ6,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************C
IF(IDDD.EQ.266)THEN
!     
CALL VMX(1,PLQ6,PLQ8,B11,1)
CALL VMX(1,B11,SQDY1,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************C
IF(IDDD.EQ.267)THEN
!     
CALL VMX(1,PLQ7,PLQ5,B11,1)
CALL VMX(1,B11,SQUZ1,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************C
IF(IDDD.EQ.268)THEN
!     
CALL VMX(1,PLQ8,PLQ6,B11,1)
CALL VMX(1,B11,SQUY1,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************C
IF(IDDD.EQ.269)THEN
!     
CALL VMX(1,PLQ5,PLQ7,B11,1)
CALL VMX(1,B11,SQDZ1,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************C
IF(IDDD.EQ.270)THEN
!
CALL VMX(1,SQDY1,DPLQ3,B11,1)
CALL VMX(1,B11,DPLQ1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************C
IF(IDDD.EQ.271)THEN
!
CALL VMX(1,SQUZ1,DPLQ2,B11,1)
CALL VMX(1,B11,DPLQ4,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************C
IF(IDDD.EQ.272)THEN
!
CALL VMX(1,SQUY1,DPLQ1,B11,1)
CALL VMX(1,B11,DPLQ3,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************C
IF(IDDD.EQ.273)THEN
!
CALL VMX(1,SQDZ1,DPLQ4,B11,1)
CALL VMX(1,B11,DPLQ2,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 12
!******************************************************************C
IF(IDDD.EQ.274)THEN
!     
CALL VMX(1,PLQ2,PLQ3,B11,1)
CALL VMX(1,B11,SQUY1,C11,1)
CALL VMX(1,C11,DPLQ7,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************C
IF(IDDD.EQ.275)THEN
!     
CALL VMX(1,PLQ3,PLQ4,B11,1)
CALL VMX(1,B11,SQUZ1,C11,1)
CALL VMX(1,C11,DPLQ6,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************C
IF(IDDD.EQ.276)THEN
!     
CALL VMX(1,PLQ4,PLQ1,B11,1)
CALL VMX(1,B11,SQDY1,C11,1)
CALL VMX(1,C11,DPLQ5,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************C
IF(IDDD.EQ.277)THEN
!     
CALL VMX(1,PLQ1,PLQ2,B11,1)
CALL VMX(1,B11,SQDZ1,C11,1)
CALL VMX(1,C11,DPLQ8,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
IF(IDDD.EQ.278)THEN
!     
CALL VMX(1,PLQ4,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ8,C11,1)
CALL VMX(1,C11,DPLQ5,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
IF(IDDD.EQ.279)THEN
!     
CALL VMX(1,PLQ1,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ7,C11,1)
CALL VMX(1,C11,DPLQ8,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
IF(IDDD.EQ.280)THEN
!     
CALL VMX(1,PLQ2,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ6,C11,1)
CALL VMX(1,C11,DPLQ7,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
IF(IDDD.EQ.281)THEN
!     
CALL VMX(1,PLQ3,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ5,C11,1)
CALL VMX(1,C11,DPLQ6,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************C
IF(IDDD.EQ.282)THEN
!     
CALL VMX(1,PLQ6,PLQ7,B11,1)
CALL VMX(1,B11,SQDY1,C11,1)
CALL VMX(1,C11,DPLQ3,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************C
IF(IDDD.EQ.283)THEN
!     
CALL VMX(1,PLQ7,PLQ8,B11,1)
CALL VMX(1,B11,SQUZ1,C11,1)
CALL VMX(1,C11,DPLQ2,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************C
IF(IDDD.EQ.284)THEN
!     
CALL VMX(1,PLQ8,PLQ5,B11,1)
CALL VMX(1,B11,SQUY1,C11,1)
CALL VMX(1,C11,DPLQ1,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************C
IF(IDDD.EQ.285)THEN
!     
CALL VMX(1,PLQ5,PLQ6,B11,1)
CALL VMX(1,B11,SQDZ1,C11,1)
CALL VMX(1,C11,DPLQ4,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************C
IF(IDDD.EQ.286)THEN
!     
CALL VMX(1,PLQ8,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ4,C11,1)
CALL VMX(1,C11,DPLQ1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************C
IF(IDDD.EQ.287)THEN
!     
CALL VMX(1,PLQ5,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ3,C11,1)
CALL VMX(1,C11,DPLQ4,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************C
IF(IDDD.EQ.288)THEN
!     
CALL VMX(1,PLQ6,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ2,C11,1)
CALL VMX(1,C11,DPLQ3,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************C
IF(IDDD.EQ.289)THEN
!     
CALL VMX(1,PLQ7,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ1,C11,1)
CALL VMX(1,C11,DPLQ2,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 13
!******************************************************************C
IF(IDDD.EQ.290)THEN
!     
CALL VMX(1,PLQ2,PLQ3,B11,1)
CALL VMX(1,B11,PLQ4,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ1,C11,1)
CALL VMX(1,C11,DPLQ2,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************C
IF(IDDD.EQ.291)THEN
!     
CALL VMX(1,PLQ3,PLQ4,B11,1)
CALL VMX(1,B11,PLQ1,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ2,C11,1)
CALL VMX(1,C11,DPLQ3,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************C
IF(IDDD.EQ.292)THEN
!     
CALL VMX(1,PLQ4,PLQ1,B11,1)
CALL VMX(1,B11,PLQ2,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ3,C11,1)
CALL VMX(1,C11,DPLQ4,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************C
IF(IDDD.EQ.293)THEN
!     
CALL VMX(1,PLQ1,PLQ2,B11,1)
CALL VMX(1,B11,PLQ3,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ4,C11,1)
CALL VMX(1,C11,DPLQ1,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
IF(IDDD.EQ.294)THEN
!     
CALL VMX(1,PLQ5,PLQ6,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ7,B11,1)
CALL VMX(1,B11,DPLQ8,C11,1)
CALL VMX(1,C11,DPLQ5,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
IF(IDDD.EQ.295)THEN
!     
CALL VMX(1,PLQ8,PLQ5,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ6,B11,1)
CALL VMX(1,B11,DPLQ7,C11,1)
CALL VMX(1,C11,DPLQ8,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
IF(IDDD.EQ.296)THEN
!     
CALL VMX(1,PLQ7,PLQ8,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ5,B11,1)
CALL VMX(1,B11,DPLQ6,C11,1)
CALL VMX(1,C11,DPLQ7,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
IF(IDDD.EQ.297)THEN
!     
CALL VMX(1,PLQ6,PLQ7,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ8,B11,1)
CALL VMX(1,B11,DPLQ5,C11,1)
CALL VMX(1,C11,DPLQ6,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************C
IF(IDDD.EQ.298)THEN
!     
CALL VMX(1,PLQ6,PLQ7,B11,1)
CALL VMX(1,B11,PLQ8,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ5,C11,1)
CALL VMX(1,C11,DPLQ6,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************C
IF(IDDD.EQ.299)THEN
!     
CALL VMX(1,PLQ7,PLQ8,B11,1)
CALL VMX(1,B11,PLQ5,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ6,C11,1)
CALL VMX(1,C11,DPLQ7,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************C
IF(IDDD.EQ.300)THEN
!     
CALL VMX(1,PLQ8,PLQ5,B11,1)
CALL VMX(1,B11,PLQ6,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ7,C11,1)
CALL VMX(1,C11,DPLQ8,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************C
IF(IDDD.EQ.301)THEN
!     
CALL VMX(1,PLQ5,PLQ6,B11,1)
CALL VMX(1,B11,PLQ7,C11,1)
CALL VMX(1,C11,PL,B11,1)
CALL VMX(1,B11,DPLQ8,C11,1)
CALL VMX(1,C11,DPLQ5,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************C
IF(IDDD.EQ.302)THEN
!     
CALL VMX(1,PLQ1,PLQ2,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ3,B11,1)
CALL VMX(1,B11,DPLQ4,C11,1)
CALL VMX(1,C11,DPLQ1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************C
IF(IDDD.EQ.303)THEN
!     
CALL VMX(1,PLQ4,PLQ1,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ2,B11,1)
CALL VMX(1,B11,DPLQ3,C11,1)
CALL VMX(1,C11,DPLQ4,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************C
IF(IDDD.EQ.304)THEN
!     
CALL VMX(1,PLQ3,PLQ4,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ1,B11,1)
CALL VMX(1,B11,DPLQ2,C11,1)
CALL VMX(1,C11,DPLQ3,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************C
IF(IDDD.EQ.305)THEN
!     
CALL VMX(1,PLQ2,PLQ3,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ4,B11,1)
CALL VMX(1,B11,DPLQ1,C11,1)
CALL VMX(1,C11,DPLQ2,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 14
!******************************************************************C
IF(IDDD.EQ.306)THEN
!     
CALL VMX(1,PLQ2,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ7,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************C
IF(IDDD.EQ.307)THEN
!     
CALL VMX(1,PLQ3,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ6,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************C
IF(IDDD.EQ.308)THEN
!     
CALL VMX(1,PLQ4,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ5,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************C
IF(IDDD.EQ.309)THEN
!     
CALL VMX(1,PLQ1,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ8,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
IF(IDDD.EQ.310)THEN
!     
CALL VMX(1,PLQ4,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ5,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
IF(IDDD.EQ.311)THEN
!     
CALL VMX(1,PLQ1,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ8,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
IF(IDDD.EQ.312)THEN
!     
CALL VMX(1,PLQ2,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ7,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
IF(IDDD.EQ.313)THEN
!     
CALL VMX(1,PLQ3,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ6,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************C
IF(IDDD.EQ.314)THEN
!     
CALL VMX(1,PLQ6,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ3,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************C
IF(IDDD.EQ.315)THEN
!     
CALL VMX(1,PLQ7,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ2,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************C
IF(IDDD.EQ.316)THEN
!     
CALL VMX(1,PLQ8,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ1,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************C
IF(IDDD.EQ.317)THEN
!     
CALL VMX(1,PLQ5,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ4,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************C
IF(IDDD.EQ.318)THEN
!     
CALL VMX(1,PLQ8,SQDY1,B11,1)
CALL VMX(1,B11,DPLQ1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************C
IF(IDDD.EQ.319)THEN
!     
CALL VMX(1,PLQ5,SQUZ1,B11,1)
CALL VMX(1,B11,DPLQ4,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************C
IF(IDDD.EQ.320)THEN
!     
CALL VMX(1,PLQ6,SQUY1,B11,1)
CALL VMX(1,B11,DPLQ3,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************C
IF(IDDD.EQ.321)THEN
!     
CALL VMX(1,PLQ7,SQDZ1,B11,1)
CALL VMX(1,B11,DPLQ2,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************C
!     PLAQUETTE OPERATORS 15
!******************************************************************C
IF(IDDD.EQ.322)THEN
!     
CALL VMX(1,PLQ2,PLQ7,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ6,A11,1)
!     
IEEE=1
!     
end if
!******************************************************************C
IF(IDDD.EQ.323)THEN
!     
CALL VMX(1,PLQ3,PLQ6,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ5,A11,1)
!     
IEEE=2
!     
end if
!******************************************************************C
IF(IDDD.EQ.324)THEN
!     
CALL VMX(1,PLQ4,PLQ5,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ8,A11,1)
!     
IEEE=3
!     
end if
!******************************************************************C
IF(IDDD.EQ.325)THEN
!     
CALL VMX(1,PLQ1,PLQ8,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ7,A11,1)
!     
IEEE=4
!     
end if
!******************************************************************C
IF(IDDD.EQ.326)THEN
!     
CALL VMX(1,PLQ1,PL,B11,1)
CALL VMX(1,B11,DPLQ4,C11,1)
CALL VMX(1,C11,DPLQ5,A11,1)
!     
IEEE=5
!     
end if
!******************************************************************C
IF(IDDD.EQ.327)THEN
!     
CALL VMX(1,PLQ2,PL,B11,1)
CALL VMX(1,B11,DPLQ1,C11,1)
CALL VMX(1,C11,DPLQ8,A11,1)
!     
IEEE=6
!     
end if
!******************************************************************C
IF(IDDD.EQ.328)THEN
!     
CALL VMX(1,PLQ3,PL,B11,1)
CALL VMX(1,B11,DPLQ2,C11,1)
CALL VMX(1,C11,DPLQ7,A11,1)
!     
IEEE=7
!     
end if
!******************************************************************C
IF(IDDD.EQ.329)THEN
!     
CALL VMX(1,PLQ4,PL,B11,1)
CALL VMX(1,B11,DPLQ3,C11,1)
CALL VMX(1,C11,DPLQ6,A11,1)
!     
IEEE=8
!     
end if
!******************************************************************C
IF(IDDD.EQ.330)THEN
!     
CALL VMX(1,PLQ6,PLQ3,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ2,A11,1)
!     
IEEE=9
!     
end if
!******************************************************************C
IF(IDDD.EQ.331)THEN
!     
CALL VMX(1,PLQ7,PLQ2,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ1,A11,1)
!     
IEEE=10
!     
end if
!******************************************************************C
IF(IDDD.EQ.332)THEN
!     
CALL VMX(1,PLQ8,PLQ1,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ4,A11,1)
!     
IEEE=11
!     
end if
!******************************************************************C
IF(IDDD.EQ.333)THEN
!     
CALL VMX(1,PLQ5,PLQ4,B11,1)
CALL VMX(1,B11,PL,C11,1)
CALL VMX(1,C11,DPLQ3,A11,1)
!     
IEEE=12
!     
end if
!******************************************************************C
IF(IDDD.EQ.334)THEN
!     
CALL VMX(1,PLQ5,PL,B11,1)
CALL VMX(1,B11,DPLQ8,C11,1)
CALL VMX(1,C11,DPLQ1,A11,1)
!     
IEEE=13
!     
end if
!******************************************************************C
IF(IDDD.EQ.335)THEN
!     
CALL VMX(1,PLQ6,PL,B11,1)
CALL VMX(1,B11,DPLQ5,C11,1)
CALL VMX(1,C11,DPLQ4,A11,1)
!     
IEEE=14
!     
end if
!******************************************************************C
IF(IDDD.EQ.336)THEN
!     
CALL VMX(1,PLQ7,PL,B11,1)
CALL VMX(1,B11,DPLQ6,C11,1)
CALL VMX(1,C11,DPLQ3,A11,1)
!     
IEEE=15
!     
end if
!******************************************************************C
IF(IDDD.EQ.337)THEN
!     
CALL VMX(1,PLQ8,PL,B11,1)
CALL VMX(1,B11,DPLQ7,C11,1)
CALL VMX(1,C11,DPLQ2,A11,1)
!     
IEEE=16
!     
end if
!******************************************************************c
end if
!******************************************************************C
IF(LCNT(IDS).LT.ICO) THEN
do IC=1, NCOL2
A11(IC)=LIN0(IC)
168 CONTINUE
end do
M2=ML
ELSE
IF(ICO.EQ.1) THEN
CALL VMX(1,A11,LIN1,C11,1)
do IC=1,NCOL2
A11(IC)=C11(IC)
169 CONTINUE
end do
end if
!
IF(ICO.EQ.2) THEN
CALL VMX(1,A11,LIN2,C11,1)
do IC=1,NCOL2
A11(IC)=C11(IC)
170 CONTINUE
end do
end if
!
IF(ICO.EQ.4) THEN
CALL VMX(1,A11,LIN4,C11,1)
do IC=1,NCOL2
A11(IC)=C11(IC)
171 CONTINUE
end do
end if
M2=ML
end if
!***********************************************************************
end if
CALL VMX(1,A11,REM11,B11,1)
!***********************************************************************
!            CALL RENORMBS(B11,UREN11)
!            DO IJ=1,NCOL2
!               AA(IJ)=UREN11(IJ)
!            ENDDO
!            JMAT=NCOL
!            CALL DETNANT(JMAT,DET,AA)
!c     
!            DNCOL=CMPLX(1.0/NCOL)
!            CDET=DET**DNCOL
!            CNORM=1.0/CDET
!c     
!            DO IC=1,NCOL2
!               B11(IC)=UREN11(IC)*CNORM
!            ENDDO
!***********************************************************************
AKT1=(0.0,0.0)
!***********************************************************************
do N1=1,NCOL
IJ=N1+NCOL*(N1-1)
AKT1=AKT1+B11(IJ)
192 CONTINUE
end do
!************ ***********************************************************     
IF(IDDD.LE.4) THEN
CSUMS(IEEE)=CSUMS(IEEE)+AKT1
CSUMSMOM(IEEE,1)=CSUMSMOM(IEEE,1)+AKT1*PF(1)
CSUMSMOM(IEEE,2)=CSUMSMOM(IEEE,2)+AKT1*PF(2)
!               CSUMSMOM(IEEE,3)=CSUMSMOM(IEEE,3)+AKT1*PF(3)
!               CSUMSMOM(IEEE,4)=CSUMSMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.4).AND.(IDDD.LT.9)) THEN
CSUM2S(IEEE)=CSUM2S(IEEE)+AKT1
CSUM2SMOM(IEEE,1)=CSUM2SMOM(IEEE,1)+AKT1*PF(1)
CSUM2SMOM(IEEE,2)=CSUM2SMOM(IEEE,2)+AKT1*PF(2)
!               CSUM2SMOM(IEEE,3)=CSUM2SMOM(IEEE,3)+AKT1*PF(3)
!               CSUM2SMOM(IEEE,4)=CSUM2SMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.8).AND.(IDDD.LT.13)) THEN
CSUM2WS(IEEE)=CSUM2WS(IEEE)+AKT1
CSUM2WSMOM(IEEE,1)=CSUM2WSMOM(IEEE,1)+AKT1*PF(1)
CSUM2WSMOM(IEEE,2)=CSUM2WSMOM(IEEE,2)+AKT1*PF(2)
!               CSUM2WSMOM(IEEE,3)=CSUM2WSMOM(IEEE,3)+AKT1*PF(3)
!               CSUM2WSMOM(IEEE,4)=CSUM2WSMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.12).AND.(IDDD.LT.17)) THEN
CSUMW(IEEE)=CSUMW(IEEE)+AKT1
CSUMWMOM(IEEE,1)=CSUMWMOM(IEEE,1)+AKT1*PF(1)
CSUMWMOM(IEEE,2)=CSUMWMOM(IEEE,2)+AKT1*PF(2)
!               CSUMWMOM(IEEE,3)=CSUMWMOM(IEEE,3)+AKT1*PF(3)
!               CSUMWMOM(IEEE,4)=CSUMWMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.16).AND.(IDDD.LT.21)) THEN
CSUM2W(IEEE)=CSUM2W(IEEE)+AKT1
CSUM2WMOM(IEEE,1)=CSUM2WMOM(IEEE,1)+AKT1*PF(1)
CSUM2WMOM(IEEE,2)=CSUM2WMOM(IEEE,2)+AKT1*PF(2)
!               CSUM2WMOM(IEEE,3)=CSUM2WMOM(IEEE,3)+AKT1*PF(3)
!               CSUM2WMOM(IEEE,4)=CSUM2WMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.20).AND.(IDDD.LT.25)) THEN
CSUM3W(IEEE)=CSUM3W(IEEE)+AKT1
CSUM3WMOM(IEEE,1)=CSUM3WMOM(IEEE,1)+AKT1*PF(1)
CSUM3WMOM(IEEE,2)=CSUM3WMOM(IEEE,2)+AKT1*PF(2)
!               CSUM3WMOM(IEEE,3)=CSUM3WMOM(IEEE,3)+AKT1*PF(3)
!               CSUM3WMOM(IEEE,4)=CSUM3WMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.24).AND.(IDDD.LT.29)) THEN
CSUMUP(IEEE)=CSUMUP(IEEE)+AKT1
CSUMUPMOM(IEEE,1)=CSUMUPMOM(IEEE,1)+AKT1*PF(1)
CSUMUPMOM(IEEE,2)=CSUMUPMOM(IEEE,2)+AKT1*PF(2)
!               CSUMUPMOM(IEEE,3)=CSUMUPMOM(IEEE,3)+AKT1*PF(3)
!               CSUMUPMOM(IEEE,4)=CSUMUPMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.28).AND.(IDDD.LT.33)) THEN
CSUMUD(IEEE)=CSUMUD(IEEE)+AKT1
CSUMUDMOM(IEEE,1)=CSUMUDMOM(IEEE,1)+AKT1*PF(1)
CSUMUDMOM(IEEE,2)=CSUMUDMOM(IEEE,2)+AKT1*PF(2)
!               CSUMUDMOM(IEEE,3)=CSUMUDMOM(IEEE,3)+AKT1*PF(3)
!               CSUMUDMOM(IEEE,4)=CSUMUDMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.32).AND.(IDDD.LT.41)) THEN
CSUMTT1(IEEE)=CSUMTT1(IEEE)+AKT1
CSUMTTMOM1(IEEE,1)=CSUMTTMOM1(IEEE,1)+AKT1*PF(1)
CSUMTTMOM1(IEEE,2)=CSUMTTMOM1(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM1(IEEE,3)=CSUMTTMOM1(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM1(IEEE,4)=CSUMTTMOM1(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.40).AND.(IDDD.LT.49)) THEN
CSUMTT2(IEEE)=CSUMTT2(IEEE)+AKT1
CSUMTTMOM2(IEEE,1)=CSUMTTMOM2(IEEE,1)+AKT1*PF(1)
CSUMTTMOM2(IEEE,2)=CSUMTTMOM2(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM2(IEEE,3)=CSUMTTMOM2(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM2(IEEE,4)=CSUMTTMOM2(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.48).AND.(IDDD.LT.57)) THEN
CSUMTT3(IEEE)=CSUMTT3(IEEE)+AKT1
CSUMTTMOM3(IEEE,1)=CSUMTTMOM3(IEEE,1)+AKT1*PF(1)
CSUMTTMOM3(IEEE,2)=CSUMTTMOM3(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM3(IEEE,3)=CSUMTTMOM3(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM3(IEEE,4)=CSUMTTMOM3(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.56).AND.(IDDD.LT.65)) THEN
CSUMTT4(IEEE)=CSUMTT4(IEEE)+AKT1
CSUMTTMOM4(IEEE,1)=CSUMTTMOM4(IEEE,1)+AKT1*PF(1)
CSUMTTMOM4(IEEE,2)=CSUMTTMOM4(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM4(IEEE,3)=CSUMTTMOM4(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM4(IEEE,4)=CSUMTTMOM4(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.64).AND.(IDDD.LT.69)) THEN
CSUMTT5(IEEE)=CSUMTT5(IEEE)+AKT1
CSUMTTMOM5(IEEE,1)=CSUMTTMOM5(IEEE,1)+AKT1*PF(1)
CSUMTTMOM5(IEEE,2)=CSUMTTMOM5(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM5(IEEE,3)=CSUMTTMOM5(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM5(IEEE,4)=CSUMTTMOM5(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.68).AND.(IDDD.LT.73)) THEN
CSUMTT6(IEEE)=CSUMTT6(IEEE)+AKT1
CSUMTTMOM6(IEEE,1)=CSUMTTMOM6(IEEE,1)+AKT1*PF(1)
CSUMTTMOM6(IEEE,2)=CSUMTTMOM6(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM6(IEEE,3)=CSUMTTMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM6(IEEE,4)=CSUMTTMOM6(IEEE,4)+AKT1*PF(4)
end if
!            
IF((IDDD.GT.72).AND.(IDDD.LT.81)) THEN
CSUMTT7(IEEE)=CSUMTT7(IEEE)+AKT1
CSUMTTMOM7(IEEE,1)=CSUMTTMOM7(IEEE,1)+AKT1*PF(1)
CSUMTTMOM7(IEEE,2)=CSUMTTMOM7(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM7(IEEE,3)=CSUMTTMOM7(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM7(IEEE,4)=CSUMTTMOM7(IEEE,4)+AKT1*PF(4)
end if
!            
IF((IDDD.GT.80).AND.(IDDD.LT.89)) THEN
CSUMTT8(IEEE)=CSUMTT8(IEEE)+AKT1
CSUMTTMOM8(IEEE,1)=CSUMTTMOM8(IEEE,1)+AKT1*PF(1)
CSUMTTMOM8(IEEE,2)=CSUMTTMOM8(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM8(IEEE,3)=CSUMTTMOM8(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM8(IEEE,4)=CSUMTTMOM8(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.88).AND.(IDDD.LT.97)) THEN
CSUMTT9(IEEE)=CSUMTT9(IEEE)+AKT1
CSUMTTMOM9(IEEE,1)=CSUMTTMOM9(IEEE,1)+AKT1*PF(1)
CSUMTTMOM9(IEEE,2)=CSUMTTMOM9(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM9(IEEE,3)=CSUMTTMOM9(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM9(IEEE,4)=CSUMTTMOM9(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.96).AND.(IDDD.LT.105)) THEN
CSUMTT10(IEEE)=CSUMTT10(IEEE)+AKT1
CSUMTTMOM10(IEEE,1)=CSUMTTMOM10(IEEE,1)+AKT1*PF(1)
CSUMTTMOM10(IEEE,2)=CSUMTTMOM10(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM10(IEEE,3)=CSUMTTMOM10(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM10(IEEE,4)=CSUMTTMOM10(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.104).AND.(IDDD.LT.113)) THEN
CSUMTT11(IEEE)=CSUMTT11(IEEE)+AKT1
CSUMTTMOM11(IEEE,1)=CSUMTTMOM11(IEEE,1)+AKT1*PF(1)
CSUMTTMOM11(IEEE,2)=CSUMTTMOM11(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM11(IEEE,3)=CSUMTTMOM11(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM11(IEEE,4)=CSUMTTMOM11(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.112).AND.(IDDD.LT.129)) THEN
CSUMTT12(IEEE)=CSUMTT12(IEEE)+AKT1
CSUMTTMOM12(IEEE,1)=CSUMTTMOM12(IEEE,1)+AKT1*PF(1)
CSUMTTMOM12(IEEE,2)=CSUMTTMOM12(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM12(IEEE,3)=CSUMTTMOM12(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM12(IEEE,4)=CSUMTTMOM12(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.128).AND.(IDDD.LT.137)) THEN
CSUMTT13(IEEE)=CSUMTT13(IEEE)+AKT1
CSUMTTMOM13(IEEE,1)=CSUMTTMOM13(IEEE,1)+AKT1*PF(1)
CSUMTTMOM13(IEEE,2)=CSUMTTMOM13(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM13(IEEE,3)=CSUMTTMOM13(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM13(IEEE,4)=CSUMTTMOM13(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.136).AND.(IDDD.LT.145)) THEN
CSUMTT14(IEEE)=CSUMTT14(IEEE)+AKT1
CSUMTTMOM14(IEEE,1)=CSUMTTMOM14(IEEE,1)+AKT1*PF(1)
CSUMTTMOM14(IEEE,2)=CSUMTTMOM14(IEEE,2)+AKT1*PF(2)
!               CSUMTTMOM14(IEEE,3)=CSUMTTMOM14(IEEE,3)+AKT1*PF(3)
!               CSUMTTMOM14(IEEE,4)=CSUMTTMOM14(IEEE,4)+AKT1*PF(4)
end if
!
IF(IDDD.EQ.145) THEN
CSUMN=CSUMN+AKT1
end if
!
IF((IDDD.GT.145).AND.(IDDD.LT.154)) THEN
CSUMPLQ(IEEE)=CSUMPLQ(IEEE)+AKT1
CSUMPLQMOM(IEEE,1)=CSUMPLQMOM(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM(IEEE,2)=CSUMPLQMOM(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM(IEEE,3)=CSUMPLQMOM(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM(IEEE,4)=CSUMPLQMOM(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.153).AND.(IDDD.LT.162)) THEN
CSUMPLQ2(IEEE)=CSUMPLQ2(IEEE)+AKT1
CSUMPLQMOM2(IEEE,1)=CSUMPLQMOM2(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM2(IEEE,2)=CSUMPLQMOM2(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM2(IEEE,3)=CSUMPLQMOM2(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM2(IEEE,4)=CSUMPLQMOM2(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.161).AND.(IDDD.LT.170)) THEN
CSUMPLQ3(IEEE)=CSUMPLQ3(IEEE)+AKT1
CSUMPLQMOM3(IEEE,1)=CSUMPLQMOM3(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM3(IEEE,2)=CSUMPLQMOM3(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM3(IEEE,3)=CSUMPLQMOM3(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM3(IEEE,4)=CSUMPLQMOM3(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.169).AND.(IDDD.LT.178)) THEN
CSUMPLQ4(IEEE)=CSUMPLQ4(IEEE)+AKT1
CSUMPLQMOM4(IEEE,1)=CSUMPLQMOM4(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM4(IEEE,2)=CSUMPLQMOM4(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM4(IEEE,3)=CSUMPLQMOM4(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM4(IEEE,4)=CSUMPLQMOM4(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.177).AND.(IDDD.LT.186)) THEN
CSUMPLQ5(IEEE)=CSUMPLQ5(IEEE)+AKT1
CSUMPLQMOM5(IEEE,1)=CSUMPLQMOM5(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM5(IEEE,2)=CSUMPLQMOM5(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM5(IEEE,3)=CSUMPLQMOM5(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM5(IEEE,4)=CSUMPLQMOM5(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.185).AND.(IDDD.LT.194)) THEN
CSUMPLQ6(IEEE)=CSUMPLQ6(IEEE)+AKT1
CSUMPLQMOM6(IEEE,1)=CSUMPLQMOM6(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM6(IEEE,2)=CSUMPLQMOM6(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.193).AND.(IDDD.LT.210)) THEN
CSUMPLQ7(IEEE)=CSUMPLQ7(IEEE)+AKT1
CSUMPLQMOM7(IEEE,1)=CSUMPLQMOM7(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM7(IEEE,2)=CSUMPLQMOM7(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.209).AND.(IDDD.LT.226)) THEN
CSUMPLQ8(IEEE)=CSUMPLQ8(IEEE)+AKT1
CSUMPLQMOM8(IEEE,1)=CSUMPLQMOM8(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM8(IEEE,2)=CSUMPLQMOM8(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.225).AND.(IDDD.LT.242)) THEN
CSUMPLQ9(IEEE)=CSUMPLQ9(IEEE)+AKT1
CSUMPLQMOM9(IEEE,1)=CSUMPLQMOM9(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM9(IEEE,2)=CSUMPLQMOM9(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.241).AND.(IDDD.LT.258)) THEN
CSUMPLQ10(IEEE)=CSUMPLQ10(IEEE)+AKT1
CSUMPLQMOM10(IEEE,1)=CSUMPLQMOM10(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM10(IEEE,2)=CSUMPLQMOM10(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.257).AND.(IDDD.LT.274)) THEN
CSUMPLQ11(IEEE)=CSUMPLQ11(IEEE)+AKT1
CSUMPLQMOM11(IEEE,1)=CSUMPLQMOM11(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM11(IEEE,2)=CSUMPLQMOM11(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.273).AND.(IDDD.LT.290)) THEN
CSUMPLQ12(IEEE)=CSUMPLQ12(IEEE)+AKT1
CSUMPLQMOM12(IEEE,1)=CSUMPLQMOM12(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM12(IEEE,2)=CSUMPLQMOM12(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.289).AND.(IDDD.LT.306)) THEN
CSUMPLQ13(IEEE)=CSUMPLQ13(IEEE)+AKT1
CSUMPLQMOM13(IEEE,1)=CSUMPLQMOM13(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM13(IEEE,2)=CSUMPLQMOM13(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.305).AND.(IDDD.LT.322)) THEN
CSUMPLQ14(IEEE)=CSUMPLQ14(IEEE)+AKT1
CSUMPLQMOM14(IEEE,1)=CSUMPLQMOM14(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM14(IEEE,2)=CSUMPLQMOM14(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!
IF((IDDD.GT.321).AND.(IDDD.LT.338)) THEN
CSUMPLQ15(IEEE)=CSUMPLQ15(IEEE)+AKT1
CSUMPLQMOM15(IEEE,1)=CSUMPLQMOM15(IEEE,1)+AKT1*PF(1)
CSUMPLQMOM15(IEEE,2)=CSUMPLQMOM15(IEEE,2)+AKT1*PF(2)
!               CSUMPLQMOM6(IEEE,3)=CSUMPLQMOM6(IEEE,3)+AKT1*PF(3)
!               CSUMPLQMOM6(IEEE,4)=CSUMPLQMOM6(IEEE,4)+AKT1*PF(4)
end if
!     
9 CONTINUE
end do
!
end do
end do
end do
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
ALINEMOM2(N4,ID,1)=(CSUMSMOM(1,1)+CSUMSMOM(2,1) &
  & +CSUMSMOM(3,1)+CSUMSMOM(4,1))*ADIV1
ALINEMOM2(N4,ID,2)=(CSUMSMOM(1,2)+CSUMSMOM(2,2) &
  & +CSUMSMOM(3,2)+CSUMSMOM(4,2))*ADIV1
!      ALINEMOM2(N4,ID,3)=(CSUMSMOM(1,3)+CSUMSMOM(2,3)
!     &+CSUMSMOM(3,3)+CSUMSMOM(4,3))*ADIV1
!      ALINEMOM2(N4,ID,4)=(CSUMSMOM(1,4)+CSUMSMOM(2,4)
!     &+CSUMSMOM(3,4)+CSUMSMOM(4,4))*ADIV1
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE3(N4,ID)=(CSUMS(1)+GIOT*CSUMS(2)-CSUMS(3)-GIOT*CSUMS(4)) &
  & *ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM3(N4,ID,1)=(CSUMSMOM(1,1)+GIOT*CSUMSMOM(2,1) &
  & -CSUMSMOM(3,1)-GIOT*CSUMSMOM(4,1))*ADIV1
ALINEMOM3(N4,ID,2)=(CSUMSMOM(1,2)+GIOT*CSUMSMOM(2,2) &
  & -CSUMSMOM(3,2)-GIOT*CSUMSMOM(4,2))*ADIV1
!      ALINEMOM3(N4,ID,3)=(CSUMSMOM(1,3)+GIOT*CSUMSMOM(2,3)
!     &-CSUMSMOM(3,3)-GIOT*CSUMSMOM(4,3))*ADIV1
!      ALINEMOM3(N4,ID,4)=(CSUMSMOM(1,4)+GIOT*CSUMSMOM(2,4)
!     &-CSUMSMOM(3,4)-GIOT*CSUMSMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE4(N4,ID)=(CSUMS(1)-CSUMS(2)+CSUMS(3)-CSUMS(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM4(N4,ID,1)=(CSUMSMOM(1,1)-CSUMSMOM(2,1) &
  & +CSUMSMOM(3,1)-CSUMSMOM(4,1))*ADIV1
ALINEMOM4(N4,ID,2)=(CSUMSMOM(1,2)-CSUMSMOM(2,2) &
  & +CSUMSMOM(3,2)-CSUMSMOM(4,2))*ADIV1
!      ALINEMOM4(N4,ID,3)=(CSUMSMOM(1,3)-CSUMSMOM(2,3)
!     &+CSUMSMOM(3,3)-CSUMSMOM(4,3))*ADIV1
!      ALINEMOM4(N4,ID,4)=(CSUMSMOM(1,4)-CSUMSMOM(2,4)
!     &+CSUMSMOM(3,4)-CSUMSMOM(4,4))*ADIV1
!**********************************************************************
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE5(N4,ID)=(CSUM2S(1)+CSUM2S(2)+CSUM2S(3)+CSUM2S(4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM5(N4,ID,1)=(CSUM2SMOM(1,1)+CSUM2SMOM(2,1) &
  & +CSUM2SMOM(3,1)+CSUM2SMOM(4,1))*ADIV1
ALINEMOM5(N4,ID,2)=(CSUM2SMOM(1,2)+CSUM2SMOM(2,2) &
  & +CSUM2SMOM(3,2)+CSUM2SMOM(4,2))*ADIV1
!      ALINEMOM5(N4,ID,3)=(CSUM2SMOM(1,3)+CSUM2SMOM(2,3)
!     &+CSUM2SMOM(3,3)+CSUM2SMOM(4,3))*ADIV1
!      ALINEMOM5(N4,ID,4)=(CSUM2SMOM(1,4)+CSUM2SMOM(2,4)
!     &+CSUM2SMOM(3,4)+CSUM2SMOM(4,4))*ADIV1
!***********************************************************************
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE6(N4,ID)=(CSUM2S(1)+GIOT*CSUM2S(2)-CSUM2S(3) &
  & -GIOT*CSUM2S(4))*ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM6(N4,ID,1)=(CSUM2SMOM(1,1)+GIOT*CSUM2SMOM(2,1) &
  & -CSUM2SMOM(3,1)-GIOT*CSUM2SMOM(4,1))*ADIV1
ALINEMOM6(N4,ID,2)=(CSUM2SMOM(1,2)+GIOT*CSUM2SMOM(2,2) &
  & -CSUM2SMOM(3,2)-GIOT*CSUM2SMOM(4,2))*ADIV1
!      ALINEMOM6(N4,ID,3)=(CSUM2SMOM(1,3)+GIOT*CSUM2SMOM(2,3)
!     &-CSUM2SMOM(3,3)-GIOT*CSUM2SMOM(4,3))*ADIV1
!      ALINEMOM6(N4,ID,4)=(CSUM2SMOM(1,4)+GIOT*CSUM2SMOM(2,4)
!     &-CSUM2SMOM(3,4)-GIOT*CSUM2SMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE7(N4,ID)=(CSUM2S(1)-CSUM2S(2)+CSUM2S(3)-CSUM2S(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM7(N4,ID,1)=(CSUM2SMOM(1,1)-CSUM2SMOM(2,1) &
  & +CSUM2SMOM(3,1)-CSUM2SMOM(4,1))*ADIV1
ALINEMOM7(N4,ID,2)=(CSUM2SMOM(1,2)-CSUM2SMOM(2,2) &
  & +CSUM2SMOM(3,2)-CSUM2SMOM(4,2))*ADIV1
!      ALINEMOM7(N4,ID,3)=(CSUM2SMOM(1,3)-CSUM2SMOM(2,3)
!     &+CSUM2SMOM(3,3)-CSUM2SMOM(4,3))*ADIV1
!      ALINEMOM7(N4,ID,4)=(CSUM2SMOM(1,4)-CSUM2SMOM(2,4)
!     &+CSUM2SMOM(3,4)-CSUM2SMOM(4,4))*ADIV1
!**********************************************************************
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE8(N4,ID)=(CSUM2WS(1)+CSUM2WS(2)+CSUM2WS(3)+CSUM2WS(4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM8(N4,ID,1)=(CSUM2WSMOM(1,1)+CSUM2WSMOM(2,1) &
  & +CSUM2WSMOM(3,1)+CSUM2WSMOM(4,1))*ADIV1
ALINEMOM8(N4,ID,2)=(CSUM2WSMOM(1,2)+CSUM2WSMOM(2,2) &
  & +CSUM2WSMOM(3,2)+CSUM2WSMOM(4,2))*ADIV1
!      ALINEMOM8(N4,ID,3)=(CSUM2WSMOM(1,3)+CSUM2WSMOM(2,3)
!     &+CSUM2WSMOM(3,3)+CSUM2WSMOM(4,3))*ADIV1
!      ALINEMOM8(N4,ID,4)=(CSUM2WSMOM(1,4)+CSUM2WSMOM(2,4)
!     &+CSUM2WSMOM(3,4)+CSUM2WSMOM(4,4))*ADIV1
!************************************************************************
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE9(N4,ID)=(CSUM2WS(1)+GIOT*CSUM2WS(2)-CSUM2WS(3) &
  & -GIOT*CSUM2WS(4))*ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM9(N4,ID,1)=(CSUM2WSMOM(1,1)+GIOT*CSUM2WSMOM(2,1) &
  & -CSUM2WSMOM(3,1)-GIOT*CSUM2WSMOM(4,1))*ADIV1
ALINEMOM9(N4,ID,2)=(CSUM2WSMOM(1,2)+GIOT*CSUM2WSMOM(2,2) &
  & -CSUM2WSMOM(3,2)-GIOT*CSUM2WSMOM(4,2))*ADIV1
!      ALINEMOM9(N4,ID,3)=(CSUM2WSMOM(1,3)+GIOT*CSUM2WSMOM(2,3)
!     &-CSUM2WSMOM(3,3)-GIOT*CSUM2WSMOM(4,3))*ADIV1
!      ALINEMOM9(N4,ID,4)=(CSUM2WSMOM(1,4)+GIOT*CSUM2WSMOM(2,4)
!     &-CSUM2WSMOM(3,4)-GIOT*CSUM2WSMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE10(N4,ID)=(CSUM2WS(1)-CSUM2WS(2)+CSUM2WS(3)-CSUM2WS(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM10(N4,ID,1)=(CSUM2WSMOM(1,1)-CSUM2WSMOM(2,1) &
  & +CSUM2WSMOM(3,1)-CSUM2WSMOM(4,1))*ADIV1
ALINEMOM10(N4,ID,2)=(CSUM2WSMOM(1,2)-CSUM2WSMOM(2,2) &
  & +CSUM2WSMOM(3,2)-CSUM2WSMOM(4,2))*ADIV1
!      ALINEMOM10(N4,ID,3)=(CSUM2WSMOM(1,3)-CSUM2WSMOM(2,3)
!     &+CSUM2WSMOM(3,3)-CSUM2WSMOM(4,3))*ADIV1
!      ALINEMOM10(N4,ID,4)=(CSUM2WSMOM(1,4)-CSUM2WSMOM(2,4)
!     &+CSUM2WSMOM(3,4)-CSUM2WSMOM(4,4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE11(N4,ID)=(CSUMW(1)+CSUMW(2)+CSUMW(3)+CSUMW(4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM11(N4,ID,1)=(CSUMWMOM(1,1)+CSUMWMOM(2,1) &
  & +CSUMWMOM(3,1)+CSUMWMOM(4,1))*ADIV1
ALINEMOM11(N4,ID,2)=(CSUMWMOM(1,2)+CSUMWMOM(2,2) &
  & +CSUMWMOM(3,2)+CSUMWMOM(4,2))*ADIV1
!      ALINEMOM11(N4,ID,3)=(CSUMWMOM(1,3)+CSUMWMOM(2,3)
!     &+CSUMWMOM(3,3)+CSUMWMOM(4,3))*ADIV1
!      ALINEMOM11(N4,ID,4)=(CSUMWMOM(1,4)+CSUMWMOM(2,4)
!     &+CSUMWMOM(3,4)+CSUMWMOM(4,4))*ADIV1
!***********************************************************************
!     J=1, Pr=-, q=0
!***********************************************************************
ALINE12(N4,ID)=(CSUMW(1)+GIOT*CSUMW(2)-CSUMW(3) &
  & -GIOT*CSUMW(4))*ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM12(N4,ID,1)=(CSUMWMOM(1,1)+GIOT*CSUMWMOM(2,1) &
  & -CSUMWMOM(3,1)-GIOT*CSUMWMOM(4,1))*ADIV1
ALINEMOM12(N4,ID,2)=(CSUMWMOM(1,2)+GIOT*CSUMWMOM(2,2) &
  & -CSUMWMOM(3,2)-GIOT*CSUMWMOM(4,2))*ADIV1
!      ALINEMOM12(N4,ID,3)=(CSUMWMOM(1,3)+GIOT*CSUMWMOM(2,3)
!     &-CSUMWMOM(3,3)-GIOT*CSUMWMOM(4,3))*ADIV1
!      ALINEMOM12(N4,ID,4)=(CSUMWMOM(1,4)+GIOT*CSUMWMOM(2,4)
!     &-CSUMWMOM(3,4)-GIOT*CSUMWMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE13(N4,ID)=(CSUMW(1)-CSUMW(2)+CSUMW(3)-CSUMW(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM13(N4,ID,1)=(CSUMWMOM(1,1)-CSUMWMOM(2,1) &
  & +CSUMWMOM(3,1)-CSUMWMOM(4,1))*ADIV1
ALINEMOM13(N4,ID,2)=(CSUMWMOM(1,2)-CSUMWMOM(2,2) &
  & +CSUMWMOM(3,2)-CSUMWMOM(4,2))*ADIV1
!      ALINEMOM13(N4,ID,3)=(CSUMWMOM(1,3)-CSUMWMOM(2,3)
!     &+CSUMWMOM(3,3)-CSUMWMOM(4,3))*ADIV1
!      ALINEMOM13(N4,ID,4)=(CSUMWMOM(1,4)-CSUMWMOM(2,4)
!     &+CSUMWMOM(3,4)-CSUMWMOM(4,4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE14(N4,ID)=(CSUM2W(1)+CSUM2W(2)+CSUM2W(3)+CSUM2W(4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM14(N4,ID,1)=(CSUM2WMOM(1,1)+CSUM2WMOM(2,1) &
  & +CSUM2WMOM(3,1)+CSUM2WMOM(4,1))*ADIV1
ALINEMOM14(N4,ID,2)=(CSUM2WMOM(1,2)+CSUM2WMOM(2,2) &
  & +CSUM2WMOM(3,2)+CSUM2WMOM(4,2))*ADIV1
!      ALINEMOM14(N4,ID,3)=(CSUM2WMOM(1,3)+CSUM2WMOM(2,3)
!     &+CSUM2WMOM(3,3)+CSUM2WMOM(4,3))*ADIV1
!      ALINEMOM14(N4,ID,4)=(CSUM2WMOM(1,4)+CSUM2WMOM(2,4)
!     &+CSUM2WMOM(3,4)+CSUM2WMOM(4,4))*ADIV1
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE15(N4,ID)=(CSUM2W(1)+GIOT*CSUM2W(2)-CSUM2W(3) &
  & -GIOT*CSUM2W(4))*ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM15(N4,ID,1)=(CSUM2WMOM(1,1)+GIOT*CSUM2WMOM(2,1) &
  & -CSUM2WMOM(3,1)-GIOT*CSUM2WMOM(4,1))*ADIV1
ALINEMOM15(N4,ID,2)=(CSUM2WMOM(1,2)+GIOT*CSUM2WMOM(2,2) &
  & -CSUM2WMOM(3,2)-GIOT*CSUM2WMOM(4,2))*ADIV1
!      ALINEMOM15(N4,ID,3)=(CSUM2WMOM(1,3)+GIOT*CSUM2WMOM(2,3)
!     &-CSUM2WMOM(3,3)-GIOT*CSUM2WMOM(4,3))*ADIV1
!      ALINEMOM15(N4,ID,4)=(CSUM2WMOM(1,4)+GIOT*CSUM2WMOM(2,4)
!     &-CSUM2WMOM(3,4)-GIOT*CSUM2WMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE16(N4,ID)=(CSUM2W(1)-CSUM2W(2)+CSUM2W(3)-CSUM2W(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM16(N4,ID,1)=(CSUM2WMOM(1,1)-CSUM2WMOM(2,1) &
  & +CSUM2WMOM(3,1)-CSUM2WMOM(4,1))*ADIV1
ALINEMOM16(N4,ID,2)=(CSUM2WMOM(1,2)-CSUM2WMOM(2,2) &
  & +CSUM2WMOM(3,2)-CSUM2WMOM(4,2))*ADIV1
!      ALINEMOM16(N4,ID,3)=(CSUM2WMOM(1,3)-CSUM2WMOM(2,3)
!     &+CSUM2WMOM(3,3)-CSUM2WMOM(4,3))*ADIV1
!      ALINEMOM16(N4,ID,4)=(CSUM2WMOM(1,4)-CSUM2WMOM(2,4)
!     &+CSUM2WMOM(3,4)-CSUM2WMOM(4,4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE17(N4,ID)=(CSUM3W(1)+CSUM3W(2)+CSUM3W(3)+CSUM3W(4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM17(N4,ID,1)=(CSUM3WMOM(1,1)+CSUM3WMOM(2,1) &
  & +CSUM3WMOM(3,1)+CSUM3WMOM(4,1))*ADIV1
ALINEMOM17(N4,ID,2)=(CSUM3WMOM(1,2)+CSUM3WMOM(2,2) &
  & +CSUM3WMOM(3,2)+CSUM3WMOM(4,2))*ADIV1
!      ALINEMOM17(N4,ID,3)=(CSUM3WMOM(1,3)+CSUM3WMOM(2,3)
!     &+CSUM3WMOM(3,3)+CSUM3WMOM(4,3))*ADIV1
!      ALINEMOM17(N4,ID,4)=(CSUM3WMOM(1,4)+CSUM3WMOM(2,4)
!     &+CSUM3WMOM(3,4)+CSUM3WMOM(4,4))*ADIV1
!***********************************************************************
!     J=1, Pr=+ q=0
!**********************************************************************
ALINE18(N4,ID)=(CSUM3W(1)+GIOT*CSUM3W(2)-CSUM3W(3)-GIOT*CSUM3W(4)) &
  & *ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM18(N4,ID,1)=(CSUM3WMOM(1,1)+GIOT*CSUM3WMOM(2,1) &
  & -CSUM3WMOM(3,1)-GIOT*CSUM3WMOM(4,1))*ADIV1
ALINEMOM18(N4,ID,2)=(CSUM3WMOM(1,2)+GIOT*CSUM3WMOM(2,2) &
  & -CSUM3WMOM(3,2)-GIOT*CSUM3WMOM(4,2))*ADIV1
!      ALINEMOM18(N4,ID,3)=(CSUM3WMOM(1,3)+GIOT*CSUM3WMOM(2,3)
!     &-CSUM3WMOM(3,3)-GIOT*CSUM3WMOM(4,3))*ADIV1
!      ALINEMOM18(N4,ID,4)=(CSUM3WMOM(1,4)+GIOT*CSUM3WMOM(2,4)
!     &-CSUM3WMOM(3,4)-GIOT*CSUM3WMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE19(N4,ID)=(CSUM3W(1)-CSUM3W(2)+CSUM3W(3)-CSUM3W(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM19(N4,ID,1)=(CSUM3WMOM(1,1)-CSUM3WMOM(2,1) &
  & +CSUM3WMOM(3,1)-CSUM3WMOM(4,1))*ADIV1
ALINEMOM19(N4,ID,2)=(CSUM3WMOM(1,2)-CSUM3WMOM(2,2) &
  & +CSUM3WMOM(3,2)-CSUM3WMOM(4,2))*ADIV1
!      ALINEMOM19(N4,ID,3)=(CSUM3WMOM(1,3)-CSUM3WMOM(2,3)
!     &+CSUM3WMOM(3,3)-CSUM3WMOM(4,3))*ADIV1
!      ALINEMOM19(N4,ID,4)=(CSUM3WMOM(1,4)-CSUM3WMOM(2,4)
!     &+CSUM3WMOM(3,4)-CSUM3WMOM(4,4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE20(N4,ID)=(CSUMUP(1)+CSUMUP(2)+CSUMUP(3)+CSUMUP(4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM20(N4,ID,1)=(CSUMUPMOM(1,1)+CSUMUPMOM(2,1) &
  & +CSUMUPMOM(3,1)+CSUMUPMOM(4,1))*ADIV1
ALINEMOM20(N4,ID,2)=(CSUMUPMOM(1,2)+CSUMUPMOM(2,2) &
  & +CSUMUPMOM(3,2)+CSUMUPMOM(4,2))*ADIV1
!      ALINEMOM20(N4,ID,3)=(CSUMUPMOM(1,3)+CSUMUPMOM(2,3)
!     &+CSUMUPMOM(3,3)+CSUMUPMOM(4,3))*ADIV1
!      ALINEMOM20(N4,ID,4)=(CSUMUPMOM(1,4)+CSUMUPMOM(2,4)
!     &+CSUMUPMOM(3,4)+CSUMUPMOM(4,4))*ADIV1
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE21(N4,ID)=(CSUMUP(1)+GIOT*CSUMUP(2)-CSUMUP(3)-GIOT*CSUMUP(4)) &
  & *ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM21(N4,ID,1)=(CSUMUPMOM(1,1)+GIOT*CSUMUPMOM(2,1) &
  & -CSUMUPMOM(3,1)-GIOT*CSUMUPMOM(4,1))*ADIV1
ALINEMOM21(N4,ID,2)=(CSUMUPMOM(1,2)+GIOT*CSUMUPMOM(2,2) &
  & -CSUMUPMOM(3,2)-GIOT*CSUMUPMOM(4,2))*ADIV1
!      ALINEMOM21(N4,ID,3)=(CSUMUPMOM(1,3)+GIOT*CSUMUPMOM(2,3)
!     &-CSUMUPMOM(3,3)-GIOT*CSUMUPMOM(4,3))*ADIV1
!      ALINEMOM21(N4,ID,4)=(CSUMUPMOM(1,4)+GIOT*CSUMUPMOM(2,4)
!     &-CSUMUPMOM(3,4)-GIOT*CSUMUPMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE22(N4,ID)=(CSUMUP(1)-CSUMUP(2)+CSUMUP(3)-CSUMUP(4))*ADIV1
!***********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!***********************************************************************
ALINEMOM22(N4,ID,1)=(CSUMUPMOM(1,1)-CSUMUPMOM(2,1) &
  & +CSUMUPMOM(3,1)-CSUMUPMOM(4,1))*ADIV1
ALINEMOM22(N4,ID,2)=(CSUMUPMOM(1,2)-CSUMUPMOM(2,2) &
  & +CSUMUPMOM(3,2)-CSUMUPMOM(4,2))*ADIV1
!      ALINEMOM22(N4,ID,3)=(CSUMUPMOM(1,3)-CSUMUPMOM(2,3)
!     &+CSUMUPMOM(3,3)-CSUMUPMOM(4,3))*ADIV1
!      ALINEMOM22(N4,ID,4)=(CSUMUPMOM(1,4)-CSUMUPMOM(2,4)
!     &+CSUMUPMOM(3,4)-CSUMUPMOM(4,4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE23(N4,ID)=(CSUMUD(1)+CSUMUD(2)+CSUMUD(3)+CSUMUD(4))*ADIV1
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM23(N4,ID,1)=(CSUMUDMOM(1,1)+CSUMUDMOM(2,1) &
  & +CSUMUDMOM(3,1)+CSUMUDMOM(4,1))*ADIV1
ALINEMOM23(N4,ID,2)=(CSUMUDMOM(1,2)+CSUMUDMOM(2,2) &
  & +CSUMUDMOM(3,2)+CSUMUDMOM(4,2))*ADIV1
!      ALINEMOM23(N4,ID,3)=(CSUMUDMOM(1,3)+CSUMUDMOM(2,3)
!     &+CSUMUDMOM(3,3)+CSUMUDMOM(4,3))*ADIV1
!      ALINEMOM23(N4,ID,4)=(CSUMUDMOM(1,4)+CSUMUDMOM(2,4)
!     &+CSUMUDMOM(3,4)+CSUMUDMOM(4,4))*ADIV1
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE24(N4,ID)=(CSUMUD(1)+GIOT*CSUMUD(2)-CSUMUD(3) &
  & -GIOT*CSUMUD(4))*ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!***********************************************************************
ALINEMOM24(N4,ID,1)=(CSUMUDMOM(1,1)+GIOT*CSUMUDMOM(2,1) &
  & -CSUMUDMOM(3,1)-GIOT*CSUMUDMOM(4,1))*ADIV1
ALINEMOM24(N4,ID,2)=(CSUMUDMOM(1,2)+GIOT*CSUMUDMOM(2,2) &
  & -CSUMUDMOM(3,2)-GIOT*CSUMUDMOM(4,2))*ADIV1
!      ALINEMOM24(N4,ID,3)=(CSUMUDMOM(1,3)+GIOT*CSUMUDMOM(2,3)
!     &-CSUMUDMOM(3,3)-GIOT*CSUMUDMOM(4,3))*ADIV1
!      ALINEMOM24(N4,ID,4)=(CSUMUDMOM(1,4)+GIOT*CSUMUDMOM(2,4)
!     &-CSUMUDMOM(3,4)-GIOT*CSUMUDMOM(4,4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE25(N4,ID)=(CSUMUD(1)-CSUMUD(2)+CSUMUD(3)-CSUMUD(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
ALINEMOM25(N4,ID,1)=(CSUMUDMOM(1,1)-CSUMUDMOM(2,1) &
  & +CSUMUDMOM(3,1)-CSUMUDMOM(4,1))*ADIV1
ALINEMOM25(N4,ID,2)=(CSUMUDMOM(1,2)-CSUMUDMOM(2,2) &
  & +CSUMUDMOM(3,2)-CSUMUDMOM(4,2))*ADIV1
!      ALINEMOM25(N4,ID,3)=(CSUMUDMOM(1,3)-CSUMUDMOM(2,3)
!     &+CSUMUDMOM(3,3)-CSUMUDMOM(4,3))*ADIV1
!      ALINEMOM25(N4,ID,4)=(CSUMUDMOM(1,4)-CSUMUDMOM(2,4)
!     &+CSUMUDMOM(3,4)-CSUMUDMOM(4,4))*ADIV1
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE26(N4,ID)=(CSUMTT1(1)+CSUMTT1(2)+CSUMTT1(3)+CSUMTT1(4)+ &
  & (CSUMTT1(7)+CSUMTT1(6)+CSUMTT1(5)+CSUMTT1(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM26(N4,ID,IK)=(CSUMTTMOM1(1,IK)+CSUMTTMOM1(2,IK) &
  & +CSUMTTMOM1(3,IK)+CSUMTTMOM1(4,IK)+CSUMTTMOM1(7,IK) &
  & +CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK) &
  & +CSUMTTMOM1(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE27(N4,ID)=(CSUMTT1(1)+CSUMTT1(2)+CSUMTT1(3)+CSUMTT1(4)- &
  & (CSUMTT1(7)+CSUMTT1(6)+CSUMTT1(5)+CSUMTT1(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM27(N4,ID,IK)=(CSUMTTMOM1(1,IK)+CSUMTTMOM1(2,IK) &
  & +CSUMTTMOM1(3,IK)+CSUMTTMOM1(4,IK)-(CSUMTTMOM1(7,IK) &
  & +CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK) &
  & +CSUMTTMOM1(8,IK)))*ADIV2
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE28(N4,ID)=(CSUMTT1(1)+GIOT*CSUMTT1(2)-CSUMTT1(3) &
  & -GIOT*CSUMTT1(4)+(CSUMTT1(6)+GIOT*CSUMTT1(7)-CSUMTT1(8) &
  & -GIOT*CSUMTT1(5)))*ADIV2
!***********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM28(N4,ID,IK)=(CSUMTTMOM1(1,IK)+GIOT*CSUMTTMOM1(2,IK) &
  & -CSUMTTMOM1(3,IK)-GIOT*CSUMTTMOM1(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE29(N4,ID)=(CSUMTT1(1)+GIOT*CSUMTT1(2)-CSUMTT1(3) &
  & -GIOT*CSUMTT1(4)-(CSUMTT1(6)+GIOT*CSUMTT1(7)-CSUMTT1(8) &
  & -GIOT*CSUMTT1(5)))*ADIV2
!***********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM29(N4,ID,IK)=(CSUMTTMOM1(6,IK)+GIOT*CSUMTTMOM1(7,IK) &
  & -CSUMTTMOM1(8,IK)-GIOT*CSUMTTMOM1(5,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE30(N4,ID)=(CSUMTT1(1)-CSUMTT1(2)+CSUMTT1(3)-CSUMTT1(4)+ &
  & CSUMTT1(7)-CSUMTT1(6)+CSUMTT1(5)-CSUMTT1(8))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM30(N4,ID,IK)=(CSUMTTMOM1(1,IK)-CSUMTTMOM1(2,IK) &
  & +CSUMTTMOM1(3,IK)-CSUMTTMOM1(4,IK)+CSUMTTMOM1(7,IK) &
  & -CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)-CSUMTTMOM1(8,IK))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE31(N4,ID)=(CSUMTT1(1)-CSUMTT1(2)+CSUMTT1(3)-CSUMTT1(4)- &
  & (CSUMTT1(7)-CSUMTT1(6)+CSUMTT1(5)-CSUMTT1(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM31(N4,ID,IK)=(CSUMTTMOM1(1,IK)-CSUMTTMOM1(2,IK) &
  & +CSUMTTMOM1(3,IK)-CSUMTTMOM1(4,IK)-(CSUMTTMOM1(7,IK) &
  & -CSUMTTMOM1(6,IK)+CSUMTTMOM1(5,IK)-CSUMTTMOM1(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE32(N4,ID)=(CSUMTT2(1)+CSUMTT2(2)+CSUMTT2(3)+CSUMTT2(4)+ &
  & CSUMTT2(5)+CSUMTT2(6)+CSUMTT2(7)+CSUMTT2(8))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM32(N4,ID,IK)=(CSUMTTMOM2(1,IK)+CSUMTTMOM2(2,IK) &
  & +CSUMTTMOM2(3,IK)+CSUMTTMOM2(4,IK)+CSUMTTMOM2(5,IK) &
  & +CSUMTTMOM2(6,IK)+CSUMTTMOM2(7,IK) &
  & +CSUMTTMOM2(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE33(N4,ID)=(CSUMTT2(1)+CSUMTT2(2)+CSUMTT2(3)+CSUMTT2(4)- &
  & (CSUMTT2(5)+CSUMTT2(6)+CSUMTT2(7)+CSUMTT2(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q
!**********************************************************************
DO IK=1,2
ALINEMOM33(N4,ID,IK)=(CSUMTTMOM2(1,IK)+CSUMTTMOM2(2,IK) &
  & +CSUMTTMOM2(3,IK)+CSUMTTMOM2(4,IK)-(CSUMTTMOM2(5,IK) &
  & +CSUMTTMOM2(6,IK)+CSUMTTMOM2(7,IK) &
  & +CSUMTTMOM2(8,IK)))*ADIV2
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE34(N4,ID)=(CSUMTT2(1)+GIOT*CSUMTT2(2)-CSUMTT2(3) &
  & -GIOT*CSUMTT2(4)+(CSUMTT2(7)+GIOT*CSUMTT2(8)-CSUMTT2(5) &
  & -GIOT*CSUMTT2(6)))*ADIV2
!***********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM34(N4,ID,IK)=(CSUMTTMOM2(1,IK)+GIOT*CSUMTTMOM2(2,IK) &
  & -CSUMTTMOM2(3,IK)-GIOT*CSUMTTMOM2(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE35(N4,ID)=(CSUMTT2(1)+GIOT*CSUMTT2(2)-CSUMTT2(3) &
  & -GIOT*CSUMTT2(4)-(CSUMTT2(7)+GIOT*CSUMTT2(8)-CSUMTT2(5) &
  & -GIOT*CSUMTT2(6)))*ADIV2
!***********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM35(N4,ID,IK)=(CSUMTTMOM2(7,IK)+GIOT*CSUMTTMOM2(8,IK) &
  & -CSUMTTMOM2(5,IK)-GIOT*CSUMTTMOM2(6,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE36(N4,ID)=(CSUMTT2(1)-CSUMTT2(2)+CSUMTT2(3)-CSUMTT2(4)+ &
  & CSUMTT2(6)-CSUMTT2(5)+CSUMTT2(8)-CSUMTT2(7))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM36(N4,ID,IK)=(CSUMTTMOM2(1,IK)-CSUMTTMOM2(2,IK) &
  & +CSUMTTMOM2(3,IK)-CSUMTTMOM2(4,IK)+CSUMTTMOM2(6,IK) &
  & -CSUMTTMOM2(5,IK)+CSUMTTMOM2(8,IK) &
  & -CSUMTTMOM2(7,IK))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE37(N4,ID)=(CSUMTT2(1)-CSUMTT2(2)+CSUMTT2(3)-CSUMTT2(4)- &
  & (CSUMTT2(6)-CSUMTT2(5)+CSUMTT2(8)-CSUMTT2(7)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM37(N4,ID,IK)=(CSUMTTMOM2(1,IK)-CSUMTTMOM2(2,IK) &
  & +CSUMTTMOM2(3,IK)-CSUMTTMOM2(4,IK)-(CSUMTTMOM2(6,IK) &
  & -CSUMTTMOM2(5,IK)+CSUMTTMOM2(8,IK) &
  & -CSUMTTMOM2(7,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE38(N4,ID)=(CSUMTT3(1)+CSUMTT3(2)+CSUMTT3(3)+CSUMTT3(4)+ &
  & CSUMTT3(5)+CSUMTT3(6)+CSUMTT3(7)+CSUMTT3(8))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM38(N4,ID,IK)=(CSUMTTMOM3(1,IK)+CSUMTTMOM3(2,IK) &
  & +CSUMTTMOM3(3,IK)+CSUMTTMOM3(4,IK)+CSUMTTMOM3(5,IK) &
  & +CSUMTTMOM3(6,IK)+CSUMTTMOM3(7,IK) &
  & +CSUMTTMOM3(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE39(N4,ID)=(CSUMTT3(1)+CSUMTT3(2)+CSUMTT3(3)+CSUMTT3(4)- &
  & (CSUMTT3(5)+CSUMTT3(6)+CSUMTT3(7)+CSUMTT3(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM39(N4,ID,IK)=(CSUMTTMOM3(1,IK)+CSUMTTMOM3(2,IK) &
  & +CSUMTTMOM3(3,IK)+CSUMTTMOM3(4,IK)-(CSUMTTMOM3(5,IK) &
  & +CSUMTTMOM3(6,IK)+CSUMTTMOM3(7,IK) &
  & +CSUMTTMOM3(8,IK)))*ADIV2
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE40(N4,ID)=(CSUMTT3(1)+GIOT*CSUMTT3(2)-CSUMTT3(3) &
  & -GIOT*CSUMTT3(4)+(CSUMTT3(6)+GIOT*CSUMTT3(7)-CSUMTT3(8) &
  & -GIOT*CSUMTT3(5)))*ADIV2
!***********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM40(N4,ID,IK)=(CSUMTTMOM3(1,IK)+GIOT*CSUMTTMOM3(2,IK) &
  & -CSUMTTMOM3(3,IK)-GIOT*CSUMTTMOM3(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE41(N4,ID)=(CSUMTT3(1)+GIOT*CSUMTT3(2)-CSUMTT3(3) &
  & -GIOT*CSUMTT3(4)-(CSUMTT3(6)+GIOT*CSUMTT3(7)-CSUMTT3(8) &
  & -GIOT*CSUMTT3(5)))*ADIV2
!***********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM41(N4,ID,IK)=(CSUMTTMOM3(6,IK)+GIOT*CSUMTTMOM3(7,IK) &
  & -CSUMTTMOM3(8,IK)-GIOT*CSUMTTMOM3(5,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE42(N4,ID)=(CSUMTT3(1)-CSUMTT3(2)+CSUMTT3(3)-CSUMTT3(4)+ &
  & (CSUMTT3(7)-CSUMTT3(6)+CSUMTT3(5)-CSUMTT3(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM42(N4,ID,IK)=(CSUMTTMOM3(1,IK)-CSUMTTMOM3(2,IK) &
  & +CSUMTTMOM3(3,IK)-CSUMTTMOM3(4,IK)+(CSUMTTMOM3(7,IK) &
  & -CSUMTTMOM3(6,IK)+CSUMTTMOM3(5,IK)-CSUMTTMOM3(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE43(N4,ID)=(CSUMTT3(1)-CSUMTT3(2)+CSUMTT3(3)-CSUMTT3(4)- &
  & (CSUMTT3(7)-CSUMTT3(6)+CSUMTT3(5)-CSUMTT3(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM43(N4,ID,IK)=(CSUMTTMOM3(1,IK)-CSUMTTMOM3(2,IK) &
  & +CSUMTTMOM3(3,IK)-CSUMTTMOM3(4,IK)-(CSUMTTMOM3(7,IK) &
  & -CSUMTTMOM3(6,IK)+CSUMTTMOM3(5,IK)-CSUMTTMOM3(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE44(N4,ID)=(CSUMTT4(1)+CSUMTT4(2)+CSUMTT4(3)+CSUMTT4(4)+ &
  & CSUMTT4(5)+CSUMTT4(6)+CSUMTT4(7)+CSUMTT4(8))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM44(N4,ID,IK)=(CSUMTTMOM4(1,IK)+CSUMTTMOM4(2,IK) &
  & +CSUMTTMOM4(3,IK)+CSUMTTMOM4(4,IK)+CSUMTTMOM4(5,IK) &
  & +CSUMTTMOM4(6,IK)+CSUMTTMOM4(7,IK) &
  & +CSUMTTMOM4(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE45(N4,ID)=(CSUMTT4(1)+CSUMTT4(2)+CSUMTT4(3)+CSUMTT4(4)- &
  & (CSUMTT4(5)+CSUMTT4(6)+CSUMTT4(7)+CSUMTT4(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM45(N4,ID,IK)=(CSUMTTMOM4(1,IK)+CSUMTTMOM4(2,IK) &
  & +CSUMTTMOM4(3,IK)+CSUMTTMOM4(4,IK)-(CSUMTTMOM4(5,IK) &
  & +CSUMTTMOM4(6,IK)+CSUMTTMOM4(7,IK) &
  & +CSUMTTMOM4(8,IK)))*ADIV2
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE46(N4,ID)=(CSUMTT4(1)+GIOT*CSUMTT4(2)-CSUMTT4(3) &
  & -GIOT*CSUMTT4(4)+(CSUMTT4(7)+GIOT*CSUMTT4(8)-CSUMTT4(5) &
  & -GIOT*CSUMTT4(6)))*ADIV2
!***********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM46(N4,ID,IK)=(CSUMTTMOM4(1,IK)+GIOT*CSUMTTMOM4(2,IK) &
  & -CSUMTTMOM4(3,IK)-GIOT*CSUMTTMOM4(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE47(N4,ID)=(CSUMTT4(1)+GIOT*CSUMTT4(2)-CSUMTT4(3) &
  & -GIOT*CSUMTT4(4)-(CSUMTT4(7)+GIOT*CSUMTT4(8)-CSUMTT4(5) &
  & -GIOT*CSUMTT4(6)))*ADIV2
!***********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM47(N4,ID,IK)=(CSUMTTMOM4(7,IK)+GIOT*CSUMTTMOM4(8,IK) &
  & -CSUMTTMOM4(5,IK)-GIOT*CSUMTTMOM4(6,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE48(N4,ID)=(CSUMTT4(1)-CSUMTT4(2)+CSUMTT4(3)-CSUMTT4(4)+ &
  & CSUMTT4(8)-CSUMTT4(7)+CSUMTT4(6)-CSUMTT4(5))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM48(N4,ID,IK)=(CSUMTTMOM4(1,IK)-CSUMTTMOM4(2,IK) &
  & +CSUMTTMOM4(3,IK)-CSUMTTMOM4(4,IK)+CSUMTTMOM4(8,IK) &
  & -CSUMTTMOM4(7,IK)+CSUMTTMOM4(6,IK)-CSUMTTMOM4(5,IK))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE49(N4,ID)=(CSUMTT4(1)-CSUMTT4(2)+CSUMTT4(3)-CSUMTT4(4)- &
  & (CSUMTT4(8)-CSUMTT4(7)+CSUMTT4(6)-CSUMTT4(5)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM49(N4,ID,IK)=(CSUMTTMOM4(1,IK)-CSUMTTMOM4(2,IK) &
  & +CSUMTTMOM4(3,IK)-CSUMTTMOM4(4,IK)-(CSUMTTMOM4(8,IK) &
  & -CSUMTTMOM4(7,IK)+CSUMTTMOM4(6,IK)-CSUMTTMOM4(5,IK)))*ADIV2
end do
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
DO IK=1,2
ALINEMOM50(N4,ID,IK)=(CSUMTTMOM5(1,IK)+CSUMTTMOM5(2,IK) &
  & +CSUMTTMOM5(3,IK)+CSUMTTMOM5(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE51(N4,ID)=(CSUMTT5(1)+GIOT*CSUMTT5(2)-CSUMTT5(3) &
  & -GIOT*CSUMTT5(4))*ADIV1
!***********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM51(N4,ID,IK)=(CSUMTTMOM5(1,IK)+GIOT*CSUMTTMOM5(2,IK) &
  & -CSUMTTMOM5(3,IK)-GIOT*CSUMTTMOM5(4,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE52(N4,ID)=(CSUMTT5(1)-CSUMTT5(2)+CSUMTT5(3)-CSUMTT5(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM52(N4,ID,IK)=(CSUMTTMOM5(1,IK)-CSUMTTMOM5(2,IK) &
  & +CSUMTTMOM5(3,IK)-CSUMTTMOM5(4,IK))*ADIV1
end do
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
DO IK=1,2
ALINEMOM53(N4,ID,IK)=(CSUMTTMOM6(1,IK)+CSUMTTMOM6(2,IK) &
  & +CSUMTTMOM6(3,IK)+CSUMTTMOM6(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE54(N4,ID)=(CSUMTT6(1)+GIOT*CSUMTT6(2)-CSUMTT6(3) &
  & -GIOT*CSUMTT6(4))*ADIV1
!***********************************************************************
!     J=1, q
!**********************************************************************
DO IK=1,2
ALINEMOM54(N4,ID,IK)=(CSUMTTMOM6(1,IK)+GIOT*CSUMTTMOM6(2,IK) &
  & -CSUMTTMOM6(3,IK)-GIOT*CSUMTTMOM6(4,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE55(N4,ID)=(CSUMTT6(1)-CSUMTT6(2)+CSUMTT6(3)-CSUMTT6(4))*ADIV1
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM55(N4,ID,IK)=(CSUMTTMOM6(1,IK)-CSUMTTMOM6(2,IK) &
  & +CSUMTTMOM6(3,IK)-CSUMTTMOM6(4,IK))*ADIV1
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 7
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-7 OPERATORS 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE56(N4,ID)=(CSUMTT7(1)+CSUMTT7(2)+CSUMTT7(3)+CSUMTT7(4)+ &
  & (CSUMTT7(5)+CSUMTT7(6)+CSUMTT7(7)+CSUMTT7(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM56(N4,ID,IK)=(CSUMTTMOM7(1,IK)+CSUMTTMOM7(2,IK) &
  & +CSUMTTMOM7(3,IK)+CSUMTTMOM7(4,IK)+(CSUMTTMOM7(5,IK) &
  & +CSUMTTMOM7(6,IK)+CSUMTTMOM7(7,IK) &
  & +CSUMTTMOM7(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE57(N4,ID)=(CSUMTT7(1)+CSUMTT7(2)+CSUMTT7(3)+CSUMTT7(4)- &
  & (CSUMTT7(5)+CSUMTT7(6)+CSUMTT7(7)+CSUMTT7(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM57(N4,ID,IK)=(CSUMTTMOM7(1,IK)+CSUMTTMOM7(2,IK) &
  & +CSUMTTMOM7(3,IK)+CSUMTTMOM7(4,IK)-(CSUMTTMOM7(5,IK) &
  & +CSUMTTMOM7(6,IK)+CSUMTTMOM7(7,IK) &
  & +CSUMTTMOM7(8,IK)))*ADIV2
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE58(N4,ID)=(CSUMTT7(1)+GIOT*CSUMTT7(2)-CSUMTT7(3) &
  & -GIOT*CSUMTT7(4)+ &
  & (CSUMTT7(7)+GIOT*CSUMTT7(8)-CSUMTT7(5)-GIOT*CSUMTT7(6)))*ADIV2
!***********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM58(N4,ID,IK)=(CSUMTTMOM7(1,IK)+GIOT*CSUMTTMOM7(2,IK) &
  & -CSUMTTMOM7(3,IK)-GIOT*CSUMTTMOM7(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE59(N4,ID)=(CSUMTT7(1)+GIOT*CSUMTT7(2)-CSUMTT7(3) &
  & -GIOT*CSUMTT7(4)- &
  & (CSUMTT7(7)+GIOT*CSUMTT7(8)-CSUMTT7(5)-GIOT*CSUMTT7(6)))*ADIV2
!***********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM59(N4,ID,IK)=(CSUMTTMOM7(7,IK)+GIOT*CSUMTTMOM7(8,IK) &
  & -CSUMTTMOM7(5,IK)-GIOT*CSUMTTMOM7(6,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE60(N4,ID)=(CSUMTT7(1)-CSUMTT7(2)+CSUMTT7(3)-CSUMTT7(4)+ &
  & (CSUMTT7(8)-CSUMTT7(7)+CSUMTT7(6)-CSUMTT7(5)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM60(N4,ID,IK)=(CSUMTTMOM7(1,IK)-CSUMTTMOM7(2,IK) &
  & +CSUMTTMOM7(3,IK)-CSUMTTMOM7(4,IK)+(CSUMTTMOM7(8,IK) &
  & -CSUMTTMOM7(7,IK)+CSUMTTMOM7(6,IK) &
  & -CSUMTTMOM7(5,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE61(N4,ID)=(CSUMTT7(1)-CSUMTT7(2)+CSUMTT7(3)-CSUMTT7(4)- &
  & (CSUMTT7(8)-CSUMTT7(7)+CSUMTT7(6)-CSUMTT7(5)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM61(N4,ID,IK)=(CSUMTTMOM7(1,IK)-CSUMTTMOM7(2,IK) &
  & +CSUMTTMOM7(3,IK)-CSUMTTMOM7(4,IK)-(CSUMTTMOM7(8,IK) &
  & -CSUMTTMOM7(7,IK)+CSUMTTMOM7(6,IK) &
  & -CSUMTTMOM7(5,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 8
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-8 OPERATORS CP=+,Pr=+
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE62(N4,ID)=(CSUMTT8(1)+CSUMTT8(2)+CSUMTT8(3)+CSUMTT8(4) &
  & +(CSUMTT8(5)+CSUMTT8(6)+CSUMTT8(7)+CSUMTT8(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM62(N4,ID,IK)=(CSUMTTMOM8(1,IK)+CSUMTTMOM8(2,IK) &
  & +CSUMTTMOM8(3,IK)+CSUMTTMOM8(4,IK)+(CSUMTTMOM8(5,IK) &
  & +CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK) &
  & +CSUMTTMOM8(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE63(N4,ID)=(CSUMTT8(1)+CSUMTT8(2)+CSUMTT8(3)+CSUMTT8(4) &
  & -(CSUMTT8(5)+CSUMTT8(6)+CSUMTT8(7)+CSUMTT8(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM63(N4,ID,IK)=(CSUMTTMOM8(1,IK)+CSUMTTMOM8(2,IK) &
  & +CSUMTTMOM8(3,IK)+CSUMTTMOM8(4,IK)-(CSUMTTMOM8(5,IK) &
  & +CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK) &
  & +CSUMTTMOM8(8,IK)))*ADIV2
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE64(N4,ID)=(CSUMTT8(1)+GIOT*CSUMTT8(2)-CSUMTT8(3) &
  & -GIOT*CSUMTT8(4))*ADIV1
!***********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM64(N4,ID,IK)=(CSUMTTMOM8(1,IK)+GIOT*CSUMTTMOM8(2,IK) &
  & -CSUMTTMOM8(3,IK)-GIOT*CSUMTTMOM8(4,IK))*ADIV1
end do
!***********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE65(N4,ID)=(CSUMTT8(7)+GIOT*CSUMTT8(6)-CSUMTT8(5) &
  & -GIOT*CSUMTT8(8))*ADIV1
!***********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM65(N4,ID,IK)=(CSUMTTMOM8(7,IK)+GIOT*CSUMTTMOM8(6,IK) &
  & -CSUMTTMOM8(5,IK)-GIOT*CSUMTTMOM8(8,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE66(N4,ID)=(CSUMTT8(1)-CSUMTT8(2)+CSUMTT8(3)-CSUMTT8(4) &
  & +(CSUMTT8(5)-CSUMTT8(6)+CSUMTT8(7)-CSUMTT8(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM66(N4,ID,IK)=(CSUMTTMOM8(1,IK)-CSUMTTMOM8(2,IK) &
  & +CSUMTTMOM8(3,IK)-CSUMTTMOM8(4,IK)+(CSUMTTMOM8(5,IK) &
  & -CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK) &
  & -CSUMTTMOM8(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE67(N4,ID)=(CSUMTT8(1)-CSUMTT8(2)+CSUMTT8(3)-CSUMTT8(4) &
  & -(CSUMTT8(5)-CSUMTT8(6)+CSUMTT8(7)-CSUMTT8(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM67(N4,ID,IK)=(CSUMTTMOM8(1,IK)-CSUMTTMOM8(2,IK) &
  & +CSUMTTMOM8(3,IK)-CSUMTTMOM8(4,IK)-(CSUMTTMOM8(5,IK) &
  & -CSUMTTMOM8(6,IK)+CSUMTTMOM8(7,IK) &
  & -CSUMTTMOM8(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 9
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-9 OPERATORS 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE68(N4,ID)=(CSUMTT9(1)+CSUMTT9(2)+CSUMTT9(3)+CSUMTT9(4)+ &
  & (CSUMTT9(5)+CSUMTT9(6)+CSUMTT9(7)+CSUMTT9(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1,2
ALINEMOM68(N4,ID,IK)=(CSUMTTMOM9(1,IK)+CSUMTTMOM9(2,IK)+ &
  & CSUMTTMOM9(3,IK)+CSUMTTMOM9(4,IK) &
  & +(CSUMTTMOM9(5,IK)+CSUMTTMOM9(6,IK)+CSUMTTMOM9(7,IK) &
  & +CSUMTTMOM9(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-9 OPERATORS CP=-,Pz=-,J=0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE69(N4,ID)=(CSUMTT9(1)+CSUMTT9(2)+CSUMTT9(3)+CSUMTT9(4) &
  & -(CSUMTT9(5)+CSUMTT9(6)+CSUMTT9(7)+CSUMTT9(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1,2
ALINEMOM69(N4,ID,IK)=(CSUMTTMOM9(1,IK)+CSUMTTMOM9(2,IK)+ &
  & CSUMTTMOM9(3,IK)+CSUMTTMOM9(4,IK) &
  & -(CSUMTTMOM9(5,IK)+CSUMTTMOM9(6,IK)+CSUMTTMOM9(7,IK) &
  & +CSUMTTMOM9(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************      
ALINE70(N4,ID)=(CSUMTT9(1)+GIOT*CSUMTT9(2)-CSUMTT9(3) &
  & -GIOT*CSUMTT9(4)+(CSUMTT9(7)+GIOT*CSUMTT9(8)-CSUMTT9(5) &
  & -GIOT*CSUMTT9(6)))*ADIV2
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************      
DO IK=1,2
ALINEMOM70(N4,ID,IK)=(CSUMTTMOM9(1,IK)+GIOT*CSUMTTMOM9(2,IK) &
  & -CSUMTTMOM9(3,IK)-GIOT*CSUMTTMOM9(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************            
ALINE71(N4,ID)=(CSUMTT9(1)+GIOT*CSUMTT9(2)-CSUMTT9(3) &
  & -GIOT*CSUMTT9(4)-(CSUMTT9(7)+GIOT*CSUMTT9(8)-CSUMTT9(5) &
  & -GIOT*CSUMTT9(6)))*ADIV2
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************            
DO IK=1,2
ALINEMOM71(N4,ID,IK)=(CSUMTTMOM9(7,IK)+GIOT*CSUMTTMOM9(8,IK) &
  & -CSUMTTMOM9(5,IK)-GIOT*CSUMTTMOM9(6,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************            
ALINE72(N4,ID)=(CSUMTT9(1)-CSUMTT9(2)+CSUMTT9(3)-CSUMTT9(4) &
  & +(CSUMTT9(7)-CSUMTT9(6)+CSUMTT9(5)-CSUMTT9(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2
!**********************************************************************            
DO IK=1,2
ALINEMOM72(N4,ID,IK)=(CSUMTTMOM9(1,IK)-CSUMTTMOM9(2,IK) &
  & +CSUMTTMOM9(3,IK)-CSUMTTMOM9(4,IK) &
  & +(CSUMTTMOM9(7,IK)-CSUMTTMOM9(6,IK)+CSUMTTMOM9(5,IK) &
  & -CSUMTTMOM9(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=-, q=0 !Here!
!**********************************************************************                  
ALINE73(N4,ID)=(CSUMTT9(1)-CSUMTT9(2)+CSUMTT9(3)-CSUMTT9(4) &
  & -(CSUMTT9(7)-CSUMTT9(6)+CSUMTT9(5)-CSUMTT9(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2
!**********************************************************************                  
DO IK=1,2
ALINEMOM73(N4,ID,IK)=(CSUMTTMOM9(1,IK)-CSUMTTMOM9(2,IK) &
  & +CSUMTTMOM9(3,IK)-CSUMTTMOM9(4,IK) &
  & -(CSUMTTMOM9(7,IK)-CSUMTTMOM9(6,IK)+CSUMTTMOM9(5,IK) &
  & -CSUMTTMOM9(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-10 OPERATORS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************                  
ALINE74(N4,ID)=(CSUMTT10(1)+CSUMTT10(2)+CSUMTT10(3)+CSUMTT10(4) &
  & +(CSUMTT10(5)+CSUMTT10(6)+CSUMTT10(7)+CSUMTT10(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2
!**********************************************************************                  
DO IK=1,2
ALINEMOM74(N4,ID,IK)=(CSUMTTMOM10(1,IK)+CSUMTTMOM10(2,IK) &
  & +CSUMTTMOM10(3,IK)+CSUMTTMOM10(4,IK) &
  & +(CSUMTTMOM10(5,IK)+CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK) &
  & +CSUMTTMOM10(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************                        
ALINE75(N4,ID)=(CSUMTT10(1)+CSUMTT10(2)+CSUMTT10(3)+CSUMTT10(4) &
  & -(CSUMTT10(5)+CSUMTT10(6)+CSUMTT10(7)+CSUMTT10(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2
!**********************************************************************                        
DO IK=1,2
ALINEMOM75(N4,ID,IK)=(CSUMTTMOM10(1,IK)+CSUMTTMOM10(2,IK) &
  & +CSUMTTMOM10(3,IK)+CSUMTTMOM10(4,IK) &
  & -(CSUMTTMOM10(5,IK)+CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK) &
  & +CSUMTTMOM10(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************                  
ALINE76(N4,ID)=(CSUMTT10(1)+GIOT*CSUMTT10(2)-CSUMTT10(3) &
  & -GIOT*CSUMTT10(4)+(CSUMTT10(7)+GIOT*CSUMTT10(6)-CSUMTT10(5) &
  & -GIOT*CSUMTT10(8)))*ADIV2
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1,2
ALINEMOM76(N4,ID,IK)=(CSUMTTMOM10(1,IK)+GIOT*CSUMTTMOM10(2,IK) &
  & -CSUMTTMOM10(3,IK)-GIOT*CSUMTTMOM10(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE77(N4,ID)=(CSUMTT10(1)+GIOT*CSUMTT10(2)-CSUMTT10(3) &
  & -GIOT*CSUMTT10(4)-(CSUMTT10(7)+GIOT*CSUMTT10(6)-CSUMTT10(5) &
  & -GIOT*CSUMTT10(8)))*ADIV2
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************
DO IK=1,2
ALINEMOM77(N4,ID,IK)=(CSUMTTMOM10(7,IK)+GIOT*CSUMTTMOM10(6,IK) &
  & -CSUMTTMOM10(5,IK)-GIOT*CSUMTTMOM10(8,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE78(N4,ID)=(CSUMTT10(1)-CSUMTT10(2)+CSUMTT10(3)-CSUMTT10(4) &
  & +(CSUMTT10(5)-CSUMTT10(6)+CSUMTT10(7)-CSUMTT10(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2
!**********************************************************************
DO IK=1,2
ALINEMOM78(N4,ID,IK)=(CSUMTTMOM10(1,IK)-CSUMTTMOM10(2,IK) &
  & +CSUMTTMOM10(3,IK)-CSUMTTMOM10(4,IK) &
  & +(CSUMTTMOM10(5,IK)-CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK) &
  & -CSUMTTMOM10(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************      
ALINE79(N4,ID)=(CSUMTT10(1)-CSUMTT10(2)+CSUMTT10(3)-CSUMTT10(4) &
  & -(CSUMTT10(5)-CSUMTT10(6)+CSUMTT10(7)-CSUMTT10(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2
!**********************************************************************      
DO IK=1,2
ALINEMOM79(N4,ID,IK)=(CSUMTTMOM10(1,IK)-CSUMTTMOM10(2,IK) &
  & +CSUMTTMOM10(3,IK)-CSUMTTMOM10(4,IK) &
  & -(CSUMTTMOM10(5,IK)-CSUMTTMOM10(6,IK)+CSUMTTMOM10(7,IK) &
  & -CSUMTTMOM10(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-11 OPERATORS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************            
ALINE80(N4,ID)=(CSUMTT11(1)+CSUMTT11(2)+CSUMTT11(3)+CSUMTT11(4) &
  & +(CSUMTT11(5)+CSUMTT11(6)+CSUMTT11(7)+CSUMTT11(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2
!**********************************************************************            
DO IK=1,2
ALINEMOM80(N4,ID,IK)=(CSUMTTMOM11(1,IK)+CSUMTTMOM11(2,IK) &
  & +CSUMTTMOM11(3,IK)+CSUMTTMOM11(4,IK) &
  & +(CSUMTTMOM11(5,IK)+CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK) &
  & +CSUMTTMOM11(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************            
ALINE81(N4,ID)=(CSUMTT11(1)+CSUMTT11(2)+CSUMTT11(3)+CSUMTT11(4) &
  & -(CSUMTT11(5)+CSUMTT11(6)+CSUMTT11(7)+CSUMTT11(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2
!**********************************************************************            
DO IK=1,2
ALINEMOM81(N4,ID,IK)=(CSUMTTMOM11(1,IK)+CSUMTTMOM11(2,IK) &
  & +CSUMTTMOM11(3,IK)+CSUMTTMOM11(4,IK) &
  & -(CSUMTTMOM11(5,IK)+CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK) &
  & +CSUMTTMOM11(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE82(N4,ID)=(CSUMTT11(1)+GIOT*CSUMTT11(2)-CSUMTT11(3) &
  & -GIOT*CSUMTT11(4))*ADIV1
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1,2
ALINEMOM82(N4,ID,IK)=(CSUMTTMOM11(1,IK)+GIOT*CSUMTTMOM11(2,IK) &
  & -CSUMTTMOM11(3,IK)-GIOT*CSUMTTMOM11(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************                        
ALINE83(N4,ID)=(CSUMTT11(7)+GIOT*CSUMTT11(6)-CSUMTT11(5) &
  & -GIOT*CSUMTT11(8))*ADIV1
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1,2
ALINEMOM83(N4,ID,IK)=(CSUMTTMOM11(7,IK)+GIOT*CSUMTTMOM11(6,IK) &
  & -CSUMTTMOM11(5,IK)-GIOT*CSUMTTMOM11(8,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE84(N4,ID)=(CSUMTT11(1)-CSUMTT11(2)+CSUMTT11(3)-CSUMTT11(4) &
  & +(CSUMTT11(5)-CSUMTT11(6)+CSUMTT11(7)-CSUMTT11(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2
!**********************************************************************                  
DO IK=1,2
ALINEMOM84(N4,ID,IK)=(CSUMTTMOM11(1,IK)-CSUMTTMOM11(2,IK) &
  & +CSUMTTMOM11(3,IK)-CSUMTTMOM11(4,IK) &
  & +(CSUMTTMOM11(5,IK)-CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK) &
  & -CSUMTTMOM11(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************                        
ALINE85(N4,ID)=(CSUMTT11(1)-CSUMTT11(2)+CSUMTT11(3)-CSUMTT11(4) &
  & -(CSUMTT11(5)-CSUMTT11(6)+CSUMTT11(7)-CSUMTT11(8)))*ADIV2
!c**********************************************************************
!     J=2, Pp=-, q=1,2
!**********************************************************************                        
DO IK=1,2
ALINEMOM85(N4,ID,IK)=(CSUMTTMOM11(1,IK)-CSUMTTMOM11(2,IK) &
  & +CSUMTTMOM11(3,IK)-CSUMTTMOM11(4,IK) &
  & -(CSUMTTMOM11(5,IK)-CSUMTTMOM11(6,IK)+CSUMTTMOM11(7,IK) &
  & -CSUMTTMOM11(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     12 THIS!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-12 OPERATORS 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE86(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4) &
  & +(CSUMTT12(5)+CSUMTT12(6)+CSUMTT12(7)+CSUMTT12(8)) &
  & +(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12)) &
  & +(CSUMTT12(13)+CSUMTT12(14)+CSUMTT12(15)+CSUMTT12(16)))*ADIV3
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM86(N4,ID,IK)=(CSUMTTMOM12(1,IK)+CSUMTTMOM12(2,IK) &
  & +CSUMTTMOM12(3,IK)+CSUMTTMOM12(4,IK) &
  & +(CSUMTTMOM12(9,IK)+CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK) &
  & +CSUMTTMOM12(12,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE87(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4) &
  & +(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12)) &
  & -(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6)) &
  & -(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM87(N4,ID,IK)=(CSUMTTMOM12(7,IK)+CSUMTTMOM12(8,IK) &
  & +CSUMTTMOM12(5,IK)+CSUMTTMOM12(6,IK)+(CSUMTTMOM12(15,IK) &
  & +CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)+CSUMTTMOM12(14,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE88(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4) &
  & -(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12)) &
  & +(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6)) &
  & -(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM88(N4,ID,IK)=(CSUMTTMOM12(1,IK)+CSUMTTMOM12(2,IK) &
  & +CSUMTTMOM12(3,IK)+CSUMTTMOM12(4,IK)-(CSUMTTMOM12(9,IK) &
  & +CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)+CSUMTTMOM12(12,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE89(N4,ID)=(CSUMTT12(1)+CSUMTT12(2)+CSUMTT12(3)+CSUMTT12(4) &
  & -(CSUMTT12(9)+CSUMTT12(10)+CSUMTT12(11)+CSUMTT12(12)) &
  & -(CSUMTT12(7)+CSUMTT12(8)+CSUMTT12(5)+CSUMTT12(6)) &
  & +(CSUMTT12(15)+CSUMTT12(16)+CSUMTT12(13)+CSUMTT12(14)))*ADIV3
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM89(N4,ID,IK)=(CSUMTTMOM12(7,IK)+CSUMTTMOM12(8,IK) &
  & +CSUMTTMOM12(5,IK)+CSUMTTMOM12(6,IK)-(CSUMTTMOM12(15,IK) &
  & +CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)+CSUMTTMOM12(14,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE90(N4,ID)=(CSUMTT12(1)+GIOT*CSUMTT12(2)-CSUMTT12(3) &
  & -GIOT*CSUMTT12(4) &
  & +(CSUMTT12(7)+GIOT*CSUMTT12(8)-CSUMTT12(5)-GIOT*CSUMTT12(6)) &
  & )*ADIV2
!**********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM90(N4,ID,IK)=(CSUMTTMOM12(1,IK)+GIOT*CSUMTTMOM12(2,IK) &
  & -CSUMTTMOM12(3,IK)-GIOT*CSUMTTMOM12(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE91(N4,ID)=(CSUMTT12(9)+GIOT*CSUMTT12(10)-CSUMTT12(11) &
  & -GIOT*CSUMTT12(12) &
  & +(CSUMTT12(15)+GIOT*CSUMTT12(16)-CSUMTT12(13)-GIOT*CSUMTT12(14)) &
  & )*ADIV2
!**********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM91(N4,ID,IK)=(CSUMTTMOM12(9,IK)+GIOT*CSUMTTMOM12(10,IK) &
  & -CSUMTTMOM12(11,IK)-GIOT*CSUMTTMOM12(12,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE92(N4,ID)=(CSUMTT12(1)+GIOT*CSUMTT12(2)-CSUMTT12(3) &
  & -GIOT*CSUMTT12(4) &
  & -(CSUMTT12(7)+GIOT*CSUMTT12(8)-CSUMTT12(5)-GIOT*CSUMTT12(6)) &
  & )*ADIV2
!**********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM92(N4,ID,IK)=(CSUMTTMOM12(7,IK)+GIOT*CSUMTTMOM12(8,IK) &
  & -CSUMTTMOM12(5,IK)-GIOT*CSUMTTMOM12(6,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE93(N4,ID)=(CSUMTT12(9)+GIOT*CSUMTT12(10)-CSUMTT12(11) &
  & -GIOT*CSUMTT12(12) &
  & -(CSUMTT12(15)+GIOT*CSUMTT12(16)-CSUMTT12(13)-GIOT*CSUMTT12(14)) &
  & )*ADIV2
!**********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM93(N4,ID,IK)=(CSUMTTMOM12(15,IK) &
  & +GIOT*CSUMTTMOM12(16,IK) &
  & -CSUMTTMOM12(13,IK)-GIOT*CSUMTTMOM12(14,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE94(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4) &
  & +(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12)) &
  & +(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6)) &
  & +(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM94(N4,ID,IK)=(CSUMTTMOM12(1,IK)-CSUMTTMOM12(2,IK) &
  & +CSUMTTMOM12(3,IK)-CSUMTTMOM12(4,IK) &
  & +(CSUMTTMOM12(9,IK)-CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK) &
  & -CSUMTTMOM12(12,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE95(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4) &
  & +(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12)) &
  & -(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6)) &
  & -(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM95(N4,ID,IK)=(CSUMTTMOM12(7,IK)-CSUMTTMOM12(8,IK) &
  & +CSUMTTMOM12(5,IK)-CSUMTTMOM12(6,IK)+(CSUMTTMOM12(15,IK) &
  & -CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)-CSUMTTMOM12(14,IK)) &
  & )*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE96(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4) &
  & -(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12)) &
  & +(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6)) &
  & -(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM96(N4,ID,IK)=(CSUMTTMOM12(1,IK)-CSUMTTMOM12(2,IK) &
  & +CSUMTTMOM12(3,IK)-CSUMTTMOM12(4,IK)-(CSUMTTMOM12(9,IK) &
  & -CSUMTTMOM12(10,IK)+CSUMTTMOM12(11,IK)-CSUMTTMOM12(12,IK)) &
  & )*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE97(N4,ID)=(CSUMTT12(1)-CSUMTT12(2)+CSUMTT12(3)-CSUMTT12(4) &
  & -(CSUMTT12(9)-CSUMTT12(10)+CSUMTT12(11)-CSUMTT12(12)) &
  & -(CSUMTT12(7)-CSUMTT12(8)+CSUMTT12(5)-CSUMTT12(6)) &
  & +(CSUMTT12(15)-CSUMTT12(16)+CSUMTT12(13)-CSUMTT12(14)))*ADIV3
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM97(N4,ID,IK)=(CSUMTTMOM12(7,IK)-CSUMTTMOM12(8,IK) &
  & +CSUMTTMOM12(5,IK)-CSUMTTMOM12(6,IK)-(CSUMTTMOM12(15,IK) &
  & -CSUMTTMOM12(16,IK)+CSUMTTMOM12(13,IK)-CSUMTTMOM12(14,IK)) &
  & )*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     13
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     TT-13 OPERATORS CP=+,Pz=+,J=0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE98(N4,ID)=(CSUMTT13(1)+CSUMTT13(2)+CSUMTT13(3)+CSUMTT13(4)+ &
  & CSUMTT13(5)+CSUMTT13(6)+CSUMTT13(7)+CSUMTT13(8))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM98(N4,ID,IK)=(CSUMTTMOM13(1,IK)+CSUMTTMOM13(2,IK) &
  & +CSUMTTMOM13(3,IK)+CSUMTTMOM13(4,IK)+CSUMTTMOM13(5,IK) &
  & +CSUMTTMOM13(6,IK)+CSUMTTMOM13(7,IK) &
  & +CSUMTTMOM13(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE99(N4,ID)=(CSUMTT13(1)+CSUMTT13(2)+CSUMTT13(3)+CSUMTT13(4) &
  & -(CSUMTT13(5)+CSUMTT13(6)+CSUMTT13(7)+CSUMTT13(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM99(N4,ID,IK)=(CSUMTTMOM13(1,IK)+CSUMTTMOM13(2,IK) &
  & +CSUMTTMOM13(3,IK)+CSUMTTMOM13(4,IK)-(CSUMTTMOM13(5,IK) &
  & +CSUMTTMOM13(6,IK)+CSUMTTMOM13(7,IK) &
  & +CSUMTTMOM13(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE100(N4,ID)=(CSUMTT13(1)+GIOT*CSUMTT13(2)-CSUMTT13(3) &
  & -GIOT*CSUMTT13(4)+(CSUMTT13(7)+GIOT*CSUMTT13(8)-CSUMTT13(5) &
  & -GIOT*CSUMTT13(6)))*ADIV2
!**********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM100(N4,ID,IK)=(CSUMTTMOM13(1,IK)+GIOT*CSUMTTMOM13(2,IK) &
  & -CSUMTTMOM13(3,IK)-GIOT*CSUMTTMOM13(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE101(N4,ID)=(CSUMTT13(1)+GIOT*CSUMTT13(2)-CSUMTT13(3) &
  & -GIOT*CSUMTT13(4) &
  & -(CSUMTT13(7)+GIOT*CSUMTT13(8)-CSUMTT13(5) &
  & -GIOT*CSUMTT13(6)))*ADIV2
!**********************************************************************
!     J=1, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM101(N4,ID,IK)=(CSUMTTMOM13(7,IK)+GIOT*CSUMTTMOM13(8,IK) &
  & -CSUMTTMOM13(5,IK)-GIOT*CSUMTTMOM13(6,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE102(N4,ID)=(CSUMTT13(1)-CSUMTT13(2)+CSUMTT13(3)-CSUMTT13(4)+ &
  & CSUMTT13(6)-CSUMTT13(5)+CSUMTT13(8)-CSUMTT13(7))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM102(N4,ID,IK)=(CSUMTTMOM13(1,IK)-CSUMTTMOM13(2,IK) &
  & +CSUMTTMOM13(3,IK)-CSUMTTMOM13(4,IK)+CSUMTTMOM13(6,IK) &
  & -CSUMTTMOM13(5,IK)+CSUMTTMOM13(8,IK) &
  & -CSUMTTMOM13(7,IK))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE103(N4,ID)=(CSUMTT13(1)-CSUMTT13(2)+CSUMTT13(3)-CSUMTT13(4) &
  & -(CSUMTT13(6)-CSUMTT13(5)+CSUMTT13(8)-CSUMTT13(7)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM103(N4,ID,IK)=(CSUMTTMOM13(1,IK)-CSUMTTMOM13(2,IK) &
  & +CSUMTTMOM13(3,IK)-CSUMTTMOM13(4,IK)-(CSUMTTMOM13(6,IK) &
  & -CSUMTTMOM13(5,IK)+CSUMTTMOM13(8,IK) &
  & -CSUMTTMOM13(7,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     14 THISSSSSSSS!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE104(N4,ID)=(CSUMTT14(1)+CSUMTT14(2)+CSUMTT14(3)+CSUMTT14(4)+ &
  & CSUMTT14(5)+CSUMTT14(6)+CSUMTT14(7)+CSUMTT14(8))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=0
!**********************************************************************
DO IK=1,2
ALINEMOM104(N4,ID,IK)=(CSUMTTMOM14(1,IK)+CSUMTTMOM14(2,IK) &
  & +CSUMTTMOM14(3,IK)+CSUMTTMOM14(4,IK)+CSUMTTMOM14(5,IK) &
  & +CSUMTTMOM14(6,IK)+CSUMTTMOM14(7,IK) &
  & +CSUMTTMOM14(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE105(N4,ID)=(CSUMTT14(1)+CSUMTT14(2)+CSUMTT14(3)+CSUMTT14(4)- &
  & (CSUMTT14(5)+CSUMTT14(6)+CSUMTT14(7)+CSUMTT14(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM105(N4,ID,IK)=(CSUMTTMOM14(1,IK)+CSUMTTMOM14(2,IK) &
  & +CSUMTTMOM14(3,IK)+CSUMTTMOM14(4,IK)-(CSUMTTMOM14(5,IK) &
  & +CSUMTTMOM14(6,IK)+CSUMTTMOM14(7,IK) &
  & +CSUMTTMOM14(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************
ALINE106(N4,ID)=(CSUMTT14(1)+GIOT*CSUMTT14(2)-CSUMTT14(3) &
  & -GIOT*CSUMTT14(4)+(CSUMTT14(6)+GIOT*CSUMTT14(7)-CSUMTT14(8) &
  & -GIOT*CSUMTT14(5)))*ADIV2
!**********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM106(N4,ID,IK)=(CSUMTTMOM14(1,IK)+GIOT*CSUMTTMOM14(2,IK) &
  & -CSUMTTMOM14(3,IK)-GIOT*CSUMTTMOM14(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************
ALINE107(N4,ID)=(CSUMTT14(1)+GIOT*CSUMTT14(2)-CSUMTT14(3) &
  & -GIOT*CSUMTT14(4)-(CSUMTT14(6)+GIOT*CSUMTT14(7)-CSUMTT14(8) &
  & -GIOT*CSUMTT14(5)))*ADIV2
!**********************************************************************
!     J=1, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM107(N4,ID,IK)=(CSUMTTMOM14(6,IK)+GIOT*CSUMTTMOM14(7,IK) &
  & -CSUMTTMOM14(8,IK)-GIOT*CSUMTTMOM14(5,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE108(N4,ID)=(CSUMTT14(1)-CSUMTT14(2)+CSUMTT14(3)-CSUMTT14(4)+ &
  & (CSUMTT14(7)-CSUMTT14(6)+CSUMTT14(5)-CSUMTT14(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM108(N4,ID,IK)=(CSUMTTMOM14(1,IK)-CSUMTTMOM14(2,IK) &
  & +CSUMTTMOM14(3,IK)-CSUMTTMOM14(4,IK)+(CSUMTTMOM14(7,IK) &
  & -CSUMTTMOM14(6,IK)+CSUMTTMOM14(5,IK)-CSUMTTMOM14(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE109(N4,ID)=(CSUMTT14(1)-CSUMTT14(2)+CSUMTT14(3)-CSUMTT14(4)- &
  & (CSUMTT14(7)-CSUMTT14(6)+CSUMTT14(5)-CSUMTT14(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2,3,4
!**********************************************************************
DO IK=1,2
ALINEMOM109(N4,ID,IK)=(CSUMTTMOM14(1,IK)-CSUMTTMOM14(2,IK) &
  & +CSUMTTMOM14(3,IK)-CSUMTTMOM14(4,IK)-(CSUMTTMOM14(7,IK) &
  & -CSUMTTMOM14(6,IK)+CSUMTTMOM14(5,IK)-CSUMTTMOM14(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATOR 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE110(N4,ID)=(CSUMPLQ(1)+CSUMPLQ(2)+CSUMPLQ(3)+CSUMPLQ(4) &
  & +CSUMPLQ(5)+CSUMPLQ(6)+CSUMPLQ(7)+CSUMPLQ(8))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM110(N4,ID,IK)=(CSUMPLQMOM(1,IK)+CSUMPLQMOM(2,IK) &
  & +CSUMPLQMOM(3,IK)+CSUMPLQMOM(4,IK)+CSUMPLQMOM(5,IK) &
  & +CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)+CSUMPLQMOM(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************      
ALINE111(N4,ID)=(CSUMPLQ(1)+CSUMPLQ(2)+CSUMPLQ(3)+CSUMPLQ(4) &
  & -CSUMPLQ(5)-CSUMPLQ(6)-CSUMPLQ(7)-CSUMPLQ(8))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM111(N4,ID,IK)=(CSUMPLQMOM(1,IK)+CSUMPLQMOM(2,IK) &
  & +CSUMPLQMOM(3,IK)+CSUMPLQMOM(4,IK)-CSUMPLQMOM(5,IK) &
  & -CSUMPLQMOM(6,IK)-CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************      
ALINE112(N4,ID)=(CSUMPLQ(1)+GIOT*CSUMPLQ(2)-CSUMPLQ(3) &
  & -GIOT*CSUMPLQ(4) &
  & +(CSUMPLQ(6)+GIOT*CSUMPLQ(5)-CSUMPLQ(8)-GIOT*CSUMPLQ(7)))*ADIV2
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************      
DO IK=1, 2
ALINEMOM112(N4,ID,IK)=(CSUMPLQMOM(1,IK)+GIOT*CSUMPLQMOM(2,IK) &
  & -CSUMPLQMOM(3,IK)-GIOT*CSUMPLQMOM(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************      
ALINE113(N4,ID)=(CSUMPLQ(1)+GIOT*CSUMPLQ(2)-CSUMPLQ(3) &
  & -GIOT*CSUMPLQ(4) &
  & -(CSUMPLQ(6)+GIOT*CSUMPLQ(5)-CSUMPLQ(8)-GIOT*CSUMPLQ(7)))*ADIV2
!**********************************************************************
!     J=1, q=0
!**********************************************************************      
DO IK=1, 2
ALINEMOM113(N4,ID,IK)=(CSUMPLQMOM(6,IK)+GIOT*CSUMPLQMOM(5,IK) &
  & -CSUMPLQMOM(8,IK)-GIOT*CSUMPLQMOM(7,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************      
ALINE114(N4,ID)=(CSUMPLQ(1)-CSUMPLQ(2)+CSUMPLQ(3)-CSUMPLQ(4) &
  & +CSUMPLQ(5)-CSUMPLQ(6)+CSUMPLQ(7)-CSUMPLQ(8))*ADIV2
!**********************************************************************
!     J=2, Pr=+, q=0
!**********************************************************************      
DO IK=1, 2
ALINEMOM114(N4,ID,IK)=(CSUMPLQMOM(1,IK)-CSUMPLQMOM(2,IK) &
  & +CSUMPLQMOM(3,IK)-CSUMPLQMOM(4,IK)+CSUMPLQMOM(5,IK) &
  & -CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************            
ALINE115(N4,ID)=(CSUMPLQ(1)-CSUMPLQ(2)+CSUMPLQ(3)-CSUMPLQ(4) &
  & -(CSUMPLQ(5)-CSUMPLQ(6)+CSUMPLQ(7)-CSUMPLQ(8)))*ADIV2
!**********************************************************************
!     J=2, Pr=-, q=0
!**********************************************************************            
DO IK=1, 2
ALINEMOM115(N4,ID,IK)=(CSUMPLQMOM(1,IK)-CSUMPLQMOM(2,IK) &
  & +CSUMPLQMOM(3,IK)-CSUMPLQMOM(4,IK)-(CSUMPLQMOM(5,IK) &
  & -CSUMPLQMOM(6,IK)+CSUMPLQMOM(7,IK)-CSUMPLQMOM(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATOR 2
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************                  
ALINE116(N4,ID)=(CSUMPLQ2(1)+CSUMPLQ2(2)+CSUMPLQ2(3)+CSUMPLQ2(4) &
  & +CSUMPLQ2(5)+CSUMPLQ2(6)+CSUMPLQ2(7)+CSUMPLQ2(8))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM116(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+CSUMPLQMOM2(2,IK) &
  & +CSUMPLQMOM2(3,IK)+CSUMPLQMOM2(4,IK)+CSUMPLQMOM2(5,IK) &
  & +CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)+CSUMPLQMOM2(8,IK))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************                  
ALINE117(N4,ID)=(CSUMPLQ2(1)+CSUMPLQ2(2)+CSUMPLQ2(3)+CSUMPLQ2(4) &
  & -(CSUMPLQ2(5)+CSUMPLQ2(6)+CSUMPLQ2(7)+CSUMPLQ2(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM117(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+CSUMPLQMOM2(2,IK) &
  & +CSUMPLQMOM2(3,IK)+CSUMPLQMOM2(4,IK)-(CSUMPLQMOM2(5,IK) &
  & +CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)+CSUMPLQMOM2(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************                        
ALINE118(N4,ID)=(CSUMPLQ2(1)+GIOT*CSUMPLQ2(2)-CSUMPLQ2(3) &
  & -GIOT*CSUMPLQ2(4) &
  & +CSUMPLQ2(6)+GIOT*CSUMPLQ2(5)-CSUMPLQ2(8)-GIOT*CSUMPLQ2(7))*ADIV2
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM118(N4,ID,IK)=(CSUMPLQMOM2(1,IK)+GIOT*CSUMPLQMOM2(2,IK) &
  & -CSUMPLQMOM2(3,IK)-GIOT*CSUMPLQMOM2(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=-, q=0
!**********************************************************************                            
ALINE119(N4,ID)=(CSUMPLQ2(1)+GIOT*CSUMPLQ2(2)-CSUMPLQ2(3) &
  & -GIOT*CSUMPLQ2(4) &
  & -(CSUMPLQ2(6)+GIOT*CSUMPLQ2(5)-CSUMPLQ2(8)-GIOT*CSUMPLQ2(7))) &
  & *ADIV2
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM119(N4,ID,IK)=(CSUMPLQMOM2(6,IK)+GIOT*CSUMPLQMOM2(5,IK) &
  & -CSUMPLQMOM2(8,IK)-GIOT*CSUMPLQMOM2(7,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************                  
ALINE120(N4,ID)=(CSUMPLQ2(1)-CSUMPLQ2(2)+CSUMPLQ2(3)-CSUMPLQ2(4) &
  & +CSUMPLQ2(5)-CSUMPLQ2(6)+CSUMPLQ2(7)-CSUMPLQ2(8))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM120(N4,ID,IK)=(CSUMPLQMOM2(1,IK)-CSUMPLQMOM2(2,IK) &
  & +CSUMPLQMOM2(3,IK)-CSUMPLQMOM2(4,IK)+CSUMPLQMOM2(5,IK) &
  & -CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)-CSUMPLQMOM2(8,IK))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************                        
ALINE121(N4,ID)=(CSUMPLQ2(1)-CSUMPLQ2(2)+CSUMPLQ2(3)-CSUMPLQ2(4) &
  & -(CSUMPLQ2(5)-CSUMPLQ2(6)+CSUMPLQ2(7)-CSUMPLQ2(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM121(N4,ID,IK)=(CSUMPLQMOM2(1,IK)-CSUMPLQMOM2(2,IK) &
  & +CSUMPLQMOM2(3,IK)-CSUMPLQMOM2(4,IK)-(CSUMPLQMOM2(5,IK) &
  & -CSUMPLQMOM2(6,IK)+CSUMPLQMOM2(7,IK)-CSUMPLQMOM2(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PLAQUETTE OPERATORS 3
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************                        
ALINE122(N4,ID)=(CSUMPLQ3(1)+CSUMPLQ3(2)+CSUMPLQ3(3) &
  & +CSUMPLQ3(4)+(CSUMPLQ3(5)+CSUMPLQ3(6)+CSUMPLQ3(7) &
  & +CSUMPLQ3(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=+, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM122(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+CSUMPLQMOM3(2,IK) &
  & +CSUMPLQMOM3(3,IK)+CSUMPLQMOM3(4,IK) &
  & +(CSUMPLQMOM3(5,IK)+CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK) &
  & +CSUMPLQMOM3(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************                        
ALINE123(N4,ID)=(CSUMPLQ3(1)+CSUMPLQ3(2)+CSUMPLQ3(3) &
  & +CSUMPLQ3(4)-(CSUMPLQ3(5)+CSUMPLQ3(6)+CSUMPLQ3(7) &
  & +CSUMPLQ3(8)))*ADIV2
!**********************************************************************
!     J=0, Pp=-, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM123(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+CSUMPLQMOM3(2,IK) &
  & +CSUMPLQMOM3(3,IK)+CSUMPLQMOM3(4,IK) &
  & -(CSUMPLQMOM3(5,IK)+CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK) &
  & +CSUMPLQMOM3(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************                            
ALINE124(N4,ID)=(CSUMPLQ3(1)+GIOT*CSUMPLQ3(2)-CSUMPLQ3(3) &
  & -GIOT*CSUMPLQ3(4))*ADIV1
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************                            
DO IK=1, 2
ALINEMOM124(N4,ID,IK)=(CSUMPLQMOM3(1,IK)+GIOT*CSUMPLQMOM3(2,IK) &
  & -CSUMPLQMOM3(3,IK)-GIOT*CSUMPLQMOM3(4,IK))*ADIV1
end do
!**********************************************************************
!     J=1, Pr=+, q=0
!**********************************************************************                            
ALINE125(N4,ID)=(CSUMPLQ3(7)+GIOT*CSUMPLQ3(6)-CSUMPLQ3(5) &
  & -GIOT*CSUMPLQ3(8))*ADIV1
!**********************************************************************
!     J=1, q=1,2
!**********************************************************************                            
DO IK=1, 2
ALINEMOM125(N4,ID,IK)=(CSUMPLQMOM3(7,IK)+GIOT*CSUMPLQMOM3(6,IK) &
  & -CSUMPLQMOM3(5,IK)-GIOT*CSUMPLQMOM3(8,IK))*ADIV1
end do
!**********************************************************************
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************                            
ALINE126(N4,ID)=(CSUMPLQ3(1)-CSUMPLQ3(2)+CSUMPLQ3(3) &
  & -CSUMPLQ3(4)+(CSUMPLQ3(5)-CSUMPLQ3(6)+CSUMPLQ3(7) &
  & -CSUMPLQ3(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=+, q=1,2
!**********************************************************************                            
DO IK=1, 2
ALINEMOM126(N4,ID,IK)=(CSUMPLQMOM3(1,IK)-CSUMPLQMOM3(2,IK) &
  & +CSUMPLQMOM3(3,IK)-CSUMPLQMOM3(4,IK) &
  & +(CSUMPLQMOM3(5,IK)-CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK) &
  & -CSUMPLQMOM3(8,IK)))*ADIV2
end do
!**********************************************************************
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************                            
ALINE127(N4,ID)=(CSUMPLQ3(1)-CSUMPLQ3(2)+CSUMPLQ3(3) &
  & -CSUMPLQ3(4)-(CSUMPLQ3(5)-CSUMPLQ3(6)+CSUMPLQ3(7) &
  & -CSUMPLQ3(8)))*ADIV2
!**********************************************************************
!     J=2, Pp=-, q=1,2
!**********************************************************************                            
DO IK=1, 2
ALINEMOM127(N4,ID,IK)=(CSUMPLQMOM3(1,IK)-CSUMPLQMOM3(2,IK) &
  & +CSUMPLQMOM3(3,IK)-CSUMPLQMOM3(4,IK) &
  & -(CSUMPLQMOM3(5,IK)-CSUMPLQMOM3(6,IK)+CSUMPLQMOM3(7,IK) &
  & -CSUMPLQMOM3(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 4
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************            
ALINE128(N4,ID)=(CSUMPLQ4(1)+CSUMPLQ4(2)+CSUMPLQ4(3) &
  & +CSUMPLQ4(4)+(CSUMPLQ4(5)+CSUMPLQ4(6)+CSUMPLQ4(7) &
  & +CSUMPLQ4(8)))*ADIV2
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM128(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+CSUMPLQMOM4(2,IK) &
  & +CSUMPLQMOM4(3,IK)+CSUMPLQMOM4(4,IK) &
  & +(CSUMPLQMOM4(5,IK)+CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK) &
  & +CSUMPLQMOM4(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************                  
ALINE129(N4,ID)=(CSUMPLQ4(1)+CSUMPLQ4(2)+CSUMPLQ4(3) &
  & +CSUMPLQ4(4)-(CSUMPLQ4(5)+CSUMPLQ4(6)+CSUMPLQ4(7) &
  & +CSUMPLQ4(8)))*ADIV2
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM129(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+CSUMPLQMOM4(2,IK) &
  & +CSUMPLQMOM4(3,IK)+CSUMPLQMOM4(4,IK) &
  & -(CSUMPLQMOM4(5,IK)+CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK) &
  & +CSUMPLQMOM4(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************                            
ALINE130(N4,ID)=(CSUMPLQ4(1)+GIOT*CSUMPLQ4(2)-CSUMPLQ4(3) &
  & -GIOT*CSUMPLQ4(4))*ADIV1
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                            
DO IK=1, 2
ALINEMOM130(N4,ID,IK)=(CSUMPLQMOM4(1,IK)+GIOT*CSUMPLQMOM4(2,IK) &
  & -CSUMPLQMOM4(3,IK)-GIOT*CSUMPLQMOM4(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************                            
ALINE131(N4,ID)=(CSUMPLQ4(7)+GIOT*CSUMPLQ4(6)-CSUMPLQ4(5) &
  & -GIOT*CSUMPLQ4(8))*ADIV1
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                            
DO IK=1, 2
ALINEMOM131(N4,ID,IK)=(CSUMPLQMOM4(7,IK)+GIOT*CSUMPLQMOM4(6,IK) &
  & -CSUMPLQMOM4(5,IK)-GIOT*CSUMPLQMOM4(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************            
ALINE132(N4,ID)=(CSUMPLQ4(1)-CSUMPLQ4(2)+CSUMPLQ4(3) &
  & -CSUMPLQ4(4)+(CSUMPLQ4(5)-CSUMPLQ4(6)+CSUMPLQ4(7) &
  & -CSUMPLQ4(8)))*ADIV2
!**********************************************************************      
!     J=2, Pp=+, q=1,2
!**********************************************************************            
DO IK=1, 2
ALINEMOM132(N4,ID,IK)=(CSUMPLQMOM4(1,IK)-CSUMPLQMOM4(2,IK) &
  & +CSUMPLQMOM4(3,IK)-CSUMPLQMOM4(4,IK) &
  & +(CSUMPLQMOM4(5,IK)-CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK) &
  & -CSUMPLQMOM4(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************                  
ALINE133(N4,ID)=(CSUMPLQ4(1)-CSUMPLQ4(2)+CSUMPLQ4(3) &
  & -CSUMPLQ4(4)-(CSUMPLQ4(5)-CSUMPLQ4(6)+CSUMPLQ4(7) &
  & -CSUMPLQ4(8)))*ADIV2
!**********************************************************************      
!     J=2, Pp=-, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM133(N4,ID,IK)=(CSUMPLQMOM4(1,IK)-CSUMPLQMOM4(2,IK) &
  & +CSUMPLQMOM4(3,IK)-CSUMPLQMOM4(4,IK) &
  & -(CSUMPLQMOM4(5,IK)-CSUMPLQMOM4(6,IK)+CSUMPLQMOM4(7,IK) &
  & -CSUMPLQMOM4(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PLAQUETTE OPERATOR 5  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************                        
ALINE134(N4,ID)=(CSUMPLQ5(1)+CSUMPLQ5(2)+CSUMPLQ5(3) &
  & +CSUMPLQ5(4)+(CSUMPLQ5(5)+CSUMPLQ5(6)+CSUMPLQ5(7) &
  & +CSUMPLQ5(8)))*ADIV2
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************                            
DO IK=1, 2
ALINEMOM134(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+CSUMPLQMOM5(2,IK) &
  & +CSUMPLQMOM5(3,IK)+CSUMPLQMOM5(4,IK) &
  & +(CSUMPLQMOM5(5,IK)+CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK) &
  & +CSUMPLQMOM5(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************                            
ALINE135(N4,ID)=(CSUMPLQ5(1)+CSUMPLQ5(2)+CSUMPLQ5(3) &
  & +CSUMPLQ5(4)-(CSUMPLQ5(5)+CSUMPLQ5(6)+CSUMPLQ5(7) &
  & +CSUMPLQ5(8)))*ADIV2
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************      
DO IK=1, 2
ALINEMOM135(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+CSUMPLQMOM5(2,IK) &
  & +CSUMPLQMOM5(3,IK)+CSUMPLQMOM5(4,IK) &
  & -(CSUMPLQMOM5(5,IK)+CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK) &
  & +CSUMPLQMOM5(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************            
ALINE136(N4,ID)=(CSUMPLQ5(1)+GIOT*CSUMPLQ5(2)-CSUMPLQ5(3) &
  & -GIOT*CSUMPLQ5(4))*ADIV1
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM136(N4,ID,IK)=(CSUMPLQMOM5(1,IK)+GIOT*CSUMPLQMOM5(2,IK) &
  & -CSUMPLQMOM5(3,IK)-GIOT*CSUMPLQMOM5(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE137(N4,ID)=(CSUMPLQ5(7)+GIOT*CSUMPLQ5(6)-CSUMPLQ5(5) &
  & -GIOT*CSUMPLQ5(8))*ADIV1
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM137(N4,ID,IK)=(CSUMPLQMOM5(7,IK)+GIOT*CSUMPLQMOM5(6,IK) &
  & -CSUMPLQMOM5(5,IK)-GIOT*CSUMPLQMOM5(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************                        
ALINE138(N4,ID)=(CSUMPLQ5(1)-CSUMPLQ5(2)+CSUMPLQ5(3) &
  & -CSUMPLQ5(4)+(CSUMPLQ5(5)-CSUMPLQ5(6)+CSUMPLQ5(7) &
  & -CSUMPLQ5(8)))*ADIV2
!**********************************************************************      
!     J=2, Pp=+, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM138(N4,ID,IK)=(CSUMPLQMOM5(1,IK)-CSUMPLQMOM5(2,IK) &
  & +CSUMPLQMOM5(3,IK)-CSUMPLQMOM5(4,IK) &
  & +(CSUMPLQMOM5(5,IK)-CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK) &
  & -CSUMPLQMOM5(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************                            
ALINE139(N4,ID)=(CSUMPLQ5(1)-CSUMPLQ5(2)+CSUMPLQ5(3) &
  & -CSUMPLQ5(4)-(CSUMPLQ5(5)-CSUMPLQ5(6)+CSUMPLQ5(7) &
  & -CSUMPLQ5(8)))*ADIV2
!**********************************************************************      
!     J=2, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM139(N4,ID,IK)=(CSUMPLQMOM5(1,IK)-CSUMPLQMOM5(2,IK) &
  & +CSUMPLQMOM5(3,IK)-CSUMPLQMOM5(4,IK) &
  & -(CSUMPLQMOM5(5,IK)-CSUMPLQMOM5(6,IK)+CSUMPLQMOM5(7,IK) &
  & -CSUMPLQMOM5(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 6
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE140(N4,ID)=(CSUMPLQ6(1)+CSUMPLQ6(2)+CSUMPLQ6(3) &
  & +CSUMPLQ6(4)+(CSUMPLQ6(5)+CSUMPLQ6(6)+CSUMPLQ6(7) &
  & +CSUMPLQ6(8)))*ADIV2
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM140(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+CSUMPLQMOM6(2,IK) &
  & +CSUMPLQMOM6(3,IK)+CSUMPLQMOM6(4,IK) &
  & +(CSUMPLQMOM6(5,IK)+CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK) &
  & +CSUMPLQMOM6(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************      
ALINE141(N4,ID)=(CSUMPLQ6(1)+CSUMPLQ6(2)+CSUMPLQ6(3) &
  & +CSUMPLQ6(4)-(CSUMPLQ6(5)+CSUMPLQ6(6)+CSUMPLQ6(7) &
  & +CSUMPLQ6(8)))*ADIV2
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************            
DO IK=1, 2
ALINEMOM141(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+CSUMPLQMOM6(2,IK) &
  & +CSUMPLQMOM6(3,IK)+CSUMPLQMOM6(4,IK) &
  & -(CSUMPLQMOM6(5,IK)+CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK) &
  & +CSUMPLQMOM6(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE142(N4,ID)=(CSUMPLQ6(1)+GIOT*CSUMPLQ6(2)-CSUMPLQ6(3) &
  & -GIOT*CSUMPLQ6(4)+(CSUMPLQ6(8)+GIOT*CSUMPLQ6(7)-CSUMPLQ6(6) &
  & -GIOT*CSUMPLQ6(5)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM142(N4,ID,IK)=(CSUMPLQMOM6(1,IK)+GIOT*CSUMPLQMOM6(2,IK) &
  & -CSUMPLQMOM6(3,IK)-GIOT*CSUMPLQMOM6(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE143(N4,ID)=(CSUMPLQ6(1)+GIOT*CSUMPLQ6(2)-CSUMPLQ6(3) &
  & -GIOT*CSUMPLQ6(4)-(CSUMPLQ6(8)+GIOT*CSUMPLQ6(7)-CSUMPLQ6(6) &
  & -GIOT*CSUMPLQ6(5)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                        
DO IK=1, 2
ALINEMOM143(N4,ID,IK)=(CSUMPLQMOM6(8,IK)+GIOT*CSUMPLQMOM6(7,IK) &
  & -CSUMPLQMOM6(6,IK)-GIOT*CSUMPLQMOM6(5,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************                        
ALINE144(N4,ID)=(CSUMPLQ6(1)-CSUMPLQ6(2)+CSUMPLQ6(3) &
  & -CSUMPLQ6(4)+(CSUMPLQ6(5)-CSUMPLQ6(6)+CSUMPLQ6(7) &
  & -CSUMPLQ6(8)))*ADIV2
!**********************************************************************      
!     J=2, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM144(N4,ID,IK)=(CSUMPLQMOM6(1,IK)-CSUMPLQMOM6(2,IK) &
  & +CSUMPLQMOM6(3,IK)-CSUMPLQMOM6(4,IK) &
  & +(CSUMPLQMOM6(5,IK)-CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK) &
  & -CSUMPLQMOM6(8,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE145(N4,ID)=(CSUMPLQ6(1)-CSUMPLQ6(2)+CSUMPLQ6(3) &
  & -CSUMPLQ6(4)-(CSUMPLQ6(5)-CSUMPLQ6(6)+CSUMPLQ6(7) &
  & -CSUMPLQ6(8)))*ADIV2
!**********************************************************************      
!     J=2, Pp=-, q=1,2
!**********************************************************************      
DO IK=1, 2
ALINEMOM145(N4,ID,IK)=(CSUMPLQMOM6(1,IK)-CSUMPLQMOM6(2,IK) &
  & +CSUMPLQMOM6(3,IK)-CSUMPLQMOM6(4,IK) &
  & -(CSUMPLQMOM6(5,IK)-CSUMPLQMOM6(6,IK)+CSUMPLQMOM6(7,IK) &
  & -CSUMPLQMOM6(8,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 7
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE146(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3) &
  & +CSUMPLQ7(4) &
  & +(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8)) &
  & +(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12)) &
  & +(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM146(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+CSUMPLQMOM7(2,IK) &
  & +CSUMPLQMOM7(3,IK)+CSUMPLQMOM7(4,IK) &
  & +(CSUMPLQMOM7(9,IK)+CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK) &
  & +CSUMPLQMOM7(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE147(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3) &
  & +CSUMPLQ7(4) &
  & -(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8)) &
  & +(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12)) &
  & -(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM147(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+CSUMPLQMOM7(6,IK) &
  & +CSUMPLQMOM7(7,IK)+CSUMPLQMOM7(8,IK) &
  & +(CSUMPLQMOM7(13,IK)+CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK) &
  & +CSUMPLQMOM7(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE148(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3) &
  & +CSUMPLQ7(4) &
  & +(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8)) &
  & -(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12)) &
  & -(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM148(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+CSUMPLQMOM7(2,IK) &
  & +CSUMPLQMOM7(3,IK)+CSUMPLQMOM7(4,IK) &
  & -(CSUMPLQMOM7(9,IK)+CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK) &
  & +CSUMPLQMOM7(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE149(N4,ID)=(CSUMPLQ7(1)+CSUMPLQ7(2)+CSUMPLQ7(3) &
  & +CSUMPLQ7(4) &
  & -(CSUMPLQ7(5)+CSUMPLQ7(6)+CSUMPLQ7(7)+CSUMPLQ7(8)) &
  & -(CSUMPLQ7(9)+CSUMPLQ7(10)+CSUMPLQ7(11)+CSUMPLQ7(12)) &
  & +(CSUMPLQ7(13)+CSUMPLQ7(14)+CSUMPLQ7(15)+CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM149(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+CSUMPLQMOM7(6,IK) &
  & +CSUMPLQMOM7(7,IK)+CSUMPLQMOM7(8,IK) &
  & -(CSUMPLQMOM7(13,IK)+CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK) &
  & +CSUMPLQMOM7(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE150(N4,ID)=(CSUMPLQ7(1)+GIOT*CSUMPLQ7(2)-CSUMPLQ7(3) &
  & -GIOT*CSUMPLQ7(4)+(CSUMPLQ7(5)+GIOT*CSUMPLQ7(6)-CSUMPLQ7(7) &
  & -GIOT*CSUMPLQ7(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM150(N4,ID,IK)=(CSUMPLQMOM7(1,IK)+GIOT*CSUMPLQMOM7(2,IK) &
  & -CSUMPLQMOM7(3,IK)-GIOT*CSUMPLQMOM7(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE151(N4,ID)=(CSUMPLQ7(1)+GIOT*CSUMPLQ7(2)-CSUMPLQ7(3) &
  & -GIOT*CSUMPLQ7(4)-(CSUMPLQ7(5)+GIOT*CSUMPLQ7(6)-CSUMPLQ7(7) &
  & -GIOT*CSUMPLQ7(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM151(N4,ID,IK)=(CSUMPLQMOM7(5,IK)+GIOT*CSUMPLQMOM7(6,IK) &
  & -CSUMPLQMOM7(7,IK)-GIOT*CSUMPLQMOM7(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE152(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3) &
  & -CSUMPLQ7(4) &
  & +(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8)) &
  & +(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12)) &
  & +(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q
!**********************************************************************
DO IK=1, 2
ALINEMOM152(N4,ID,IK)=(CSUMPLQMOM7(1,IK)-CSUMPLQMOM7(2,IK) &
  & +CSUMPLQMOM7(3,IK)-CSUMPLQMOM7(4,IK) &
  & +(CSUMPLQMOM7(9,IK)-CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK) &
  & -CSUMPLQMOM7(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE153(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3) &
  & -CSUMPLQ7(4) &
  & -(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8)) &
  & +(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12)) &
  & -(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM153(N4,ID,IK)=(CSUMPLQMOM7(5,IK)-CSUMPLQMOM7(6,IK) &
  & +CSUMPLQMOM7(7,IK)-CSUMPLQMOM7(8,IK) &
  & +(CSUMPLQMOM7(13,IK)-CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK) &
  & -CSUMPLQMOM7(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE154(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3) &
  & -CSUMPLQ7(4) &
  & +(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8)) &
  & -(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12)) &
  & -(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM154(N4,ID,IK)=(CSUMPLQMOM7(1,IK)-CSUMPLQMOM7(2,IK) &
  & +CSUMPLQMOM7(3,IK)-CSUMPLQMOM7(4,IK) &
  & -(CSUMPLQMOM7(9,IK)-CSUMPLQMOM7(10,IK)+CSUMPLQMOM7(11,IK) &
  & -CSUMPLQMOM7(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE155(N4,ID)=(CSUMPLQ7(1)-CSUMPLQ7(2)+CSUMPLQ7(3) &
  & -CSUMPLQ7(4) &
  & -(CSUMPLQ7(5)-CSUMPLQ7(6)+CSUMPLQ7(7)-CSUMPLQ7(8)) &
  & -(CSUMPLQ7(9)-CSUMPLQ7(10)+CSUMPLQ7(11)-CSUMPLQ7(12)) &
  & +(CSUMPLQ7(13)-CSUMPLQ7(14)+CSUMPLQ7(15)-CSUMPLQ7(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM155(N4,ID,IK)=(CSUMPLQMOM7(5,IK)-CSUMPLQMOM7(6,IK) &
  & +CSUMPLQMOM7(7,IK)-CSUMPLQMOM7(8,IK) &
  & -(CSUMPLQMOM7(13,IK)-CSUMPLQMOM7(14,IK)+CSUMPLQMOM7(15,IK) &
  & -CSUMPLQMOM7(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 8
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE156(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3) &
  & +CSUMPLQ8(4) &
  & +(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8)) &
  & +(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12)) &
  & +(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM156(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+CSUMPLQMOM8(2,IK) &
  & +CSUMPLQMOM8(3,IK)+CSUMPLQMOM8(4,IK) &
  & +(CSUMPLQMOM8(9,IK)+CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK) &
  & +CSUMPLQMOM8(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE157(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3) &
  & +CSUMPLQ8(4) &
  & -(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8)) &
  & +(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12)) &
  & -(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM157(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+CSUMPLQMOM8(6,IK) &
  & +CSUMPLQMOM8(7,IK)+CSUMPLQMOM8(8,IK) &
  & +(CSUMPLQMOM8(13,IK)+CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK) &
  & +CSUMPLQMOM8(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE158(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3) &
  & +CSUMPLQ8(4) &
  & +(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8)) &
  & -(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12)) &
  & -(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM158(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+CSUMPLQMOM8(2,IK) &
  & +CSUMPLQMOM8(3,IK)+CSUMPLQMOM8(4,IK) &
  & -(CSUMPLQMOM8(9,IK)+CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK) &
  & +CSUMPLQMOM8(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE159(N4,ID)=(CSUMPLQ8(1)+CSUMPLQ8(2)+CSUMPLQ8(3) &
  & +CSUMPLQ8(4) &
  & -(CSUMPLQ8(5)+CSUMPLQ8(6)+CSUMPLQ8(7)+CSUMPLQ8(8)) &
  & -(CSUMPLQ8(9)+CSUMPLQ8(10)+CSUMPLQ8(11)+CSUMPLQ8(12)) &
  & +(CSUMPLQ8(13)+CSUMPLQ8(14)+CSUMPLQ8(15)+CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM159(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+CSUMPLQMOM8(6,IK) &
  & +CSUMPLQMOM8(7,IK)+CSUMPLQMOM8(8,IK) &
  & -(CSUMPLQMOM8(13,IK)+CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK) &
  & +CSUMPLQMOM8(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE160(N4,ID)=(CSUMPLQ8(1)+GIOT*CSUMPLQ8(2)-CSUMPLQ8(3) &
  & -GIOT*CSUMPLQ8(4)+(CSUMPLQ8(5)+GIOT*CSUMPLQ8(6)-CSUMPLQ8(7) &
  & -GIOT*CSUMPLQ8(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM160(N4,ID,IK)=(CSUMPLQMOM8(1,IK)+GIOT*CSUMPLQMOM8(2,IK) &
  & -CSUMPLQMOM8(3,IK)-GIOT*CSUMPLQMOM8(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE161(N4,ID)=(CSUMPLQ8(1)+GIOT*CSUMPLQ8(2)-CSUMPLQ8(3) &
  & -GIOT*CSUMPLQ8(4)-(CSUMPLQ8(5)+GIOT*CSUMPLQ8(6)-CSUMPLQ8(7) &
  & -GIOT*CSUMPLQ8(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM161(N4,ID,IK)=(CSUMPLQMOM8(5,IK)+GIOT*CSUMPLQMOM8(6,IK) &
  & -CSUMPLQMOM8(7,IK)-GIOT*CSUMPLQMOM8(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE162(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3) &
  & -CSUMPLQ8(4) &
  & +(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8)) &
  & +(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12)) &
  & +(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM162(N4,ID,IK)=(CSUMPLQMOM8(1,IK)-CSUMPLQMOM8(2,IK) &
  & +CSUMPLQMOM8(3,IK)-CSUMPLQMOM8(4,IK) &
  & +(CSUMPLQMOM8(9,IK)-CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK) &
  & -CSUMPLQMOM8(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE163(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3) &
  & -CSUMPLQ8(4) &
  & -(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8)) &
  & +(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12)) &
  & -(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM163(N4,ID,IK)=(CSUMPLQMOM8(5,IK)-CSUMPLQMOM8(6,IK) &
  & +CSUMPLQMOM8(7,IK)-CSUMPLQMOM8(8,IK) &
  & +(CSUMPLQMOM8(13,IK)-CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK) &
  & -CSUMPLQMOM8(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE164(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3) &
  & -CSUMPLQ8(4) &
  & +(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8)) &
  & -(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12)) &
  & -(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM164(N4,ID,IK)=(CSUMPLQMOM8(1,IK)-CSUMPLQMOM8(2,IK) &
  & +CSUMPLQMOM8(3,IK)-CSUMPLQMOM8(4,IK) &
  & -(CSUMPLQMOM8(9,IK)-CSUMPLQMOM8(10,IK)+CSUMPLQMOM8(11,IK) &
  & -CSUMPLQMOM8(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE165(N4,ID)=(CSUMPLQ8(1)-CSUMPLQ8(2)+CSUMPLQ8(3) &
  & -CSUMPLQ8(4) &
  & -(CSUMPLQ8(5)-CSUMPLQ8(6)+CSUMPLQ8(7)-CSUMPLQ8(8)) &
  & -(CSUMPLQ8(9)-CSUMPLQ8(10)+CSUMPLQ8(11)-CSUMPLQ8(12)) &
  & +(CSUMPLQ8(13)-CSUMPLQ8(14)+CSUMPLQ8(15)-CSUMPLQ8(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM165(N4,ID,IK)=(CSUMPLQMOM8(5,IK)-CSUMPLQMOM8(6,IK) &
  & +CSUMPLQMOM8(7,IK)-CSUMPLQMOM8(8,IK) &
  & -(CSUMPLQMOM8(13,IK)-CSUMPLQMOM8(14,IK)+CSUMPLQMOM8(15,IK) &
  & -CSUMPLQMOM8(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 9
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE166(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3) &
  & +CSUMPLQ9(4) &
  & +(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8)) &
  & +(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12)) &
  & +(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM166(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+CSUMPLQMOM9(2,IK) &
  & +CSUMPLQMOM9(3,IK)+CSUMPLQMOM9(4,IK) &
  & +(CSUMPLQMOM9(9,IK)+CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK) &
  & +CSUMPLQMOM9(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE167(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3) &
  & +CSUMPLQ9(4) &
  & -(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8)) &
  & +(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12)) &
  & -(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM167(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+CSUMPLQMOM9(6,IK) &
  & +CSUMPLQMOM9(7,IK)+CSUMPLQMOM9(8,IK) &
  & +(CSUMPLQMOM9(13,IK)+CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK) &
  & +CSUMPLQMOM9(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE168(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3) &
  & +CSUMPLQ9(4) &
  & +(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8)) &
  & -(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12)) &
  & -(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM168(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+CSUMPLQMOM9(2,IK) &
  & +CSUMPLQMOM9(3,IK)+CSUMPLQMOM9(4,IK) &
  & -(CSUMPLQMOM9(9,IK)+CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK) &
  & +CSUMPLQMOM9(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE169(N4,ID)=(CSUMPLQ9(1)+CSUMPLQ9(2)+CSUMPLQ9(3) &
  & +CSUMPLQ9(4) &
  & -(CSUMPLQ9(5)+CSUMPLQ9(6)+CSUMPLQ9(7)+CSUMPLQ9(8)) &
  & -(CSUMPLQ9(9)+CSUMPLQ9(10)+CSUMPLQ9(11)+CSUMPLQ9(12)) &
  & +(CSUMPLQ9(13)+CSUMPLQ9(14)+CSUMPLQ9(15)+CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM169(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+CSUMPLQMOM9(6,IK) &
  & +CSUMPLQMOM9(7,IK)+CSUMPLQMOM9(8,IK) &
  & -(CSUMPLQMOM9(13,IK)+CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK) &
  & +CSUMPLQMOM9(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE170(N4,ID)=(CSUMPLQ9(1)+GIOT*CSUMPLQ9(2)-CSUMPLQ9(3) &
  & -GIOT*CSUMPLQ9(4)+(CSUMPLQ9(5)+GIOT*CSUMPLQ9(6)-CSUMPLQ9(7) &
  & -GIOT*CSUMPLQ9(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM170(N4,ID,IK)=(CSUMPLQMOM9(1,IK)+GIOT*CSUMPLQMOM9(2,IK) &
  & -CSUMPLQMOM9(3,IK)-GIOT*CSUMPLQMOM9(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE171(N4,ID)=(CSUMPLQ9(1)+GIOT*CSUMPLQ9(2)-CSUMPLQ9(3) &
  & -GIOT*CSUMPLQ9(4)-(CSUMPLQ9(5)+GIOT*CSUMPLQ9(6)-CSUMPLQ9(7) &
  & -GIOT*CSUMPLQ9(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM171(N4,ID,IK)=(CSUMPLQMOM9(5,IK)+GIOT*CSUMPLQMOM9(6,IK) &
  & -CSUMPLQMOM9(7,IK)-GIOT*CSUMPLQMOM9(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE172(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3) &
  & -CSUMPLQ9(4) &
  & +(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8)) &
  & +(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12)) &
  & +(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM172(N4,ID,IK)=(CSUMPLQMOM9(1,IK)-CSUMPLQMOM9(2,IK) &
  & +CSUMPLQMOM9(3,IK)-CSUMPLQMOM9(4,IK) &
  & +(CSUMPLQMOM9(9,IK)-CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK) &
  & -CSUMPLQMOM9(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE173(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3) &
  & -CSUMPLQ9(4) &
  & -(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8)) &
  & +(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12)) &
  & -(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM173(N4,ID,IK)=(CSUMPLQMOM9(5,IK)-CSUMPLQMOM9(6,IK) &
  & +CSUMPLQMOM9(7,IK)-CSUMPLQMOM9(8,IK) &
  & +(CSUMPLQMOM9(13,IK)-CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK) &
  & -CSUMPLQMOM9(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE174(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3) &
  & -CSUMPLQ9(4) &
  & +(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8)) &
  & -(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12)) &
  & -(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM174(N4,ID,IK)=(CSUMPLQMOM9(1,IK)-CSUMPLQMOM9(2,IK) &
  & +CSUMPLQMOM9(3,IK)-CSUMPLQMOM9(4,IK) &
  & -(CSUMPLQMOM9(9,IK)-CSUMPLQMOM9(10,IK)+CSUMPLQMOM9(11,IK) &
  & -CSUMPLQMOM9(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE175(N4,ID)=(CSUMPLQ9(1)-CSUMPLQ9(2)+CSUMPLQ9(3) &
  & -CSUMPLQ9(4) &
  & -(CSUMPLQ9(5)-CSUMPLQ9(6)+CSUMPLQ9(7)-CSUMPLQ9(8)) &
  & -(CSUMPLQ9(9)-CSUMPLQ9(10)+CSUMPLQ9(11)-CSUMPLQ9(12)) &
  & +(CSUMPLQ9(13)-CSUMPLQ9(14)+CSUMPLQ9(15)-CSUMPLQ9(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM175(N4,ID,IK)=(CSUMPLQMOM9(5,IK)-CSUMPLQMOM9(6,IK) &
  & +CSUMPLQMOM9(7,IK)-CSUMPLQMOM9(8,IK) &
  & -(CSUMPLQMOM9(13,IK)-CSUMPLQMOM9(14,IK)+CSUMPLQMOM9(15,IK) &
  & -CSUMPLQMOM9(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 10
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE176(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3) &
  & +CSUMPLQ10(4) &
  & +(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8)) &
  & +(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12)) &
  & +(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM176(N4,ID,IK)=(CSUMPLQMOM10(1,IK)+CSUMPLQMOM10(2,IK) &
  & +CSUMPLQMOM10(3,IK)+CSUMPLQMOM10(4,IK) &
  & +(CSUMPLQMOM10(9,IK)+CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK) &
  & +CSUMPLQMOM10(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE177(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3) &
  & +CSUMPLQ10(4) &
  & -(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8)) &
  & +(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12)) &
  & -(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM177(N4,ID,IK)=(CSUMPLQMOM10(5,IK)+CSUMPLQMOM10(6,IK) &
  & +CSUMPLQMOM10(7,IK)+CSUMPLQMOM10(8,IK) &
  & +(CSUMPLQMOM10(13,IK)+CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK) &
  & +CSUMPLQMOM10(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE178(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3) &
  & +CSUMPLQ10(4) &
  & +(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8)) &
  & -(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12)) &
  & -(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM178(N4,ID,IK)=(CSUMPLQMOM10(1,IK)+CSUMPLQMOM10(2,IK) &
  & +CSUMPLQMOM10(3,IK)+CSUMPLQMOM10(4,IK) &
  & -(CSUMPLQMOM10(9,IK)+CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK) &
  & +CSUMPLQMOM10(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE179(N4,ID)=(CSUMPLQ10(1)+CSUMPLQ10(2)+CSUMPLQ10(3) &
  & +CSUMPLQ10(4) &
  & -(CSUMPLQ10(5)+CSUMPLQ10(6)+CSUMPLQ10(7)+CSUMPLQ10(8)) &
  & -(CSUMPLQ10(9)+CSUMPLQ10(10)+CSUMPLQ10(11)+CSUMPLQ10(12)) &
  & +(CSUMPLQ10(13)+CSUMPLQ10(14)+CSUMPLQ10(15)+CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM179(N4,ID,IK)=(CSUMPLQMOM10(5,IK)+CSUMPLQMOM10(6,IK) &
  & +CSUMPLQMOM10(7,IK)+CSUMPLQMOM10(8,IK) &
  & -(CSUMPLQMOM10(13,IK)+CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK) &
  & +CSUMPLQMOM10(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE180(N4,ID)=(CSUMPLQ10(1)+GIOT*CSUMPLQ10(2)-CSUMPLQ10(3) &
  & -GIOT*CSUMPLQ10(4)+(CSUMPLQ10(5)+GIOT*CSUMPLQ10(6)-CSUMPLQ10(7) &
  & -GIOT*CSUMPLQ10(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM180(N4,ID,IK)=(CSUMPLQMOM10(1,IK) &
  & +GIOT*CSUMPLQMOM10(2,IK) &
  & -CSUMPLQMOM10(3,IK)-GIOT*CSUMPLQMOM10(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE181(N4,ID)=(CSUMPLQ10(1)+GIOT*CSUMPLQ10(2)-CSUMPLQ10(3) &
  & -GIOT*CSUMPLQ10(4)-(CSUMPLQ10(5)+GIOT*CSUMPLQ10(6)-CSUMPLQ10(7) &
  & -GIOT*CSUMPLQ10(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM181(N4,ID,IK)=(CSUMPLQMOM10(5,IK) &
  & +GIOT*CSUMPLQMOM10(6,IK) &
  & -CSUMPLQMOM10(7,IK)-GIOT*CSUMPLQMOM10(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE182(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3) &
  & -CSUMPLQ10(4) &
  & +(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8)) &
  & +(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12)) &
  & +(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM182(N4,ID,IK)=(CSUMPLQMOM10(1,IK)-CSUMPLQMOM10(2,IK) &
  & +CSUMPLQMOM10(3,IK)-CSUMPLQMOM10(4,IK) &
  & +(CSUMPLQMOM10(9,IK)-CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK) &
  & -CSUMPLQMOM10(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE183(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3) &
  & -CSUMPLQ10(4) &
  & -(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8)) &
  & +(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12)) &
  & -(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM183(N4,ID,IK)=(CSUMPLQMOM10(5,IK)-CSUMPLQMOM10(6,IK) &
  & +CSUMPLQMOM10(7,IK)-CSUMPLQMOM10(8,IK) &
  & +(CSUMPLQMOM10(13,IK)-CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK) &
  & -CSUMPLQMOM10(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE184(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3) &
  & -CSUMPLQ10(4) &
  & +(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8)) &
  & -(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12)) &
  & -(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM184(N4,ID,IK)=(CSUMPLQMOM10(1,IK)-CSUMPLQMOM10(2,IK) &
  & +CSUMPLQMOM10(3,IK)-CSUMPLQMOM10(4,IK) &
  & -(CSUMPLQMOM10(9,IK)-CSUMPLQMOM10(10,IK)+CSUMPLQMOM10(11,IK) &
  & -CSUMPLQMOM10(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE185(N4,ID)=(CSUMPLQ10(1)-CSUMPLQ10(2)+CSUMPLQ10(3) &
  & -CSUMPLQ10(4) &
  & -(CSUMPLQ10(5)-CSUMPLQ10(6)+CSUMPLQ10(7)-CSUMPLQ10(8)) &
  & -(CSUMPLQ10(9)-CSUMPLQ10(10)+CSUMPLQ10(11)-CSUMPLQ10(12)) &
  & +(CSUMPLQ10(13)-CSUMPLQ10(14)+CSUMPLQ10(15)-CSUMPLQ10(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM185(N4,ID,IK)=(CSUMPLQMOM10(5,IK)-CSUMPLQMOM10(6,IK) &
  & +CSUMPLQMOM10(7,IK)-CSUMPLQMOM10(8,IK) &
  & -(CSUMPLQMOM10(13,IK)-CSUMPLQMOM10(14,IK)+CSUMPLQMOM10(15,IK) &
  & -CSUMPLQMOM10(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 11
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE186(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3) &
  & +CSUMPLQ11(4) &
  & +(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8)) &
  & +(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12)) &
  & +(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM186(N4,ID,IK)=(CSUMPLQMOM11(1,IK)+CSUMPLQMOM11(2,IK) &
  & +CSUMPLQMOM11(3,IK)+CSUMPLQMOM11(4,IK) &
  & +(CSUMPLQMOM11(9,IK)+CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK) &
  & +CSUMPLQMOM11(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE187(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3) &
  & +CSUMPLQ11(4) &
  & -(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8)) &
  & +(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12)) &
  & -(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM187(N4,ID,IK)=(CSUMPLQMOM11(5,IK)+CSUMPLQMOM11(6,IK) &
  & +CSUMPLQMOM11(7,IK)+CSUMPLQMOM11(8,IK) &
  & +(CSUMPLQMOM11(13,IK)+CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK) &
  & +CSUMPLQMOM11(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE188(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3) &
  & +CSUMPLQ11(4) &
  & +(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8)) &
  & -(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12)) &
  & -(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM188(N4,ID,IK)=(CSUMPLQMOM11(1,IK)+CSUMPLQMOM11(2,IK) &
  & +CSUMPLQMOM11(3,IK)+CSUMPLQMOM11(4,IK) &
  & -(CSUMPLQMOM11(9,IK)+CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK) &
  & +CSUMPLQMOM11(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE189(N4,ID)=(CSUMPLQ11(1)+CSUMPLQ11(2)+CSUMPLQ11(3) &
  & +CSUMPLQ11(4) &
  & -(CSUMPLQ11(5)+CSUMPLQ11(6)+CSUMPLQ11(7)+CSUMPLQ11(8)) &
  & -(CSUMPLQ11(9)+CSUMPLQ11(10)+CSUMPLQ11(11)+CSUMPLQ11(12)) &
  & +(CSUMPLQ11(13)+CSUMPLQ11(14)+CSUMPLQ11(15)+CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM189(N4,ID,IK)=(CSUMPLQMOM11(5,IK)+CSUMPLQMOM11(6,IK) &
  & +CSUMPLQMOM11(7,IK)+CSUMPLQMOM11(8,IK) &
  & -(CSUMPLQMOM11(13,IK)+CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK) &
  & +CSUMPLQMOM11(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE190(N4,ID)=(CSUMPLQ11(1)+GIOT*CSUMPLQ11(2)-CSUMPLQ11(3) &
  & -GIOT*CSUMPLQ11(4)+(CSUMPLQ11(5)+GIOT*CSUMPLQ11(6)-CSUMPLQ11(7) &
  & -GIOT*CSUMPLQ11(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM190(N4,ID,IK)=(CSUMPLQMOM11(1,IK) &
  & +GIOT*CSUMPLQMOM11(2,IK) &
  & -CSUMPLQMOM11(3,IK)-GIOT*CSUMPLQMOM11(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE191(N4,ID)=(CSUMPLQ11(1)+GIOT*CSUMPLQ11(2)-CSUMPLQ11(3) &
  & -GIOT*CSUMPLQ11(4)-(CSUMPLQ11(5)+GIOT*CSUMPLQ11(6)-CSUMPLQ11(7) &
  & -GIOT*CSUMPLQ11(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM191(N4,ID,IK)=(CSUMPLQMOM11(5,IK) &
  & +GIOT*CSUMPLQMOM11(6,IK) &
  & -CSUMPLQMOM11(7,IK)-GIOT*CSUMPLQMOM11(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE192(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3) &
  & -CSUMPLQ11(4) &
  & +(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8)) &
  & +(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12)) &
  & +(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM192(N4,ID,IK)=(CSUMPLQMOM11(1,IK)-CSUMPLQMOM11(2,IK) &
  & +CSUMPLQMOM11(3,IK)-CSUMPLQMOM11(4,IK) &
  & +(CSUMPLQMOM11(9,IK)-CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK) &
  & -CSUMPLQMOM11(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE193(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3) &
  & -CSUMPLQ11(4) &
  & -(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8)) &
  & +(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12)) &
  & -(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM193(N4,ID,IK)=(CSUMPLQMOM11(5,IK)-CSUMPLQMOM11(6,IK) &
  & +CSUMPLQMOM11(7,IK)-CSUMPLQMOM11(8,IK) &
  & +(CSUMPLQMOM11(13,IK)-CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK) &
  & -CSUMPLQMOM11(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE194(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3) &
  & -CSUMPLQ11(4) &
  & +(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8)) &
  & -(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12)) &
  & -(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM194(N4,ID,IK)=(CSUMPLQMOM11(1,IK)-CSUMPLQMOM11(2,IK) &
  & +CSUMPLQMOM11(3,IK)-CSUMPLQMOM11(4,IK) &
  & -(CSUMPLQMOM11(9,IK)-CSUMPLQMOM11(10,IK)+CSUMPLQMOM11(11,IK) &
  & -CSUMPLQMOM11(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE195(N4,ID)=(CSUMPLQ11(1)-CSUMPLQ11(2)+CSUMPLQ11(3) &
  & -CSUMPLQ11(4) &
  & -(CSUMPLQ11(5)-CSUMPLQ11(6)+CSUMPLQ11(7)-CSUMPLQ11(8)) &
  & -(CSUMPLQ11(9)-CSUMPLQ11(10)+CSUMPLQ11(11)-CSUMPLQ11(12)) &
  & +(CSUMPLQ11(13)-CSUMPLQ11(14)+CSUMPLQ11(15)-CSUMPLQ11(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM195(N4,ID,IK)=(CSUMPLQMOM11(5,IK)-CSUMPLQMOM11(6,IK) &
  & +CSUMPLQMOM11(7,IK)-CSUMPLQMOM11(8,IK) &
  & -(CSUMPLQMOM11(13,IK)-CSUMPLQMOM11(14,IK)+CSUMPLQMOM11(15,IK) &
  & -CSUMPLQMOM11(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 12
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE196(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3) &
  & +CSUMPLQ12(4) &
  & +(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8)) &
  & +(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12)) &
  & +(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM196(N4,ID,IK)=(CSUMPLQMOM12(1,IK)+CSUMPLQMOM12(2,IK) &
  & +CSUMPLQMOM12(3,IK)+CSUMPLQMOM12(4,IK) &
  & +(CSUMPLQMOM12(9,IK)+CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK) &
  & +CSUMPLQMOM12(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE197(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3) &
  & +CSUMPLQ12(4) &
  & -(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8)) &
  & +(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12)) &
  & -(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM197(N4,ID,IK)=(CSUMPLQMOM12(5,IK)+CSUMPLQMOM12(6,IK) &
  & +CSUMPLQMOM12(7,IK)+CSUMPLQMOM12(8,IK) &
  & +(CSUMPLQMOM12(13,IK)+CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK) &
  & +CSUMPLQMOM12(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE198(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3) &
  & +CSUMPLQ12(4) &
  & +(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8)) &
  & -(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12)) &
  & -(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM198(N4,ID,IK)=(CSUMPLQMOM12(1,IK)+CSUMPLQMOM12(2,IK) &
  & +CSUMPLQMOM12(3,IK)+CSUMPLQMOM12(4,IK) &
  & -(CSUMPLQMOM12(9,IK)+CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK) &
  & +CSUMPLQMOM12(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE199(N4,ID)=(CSUMPLQ12(1)+CSUMPLQ12(2)+CSUMPLQ12(3) &
  & +CSUMPLQ12(4) &
  & -(CSUMPLQ12(5)+CSUMPLQ12(6)+CSUMPLQ12(7)+CSUMPLQ12(8)) &
  & -(CSUMPLQ12(9)+CSUMPLQ12(10)+CSUMPLQ12(11)+CSUMPLQ12(12)) &
  & +(CSUMPLQ12(13)+CSUMPLQ12(14)+CSUMPLQ12(15)+CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM199(N4,ID,IK)=(CSUMPLQMOM12(5,IK)+CSUMPLQMOM12(6,IK) &
  & +CSUMPLQMOM12(7,IK)+CSUMPLQMOM12(8,IK) &
  & -(CSUMPLQMOM12(13,IK)+CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK) &
  & +CSUMPLQMOM12(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE200(N4,ID)=(CSUMPLQ12(1)+GIOT*CSUMPLQ12(2)-CSUMPLQ12(3) &
  & -GIOT*CSUMPLQ12(4)+(CSUMPLQ12(5)+GIOT*CSUMPLQ12(6)-CSUMPLQ12(7) &
  & -GIOT*CSUMPLQ12(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM200(N4,ID,IK)=(CSUMPLQMOM12(1,IK) &
  & +GIOT*CSUMPLQMOM12(2,IK) &
  & -CSUMPLQMOM12(3,IK)-GIOT*CSUMPLQMOM12(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE201(N4,ID)=(CSUMPLQ12(1)+GIOT*CSUMPLQ12(2)-CSUMPLQ12(3) &
  & -GIOT*CSUMPLQ12(4)-(CSUMPLQ12(5)+GIOT*CSUMPLQ12(6)-CSUMPLQ12(7) &
  & -GIOT*CSUMPLQ12(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM201(N4,ID,IK)=(CSUMPLQMOM12(5,IK) &
  & +GIOT*CSUMPLQMOM12(6,IK) &
  & -CSUMPLQMOM12(7,IK)-GIOT*CSUMPLQMOM12(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE202(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3) &
  & -CSUMPLQ12(4) &
  & +(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8)) &
  & +(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12)) &
  & +(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM202(N4,ID,IK)=(CSUMPLQMOM12(1,IK)-CSUMPLQMOM12(2,IK) &
  & +CSUMPLQMOM12(3,IK)-CSUMPLQMOM12(4,IK) &
  & +(CSUMPLQMOM12(9,IK)-CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK) &
  & -CSUMPLQMOM12(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE203(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3) &
  & -CSUMPLQ12(4) &
  & -(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8)) &
  & +(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12)) &
  & -(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM203(N4,ID,IK)=(CSUMPLQMOM12(5,IK)-CSUMPLQMOM12(6,IK) &
  & +CSUMPLQMOM12(7,IK)-CSUMPLQMOM12(8,IK) &
  & +(CSUMPLQMOM12(13,IK)-CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK) &
  & -CSUMPLQMOM12(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE204(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3) &
  & -CSUMPLQ12(4) &
  & +(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8)) &
  & -(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12)) &
  & -(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM204(N4,ID,IK)=(CSUMPLQMOM12(1,IK)-CSUMPLQMOM12(2,IK) &
  & +CSUMPLQMOM12(3,IK)-CSUMPLQMOM12(4,IK) &
  & -(CSUMPLQMOM12(9,IK)-CSUMPLQMOM12(10,IK)+CSUMPLQMOM12(11,IK) &
  & -CSUMPLQMOM12(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE205(N4,ID)=(CSUMPLQ12(1)-CSUMPLQ12(2)+CSUMPLQ12(3) &
  & -CSUMPLQ12(4) &
  & -(CSUMPLQ12(5)-CSUMPLQ12(6)+CSUMPLQ12(7)-CSUMPLQ12(8)) &
  & -(CSUMPLQ12(9)-CSUMPLQ12(10)+CSUMPLQ12(11)-CSUMPLQ12(12)) &
  & +(CSUMPLQ12(13)-CSUMPLQ12(14)+CSUMPLQ12(15)-CSUMPLQ12(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM205(N4,ID,IK)=(CSUMPLQMOM12(5,IK)-CSUMPLQMOM12(6,IK) &
  & +CSUMPLQMOM12(7,IK)-CSUMPLQMOM12(8,IK) &
  & -(CSUMPLQMOM12(13,IK)-CSUMPLQMOM12(14,IK)+CSUMPLQMOM12(15,IK) &
  & -CSUMPLQMOM12(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 13
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE206(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3) &
  & +CSUMPLQ13(4) &
  & +(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8)) &
  & +(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12)) &
  & +(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM206(N4,ID,IK)=(CSUMPLQMOM13(1,IK)+CSUMPLQMOM13(2,IK) &
  & +CSUMPLQMOM13(3,IK)+CSUMPLQMOM13(4,IK) &
  & +(CSUMPLQMOM13(9,IK)+CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK) &
  & +CSUMPLQMOM13(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE207(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3) &
  & +CSUMPLQ13(4) &
  & -(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8)) &
  & +(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12)) &
  & -(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM207(N4,ID,IK)=(CSUMPLQMOM13(5,IK)+CSUMPLQMOM13(6,IK) &
  & +CSUMPLQMOM13(7,IK)+CSUMPLQMOM13(8,IK) &
  & +(CSUMPLQMOM13(13,IK)+CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK) &
  & +CSUMPLQMOM13(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE208(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3) &
  & +CSUMPLQ13(4) &
  & +(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8)) &
  & -(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12)) &
  & -(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM208(N4,ID,IK)=(CSUMPLQMOM13(1,IK)+CSUMPLQMOM13(2,IK) &
  & +CSUMPLQMOM13(3,IK)+CSUMPLQMOM13(4,IK) &
  & -(CSUMPLQMOM13(9,IK)+CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK) &
  & +CSUMPLQMOM13(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE209(N4,ID)=(CSUMPLQ13(1)+CSUMPLQ13(2)+CSUMPLQ13(3) &
  & +CSUMPLQ13(4) &
  & -(CSUMPLQ13(5)+CSUMPLQ13(6)+CSUMPLQ13(7)+CSUMPLQ13(8)) &
  & -(CSUMPLQ13(9)+CSUMPLQ13(10)+CSUMPLQ13(11)+CSUMPLQ13(12)) &
  & +(CSUMPLQ13(13)+CSUMPLQ13(14)+CSUMPLQ13(15)+CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM209(N4,ID,IK)=(CSUMPLQMOM13(5,IK)+CSUMPLQMOM13(6,IK) &
  & +CSUMPLQMOM13(7,IK)+CSUMPLQMOM13(8,IK) &
  & -(CSUMPLQMOM13(13,IK)+CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK) &
  & +CSUMPLQMOM13(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE210(N4,ID)=(CSUMPLQ13(1)+GIOT*CSUMPLQ13(2)-CSUMPLQ13(3) &
  & -GIOT*CSUMPLQ13(4)+(CSUMPLQ13(5)+GIOT*CSUMPLQ13(6)-CSUMPLQ13(7) &
  & -GIOT*CSUMPLQ13(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM210(N4,ID,IK)=(CSUMPLQMOM13(1,IK) &
  & +GIOT*CSUMPLQMOM13(2,IK) &
  & -CSUMPLQMOM13(3,IK)-GIOT*CSUMPLQMOM13(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE211(N4,ID)=(CSUMPLQ13(1)+GIOT*CSUMPLQ13(2)-CSUMPLQ13(3) &
  & -GIOT*CSUMPLQ13(4)-(CSUMPLQ13(5)+GIOT*CSUMPLQ13(6)-CSUMPLQ13(7) &
  & -GIOT*CSUMPLQ13(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM211(N4,ID,IK)=(CSUMPLQMOM13(5,IK) &
  & +GIOT*CSUMPLQMOM13(6,IK) &
  & -CSUMPLQMOM13(7,IK)-GIOT*CSUMPLQMOM13(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE212(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3) &
  & -CSUMPLQ13(4) &
  & +(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8)) &
  & +(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12)) &
  & +(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM212(N4,ID,IK)=(CSUMPLQMOM13(1,IK)-CSUMPLQMOM13(2,IK) &
  & +CSUMPLQMOM13(3,IK)-CSUMPLQMOM13(4,IK) &
  & +(CSUMPLQMOM13(9,IK)-CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK) &
  & -CSUMPLQMOM13(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE213(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3) &
  & -CSUMPLQ13(4) &
  & -(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8)) &
  & +(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12)) &
  & -(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM213(N4,ID,IK)=(CSUMPLQMOM13(5,IK)-CSUMPLQMOM13(6,IK) &
  & +CSUMPLQMOM13(7,IK)-CSUMPLQMOM13(8,IK) &
  & +(CSUMPLQMOM13(13,IK)-CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK) &
  & -CSUMPLQMOM13(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE214(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3) &
  & -CSUMPLQ13(4) &
  & +(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8)) &
  & -(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12)) &
  & -(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM214(N4,ID,IK)=(CSUMPLQMOM13(1,IK)-CSUMPLQMOM13(2,IK) &
  & +CSUMPLQMOM13(3,IK)-CSUMPLQMOM13(4,IK) &
  & -(CSUMPLQMOM13(9,IK)-CSUMPLQMOM13(10,IK)+CSUMPLQMOM13(11,IK) &
  & -CSUMPLQMOM13(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE215(N4,ID)=(CSUMPLQ13(1)-CSUMPLQ13(2)+CSUMPLQ13(3) &
  & -CSUMPLQ13(4) &
  & -(CSUMPLQ13(5)-CSUMPLQ13(6)+CSUMPLQ13(7)-CSUMPLQ13(8)) &
  & -(CSUMPLQ13(9)-CSUMPLQ13(10)+CSUMPLQ13(11)-CSUMPLQ13(12)) &
  & +(CSUMPLQ13(13)-CSUMPLQ13(14)+CSUMPLQ13(15)-CSUMPLQ13(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM215(N4,ID,IK)=(CSUMPLQMOM13(5,IK)-CSUMPLQMOM13(6,IK) &
  & +CSUMPLQMOM13(7,IK)-CSUMPLQMOM13(8,IK) &
  & -(CSUMPLQMOM13(13,IK)-CSUMPLQMOM13(14,IK)+CSUMPLQMOM13(15,IK) &
  & -CSUMPLQMOM13(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 14
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE216(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3) &
  & +CSUMPLQ14(4) &
  & +(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8)) &
  & +(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12)) &
  & +(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM216(N4,ID,IK)=(CSUMPLQMOM14(1,IK)+CSUMPLQMOM14(2,IK) &
  & +CSUMPLQMOM14(3,IK)+CSUMPLQMOM14(4,IK) &
  & +(CSUMPLQMOM14(9,IK)+CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK) &
  & +CSUMPLQMOM14(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE217(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3) &
  & +CSUMPLQ14(4) &
  & -(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8)) &
  & +(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12)) &
  & -(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM217(N4,ID,IK)=(CSUMPLQMOM14(5,IK)+CSUMPLQMOM14(6,IK) &
  & +CSUMPLQMOM14(7,IK)+CSUMPLQMOM14(8,IK) &
  & +(CSUMPLQMOM14(13,IK)+CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK) &
  & +CSUMPLQMOM14(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE218(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3) &
  & +CSUMPLQ14(4) &
  & +(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8)) &
  & -(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12)) &
  & -(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM218(N4,ID,IK)=(CSUMPLQMOM14(1,IK)+CSUMPLQMOM14(2,IK) &
  & +CSUMPLQMOM14(3,IK)+CSUMPLQMOM14(4,IK) &
  & -(CSUMPLQMOM14(9,IK)+CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK) &
  & +CSUMPLQMOM14(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE219(N4,ID)=(CSUMPLQ14(1)+CSUMPLQ14(2)+CSUMPLQ14(3) &
  & +CSUMPLQ14(4) &
  & -(CSUMPLQ14(5)+CSUMPLQ14(6)+CSUMPLQ14(7)+CSUMPLQ14(8)) &
  & -(CSUMPLQ14(9)+CSUMPLQ14(10)+CSUMPLQ14(11)+CSUMPLQ14(12)) &
  & +(CSUMPLQ14(13)+CSUMPLQ14(14)+CSUMPLQ14(15)+CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM219(N4,ID,IK)=(CSUMPLQMOM14(5,IK)+CSUMPLQMOM14(6,IK) &
  & +CSUMPLQMOM14(7,IK)+CSUMPLQMOM14(8,IK) &
  & -(CSUMPLQMOM14(13,IK)+CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK) &
  & +CSUMPLQMOM14(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE220(N4,ID)=(CSUMPLQ14(1)+GIOT*CSUMPLQ14(2)-CSUMPLQ14(3) &
  & -GIOT*CSUMPLQ14(4)+(CSUMPLQ14(5)+GIOT*CSUMPLQ14(6)-CSUMPLQ14(7) &
  & -GIOT*CSUMPLQ14(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM220(N4,ID,IK)=(CSUMPLQMOM14(1,IK) &
  & +GIOT*CSUMPLQMOM14(2,IK) &
  & -CSUMPLQMOM14(3,IK)-GIOT*CSUMPLQMOM14(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE221(N4,ID)=(CSUMPLQ14(1)+GIOT*CSUMPLQ14(2)-CSUMPLQ14(3) &
  & -GIOT*CSUMPLQ14(4)-(CSUMPLQ14(5)+GIOT*CSUMPLQ14(6)-CSUMPLQ14(7) &
  & -GIOT*CSUMPLQ14(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM221(N4,ID,IK)=(CSUMPLQMOM14(5,IK) &
  & +GIOT*CSUMPLQMOM14(6,IK) &
  & -CSUMPLQMOM14(7,IK)-GIOT*CSUMPLQMOM14(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE222(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3) &
  & -CSUMPLQ14(4) &
  & +(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8)) &
  & +(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12)) &
  & +(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM222(N4,ID,IK)=(CSUMPLQMOM14(1,IK)-CSUMPLQMOM14(2,IK) &
  & +CSUMPLQMOM14(3,IK)-CSUMPLQMOM14(4,IK) &
  & +(CSUMPLQMOM14(9,IK)-CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK) &
  & -CSUMPLQMOM14(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE223(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3) &
  & -CSUMPLQ14(4) &
  & -(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8)) &
  & +(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12)) &
  & -(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM223(N4,ID,IK)=(CSUMPLQMOM14(5,IK)-CSUMPLQMOM14(6,IK) &
  & +CSUMPLQMOM14(7,IK)-CSUMPLQMOM14(8,IK) &
  & +(CSUMPLQMOM14(13,IK)-CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK) &
  & -CSUMPLQMOM14(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE224(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3) &
  & -CSUMPLQ14(4) &
  & +(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8)) &
  & -(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12)) &
  & -(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM224(N4,ID,IK)=(CSUMPLQMOM14(1,IK)-CSUMPLQMOM14(2,IK) &
  & +CSUMPLQMOM14(3,IK)-CSUMPLQMOM14(4,IK) &
  & -(CSUMPLQMOM14(9,IK)-CSUMPLQMOM14(10,IK)+CSUMPLQMOM14(11,IK) &
  & -CSUMPLQMOM14(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE225(N4,ID)=(CSUMPLQ14(1)-CSUMPLQ14(2)+CSUMPLQ14(3) &
  & -CSUMPLQ14(4) &
  & -(CSUMPLQ14(5)-CSUMPLQ14(6)+CSUMPLQ14(7)-CSUMPLQ14(8)) &
  & -(CSUMPLQ14(9)-CSUMPLQ14(10)+CSUMPLQ14(11)-CSUMPLQ14(12)) &
  & +(CSUMPLQ14(13)-CSUMPLQ14(14)+CSUMPLQ14(15)-CSUMPLQ14(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM225(N4,ID,IK)=(CSUMPLQMOM14(5,IK)-CSUMPLQMOM14(6,IK) &
  & +CSUMPLQMOM14(7,IK)-CSUMPLQMOM14(8,IK) &
  & -(CSUMPLQMOM14(13,IK)-CSUMPLQMOM14(14,IK)+CSUMPLQMOM14(15,IK) &
  & -CSUMPLQMOM14(16,IK)))*ADIV2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PLAQUETTE OPERATORS 15
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!**********************************************************************      
!     J=0, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE226(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3) &
  & +CSUMPLQ15(4) &
  & +(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8)) &
  & +(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12)) &
  & +(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM226(N4,ID,IK)=(CSUMPLQMOM15(1,IK)+CSUMPLQMOM15(2,IK) &
  & +CSUMPLQMOM15(3,IK)+CSUMPLQMOM15(4,IK) &
  & +(CSUMPLQMOM15(9,IK)+CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK) &
  & +CSUMPLQMOM15(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE227(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3) &
  & +CSUMPLQ15(4) &
  & -(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8)) &
  & +(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12)) &
  & -(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=+, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM227(N4,ID,IK)=(CSUMPLQMOM15(5,IK)+CSUMPLQMOM15(6,IK) &
  & +CSUMPLQMOM15(7,IK)+CSUMPLQMOM15(8,IK) &
  & +(CSUMPLQMOM15(13,IK)+CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK) &
  & +CSUMPLQMOM15(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE228(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3) &
  & +CSUMPLQ15(4) &
  & +(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8)) &
  & -(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12)) &
  & -(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM228(N4,ID,IK)=(CSUMPLQMOM15(1,IK)+CSUMPLQMOM15(2,IK) &
  & +CSUMPLQMOM15(3,IK)+CSUMPLQMOM15(4,IK) &
  & -(CSUMPLQMOM15(9,IK)+CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK) &
  & +CSUMPLQMOM15(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=0, Pp=-, Pr=-, q=0
!**********************************************************************
ALINE229(N4,ID)=(CSUMPLQ15(1)+CSUMPLQ15(2)+CSUMPLQ15(3) &
  & +CSUMPLQ15(4) &
  & -(CSUMPLQ15(5)+CSUMPLQ15(6)+CSUMPLQ15(7)+CSUMPLQ15(8)) &
  & -(CSUMPLQ15(9)+CSUMPLQ15(10)+CSUMPLQ15(11)+CSUMPLQ15(12)) &
  & +(CSUMPLQ15(13)+CSUMPLQ15(14)+CSUMPLQ15(15)+CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=0, Pp=-, q=1,2
!**********************************************************************
DO IK=1, 2
ALINEMOM229(N4,ID,IK)=(CSUMPLQMOM15(5,IK)+CSUMPLQMOM15(6,IK) &
  & +CSUMPLQMOM15(7,IK)+CSUMPLQMOM15(8,IK) &
  & -(CSUMPLQMOM15(13,IK)+CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK) &
  & +CSUMPLQMOM15(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=1, Pr=+, q=0
!**********************************************************************            
ALINE230(N4,ID)=(CSUMPLQ15(1)+GIOT*CSUMPLQ15(2)-CSUMPLQ15(3) &
  & -GIOT*CSUMPLQ15(4)+(CSUMPLQ15(5)+GIOT*CSUMPLQ15(6)-CSUMPLQ15(7) &
  & -GIOT*CSUMPLQ15(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2
!**********************************************************************                  
DO IK=1, 2
ALINEMOM230(N4,ID,IK)=(CSUMPLQMOM15(1,IK) &
  & +GIOT*CSUMPLQMOM15(2,IK) &
  & -CSUMPLQMOM15(3,IK)-GIOT*CSUMPLQMOM15(4,IK))*ADIV1
end do
!**********************************************************************      
!     J=1, Pr=-, q=0
!**********************************************************************                  
ALINE231(N4,ID)=(CSUMPLQ15(1)+GIOT*CSUMPLQ15(2)-CSUMPLQ15(3) &
  & -GIOT*CSUMPLQ15(4)-(CSUMPLQ15(5)+GIOT*CSUMPLQ15(6)-CSUMPLQ15(7) &
  & -GIOT*CSUMPLQ15(8)))*ADIV2
!**********************************************************************      
!     J=1, q=1,2  Here!
!**********************************************************************                        
DO IK=1, 2
ALINEMOM231(N4,ID,IK)=(CSUMPLQMOM15(5,IK) &
  & +GIOT*CSUMPLQMOM15(6,IK) &
  & -CSUMPLQMOM15(7,IK)-GIOT*CSUMPLQMOM15(8,IK))*ADIV1
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=+, q=0
!**********************************************************************
ALINE232(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3) &
  & -CSUMPLQ15(4) &
  & +(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8)) &
  & +(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12)) &
  & +(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM232(N4,ID,IK)=(CSUMPLQMOM15(1,IK)-CSUMPLQMOM15(2,IK) &
  & +CSUMPLQMOM15(3,IK)-CSUMPLQMOM15(4,IK) &
  & +(CSUMPLQMOM15(9,IK)-CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK) &
  & -CSUMPLQMOM15(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=+, Pr=-, q=0
!**********************************************************************
ALINE233(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3) &
  & -CSUMPLQ15(4) &
  & -(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8)) &
  & +(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12)) &
  & -(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=+, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM233(N4,ID,IK)=(CSUMPLQMOM15(5,IK)-CSUMPLQMOM15(6,IK) &
  & +CSUMPLQMOM15(7,IK)-CSUMPLQMOM15(8,IK) &
  & +(CSUMPLQMOM15(13,IK)-CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK) &
  & -CSUMPLQMOM15(16,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=+, q=0
!**********************************************************************
ALINE234(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3) &
  & -CSUMPLQ15(4) &
  & +(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8)) &
  & -(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12)) &
  & -(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM234(N4,ID,IK)=(CSUMPLQMOM15(1,IK)-CSUMPLQMOM15(2,IK) &
  & +CSUMPLQMOM15(3,IK)-CSUMPLQMOM15(4,IK) &
  & -(CSUMPLQMOM15(9,IK)-CSUMPLQMOM15(10,IK)+CSUMPLQMOM15(11,IK) &
  & -CSUMPLQMOM15(12,IK)))*ADIV2
end do
!**********************************************************************      
!     J=2, Pp=-, Pr=-, q=0  test
!**********************************************************************
ALINE235(N4,ID)=(CSUMPLQ15(1)-CSUMPLQ15(2)+CSUMPLQ15(3) &
  & -CSUMPLQ15(4) &
  & -(CSUMPLQ15(5)-CSUMPLQ15(6)+CSUMPLQ15(7)-CSUMPLQ15(8)) &
  & -(CSUMPLQ15(9)-CSUMPLQ15(10)+CSUMPLQ15(11)-CSUMPLQ15(12)) &
  & +(CSUMPLQ15(13)-CSUMPLQ15(14)+CSUMPLQ15(15)-CSUMPLQ15(16)))*ADIV3
!**********************************************************************      
!     J=2, Pp=-, q=0
!**********************************************************************
DO IK=1, 2
ALINEMOM235(N4,ID,IK)=(CSUMPLQMOM15(5,IK)-CSUMPLQMOM15(6,IK) &
  & +CSUMPLQMOM15(7,IK)-CSUMPLQMOM15(8,IK) &
  & -(CSUMPLQMOM15(13,IK)-CSUMPLQMOM15(14,IK)+CSUMPLQMOM15(15,IK) &
  & -CSUMPLQMOM15(16,IK)))*ADIV2
end do
!**********************************************************************
2 CONTINUE
end do
!**********************************************************************
RETURN
end subroutine THERML1
!**********************************************************************
!**********************************************************************
!     HERE I CREATE THE DIAGONAL OPERATORS
!**********************************************************************
!********************************************************************
!   HERE, WE CREATE HERMITIAN CORRELATORS!!!!                  
!   COR(I,J)=F(I,T+DT)*F^+(J,T)+F^+(J,T+DT)*F(I,T)             
!   COR(J,I)=F(J,T+DT)*F^+(I,T)+F^+(I,T+DT)*F(J,T)             
!   COR^+(J,I)=F^+(J,T+DT)*F(I,T)+F(I,T+DT)*F^+(J,T)           
!   => COR(I,J)=COR^+(J,I)!!                                   
!*********************************************************************
SUBROUTINE POT_OLD
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
PARAMETER(NUMBIN=2,LMAX=LX4/2+1)
PARAMETER(LMAXIR=3)
!******************************************************************************
!******************************************************************************     
PARAMETER(I0PP=38,I0PM=10,I0MP=15,I0MM=22) !TBC!
PARAMETER(I1P=34,I1M=32) !TBC!
PARAMETER(I2PP=27,I2PM=20,I2MP=25,I2MM=12) !TBC!
!----------------------------------------------
PARAMETER(I0PMOM=47,I0MMOM=37) !TBC!
PARAMETER(I1MOM=66) !TBC!
PARAMETER(I2PMOM=47,I2MMOM=37) !TBC!
!-----------------------------------
PARAMETER(NOPJ0PP=I0PP*IBLOK,NOPJ0PM=I0PM*IBLOK)
PARAMETER(NOPJ0MP=I0MP*IBLOK,NOPJ0MM=I0MM*IBLOK)
PARAMETER(NOPJ1P=I1P*IBLOK,NOPJ1M=I1M*IBLOK)
PARAMETER(NOPJ2PP=I2PP*IBLOK,NOPJ2PM=I2PM*IBLOK)
PARAMETER(NOPJ2MP=I2MP*IBLOK,NOPJ2MM=I2MM*IBLOK)
!-----------------------------------------------------
PARAMETER(NOPJ0PMOM=I0PMOM*IBLOK,NOPJ0MMOM=I0MMOM*IBLOK)
PARAMETER(NOPJ1MOM=I1MOM*IBLOK)
PARAMETER(NOPJ2PMOM=I2PMOM*IBLOK,NOPJ2MMOM=I2MMOM*IBLOK)
!-------------------------------------------------------------
PARAMETER(NOPFULJ0=NOPJ0PP+NOPJ0PM+NOPJ0MP+NOPJ0MM)
PARAMETER(NOPFULJ1=NOPJ1P+NOPJ1M)
PARAMETER(NOPFULJ2=NOPJ2PP+NOPJ2PM+NOPJ2MP+NOPJ2MM)
PARAMETER(NTOTAL=NOPFULJ0+NOPFULJ1+NOPFULJ2)
!-------------------------------------------------
PARAMETER(NOPFULJ0MOM=NOPJ0PMOM+NOPJ0MMOM)
PARAMETER(NOPFULJ1MOM=NOPJ1MOM)
PARAMETER(NOPFULJ2MOM=NOPJ2PMOM+NOPJ2MMOM)
PARAMETER(NTOTALMOM=NOPFULJ0MOM+NOPFULJ1MOM+NOPFULJ2MOM)
!----------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************  
!     
COMMON/PBLOCKP/ACORLP(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLP(NUMBIN,NTOTAL)
!
COMMON/PBLOCKMOMP/ACORLPMOM(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM) &
  & ,AVACLPMOM(2,NUMBIN,NTOTALMOM)
!
COMMON/PBLOCKJ0PP/ACORLJ0PP(NUMBIN,LMAX,NOPJ0PP,NOPJ0PP) &
  & ,AVACLJ0PP(NUMBIN,NOPJ0PP)
COMMON/PBLOCKJ0PM/ACORLJ0PM(NUMBIN,LMAX,NOPJ0PM,NOPJ0PM) &
  & ,AVACLJ0PM(NUMBIN,NOPJ0PM)
COMMON/PBLOCKJ0MP/ACORLJ0MP(NUMBIN,LMAX,NOPJ0MP,NOPJ0MP) &
  & ,AVACLJ0MP(NUMBIN,NOPJ0MP)
COMMON/PBLOCKJ0MM/ACORLJ0MM(NUMBIN,LMAX,NOPJ0MM,NOPJ0MM) &
  & ,AVACLJ0MM(NUMBIN,NOPJ0MM)
COMMON/PBLOCKJ1P/ACORLJ1P(NUMBIN,LMAX,NOPJ1P,NOPJ1P) &
  & ,AVACLJ1P(NUMBIN,NOPJ1P)
COMMON/PBLOCKJ1M/ACORLJ1M(NUMBIN,LMAX,NOPJ1M,NOPJ1M) &
  & ,AVACLJ1M(NUMBIN,NOPJ1M)
COMMON/PBLOCKJ2PP/ACORLJ2PP(NUMBIN,LMAX,NOPJ2PP,NOPJ2PP) &
  & ,AVACLJ2PP(NUMBIN,NOPJ2PP)
COMMON/PBLOCKJ2PM/ACORLJ2PM(NUMBIN,LMAX,NOPJ2PM,NOPJ2PM) &
  & ,AVACLJ2PM(NUMBIN,NOPJ2PM)
COMMON/PBLOCKJ2MP/ACORLJ2MP(NUMBIN,LMAX,NOPJ2MP,NOPJ2MP) &
  & ,AVACLJ2MP(NUMBIN,NOPJ2MP)
COMMON/PBLOCKJ2MM/ACORLJ2MM(NUMBIN,LMAX,NOPJ2MM,NOPJ2MM) &
  & ,AVACLJ2MM(NUMBIN,NOPJ2MM)
!
COMMON/PBLOCKJ0/ACORLJ0(NUMBIN,LMAX,NOPFULJ0,NOPFULJ0) &
  & ,AVACLJ0(NUMBIN,NOPFULJ0)
COMMON/PBLOCKJ1/ACORLJ1(NUMBIN,LMAX,NOPFULJ1,NOPFULJ1) &
  & ,AVACLJ1(NUMBIN,NOPFULJ1)
COMMON/PBLOCKJ2/ACORLJ2(NUMBIN,LMAX,NOPFULJ2,NOPFULJ2) &
  & ,AVACLJ2(NUMBIN,NOPFULJ2)
!
COMMON/PBLOCKJALL/ACORLALL(NUMBIN,LMAXIR,NTOTAL,NTOTAL) &
  & ,AVACLALL(NUMBIN,NTOTAL)
!
COMMON/PBLOCKMOMJ0P/ &
  & ACORLPMOMJ0P(2,NUMBIN,LMAX,NOPJ0PMOM,NOPJ0PMOM) &
  & ,AVACLMOMJ0P(2,NUMBIN,NOPJ0PMOM)
COMMON/PBLOCKMOMJ0M/ &
  & ACORLPMOMJ0M(2,NUMBIN,LMAX,NOPJ0MMOM,NOPJ0MMOM) &
  & ,AVACLMOMJ0M(2,NUMBIN,NOPJ0MMOM)
COMMON/PBLOCKMOMJ1/ &
  & ACORLPMOMJ1(2,NUMBIN,LMAX,NOPJ1MOM,NOPJ1MOM) &
  & ,AVACLMOMJ1(2,NUMBIN,NOPJ1MOM)
COMMON/PBLOCKMOMJ2P/ &
  & ACORLPMOMJ2P(2,NUMBIN,LMAX,NOPJ2PMOM,NOPJ2PMOM) &
  & ,AVACLMOMJ2P(2,NUMBIN,NOPJ2PMOM)
COMMON/PBLOCKMOMJ2M/ &
  & ACORLPMOMJ2M(2,NUMBIN,LMAX,NOPJ2MMOM,NOPJ2MMOM) &
  & ,AVACLMOMJ2M(2,NUMBIN,NOPJ2MMOM)
!
COMMON/PBLOCKMOMJ0/ &
  & ACORLPMOMJ0(2,NUMBIN,LMAX,NOPFULJ0MOM,NOPFULJ0MOM) &
  & ,AVACLMOMJ0(2,NUMBIN,NOPFULJ0MOM)
COMMON/PBLOCKMOMJ2/ &
  & ACORLPMOMJ2(2,NUMBIN,LMAX,NOPFULJ2MOM,NOPFULJ2MOM) &
  & ,AVACLMOMJ2(2,NUMBIN,NOPFULJ2MOM)
!
COMMON/PBLOCKMOMJALL/ &
  & ACORLPMOMJALL(2,NUMBIN,LMAXIR,NTOTALMOM,NTOTALMOM) &
  & ,AVACLMOMALL(2,NUMBIN,NTOTALMOM)
!      
!***************************************************************
COMMON/LINES/ALINE1(LX4,IBLOK),ALINE2(LX4,IBLOK) &
  & ,ALINE3(LX4,IBLOK),ALINE4(LX4,IBLOK),ALINE5(LX4,IBLOK) &
  & ,ALINE6(LX4,IBLOK),ALINE7(LX4,IBLOK),ALINE8(LX4,IBLOK) &
  & ,ALINE9(LX4,IBLOK),ALINE10(LX4,IBLOK),ALINE11(LX4,IBLOK) &
  & ,ALINE12(LX4,IBLOK),ALINE13(LX4,IBLOK),ALINE14(LX4,IBLOK) &
  & ,ALINE15(LX4,IBLOK),ALINE16(LX4,IBLOK),ALINE17(LX4,IBLOK) &
  & ,ALINE18(LX4,IBLOK),ALINE19(LX4,IBLOK),ALINE20(LX4,IBLOK) &
  & ,ALINE21(LX4,IBLOK),ALINE22(LX4,IBLOK),ALINE23(LX4,IBLOK) &
  & ,ALINE24(LX4,IBLOK),ALINE25(LX4,IBLOK),ALINE26(LX4,IBLOK) &
  & ,ALINE27(LX4,IBLOK),ALINE28(LX4,IBLOK),ALINE29(LX4,IBLOK) &
  & ,ALINE30(LX4,IBLOK),ALINE31(LX4,IBLOK),ALINE32(LX4,IBLOK) &
  & ,ALINE33(LX4,IBLOK),ALINE34(LX4,IBLOK),ALINE35(LX4,IBLOK) &
  & ,ALINE36(LX4,IBLOK),ALINE37(LX4,IBLOK),ALINE38(LX4,IBLOK) &
  & ,ALINE39(LX4,IBLOK),ALINE40(LX4,IBLOK),ALINE41(LX4,IBLOK) &
  & ,ALINE42(LX4,IBLOK),ALINE43(LX4,IBLOK),ALINE44(LX4,IBLOK) &
  & ,ALINE45(LX4,IBLOK),ALINE46(LX4,IBLOK),ALINE47(LX4,IBLOK) &
  & ,ALINE48(LX4,IBLOK),ALINE49(LX4,IBLOK),ALINE50(LX4,IBLOK) &
  & ,ALINE51(LX4,IBLOK),ALINE52(LX4,IBLOK),ALINE53(LX4,IBLOK) &
  & ,ALINE54(LX4,IBLOK),ALINE55(LX4,IBLOK),ALINE56(LX4,IBLOK) &
  & ,ALINE57(LX4,IBLOK),ALINE58(LX4,IBLOK),ALINE59(LX4,IBLOK) &
  & ,ALINE60(LX4,IBLOK),ALINE61(LX4,IBLOK),ALINE62(LX4,IBLOK) &
  & ,ALINE63(LX4,IBLOK),ALINE64(LX4,IBLOK),ALINE65(LX4,IBLOK) &
  & ,ALINE66(LX4,IBLOK),ALINE67(LX4,IBLOK),ALINE68(LX4,IBLOK) &
  & ,ALINE69(LX4,IBLOK),ALINE70(LX4,IBLOK),ALINE71(LX4,IBLOK) &
  & ,ALINE72(LX4,IBLOK),ALINE73(LX4,IBLOK),ALINE74(LX4,IBLOK) &
  & ,ALINE75(LX4,IBLOK),ALINE76(LX4,IBLOK),ALINE77(LX4,IBLOK) &
  & ,ALINE78(LX4,IBLOK),ALINE79(LX4,IBLOK),ALINE80(LX4,IBLOK) &
  & ,ALINE81(LX4,IBLOK),ALINE82(LX4,IBLOK),ALINE83(LX4,IBLOK) &
  & ,ALINE84(LX4,IBLOK),ALINE85(LX4,IBLOK),ALINE86(LX4,IBLOK) &
  & ,ALINE87(LX4,IBLOK),ALINE88(LX4,IBLOK),ALINE89(LX4,IBLOK) &
  & ,ALINE90(LX4,IBLOK),ALINE91(LX4,IBLOK),ALINE92(LX4,IBLOK) &
  & ,ALINE93(LX4,IBLOK),ALINE94(LX4,IBLOK),ALINE95(LX4,IBLOK) &
  & ,ALINE96(LX4,IBLOK),ALINE97(LX4,IBLOK),ALINE98(LX4,IBLOK) &
  & ,ALINE99(LX4,IBLOK),ALINE100(LX4,IBLOK),ALINE101(LX4,IBLOK) &
  & ,ALINE102(LX4,IBLOK),ALINE103(LX4,IBLOK),ALINE104(LX4,IBLOK) &
  & ,ALINE105(LX4,IBLOK),ALINE106(LX4,IBLOK),ALINE107(LX4,IBLOK) &
  & ,ALINE108(LX4,IBLOK),ALINE109(LX4,IBLOK),ALINE110(LX4,IBLOK) &
  & ,ALINE111(LX4,IBLOK),ALINE112(LX4,IBLOK),ALINE113(LX4,IBLOK) &
  & ,ALINE114(LX4,IBLOK),ALINE115(LX4,IBLOK),ALINE116(LX4,IBLOK) &
  & ,ALINE117(LX4,IBLOK),ALINE118(LX4,IBLOK),ALINE119(LX4,IBLOK) &
  & ,ALINE120(LX4,IBLOK),ALINE121(LX4,IBLOK),ALINE122(LX4,IBLOK) &
  & ,ALINE123(LX4,IBLOK),ALINE124(LX4,IBLOK),ALINE125(LX4,IBLOK) &
  & ,ALINE126(LX4,IBLOK),ALINE127(LX4,IBLOK),ALINE128(LX4,IBLOK) &
  & ,ALINE129(LX4,IBLOK),ALINE130(LX4,IBLOK),ALINE131(LX4,IBLOK) &
  & ,ALINE132(LX4,IBLOK),ALINE133(LX4,IBLOK),ALINE134(LX4,IBLOK) &
  & ,ALINE135(LX4,IBLOK),ALINE136(LX4,IBLOK),ALINE137(LX4,IBLOK) &
  & ,ALINE138(LX4,IBLOK),ALINE139(LX4,IBLOK),ALINE140(LX4,IBLOK) &
  & ,ALINE141(LX4,IBLOK),ALINE142(LX4,IBLOK),ALINE143(LX4,IBLOK) &
  & ,ALINE144(LX4,IBLOK),ALINE145(LX4,IBLOK),ALINE146(LX4,IBLOK) &
  & ,ALINE147(LX4,IBLOK),ALINE148(LX4,IBLOK),ALINE149(LX4,IBLOK) &
  & ,ALINE150(LX4,IBLOK),ALINE151(LX4,IBLOK),ALINE152(LX4,IBLOK) &
  & ,ALINE153(LX4,IBLOK),ALINE154(LX4,IBLOK),ALINE155(LX4,IBLOK) &
  & ,ALINE156(LX4,IBLOK),ALINE157(LX4,IBLOK),ALINE158(LX4,IBLOK) &
  & ,ALINE159(LX4,IBLOK),ALINE160(LX4,IBLOK),ALINE161(LX4,IBLOK) &
  & ,ALINE162(LX4,IBLOK),ALINE163(LX4,IBLOK),ALINE164(LX4,IBLOK) &
  & ,ALINE165(LX4,IBLOK),ALINE166(LX4,IBLOK),ALINE167(LX4,IBLOK) &
  & ,ALINE168(LX4,IBLOK),ALINE169(LX4,IBLOK),ALINE170(LX4,IBLOK) &
  & ,ALINE171(LX4,IBLOK),ALINE172(LX4,IBLOK),ALINE173(LX4,IBLOK) &
  & ,ALINE174(LX4,IBLOK),ALINE175(LX4,IBLOK),ALINE176(LX4,IBLOK) &
  & ,ALINE177(LX4,IBLOK),ALINE178(LX4,IBLOK),ALINE179(LX4,IBLOK) &
  & ,ALINE180(LX4,IBLOK),ALINE181(LX4,IBLOK),ALINE182(LX4,IBLOK) &
  & ,ALINE183(LX4,IBLOK),ALINE184(LX4,IBLOK),ALINE185(LX4,IBLOK) &
  & ,ALINE186(LX4,IBLOK),ALINE187(LX4,IBLOK),ALINE188(LX4,IBLOK) &
  & ,ALINE189(LX4,IBLOK),ALINE190(LX4,IBLOK),ALINE191(LX4,IBLOK) &
  & ,ALINE192(LX4,IBLOK),ALINE193(LX4,IBLOK),ALINE194(LX4,IBLOK) &
  & ,ALINE195(LX4,IBLOK),ALINE196(LX4,IBLOK),ALINE197(LX4,IBLOK) &
  & ,ALINE198(LX4,IBLOK),ALINE199(LX4,IBLOK),ALINE200(LX4,IBLOK) &
  & ,ALINE201(LX4,IBLOK),ALINE202(LX4,IBLOK),ALINE203(LX4,IBLOK) &
  & ,ALINE204(LX4,IBLOK),ALINE205(LX4,IBLOK),ALINE206(LX4,IBLOK) &
  & ,ALINE207(LX4,IBLOK),ALINE208(LX4,IBLOK),ALINE209(LX4,IBLOK) &
  & ,ALINE210(LX4,IBLOK),ALINE211(LX4,IBLOK),ALINE212(LX4,IBLOK) &
  & ,ALINE213(LX4,IBLOK),ALINE214(LX4,IBLOK),ALINE215(LX4,IBLOK) &
  & ,ALINE216(LX4,IBLOK),ALINE217(LX4,IBLOK),ALINE218(LX4,IBLOK) &
  & ,ALINE219(LX4,IBLOK),ALINE220(LX4,IBLOK),ALINE221(LX4,IBLOK) &
  & ,ALINE222(LX4,IBLOK),ALINE223(LX4,IBLOK),ALINE224(LX4,IBLOK) &
  & ,ALINE225(LX4,IBLOK),ALINE226(LX4,IBLOK),ALINE227(LX4,IBLOK) &
  & ,ALINE228(LX4,IBLOK),ALINE229(LX4,IBLOK),ALINE230(LX4,IBLOK) &
  & ,ALINE231(LX4,IBLOK),ALINE232(LX4,IBLOK),ALINE233(LX4,IBLOK) &
  & ,ALINE234(LX4,IBLOK),ALINE235(LX4,IBLOK)
!********************************************************
COMMON/LINESMOM/ALINEMOM2(LX4,IBLOK,2),ALINEMOM3(LX4,IBLOK,2) &
  & ,ALINEMOM4(LX4,IBLOK,2),ALINEMOM5(LX4,IBLOK,2) &
  & ,ALINEMOM6(LX4,IBLOK,2),ALINEMOM7(LX4,IBLOK,2) &
  & ,ALINEMOM8(LX4,IBLOK,2),ALINEMOM9(LX4,IBLOK,2) &
  & ,ALINEMOM10(LX4,IBLOK,2),ALINEMOM11(LX4,IBLOK,2) &
  & ,ALINEMOM12(LX4,IBLOK,2),ALINEMOM13(LX4,IBLOK,2) &
  & ,ALINEMOM14(LX4,IBLOK,2),ALINEMOM15(LX4,IBLOK,2) &
  & ,ALINEMOM16(LX4,IBLOK,2),ALINEMOM17(LX4,IBLOK,2) &
  & ,ALINEMOM18(LX4,IBLOK,2),ALINEMOM19(LX4,IBLOK,2) &
  & ,ALINEMOM20(LX4,IBLOK,2),ALINEMOM21(LX4,IBLOK,2) &
  & ,ALINEMOM22(LX4,IBLOK,2),ALINEMOM23(LX4,IBLOK,2) &
  & ,ALINEMOM24(LX4,IBLOK,2),ALINEMOM25(LX4,IBLOK,2) &
  & ,ALINEMOM26(LX4,IBLOK,2),ALINEMOM27(LX4,IBLOK,2) &
  & ,ALINEMOM28(LX4,IBLOK,2),ALINEMOM29(LX4,IBLOK,2) &
  & ,ALINEMOM30(LX4,IBLOK,2),ALINEMOM31(LX4,IBLOK,2) &
  & ,ALINEMOM32(LX4,IBLOK,2),ALINEMOM33(LX4,IBLOK,2) &
  & ,ALINEMOM34(LX4,IBLOK,2),ALINEMOM35(LX4,IBLOK,2) &
  & ,ALINEMOM36(LX4,IBLOK,2),ALINEMOM37(LX4,IBLOK,2) &
  & ,ALINEMOM38(LX4,IBLOK,2),ALINEMOM39(LX4,IBLOK,2) &
  & ,ALINEMOM40(LX4,IBLOK,2),ALINEMOM41(LX4,IBLOK,2) &
  & ,ALINEMOM42(LX4,IBLOK,2),ALINEMOM43(LX4,IBLOK,2) &
  & ,ALINEMOM44(LX4,IBLOK,2),ALINEMOM45(LX4,IBLOK,2) &
  & ,ALINEMOM46(LX4,IBLOK,2),ALINEMOM47(LX4,IBLOK,2) &
  & ,ALINEMOM48(LX4,IBLOK,2),ALINEMOM49(LX4,IBLOK,2) &
  & ,ALINEMOM50(LX4,IBLOK,2),ALINEMOM51(LX4,IBLOK,2) &
  & ,ALINEMOM52(LX4,IBLOK,2),ALINEMOM53(LX4,IBLOK,2) &
  & ,ALINEMOM54(LX4,IBLOK,2),ALINEMOM55(LX4,IBLOK,2) &
  & ,ALINEMOM56(LX4,IBLOK,2),ALINEMOM57(LX4,IBLOK,2) &
  & ,ALINEMOM58(LX4,IBLOK,2),ALINEMOM59(LX4,IBLOK,2) &
  & ,ALINEMOM60(LX4,IBLOK,2),ALINEMOM61(LX4,IBLOK,2) &
  & ,ALINEMOM62(LX4,IBLOK,2),ALINEMOM63(LX4,IBLOK,2) &
  & ,ALINEMOM64(LX4,IBLOK,2) &
  & ,ALINEMOM65(LX4,IBLOK,2),ALINEMOM66(LX4,IBLOK,2) &
  & ,ALINEMOM67(LX4,IBLOK,2),ALINEMOM68(LX4,IBLOK,2) &
  & ,ALINEMOM69(LX4,IBLOK,2),ALINEMOM70(LX4,IBLOK,2) &
  & ,ALINEMOM71(LX4,IBLOK,2),ALINEMOM72(LX4,IBLOK,2) &
  & ,ALINEMOM73(LX4,IBLOK,2),ALINEMOM74(LX4,IBLOK,2) &
  & ,ALINEMOM75(LX4,IBLOK,2),ALINEMOM76(LX4,IBLOK,2) &
  & ,ALINEMOM77(LX4,IBLOK,2),ALINEMOM78(LX4,IBLOK,2) &
  & ,ALINEMOM79(LX4,IBLOK,2),ALINEMOM80(LX4,IBLOK,2) &
  & ,ALINEMOM81(LX4,IBLOK,2),ALINEMOM82(LX4,IBLOK,2) &
  & ,ALINEMOM83(LX4,IBLOK,2),ALINEMOM84(LX4,IBLOK,2) &
  & ,ALINEMOM85(LX4,IBLOK,2),ALINEMOM86(LX4,IBLOK,2) &
  & ,ALINEMOM87(LX4,IBLOK,2),ALINEMOM88(LX4,IBLOK,2) &
  & ,ALINEMOM89(LX4,IBLOK,2),ALINEMOM90(LX4,IBLOK,2) &
  & ,ALINEMOM91(LX4,IBLOK,2),ALINEMOM92(LX4,IBLOK,2) &
  & ,ALINEMOM93(LX4,IBLOK,2),ALINEMOM94(LX4,IBLOK,2) &
  & ,ALINEMOM95(LX4,IBLOK,2),ALINEMOM96(LX4,IBLOK,2) &
  & ,ALINEMOM97(LX4,IBLOK,2),ALINEMOM98(LX4,IBLOK,2) &
  & ,ALINEMOM99(LX4,IBLOK,2),ALINEMOM100(LX4,IBLOK,2) &
  & ,ALINEMOM101(LX4,IBLOK,2),ALINEMOM102(LX4,IBLOK,2) &
  & ,ALINEMOM103(LX4,IBLOK,2),ALINEMOM104(LX4,IBLOK,2) &
  & ,ALINEMOM105(LX4,IBLOK,2),ALINEMOM106(LX4,IBLOK,2) &
  & ,ALINEMOM107(LX4,IBLOK,2),ALINEMOM108(LX4,IBLOK,2) &
  & ,ALINEMOM109(LX4,IBLOK,2),ALINEMOM110(LX4,IBLOK,2) &
  & ,ALINEMOM111(LX4,IBLOK,2),ALINEMOM112(LX4,IBLOK,2) &
  & ,ALINEMOM113(LX4,IBLOK,2),ALINEMOM114(LX4,IBLOK,2) &
  & ,ALINEMOM115(LX4,IBLOK,2),ALINEMOM116(LX4,IBLOK,2) &
  & ,ALINEMOM117(LX4,IBLOK,2),ALINEMOM118(LX4,IBLOK,2) &
  & ,ALINEMOM119(LX4,IBLOK,2),ALINEMOM120(LX4,IBLOK,2) &
  & ,ALINEMOM121(LX4,IBLOK,2),ALINEMOM122(LX4,IBLOK,2) &
  & ,ALINEMOM123(LX4,IBLOK,2),ALINEMOM124(LX4,IBLOK,2) &
  & ,ALINEMOM125(LX4,IBLOK,2),ALINEMOM126(LX4,IBLOK,2) &
  & ,ALINEMOM127(LX4,IBLOK,2),ALINEMOM128(LX4,IBLOK,2) &
  & ,ALINEMOM129(LX4,IBLOK,2),ALINEMOM130(LX4,IBLOK,2) &
  & ,ALINEMOM131(LX4,IBLOK,2),ALINEMOM132(LX4,IBLOK,2) &
  & ,ALINEMOM133(LX4,IBLOK,2),ALINEMOM134(LX4,IBLOK,2) &
  & ,ALINEMOM135(LX4,IBLOK,2),ALINEMOM136(LX4,IBLOK,2) &
  & ,ALINEMOM137(LX4,IBLOK,2),ALINEMOM138(LX4,IBLOK,2) &
  & ,ALINEMOM139(LX4,IBLOK,2),ALINEMOM140(LX4,IBLOK,2) &
  & ,ALINEMOM141(LX4,IBLOK,2),ALINEMOM142(LX4,IBLOK,2) &
  & ,ALINEMOM143(LX4,IBLOK,2),ALINEMOM144(LX4,IBLOK,2) &
  & ,ALINEMOM145(LX4,IBLOK,2),ALINEMOM146(LX4,IBLOK,2) &
  & ,ALINEMOM147(LX4,IBLOK,2),ALINEMOM148(LX4,IBLOK,2) &
  & ,ALINEMOM149(LX4,IBLOK,2),ALINEMOM150(LX4,IBLOK,2) &
  & ,ALINEMOM151(LX4,IBLOK,2),ALINEMOM152(LX4,IBLOK,2) &
  & ,ALINEMOM153(LX4,IBLOK,2),ALINEMOM154(LX4,IBLOK,2) &
  & ,ALINEMOM155(LX4,IBLOK,2),ALINEMOM156(LX4,IBLOK,2) &
  & ,ALINEMOM157(LX4,IBLOK,2),ALINEMOM158(LX4,IBLOK,2) &
  & ,ALINEMOM159(LX4,IBLOK,2),ALINEMOM160(LX4,IBLOK,2) &
  & ,ALINEMOM161(LX4,IBLOK,2),ALINEMOM162(LX4,IBLOK,2) &
  & ,ALINEMOM163(LX4,IBLOK,2),ALINEMOM164(LX4,IBLOK,2) &
  & ,ALINEMOM165(LX4,IBLOK,2),ALINEMOM166(LX4,IBLOK,2) &
  & ,ALINEMOM167(LX4,IBLOK,2),ALINEMOM168(LX4,IBLOK,2) &
  & ,ALINEMOM169(LX4,IBLOK,2),ALINEMOM170(LX4,IBLOK,2) &
  & ,ALINEMOM171(LX4,IBLOK,2),ALINEMOM172(LX4,IBLOK,2) &
  & ,ALINEMOM173(LX4,IBLOK,2),ALINEMOM174(LX4,IBLOK,2) &
  & ,ALINEMOM175(LX4,IBLOK,2),ALINEMOM176(LX4,IBLOK,2) &
  & ,ALINEMOM177(LX4,IBLOK,2),ALINEMOM178(LX4,IBLOK,2) &
  & ,ALINEMOM179(LX4,IBLOK,2),ALINEMOM180(LX4,IBLOK,2) &
  & ,ALINEMOM181(LX4,IBLOK,2),ALINEMOM182(LX4,IBLOK,2) &
  & ,ALINEMOM183(LX4,IBLOK,2),ALINEMOM184(LX4,IBLOK,2) &
  & ,ALINEMOM185(LX4,IBLOK,2),ALINEMOM186(LX4,IBLOK,2) &
  & ,ALINEMOM187(LX4,IBLOK,2),ALINEMOM188(LX4,IBLOK,2) &
  & ,ALINEMOM189(LX4,IBLOK,2),ALINEMOM190(LX4,IBLOK,2) &
  & ,ALINEMOM191(LX4,IBLOK,2),ALINEMOM192(LX4,IBLOK,2) &
  & ,ALINEMOM193(LX4,IBLOK,2),ALINEMOM194(LX4,IBLOK,2) &
  & ,ALINEMOM195(LX4,IBLOK,2),ALINEMOM196(LX4,IBLOK,2) &
  & ,ALINEMOM197(LX4,IBLOK,2),ALINEMOM198(LX4,IBLOK,2) &
  & ,ALINEMOM199(LX4,IBLOK,2),ALINEMOM200(LX4,IBLOK,2) &
  & ,ALINEMOM201(LX4,IBLOK,2),ALINEMOM202(LX4,IBLOK,2) &
  & ,ALINEMOM203(LX4,IBLOK,2),ALINEMOM204(LX4,IBLOK,2) &
  & ,ALINEMOM205(LX4,IBLOK,2),ALINEMOM206(LX4,IBLOK,2) &
  & ,ALINEMOM207(LX4,IBLOK,2),ALINEMOM208(LX4,IBLOK,2) &
  & ,ALINEMOM209(LX4,IBLOK,2),ALINEMOM210(LX4,IBLOK,2) &
  & ,ALINEMOM211(LX4,IBLOK,2),ALINEMOM212(LX4,IBLOK,2) &
  & ,ALINEMOM213(LX4,IBLOK,2),ALINEMOM214(LX4,IBLOK,2) &
  & ,ALINEMOM215(LX4,IBLOK,2),ALINEMOM216(LX4,IBLOK,2) &
  & ,ALINEMOM217(LX4,IBLOK,2),ALINEMOM218(LX4,IBLOK,2) &
  & ,ALINEMOM219(LX4,IBLOK,2),ALINEMOM220(LX4,IBLOK,2) &
  & ,ALINEMOM221(LX4,IBLOK,2),ALINEMOM222(LX4,IBLOK,2) &
  & ,ALINEMOM223(LX4,IBLOK,2),ALINEMOM224(LX4,IBLOK,2) &
  & ,ALINEMOM225(LX4,IBLOK,2),ALINEMOM226(LX4,IBLOK,2) &
  & ,ALINEMOM227(LX4,IBLOK,2),ALINEMOM228(LX4,IBLOK,2) &
  & ,ALINEMOM229(LX4,IBLOK,2),ALINEMOM230(LX4,IBLOK,2) &
  & ,ALINEMOM231(LX4,IBLOK,2),ALINEMOM232(LX4,IBLOK,2) &
  & ,ALINEMOM233(LX4,IBLOK,2),ALINEMOM234(LX4,IBLOK,2) &
  & ,ALINEMOM235(LX4,IBLOK,2)
!***************************************************************
COMMON/ITEM/ITR,IBIN,ITOT
!***********************************************************
complex(real64) ALLP
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
complex(real64) ALLPJ0PP, ALLPJ0PM, ALLPJ0MP, ALLPJ0MM
complex(real64) ALLPJ2PP, ALLPJ2PM, ALLPJ2MP, ALLPJ2MM
complex(real64) ALLPJ1P, ALLPJ1M
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
complex(real64) ALLPJ0, ALLPJ1, ALLPJ2, ALLPJALL
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
complex(real64) ALINE1,ALINE2,ALINE3,ALINE4,ALINE5,ALINE6,ALINE7
complex(real64) ALINE8,ALINE9,ALINE10,ALINE11,ALINE12,ALINE13
complex(real64) ALINE14,ALINE15,ALINE16,ALINE17,ALINE18,ALINE19
complex(real64) ALINE20,ALINE21,ALINE22,ALINE23,ALINE24,ALINE25
complex(real64) ALINE26,ALINE27,ALINE28,ALINE29,ALINE30
complex(real64) ALINE31,ALINE32,ALINE33,ALINE34,ALINE35
complex(real64) ALINE36,ALINE37,ALINE38,ALINE39,ALINE40
complex(real64) ALINE41,ALINE42,ALINE43,ALINE44,ALINE45
complex(real64) ALINE46,ALINE47,ALINE48,ALINE49,ALINE50
complex(real64) ALINE51,ALINE52,ALINE53,ALINE54,ALINE55
complex(real64) ALINE56,ALINE57,ALINE58,ALINE59,ALINE60
complex(real64) ALINE61,ALINE62,ALINE63,ALINE64,ALINE65
complex(real64) ALINE66,ALINE67,ALINE68,ALINE69,ALINE70
complex(real64) ALINE71,ALINE72,ALINE73,ALINE74,ALINE75
complex(real64) ALINE76,ALINE77,ALINE78,ALINE79,ALINE80
complex(real64) ALINE81,ALINE82,ALINE83,ALINE84,ALINE85
complex(real64) ALINE86,ALINE87,ALINE88,ALINE89,ALINE90
complex(real64) ALINE91,ALINE92,ALINE93,ALINE94,ALINE95
complex(real64) ALINE96,ALINE97,ALINE98,ALINE99,ALINE100
complex(real64) ALINE101,ALINE102,ALINE103,ALINE104,ALINE105
complex(real64) ALINE106,ALINE107,ALINE108,ALINE109,ALINE110
complex(real64) ALINE111,ALINE112,ALINE113,ALINE114,ALINE115
complex(real64) ALINE116,ALINE117,ALINE118,ALINE119,ALINE120
complex(real64) ALINE121,ALINE122,ALINE123,ALINE124,ALINE125
complex(real64) ALINE126,ALINE127,ALINE128,ALINE129,ALINE130
complex(real64) ALINE131,ALINE132,ALINE133,ALINE134,ALINE135
complex(real64) ALINE136,ALINE137,ALINE138,ALINE139,ALINE140
complex(real64) ALINE141,ALINE142,ALINE143,ALINE144,ALINE145
complex(real64) ALINE146,ALINE147,ALINE148,ALINE149,ALINE150
complex(real64) ALINE151,ALINE152,ALINE153,ALINE154,ALINE155
complex(real64) ALINE156,ALINE157,ALINE158,ALINE159,ALINE160
complex(real64) ALINE161,ALINE162,ALINE163,ALINE164,ALINE165
complex(real64) ALINE166,ALINE167,ALINE168,ALINE169,ALINE170
complex(real64) ALINE171,ALINE172,ALINE173,ALINE174,ALINE175
complex(real64) ALINE176,ALINE177,ALINE178,ALINE179,ALINE180
complex(real64) ALINE181,ALINE182,ALINE183,ALINE184,ALINE185
complex(real64) ALINE186,ALINE187,ALINE188,ALINE189,ALINE190
complex(real64) ALINE191,ALINE192,ALINE193,ALINE194,ALINE195
complex(real64) ALINE196,ALINE197,ALINE198,ALINE199,ALINE200
complex(real64) ALINE201,ALINE202,ALINE203,ALINE204,ALINE205
complex(real64) ALINE206,ALINE207,ALINE208,ALINE209,ALINE210
complex(real64) ALINE211,ALINE212,ALINE213,ALINE214,ALINE215
complex(real64) ALINE216,ALINE217,ALINE218,ALINE219,ALINE220
complex(real64) ALINE221,ALINE222,ALINE223,ALINE224,ALINE225
complex(real64) ALINE226,ALINE227,ALINE228,ALINE229,ALINE230
complex(real64) ALINE231,ALINE232,ALINE233,ALINE234,ALINE235
!***********************************************************
complex(real64) ALLPMOM
complex(real64) ALLPMOMJ0P, ALLPMOMJ0M
complex(real64) ALLPMOMJ1
complex(real64) ALLPMOMJ2P, ALLPMOMJ2M
complex(real64) ALLPMOMJ0, ALLPMOMJ2
complex(real64) ALLPMOMJALL
!***********************************************************
complex(real64) ALINEMOM2,ALINEMOM3,ALINEMOM4
complex(real64) ALINEMOM5,ALINEMOM6,ALINEMOM7,ALINEMOM8
complex(real64) ALINEMOM9,ALINEMOM10,ALINEMOM11,ALINEMOM12
complex(real64) ALINEMOM13,ALINEMOM14,ALINEMOM15,ALINEMOM16
complex(real64) ALINEMOM17,ALINEMOM18,ALINEMOM19,ALINEMOM20
complex(real64) ALINEMOM21,ALINEMOM22,ALINEMOM23,ALINEMOM24
complex(real64) ALINEMOM25,ALINEMOM26,ALINEMOM27,ALINEMOM28
complex(real64) ALINEMOM29,ALINEMOM30,ALINEMOM31,ALINEMOM32
complex(real64) ALINEMOM33,ALINEMOM34,ALINEMOM35,ALINEMOM36
complex(real64) ALINEMOM37,ALINEMOM38,ALINEMOM39,ALINEMOM40
complex(real64) ALINEMOM41,ALINEMOM42,ALINEMOM43,ALINEMOM44
complex(real64) ALINEMOM45,ALINEMOM46,ALINEMOM47,ALINEMOM48
complex(real64) ALINEMOM49,ALINEMOM50,ALINEMOM51,ALINEMOM52
complex(real64) ALINEMOM53,ALINEMOM54,ALINEMOM55,ALINEMOM56
complex(real64) ALINEMOM57,ALINEMOM58,ALINEMOM59,ALINEMOM60
complex(real64) ALINEMOM61,ALINEMOM62,ALINEMOM63,ALINEMOM64
complex(real64) ALINEMOM65,ALINEMOM66,ALINEMOM67,ALINEMOM68
complex(real64) ALINEMOM69,ALINEMOM70,ALINEMOM71,ALINEMOM72
complex(real64) ALINEMOM73,ALINEMOM74,ALINEMOM75,ALINEMOM76
complex(real64) ALINEMOM77,ALINEMOM78,ALINEMOM79,ALINEMOM80
complex(real64) ALINEMOM81,ALINEMOM82,ALINEMOM83,ALINEMOM84
complex(real64) ALINEMOM85,ALINEMOM86,ALINEMOM87,ALINEMOM88
complex(real64) ALINEMOM89,ALINEMOM90,ALINEMOM91,ALINEMOM92
complex(real64) ALINEMOM93,ALINEMOM94,ALINEMOM95,ALINEMOM96
complex(real64) ALINEMOM97,ALINEMOM98,ALINEMOM99,ALINEMOM100
complex(real64) ALINEMOM101,ALINEMOM102,ALINEMOM103,ALINEMOM104
complex(real64) ALINEMOM105,ALINEMOM106,ALINEMOM107,ALINEMOM108
complex(real64) ALINEMOM109,ALINEMOM110
complex(real64) ALINEMOM111,ALINEMOM112,ALINEMOM113,ALINEMOM114
complex(real64) ALINEMOM115,ALINEMOM116,ALINEMOM117,ALINEMOM118
complex(real64) ALINEMOM119,ALINEMOM120,ALINEMOM121,ALINEMOM122
complex(real64) ALINEMOM123,ALINEMOM124,ALINEMOM125,ALINEMOM126
complex(real64) ALINEMOM127,ALINEMOM128,ALINEMOM129,ALINEMOM130
complex(real64) ALINEMOM131,ALINEMOM132,ALINEMOM133,ALINEMOM134
complex(real64) ALINEMOM135,ALINEMOM136,ALINEMOM137,ALINEMOM138
complex(real64) ALINEMOM139,ALINEMOM140,ALINEMOM141,ALINEMOM142
complex(real64) ALINEMOM143,ALINEMOM144,ALINEMOM145,ALINEMOM146
complex(real64) ALINEMOM147,ALINEMOM148,ALINEMOM149,ALINEMOM150
complex(real64) ALINEMOM151,ALINEMOM152,ALINEMOM153,ALINEMOM154
complex(real64) ALINEMOM155,ALINEMOM156,ALINEMOM157,ALINEMOM158
complex(real64) ALINEMOM159,ALINEMOM160,ALINEMOM161,ALINEMOM162
complex(real64) ALINEMOM163,ALINEMOM164,ALINEMOM165,ALINEMOM166
complex(real64) ALINEMOM167,ALINEMOM168,ALINEMOM169,ALINEMOM170
complex(real64) ALINEMOM171,ALINEMOM172,ALINEMOM173,ALINEMOM174
complex(real64) ALINEMOM175,ALINEMOM176,ALINEMOM177,ALINEMOM178
complex(real64) ALINEMOM179,ALINEMOM180,ALINEMOM181,ALINEMOM182
complex(real64) ALINEMOM183,ALINEMOM184,ALINEMOM185,ALINEMOM186
complex(real64) ALINEMOM187,ALINEMOM188,ALINEMOM189,ALINEMOM190
complex(real64) ALINEMOM191,ALINEMOM192,ALINEMOM193,ALINEMOM194
complex(real64) ALINEMOM195,ALINEMOM196,ALINEMOM197,ALINEMOM198
complex(real64) ALINEMOM199,ALINEMOM200,ALINEMOM201,ALINEMOM202
complex(real64) ALINEMOM203,ALINEMOM204,ALINEMOM205,ALINEMOM206
complex(real64) ALINEMOM207,ALINEMOM208,ALINEMOM209,ALINEMOM210
complex(real64) ALINEMOM211,ALINEMOM212,ALINEMOM213,ALINEMOM214
complex(real64) ALINEMOM215,ALINEMOM216,ALINEMOM217,ALINEMOM218
complex(real64) ALINEMOM219,ALINEMOM220,ALINEMOM221,ALINEMOM222
complex(real64) ALINEMOM223,ALINEMOM224,ALINEMOM225,ALINEMOM226
complex(real64) ALINEMOM227,ALINEMOM228,ALINEMOM229,ALINEMOM230
complex(real64) ALINEMOM231,ALINEMOM232,ALINEMOM233,ALINEMOM234
complex(real64) ALINEMOM235
!***********************************************************
complex(real64) ACORLP
complex(real64) AVACLP
!***********************************************************
complex(real64) ACORLPMOM
complex(real64) AVACLPMOM
!***********************************************************
complex(real64) ACORLPMOMJ0P, ACORLPMOMJ0M
complex(real64) AVACLMOMJ0P, AVACLMOMJ0M
!***********************************************************
complex(real64) ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
complex(real64) AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
complex(real64) ACORLPMOMJALL
complex(real64) AVACLMOMALL
!***********************************************************
complex(real64) ACORLPMOMJ2P, ACORLPMOMJ2M
complex(real64) AVACLMOMJ2P, AVACLMOMJ2M
!***********************************************************
complex(real64) ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
complex(real64) AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
!***********************************************************
complex(real64) ACORLJ1P,ACORLJ1M
complex(real64) AVACLJ1P,AVACLJ1M
!***********************************************************
complex(real64) ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
complex(real64) AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
!***********************************************************
complex(real64) ACORLJ0, ACORLJ1, ACORLJ2, ACORLALL
complex(real64) AVACLJ0, AVACLJ1, AVACLJ2, AVACLALL
!***********************************************************      
complex(real64) ALINEP(LX4,NTOTAL)
complex(real64) ALINEPMOM(2,LX4,NTOTALMOM)
!***********************************************************
!     NEW CORRELATORS
!***********************************************************
complex(real64) ALINEJ0PP(LX4,NOPJ0PP)
complex(real64) ALINEJ0PM(LX4,NOPJ0PM)
complex(real64) ALINEJ0MP(LX4,NOPJ0MP)
complex(real64) ALINEJ0MM(LX4,NOPJ0MM)
complex(real64) ALINEJ2PP(LX4,NOPJ2PP)
complex(real64) ALINEJ2PM(LX4,NOPJ2PM)
complex(real64) ALINEJ2MP(LX4,NOPJ2MP)
complex(real64) ALINEJ2MM(LX4,NOPJ2MM)
complex(real64) ALINEJ1P(LX4,NOPJ1P)
complex(real64) ALINEJ1M(LX4,NOPJ1M)
!***********************************************************
complex(real64) ALINEJ0(LX4,NOPFULJ0)
complex(real64) ALINEJ2(LX4,NOPFULJ2)
complex(real64) ALINEJ1(LX4,NOPFULJ1)
complex(real64) ALINEJALL(LX4,NTOTAL)
!***********************************************************
complex(real64) ALINEMOMJ0P(2,LX4,NOPJ0PMOM)
complex(real64) ALINEMOMJ0M(2,LX4,NOPJ0MMOM)
complex(real64) ALINEMOMJ2P(2,LX4,NOPJ2PMOM)
complex(real64) ALINEMOMJ2M(2,LX4,NOPJ2MMOM)
complex(real64) ALINEMOMJ1(2,LX4,NOPJ1MOM)
!***********************************************************
complex(real64) ALINEMOMJ0(2,LX4,NOPFULJ0MOM)
complex(real64) ALINEMOMJ2(2,LX4,NOPFULJ2MOM)
!***********************************************************
complex(real64) ALINEMOMJALL(2,LX4,NTOTALMOM)
!************************************************************
JBIN=(ITR-1)/IBIN+1
!************************************************************
do N4=1,LX4
!
do IDL=1,IBLOK
ALINEP(N4,IDL)=ALINE1(N4,IDL)
ALINEJ0PP(N4,IDL)=ALINE1(N4,IDL)
ALINEJ0PM(N4,IDL)=ALINE87(N4,IDL)
ALINEJ0MP(N4,IDL)=ALINE63(N4,IDL)
ALINEJ0MM(N4,IDL)=ALINE27(N4,IDL)
!
ALINEJ1P(N4,IDL)=ALINE3(N4,IDL)
ALINEJ1M(N4,IDL)=ALINE9(N4,IDL)
!
ALINEJ2PP(N4,IDL)=ALINE4(N4,IDL)
ALINEJ2PM(N4,IDL)=ALINE30(N4,IDL)
ALINEJ2MP(N4,IDL)=ALINE31(N4,IDL)
ALINEJ2MM(N4,IDL)=ALINE97(N4,IDL)
!            
3 CONTINUE
end do
!
ID=IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE2(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE2(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE147(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE81(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE33(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE6(N4,IDL)
ALINEJ1M(N4,ID)=ALINE12(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE7(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE36(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE37(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE155(N4,IDL)
!            
4 CONTINUE
end do
!
ID=2*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE3(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE5(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE157(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE88(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE39(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE18(N4,IDL)
ALINEJ1M(N4,ID)=ALINE15(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE10(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE42(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE43(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE165(N4,IDL)
!
5 CONTINUE
end do
!
ID=3*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE4(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE8(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE167(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE123(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE45(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE21(N4,IDL)
ALINEJ1M(N4,ID)=ALINE24(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE13(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE48(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE49(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE175(N4,IDL)
!            
6 CONTINUE
end do
!
ID=4*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE5(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE11(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE177(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE129(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE57(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE28(N4,IDL)
ALINEJ1M(N4,ID)=ALINE29(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE16(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE60(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE61(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE185(N4,IDL)
!
7 CONTINUE
end do
!
ID=5*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE6(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE14(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE187(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE135(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE69(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE34(N4,IDL)
ALINEJ1M(N4,ID)=ALINE35(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE19(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE95(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE67(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE195(N4,IDL)
!
8 CONTINUE
end do
!
ID=6*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE7(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE17(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE197(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE148(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE75(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE40(N4,IDL)
ALINEJ1M(N4,ID)=ALINE41(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE22(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE102(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE85(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE205(N4,IDL)
!
9 CONTINUE
end do
!
ID=7*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE8(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE20(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE207(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE158(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE89(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE46(N4,IDL)
ALINEJ1M(N4,ID)=ALINE47(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE25(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE108(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE96(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE215(N4,IDL)
!
10 CONTINUE
end do
!     
ID=8*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE9(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE23(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE217(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE168(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE99(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE51(N4,IDL)
ALINEJ1M(N4,ID)=ALINE54(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE52(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE114(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE103(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE225(N4,IDL)
!
11 CONTINUE
end do
!     
ID=9*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE10(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE26(N4,IDL)
ALINEJ0PM(N4,ID)=ALINE227(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE178(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE105(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE58(N4,IDL)
ALINEJ1M(N4,ID)=ALINE59(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE55(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE120(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE109(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE235(N4,IDL)
!
12 CONTINUE
end do
!     
ID=10*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE11(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE32(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE188(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE111(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE64(N4,IDL)
ALINEJ1M(N4,ID)=ALINE71(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE66(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE144(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE115(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE73(N4,IDL)
!
13 CONTINUE
end do
!
ID=11*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE12(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE38(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE198(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE117(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE65(N4,IDL)
ALINEJ1M(N4,ID)=ALINE77(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE72(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE153(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE121(N4,IDL)
ALINEJ2MM(N4,ID)=ALINE79(N4,IDL)
!
14 CONTINUE
end do
!
ID=12*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE13(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE44(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE208(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE141(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE70(N4,IDL)
ALINEJ1M(N4,ID)=ALINE82(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE78(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE163(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE127(N4,IDL)
!
15 CONTINUE
end do
!
ID=13*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE14(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE50(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE218(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE149(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE76(N4,IDL)
ALINEJ1M(N4,ID)=ALINE83(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE84(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE173(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE133(N4,IDL)
!
16 CONTINUE
end do
!
ID=14*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE15(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE53(N4,IDL)
ALINEJ0MP(N4,ID)=ALINE228(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE159(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE90(N4,IDL)
ALINEJ1M(N4,ID)=ALINE92(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE94(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE183(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE139(N4,IDL)
!
17 CONTINUE
end do
!
ID=15*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE16(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE56(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE169(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE91(N4,IDL)
ALINEJ1M(N4,ID)=ALINE93(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE126(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE193(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE145(N4,IDL)
!
18 CONTINUE
end do
!
ID=16*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE17(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE62(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE179(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE100(N4,IDL)
ALINEJ1M(N4,ID)=ALINE101(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE132(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE203(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE154(N4,IDL)
!
19 CONTINUE
end do
!
ID=17*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE18(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE68(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE189(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE106(N4,IDL)
ALINEJ1M(N4,ID)=ALINE107(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE138(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE213(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE164(N4,IDL)
!
20 CONTINUE
end do
!
ID=18*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE19(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE74(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE199(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE112(N4,IDL)
ALINEJ1M(N4,ID)=ALINE113(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE152(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE223(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE174(N4,IDL)
!  
21 CONTINUE
end do
!
ID=19*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE20(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE80(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE209(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE118(N4,IDL)
ALINEJ1M(N4,ID)=ALINE136(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE162(N4,IDL)
ALINEJ2PM(N4,ID)=ALINE233(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE184(N4,IDL)
!
22 CONTINUE
end do
!
ID=20*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE21(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE86(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE219(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE124(N4,IDL)
ALINEJ1M(N4,ID)=ALINE137(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE172(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE194(N4,IDL)
!            
23 CONTINUE
end do
!
ID=21*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE22(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE98(N4,IDL)
ALINEJ0MM(N4,ID)=ALINE229(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE125(N4,IDL)
ALINEJ1M(N4,ID)=ALINE143(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE182(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE204(N4,IDL)
!
24 CONTINUE
end do
!
ID=22*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE23(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE104(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE130(N4,IDL)
ALINEJ1M(N4,ID)=ALINE151(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE192(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE214(N4,IDL)
!
25 CONTINUE
end do
!
ID=23*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE24(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE110(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE131(N4,IDL)
ALINEJ1M(N4,ID)=ALINE161(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE202(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE224(N4,IDL)
!
26 CONTINUE
end do
!
ID=24*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE25(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE116(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE142(N4,IDL)
ALINEJ1M(N4,ID)=ALINE171(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE212(N4,IDL)
ALINEJ2MP(N4,ID)=ALINE234(N4,IDL)
!
27 CONTINUE
end do
!
ID=25*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE26(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE122(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE150(N4,IDL)
ALINEJ1M(N4,ID)=ALINE181(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE222(N4,IDL)
!
28 CONTINUE
end do
!
ID=26*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE27(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE128(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE160(N4,IDL)
ALINEJ1M(N4,ID)=ALINE191(N4,IDL)
!
ALINEJ2PP(N4,ID)=ALINE232(N4,IDL)
!
29 CONTINUE
end do
!
ID=27*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE28(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE134(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE170(N4,IDL)
ALINEJ1M(N4,ID)=ALINE201(N4,IDL)
!
31 CONTINUE
end do
!     
ID=28*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE29(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE140(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE180(N4,IDL)
ALINEJ1M(N4,ID)=ALINE211(N4,IDL)
!
32 CONTINUE
end do
!
ID=29*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE30(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE146(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE190(N4,IDL)
ALINEJ1M(N4,ID)=ALINE221(N4,IDL)
!
33 CONTINUE
end do
!
ID=30*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE31(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE156(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE200(N4,IDL)
ALINEJ1M(N4,ID)=ALINE231(N4,IDL)
!
34 CONTINUE
end do
!
ID=31*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE32(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE166(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE210(N4,IDL)
ALINEJ1M(N4,ID)=ALINE119(N4,IDL)
!     
35 CONTINUE
end do
!
ID=32*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE33(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE176(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE220(N4,IDL)
!
36 CONTINUE
end do
!
ID=33*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE34(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE186(N4,IDL)
!
ALINEJ1P(N4,ID)=ALINE230(N4,IDL)
!
37 CONTINUE
end do
!     
ID=34*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE35(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE196(N4,IDL)
!
!
38 CONTINUE
end do
!
ID=35*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE36(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE206(N4,IDL)
39 CONTINUE
end do
!     
ID=36*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE37(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE216(N4,IDL)
40 CONTINUE
end do
!
ID=37*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE38(N4,IDL)
ALINEJ0PP(N4,ID)=ALINE226(N4,IDL)
41 CONTINUE
end do
!
ID=38*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE39(N4,IDL)
42 CONTINUE
end do
!
ID=39*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE40(N4,IDL)
43 CONTINUE
end do
!
ID=40*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE41(N4,IDL)
44 CONTINUE
end do
!
ID=41*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE42(N4,IDL)
45 CONTINUE
end do
!
ID=42*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE43(N4,IDL)
46 CONTINUE
end do
!
ID=43*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE44(N4,IDL)
47 CONTINUE
end do
!
ID=44*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE45(N4,IDL)
48 CONTINUE
end do
!     
ID=45*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE46(N4,IDL)
49 CONTINUE
end do
!
ID=46*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE47(N4,IDL)
50 CONTINUE
end do
!
ID=47*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE48(N4,IDL)
51 CONTINUE
end do
!
ID=48*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE49(N4,IDL)
52 CONTINUE
end do
!
ID=49*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE50(N4,IDL)
53 CONTINUE
end do
!
ID=50*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE51(N4,IDL)
54 CONTINUE
end do
!
ID=51*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE52(N4,IDL)
55 CONTINUE
end do
!
ID=52*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE53(N4,IDL)
56 CONTINUE
end do
!
ID=53*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE54(N4,IDL)
57 CONTINUE
end do
!
ID=54*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE55(N4,IDL)
58 CONTINUE
end do
!
ID=55*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE56(N4,IDL)
59 CONTINUE
end do
!
ID=56*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE57(N4,IDL)
60 CONTINUE
end do
!
ID=57*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE58(N4,IDL)
61 CONTINUE
end do
!
ID=58*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE59(N4,IDL)
62 CONTINUE
end do
!
ID=59*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE60(N4,IDL)
63 CONTINUE
end do
!
ID=60*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE61(N4,IDL)
64 CONTINUE
end do
!
ID=61*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE62(N4,IDL)
65 CONTINUE
end do
!
ID=62*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE63(N4,IDL)
66 CONTINUE
end do
!
ID=63*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE64(N4,IDL)
67 CONTINUE
end do
!
ID=64*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE65(N4,IDL)
68 CONTINUE
end do
!
ID=65*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE66(N4,IDL)
69 CONTINUE
end do
!
ID=66*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE67(N4,IDL)
70 CONTINUE
end do
!
ID=67*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE68(N4,IDL)
71 CONTINUE
end do
!
ID=68*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE69(N4,IDL)
72 CONTINUE
end do
!
ID=69*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE70(N4,IDL)
73 CONTINUE
end do
!
ID=70*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE71(N4,IDL)
74 CONTINUE
end do
!
ID=71*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE72(N4,IDL)
75 CONTINUE
end do
!
ID=72*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE73(N4,IDL)
76 CONTINUE
end do
!
ID=73*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE74(N4,IDL)
77 CONTINUE
end do
!
ID=74*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE75(N4,IDL)
78 CONTINUE
end do
!
ID=75*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE76(N4,IDL)
79 CONTINUE
end do
!
ID=76*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE77(N4,IDL)
80 CONTINUE
end do
!
ID=77*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE78(N4,IDL)
81 CONTINUE
end do
!     
ID=78*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE79(N4,IDL)
82 CONTINUE
end do
!
ID=79*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE80(N4,IDL)
83 CONTINUE
end do
!
ID=80*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE81(N4,IDL)
84 CONTINUE
end do
!
ID=81*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE82(N4,IDL)
85 CONTINUE
end do
!
ID=82*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE83(N4,IDL)
86 CONTINUE
end do
!
ID=83*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE84(N4,IDL)
87 CONTINUE
end do
!
ID=84*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE85(N4,IDL)
88 CONTINUE
end do
!
ID=85*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE86(N4,IDL)
89 CONTINUE
end do
!
ID=86*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE87(N4,IDL)
90 CONTINUE
end do
!
ID=87*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE88(N4,IDL)
91 CONTINUE
end do
!
ID=88*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE89(N4,IDL)
92 CONTINUE
end do
!
ID=89*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE90(N4,IDL)
93 CONTINUE
end do
!
ID=90*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE91(N4,IDL)
94 CONTINUE
end do
!
ID=91*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE92(N4,IDL)
95 CONTINUE
end do
!
ID=92*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE93(N4,IDL)
96 CONTINUE
end do
!
ID=93*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE94(N4,IDL)
97 CONTINUE
end do
!
ID=94*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE95(N4,IDL)
98 CONTINUE
end do
!
ID=95*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE96(N4,IDL)
99 CONTINUE
end do
!
ID=96*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE97(N4,IDL)
300 CONTINUE
end do
!
ID=97*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE98(N4,IDL)
301 CONTINUE
end do
!
ID=98*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE99(N4,IDL)
302 CONTINUE
end do
!
ID=99*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE100(N4,IDL)
303 CONTINUE
end do
!
ID=100*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE101(N4,IDL)
304 CONTINUE
end do
!
ID=101*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE102(N4,IDL)
305 CONTINUE
end do
!
ID=102*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE103(N4,IDL)
306 CONTINUE
end do
!
ID=103*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE104(N4,IDL)
307 CONTINUE
end do
!
ID=104*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE105(N4,IDL)
308 CONTINUE
end do
!
ID=105*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE106(N4,IDL)
309 CONTINUE
end do
!
ID=106*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE107(N4,IDL)
310 CONTINUE
end do
!
ID=107*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE108(N4,IDL)
311 CONTINUE
end do
!
ID=108*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE109(N4,IDL)
312 CONTINUE
end do
!
ID=109*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE110(N4,IDL)
313 CONTINUE
end do
!
ID=110*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE111(N4,IDL)
314 CONTINUE
end do
!
ID=111*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE112(N4,IDL)
315 CONTINUE
end do
!
ID=112*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE113(N4,IDL)
316 CONTINUE
end do
!
ID=113*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE114(N4,IDL)
317 CONTINUE
end do
!
ID=114*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE115(N4,IDL)
318 CONTINUE
end do
!
ID=115*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE116(N4,IDL)
319 CONTINUE
end do
!
ID=116*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE117(N4,IDL)
320 CONTINUE
end do
!
ID=117*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE118(N4,IDL)
321 CONTINUE
end do
!     
ID=118*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE119(N4,IDL)
322 CONTINUE
end do
!
ID=119*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE120(N4,IDL)
323 CONTINUE
end do
!
ID=120*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE121(N4,IDL)
324 CONTINUE
end do
!
ID=121*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE122(N4,IDL)
325 CONTINUE
end do
!
ID=122*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE123(N4,IDL)
326 CONTINUE
end do
!
ID=123*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE124(N4,IDL)
327 CONTINUE
end do
!     
ID=124*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE125(N4,IDL)
328 CONTINUE
end do
!
ID=125*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE126(N4,IDL)
329 CONTINUE
end do
!     
ID=126*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE127(N4,IDL)
330 CONTINUE
end do
!
ID=127*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE128(N4,IDL)
331 CONTINUE
end do
!
ID=128*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE129(N4,IDL)
332 CONTINUE
end do
!
ID=129*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE130(N4,IDL)
333 CONTINUE
end do
!
ID=130*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE131(N4,IDL)
334 CONTINUE
end do
!
ID=131*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE132(N4,IDL)
335 CONTINUE
end do
!
ID=132*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE133(N4,IDL)
336 CONTINUE
end do
!
ID=133*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE134(N4,IDL)
337 CONTINUE
end do
!
ID=134*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE135(N4,IDL)
338 CONTINUE
end do
!     
ID=135*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE136(N4,IDL)
339 CONTINUE
end do
!
ID=136*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE137(N4,IDL)
340 CONTINUE
end do
!
ID=137*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE138(N4,IDL)
341 CONTINUE
end do
!
ID=138*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE139(N4,IDL)
342 CONTINUE
end do
!
ID=139*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE140(N4,IDL)
343 CONTINUE
end do
!
ID=140*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE141(N4,IDL)
344 CONTINUE
end do
!
ID=141*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE142(N4,IDL)
345 CONTINUE
end do
!
ID=142*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE143(N4,IDL)
346 CONTINUE
end do
!
ID=143*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE144(N4,IDL)
347 CONTINUE
end do
!
ID=144*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE145(N4,IDL)
348 CONTINUE
end do
!
ID=145*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE146(N4,IDL)
349 CONTINUE
end do
!
ID=146*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE147(N4,IDL)
350 CONTINUE
end do
!
ID=147*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE148(N4,IDL)
351 CONTINUE
end do
!
ID=148*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE149(N4,IDL)
352 CONTINUE
end do
!
ID=149*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE150(N4,IDL)
353 CONTINUE
end do
!
ID=150*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE151(N4,IDL)
354 CONTINUE
end do
!
ID=151*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE152(N4,IDL)
355 CONTINUE
end do
!
ID=152*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE153(N4,IDL)
356 CONTINUE
end do
!
ID=153*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE154(N4,IDL)
357 CONTINUE
end do
!
ID=154*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE155(N4,IDL)
358 CONTINUE
end do
!
ID=155*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE156(N4,IDL)
359 CONTINUE
end do
!
ID=156*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE157(N4,IDL)
360 CONTINUE
end do
!
ID=157*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE158(N4,IDL)
361 CONTINUE
end do
!
ID=158*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE159(N4,IDL)
362 CONTINUE
end do
!
ID=159*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE160(N4,IDL)
363 CONTINUE
end do
!
ID=160*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE161(N4,IDL)
364 CONTINUE
end do
!
ID=161*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE162(N4,IDL)
365 CONTINUE
end do
!
ID=162*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE163(N4,IDL)
366 CONTINUE
end do
!
ID=163*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE164(N4,IDL)
367 CONTINUE
end do
!
ID=164*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE165(N4,IDL)
368 CONTINUE
end do
!
ID=165*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE166(N4,IDL)
369 CONTINUE
end do
!
ID=166*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE167(N4,IDL)
370 CONTINUE
end do
!
ID=167*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE168(N4,IDL)
371 CONTINUE
end do
!
ID=168*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE169(N4,IDL)
372 CONTINUE
end do
!
ID=169*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE170(N4,IDL)
373 CONTINUE
end do
!
ID=170*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE171(N4,IDL)
374 CONTINUE
end do
!
ID=171*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE172(N4,IDL)
375 CONTINUE
end do
!
ID=172*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE173(N4,IDL)
376 CONTINUE
end do
!
ID=173*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE174(N4,IDL)
377 CONTINUE
end do
!
ID=174*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE175(N4,IDL)
378 CONTINUE
end do
!
ID=175*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE176(N4,IDL)
379 CONTINUE
end do
!
ID=176*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE177(N4,IDL)
380 CONTINUE
end do
!
ID=177*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE178(N4,IDL)
381 CONTINUE
end do
!     
ID=178*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE179(N4,IDL)
382 CONTINUE
end do
!
ID=179*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE180(N4,IDL)
383 CONTINUE
end do
!
ID=180*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE181(N4,IDL)
384 CONTINUE
end do
!
ID=181*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE182(N4,IDL)
385 CONTINUE
end do
!
ID=182*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE183(N4,IDL)
386 CONTINUE
end do
!
ID=183*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE184(N4,IDL)
387 CONTINUE
end do
!
ID=184*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE185(N4,IDL)
388 CONTINUE
end do
!
ID=185*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE186(N4,IDL)
389 CONTINUE
end do
!
ID=186*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE187(N4,IDL)
390 CONTINUE
end do
!
ID=187*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE188(N4,IDL)
391 CONTINUE
end do
!
ID=188*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE189(N4,IDL)
392 CONTINUE
end do
!
ID=189*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE190(N4,IDL)
393 CONTINUE
end do
!
ID=190*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE191(N4,IDL)
394 CONTINUE
end do
!
ID=191*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE192(N4,IDL)
395 CONTINUE
end do
!
ID=192*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE193(N4,IDL)
396 CONTINUE
end do
!
ID=193*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE194(N4,IDL)
397 CONTINUE
end do
!
ID=194*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE195(N4,IDL)
398 CONTINUE
end do
!
ID=195*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE196(N4,IDL)
399 CONTINUE
end do
!
ID=196*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE197(N4,IDL)
400 CONTINUE
end do
!
ID=197*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE198(N4,IDL)
401 CONTINUE
end do
!
ID=198*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE199(N4,IDL)
402 CONTINUE
end do
!
ID=199*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE200(N4,IDL)
403 CONTINUE
end do
!
ID=200*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE201(N4,IDL)
404 CONTINUE
end do
!
ID=201*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE202(N4,IDL)
405 CONTINUE
end do
!
ID=202*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE203(N4,IDL)
406 CONTINUE
end do
!
ID=203*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE204(N4,IDL)
407 CONTINUE
end do
!
ID=204*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE205(N4,IDL)
408 CONTINUE
end do
!
ID=205*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE206(N4,IDL)
409 CONTINUE
end do
!
ID=206*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE207(N4,IDL)
410 CONTINUE
end do
!
ID=207*IBLOK
do IDL=1,IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE208(N4,IDL)
411 CONTINUE
end do
!
ID=208*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE209(N4,IDL)
412 CONTINUE
end do
!
ID=209*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE210(N4,IDL)
413 CONTINUE
end do
!
ID=210*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE211(N4,IDL)
414 CONTINUE
end do
!
ID=211*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE212(N4,IDL)
415 CONTINUE
end do
!
ID=212*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE213(N4,IDL)
416 CONTINUE
end do
!
ID=213*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE214(N4,IDL)
417 CONTINUE
end do
!
ID=214*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE215(N4,IDL)
418 CONTINUE
end do
!
ID=215*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE216(N4,IDL)
419 CONTINUE
end do
!
ID=216*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE217(N4,IDL)
420 CONTINUE
end do
!
ID=217*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE218(N4,IDL)
421 CONTINUE
end do
!     
ID=218*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE219(N4,IDL)
422 CONTINUE
end do
!
ID=219*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE220(N4,IDL)
423 CONTINUE
end do
!
ID=220*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE221(N4,IDL)
424 CONTINUE
end do
!
ID=221*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE222(N4,IDL)
425 CONTINUE
end do
!
ID=222*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE223(N4,IDL)
426 CONTINUE
end do
!
ID=223*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE224(N4,IDL)
427 CONTINUE
end do
!     
ID=224*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE225(N4,IDL)
428 CONTINUE
end do
!
ID=225*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE226(N4,IDL)
429 CONTINUE
end do
!     
ID=226*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE227(N4,IDL)
430 CONTINUE
end do
!
ID=227*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE228(N4,IDL)
431 CONTINUE
end do
!
ID=228*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE229(N4,IDL)
432 CONTINUE
end do
!
ID=229*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE230(N4,IDL)
433 CONTINUE
end do
!
ID=230*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE231(N4,IDL)
434 CONTINUE
end do
!
ID=231*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE232(N4,IDL)
435 CONTINUE
end do
!
ID=232*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE233(N4,IDL)
436 CONTINUE
end do
!
ID=233*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE234(N4,IDL)
437 CONTINUE
end do
!
ID=234*IBLOK
do IDL=1, IBLOK
ID=ID+1
ALINEP(N4,ID)=ALINE235(N4,IDL)
438 CONTINUE
end do
!          
2 CONTINUE
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCC FULL J LINES
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
DO N4=1, LX4
!         
DO IDL=1, NOPJ0PP
ALINEJ0(N4,IDL)=ALINEJ0PP(N4,IDL)
end do
!
ID=NOPJ0PP
DO IDL=1, NOPJ0PM
ID=ID+1
ALINEJ0(N4,ID)=ALINEJ0PM(N4,IDL)
end do
!
ID=NOPJ0PP+NOPJ0PM
DO IDL=1, NOPJ0MP
ID=ID+1
ALINEJ0(N4,ID)=ALINEJ0MP(N4,IDL)
end do
!
ID=NOPJ0PP+NOPJ0PM+NOPJ0MP
DO IDL=1, NOPJ0MM
ID=ID+1
ALINEJ0(N4,ID)=ALINEJ0MM(N4,IDL)
end do
!
!     *************************************************************
!
DO IDL=1, NOPJ1P
ALINEJ1(N4,IDL)=ALINEJ1P(N4,IDL)
end do
!
ID=NOPJ1P
DO IDL=1, NOPJ1M
ID=ID+1
ALINEJ1(N4,ID)=ALINEJ1M(N4,IDL)
end do
!
!     *************************************************************
!
!         
DO IDL=1, NOPJ2PP
ALINEJ2(N4,IDL)=ALINEJ2PP(N4,IDL)
end do
!
ID=NOPJ2PP
DO IDL=1, NOPJ2PM
ID=ID+1
ALINEJ2(N4,ID)=ALINEJ2PM(N4,IDL)
end do
!
ID=NOPJ2PP+NOPJ2PM
DO IDL=1, NOPJ2MP
ID=ID+1
ALINEJ2(N4,ID)=ALINEJ2MP(N4,IDL)
end do
!
ID=NOPJ2PP+NOPJ2PM+NOPJ2MP
DO IDL=1, NOPJ2MM
ID=ID+1
ALINEJ2(N4,ID)=ALINEJ2MM(N4,IDL)
end do
!     
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCC FULL J LINES
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
DO N4=1, LX4
!
DO IDL=1, NOPFULJ0
ALINEJALL(N4,IDL)=ALINEJ0(N4,IDL)
end do
!
ID=NOPFULJ0
DO IDL=1, NOPFULJ1
ID=ID+1
ALINEJALL(N4,ID)=ALINEJ1(N4,IDL)
end do
!
ID=NOPFULJ0+NOPFULJ1
DO IDL=1, NOPFULJ2
ID=ID+1
ALINEJALL(N4,ID)=ALINEJ2(N4,IDL)
end do
!         
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
do N4=1,LX4
!     
do IDL=1,IBLOK
DO IJ=1,2
ALINEPMOM(IJ,N4,IDL)=ALINEMOM2(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,IDL)=ALINEMOM2(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,IDL)=ALINEMOM27(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,IDL)=ALINEMOM3(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,IDL)=ALINEMOM4(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,IDL)=ALINEMOM31(N4,IDL,IJ)
!     
end do
103 CONTINUE
end do
!
ID=IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM3(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM5(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM33(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM6(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM7(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM37(N4,IDL,IJ)
! 
end do
104 CONTINUE
end do
!
ID=2*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM4(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM8(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM39(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM9(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM10(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM43(N4,IDL,IJ)
! 
end do
105 CONTINUE
end do
!     
ID=3*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM5(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM11(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM45(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM12(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM13(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM49(N4,IDL,IJ)
! 
end do
106 CONTINUE
end do
!
ID=4*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM6(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM14(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM57(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM15(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM16(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM61(N4,IDL,IJ)
! 
end do
107 CONTINUE
end do
!
ID=5*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM7(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM17(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM63(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM18(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM19(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM67(N4,IDL,IJ)
! 
end do
108 CONTINUE
end do
!
ID=6*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM8(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM20(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM69(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM21(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM22(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM73(N4,IDL,IJ)
! 
end do
109 CONTINUE
end do
!
ID=7*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM9(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM23(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM75(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM24(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM25(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM79(N4,IDL,IJ)
! 
end do
110 CONTINUE
end do
!
ID=8*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM10(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM26(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM81(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM28(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM30(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM85(N4,IDL,IJ)
!               
end do
111 CONTINUE
end do
!
ID=9*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM11(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM32(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM88(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM29(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM36(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM96(N4,IDL,IJ)
! 
end do
112 CONTINUE
end do
!
ID=10*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM12(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM38(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM89(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM34(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM42(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM97(N4,IDL,IJ)
!
end do
113 CONTINUE
end do
!
ID=11*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM13(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM44(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM105(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM35(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM48(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM103(N4,IDL,IJ)
!
end do
114 CONTINUE
end do
!
ID=12*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM14(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM50(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM111(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM40(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM52(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM109(N4,IDL,IJ)
!
end do
115 CONTINUE
end do
!
ID=13*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM15(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM53(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM117(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM41(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM55(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM115(N4,IDL,IJ)
!
end do
116 CONTINUE
end do
!            
ID=14*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM16(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM56(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM123(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM46(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM60(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM121(N4,IDL,IJ)
!
end do
117 CONTINUE
end do
!
ID=15*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM17(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM62(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM129(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM47(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM66(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM127(N4,IDL,IJ)
!
end do
118 CONTINUE
end do
!
ID=16*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM18(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM68(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM135(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM51(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM72(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM133(N4,IDL,IJ)
!
end do
119 CONTINUE
end do
!
ID=17*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM19(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM74(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM141(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM54(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM78(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM139(N4,IDL,IJ)
!
end do
120 CONTINUE
end do
!
ID=18*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM20(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM80(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM148(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM58(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM84(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM145(N4,IDL,IJ)
!               
end do
121 CONTINUE
end do
!
ID=19*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM21(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM86(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM149(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM59(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM94(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM154(N4,IDL,IJ)
!
end do
122 CONTINUE
end do
!
ID=20*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM22(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM87(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM158(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM64(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM95(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM155(N4,IDL,IJ)
!
end do
123 CONTINUE
end do
!
ID=21*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM23(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM98(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM159(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM65(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM102(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM164(N4,IDL,IJ)
!               
end do
724 CONTINUE
end do
!
ID=22*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM24(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM227(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM168(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM70(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM108(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM165(N4,IDL,IJ)
!
end do
124 CONTINUE
end do
!
ID=23*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM25(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM104(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM169(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM71(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM114(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM174(N4,IDL,IJ)
!
end do
125 CONTINUE
end do
!
ID=24*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM26(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM110(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM178(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM76(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM120(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM175(N4,IDL,IJ)
!
end do
126 CONTINUE
end do
!
ID=25*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM27(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM116(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM179(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM77(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM126(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM184(N4,IDL,IJ)
!
end do
127 CONTINUE
end do
!
ID=26*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM28(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM122(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM188(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM82(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM132(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM185(N4,IDL,IJ)
!
end do
128 CONTINUE
end do
!
ID=27*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM29(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM128(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM189(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM83(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM138(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM194(N4,IDL,IJ)
!
end do
129 CONTINUE
end do
!
ID=28*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM30(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM134(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM198(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM90(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM144(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM195(N4,IDL,IJ)
!
end do
130 CONTINUE
end do
!
ID=29*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM31(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM140(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM199(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM91(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM152(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM204(N4,IDL,IJ)
!              
end do
131 CONTINUE
end do
!
ID=30*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM32(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM146(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM208(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM92(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM153(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM205(N4,IDL,IJ)
!
end do
132 CONTINUE
end do
!
ID=31*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM33(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM147(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM209(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM93(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM162(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM214(N4,IDL,IJ)
!
end do
133 CONTINUE
end do
!
ID=32*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM34(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM156(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM218(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM100(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM163(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM215(N4,IDL,IJ)
!
end do
134 CONTINUE
end do
!
ID=33*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM35(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM157(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM219(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM101(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM172(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM224(N4,IDL,IJ)
!
end do
135 CONTINUE
end do
!
ID=34*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM36(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM166(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM228(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM106(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM173(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM225(N4,IDL,IJ)
!
end do
136 CONTINUE
end do
!
ID=35*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM37(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM167(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM229(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM107(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM182(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM234(N4,IDL,IJ)
!
end do
137 CONTINUE
end do
!
ID=36*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM38(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM176(N4,IDL,IJ)
ALINEMOMJ0M(IJ,N4,ID)=ALINEMOM99(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM112(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM183(N4,IDL,IJ)
ALINEMOMJ2M(IJ,N4,ID)=ALINEMOM235(N4,IDL,IJ)
!
end do
138 CONTINUE
end do
!
ID=37*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM39(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM177(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM113(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM192(N4,IDL,IJ)
!
end do
139 CONTINUE
end do
!
ID=38*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM40(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM186(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM118(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM193(N4,IDL,IJ)
!
end do
140 CONTINUE
end do
!
ID=39*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM41(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM187(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM119(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM202(N4,IDL,IJ)
!
end do
141 CONTINUE
end do
!
ID=40*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM42(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM196(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM124(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM203(N4,IDL,IJ)
!
end do
142 CONTINUE
end do
!
ID=41*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM43(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM197(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM125(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM212(N4,IDL,IJ)
!
end do
143 CONTINUE
end do
!
ID=42*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM44(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM206(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM130(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM213(N4,IDL,IJ)
!
end do
144 CONTINUE
end do
!
ID=43*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM45(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM207(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM131(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM222(N4,IDL,IJ)
!
end do
145 CONTINUE
end do
!
ID=44*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM46(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM216(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM136(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM223(N4,IDL,IJ)
!
end do
146 CONTINUE
end do
!
ID=45*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM47(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM217(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM137(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM232(N4,IDL,IJ)
!
end do
147 CONTINUE
end do
!
ID=46*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM48(N4,IDL,IJ)
!
ALINEMOMJ0P(IJ,N4,ID)=ALINEMOM226(N4,IDL,IJ)
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM142(N4,IDL,IJ)
ALINEMOMJ2P(IJ,N4,ID)=ALINEMOM233(N4,IDL,IJ)
!
end do
148 CONTINUE
end do
!
ID=47*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM49(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM143(N4,IDL,IJ)
!
end do
149 CONTINUE
end do
!
ID=48*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM50(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM150(N4,IDL,IJ)
!
end do
150 CONTINUE
end do
!
ID=49*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM51(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM151(N4,IDL,IJ)
!
end do
151 CONTINUE
end do
!
ID=50*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM52(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM160(N4,IDL,IJ)
!
end do
152 CONTINUE
end do
!
ID=51*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM53(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM161(N4,IDL,IJ)
!
end do
153 CONTINUE
end do
!
ID=52*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM54(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM170(N4,IDL,IJ)
!
end do
154 CONTINUE
end do
!
ID=53*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM55(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM171(N4,IDL,IJ)
!
end do
155 CONTINUE
end do
!
ID=54*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM56(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM180(N4,IDL,IJ)
!
end do
156 CONTINUE
end do
!
ID=55*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM57(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM181(N4,IDL,IJ)
!
end do
157 CONTINUE
end do
!
ID=56*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM58(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM190(N4,IDL,IJ)
!
end do
158 CONTINUE
end do
!
ID=57*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM59(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM191(N4,IDL,IJ)
!
end do
159 CONTINUE
end do
!     
ID=58*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM60(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM200(N4,IDL,IJ)
!
end do
160 CONTINUE
end do
!
ID=59*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM61(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM201(N4,IDL,IJ)
!
end do
161 CONTINUE
end do
!
ID=60*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM62(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM210(N4,IDL,IJ)
!
end do
162 CONTINUE
end do
!
ID=61*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM63(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM211(N4,IDL,IJ)
!
end do
163 CONTINUE
end do
!
ID=62*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM64(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM220(N4,IDL,IJ)
!
end do
164 CONTINUE
end do
!
ID=63*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM65(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM221(N4,IDL,IJ)
!
end do
165 CONTINUE
end do
!
ID=64*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM66(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM230(N4,IDL,IJ)
!
end do
166 CONTINUE
end do
!
ID=65*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM67(N4,IDL,IJ)
!
ALINEMOMJ1(IJ,N4,ID)=ALINEMOM231(N4,IDL,IJ)
!
end do
167 CONTINUE
end do
!
ID=66*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM68(N4,IDL,IJ)
end do
168 CONTINUE
end do
!
ID=67*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM69(N4,IDL,IJ)
end do
169 CONTINUE
end do
!
ID=68*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM70(N4,IDL,IJ)
end do
170 CONTINUE
end do
!
ID=69*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM71(N4,IDL,IJ)
end do
171 CONTINUE
end do
!
ID=70*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM72(N4,IDL,IJ)
end do
172 CONTINUE
end do
!
ID=71*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM73(N4,IDL,IJ)
end do
173 CONTINUE
end do
!
ID=72*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM74(N4,IDL,IJ)
end do
174 CONTINUE
end do
!
ID=73*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM75(N4,IDL,IJ)
end do
175 CONTINUE
end do
!
ID=74*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM76(N4,IDL,IJ)
end do
176 CONTINUE
end do
!
ID=75*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM77(N4,IDL,IJ)
end do
177 CONTINUE
end do
!
ID=76*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM78(N4,IDL,IJ)
end do
178 CONTINUE
end do
!
ID=77*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM79(N4,IDL,IJ)
end do
179 CONTINUE
end do
!
ID=78*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM80(N4,IDL,IJ)
end do
180 CONTINUE
end do
!
ID=79*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM81(N4,IDL,IJ)
end do
181 CONTINUE
end do
!
ID=80*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM82(N4,IDL,IJ)
end do
182 CONTINUE
end do
!
ID=81*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM83(N4,IDL,IJ)
end do
183 CONTINUE
end do
!
ID=81*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM84(N4,IDL,IJ)
end do
184 CONTINUE
end do
!
ID=83*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM85(N4,IDL,IJ)
end do
185 CONTINUE
end do
!
ID=84*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM86(N4,IDL,IJ)
end do
186 CONTINUE
end do
!
ID=85*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM87(N4,IDL,IJ)
end do
187 CONTINUE
end do
!
ID=86*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM88(N4,IDL,IJ)
end do
188 CONTINUE
end do
!
ID=87*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM89(N4,IDL,IJ)
end do
189 CONTINUE
end do
!
ID=88*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM90(N4,IDL,IJ)
end do
190 CONTINUE
end do
!
ID=89*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM91(N4,IDL,IJ)
end do
191 CONTINUE
end do
!
ID=90*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM92(N4,IDL,IJ)
end do
192 CONTINUE
end do
!
ID=91*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM93(N4,IDL,IJ)
end do
193 CONTINUE
end do
!
ID=92*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM94(N4,IDL,IJ)
end do
194 CONTINUE
end do
!
ID=93*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM95(N4,IDL,IJ)
end do
195 CONTINUE
end do
!
ID=94*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM96(N4,IDL,IJ)
end do
196 CONTINUE
end do
!
ID=95*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM97(N4,IDL,IJ)
end do
197 CONTINUE
end do
!
ID=96*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM98(N4,IDL,IJ)
end do
198 CONTINUE
end do
!
ID=97*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM99(N4,IDL,IJ)
end do
199 CONTINUE
end do
!
ID=98*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM100(N4,IDL,IJ)
end do
200 CONTINUE
end do
!
ID=99*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM101(N4,IDL,IJ)
end do
201 CONTINUE
end do
!
ID=100*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM102(N4,IDL,IJ)
end do
202 CONTINUE
end do
!
ID=101*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM103(N4,IDL,IJ)
end do
203 CONTINUE
end do
!
ID=102*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM104(N4,IDL,IJ)
end do
204 CONTINUE
end do
!
ID=103*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM105(N4,IDL,IJ)
end do
205 CONTINUE
end do
!
ID=104*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM106(N4,IDL,IJ)
end do
206 CONTINUE
end do
!
ID=105*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM107(N4,IDL,IJ)
end do
207 CONTINUE
end do
!
ID=106*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM108(N4,IDL,IJ)
end do
208 CONTINUE
end do
!
ID=107*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM109(N4,IDL,IJ)
end do
209 CONTINUE
end do
!
ID=108*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM110(N4,IDL,IJ)
end do
210 CONTINUE
end do
!
ID=109*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM111(N4,IDL,IJ)
end do
211 CONTINUE
end do
!
ID=110*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM112(N4,IDL,IJ)
end do
212 CONTINUE
end do
!
ID=111*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM113(N4,IDL,IJ)
end do
213 CONTINUE
end do
!
ID=112*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM114(N4,IDL,IJ)
end do
214 CONTINUE
end do
!
ID=113*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM115(N4,IDL,IJ)
end do
215 CONTINUE
end do
!
ID=114*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM116(N4,IDL,IJ)
end do
216 CONTINUE
end do
!
ID=115*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM117(N4,IDL,IJ)
end do
217 CONTINUE
end do
!
ID=116*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM118(N4,IDL,IJ)
end do
218 CONTINUE
end do
!
ID=117*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM119(N4,IDL,IJ)
end do
219 CONTINUE
end do
!
ID=118*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM120(N4,IDL,IJ)
end do
220 CONTINUE
end do
!
ID=119*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM121(N4,IDL,IJ)
end do
221 CONTINUE
end do
!
ID=120*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM122(N4,IDL,IJ)
end do
222 CONTINUE
end do
!
ID=121*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM123(N4,IDL,IJ)
end do
223 CONTINUE
end do
!
ID=122*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM124(N4,IDL,IJ)
end do
224 CONTINUE
end do
!
ID=123*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM125(N4,IDL,IJ)
end do
225 CONTINUE
end do
!
ID=124*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM126(N4,IDL,IJ)
end do
226 CONTINUE
end do
!
ID=125*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM127(N4,IDL,IJ)
end do
227 CONTINUE
end do
!
ID=126*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM128(N4,IDL,IJ)
end do
228 CONTINUE
end do
!
ID=127*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM129(N4,IDL,IJ)
end do
229 CONTINUE
end do
!
ID=128*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM130(N4,IDL,IJ)
end do
230 CONTINUE
end do
!
ID=129*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM131(N4,IDL,IJ)
end do
231 CONTINUE
end do
!
ID=130*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM132(N4,IDL,IJ)
end do
232 CONTINUE
end do
!
ID=131*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM133(N4,IDL,IJ)
end do
233 CONTINUE
end do
!
ID=132*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM134(N4,IDL,IJ)
end do
234 CONTINUE
end do
!
ID=133*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM135(N4,IDL,IJ)
end do
235 CONTINUE
end do
!
ID=134*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM136(N4,IDL,IJ)
end do
236 CONTINUE
end do
!
ID=135*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM137(N4,IDL,IJ)
end do
237 CONTINUE
end do
!
ID=136*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM138(N4,IDL,IJ)
end do
238 CONTINUE
end do
!
ID=137*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM139(N4,IDL,IJ)
end do
239 CONTINUE
end do
!
ID=138*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM140(N4,IDL,IJ)
end do
240 CONTINUE
end do
!
ID=139*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM141(N4,IDL,IJ)
end do
241 CONTINUE
end do
!
ID=140*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM142(N4,IDL,IJ)
end do
242 CONTINUE
end do
!
ID=141*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM143(N4,IDL,IJ)
end do
243 CONTINUE
end do
!
ID=142*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM144(N4,IDL,IJ)
end do
244 CONTINUE
end do
!
ID=143*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM145(N4,IDL,IJ)
end do
245 CONTINUE
end do
!CCCCCCCCCCCCCCCCCC
ID=144*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM146(N4,IDL,IJ)
end do
246 CONTINUE
end do
!
ID=145*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM147(N4,IDL,IJ)
end do
247 CONTINUE
end do
!
ID=146*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM148(N4,IDL,IJ)
end do
248 CONTINUE
end do
!
ID=147*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM149(N4,IDL,IJ)
end do
249 CONTINUE
end do
!
ID=148*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM150(N4,IDL,IJ)
end do
250 CONTINUE
end do
!
ID=149*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM151(N4,IDL,IJ)
end do
251 CONTINUE
end do
!
ID=150*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM152(N4,IDL,IJ)
end do
252 CONTINUE
end do
!
ID=151*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM153(N4,IDL,IJ)
end do
253 CONTINUE
end do
!
ID=152*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM154(N4,IDL,IJ)
end do
254 CONTINUE
end do
!
ID=153*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM155(N4,IDL,IJ)
end do
255 CONTINUE
end do
!
ID=154*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM156(N4,IDL,IJ)
end do
256 CONTINUE
end do
!
ID=155*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM157(N4,IDL,IJ)
end do
257 CONTINUE
end do
!
ID=156*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM158(N4,IDL,IJ)
end do
258 CONTINUE
end do
!
ID=157*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM159(N4,IDL,IJ)
end do
259 CONTINUE
end do
!     
ID=158*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM160(N4,IDL,IJ)
end do
260 CONTINUE
end do
!
ID=159*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM161(N4,IDL,IJ)
end do
261 CONTINUE
end do
!
ID=160*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM162(N4,IDL,IJ)
end do
262 CONTINUE
end do
!
ID=161*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM163(N4,IDL,IJ)
end do
263 CONTINUE
end do
!
ID=162*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM164(N4,IDL,IJ)
end do
264 CONTINUE
end do
!
ID=163*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM165(N4,IDL,IJ)
end do
265 CONTINUE
end do
!
ID=164*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM166(N4,IDL,IJ)
end do
266 CONTINUE
end do
!
ID=165*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM167(N4,IDL,IJ)
end do
267 CONTINUE
end do
!
ID=166*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM168(N4,IDL,IJ)
end do
268 CONTINUE
end do
!
ID=167*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM169(N4,IDL,IJ)
end do
269 CONTINUE
end do
!
ID=168*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM170(N4,IDL,IJ)
end do
270 CONTINUE
end do
!
ID=169*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM171(N4,IDL,IJ)
end do
271 CONTINUE
end do
!
ID=170*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM172(N4,IDL,IJ)
end do
272 CONTINUE
end do
!
ID=171*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM173(N4,IDL,IJ)
end do
273 CONTINUE
end do
!
ID=172*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM174(N4,IDL,IJ)
end do
274 CONTINUE
end do
!
ID=173*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM175(N4,IDL,IJ)
end do
275 CONTINUE
end do
!
ID=174*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM176(N4,IDL,IJ)
end do
276 CONTINUE
end do
!
ID=175*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM177(N4,IDL,IJ)
end do
277 CONTINUE
end do
!
ID=176*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM178(N4,IDL,IJ)
end do
278 CONTINUE
end do
!
ID=177*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM179(N4,IDL,IJ)
end do
279 CONTINUE
end do
!
ID=178*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM180(N4,IDL,IJ)
end do
280 CONTINUE
end do
!
ID=179*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM181(N4,IDL,IJ)
end do
281 CONTINUE
end do
!
ID=180*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM182(N4,IDL,IJ)
end do
282 CONTINUE
end do
!
ID=181*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM183(N4,IDL,IJ)
end do
283 CONTINUE
end do
!
ID=181*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM184(N4,IDL,IJ)
end do
284 CONTINUE
end do
!
ID=183*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM185(N4,IDL,IJ)
end do
285 CONTINUE
end do
!
ID=184*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM186(N4,IDL,IJ)
end do
286 CONTINUE
end do
!
ID=185*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM187(N4,IDL,IJ)
end do
287 CONTINUE
end do
!
ID=186*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM188(N4,IDL,IJ)
end do
288 CONTINUE
end do
!
ID=187*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM189(N4,IDL,IJ)
end do
289 CONTINUE
end do
!
ID=188*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM190(N4,IDL,IJ)
end do
290 CONTINUE
end do
!
ID=189*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM191(N4,IDL,IJ)
end do
291 CONTINUE
end do
!
ID=190*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM192(N4,IDL,IJ)
end do
292 CONTINUE
end do
!
ID=191*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM193(N4,IDL,IJ)
end do
293 CONTINUE
end do
!
ID=192*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM194(N4,IDL,IJ)
end do
294 CONTINUE
end do
!
ID=193*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM195(N4,IDL,IJ)
end do
295 CONTINUE
end do
!
ID=194*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM196(N4,IDL,IJ)
end do
296 CONTINUE
end do
!
ID=195*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM197(N4,IDL,IJ)
end do
297 CONTINUE
end do
!
ID=196*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM198(N4,IDL,IJ)
end do
298 CONTINUE
end do
!
ID=197*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM199(N4,IDL,IJ)
end do
299 CONTINUE
end do
!
ID=198*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM200(N4,IDL,IJ)
end do
500 CONTINUE
end do
!
ID=199*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM201(N4,IDL,IJ)
end do
501 CONTINUE
end do
!
ID=200*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM202(N4,IDL,IJ)
end do
502 CONTINUE
end do
!
ID=201*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM203(N4,IDL,IJ)
end do
503 CONTINUE
end do
!
ID=202*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM204(N4,IDL,IJ)
end do
504 CONTINUE
end do
!
ID=203*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM205(N4,IDL,IJ)
end do
505 CONTINUE
end do
!
ID=204*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM206(N4,IDL,IJ)
end do
506 CONTINUE
end do
!
ID=205*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM207(N4,IDL,IJ)
end do
507 CONTINUE
end do
!
ID=206*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM208(N4,IDL,IJ)
end do
508 CONTINUE
end do
!
ID=207*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM209(N4,IDL,IJ)
end do
509 CONTINUE
end do
!
ID=208*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM210(N4,IDL,IJ)
end do
510 CONTINUE
end do
!
ID=209*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM211(N4,IDL,IJ)
end do
511 CONTINUE
end do
!
ID=210*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM212(N4,IDL,IJ)
end do
512 CONTINUE
end do
!
ID=211*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM213(N4,IDL,IJ)
end do
513 CONTINUE
end do
!
ID=212*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM214(N4,IDL,IJ)
end do
514 CONTINUE
end do
!
ID=213*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM215(N4,IDL,IJ)
end do
515 CONTINUE
end do
!
ID=214*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM216(N4,IDL,IJ)
end do
516 CONTINUE
end do
!
ID=215*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM217(N4,IDL,IJ)
end do
517 CONTINUE
end do
!
ID=216*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM218(N4,IDL,IJ)
end do
518 CONTINUE
end do
!
ID=217*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM219(N4,IDL,IJ)
end do
519 CONTINUE
end do
!
ID=218*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM220(N4,IDL,IJ)
end do
520 CONTINUE
end do
!
ID=219*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM221(N4,IDL,IJ)
end do
521 CONTINUE
end do
!
ID=220*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM222(N4,IDL,IJ)
end do
522 CONTINUE
end do
!
ID=221*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM223(N4,IDL,IJ)
end do
523 CONTINUE
end do
!
ID=222*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM224(N4,IDL,IJ)
end do
524 CONTINUE
end do
!
ID=223*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM225(N4,IDL,IJ)
end do
525 CONTINUE
end do
!
ID=224*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM226(N4,IDL,IJ)
end do
526 CONTINUE
end do
!
ID=225*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM227(N4,IDL,IJ)
end do
527 CONTINUE
end do
!
ID=226*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM228(N4,IDL,IJ)
end do
528 CONTINUE
end do
!
ID=227*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM229(N4,IDL,IJ)
end do
529 CONTINUE
end do
!
ID=228*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM230(N4,IDL,IJ)
end do
530 CONTINUE
end do
!
ID=229*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM231(N4,IDL,IJ)
end do
531 CONTINUE
end do
!
ID=230*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM232(N4,IDL,IJ)
end do
532 CONTINUE
end do
!
ID=231*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM233(N4,IDL,IJ)
end do
533 CONTINUE
end do
!
ID=232*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM234(N4,IDL,IJ)
end do
534 CONTINUE
end do
!
ID=233*IBLOK
do IDL=1,IBLOK
ID=ID+1
DO IJ=1,2
ALINEPMOM(IJ,N4,ID)=ALINEMOM235(N4,IDL,IJ)
end do
535 CONTINUE
end do
!
!
102 CONTINUE
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCC FULL J MOM-LINES
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
DO N4=1, LX4
!         
DO IDL=1, NOPJ0PMOM
DO IJ=1, 2
ALINEMOMJ0(IJ,N4,IDL)=ALINEMOMJ0P(IJ,N4,IDL)
end do
end do
!
ID=NOPJ0PMOM
DO IDL=1, NOPJ0MMOM
ID=ID+1
DO IJ=1, 2
ALINEMOMJ0(IJ,N4,ID)=ALINEMOMJ0M(IJ,N4,IDL)
end do
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
DO IDL=1, NOPJ2PMOM
DO IJ=1, 2
ALINEMOMJ2(IJ,N4,IDL)=ALINEMOMJ2P(IJ,N4,IDL)
end do
end do
!
ID=NOPJ2PMOM
DO IDL=1, NOPJ2MMOM
ID=ID+1
DO IJ=1, 2
ALINEMOMJ2(IJ,N4,ID)=ALINEMOMJ2M(IJ,N4,IDL)
end do
end do
!         
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCC FULL J ALL MOM-LINES
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
DO N4=1, LX4
!
DO IDL=1, NOPFULJ0MOM
DO IJ=1, 2
ALINEMOMJALL(IJ,N4,IDL)=ALINEMOMJ0(IJ,N4,IDL)
end do
end do
!
ID=NOPFULJ0MOM
DO IDL=1, NOPFULJ1MOM
ID=ID+1
DO IJ=1, 2
ALINEMOMJALL(IJ,N4,ID)=ALINEMOMJ1(IJ,N4,IDL)
end do
end do
!
ID=NOPFULJ0MOM+NOPFULJ1MOM
DO IDL=1, NOPFULJ2MOM
ID=ID+1
DO IJ=1, 2
ALINEMOMJALL(IJ,N4,ID)=ALINEMOMJ2(IJ,N4,IDL)
end do
end do
!         
end do
!***************************************************************
do N4=1,LX4
!***************************************************************     
do ID=1, NTOTAL
ALLP=ALINEP(N4,ID)
AVACLP(JBIN,ID)=AVACLP(JBIN,ID)+ALLP
606 CONTINUE
end do
!
do ID=1, NOPJ0PP
ALLPJ0PP=ALINEJ0PP(N4,ID)
AVACLJ0PP(JBIN,ID)=AVACLJ0PP(JBIN,ID)+ALLPJ0PP
607 CONTINUE
end do
!
do ID=1, NOPJ0PM
ALLPJ0PM=ALINEJ0PM(N4,ID)
AVACLJ0PM(JBIN,ID)=AVACLJ0PM(JBIN,ID)+ALLPJ0PM
608 CONTINUE
end do
!
do ID=1, NOPJ0MP
ALLPJ0MP=ALINEJ0MP(N4,ID)
AVACLJ0MP(JBIN,ID)=AVACLJ0MP(JBIN,ID)+ALLPJ0MP
609 CONTINUE
end do
!
do ID=1, NOPJ0MM
ALLPJ0MM=ALINEJ0MM(N4,ID)
AVACLJ0MM(JBIN,ID)=AVACLJ0MM(JBIN,ID)+ALLPJ0MM
610 CONTINUE
end do
!
do ID=1, NOPJ1P
ALLPJ1P=ALINEJ1P(N4,ID)
AVACLJ1P(JBIN,ID)=AVACLJ1P(JBIN,ID)+ALLPJ1P
611 CONTINUE
end do
!
do ID=1, NOPJ1M
ALLPJ1M=ALINEJ1M(N4,ID)
AVACLJ1M(JBIN,ID)=AVACLJ1M(JBIN,ID)+ALLPJ1M
612 CONTINUE
end do
!
do ID=1, NOPJ2PP
ALLPJ2PP=ALINEJ2PP(N4,ID)
AVACLJ2PP(JBIN,ID)=AVACLJ2PP(JBIN,ID)+ALLPJ2PP
613 CONTINUE
end do
!
do ID=1, NOPJ2PM
ALLPJ2PM=ALINEJ2PM(N4,ID)
AVACLJ2PM(JBIN,ID)=AVACLJ2PM(JBIN,ID)+ALLPJ2PM
614 CONTINUE
end do
!
do ID=1, NOPJ2MP
ALLPJ2MP=ALINEJ2MP(N4,ID)
AVACLJ2MP(JBIN,ID)=AVACLJ2MP(JBIN,ID)+ALLPJ2MP
615 CONTINUE
end do
!
do ID=1, NOPJ2MM
ALLPJ2MM=ALINEJ2MM(N4,ID)
AVACLJ2MM(JBIN,ID)=AVACLJ2MM(JBIN,ID)+ALLPJ2MM
616 CONTINUE
end do
!
do ID=1, NOPFULJ0
ALLPJ0=ALINEJ0(N4,ID)
AVACLJ0(JBIN,ID)=AVACLJ0(JBIN,ID)+ALLPJ0
617 CONTINUE
end do
!
do ID=1, NOPFULJ1
ALLPJ1=ALINEJ1(N4,ID)
AVACLJ1(JBIN,ID)=AVACLJ1(JBIN,ID)+ALLPJ1
618 CONTINUE
end do
!
do ID=1, NOPFULJ2
ALLPJ2=ALINEJ2(N4,ID)
AVACLJ2(JBIN,ID)=AVACLJ2(JBIN,ID)+ALLPJ2
619 CONTINUE
end do
!
do ID=1, NTOTAL
ALLPJALL=ALINEJALL(N4,ID)
AVACLALL(JBIN,ID)=AVACLALL(JBIN,ID)+ALLPJALL
620 CONTINUE
end do
!         
do ID=1, NTOTALMOM
DO IJ=1,2
ALLPMOM=ALINEPMOM(IJ,N4,ID)
AVACLPMOM(IJ,JBIN,ID)=AVACLPMOM(IJ,JBIN,ID)+ALLPMOM
end do
807 CONTINUE
end do
!
do ID=1, NOPJ0PMOM
DO IJ=1,2
ALLPMOMJ0P=ALINEMOMJ0P(IJ,N4,ID)
AVACLMOMJ0P(IJ,JBIN,ID)=AVACLMOMJ0P(IJ,JBIN,ID) &
  & +ALLPMOMJ0P
end do
808 CONTINUE
end do
!
do ID=1, NOPJ0MMOM
DO IJ=1,2
ALLPMOMJ0M=ALINEMOMJ0M(IJ,N4,ID)
AVACLMOMJ0M(IJ,JBIN,ID)=AVACLMOMJ0M(IJ,JBIN,ID) &
  & +ALLPMOMJ0M
end do
809 CONTINUE
end do
!
do ID=1, NOPJ1MOM
DO IJ=1,2
ALLPMOMJ1=ALINEMOMJ1(IJ,N4,ID)
AVACLMOMJ1(IJ,JBIN,ID)=AVACLMOMJ1(IJ,JBIN,ID) &
  & +ALLPMOMJ1
end do
810 CONTINUE
end do
!
do ID=1, NOPJ2PMOM
DO IJ=1,2
ALLPMOMJ2P=ALINEMOMJ2P(IJ,N4,ID)
AVACLMOMJ2P(IJ,JBIN,ID)=AVACLMOMJ2P(IJ,JBIN,ID) &
  & +ALLPMOMJ2P
end do
811 CONTINUE
end do
!
do ID=1, NOPJ2MMOM
DO IJ=1,2
ALLPMOMJ2M=ALINEMOMJ2M(IJ,N4,ID)
AVACLMOMJ2M(IJ,JBIN,ID)=AVACLMOMJ2M(IJ,JBIN,ID) &
  & +ALLPMOMJ2M
end do
812 CONTINUE
end do
!
do ID=1, NOPFULJ0MOM
DO IJ=1,2
ALLPMOMJ0=ALINEMOMJ0(IJ,N4,ID)
AVACLMOMJ0(IJ,JBIN,ID)=AVACLMOMJ0(IJ,JBIN,ID) &
  & +ALLPMOMJ0
end do
813 CONTINUE
end do
!
do ID=1, NOPFULJ2MOM
DO IJ=1,2
ALLPMOMJ2=ALINEMOMJ2(IJ,N4,ID)
AVACLMOMJ2(IJ,JBIN,ID)=AVACLMOMJ2(IJ,JBIN,ID) &
  & +ALLPMOMJ2
end do
814 CONTINUE
end do
!
do ID=1, NTOTALMOM
DO IJ=1,2
ALLPMOMJALL=ALINEMOMJALL(IJ,N4,ID)
AVACLMOMALL(IJ,JBIN,ID)=AVACLMOMALL(IJ,JBIN,ID) &
  & +ALLPMOMJALL
end do
815 CONTINUE
end do
!         
605 CONTINUE
end do
!***************************************************************
!     FULL CORRELATOR
!***************************************************************      
INO=0
DO ID2=1, NTOTAL
DO ID1=1, NTOTAL
!
DO N4=1,LX4
DO NT=1,LMAXIR
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLP(JBIN,NT,ID1,ID2)=ACORLP(JBIN,NT,ID1,ID2) &
  & +(ALINEP(N4,ID1)*CONJG(ALINEP(N4X,ID2)) &
  & +CONJG(ALINEP(N4,ID2))*ALINEP(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=0, Pp=+, Pr=+, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ0PP
DO ID1=1,NOPJ0PP
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ0PP(JBIN,NT,ID1,ID2)=ACORLJ0PP(JBIN,NT,ID1,ID2) &
  & +(ALINEJ0PP(N4,ID1)*CONJG(ALINEJ0PP(N4X,ID2)) &
  & +CONJG(ALINEJ0PP(N4,ID2))*ALINEJ0PP(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=0, Pp=+, Pr=-, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ0PM
DO ID1=1,NOPJ0PM
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ0PM(JBIN,NT,ID1,ID2)=ACORLJ0PM(JBIN,NT,ID1,ID2) &
  & +(ALINEJ0PM(N4,ID1)*CONJG(ALINEJ0PM(N4X,ID2)) &
  & +CONJG(ALINEJ0PM(N4,ID2))*ALINEJ0PM(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=0, Pp=-, Pr=+, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ0MP
DO ID1=1,NOPJ0MP
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ0MP(JBIN,NT,ID1,ID2)=ACORLJ0MP(JBIN,NT,ID1,ID2) &
  & +(ALINEJ0MP(N4,ID1)*CONJG(ALINEJ0MP(N4X,ID2)) &
  & +CONJG(ALINEJ0MP(N4,ID2))*ALINEJ0MP(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=0, Pp=-, Pr=-, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ0MM
DO ID1=1,NOPJ0MM
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ0MM(JBIN,NT,ID1,ID2)=ACORLJ0MM(JBIN,NT,ID1,ID2) &
  & +(ALINEJ0MM(N4,ID1)*CONJG(ALINEJ0MM(N4X,ID2)) &
  & +CONJG(ALINEJ0MM(N4,ID2))*ALINEJ0MM(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=1, Pr=+, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ1P
DO ID1=1,NOPJ1P
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ1P(JBIN,NT,ID1,ID2)=ACORLJ1P(JBIN,NT,ID1,ID2) &
  & +(ALINEJ1P(N4,ID1)*CONJG(ALINEJ1P(N4X,ID2)) &
  & +CONJG(ALINEJ1P(N4,ID2))*ALINEJ1P(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=1, Pr=-, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ1M
DO ID1=1,NOPJ1M
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ1M(JBIN,NT,ID1,ID2)=ACORLJ1M(JBIN,NT,ID1,ID2) &
  & +(ALINEJ1M(N4,ID1)*CONJG(ALINEJ1M(N4X,ID2)) &
  & +CONJG(ALINEJ1M(N4,ID2))*ALINEJ1M(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=2, Pp=+, Pr=+, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ2PP
DO ID1=1,NOPJ2PP
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ2PP(JBIN,NT,ID1,ID2)=ACORLJ2PP(JBIN,NT,ID1,ID2) &
  & +(ALINEJ2PP(N4,ID1)*CONJG(ALINEJ2PP(N4X,ID2)) &
  & +CONJG(ALINEJ2PP(N4,ID2))*ALINEJ2PP(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=2, Pp=+, Pr=-, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ2PM
DO ID1=1,NOPJ2PM
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ2PM(JBIN,NT,ID1,ID2)=ACORLJ2PM(JBIN,NT,ID1,ID2) &
  & +(ALINEJ2PM(N4,ID1)*CONJG(ALINEJ2PM(N4X,ID2)) &
  & +CONJG(ALINEJ2PM(N4,ID2))*ALINEJ2PM(N4X,ID1))*0.5
!*******************************************
end do
end do
!     
end do
end do
!***************************************************************
!     J=2, Pp=-, Pr=+, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ2MP
DO ID1=1,NOPJ2MP
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ2MP(JBIN,NT,ID1,ID2)=ACORLJ2MP(JBIN,NT,ID1,ID2) &
  & +(ALINEJ2MP(N4,ID1)*CONJG(ALINEJ2MP(N4X,ID2)) &
  & +CONJG(ALINEJ2MP(N4,ID2))*ALINEJ2MP(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=2, Pp=-, Pr=-, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPJ2MM
DO ID1=1,NOPJ2MM
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ2MM(JBIN,NT,ID1,ID2)=ACORLJ2MM(JBIN,NT,ID1,ID2) &
  & +(ALINEJ2MM(N4,ID1)*CONJG(ALINEJ2MM(N4X,ID2)) &
  & +CONJG(ALINEJ2MM(N4,ID2))*ALINEJ2MM(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=0 FULL, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPFULJ0
DO ID1=1,NOPFULJ0
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ0(JBIN,NT,ID1,ID2)=ACORLJ0(JBIN,NT,ID1,ID2) &
  & +(ALINEJ0(N4,ID1)*CONJG(ALINEJ0(N4X,ID2)) &
  & +CONJG(ALINEJ0(N4,ID2))*ALINEJ0(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=1 FULL, q=0 correlator
!***************************************************************      
INO=0
DO ID2=1,NOPFULJ1
DO ID1=1,NOPFULJ1
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ1(JBIN,NT,ID1,ID2)=ACORLJ1(JBIN,NT,ID1,ID2) &
  & +(ALINEJ1(N4,ID1)*CONJG(ALINEJ1(N4X,ID2)) &
  & +CONJG(ALINEJ1(N4,ID2))*ALINEJ1(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=2 FULL, q=0 correlator
!***************************************************************
INO=0
DO ID2=1,NOPFULJ2
DO ID1=1,NOPFULJ2
!
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLJ2(JBIN,NT,ID1,ID2)=ACORLJ2(JBIN,NT,ID1,ID2) &
  & +(ALINEJ2(N4,ID1)*CONJG(ALINEJ2(N4X,ID2)) &
  & +CONJG(ALINEJ2(N4,ID2))*ALINEJ2(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!     J=ALL FULL, q=0 correlator
!***************************************************************      
INO=0
DO ID2=1,NTOTAL
DO ID1=1,NTOTAL
!
DO N4=1,LX4
DO NT=1,LMAXIR
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
ACORLALL(JBIN,NT,ID1,ID2)=ACORLALL(JBIN,NT,ID1,ID2) &
  & +(ALINEJALL(N4,ID1)*CONJG(ALINEJALL(N4X,ID2)) &
  & +CONJG(ALINEJALL(N4,ID2))*ALINEJALL(N4X,ID1))*0.5
!*******************************************
end do
end do
!
end do
end do
!***************************************************************
!***************************************************************      
!***************************************************************
!CCCCC   Momentum
!***************************************************************
!***************************************************************      
!***************************************************************      
!     ID12=0
INO=0
DO ID2=1, NTOTALMOM
DO ID1=1, NTOTALMOM
!     ID12=ID12+1
!     
DO N4=1,LX4
DO NT=1,LMAXIR
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOM(IJ,JBIN,NT,ID1,ID2)=ACORLPMOM(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEPMOM(IJ,N4,ID1)*CONJG(ALINEPMOM(IJ,N4X,ID2)) &
  & +CONJG(ALINEPMOM(IJ,N4,ID2))*ALINEPMOM(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!
end do
end do
!****************************************************************
!     J=0, Pp=+, q=1,2
!****************************************************************
INO=0
DO ID2=1,NOPJ0PMOM
DO ID1=1,NOPJ0PMOM
!     
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJ0P(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ0P(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJ0P(IJ,N4,ID1)*CONJG(ALINEMOMJ0P(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJ0P(IJ,N4,ID2))*ALINEMOMJ0P(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!               
end do
end do
!****************************************************************
!     J=0, Pp=-, q=1,2
!****************************************************************
INO=0
DO ID2=1,NOPJ0MMOM
DO ID1=1,NOPJ0MMOM
!     
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJ0M(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ0M(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJ0M(IJ,N4,ID1)*CONJG(ALINEMOMJ0M(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJ0M(IJ,N4,ID2))*ALINEMOMJ0M(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!
end do
end do
!****************************************************************
!     J=1, q=1,2
!****************************************************************
INO=0
DO ID2=1,NOPJ1MOM
DO ID1=1,NOPJ1MOM
!     
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJ1(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ1(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJ1(IJ,N4,ID1)*CONJG(ALINEMOMJ1(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJ1(IJ,N4,ID2))*ALINEMOMJ1(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!
end do
end do
!****************************************************************
!     J=2, Pp=+, q=1,2
!****************************************************************
INO=0
DO ID2=1,NOPJ2PMOM
DO ID1=1,NOPJ2PMOM
!     
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJ2P(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ2P(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJ2P(IJ,N4,ID1)*CONJG(ALINEMOMJ2P(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJ2P(IJ,N4,ID2))*ALINEMOMJ2P(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!                
end do
end do
!****************************************************************
!     J=2, Pp=-, q=1,2
!****************************************************************
INO=0
DO ID2=1,NOPJ2MMOM
DO ID1=1,NOPJ2MMOM
!     
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJ2M(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ2M(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJ2M(IJ,N4,ID1)*CONJG(ALINEMOMJ2M(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJ2M(IJ,N4,ID2))*ALINEMOMJ2M(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!
end do
end do
!****************************************************************
!     J=0, q=1,2
!****************************************************************
INO=0
DO ID2=1,NOPFULJ0MOM
DO ID1=1,NOPFULJ0MOM
!     
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJ0(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ0(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJ0(IJ,N4,ID1)*CONJG(ALINEMOMJ0(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJ0(IJ,N4,ID2))*ALINEMOMJ0(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!
end do
end do
!****************************************************************
!     J=2, q=1,2
!****************************************************************
INO=0
DO ID2=1,NOPFULJ2MOM
DO ID1=1,NOPFULJ2MOM
!     
DO N4=1,LX4
DO NT=1,LMAX
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJ2(IJ,JBIN,NT,ID1,ID2)=ACORLPMOMJ2(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJ2(IJ,N4,ID1)*CONJG(ALINEMOMJ2(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJ2(IJ,N4,ID2))*ALINEMOMJ2(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!  
end do
end do
!****************************************************************
!     J=ALL, q=1,2
!****************************************************************
INO=0
DO ID2=1, NTOTALMOM
DO ID1=1, NTOTALMOM
!     
DO N4=1,LX4
DO NT=1,LMAXIR
N4X=N4+NT-1
IF(N4X.GT.LX4) N4X=N4X-LX4
!*******************************************
DO IJ=1,2
!*******************************************
ACORLPMOMJALL(IJ,JBIN,NT,ID1,ID2)= &
  & ACORLPMOMJALL(IJ,JBIN,NT,ID1,ID2) &
  & +(ALINEMOMJALL(IJ,N4,ID1)*CONJG(ALINEMOMJALL(IJ,N4X,ID2)) &
  & +CONJG(ALINEMOMJALL(IJ,N4,ID2))*ALINEMOMJALL(IJ,N4X,ID1))*0.5
!*******************************************
end do
!*******************************************
end do
end do
!
end do
end do
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
RETURN
end subroutine POT_OLD
!*********************************************************************
!*********************************************************************
!  Here we set up blocked link pointers...
!*********************************************************************
!*********************************************************************
SUBROUTINE SETUPB
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52,IBLOK=5)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
!
COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
DIMENSION LB(IBLOK+1),LB1(IBLOK+1),LB12(IBLOK+1)
!
LX12=LX1*LX2
LX123=LX12*LX3
LB(1)=1
LB1(1)=LX1
LB12(1)=LX12
do IBL=2,IBLOK+1
LB(IBL)=2*LB(IBL-1)
LB1(IBL)=2*LB1(IBL-1)
LB12(IBL)=2*LB12(IBL-1)
1 CONTINUE
end do
!
NN=0
DO L3=1,LX3
DO L2=1,LX2
DO L1=1,LX1
NN=NN+1
do ID=1,IBLOK+1
!     
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
!     
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
!     
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
!     
3 CONTINUE
end do
end do
end do
end do
!     
RETURN
end subroutine SETUPB
!*********************************************************************
!*********************************************************************
!                        SUBROUTINE ACTION                           *
!*********************************************************************
!*********************************************************************
SUBROUTINE ACTION(IPR,ITER,TOTACT)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(NITER=2,NUMBIN=2)
PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
COMMON/ASTORE/ACTN(NUMBIN,6),PLAQ(NITER,6)
COMMON/COMMUNICATE/AVACS(NUMBIN),AVACT(NUMBIN)
COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
DIMENSION DUM11(NCOL2) &
  & ,A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2)
complex(real32) U11,A11,B11,C11,D11,DUM11,ACT
!
DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN)
DIMENSION VAL(NUMBIN),AV(NUMBIN)
!
IBIN=NITER/NUMBIN
!
IF(IPR.EQ.0)THEN
JBIN=(ITER-1)/IBIN+1
!
IF(ITER.EQ.1)THEN
do NB=1,NUMBIN
AVAC(NB)=0.0d0
AVACSQ(NB)=0.0d0
AVACS(NB)=0.0d0
AVACT(NB)=0.0d0
3 CONTINUE
end do
end if
!
VACTS=0.0
VACTT=0.0
!
DO NN=1,LSIZE
DO MU=1,3
M1=NN
!     
do IJ=1,NCOL2
A11(IJ)=U11(IJ,M1,MU)
6 CONTINUE
end do
M2=IUP(M1,MU)
!     
do NU=MU+1,4
IPLAQ=6-NU-MU+5*(NU/4)
!     
do IJ=1,NCOL2
B11(IJ)=U11(IJ,M2,NU)
7 CONTINUE
end do
CALL VMX(1,A11,B11,C11,1)
M3=IUP(M1,NU)
do IJ=1,NCOL2
D11(IJ)=U11(IJ,M3,MU)
8 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
CALL VMX(1,C11,D11,B11,1)
do IJ=1,NCOL2
D11(IJ)=U11(IJ,M1,NU)
9 CONTINUE
end do
CALL HERM(1,D11,DUM11,1)
CALL TRVMX(1,B11,D11,ACT,1)
ANN=1.0/NCOL
ACT=ANN*REAL(ACT)
ACTN(JBIN,IPLAQ)=ACTN(JBIN,IPLAQ)+ACT
IF(NU.NE.4)THEN
VACTS=VACTS+ACT
ELSE
VACTT=VACTT+ACT
end if
!     
11 CONTINUE
end do
end do
end do
!
VACTS=VACTS/(3.0*LSIZE)
VACTT=VACTT/(3.0*LSIZE)
TOTACT = 0.5*(VACTS+VACTT)
write (91,*) VACTS, VACTT, TOTACT
!
AVACS(JBIN)=AVACS(JBIN)+VACTS
AVACT(JBIN)=AVACT(JBIN)+VACTT
!
IF(ITER.EQ.NITER)THEN
DO IP=1,6
DO JB=1,NUMBIN
ACTN(JB,IP)=ACTN(JB,IP)/(IBIN*LSIZE)
end do
end do
end if
!
end if
!
IF(IPR.EQ.1)THEN
!
WRITE(6,80)
80 FORMAT(' *****************************************************')
WRITE(6,100)
100 FORMAT('                 ACTION            ')
WRITE(6,80)
WRITE(6,81)
81 FORMAT(' *  ')
do NB=1,NUMBIN
AV(NB)=AVACS(NB)/IBIN
20 CONTINUE
end do
CALL JACKM(NUMBIN,AV,AVRS,ERS)
do NB=1,NUMBIN
AV(NB)=AVACT(NB)/IBIN
22 CONTINUE
end do
CALL JACKM(NUMBIN,AV,AVRT,ERT)
WRITE(6,110) AVRS,ERS
110 FORMAT('   AVER,ERR   SPACE ACTION =',2F12.8)
WRITE(6,112) AVRT,ERT
112 FORMAT('   AVER,ERR   TIME  ACTION =',2F12.8)
WRITE(6,81)
!
end if
!
RETURN
end subroutine ACTION
!***********************************************************************
!***********************************************************************
!                          THERMAL LINES
!***********************************************************************
!***********************************************************************
SUBROUTINE POLY(IPR,ITER)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(NITER=2,NUMBIN=2)
PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
COMMON/POLYT/TLINE(NITER),SQTLINE(NITER)
COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
COMMON/COMMUNICATE2/AVAC(NUMBIN)
COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2)
DIMENSION AV1(NUMBIN),AV2(NUMBIN)
!
complex(real32) U11,A11,B11,C11,AKT1,TLINE,AVAC,CSUM
!
IBIN=NITER/NUMBIN
IF(IPR.EQ.1)go to 100
!
JBIN=(ITER-1)/IBIN+1
!
IF(ITER.EQ.1)THEN
do NB=1,NUMBIN
AVAC(NB)=(0.0,0.0)
3 CONTINUE
end do
end if
!
LX12=LX1*LX2
LX123=LX1*LX2*LX3
!
AKT1=(0.0,0.0)
AKT2=0.0
DO N4=1,LX4
DO N3=1,LX3
DO N2=1,LX2
M2=1+LX1*(N2-1)+LX12*(N3-1)+LX123*(N4-1)
!
do IC=1,NCOL2
A11(IC)=(0.0,0.0)
2 CONTINUE
end do
do NC=1,NCOL
IC=NC+NCOL*(NC-1)
A11(IC)=(1.0,0.0)
4 CONTINUE
end do
!
do N1=1,LX1
do IC=1,NCOL2
B11(IC)=U11(IC,M2,1)
6 CONTINUE
end do
CALL VMX(1,A11,B11,C11,1)
M3=IUP(M2,1)
M2=M3
do IC=1,NCOL2
A11(IC)=C11(IC)
8 CONTINUE
end do
10 CONTINUE
end do
!
CSUM=(0.0,0.0)
do NC=1,NCOL
IC=NC+NCOL*(NC-1)
CSUM=CSUM+A11(IC)
12 CONTINUE
end do
CSUM=CSUM/NCOL
AKT1=AKT1+CSUM
AKT2=AKT2+real(CSUM*CONJG(CSUM))
!
end do
end do
end do
!
AKT1=AKT1/(LX2*LX3)
AKT2=AKT2/(LX2*LX3)
TLINE(ITER)=AKT1
SQTLINE(ITER)=AKT2
AVAC(JBIN)=AVAC(JBIN)+AKT1
!
go to 200
100 CONTINUE
!
WRITE(6,80)
80 FORMAT(' *****************************************************')
WRITE(6,108)
108 FORMAT('     THERMAL LINES           ')
WRITE(6,80)
WRITE(6,81)
81 FORMAT(' *  ')
do NB=1,NUMBIN
AV1(NB)=real(AVAC(NB))/IBIN
AV2(NB)=aimag(AVAC(NB))/IBIN
WRITE(6,109) NB,AV1(NB),AV2(NB)
109 FORMAT(' BIN =',I4,'   REAL,IMAG  POLY =',2F9.4)
20 CONTINUE
end do
CALL JACKM(NUMBIN,AV1,AVR,ERR)
CALL JACKM(NUMBIN,AV2,AVI,ERI)
WRITE(6,81)
WRITE(6,110) AVR,ERR
110 FORMAT('   AV,ER  REAL POLY =',2F12.8)
WRITE(6,111) AVI,ERI
111 FORMAT('   AV,ER  IMAG POLY =',2F12.8)
WRITE(6,81)
!
200 RETURN
end subroutine POLY
!**********************************************************
!  THIS ROUTINE REIMPOSES THE UNITARITY CONSTRAINTS ON
!  OUR SU(3) MATRICES
!**********************************************************
SUBROUTINE RENORM
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
COMMON/ARRAYS/U11(NCOL,NCOL,LSIZE,4)
complex(real32) ADUM(NCOL,NCOL),U11,CSUM
!
complex(real32) aa(ncol,ncol)
complex(real32) det,cnorm,cdet,dncol
!
DO MU=1,4
DO NN=1,LSIZE
!
DO J=1,NCOL
DO I=1,NCOL
ADUM(I,J)=U11(I,J,NN,MU)
end do
end do
!
do N2=1,NCOL
!
do N3=1,N2-1
!
CSUM=(0.0,0.0)
do N1=1,NCOL
CSUM=CSUM+ADUM(N1,N2)*CONJG(ADUM(N1,N3))
5 CONTINUE
end do
do N1=1,NCOL
ADUM(N1,N2)=ADUM(N1,N2)-CSUM*ADUM(N1,N3)
6 CONTINUE
end do
!     
10 CONTINUE
end do
!
SUM=0.0
do N1=1,NCOL
SUM=SUM+ADUM(N1,N2)*CONJG(ADUM(N1,N2))
7 CONTINUE
end do
ANORM=1.0/SQRT(SUM)
do N1=1,NCOL
ADUM(N1,N2)=ADUM(N1,N2)*ANORM
8 CONTINUE
end do
!
20 CONTINUE
end do
!
DO J=1,NCOL
DO I=1,NCOL
U11(I,J,NN,MU)=ADUM(I,J)
end do
end do
!
do j=1,ncol
do i=1,ncol
aa(i,j)=adum(i,j)
end do
end do
jmat=ncol
CALL DETNANT(jmat,det,aa)
!
dncol=cmplx(1.0/ncol)
cdet=det**dncol
cnorm=1.0/cdet
do j=1,ncol
do i=1,ncol
aa(i,j)=adum(i,j)*cnorm
u11(i,j,nn,mu)=aa(i,j)
end do
end do
!
end do
end do
RETURN
end subroutine RENORM
!********************************************************
!********************************************************
SUBROUTINE DETNANT(NUMOP,DET,AA)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
!
complex(real32) AA(NUMOP,NUMOP),B(NUMOP,NUMOP),det
!
NP=NUMOP
!
DO I=1,NP
DO J=1,NP
B(I,J)=AA(I,J)
end do
end do
CALL DETMAT(B,NP,DET)
!     
RETURN
end subroutine DETNANT
!***************************************************************
SUBROUTINE DETMAT(A,NP,DET)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
!     
DIMENSION INDX(NP)
complex(real32) A(NP,NP),DET,DD
!     
N=NP
CALL LUDCMP(A,N,NP,INDX,D)
DD=D
do J=1,NP
DD=DD*A(J,J)
13 CONTINUE
end do
DET=DD
!     
RETURN
end subroutine DETMAT
!***************************************************************
!***************************************************************
SUBROUTINE LUDCMP(A,N,NP,INDX,D)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(NCOL=2)
PARAMETER(NMAX=NCOL,TINY=1.0E-20)
!     
complex(real32) A(NP,NP),SUM,CDUM
DIMENSION INDX(N),VV(NMAX)
!     
D=1.
do I=1,N
AAMAX=0.
do J=1,N
IF(CABS(A(I,J)).GT.AAMAX) AAMAX=CABS(A(I,J))
11 CONTINUE
end do
!     IF(AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX IN LUDCMP'
IF(AAMAX.EQ.0.) return
VV(I)=1./AAMAX
12 CONTINUE
end do
!     
do J=1,N
do I=1,J-1
SUM=A(I,J)
do K=1,I-1
SUM=SUM-A(I,K)*A(K,J)
13 CONTINUE
end do
A(I,J)=SUM
14 CONTINUE
end do
!     
AAMAX=0.
do I=J,N
SUM=A(I,J)
do K=1,J-1
SUM=SUM-A(I,K)*A(K,J)
15 CONTINUE
end do
A(I,J)=SUM
DUM=VV(I)*CABS(SUM)
IF(DUM.GE.AAMAX)THEN
IMAX=I
AAMAX=DUM
end if
16 CONTINUE
end do
!     
IF(J.NE.IMAX)THEN
do K=1,N
CDUM=A(IMAX,K)
A(IMAX,K)=A(J,K)
A(J,K)=CDUM
17 CONTINUE
end do
D=-D
VV(IMAX)=VV(J)
end if
!     
INDX(J)=IMAX
IF(CABS(A(J,J)).EQ.0.) A(J,J)=TINY
IF(J.NE.N)THEN
CDUM=1./A(J,J)
do I=J+1,N
A(I,J)=A(I,J)*CDUM
18 CONTINUE
end do
end if
!     
19 CONTINUE
end do
!     
RETURN
end subroutine LUDCMP
!****************************************************************
!****************************************************************
!     HERE WE SET UP  LINK POINTERS FOR FULL LATTICE            *
!****************************************************************
!****************************************************************
SUBROUTINE SETUP
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(LX1=26,LX2=26,LX3=26,LX4=52)
PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
!
COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
!
LX12=LX1*LX2
LX123=LX12*LX3
LX1234=LX123*LX4
!
NN=0
DO L4=1,LX4
DO L3=1,LX3
DO L2=1,LX2
DO L1=1,LX1
NN=NN+1
!     
NU1=NN+1
IF((L1+1).GT.LX1) NU1=NU1-LX1
IUP(NN,1)=NU1
ND1=NN-1
IF((L1-1).LT.1) ND1=ND1+LX1
IDN(NN,1)=ND1
!     
NU2=NN+LX1
IF((L2+1).GT.LX2) NU2=NU2-LX12
IUP(NN,2)=NU2
ND2=NN-LX1
IF((L2-1).LT.1) ND2=ND2+LX12
IDN(NN,2)=ND2
!     
NU3=NN+LX12
IF((L3+1).GT.LX3) NU3=NU3-LX123
IUP(NN,3)=NU3
ND3=NN-LX12
IF((L3-1).LT.1) ND3=ND3+LX123
IDN(NN,3)=ND3
!     
NU4=NN+LX123
IF((L4+1).GT.LX4) NU4=NU4-LX1234
IUP(NN,4)=NU4
ND4=NN-LX123
IF((L4-1).LT.1) ND4=ND4+LX1234
IDN(NN,4)=ND4
!     
end do
end do
end do
end do
!
RETURN
end subroutine SETUP
!***********************************************************************
!                  jack-knife errors for av(i=1,..,num)                *
!***********************************************************************
SUBROUTINE JACKM(NUM,AV,AVER,ERR)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
!
DIMENSION AV(*)
!
SUM=0.0d0
do N=1,NUM
SUM=SUM+AV(N)
10 CONTINUE
end do
SUM=SUM/NUM
!
ESUM=0.0d0
do N=1,NUM
DIF=SUM-AV(N)
ESUM=ESUM+DIF*DIF
12 CONTINUE
end do
ESUM=ESUM/NUM
!
AVER=SUM
ERR=SQRT(ESUM*NUM)/(NUM-1)
!
RETURN
end subroutine JACKM
!***********************************************************************
!                 jack-knife errors for ratio avu to avd               *
!***********************************************************************
SUBROUTINE JACKMR(NUM,NMAX,AVU,AVD,AVER,ERR)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
!
DIMENSION AVU(NMAX),AVD(NMAX)
DIMENSION DIFU(NMAX),DIFD(NMAX)
!
AVER=0.0d0
ERR=0.0d0
!
SUMU=0.0d0
SUMD=0.0d0
do N=1,NUM
SUMU=SUMU+AVU(N)
SUMD=SUMD+AVD(N)
10 CONTINUE
end do
do N=1,NUM
DIFU(N)=SUMU-AVU(N)
DIFD(N)=SUMD-AVD(N)
12 CONTINUE
end do
do N=1,NUM
IF(DIFD(N).EQ.0.0d0)THEN
ERR=99.0d0
go to 99
end if
14 CONTINUE
end do
!
ASUM=0.0d0
ESUM=0.0d0
do N=1,NUM
DIF=DIFU(N)/DIFD(N)
ASUM=ASUM+DIF
ESUM=ESUM+DIF*DIF
16 CONTINUE
end do
ASUM=ASUM/NUM
ESUM=ESUM/NUM
ESUM=ESUM-ASUM*ASUM
!
IF(SUMD.NE.0.0) AVER=SUMU/SUMD
IF(ESUM.GT.0.0) ERR=SQRT(ESUM*NUM)
!
99 CONTINUE
RETURN
end subroutine JACKMR
!***********************************************************************
!                           fitting masses                             *
!***********************************************************************
SUBROUTINE JACK(NBIN,ISUB,LXI,NUMBIN,LMAX,COR,VAC)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
!
COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
!
DIMENSION COR(NUMBIN,LMAX),VAC(NUMBIN)
DIMENSION CORR(200),AW(200)
!
do I4=1,LMAX
ACOR(I4)=0.0d0
SCOR(I4)=0.0d0
AMM(I4)=0.0d0
SMM(I4)=0.0d0
1 CONTINUE
end do
!
do IB=1,NBIN
!
IF(ISUB.EQ.0)THEN
VACC=0.0d0
do JB=1,NBIN
IF(JB.EQ.IB)go to 4
VACC=VACC+VAC(JB)
4 CONTINUE
end do
end if
!
do NT=1,LMAX
CORR(NT)=0.0d0
do JB=1,NBIN
IF(JB.EQ.IB)go to 7
CORR(NT)=CORR(NT)+COR(JB,NT)
7 CONTINUE
end do
IF(ISUB.EQ.0) CORR(NT)=CORR(NT)-VACC*VACC/(NBIN-1)
6 CONTINUE
end do
ANORM=CORR(1)
IF(ANORM.EQ.0.0)go to 99
do NT=1,LMAX
CORR(NT)=CORR(NT)/ANORM
8 CONTINUE
end do
CALL FITT(LXI,CORR,AW,NTMAX)
do NT=1,LMAX
ACOR(NT)=ACOR(NT)+CORR(NT)
SCOR(NT)=SCOR(NT)+CORR(NT)**2
AMM(NT)=AMM(NT)+AW(NT)
SMM(NT)=SMM(NT)+AW(NT)**2
9 CONTINUE
end do
2 CONTINUE
end do
!  C
do NT=1,LMAX
ACOR(NT)=ACOR(NT)/NBIN
SCOR(NT)=SCOR(NT)/NBIN
SCOR(NT)=(SCOR(NT)-ACOR(NT)*ACOR(NT))*NBIN
IF(SCOR(NT).GT.0.0) SCOR(NT)=SQRT(SCOR(NT))
AMM(NT)=AMM(NT)/NBIN
SMM(NT)=SMM(NT)/NBIN
SMM(NT)=(SMM(NT)-AMM(NT)*AMM(NT))*NBIN
IF(SMM(NT).GT.0.0) SMM(NT)=SQRT(SMM(NT))
10 CONTINUE
end do
!     
IF(ISUB.EQ.0)THEN
VACC=0.0d0
do JB=1,NBIN
VACC=VACC+VAC(JB)
22 CONTINUE
end do
end if
do NT=1,LMAX
CORR(NT)=0.0d0
do JB=1,NBIN
CORR(NT)=CORR(NT)+COR(JB,NT)
26 CONTINUE
end do
IF(ISUB.EQ.0) CORR(NT)=CORR(NT)-VACC*VACC/NBIN
24 CONTINUE
end do
ANORM=CORR(1)
IF(ANORM.EQ.0.0)go to 99
do NT=1,LMAX
CORR(NT)=CORR(NT)/ANORM
28 CONTINUE
end do
CALL FITT(LXI,CORR,AW,NTMAX)
do NT=1,LMAX
ACOR(NT)=CORR(NT)
AMM(NT)=AW(NT)
29 CONTINUE
end do
!     C
99 RETURN
end subroutine JACK
!*********************************************************************
! effective energies using a local cosh fit - lattice length LT.
! AW(timedif+1) is input corrln function;
! BW(nt) is eff energy from  time differences nt-1 to nt.
!*********************************************************************
SUBROUTINE FITT(LT,AW,BW,NTMAX)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
!
DIMENSION AW(200),BW(200)
!
do I=1,LT/2+1
BW(I)=0.0d0
1 CONTINUE
end do
!
NTMAX=LT/2+1
do I4=1,LT/2
!  C
IF(AW(I4).LE.0.000001.OR.AW(I4+1).LE.0.000001)THEN
NTMAX=I4-1
go to 99
end if
AMSM=0.0d0
FTM=(AW(I4)/AW(I4+1))
IF(FTM.GT.1.0)THEN
AML=DLOG(FTM)
AMU=DLOG(2.0d0*FTM)
do NS=1,20
AMS=(AML+AMU)/2
FTS=(EXP(-AMS*(I4-1))+EXP(-(LT-I4+1)*AMS)) &
  & /(EXP(-AMS*(I4))+EXP(-(LT-I4)*AMS))
IF(FTS.LT.FTM)THEN
AML=AMS
ELSE
AMU=AMS
end if
13 CONTINUE
end do
AMSM=(AML+AMU)/2
ELSE
NTMAX=I4-1
go to 99
end if
BW(I4)=AMSM
!     
2 CONTINUE
end do
!
99 CONTINUE
RETURN
end subroutine FITT
!***********************************************************************
!                 VECTOR MATRIX MULTIPLY ... 5*5 COMPLEX               *
!***********************************************************************
SUBROUTINE VMX(NNN1,A,B,C,NNN2)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
complex(real32) A(NCOL,NCOL),B(NCOL,NCOL),C(NCOL,NCOL),CSUM
!
DO J=1,NCOL
DO I=1,NCOL
CSUM=(0.0,0.0)
do K=1,NCOL
CSUM=CSUM+A(I,K)*B(K,J)
2 CONTINUE
end do
C(I,J)=CSUM
end do
end do
!
RETURN
end subroutine VMX
!***********************************************************************
!***********************************************************************
!                   TRACE PRODUCT .. Ncol*Ncol COMPLEX                 *
!***********************************************************************
!***********************************************************************
SUBROUTINE TRVMX(NNN1,A,B,CC,NNN2)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
complex(real32) A(NCOL,NCOL),B(NCOL,NCOL),CC
!
CC=(0.0,0.0)
DO I=1,NCOL
DO K=1,NCOL
CC=CC+A(I,K)*B(K,I)
end do
end do
!
RETURN
end subroutine TRVMX
!***********************************************************************
!***********************************************************************
!                           HERMITIAN CONJUGATE                        *
!***********************************************************************
!***********************************************************************
SUBROUTINE HERM(NNN1,A11,DUM11,NNN2)
use iso_fortran_env, only : real32, real64, int32
implicit real(real64) (A-H,O-Z)
PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
!
complex(real32) A11(NCOL,NCOL),DUM11(NCOL,NCOL)
!
DO I=1,NCOL
DO J=1,NCOL
DUM11(I,J)=CONJG(A11(J,I))
end do
end do
DO I=1,NCOL
DO J=1,NCOL
A11(I,J)=DUM11(I,J)
end do
end do
!
RETURN
end subroutine HERM
