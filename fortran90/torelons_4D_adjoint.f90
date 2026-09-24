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
PARAMETER(ICMIN=134764,ICMAX=NITER+ICMIN-1)
PARAMETER(IBLOK=5,IBING=109,NUMBIN=2)
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
complex(real32) :: U11
real(real64) :: rndnum
!
character(len=43) :: homepath
character(len=24) :: conf_directory1
character(len=35) :: conf_directory2
character(len=102) :: directory
character(len=120) :: file_name
character(len=256) :: copy_file
character(len=256) :: trnsf_file
character(len=7) :: confnum
character(len=6) :: name_of_file
!
!homepath=trim('/home/dp208/dp208/dc-athe1/AXIONS/NF2/b2.3/')
!conf_directory1=trim('m-1.0/26x26x26x52/confs/')
!conf_directory2=trim('run1_52x26x26x26nc2rADJnf2b2.300000')
!directory=trim(homepath//conf_directory1//conf_directory2)
homepath="/gpfs/scratch/ehpc598/torelons/"
conf_directory1="cnfg/"
conf_directory2="run1_52x26x26x26nc2rADJnf2b2.300000"
directory=trim(homepath) // trim(conf_directory1) //trim(conf_directory2)
!
CALL SETUP
ISEED=3591
call rluxgo(3, iseed, 0, 0)
!
ITOT=NITER
!
WRITE(*,*) "                                                "
WRITE(6,90)
90 FORMAT(" *******************************************************")
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
91 FORMAT(" *")
WRITE(6,92) LX1,LX2,LX3,LX4
92 FORMAT("[Info][Lattice Size]  ","             LX   = ",4I4)
WRITE(6,193) BETAG
193 FORMAT("[Info][Value of Beta] ","           beta   =",F8.4)
WRITE(6,194) NCOL
194 FORMAT("[Info][Colors]        "," Number of Colors = ",I6)
WRITE(6,195) IHEAT
195 FORMAT("[Info][Thermalisation]","   Therml. sweeps = ",I6)
WRITE(6,91)
WRITE(6,90)
WRITE(6,91)
WRITE(6,196) ITOT
196 FORMAT("[Info][Measurements]","   Number of MC Iterations = ",I6)
WRITE(6,197) ICALLG
197 FORMAT("[Info][Measurements]","   Sweeps per measurements = ",I6)
WRITE(6,198) IBING
198 FORMAT("[Info][Measurements]","      Measurements per bin = ",I6)
WRITE(6,91)
WRITE(6,90)
200 FORMAT("[Info][Number of Operators q=0]","  Total Number = ",I4)
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
205 FORMAT("[Info][Number of Operators q=0]","     J   P   R      #")
206 FORMAT("[Info][Number of Operators q=0]","     ----------------")
210 FORMAT("[Info][Number of Operators q=0]","     0   +   +   ",I4)
211 FORMAT("[Info][Number of Operators q=0]","     0   +   -   ",I4)
212 FORMAT("[Info][Number of Operators q=0]","     0   -   +   ",I4)
213 FORMAT("[Info][Number of Operators q=0]","     0   -   -   ",I4)
214 FORMAT("[Info][Number of Operators q=0]","     1   #   +   ",I4)
216 FORMAT("[Info][Number of Operators q=0]","     1   #   -   ",I4)
218 FORMAT("[Info][Number of Operators q=0]","     2   +   +   ",I4)
219 FORMAT("[Info][Number of Operators q=0]","     2   +   -   ",I4)
220 FORMAT("[Info][Number of Operators q=0]","     2   -   +   ",I4)
221 FORMAT("[Info][Number of Operators q=0]","     2   -   -   ",I4)
207 FORMAT("[Info][Number of Operators q=0]","     ----------------")
WRITE(6,91)
WRITE(6,90)
300 FORMAT("[Info][Number of Operators q=1,2]","  Total Number = ",I4)
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
305 FORMAT("[Info][Number of Operators q=1,2]","     J   P   #")
306 FORMAT("[Info][Number of Operators q=1,2]","     ---------------")
310 FORMAT("[Info][Number of Operators q=1,2]","     0   + ",I4)
312 FORMAT("[Info][Number of Operators q=1,2]","     0   - ",I4)
314 FORMAT("[Info][Number of Operators q=1,2]","     1     ",I4)
318 FORMAT("[Info][Number of Operators q=1,2]","     2   + ",I4)
320 FORMAT("[Info][Number of Operators q=1,2]","     2   - ",I4)
307 FORMAT("[Info][Number of Operators q=1,2]","     ---------------")
WRITE(6,91)
WRITE(6,90)
!
IFILE = ICMIN
do ITER=1,NITER
IFILE = IFILE + 16
!
ist=IFILE
write(confnum, "(i0)") ist
file_name=trim(directory)//"m1.000000n"//trim(confnum)
!
write(6,*) "------------------------------"
write(6,*) "Execution of external commands"
write(6,*) "                              "
write(6,*) "Commands to be executed"
write(6,*) "------------------------------"

copy_file=trim("cp "//file_name//" ./conf")
write(6,*) copy_file
write(6,*) "                              "
write(6,*) "                              "
write(6,*) "*************************************************"
call execute_command_line(copy_file, WAIT=.true.)
write(6,*) "*************************************************"
write(6,*) "                              "
write(6,*) "                              "
name_of_file=trim("./conf")
WRITE(6,*) "Config to be read: ", file_name
!
CALL cpu_time(t1)
CALL READ_GF(name_of_file)
CALL cpu_time(t2)
WRITE(6,902) real(t2-t1, kind=real32)
902 FORMAT("[Info][Time]","    Configuration Read time:", F8.4)
!
CALL ACTION(0,ITER,TOTACT)
511 format("[FM][0]Check plaq = ", f8.6)
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
integer(int32) :: nc_read, nx_read, ny_read, nz_read, nt_read
integer :: t, x, y, z, dir, dir_target, iun, iq
real(real64) :: plaquette_read

real(real64) :: quaternion(4)
complex(real32), dimension(2, 2) :: quaternion_to_matrix

COMMON/ARRAYS/U11(NCOL,NCOL,LX1,LX2,LX3,LX4,4)

complex(real32) :: U11

open(newunit=iun, file=filename, access="stream", form="unformatted", &
     status="old", action="read")
nc_read      = read_be_int32(iun)
nt_read      = read_be_int32(iun)
nx_read      = read_be_int32(iun)
ny_read      = read_be_int32(iun)
nz_read      = read_be_int32(iun)
plaquette_read = read_be_real64(iun)

WRITE(6,101) plaquette_read
101 FORMAT("[I/O][Plaq]", "Plaquette value:",F8.6)
WRITE(6,102) nc_read
102 FORMAT("[I/O][Ncol]", "Number of Colors:",I2.1)
WRITE(6,103) nt_read, nx_read, ny_read, nz_read
103 FORMAT("[I/O][Dim]", "T x X x Y x Z=",I3.2, I3.2, I3.2, I3.2)
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

complex(real32) :: U11
complex(real64) :: UR11(LX4,LX3,LX2,LX1,4,NCOL,NCOL)
character(len=6) :: filename
integer :: iun
real(real64) :: rpart, ipart
!
open(newunit=iun, file=trim(filename), access="stream", &
     form="unformatted", status="old", action="read")
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
WRITE(6,*) "Configuration read"
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
complex(real32) :: U11,A11,B11,C11,D11,DUM11,ACT
!
DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN),AVACS(NUMBIN),AVACT(NUMBIN)
DIMENSION VAL(NUMBIN),AV(NUMBIN)
dimension icoordvect(4)
complex(real32) :: cpol1
complex(real64) :: cpol(4)
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
888 format("Poly(" ,i1, ")  =  (",f16.12,",",f16.12,")")
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
complex(real64) :: ACORLP
complex(real64) :: AVACLP
complex(real64) :: ACORLJ0PP
complex(real64) :: AVACLJ0PP
complex(real64) :: ACORLJ0PM
complex(real64) :: AVACLJ0PM
complex(real64) :: ACORLJ0MP
complex(real64) :: AVACLJ0MP
complex(real64) :: ACORLJ0MM
complex(real64) :: AVACLJ0MM
complex(real64) :: ACORLJ1P
complex(real64) :: AVACLJ1P
complex(real64) :: ACORLJ1M
complex(real64) :: AVACLJ1M
complex(real64) :: ACORLJ2PP
complex(real64) :: AVACLJ2PP
complex(real64) :: ACORLJ2PM
complex(real64) :: AVACLJ2PM
complex(real64) :: ACORLJ2MP
complex(real64) :: AVACLJ2MP
complex(real64) :: ACORLJ2MM
complex(real64) :: AVACLJ2MM
!
complex(real64) :: ACORLJ0
complex(real64) :: AVACLJ0
complex(real64) :: ACORLJ1
complex(real64) :: AVACLJ1
complex(real64) :: ACORLJ2
complex(real64) :: AVACLJ2
complex(real64) :: ACORLALL
complex(real64) :: AVACLALL
!
complex(real32) :: TLINE
!
OPEN (11,FILE="ACORLQ0.DAT")
OPEN (12,FILE="ACORLJ0PPQ0.DAT")
OPEN (13,FILE="ACORLJ0PMQ0.DAT")
OPEN (14,FILE="ACORLJ0MPQ0.DAT")
OPEN (15,FILE="ACORLJ0MMQ0.DAT")
OPEN (16,FILE="ACORLJ1PQ0.DAT")
OPEN (17,FILE="ACORLJ1MQ0.DAT")
OPEN (18,FILE="ACORLJ2PPQ0.DAT")
OPEN (19,FILE="ACORLJ2PMQ0.DAT")
OPEN (20,FILE="ACORLJ2MPQ0.DAT")
OPEN (21,FILE="ACORLJ2MMQ0.DAT")
OPEN (22,FILE="ACORLJ0Q0.DAT")
OPEN (23,FILE="ACORLJ1Q0.DAT")
OPEN (24,FILE="ACORLJ2Q0.DAT")
OPEN (25,FILE="ACORLJALL.DAT")
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
complex(real64) :: ACORLPMOM
complex(real64) :: AVACLPMOM
complex(real64) :: ACORLPMOMJ0P
complex(real64) :: AVACLMOMJ0P
complex(real64) :: ACORLPMOMJ0M
complex(real64) :: AVACLMOMJ0M
complex(real64) :: ACORLPMOMJ1
complex(real64) :: AVACLMOMJ1
complex(real64) :: ACORLPMOMJ2P
complex(real64) :: AVACLMOMJ2P
complex(real64) :: ACORLPMOMJ2M
complex(real64) :: AVACLMOMJ2M
complex(real64) :: ACORLPMOMJ0
complex(real64) :: AVACLMOMJ0
complex(real64) :: ACORLPMOMJ2
complex(real64) :: AVACLMOMJ2
complex(real64) :: ACORLPMOMJALL
complex(real64) :: AVACLMOMALL
!*****************************************************************
OPEN (17, FILE="ACORLPMOM1.DAT")
OPEN (18, FILE="ACORLPMOM2.DAT")
OPEN (19, FILE="ACORLJ0PMOMQ1.DAT")
OPEN (20, FILE="ACORLJ0PMOMQ2.DAT")
OPEN (21, FILE="ACORLJ0MMOMQ1.DAT")
OPEN (22, FILE="ACORLJ0MMOMQ2.DAT")
OPEN (23, FILE="ACORLJ1MOMQ1.DAT")
OPEN (24, FILE="ACORLJ1MOMQ2.DAT")
OPEN (25, FILE="ACORLJ2PMOMQ1.DAT")
OPEN (26, FILE="ACORLJ2PMOMQ2.DAT")
OPEN (27, FILE="ACORLJ2MMOMQ1.DAT")
OPEN (28, FILE="ACORLJ2MMOMQ2.DAT")
OPEN (29, FILE="ACORLJ0MOMQ1.DAT")
OPEN (30, FILE="ACORLJ0MOMQ2.DAT")
OPEN (33, FILE="ACORLJ2MOMQ1.DAT")
OPEN (34, FILE="ACORLJ2MOMQ2.DAT")
OPEN (35, FILE="ACORLJALLMOMQ1.DAT")
OPEN (36, FILE="ACORLJALLMOMQ2.DAT")
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
IF(ITER==(ITERG*ICALLG).AND.ITERG<=NTOTG) THEN
PAR=1.0d0
!
CALL SETUPB
!
!        WRITE(*,*) "CALLING BLOCK"
CALL BLOCK
WRITE(*,*) "Number of iteration:", ITER
end if
!
IF(ITER==NTOTG*ICALLG) THEN
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
complex(real64) :: ACORLP,AVACLP
!
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE="DIAGNALFULLQ0.DAT")
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
80 FORMAT("****************************************************")
WRITE(23,181)
181 FORMAT("****************** OPERATORS WITH J = 0 ************")
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(" *")
WRITE(23,91)
91 FORMAT(" AVERAGE  X,Y  LINES ")
WRITE(23,81)
82 FORMAT(" ******************* ")
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
101 FORMAT("OP =",I4,"  AV,ER R,I LINES =",4F9.4)
10 CONTINUE
end do
!
WRITE(23,80)
WRITE(23,81)
WRITE(23,92)
92 FORMAT(" AVERAGE   Y   LINES ")
WRITE(23,81)

do ID=1, NTOTAL
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(" OPERATOR = ",I4)
WRITE(23,82)
DO NT=1,LMAXIR
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLP(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1))<=EPS)goto 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAXIR,CORP,VAC)
!         DO 26 NT=1,MAXDTLS
do NT=1,LMAXIR-1
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT("  DT=",I3,"   AV,ER COR = ",2F8.4,"    E=",2F8.4)
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
complex(real64) :: ACORLPMOM,AVACLPMOM
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (24, FILE="DIAGNALFULLMOMQ0.DAT")
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
80 FORMAT("****************************************************")
181 FORMAT("************ OPERATORS ************")
WRITE(24,80)
WRITE(24,80)
WRITE(24,81)
WRITE(24,80)
WRITE(24,181)
WRITE(24,80)
WRITE(24,81)
81 FORMAT(" *")
WRITE(24,91)
91 FORMAT(" AVERAGE  X,Y  LINES ")
WRITE(24,81)
82 FORMAT(" ******************* ")
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
101 FORMAT("OP =",I4,"  AV,ER R,I LINES =",4F9.4)
10 CONTINUE
end do
!
WRITE(24,80)
WRITE(24,81)
WRITE(24,92)
92 FORMAT(" AVERAGE   Y   LINES ")
WRITE(24,81)

do ID=1, NTOTALMOM
WRITE(24,81)
WRITE(24,82)
WRITE(24,94) ID
94 FORMAT(" OPERATOR = ",I4)
WRITE(24,82)
DO NT=1,LMAXIR
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOM(IJ,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1))<=EPS)goto 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAXIR,CORP,VAC)
!
!         DO 26 NT=1,MAXDTLS
do NT=1, LMAXIR-1
WRITE(24,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT("  DT=",I3,"   AV,ER COR = ",2F8.4,"    E=",2F8.4)
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
complex(real64) :: ACORLP,AVACLP
!***********************************************************
complex(real64) :: ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
complex(real64) :: AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
!***********************************************************
complex(real64) :: ACORLJ1P,ACORLJ1M
complex(real64) :: AVACLJ1P,AVACLJ1M
!***********************************************************
complex(real64) :: ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
complex(real64) :: AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
!***********************************************************
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE="DIAGNALQ0IND.DAT")
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
80 FORMAT("****************************************************")
WRITE(23,181)
181 FORMAT("****************** OPERATORS WITH Q = 0 ************")
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(" *")
WRITE(23,91)
91 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=+, Pr=+ ")
WRITE(23,81)
82 FORMAT(" ******************* ")
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
101 FORMAT("OP =",I4,"  AV,ER R,I LINES =",4F9.4)
10 CONTINUE
end do
WRITE(23,81)
WRITE(23,201)
201 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=+, Pr=- ")
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
202 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=-, Pr=+ ")
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
203 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=-, Pr=- ")
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
204 FORMAT(" AVERAGE  X,Y  LINES, J=1, Pr=+ ")
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
205 FORMAT(" AVERAGE  X,Y  LINES, J=1, Pr=- ")
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
206 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=+, Pr=+ ")
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
208 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=+, Pr=- ")
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
209 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=-, Pr=+ ")
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
210 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=-, Pr=- ")
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
220 FORMAT(" DIAGONAL TORELONS J=0, Pp=+, Pr=+, q=0 ")
WRITE(23,81)

do ID=1, NOPJ0PP
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(" OPERATOR = ",I4)
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLJ0PP(IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1))<=EPS)goto 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT("  DT=",I3,"   AV,ER COR = ",2F8.4,"    E=",2F8.4)
end do
20 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,211)
211 FORMAT(" DIAGONAL TORELONS J=0, Pp=+, Pr=-, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 21
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
21 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,212)
212 FORMAT(" DIAGONAL TORELONS J=0, Pp=-, Pr=+, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 22
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
22 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,213)
213 FORMAT(" DIAGONAL TORELONS J=0, Pp=-, Pr=-, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 23
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
23 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,214)
214 FORMAT(" DIAGONAL TORELONS J=1, Pr=+, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 24
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
24 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,215)
215 FORMAT(" DIAGONAL TORELONS J=1, Pr=-, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 25
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
25 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,216)
216 FORMAT(" DIAGONAL TORELONS J=2, Pp=+, Pr=+, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 26
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
26 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,217)
217 FORMAT(" DIAGONAL TORELONS J=2, Pp=+, Pr=-, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 27
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
27 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,218)
218 FORMAT(" DIAGONAL TORELONS J=2, Pp=-, Pr=+, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 28
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
28 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,219)
219 FORMAT(" DIAGONAL TORELONS J=2, Pp=-, Pr=-, q=0 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 29
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
complex(real64) :: ACORLP,AVACLP
!***********************************************************
complex(real64) :: ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
complex(real64) :: AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
!***********************************************************
complex(real64) :: ACORLJ1P,ACORLJ1M
complex(real64) :: AVACLJ1P,AVACLJ1M
!***********************************************************
complex(real64) :: ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
complex(real64) :: AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
!***********************************************************
complex(real64) :: ACORLPMOMJ0P, ACORLPMOMJ0M
complex(real64) :: AVACLMOMJ0P, AVACLMOMJ0M
!***********************************************************
complex(real64) :: ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
complex(real64) :: AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
!***********************************************************
complex(real64) :: ACORLPMOMJ2P, ACORLPMOMJ2M
complex(real64) :: AVACLMOMJ2P, AVACLMOMJ2M
!***********************************************************
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE="DIAGNALQ1IND.DAT")
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
80 FORMAT("****************************************************")
WRITE(23,181)
181 FORMAT("****************** OPERATORS WITH Q = 1 ************")
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(" *")
WRITE(23,91)
91 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=+")
WRITE(23,81)
82 FORMAT(" ******************* ")
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
101 FORMAT("OP =",I4,"  AV,ER R,I LINES =",4F9.4)
10 CONTINUE
end do
WRITE(23,81)
WRITE(23,202)
202 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=- ")
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
204 FORMAT(" AVERAGE  X,Y  LINES, J=1 ")
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
206 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=+ ")
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
207 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=- ")
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
220 FORMAT(" DIAGONAL TORELONS J=0, Pp=+, q=1 ")
WRITE(23,81)

do ID=1, NOPJ0PMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(" OPERATOR = ",I4)
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ0P(1,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1))<=EPS)goto 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT("  DT=",I3,"   AV,ER COR = ",2F8.4,"    E=",2F8.4)
end do
20 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,211)
211 FORMAT(" DIAGONAL TORELONS J=0, Pp=-, q=1 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 21
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
21 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,212)
212 FORMAT(" DIAGONAL TORELONS J=1, q=1 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 22
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
22 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,213)
213 FORMAT(" DIAGONAL TORELONS J=2, Pp=+, q=1 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 23
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
23 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,214)
214 FORMAT(" DIAGONAL TORELONS J=2, Pr=-, q=1 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 24
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
complex(real64) :: ACORLP,AVACLP
!***********************************************************
complex(real64) :: ACORLJ0PP,ACORLJ0PM,ACORLJ0MP,ACORLJ0MM
complex(real64) :: AVACLJ0PP,AVACLJ0PM,AVACLJ0MP,AVACLJ0MM
!***********************************************************
complex(real64) :: ACORLJ1P,ACORLJ1M
complex(real64) :: AVACLJ1P,AVACLJ1M
!***********************************************************
complex(real64) :: ACORLJ2PP,ACORLJ2PM,ACORLJ2MP,ACORLJ2MM
complex(real64) :: AVACLJ2PP,AVACLJ2PM,AVACLJ2MP,AVACLJ2MM
!***********************************************************
complex(real64) :: ACORLPMOMJ0P, ACORLPMOMJ0M
complex(real64) :: AVACLMOMJ0P, AVACLMOMJ0M
!***********************************************************
complex(real64) :: ACORLPMOMJ0, ACORLPMOMJ1, ACORLPMOMJ2
complex(real64) :: AVACLMOMJ0, AVACLMOMJ1, AVACLMOMJ2
!***********************************************************
complex(real64) :: ACORLPMOMJ2P, ACORLPMOMJ2M
complex(real64) :: AVACLMOMJ2P, AVACLMOMJ2M
!***********************************************************
DIMENSION CORP(NUMBIN,LMAX)
DIMENSION VAC(NUMBIN),LS(3),AV(NUMBIN)
!
OPEN (23, FILE="DIAGNALQ2IND.DAT")
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
80 FORMAT("****************************************************")
WRITE(23,181)
181 FORMAT("****************** OPERATORS WITH Q = 2 ************")
WRITE(23,80)
WRITE(23,81)
WRITE(23,80)
WRITE(23,80)
WRITE(23,81)
81 FORMAT(" *")
WRITE(23,91)
91 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=+")
WRITE(23,81)
82 FORMAT(" ******************* ")
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
101 FORMAT("OP =",I4,"  AV,ER R,I LINES =",4F9.4)
10 CONTINUE
end do
WRITE(23,81)
WRITE(23,202)
202 FORMAT(" AVERAGE  X,Y  LINES, J=0, Pp=- ")
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
204 FORMAT(" AVERAGE  X,Y  LINES, J=1 ")
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
206 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=+ ")
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
207 FORMAT(" AVERAGE  X,Y  LINES, J=2, Pp=- ")
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
220 FORMAT(" DIAGONAL TORELONS J=0, Pp=+, q=2 ")
WRITE(23,81)

do ID=1, NOPJ0PMOM
WRITE(23,81)
WRITE(23,82)
WRITE(23,94) ID
94 FORMAT(" OPERATOR = ",I4)
WRITE(23,82)
DO NT=1,LMAX
DO IB=1,NBIN
CORP(IB,NT)=real(ACORLPMOMJ0P(2,IB,NT,ID,ID))
end do
end do
IF(ABS(CORP(1,1))<=EPS)goto 20
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102 FORMAT("  DT=",I3,"   AV,ER COR = ",2F8.4,"    E=",2F8.4)
end do
20 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,211)
211 FORMAT(" DIAGONAL TORELONS J=0, Pp=-, q=2 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 21
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
21 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,212)
212 FORMAT(" DIAGONAL TORELONS J=1, q=2 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 22
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
22 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,213)
213 FORMAT(" DIAGONAL TORELONS J=2, Pp=+, q=2 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 23
CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,CORP,VAC)
DO NT=1,MAXDTLS
WRITE(23,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
end do
23 CONTINUE
end do
!*********************************************************************
WRITE(23,81)
WRITE(23,214)
214 FORMAT(" DIAGONAL TORELONS J=2, Pr=-, q=2 ")
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
IF(ABS(CORP(1,1))<=EPS)goto 24
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
use loop_builder, only : THERML1
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
complex(real32) :: U11,UB11,A11,B11,C11,UC11
!
CALL cpu_time(t1)
do I4=1,LX4
!
do IBL=1,IBLOK
!
IF(IBL==1)THEN
DO MU=1,3
DO NN=1,LSIZEB
DO IJ=1,NCOL2
UB11(IJ,NN,MU,1)=U11(IJ,NN,I4,MU)
end do
end do
end do
goto 11
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
901 FORMAT("[Info][Time]","     Thermal Line time:", F8.4)
!
CALL cpu_time(t1)
CALL POT
CALL cpu_time(t2)
WRITE(6,902) real(t2-t1, kind=real32)
902 FORMAT("[Info][Time]","     Correlation Creation time:", F8.4)
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
complex(real32) :: A11,B11,C11,UINT11,UREN11,UC11,DUM11,UCC11,UDD
complex(real32) :: UUC11
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
IF(IDIAG==1) CALL DIAG(MU)
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
IF(MU==NU)goto 40
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
IF(IDIAG/=1)goto 50
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
complex(real32) :: ADUM(NCOL,NCOL),U11,CSUM
!
complex(real32) :: AA(NCOL,NCOL)
complex(real32) :: DET,CNORM,CDET,DNCOL
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
complex(real32) :: C11(NCOL2),F11(NCOL2),DUM11(NCOL2)
complex(real32) :: A11(NCOL2),D11(NCOL2),E11(NCOL2)
complex(real32) :: AA11(NCOL2),DD11(NCOL2),EE11(NCOL2)
complex(real32) :: UDD,UC11
!
NU=MU+1
IF(MU==3)NU=1
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
complex(real32) :: UB1,UU1,CSUM,ADUM(NCOL,NCOL),UREN11
complex(real32) :: A11(NCOL2),B11(NCOL2),C11(NCOL2)
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
complex(real32) :: B11(NCOL2),C11(NCOL2)
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
complex(real32) :: B11(NCOL,NCOL),C11(NCOL,NCOL)
complex(real32) :: F11,F12,A11(NCOL,NCOL),S11(NCOL,NCOL)
complex(real32) :: T11(NCOL,NCOL)
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
IF((L1+LB(ID))>(LX1*1)) NU1=NU1-LX1
IF((L1+LB(ID))>(LX1*2)) NU1=NU1-LX1
IF((L1+LB(ID))>(LX1*3)) NU1=NU1-LX1
IF((L1+LB(ID))>(LX1*4)) NU1=NU1-LX1
IUPB(NN,1,ID)=NU1
ND1=NN-LB(ID)
IF((L1-LB(ID))<(1-LX1*0)) ND1=ND1+LX1
IF((L1-LB(ID))<(1-LX1*1)) ND1=ND1+LX1
IF((L1-LB(ID))<(1-LX1*2)) ND1=ND1+LX1
IF((L1-LB(ID))<(1-LX1*3)) ND1=ND1+LX1
IDNB(NN,1,ID)=ND1
!
NU2=NN+LB1(ID)
IF((L2+LB(ID))>(LX2*1)) NU2=NU2-LX12
IF((L2+LB(ID))>(LX2*2)) NU2=NU2-LX12
IF((L2+LB(ID))>(LX2*3)) NU2=NU2-LX12
IF((L2+LB(ID))>(LX2*4)) NU2=NU2-LX12
IUPB(NN,2,ID)=NU2
ND2=NN-LB1(ID)
IF((L2-LB(ID))<(1-LX2*0)) ND2=ND2+LX12
IF((L2-LB(ID))<(1-LX2*1)) ND2=ND2+LX12
IF((L2-LB(ID))<(1-LX2*2)) ND2=ND2+LX12
IF((L2-LB(ID))<(1-LX2*3)) ND2=ND2+LX12
IDNB(NN,2,ID)=ND2
!
NU3=NN+LB12(ID)
IF((L3+LB(ID))>(LX3*1)) NU3=NU3-LX123
IF((L3+LB(ID))>(LX3*2)) NU3=NU3-LX123
IF((L3+LB(ID))>(LX3*3)) NU3=NU3-LX123
IF((L3+LB(ID))>(LX3*4)) NU3=NU3-LX123
IUPB(NN,3,ID)=NU3
ND3=NN-LB12(ID)
IF((L3-LB(ID))<(1-LX3*0)) ND3=ND3+LX123
IF((L3-LB(ID))<(1-LX3*1)) ND3=ND3+LX123
IF((L3-LB(ID))<(1-LX3*2)) ND3=ND3+LX123
IF((L3-LB(ID))<(1-LX3*3)) ND3=ND3+LX123
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
complex(real32) :: U11,A11,B11,C11,D11,DUM11,ACT
!
DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN)
DIMENSION VAL(NUMBIN),AV(NUMBIN)
!
IBIN=NITER/NUMBIN
!
IF(IPR==0)THEN
JBIN=(ITER-1)/IBIN+1
!
IF(ITER==1)THEN
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
IF(NU/=4)THEN
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
IF(ITER==NITER)THEN
DO IP=1,6
DO JB=1,NUMBIN
ACTN(JB,IP)=ACTN(JB,IP)/(IBIN*LSIZE)
end do
end do
end if
!
end if
!
IF(IPR==1)THEN
!
WRITE(6,80)
80 FORMAT(" *****************************************************")
WRITE(6,100)
100 FORMAT("                 ACTION            ")
WRITE(6,80)
WRITE(6,81)
81 FORMAT(" *  ")
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
110 FORMAT("   AVER,ERR   SPACE ACTION =",2F12.8)
WRITE(6,112) AVRT,ERT
112 FORMAT("   AVER,ERR   TIME  ACTION =",2F12.8)
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
complex(real32) :: U11,A11,B11,C11,AKT1,TLINE,AVAC,CSUM
!
IBIN=NITER/NUMBIN
IF(IPR==1)goto 100
!
JBIN=(ITER-1)/IBIN+1
!
IF(ITER==1)THEN
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
goto 200
100 CONTINUE
!
WRITE(6,80)
80 FORMAT(" *****************************************************")
WRITE(6,108)
108 FORMAT("     THERMAL LINES           ")
WRITE(6,80)
WRITE(6,81)
81 FORMAT(" *  ")
do NB=1,NUMBIN
AV1(NB)=real(AVAC(NB))/IBIN
AV2(NB)=aimag(AVAC(NB))/IBIN
WRITE(6,109) NB,AV1(NB),AV2(NB)
109 FORMAT(" BIN =",I4,"   REAL,IMAG  POLY =",2F9.4)
20 CONTINUE
end do
CALL JACKM(NUMBIN,AV1,AVR,ERR)
CALL JACKM(NUMBIN,AV2,AVI,ERI)
WRITE(6,81)
WRITE(6,110) AVR,ERR
110 FORMAT("   AV,ER  REAL POLY =",2F12.8)
WRITE(6,111) AVI,ERI
111 FORMAT("   AV,ER  IMAG POLY =",2F12.8)
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
complex(real32) :: ADUM(NCOL,NCOL),U11,CSUM
!
complex(real32) :: aa(ncol,ncol)
complex(real32) :: det,cnorm,cdet,dncol
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
complex(real32) :: AA(NUMOP,NUMOP),B(NUMOP,NUMOP),det
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
complex(real32) :: A(NP,NP),DET,DD
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
complex(real32) :: A(NP,NP),SUM,CDUM
DIMENSION INDX(N),VV(NMAX)
!
D=1.
do I=1,N
AAMAX=0.
do J=1,N
IF(CABS(A(I,J))>AAMAX) AAMAX=CABS(A(I,J))
11 CONTINUE
end do
!     IF(AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX IN LUDCMP'
IF(AAMAX==0.) return
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
IF(DUM>=AAMAX)THEN
IMAX=I
AAMAX=DUM
end if
16 CONTINUE
end do
!
IF(J/=IMAX)THEN
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
IF(CABS(A(J,J))==0.) A(J,J)=TINY
IF(J/=N)THEN
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
IF((L1+1)>LX1) NU1=NU1-LX1
IUP(NN,1)=NU1
ND1=NN-1
IF((L1-1)<1) ND1=ND1+LX1
IDN(NN,1)=ND1
!
NU2=NN+LX1
IF((L2+1)>LX2) NU2=NU2-LX12
IUP(NN,2)=NU2
ND2=NN-LX1
IF((L2-1)<1) ND2=ND2+LX12
IDN(NN,2)=ND2
!
NU3=NN+LX12
IF((L3+1)>LX3) NU3=NU3-LX123
IUP(NN,3)=NU3
ND3=NN-LX12
IF((L3-1)<1) ND3=ND3+LX123
IDN(NN,3)=ND3
!
NU4=NN+LX123
IF((L4+1)>LX4) NU4=NU4-LX1234
IUP(NN,4)=NU4
ND4=NN-LX123
IF((L4-1)<1) ND4=ND4+LX1234
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
IF(DIFD(N)==0.0d0)THEN
ERR=99.0d0
goto 99
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
IF(SUMD/=0.0) AVER=SUMU/SUMD
IF(ESUM>0.0) ERR=SQRT(ESUM*NUM)
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
IF(ISUB==0)THEN
VACC=0.0d0
do JB=1,NBIN
IF(JB==IB)goto 4
VACC=VACC+VAC(JB)
4 CONTINUE
end do
end if
!
do NT=1,LMAX
CORR(NT)=0.0d0
do JB=1,NBIN
IF(JB==IB)goto 7
CORR(NT)=CORR(NT)+COR(JB,NT)
7 CONTINUE
end do
IF(ISUB==0) CORR(NT)=CORR(NT)-VACC*VACC/(NBIN-1)
6 CONTINUE
end do
ANORM=CORR(1)
IF(ANORM==0.0)goto 99
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
IF(SCOR(NT)>0.0) SCOR(NT)=SQRT(SCOR(NT))
AMM(NT)=AMM(NT)/NBIN
SMM(NT)=SMM(NT)/NBIN
SMM(NT)=(SMM(NT)-AMM(NT)*AMM(NT))*NBIN
IF(SMM(NT)>0.0) SMM(NT)=SQRT(SMM(NT))
10 CONTINUE
end do
!
IF(ISUB==0)THEN
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
IF(ISUB==0) CORR(NT)=CORR(NT)-VACC*VACC/NBIN
24 CONTINUE
end do
ANORM=CORR(1)
IF(ANORM==0.0)goto 99
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
IF(AW(I4)<=0.000001.OR.AW(I4+1)<=0.000001)THEN
NTMAX=I4-1
goto 99
end if
AMSM=0.0d0
FTM=(AW(I4)/AW(I4+1))
IF(FTM>1.0)THEN
AML=DLOG(FTM)
AMU=DLOG(2.0d0*FTM)
do NS=1,20
AMS=(AML+AMU)/2
FTS=(EXP(-AMS*(I4-1))+EXP(-(LT-I4+1)*AMS)) &
  & /(EXP(-AMS*(I4))+EXP(-(LT-I4)*AMS))
IF(FTS<FTM)THEN
AML=AMS
ELSE
AMU=AMS
end if
13 CONTINUE
end do
AMSM=(AML+AMU)/2
ELSE
NTMAX=I4-1
goto 99
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
complex(real32) :: A(NCOL,NCOL),B(NCOL,NCOL),C(NCOL,NCOL),CSUM
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
complex(real32) :: A(NCOL,NCOL),B(NCOL,NCOL),CC
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
complex(real32) :: A11(NCOL,NCOL),DUM11(NCOL,NCOL)
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
