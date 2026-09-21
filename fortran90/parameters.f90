module parameters
    use iso_fortran_env, only : int8, int32, int64, real32, real64
    implicit none

    ! Declare lattice variables
    integer :: NCOL, LX1, LX2, LX3, LX4

    ! Declare dependant variables
    integer :: SLICE_VOLUME ! Number of sites in one slice, renamed from LSIZEB
    integer :: LATTICE_VOLUME ! Total number of sites, renamed from LSIZE
    integer :: NCOL2 ! Square number of colours
    integer :: MAX_DELTA_T ! Maximum time extent of correlation functions, renamed from MAXDTLS

    !! TODO: implement all other global parameters
end module