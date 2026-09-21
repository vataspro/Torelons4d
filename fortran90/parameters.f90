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

    contains

    ! Read all parameters from file and set dependant parameters
    subroutine initialise_parameters(parameter_filename)
        implicit none
        character(len=*) :: parameter_filename

        ! Define variables to read in from parameter file
        namelist /params/ NCOL, LX1, LX2, LX3, LX4

        ! Read parameters from parameter file
        open(10, file=trim(parameter_filename))
        read(10, params)
        close(10)

        ! Set dependant parameters
        SLICE_VOLUME = LX1 * LX2 * LX3
        LATTICE_VOLUME = SLICE_VOLUME * LX4
        NCOL2 = NCOL * NCOL
        MAX_DELTA_T = LX4 / 2
    end subroutine
end module