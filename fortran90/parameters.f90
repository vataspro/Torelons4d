module parameters
    use iso_fortran_env, only : int8, int32, int64, real32, real64
    implicit none

    !! Define pi
    real(real64), parameter :: PI = 4.0d0 * atan(1.0d0)

    !! Declare input parameters
    ! Configuration file path and file name (without configuration number)
    character(len=256) :: FILEPATH, FILENAME

    ! Lattice parameters
    integer :: NCOL, LX1, LX2, LX3, LX4

    ! File access parameters
    integer :: NCONFIG ! Number of configurations read, renamed from NITER
    integer :: CONFIG_STEP ! Number of steps between configurations, renamed from ICALLG
    integer :: CONFIG_START ! Index of first configuration to be read, renamed from ICMIN

    ! Blocking parameters
    integer :: MAX_BLOCKING_LEVEL ! Highest blocking level inclusive, renamed from IBLOK

    ! Smearing parameters
    real(real64) :: STAPLE_WEIGHT ! Weighting of staples in smearing procedure, renamed from PARBS
    real(real64) :: DIAGONAL_STAPLE_WEIGHT ! Weighting of diagonal staples in smearing procedure, renamed from PARBDS
    real(real64) :: TOL_SVD ! Tolerance on SVD-unitarisation algorithm. Set to between 1e-6 - 1e-10.

    !! Declare dependant parameters
    ! Lattice parameters
    integer :: SLICE_VOLUME ! Number of sites in one slice, renamed from LSIZEB
    integer :: LATTICE_VOLUME ! Total number of sites, renamed from LSIZE
    integer :: NCOL2 ! Square number of colours
    integer :: MAX_DELTA_T ! Maximum time extent of correlation functions, renamed from MAXDTLS

    ! File access parameters
    integer :: CONFIG_STOP ! Index of last configuration to be read, renamed from ICMAX

    !! TODO: implement all other global parameters

    contains

    ! Read all parameters from file and set dependant parameters
    subroutine initialise_parameters(parameter_filename)
        implicit none
        character(len=*), intent(in) :: parameter_filename

        ! Define variables to read in from parameter file
        namelist /params/ FILEPATH, FILENAME, NCOL, LX1, LX2, LX3, LX4, NCONFIG, CONFIG_STEP, CONFIG_START, &
        MAX_BLOCKING_LEVEL, STAPLE_WEIGHT, DIAGONAL_STAPLE_WEIGHT, TOL_SVD

        ! Read parameters from parameter file
        open(10, file=trim(parameter_filename))
        read(10, params)
        close(10)

        ! Set dependant parameters
        SLICE_VOLUME = LX1 * LX2 * LX3
        LATTICE_VOLUME = SLICE_VOLUME * LX4
        NCOL2 = NCOL * NCOL
        MAX_DELTA_T = LX4 / 2
        CONFIG_STOP = CONFIG_START + (NCONFIG - 1) * CONFIG_STEP
    end subroutine initialise_parameters
end module parameters
