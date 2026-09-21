module read_field_config
    use iso_fortran_env, only : int8, int32, int64, real32, real64
    use parameters
    implicit none

    ! Field configuration container
    complex(real32), allocatable :: U11(:,:,:,:,:,:,:)
end module