module axis_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use axis_module
  implicit none
contains

  subroutine get_axis0_fortran(nat,at_ptr,coord_ptr,rot_ptr,avmom_ptr,evec_ptr) bind(C,name="get_axis0_fortran")
    use,intrinsic :: iso_c_binding
    implicit none
    integer(c_int),value :: nat
    type(c_ptr),value :: at_ptr
    type(c_ptr),value :: coord_ptr

    type(c_ptr),value :: rot_ptr
    type(c_ptr),value :: evec_ptr
    type(c_ptr),value :: avmom_ptr

    ! Fortran pointer views of the incoming C buffers
    integer(c_int),pointer :: at(:)
    real(c_double),pointer :: coord(:)     ! length 3*natoms, flat
    real(c_double),pointer :: rot(:)     ! length 3
    real(c_double),pointer :: evec(:,:)    ! length 3, 3
    real(c_double),pointer :: avmom(:)    ! length 1
    
    real(wp) :: rot_f(3), evec_f(3,3), avmom_f

    call c_f_pointer(at_ptr,at, [nat])
    call c_f_pointer(coord_ptr,coord, [3*nat])
    call c_f_pointer(rot_ptr,rot, [3])
    call c_f_pointer(evec_ptr,evec, [3, 3])
    call c_f_pointer(avmom_ptr,avmom, [1])

    call axis_0(nat,at,coord,rot_f,avmom_f,evec_f)

    avmom = real(avmom_f,kind=c_double)
    rot = real(rot_f,kind=c_double)
    evec = real(evec_f,kind=c_double)

  end subroutine get_axis0_fortran

end module axis_exposed
