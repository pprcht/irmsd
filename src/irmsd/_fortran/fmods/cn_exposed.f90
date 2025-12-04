module cn_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use strucrd
  implicit none
contains

  subroutine get_cn_fortran(natoms,types_ptr,coords_ptr,cn_ptr) bind(C,name="get_cn_fortran")
    use,intrinsic :: iso_c_binding
    implicit none
    integer(c_int),value :: natoms
    type(c_ptr),value :: types_ptr
    type(c_ptr),value :: coords_ptr
    type(c_ptr),value :: cn_ptr

    ! Fortran pointer views of the incoming C buffers
    real(c_double),pointer :: cn(:)         ! length natoms
    
    integer :: i
    type(coord) :: mol
    real(wp),allocatable :: cn_f(:)

    ! Map raw C pointers to Fortran pointers with explicit shapes
    call c_f_pointer(cn_ptr,cn, [natoms])

    !>--- add to mol object
    call mol%C_to_mol(natoms,types_ptr,coords_ptr,.true.)

    !>--- get CN
    call mol%get_CN(cn_f,cn_type="cov")

    !>-- typecast to C
    do i=1,natoms
       cn(i) = cn_f(i)
    enddo

  end subroutine get_cn_fortran

end module cn_exposed

