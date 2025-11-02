module xyz_bridge
  use, intrinsic :: iso_c_binding
  implicit none
contains

  !> C-entry subroutine. Interoperable signature via c_ptr + c_f_pointer.
  !!
  !! Arguments (C-side / Python):
  !!   natoms      : integer(c_int), by value
  !!   types_ptr   : pointer to int32 array of length natoms
  !!   coords_ptr  : pointer to float64 array of length 3*natoms (x1,y1,z1,x2,y2,z2,...)
  !!   mat_ptr     : pointer to float64 array of length 9 (3x3), column-major expected
  !!
  !! Effects:
  !!   - overwrites coords in-place (dummy example: translates all atoms by +0.1 in each coord)
  !!   - writes a 3x3 matrix (dummy example: identity * real(natoms))
  !!
  subroutine xyz_to_fortran(natoms, types_ptr, coords_ptr, mat_ptr) bind(C, name="xyz_to_fortran")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: natoms
    type(c_ptr),   value :: types_ptr
    type(c_ptr),   value :: coords_ptr
    type(c_ptr),   value :: mat_ptr

    ! Fortran pointer views of the incoming C buffers
    integer(c_int),  pointer :: types(:)
    real(c_double),  pointer :: coords(:)     ! length 3*natoms, flat
    real(c_double),  pointer :: mat(:,:)      ! shape (3,3), column-major

    integer :: i
    real(c_double) :: shift

    ! Map raw C pointers to Fortran pointers with explicit shapes
    call c_f_pointer(types_ptr,  types,  [natoms])
    call c_f_pointer(coords_ptr, coords, [3*natoms])
    call c_f_pointer(mat_ptr,    mat,    [3,3])

    ! --- Dummy computation on coordinates: translate everything by +0.1 Å
    shift = 0.1_c_double
    do i = 1, 3*natoms
      coords(i) = coords(i) + shift
    end do

    ! --- Dummy 3x3 output matrix: natoms * Identity
    mat(:,:) = 0.0_c_double
    mat(1,1) = real(natoms, c_double)
    mat(2,2) = real(natoms, c_double)
    mat(3,3) = real(natoms, c_double)
  end subroutine xyz_to_fortran

end module xyz_bridge

