module xyz_bridge
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use strucrd
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
  subroutine xyz_to_fortran(natoms,types_ptr,coords_ptr,mat_ptr) bind(C,name="xyz_to_fortran")
    use,intrinsic :: iso_c_binding
    implicit none
    integer(c_int),value :: natoms
    type(c_ptr),value :: types_ptr
    type(c_ptr),value :: coords_ptr
    type(c_ptr),value :: mat_ptr

    ! Fortran pointer views of the incoming C buffers
    integer(c_int),pointer :: types(:)
    real(c_double),pointer :: coords(:)     ! length 3*natoms, flat
    real(c_double),pointer :: mat(:,:)      ! shape (3,3), column-major

    integer :: i
    real(c_double) :: shift
    type(coord) :: mol

    ! Map raw C pointers to Fortran pointers with explicit shapes
    call c_f_pointer(types_ptr,types, [natoms])
    call c_f_pointer(coords_ptr,coords, [3*natoms])
    call c_f_pointer(mat_ptr,mat, [3,3])

    !>--- add to mol (for tests)
    call mol%C_to_mol(natoms,types,coords,.true.)
    write (*,*) 'Hello from Fortran. These are your coords:'
    do i = 1,mol%nat
      write (*,'(i5,3f20.10)') mol%at(i),mol%xyz(1:3,i)*0.52917726_wp
    end do
    write (*,*)

    ! --- Dummy computation on coordinates: translate everything by +1.0 Å
    do i = 1,mol%nat
      mol%xyz(:,i) = mol%xyz(:,i)+1.0_wp/0.52917726_wp
    end do
    call mol%mol_to_C(types,coords,.true.)

    ! --- Dummy 3x3 output matrix: natoms * Identity
    mat(:,:) = 0.0_c_double
    mat(1,1) = real(natoms,c_double)
    mat(2,2) = real(natoms,c_double)
    mat(3,3) = real(natoms,c_double)
  end subroutine xyz_to_fortran

!###############################################################################!
!===============================================================================!
!###############################################################################!

  subroutine xyz_to_fortran_pair(n1,types1_ptr,coords1_ptr, &
                                 n2,types2_ptr,coords2_ptr, &
                                 mat1_ptr,mat2_ptr) &
    bind(C,name="xyz_to_fortran_pair")
    use,intrinsic :: iso_c_binding
    implicit none

    ! sizes
    integer(c_int),value :: n1,n2

    ! incoming pointers
    type(c_ptr),value :: types1_ptr,coords1_ptr
    type(c_ptr),value :: types2_ptr,coords2_ptr
    type(c_ptr),value :: mat1_ptr,mat2_ptr

    ! Fortran views
    integer(c_int),pointer :: types1(:),types2(:)
    real(c_double),pointer :: coords1(:),coords2(:)      ! length 3*n1 / 3*n2
    real(c_double),pointer :: mat1(:,:),mat2(:,:)        ! (3,3) each

    integer :: i
    real(c_double) :: s1,s2

    ! Map raw pointers
    call c_f_pointer(types1_ptr,types1, [n1])
    call c_f_pointer(coords1_ptr,coords1, [3*n1])
    call c_f_pointer(types2_ptr,types2, [n2])
    call c_f_pointer(coords2_ptr,coords2, [3*n2])
    call c_f_pointer(mat1_ptr,mat1, [3,3])
    call c_f_pointer(mat2_ptr,mat2, [3,3])

    ! Dummy operations:
    ! - translate set #1 by +0.1 Å, set #2 by +0.2 Å
    s1 = 0.1_c_double
    s2 = 0.2_c_double

    do i = 1,3*n1
      coords1(i) = coords1(i)+s1
    end do
    do i = 1,3*n2
      coords2(i) = coords2(i)+s2
    end do

    ! Matrices: n1*I3 and n2*I3
    mat1(:,:) = 0.0_c_double
    mat2(:,:) = 0.0_c_double
    mat1(1,1) = real(n1,c_double); mat1(2,2) = real(n1,c_double); mat1(3,3) = real(n1,c_double)
    mat2(1,1) = real(n2,c_double); mat2(2,2) = real(n2,c_double); mat2(3,3) = real(n2,c_double)
  end subroutine xyz_to_fortran_pair

end module xyz_bridge

