module rmsd_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use crest_parameters,only:autoaa
  use irmsd_module
  use strucrd
  implicit none
contains

  subroutine get_quaternion_rmsd_fortran(natoms1,types1_ptr,coords1_ptr, &
                                         natoms2,types2_ptr,coords2_ptr, &
                                         rmsd_c,Umat_ptr) &
    bind(C,name="get_quaternion_rmsd_fortran")
    use,intrinsic :: iso_c_binding
    implicit none
    !> IN-/OUTPUTS
    integer(c_int),value :: natoms1,natoms2
    type(c_ptr),value :: types1_ptr,coords1_ptr
    type(c_ptr),value :: types2_ptr,coords2_ptr
    type(c_ptr),value :: Umat_ptr
    real(c_double),intent(out) :: rmsd_c
    !> LOCAL
    real(c_double),pointer :: Umat(:,:)        ! (3,3) each

    real(wp) :: rmsdval,rotmat(3,3)
    type(coord) :: ref
    type(coord) :: mol

    call ref%C_to_mol(natoms1,types1_ptr,coords1_ptr,.true.)
    call mol%C_to_mol(natoms2,types2_ptr,coords2_ptr,.true.)

    
    !> the quaternion rmsd, converted to angström
    rmsdval = rmsd(ref,mol,rotmat=rotmat)*autoaa

    call c_f_pointer(Umat_ptr,Umat, [3,3])
    Umat(:,:) = real(rotmat(:,:),c_double)
    rmsd_c = real(rmsdval,c_double)
  end subroutine get_quaternion_rmsd_fortran

end module rmsd_exposed
