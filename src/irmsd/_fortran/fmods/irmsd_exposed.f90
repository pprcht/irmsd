module irmsd_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use strucrd
  use crest_parameters
  use axis_module
  use irmsd_module
  use canonical_mod
  implicit none
contains

  subroutine get_irmsd_fortran(natoms1,types1_ptr,coords1_ptr, &
                               natoms2,types2_ptr,coords2_ptr, &
                               iinversion_c,rmsd_c,types_out1_ptr,coords_out1_ptr, &
                               types_out2_ptr,coords_out2_ptr) &
    bind(C,name="get_irmsd_fortran")
    use,intrinsic :: iso_c_binding
    implicit none
    !> IN-/OUTPUTS
    integer(c_int),value :: natoms1,natoms2,iinversion_c
    type(c_ptr),value :: types1_ptr,coords1_ptr
    type(c_ptr),value :: types2_ptr,coords2_ptr
    type(c_ptr),value :: types_out1_ptr
    type(c_ptr),value :: coords_out1_ptr
    type(c_ptr),value :: types_out2_ptr
    type(c_ptr),value :: coords_out2_ptr
    real(c_double),intent(out) :: rmsd_c

    integer :: iinversion
    integer(c_int),pointer :: types_out1_c(:)
    real(c_double),pointer :: coords_out1_c(:)
    integer(c_int),pointer :: types_out2_c(:)
    real(c_double),pointer :: coords_out2_c(:)
    type(coord) :: mol,ref
    real(wp) :: rmsdval,tmpd(3),tmpdist
    integer :: i
    type(rmsd_cache) :: rcache
    type(canonical_sorter) :: canmol
    type(canonical_sorter) :: canref
    logical :: mirror

    logical,parameter :: debug = .false.

    call ref%C_to_mol(natoms1,types1_ptr,coords1_ptr,.true.)
    call mol%C_to_mol(natoms2,types2_ptr,coords2_ptr,.true.)

    if (natoms1 /= natoms2) then
      error stop 'both molecules need to have the same number of atoms'
    end if

    iinversion = iinversion_c

    !> move ref to CMA and align rotational axes
    call axis(ref%nat,ref%at,ref%xyz)

    !> allocate memory
    call rcache%allocate(ref%nat)

    !> canonical atom ranks
    call canref%init(ref,invtype='apsp+',heavy=.false.)
    rcache%stereocheck = .not. (canref%hasstereo(ref))
    call canref%shrink()
    if (debug) write (stdout,*) 'false enantiomers possible?: ',rcache%stereocheck
    select case (iinversion)
    case (0)  !> whatever rcache%stereocheck says
      mirror = .true.
    case (1)  !> force on
      mirror = .true.
      rcache%stereocheck = .true.
    case (2) !> force off
      mirror = .false.
      rcache%stereocheck = .false.
    end select

    if (debug) write (stdout,*) 'allow inversion?:            ',mirror

    call canmol%init(mol,invtype='apsp+',heavy=.false.)
    call canmol%shrink()

    !> check if we can work with the determined ranks
    if (checkranks(ref%nat,canref%rank,canmol%rank)) then
      if (debug) write (stdout,*) 'using canonical atom identities as rank backend'
      rcache%rank(:,1) = canref%rank(:)
      rcache%rank(:,2) = canmol%rank(:)
      if (debug) then
        write (stdout,*) 'iRMSD ranks:'
        write (stdout,*) 'atom',' rank(ref)',' rank(mol)'
        do i = 1,ref%nat
          write (stdout,*) i,rcache%rank(i,1),rcache%rank(i,2)
        end do
        write (stdout,*)
      end if
    else
      !> if not, fall back to atom types
      if (debug) write (stdout,*) 'using atom types as rank backend'
      call fallbackranks(ref,mol,ref%nat,rcache%rank)
    end if

    call min_rmsd(ref,mol,rcache=rcache,rmsdout=rmsdval,align=.true.)

    if (debug) then
      do i = 1,mol%nat
        tmpd(:) = (mol%xyz(:,i)-ref%xyz(:,i))**2
        tmpdist = sqrt(sum(tmpd(:)))*autoaa
        if (tmpdist > 0.01_wp) then
          write (stdout,*) i,mol%at(i),tmpdist
        end if
      end do
    end if

    rmsdval = rmsdval*autoaa
    rmsd_c = real(rmsdval,c_double)

    call c_f_pointer(types_out1_ptr,types_out1_c, [natoms1])
    call c_f_pointer(coords_out1_ptr,coords_out1_c, [3*natoms1])

    call c_f_pointer(types_out2_ptr,types_out2_c, [natoms2])
    call c_f_pointer(coords_out2_ptr,coords_out2_c, [3*natoms2])

    call ref%mol_to_C(types_out1_c,coords_out1_c,.true.)
    call mol%mol_to_C(types_out2_c,coords_out2_c,.true.)

  end subroutine get_irmsd_fortran
end module irmsd_exposed
