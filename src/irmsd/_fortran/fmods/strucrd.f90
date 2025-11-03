!*******************************************************************
!* This module is a drop-in replacement for CREST's strucrd module
!* which primarily exports the "coord" type. We do not want the
!* entire source of the original file.
!******************************************************************

module strucrd
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use,intrinsic :: iso_c_binding
  implicit none

!>--- private module variables and parameters
  private

!>--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp

  !>--- selected public exports of the module
  public :: coord

!=========================================================================================!
  !coord class. contains a single structure
  !by convention coordinates are in atomic units (Bohr) for a single structure!
  type :: coord
    !********************************************!
    !> data that's typically used in coord type <!
    !********************************************!
    !>-- number of atoms
    integer :: nat = 0
    !>-- atom types as integer, dimension will be at(nat)
    integer,allocatable  :: at(:)
    !>-- atomic coordinates, by convention in Bohrs
    real(wp),allocatable :: xyz(:,:)

    !**************************************!
    !> (optional) data, often not present <!
    !**************************************!
    !>-- energy
    real(wp) :: energy = 0.0_wp
    !>-- molecular charge
    integer :: chrg = 0
    !>-- multiplicity information
    integer :: uhf = 0
    !>-- number of bonds
    integer :: nbd = 0
    !>-- bond info
    integer,allocatable :: bond(:,:)
    !>-- lattice vectors
    real(wp),allocatable :: lat(:,:)

    !>-- atomic charges
    real(wp),allocatable :: qat(:)

  contains
    procedure :: deallocate => deallocate_coord !> clear memory space
!    procedure :: open => opencoord              !> read an coord file
!    procedure :: write => writecoord            !> write
!    procedure :: dist => coord_getdistance      !> calculate distance between two atoms
!    procedure :: angle => coord_getangle        !> calculate angle between three atoms
!    procedure :: dihedral => coord_getdihedral  !> calculate dihedral angle between four atoms
!    procedure :: get_CN => coord_get_CN         !> calculate coordination number
!    procedure :: get_z => coord_get_z           !> calculate nuclear charge
!    procedure :: cn_to_bond => coord_cn_to_bond !> generate neighbour matrix from CN
    procedure :: swap => atswp                  !> swap two atoms coordinates and their at() entries
    procedure :: C_to_mol
    procedure :: mol_to_C
  end type coord
!=========================================================================================!

!==============================================================================!
contains   !> MODULE PROCEDURES START HERE
!==============================================================================!

  subroutine deallocate_coord(self)
    !**************************************
    !* deallocate data of an coord object
    !**************************************
    implicit none
    class(coord) :: self
    self%nat = 0
    self%energy = 0.0_wp
    self%chrg = 0
    self%uhf = 0
    self%nbd = 0
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%bond)) deallocate (self%bond)
    if (allocated(self%lat)) deallocate (self%lat)
    if (allocated(self%qat)) deallocate (self%qat)
    return
  end subroutine deallocate_coord

  subroutine C_to_mol(self,natoms_c,at_c,xyz_c,convert_to_Bohr)
    !***************************************************
    !* Pass number of atoms and coordinats from C types
    !* and allocate coord object in fortran types
    !***************************************************
    implicit none
    class(coord) :: self
    integer(c_int),value :: natoms_c
    integer(c_int),pointer :: at_c(:)
    real(c_double),pointer :: xyz_c(:)
    logical,intent(in) :: convert_to_Bohr
    integer :: i,j,k
    real(wp) :: convert
    call self%deallocate()
    if (convert_to_Bohr) then
      convert = 1.0_wp/bohr  !> Input in Ang, convert to Bohr
    else
      convert = 1.0_wp  !> Input already in Bohr
    end if
    self%nat = int(natoms_c)
    allocate (self%at(self%nat),source=0)
    allocate (self%xyz(3,self%nat),source=0.0_wp)
    k = 0
    do i = 1,self%nat
      self%at(i) = int(at_c(i))
      do j = 1,3
        k = k+1
        self%xyz(j,i) = real(xyz_c(k),wp) * convert
      end do
    end do
  end subroutine C_to_mol

  subroutine mol_to_C(self,at_c,xyz_c,convert_to_Ang)
    implicit none
    class(coord)           :: self
    integer(c_int),pointer :: at_c(:)
    real(c_double),pointer :: xyz_c(:)
    logical,intent(in) :: convert_to_Ang
    integer :: i,j,k
    real(wp) :: convert
    if (convert_to_Ang) then
      convert = bohr  !> Output in Ang, convert from Bohr
    else
      convert = 1.0_wp  !> Output in Bohr
    end if
    !> Copy Z numbers
    do i = 1,self%nat
      at_c(i) = int(self%at(i),c_int)
    end do
    !> Pack coordinates as a flat array (x,y,z for atom 1, then atom 2, ...)
    k = 0
    do i = 1,self%nat
      do j = 1,3
        k = k+1
        xyz_c(k) = self%xyz(j,i) * convert
      end do
    end do
  end subroutine mol_to_C

  subroutine atswp(self,ati,atj)
    !********************************
    !* swap atom ati with atj in mol
    !********************************
    implicit none
    class(coord),intent(inout) :: self
    integer,intent(in) :: ati,atj
    real(wp) :: xyztmp(3)
    integer :: attmp
    xyztmp(1:3) = self%xyz(1:3,ati)
    attmp = self%at(ati)
    self%xyz(1:3,ati) = self%xyz(1:3,atj)
    self%at(ati) = self%at(atj)
    self%xyz(1:3,atj) = xyztmp(1:3)
    self%at(atj) = attmp
  end subroutine atswp

end module strucrd
