!*******************************************************************
!* This module is a drop-in replacement for CREST's strucrd module
!* which primarily exports the "coord" type. We do not want the
!* entire source of the original file.
!******************************************************************

module strucrd
  use iso_fortran_env,only:wp => real64
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

!  contains
!    procedure :: deallocate => deallocate_coord !> clear memory space
!    procedure :: open => opencoord              !> read an coord file
!    procedure :: write => writecoord            !> write
!    procedure :: dist => coord_getdistance      !> calculate distance between two atoms
!    procedure :: angle => coord_getangle        !> calculate angle between three atoms
!    procedure :: dihedral => coord_getdihedral  !> calculate dihedral angle between four atoms
!    procedure :: get_CN => coord_get_CN         !> calculate coordination number
!    procedure :: get_z => coord_get_z           !> calculate nuclear charge
!    procedure :: cn_to_bond => coord_cn_to_bond !> generate neighbour matrix from CN
  end type coord
!=========================================================================================!

!==============================================================================!
contains   !> MODULE PROCEDURES START HERE
!==============================================================================!

end module strucrd
