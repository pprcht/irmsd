!*******************************************************************
!* This module is a drop-in replacement for CREST's strucrd module
!* which primarily exports the "coord" type. We do not want the
!* entire source of the original file.
!******************************************************************

module strucrd
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use,intrinsic :: iso_c_binding
  use crest_cn_module,only:calculate_CN
  use miscdata,only:PSE
  implicit none

!>--- private module variables and parameters
  private

!>--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp

  !>--- selected public exports of the module
  public :: coord, C_to_ensemble
  public :: i2e

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
    procedure :: get_CN => coord_get_CN         !> calculate coordination number
    procedure :: get_z => coord_get_z           !> calculate nuclear charge
    procedure :: cn_to_bond => coord_cn_to_bond !> generate neighbour matrix from CN
    procedure :: swap => atswp                  !> swap two atoms coordinates and their at() entries
    procedure :: C_to_mol
    procedure :: mol_to_C
    procedure :: append
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

  subroutine C_to_mol(self,natoms_c,at_ptr,xyz_ptr,convert_to_Bohr)
    !***************************************************
    !* Pass number of atoms and coordinats as ptrs from C
    !* and allocate coord object in fortran types
    !***************************************************
    implicit none
    class(coord) :: self
    integer(c_int),value :: natoms_c

    type(c_ptr),value :: at_ptr
    type(c_ptr),value :: xyz_ptr
    logical,intent(in) :: convert_to_Bohr

    integer(c_int),pointer :: at_c(:)
    real(c_double),pointer :: xyz_c(:)
    integer :: i,j,k
    real(wp) :: convert

    call c_f_pointer(at_ptr,at_c, [natoms_c])
    call c_f_pointer(xyz_ptr,xyz_c, [3*natoms_c])

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
        self%xyz(j,i) = real(xyz_c(k),wp)*convert
      end do
    end do
  end subroutine C_to_mol

  subroutine C_to_ensemble(structures,nstructures_c,many_natoms_ptr, &
  &                        many_at_ptr,many_xyz_ptr,convert_to_Bohr_c)
    !***************************************************
    !* Retrieve a long list of coordinates and atom types
    !* and return an array of coord objectas
    !***************************************************
    implicit none
    !> IN-/OUTPUTS
    type(coord),allocatable,intent(out) :: structures(:)
    integer(c_int),value :: nstructures_c
    type(c_ptr),value :: many_natoms_ptr
    type(c_ptr),value :: many_at_ptr
    type(c_ptr),value :: many_xyz_ptr
    logical(c_bool),value,intent(in) :: convert_to_Bohr_c
    !> LOCAL
    integer(c_int),pointer :: many_natoms_c(:)
    integer(c_int),pointer :: many_at_c(:)
    real(c_double),pointer :: many_xyz_c(:)
    logical :: convert_to_Bohr
    integer :: i,j,k,fulllength,nat,k1,k2
    real(wp) :: convert

    convert_to_Bohr = convert_to_Bohr_c
    call c_f_pointer(many_natoms_ptr,many_natoms_c, [nstructures_c])
    fulllength = 0
    do i = 1,nstructures_c
      fulllength = fulllength+many_natoms_c(i)
    end do
    call c_f_pointer(many_at_ptr,many_at_c, [fulllength])
    call c_f_pointer(many_xyz_ptr,many_xyz_c, [3*fulllength])

    allocate (structures(nstructures_c))
    if (convert_to_Bohr) then
      convert = 1.0_wp/bohr  !> Input in Ang, convert to Bohr
    else
      convert = 1.0_wp  !> Input already in Bohr
    end if
    k1 = 0
    k2 = 0
    do i = 1,nstructures_c
      nat = many_natoms_c(i)
      allocate (structures(i)%xyz(3,nat),source=0.0_wp)
      allocate (structures(i)%at(nat),source=0)
      structures(i)%nat = nat
      do j = 1,nat
         k1=k1+1
         structures(i)%at(j) = many_at_c(k1)
         do k=1,3
           k2=k2+1
           structures(i)%xyz(k,j) = many_xyz_c(k2)
         enddo
      end do
    end do
  end subroutine C_to_ensemble

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
        xyz_c(k) = self%xyz(j,i)*convert
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

!==================================================================!

  subroutine coord_get_CN(self,cn,cn_type,cn_thr,dcndr)
    implicit none
    class(coord) :: self
    real(wp),intent(out),allocatable :: cn(:)
    real(wp),intent(in),optional :: cn_thr
    character(len=*),intent(in),optional :: cn_type
    real(wp),intent(out),optional :: dcndr(3,self%nat,self%nat)
    if (self%nat <= 0) return
    if (.not.allocated(self%xyz).or..not.allocated(self%at)) return
    allocate (cn(self%nat),source=0.0_wp)
    call calculate_CN(self%nat,self%at,self%xyz,cn, &
    & cntype=cn_type,cnthr=cn_thr,dcndr=dcndr)
  end subroutine coord_get_CN

!==================================================================!

  subroutine coord_get_z(self,z)
    implicit none
    class(coord) :: self
    real(wp),intent(out),allocatable :: z(:)
    integer :: i
    if (self%nat <= 0) return
    if (.not.allocated(self%xyz).or..not.allocated(self%at)) return
    allocate (z(self%nat),source=0.0_wp)
    do i = 1,self%nat
      z(i) = real(self%at(i),wp)-real(ncore(self%at(i)))
      if (self%at(i) > 57.and.self%at(i) < 72) z(i) = 3.0_wp
    end do
  end subroutine coord_get_z

!==================================================================!

  subroutine coord_cn_to_bond(self,cn,bond,cn_type,cn_thr)
    implicit none
    class(coord) :: self
    real(wp),intent(out),allocatable :: cn(:)
    real(wp),intent(out),allocatable,optional :: bond(:,:)
    real(wp),intent(in),optional :: cn_thr
    character(len=*),intent(in),optional :: cn_type
    if (self%nat <= 0) return
    if (.not.allocated(self%xyz).or..not.allocated(self%at)) return
    allocate (cn(self%nat),source=0.0_wp)
    call calculate_CN(self%nat,self%at,self%xyz,cn, &
    & cntype=cn_type,cnthr=cn_thr,bond=bond)
  end subroutine coord_cn_to_bond

!=============================================================!

  pure elemental integer function ncore(at)
    integer,intent(in) :: at
    if (at .le. 2) then
      ncore = 0
    elseif (at .le. 10) then
      ncore = 2
    elseif (at .le. 18) then
      ncore = 10
    elseif (at .le. 29) then   !zn
      ncore = 18
    elseif (at .le. 36) then
      ncore = 28
    elseif (at .le. 47) then
      ncore = 36
    elseif (at .le. 54) then
      ncore = 46
    elseif (at .le. 71) then
      ncore = 54
    elseif (at .le. 79) then
      ncore = 68
    elseif (at .le. 86) then
      ncore = 78
    end if
  end function ncore

!==============================================================================!

  subroutine append(self,iunit,filename)
    implicit none
    class(coord) :: self
    integer,intent(in) :: iunit
    character(len=*),intent(in),optional :: filename
    integer :: i,iunit_loc
    if (present(filename)) then
      open (newunit=iunit_loc,file=filename)
    else
      iunit_loc = iunit
    end if
    write (iunit_loc,'(1x,i0)') self%nat
    !write (iunit_loc,'(1x,a,F20.10)') 'Energy=',self%energy
    write (iunit_loc,*)
    do i = 1,self%nat
      write (iunit_loc,'(a,1x,3F25.10)') PSE(self%at(i)),self%xyz(1:3,i)*bohr
    end do
    if (present(filename)) close (iunit_loc)
  end subroutine append

!============================================================================!

  function i2e(iat,spec) result(esym)
    integer,intent(in) :: iat
    character(len=*),intent(in),optional :: spec  !> this is ignored, legacy arg
    character(len=:),allocatable :: esym
    if(present(spec)) continue !> again, some legacy stuff
    esym = PSE(iat)
  end function i2e

end module strucrd
