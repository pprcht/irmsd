module sorter_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use strucrd
  use crest_parameters
  use axis_module
  use irmsd_module
  use canonical_mod
  implicit none
contains

  subroutine sorter_exposed_xyz_fortran( &
    &                     nat,nall,xyzall_ptr,atall_ptr, &
    &                     groups_ptr,rthresh,iinversion,allcanon_c,printlvl &
    &                   ) bind(C,name="sorter_exposed_xyz_fortran")
    use,intrinsic :: iso_c_binding
    implicit none

    ! C-facing arguments
    integer(c_int),value :: nat
    integer(c_int),value :: nall
    type(c_ptr),value :: xyzall_ptr
    type(c_ptr),value :: atall_ptr
    type(c_ptr),value :: groups_ptr
    real(c_double),value :: rthresh
    integer(c_int),value :: iinversion
    logical(c_bool),value :: allcanon_c
    integer(c_int),value :: printlvl

    ! Fortran pointer views of C buffers
    real(c_double),pointer :: xyzall(:)
    integer(c_int),pointer :: atall(:)
    integer(c_int),pointer :: groups(:)

    ! Local Fortran variables
    type(coord),allocatable :: structures(:)
    integer :: i,j,k1,k2,l
    logical :: allcanon

    ! Map raw C pointers to Fortran pointers with shapes
    call c_f_pointer(xyzall_ptr,xyzall, [3*nat*nall])
    call c_f_pointer(atall_ptr,atall, [nat*nall])
    call c_f_pointer(groups_ptr,groups, [nall])
    allcanon = allcanon_c

    write (*,*) "Hello from Fortran!"
    write (*,*) "We have this info:"
    write (*,'(a,i0)') "xyzall with length ",size(xyzall,1)
    write (*,'(a,i0)') "atall with length ",size(atall,1)
    write (*,'(a,i0)') "groups with length ",size(groups,1)
    write(*,*) "rthr",rthresh
    write(*,*) "iinversion",iinversion
    write(*,*) "allcanon",allcanon
    write(*,*) "printlvl",printlvl

    ! Call original Fortran routine
    allocate (structures(nall))
    k1 = 0
    k2 = 0
    do i = 1,nall
      structures(i)%nat = nat 
      allocate (structures(i)%at(nat),source=0)
      allocate (structures(i)%xyz(3,nat),source=0.0_wp)
      do j = 1,nat
        k1 = k1+1
        structures(i)%at(j) = atall(k1)
        do l = 1,3
          k2 = k2+1
          structures(i)%xyz(l,j) = xyzall(k2)*aatoau  !> Angström to BOHR
        end do
      end do
      call structures(i)%append(6)
    end do

    stop "something is broken in the following call, stopping here"
    call cregen_irmsd_sort(nall,structures,groups,rthresh,iinversion,allcanon,printlvl)

    write(*,*)
    write(*,*) "groups:"
    write(*,*) groups(1:nall)
    stop
    !> coordinates for each structure have been aligned and sorted,
    !> so we need to pass them back

    k1 = 0
    k2 = 0
    do i = 1,nall
      do j = 1,nat
        k1 = k1+1
        atall(k1) = structures(i)%at(j)
        do l = 1,3
          k2 = k2+1
          xyzall(k2) = structures(i)%xyz(l,j)*autoaa
        end do
      end do
    end do

    !> free memory
    deallocate (structures)

    !> group assignment (== unique conformers) on "groups" array

  end subroutine sorter_exposed_xyz_fortran

  subroutine cregen_irmsd_sort(nall,structures,groups,rthresh,iinversion,allcanon,printlvl)
!*******************************************************
!* Proof-of-concept routine to analyze an
!* ensemble only via the iRMSD procedure.
!* Conformers are identified by the rthr threshold only
!*******************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nall
    type(coord),intent(inout),target :: structures(nall)
    integer,intent(inout) :: groups(nall)
    real(wp),intent(in) :: RTHRESH
    integer,intent(in) :: iinversion
    logical,intent(in),optional :: allcanon
    integer,intent(in),optional :: printlvl

    !> LOCAL
    integer :: i,j,ii,jj,T,Tn,nallpairs,cc,nat
    integer :: gcount
    integer :: prlvl,iunit
    type(rmsd_cache),allocatable :: rcaches(:)
    type(coord),allocatable,target :: workmols(:)
    type(canonical_sorter),allocatable :: sorters(:)
    real(wp),allocatable :: rmsds(:)
    type(coord),pointer :: ref,mol
    type(coord) :: molloc
    real(wp) :: rmsdval,runtime,RTHR
    logical :: stereocheck,individual_IDs

    logical,parameter :: debug = .true.

!>--- handle optional arguments
    if (present(allcanon)) then
      individual_IDs = allcanon
    else
      individual_IDs = .false.
    end if
    if (present(printlvl)) then
      prlvl = printlvl
    else
      prlvl = 1
    end if

!>--- set up parallelization
!     ...
    T = 1 !> doing it serial for now

!>--- set up parameters (note we are working with BOHR internally)
    RTHR = RTHRESH*aatoau

!>--- print some sorting data
    if (prlvl > 0) then
      write (stdout,'(a)') 'Info for iRMSD sorting:'
      write (stdout,'(2x,a,i9)') 'number of structures  :',nall
      write (stdout,'(2x,a,f9.5,a)') 'RTHR (RMSD threshold) :',RTHR*autoaa,' Å'
      write (stdout,'(2x,a,i9)') 'OpenMP threads        :',T
      write (stdout,'(2x,a,l9)') 'Individual atom IDs?  :',individual_IDs
      write (stdout,'(2x,a)',advance='no') 'False rotamer check?  :'
      select case (iinversion)
      case (0)
        write (stdout,'(a9)') 'auto'
      case (1)
        write (stdout,'(a9)') 'on'
      case (2)
        write (stdout,'(a9)') 'off'
      end select
      write (stdout,*)
    end if

!>--- Set up atom identities (either for all, or just the first structure)
    if (individual_IDs) then
      allocate (sorters(nall))
    else
      allocate (sorters(1))
    end if
    ref => structures(1)
    ! !$omp parallel &
    ! !$omp shared(sorters, structures, stereocheck) &
    ! !$omp private(mol,ii)
    ! !$omp do schedule(dynamic)
    do ii = 1,nall
      mol => structures(ii)
      call axis(mol%nat,mol%at,mol%xyz)
      if (individual_IDs.or.ii == 1) then
        call sorters(ii)%init(mol,invtype='apsp+',heavy=.false.)
      end if
      if (ii == 1) then
        stereocheck = .not. (sorters(ii)%hasstereo(ref))
      end if
      if (individual_IDs.or.ii == 1) then
        call sorters(ii)%shrink()
      end if
    end do
    ! !$omp end do
    ! !$omp end parallel

    !>--- allow user to set inversion check (false rotamers)
    select case (iinversion)
    case (0)
      continue
    case (1)
      stereocheck = .true.
    case (2)
      stereocheck = .false.
    end select
    if (prlvl > 1) then
      write (stdout,'(a,l2)') 'CREGEN> Check for false rotamers (geometry inversion)? -->',stereocheck
    end if

!>--- allocate work cache
    if (prlvl > 0) then
      write (stdout,'(a)',advance='no') 'CREGEN> Allocating iRMSD work cache ... '
      flush (stdout)
    end if
    allocate (rcaches(T))
    ref => structures(1)
    nat = ref%nat
    allocate (workmols(T))
    do i = 1,T
      mol => workmols(i)
      allocate (mol%at(ref%nat))
      allocate (mol%xyz(3,ref%nat))
      nullify (mol)
      call rcaches(i)%allocate(ref%nat)
      rcaches(i)%stereocheck = stereocheck
    end do
    if (prlvl > 0) then
      write (stdout,'(a)') 'done.'
      write (stdout,*)
    end if

!>--- run the checks
    gcount = maxval(groups(:))
    do ii = 1,nall
!>--- find next unassigned conformer and assign a new group
      if (groups(ii) .ne. 0) cycle
      gcount = gcount+1
      groups(ii) = gcount

!>--- Then, cross-check all other unassigned conformers
      cc = 1  !> again, serial implementation for now
      ! !$omp parallel &
      ! !$omp shared(nall, nat, groups, individual_IDs, sorters, rcaches) &
      ! !$omp shared(workmols, structures, ii) &
      ! !$omp private(jj,rmsdval,cc)
      ! !$omp do schedule(dynamic)
      do jj = ii+1,nall
        !cc = omp_get_thread_num()+1
        if (groups(jj) .ne. 0) cycle
        if (individual_IDs) then
          rcaches(cc)%rank(1:nat,1) = sorters(ii)%rank(1:nat)
          rcaches(cc)%rank(1:nat,2) = sorters(jj)%rank(1:nat)
        else
          rcaches(cc)%rank(1:nat,1) = sorters(1)%rank(1:nat)
          rcaches(cc)%rank(1:nat,2) = sorters(1)%rank(1:nat)
        end if
        workmols(cc)%nat = structures(jj)%nat
        workmols(cc)%at(:) = structures(jj)%at(:)
        workmols(cc)%xyz(:,:) = structures(jj)%xyz(:,:)
        call min_rmsd(structures(ii),workmols(cc), &
        &        rcache=rcaches(cc),rmsdout=rmsdval)
        if (rmsdval < RTHR) groups(jj) = gcount
      end do
      ! !$omp end do
      ! !$omp end parallel
    end do

    if (debug) then
      write (*,*) 'assigned groups, and count'
      do ii = 1,maxval(groups(:))
        write (*,*) ii,count(groups(:) == ii)
      end do
    end if

  end subroutine cregen_irmsd_sort

end module sorter_exposed
