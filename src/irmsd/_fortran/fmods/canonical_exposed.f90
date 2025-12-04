module canonical_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use canonical_mod
  use strucrd
  implicit none
contains

  subroutine get_canonical_sorter_fortran(natoms,types_ptr,coord_ptr, &
  &    wbo_ptr,invtype_ptr,heavy,rank_ptr) bind(C,name="get_canonical_sorter_fortran")
    implicit none
    integer(c_int),value :: natoms
    type(c_ptr),value :: types_ptr
    type(c_ptr),value :: coord_ptr
    type(c_ptr),value :: wbo_ptr
    type(c_ptr),value :: rank_ptr
    character(kind=c_char),dimension(*) :: invtype_ptr
    logical(c_bool),value :: heavy
    logical :: heavy_f

    ! Fortran pointer views of the incoming C buffers
    real(c_double),pointer :: wbo(:,:)      ! length natoms x natoms
    integer(c_int),pointer :: rank(:)       ! length natoms

    type(coord) :: mol
    type(canonical_sorter) :: canonical

    integer :: n
    character(:),allocatable :: invtype_f

    heavy_f = heavy

    ! TODO: refactor string handling into a utility function?
    ! NOTE: Not sure whether this is the best way to handle it.
    ! Find the length of the C string (up to null character)
    n = 0
    do while (invtype_ptr(n+1) /= c_null_char)
      n = n+1
    end do

    ! Allocate a normal Fortran CHARACTER string
    allocate (character(len=n) :: invtype_f)

    ! Copy the characters
    invtype_f = transfer(invtype_ptr(1:n),invtype_f)

    call mol%C_to_mol(natoms,types_ptr,coord_ptr,.true.) ! last arguments indicates convert to bohr
    if (c_associated(wbo_ptr)) then
      call c_f_pointer(wbo_ptr,wbo, [natoms,natoms])
      call canonical%init(mol,wbo,invtype_f,heavy_f)
    else
      call canonical%init(mol,invtype=invtype_f,heavy=heavy_f)
    end if

    call c_f_pointer(rank_ptr,rank,shape(canonical%rank))

    rank = canonical%rank

  end subroutine get_canonical_sorter_fortran

  subroutine get_ids_from_connect_fortran(natoms,types_ptr,connect_ptr,heavy,rank_ptr) bind(C,name="get_ids_from_connect_fortran")
    implicit none
    integer(c_int),value :: natoms
    type(c_ptr),value :: types_ptr
    type(c_ptr),value :: connect_ptr
    type(c_ptr),value :: rank_ptr
    logical(c_bool),value :: heavy
    logical :: heavy_f

    ! Fortran pointer views of the incoming C buffers
    integer(c_int),pointer :: types_c(:)        ! length natoms
    integer(c_int),pointer :: connect_c(:,:)    ! length natoms x natoms
    integer(c_int),pointer :: rank(:)           ! length natoms

    type(canonical_sorter) :: canonical

    real(wp),allocatable :: tmpwbo(:,:)

    heavy_f = heavy
    call c_f_pointer(types_ptr,types_c,[natoms]) 
    call c_f_pointer(connect_ptr,connect_c, [natoms,natoms]) 

    !>-- blowup
    allocate(tmpwbo(natoms,natoms))
    tmpwbo(:,:) = real(connect_c(:,:),wp)
    call canonical%init_connect(types_c,tmpwbo,invtype="asap+",heavy=heavy_f)

    call c_f_pointer(rank_ptr,rank,shape(canonical%rank))

    rank = canonical%rank

  end subroutine get_ids_from_connect_fortran

end module canonical_exposed
