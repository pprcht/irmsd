module symmetry_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use symmetry_i,only:schoenflies,schoenflies_elements,getsym_params
  implicit none
  private

  public :: get_symmetry_fortran
  public :: get_symmetry_elements_fortran

  real(wp),parameter :: bohr = 0.52917726_wp

contains

  subroutine get_symmetry_fortran(nat,at_ptr,coord_ptr,thr,thr_primary, &
  &          maxorder,maxcycles,sym_ptr,symlen) bind(C,name="get_symmetry_fortran")
    !***********************************
    !* Schoenflies symbol, null-terminated into caller buffer sym_ptr(symlen).
    !* coord_ptr(3*nat) in Angstrom; thr, thr_primary (atom pairing) in Bohr.
    !***********************************
    implicit none
    integer(c_int),value :: nat
    type(c_ptr),value :: at_ptr
    type(c_ptr),value :: coord_ptr
    real(c_double),value :: thr,thr_primary
    integer(c_int),value :: maxorder,maxcycles
    type(c_ptr),value :: sym_ptr
    integer(c_int),value :: symlen

    integer(c_int),pointer :: at(:)
    real(c_double),pointer :: coord(:)
    real(wp),allocatable :: xyz(:,:)
    character(len=8) :: symbol

    call c_f_pointer(at_ptr,at, [nat])
    call c_f_pointer(coord_ptr,coord, [3*nat])
    xyz = reshape(real(coord,wp)/bohr, [3,nat])

    symbol = ''
    call schoenflies(int(nat),int(at),xyz,symbol, &
    &    params(thr,thr_primary,maxorder,maxcycles))
    call symbol_to_c(symbol,sym_ptr,symlen)
  end subroutine get_symmetry_fortran

  subroutine get_symmetry_elements_fortran(nat,at_ptr,coord_ptr,thr,thr_primary, &
  &          maxorder,maxcycles,sym_ptr,symlen, &
  &          maxel,nel,type_ptr,order_ptr,rmat_ptr,tvec_ptr,axis_ptr,point_ptr, &
  &          maxdev_ptr,perm_ptr) bind(C,name="get_symmetry_elements_fortran")
    !***********************************
    !* As get_symmetry_fortran, plus the first min(nel,maxel) of nel elements
    !* in caller buffers, C order: type/order/maxdev (maxel), perm (maxel,nat)
    !* 0-based, rmat (maxel,3,3) as x' = R x + t, tvec/axis/point (maxel,3).
    !* tvec, point, maxdev in Angstrom.
    !***********************************
    implicit none
    integer(c_int),value :: nat
    type(c_ptr),value :: at_ptr
    type(c_ptr),value :: coord_ptr
    real(c_double),value :: thr,thr_primary
    integer(c_int),value :: maxorder,maxcycles
    type(c_ptr),value :: sym_ptr
    integer(c_int),value :: symlen
    integer(c_int),value :: maxel
    integer(c_int),intent(out) :: nel
    type(c_ptr),value :: type_ptr,order_ptr,rmat_ptr,tvec_ptr
    type(c_ptr),value :: axis_ptr,point_ptr,maxdev_ptr,perm_ptr

    integer(c_int),pointer :: at(:)
    real(c_double),pointer :: coord(:)
    integer(c_int),pointer :: eltype_c(:),elorder_c(:),perm_c(:,:)
    real(c_double),pointer :: rmat_c(:,:,:),tvec_c(:,:),axis_c(:,:)
    real(c_double),pointer :: point_c(:,:),maxdev_c(:)

    real(wp),allocatable :: xyz(:,:)
    character(len=8) :: symbol
    integer :: n,k
    integer,allocatable  :: eltype(:),elorder(:),perm(:,:)
    real(wp),allocatable :: rmat(:,:,:),tvec(:,:),axis(:,:),point(:,:),maxdev(:)

    call c_f_pointer(at_ptr,at, [nat])
    call c_f_pointer(coord_ptr,coord, [3*nat])
    xyz = reshape(real(coord,wp)/bohr, [3,nat])

    symbol = ''
    call schoenflies_elements(int(nat),int(at),xyz,symbol,n,eltype,elorder, &
    &                         rmat,tvec,axis,point,maxdev,perm, &
    &                         params(thr,thr_primary,maxorder,maxcycles))
    call symbol_to_c(symbol,sym_ptr,symlen)
    nel = int(n,c_int)

    if (maxel < 1) return
    call c_f_pointer(type_ptr,eltype_c, [maxel])
    call c_f_pointer(order_ptr,elorder_c, [maxel])
    call c_f_pointer(rmat_ptr,rmat_c, [3,3,maxel])
    call c_f_pointer(tvec_ptr,tvec_c, [3,maxel])
    call c_f_pointer(axis_ptr,axis_c, [3,maxel])
    call c_f_pointer(point_ptr,point_c, [3,maxel])
    call c_f_pointer(maxdev_ptr,maxdev_c, [maxel])
    call c_f_pointer(perm_ptr,perm_c, [nat,maxel])

    do k = 1,min(n,int(maxel))
      eltype_c(k) = eltype(k)
      elorder_c(k) = elorder(k)
!>-- transpose so that the caller reads R[k,i,j] = R_ij
      rmat_c(:,:,k) = transpose(rmat(:,:,k))
      tvec_c(:,k) = tvec(:,k)*bohr
      axis_c(:,k) = axis(:,k)
      point_c(:,k) = point(:,k)*bohr
      maxdev_c(k) = maxdev(k)*bohr
      perm_c(:,k) = perm(:,k)-1
    end do
  end subroutine get_symmetry_elements_fortran

  function params(thr,thr_primary,maxorder,maxcycles) result(paramar)
    !***********************************
    !* getsym_params with the user-adjustable entries overridden.
    !***********************************
    real(c_double),intent(in) :: thr,thr_primary
    integer(c_int),intent(in) :: maxorder,maxcycles
    real(wp) :: paramar(11)
    paramar = getsym_params(real(thr,wp))
    paramar(2) = real(maxorder,wp)
    paramar(3) = real(maxcycles,wp)
    paramar(5) = real(thr_primary,wp)
  end function params

  subroutine symbol_to_c(symbol,sym_ptr,symlen)
    !***********************************
    !* Copy a Fortran string into a null-terminated C buffer, truncated.
    !***********************************
    character(len=*),intent(in) :: symbol
    type(c_ptr),value :: sym_ptr
    integer(c_int),value :: symlen
    character(kind=c_char),pointer :: sym(:)
    integer :: i,n

    call c_f_pointer(sym_ptr,sym, [symlen])
    n = min(len_trim(symbol),int(symlen)-1)
    do i = 1,n
      sym(i) = symbol(i:i)
    end do
    sym(n+1) = c_null_char
  end subroutine symbol_to_c

end module symmetry_exposed
