module linalg_mod
  use iso_c_binding
  implicit none
contains
  subroutine saxpy(n, a, x, y) bind(C, name="saxpy")
    integer(c_int), value :: n
    real(c_double), value :: a
    real(c_double) :: x(*), y(*)
    integer :: i
    do i = 1, n
      y(i) = y(i) + a * x(i)
    end do
  end subroutine saxpy

  subroutine scal(n, a, x) bind(C, name="scal")
    integer(c_int), value :: n
    real(c_double), value :: a
    real(c_double) :: x(*)
    integer :: i
    do i = 1, n
      x(i) = a * x(i)
    end do
  end subroutine scal

  subroutine dot(n, x, y, out) bind(C, name="dot")
    integer(c_int), value :: n
    real(c_double) :: x(*), y(*), out(*)
    integer :: i
    real(c_double) :: s
    s = 0.0_c_double
    do i = 1, n
      s = s + x(i)*y(i)
    end do
    out(1) = s
  end subroutine dot
end module linalg_mod


