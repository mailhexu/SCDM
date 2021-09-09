#include "abi_common.h"
module m_test_utils
  use defs_basis
  implicit none
  private

  public:: assertion
  public:: all_close_real
  public:: all_close_integer
contains
  subroutine assertion(condition, message)
    logical, intent(in) :: condition
    character (*), intent(in) :: message
    if (.not. condition) then
       ABI_ERROR(message)
    end if
  end subroutine assertion

  function all_close_real(a, b, tol) result (ret)
    real(dp) :: a(:), b(:)
    real(dp) :: tol
    logical :: ret
    ret=sum(abs(a-b))< tol
  end function all_close_real

  function all_close_integer(a, b) result (ret)
    integer :: a(:), b(:)
    logical :: ret
    ret=(sum(abs(a-b))==0)
  end function all_close_integer


end module m_test_utils
