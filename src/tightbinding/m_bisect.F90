#include "abi_common.h"

module m_bisect
  use defs_basis, only: dp
  implicit none
  private

  type, public :: bisect_solver
     integer :: counter
     integer :: max_counter=60
     real(dp) :: a, b
     real(dp) :: fa=0.0_dp, fb=0.0_dp
   contains
     procedure :: initialize
     procedure :: nextx
     procedure :: too_many_steps
  end type bisect_solver

contains

  subroutine initialize(self, a, b, fa, fb, max_counter)
    class(bisect_solver), intent(inout) :: self
    real(dp), intent(in) :: a, b, fa, fb
    integer, optional, intent(in) :: max_counter
    self%a=a
    self%b=b
    self%fa=fa
    self%fb=fb
    if(fa*fb>0) then
       print *, "Error: The upper and lower bound of bisect solver has the same sign. Exit"
       stop
    end if
    self%counter=0
    if (present(max_counter)) self%max_counter=max_counter
  end subroutine initialize

  subroutine nextx(self, x, fx)
    class(bisect_solver), intent(inout) :: self
    real(dp),  intent(inout) :: x
    real(dp), intent(in) :: fx
    self%counter=self%counter+1
    if ( (fx<0) .eqv. (self%fa<0) ) then
       self%a = x
       self%fa = fx
    else
       self%b = x
       self%fb = fx
    end if
    x = ( self%a + self%b ) / 2.0D+00
  end subroutine nextx

  function too_many_steps(self)
    class(bisect_solver), intent(inout) :: self
    logical :: too_many_steps
    too_many_steps= (self%counter>self%max_counter)
  end function too_many_steps

end module m_bisect

