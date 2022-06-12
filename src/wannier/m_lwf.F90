
!===============================================================
! m_lwf: constructing of lattice wannier function
! 
!> @description: Select Column Density matrix method for
!   generating Wannier functions.
!===============================================================

#include "abi_common.h"

module m_lwf
  use m_wannier_builder, only: WannierBuilder_witheigen_t
  use defs_basis, only: dp
  implicit none

  public :: lwf
  private

  type :: ifc

  end type ifc

  !===============================================================
  ! scdmk type:
  !> @ description: the class for scdmk method.
  !===============================================================
  type :: lwf
     complex(dp), allocatable :: ifc(:,:,:)
     real(dp), allocatable :: kpts(:,:)
     integer :: nkpts
     type(WannierBuilder_witheigen_t) :: myscdmk
   contains
     procedure :: initialize
     procedure :: finalize
     !procedure :: read_ifc
     procedure :: build_kmesh
     procedure :: get_eigen
     procedure :: build_lwf
  end type lwf
  
contains
  subroutine initialize(self)
    class(lwf),intent(inout) :: self
  end subroutine initialize

  subroutine finalize(self)
    class(lwf), intent(inout) :: self
  end subroutine finalize


  subroutine write_lwf(self)
    class(lwf) , intent(inout):: self
  end subroutine write_lwf

  subroutine build_kmesh(self)
    class(lwf) , intent(inout):: self
  end subroutine build_kmesh

  subroutine get_eigen(self)
      class(lwf) , intent(inout):: self
  end subroutine get_eigen

  subroutine build_lwf(self)
    class(lwf),intent(inout) :: self
    
  end subroutine build_lwf

end module m_lwf
