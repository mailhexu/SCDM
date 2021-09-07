
!===============================================================
! m_lwf: constructing of lattice wannier function
! 
!> @description: Select Column Density matrix method for
!   generating Wannier functions.
!===============================================================
module m_lwf
  use m_scdm, only: scdmk
  implicit none
  integer, parameter :: dp = 8
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
     type(scdmk) :: myscdmk
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: read_ifc
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
    class(lwf) :: self
  end subroutine write_lwf
end module m_lwf
