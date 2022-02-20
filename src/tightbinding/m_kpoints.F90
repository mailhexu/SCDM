#include "abi_common.h"
module m_kpoints
  use defs_basis
  use m_test_utils
  implicit none
  private

  !===============================================================
  ! a object containing the kpoints and the weights.
  !> @
  !===============================================================
  type, public ::  kpoints
     integer :: nk=0
     real(dp), allocatable :: kpts(:,:)
     real(dp), allocatable :: kweights(:)
   contains
     procedure, public :: monkhorst_pack => kpoints_monkhorst_pack
     procedure, public :: finalize
  end type kpoints

  public :: Monkhorst_Pack
  public :: test_kpoints
  public :: build_Rgrid
contains

  ! -4/2-> -1, 4/2->2. Then (-1,0,1,2)
  pure function div2(x) result(y)
      integer, intent(in) :: x
      integer :: y
      if(mod(x, 2)==0 .and. x<0) then 
          y=x/2+1
      else
          y=x/2
      endif 
  end function div2

  subroutine build_Rgrid(kmesh, Rlist)
    integer, intent(in) :: kmesh(3)
    integer, allocatable, intent(out) :: Rlist(:,:)
    integer :: n, i1, i2, i3, i

    n=kmesh(1) * kmesh(2) *kmesh(3)
    ABI_MALLOC(Rlist, (3,n))
    i=0

    ! Note that C/Fortran integer division is "truncate towards 0" division,
    ! whereas Python one is "floor" division.
    ! For C/Fortran, the behavior for even and odd numbers is 
    !  not consistent and need special treatment in div2.
    do i3 = div2(-kmesh(3)), div2(kmesh(3))
      do i2 = div2(-kmesh(2)), div2(kmesh(2))
        do i1 = div2(-kmesh(1)), div2(kmesh(1))
             i=i+1
             Rlist(:, i) = [i1, i2, i3]
          end do
       end do
    end do
  end subroutine build_Rgrid


  !===============================================================
  ! Generate Monkhorst pack k-grid
  !> @ kmesh: an integer array of size 3
  !> @ kpts: output. size:(3,nkpts).
  !===============================================================
  subroutine Monkhorst_Pack(kmesh, kpts)
    integer, intent(in) :: kmesh(3)
    !real(dp) :: shift(3)
    real(dp), allocatable, intent(inout) :: kpts(:,:)
    integer :: i, nk, ix, iy, iz
    nk=product(kmesh)
    ABI_MALLOC(kpts,(3, nk))
    i=0
    do ix=0, kmesh(1)-1
       do iy=0, kmesh(2)-1
          do iz=0, kmesh(3)-1
             i=i+1
             kpts(:, i)=[ix, iy, iz]
          end do
       end do
    end do
    do i =1, nk
       kpts(:, i)=(kpts(:,i)+0.5)/kmesh -0.5
    end do
  end subroutine Monkhorst_Pack


  !===============================================================
  ! Generate monkhorst pack grid into a kpoints type.
  !> @ kmesh: a triplet defining the k mesh.
  !===============================================================
  subroutine kpoints_monkhorst_pack(self, kmesh)
    class(kpoints), intent(inout) :: self
    integer :: kmesh(3)
    call Monkhorst_Pack(kmesh, self%kpts)
    self%nk=size(self%kpts, 2)
    ABI_MALLOC(self%kweights, (self%nk))
    self%kweights(:)=1.0/self%nk
  end subroutine kpoints_monkhorst_pack

  !===============================================================
  ! Finalize kpoints object
  !===============================================================
  subroutine finalize(self)
    class(kpoints), intent(inout) :: self
    ABI_SFREE(self%kpts)
    ABI_SFREE(self%kweights)
    self%nk=0
  end subroutine finalize

  !===============================================================
  ! Unit test to the monkhorst function
  !===============================================================

  subroutine test_kpoints()
    type(kpoints) :: ks
    integer :: kmesh(3)=[3,3,3], i

    integer, allocatable :: Rlist(:,:)
    call build_Rgrid(kmesh, Rlist)
    do i = 1, 27 
       print *, Rlist(:, i)
    end do

    call ks%monkhorst_pack(kmesh)
    call assertion(ks%nk==27, " Wrong number of k points.")
    call assertion(all_close_real(ks%kpts(:,1), [-1.0d0/3, -1.0d0/3, -1.0d0/3], 1d-6), "wrong k-points generated")
    call assertion(all_close_real(ks%kpts(:,27), [1.0d0/3, 1.0d0/3, 1.0d0/3], 1d-6), "wrong k-points generated")
    call assertion(abs(ks%kweights(1)-3.7037037312984467E-002)< 1d-6, "wrong k-weights")
    call ks%finalize()
  end subroutine test_kpoints

end module m_kpoints

