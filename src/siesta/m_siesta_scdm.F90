
#include "../common/abi_common.h"

module m_siesta_scdm
    use defs_basis, only: dp, fnlen
    use m_kpoints, only: kpoints, build_Rgrid
    use m_test_utils, only: assertion
    use m_scdm, only: scdmk_t
    implicit none
    private

    type :: siesta_scdm
        type(scdmk_t) :: objscdmk
    contains
        procedure :: initialize
        procedure :: finalize
    end type siesta_scdm

contains

    subroutine initialize(self)
        implicit none
        class(siesta_scdm), intent(in) :: self
    end subroutine initialize

    subroutine finalize(self)
        implicit none
        class(siesta_scdm), intent(in) :: self
    end subroutine finalize

    subroutine test_siesta_scdm(self)
        implicit none
        type(siesta_scdm) :: self

        integer :: kmesh(3), nspin
        integer, allocatable :: Rlist(:,:)
        real(dp) :: ne, beta
        integer :: ncell, i
        character(fnlen) :: fname
        integer, parameter :: nbasis = 14, nkpt = 7*7*7
        real(dp), allocatable :: evals(:, :)  ! nbasis, nkpt
        complex(dp), allocatable :: evecs(:, :, :) ! nbasis, nband, nkpt
        integer :: exclude_bands(0)



        ABI_SFREE(Rlist)
        ABI_FREE(evals)
        ABI_FREE(evecs)
    end subroutine test_siesta_scdm

end module m_siesta_scdm
