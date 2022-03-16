
#include "../common/abi_common.h"

module m_tb_scdm
    use defs_basis, only: dp, fnlen
    use m_tight_binding, only: TBHam
    use m_kpoints, only: kpoints, build_Rgrid
    use m_test_utils, only: assertion
    use m_scdm, only: scdmk_t
    implicit none
    private

    type :: tb_scdm
        type(TBHam) :: ham
        type(kpoints) :: kpts
        type(scdmk_t) :: objscdmk
    contains
        procedure :: finalize
    end type tb_scdm

    public :: tb_scdm
    public :: run_scdm
    public :: test_tb_scdm
contains

    subroutine finalize(self)
        implicit none
        class(tb_scdm), intent(in) :: self

    end subroutine finalize

    subroutine run_scdm(ham, kpts, scdmk)
        implicit none
        class(TBHam), intent(inout) :: ham
        type(kpoints), intent(inout) :: kpts
        type(scdmk_t), intent(inout) :: scdmk

    end subroutine run_scdm

    subroutine test_tb_scdm()
        implicit none
        type(TBHam) :: ham
        type(kpoints) :: kpts
        type(scdmk_t) :: scdmk

        integer :: kmesh(3), nspin
        integer, allocatable :: Rlist(:,:)
        real(dp) :: ne, beta
        integer :: ncell, i
        character(fnlen) :: fname
        integer, parameter :: nbasis = 14, nkpt = 7*7*7
        real(dp), allocatable :: evals(:, :)  ! nbasis, nkpt
        complex(dp), allocatable :: evecs(:, :, :) ! nbasis, nband, nkpt
        integer :: exclude_bands(0)

        print *, "Reading TB Hamiltonian"
        fname = "/home/hexu/projects/computeU/test/SrTiO3/STO.nc"
        call ham%read_from_netcdf(fname)
        kmesh = [7, 7, 7]
        call kpts%monkhorst_pack(kmesh)
        call build_Rgrid(kmesh, Rlist)

        ABI_MALLOC(evals, (nbasis, nkpt))
        ABI_MALLOC(evecs, (nbasis, nbasis, nkpt))
        call ham%solve_all(kpts%kpts, evals, evecs)
        print *, "Initialize scdmk for tb."
        call scdmk%initialize(evals=evals, psi=evecs, &
            & kpts=kpts%kpts, kweights=kpts%kweights, &
            & Rlist= Rlist, &
        &  nwann=5,  &
        & disentangle_func_type=3, mu=6.0, sigma=3.5, &
        & exclude_bands=exclude_bands )

        print *, "E{gamma}", evals(:,1)
        print *, "Setting anchor point"
        !call scdmk%set_anchor(anchor_kpt=[0.0_dp, 0.0_dp, 0.0_dp], anchor_ibands=[1,2,3, 4, 5])
        call scdmk%set_anchor(anchor_kpt=[0.0_dp, 0.0_dp, 0.0_dp]) 
        call scdmk%run_all()

        ABI_SFREE(Rlist)
        ABI_FREE(evals)
        ABI_FREE(evecs)
    end subroutine test_tb_scdm

end module m_tb_scdm
