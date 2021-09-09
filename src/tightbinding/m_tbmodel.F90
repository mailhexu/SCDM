#include "abi_common.h"

module m_tb_model
  use defs_basis , only : dp,dpc, j_dpc, two_pi, fnlen
  use m_mathfuncs, only : eigensolver
  use m_test_utils, only: assertion, all_close_real
  use m_tight_binding,only: TBHam
  use m_kpoints, only : kpoints
  use m_occ, only: find_efermi, get_density, get_density_matrix
  implicit none

  private

  !===============================================================
  ! Tight binding Hamiltonian
  !===============================================================
  type, public :: TBModel
     type(TBHam), pointer :: Htb
     type(kpoints) :: kpts
     integer :: nspin
     real(dp) :: ne
     real(dp) :: efermi
     real(dp) :: beta
     real(dp):: tol_fermi=1.0d-8
     real(dp), allocatable :: occs(:,:)
     real(dp), allocatable :: evals(:,:)
     complex(dp), allocatable :: evecs(:,:,:)
     real(dp), allocatable :: den(:)
     complex(dp), allocatable :: den_mat(:,:) ! density matrix
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_kmesh
     procedure :: set_tb_ham
     procedure :: run
  end type TBModel

  public :: test_tbmodel
  public :: test_tbmodel_supercell
  public :: test_tbmodel_netcdf
  public :: test_tbmodel_w90_hr

contains

  !===============================================================
  !
  !> @ tbham_fname: the file 
  !===============================================================
  subroutine initialize(self, tb_ham, kmesh, ne,nspin, beta, tol_fermi)
    class(TBModel), intent(inout) :: self
    type(TBHam), target :: tb_ham
    integer, intent(in) :: kmesh(:), nspin
    real(dp), intent(in) :: ne, beta
    real(dp), optional :: tol_fermi
    call self%set_tb_ham(tb_ham)
    call self%kpts%monkhorst_pack(kmesh)
    self%nspin=nspin
    self%ne=ne
    self%beta=beta
    if (present(tol_fermi)) self%tol_fermi=tol_fermi
    ABI_MALLOC(self%evals, (self%Htb%nwann, self%kpts%nk ))
    ABI_MALLOC(self%evecs, (self%Htb%nwann, self%Htb%nwann, self%kpts%nk ))
    ABI_MALLOC(self%den, (self%Htb%nwann))
    ABI_MALLOC(self%den_mat, (self%Htb%nwann, self%Htb%nwann))
    ABI_MALLOC(self%occs, (self%Htb%nwann, self%kpts%nk))
  end subroutine initialize

  subroutine finalize(self)
    class(TBModel), intent(inout) :: self
    call self%Htb%finalize()
    nullify(self%Htb)
    call self%kpts%finalize
    ABI_SFREE(self%occs)
    ABI_SFREE(self%evals)
    ABI_SFREE(self%evecs)
    ABI_SFREE(self%den)
    ABI_SFREE(self%den_mat)
  end subroutine finalize

  subroutine set_tb_ham(self, tb_ham)
    class(TBModel), intent(inout) :: self
    type(TBHam), target :: tb_ham
    self%Htb=>tb_ham
  end subroutine set_tb_ham

  subroutine set_kmesh(self, kmesh)
    class(TBModel), intent(inout) :: self
    integer, intent(in) ::  kmesh(:)
    call self%kpts%monkhorst_pack(kmesh)
  end subroutine set_kmesh

  !===============================================================
  ! 
  !> @
  !===============================================================
  subroutine run(self)
    class(TBModel), intent(inout) :: self
    call self%Htb%solve_all(kpts=self%kpts%kpts, evals=self%evals, evecs=self%evecs)
    call find_efermi(evals=self%evals, kweights=self%kpts%kweights, beta=self%beta, occs=self%occs, &
         & nspin=self%nspin, nel=self%ne, thr=self%tol_fermi, mu=self%efermi)
    call get_density(self%evecs, self%kpts%kweights,self%occs, self%den)
  end subroutine run


!############################ Tests

  subroutine test_tbmodel()
    use m_tight_binding, only: build_1d_chain

    type(TBHam) :: tb_ham
    type(TBModel) :: tb_model
    integer :: kmesh(3), nspin
    real(dp) :: ne, beta

    call build_1d_chain(tb_ham, on_site=cmplx(0.0_dp,0.0_dp, kind=dp), &
         & hopping=cmplx(1.0_dp,0.0_dp, kind=dp))
    kmesh=[18,1,1]
    !ne=1.0d-6
    ne=1.4
    nspin=1
    beta=0.01
    call tb_model%initialize( tb_ham, kmesh, ne,nspin, beta)

    call tb_model%run()
    call  assertion(all_close_real(tb_model%den, [ne], 1d-4 ), "wrong density " )
  end subroutine test_tbmodel


  subroutine test_tbmodel_supercell()
    use m_tight_binding, only: build_1d_chain

    type(TBHam) :: tb_ham, sc_tb_ham
    type(TBModel) :: tb_model
    integer :: kmesh(3), nspin
    real(dp) :: ne, beta
    integer :: ncell, i
    ncell=3

    call build_1d_chain(tb_ham, on_site=cmplx(0.0_dp,0.0_dp, kind=dp), &
         & hopping=cmplx(1.0_dp,0.0_dp, kind=dp))
    call tb_ham%fill_supercell_with_matrix([ncell,0,0,0,1,0,0,0,1], sc_tb_ham)
    kmesh=[31,1,1]
    !ne=1.0d-6
    ne=0.1*ncell
    nspin=1
    beta=0.01
    call tb_model%initialize( sc_tb_ham, kmesh, ne,nspin, beta)
    call tb_model%run()
    call  assertion(all_close_real(tb_model%den, [(ne/ncell, i=1, ncell)], 1d-4 ), "wrong density " )
  end subroutine test_tbmodel_supercell


  subroutine test_tbmodel_netcdf()
    use m_tight_binding, only: build_1d_chain

    character(fnlen) :: fname
    type(TBHam) :: tb_ham, sc_tb_ham
    type(TBModel) :: tb_model
    integer :: kmesh(3), nspin
    real(dp) :: ne, beta
    integer :: ncell, i
    ncell=1
    fname="/home/hexu/projects/computeU/test/SrTiO3/STO.nc"

    call tb_ham%read_from_netcdf(fname)
    !call tb_ham%fill_supercell_with_matrix([ncell,0,0,0,1,0,0,0,1], sc_tb_ham)
    kmesh=[8,8,8]
    !ne=1.0d-6
    ne=18*ncell
    nspin=1
    beta=0.01
    call tb_model%initialize( tb_ham, kmesh, ne,nspin, beta)
    call tb_model%run()
    print*, sum(tb_model%den)
    print*, tb_model%den
    !call  assertion(all_close_real(tb_model%den, [(ne/ncell, i=1, ncell)], 1d-4 ), "wrong density " )
  end subroutine test_tbmodel_netcdf

  subroutine test_tbmodel_w90_hr()
    use m_tight_binding, only: build_1d_chain

    character(fnlen) :: fname
    type(TBHam) :: tb_ham, sc_tb_ham
    type(TBModel) :: tb_model
    integer :: kmesh(3), nspin
    real(dp) :: ne, beta
    integer :: ncell, i
    ncell=1
    !fname="/home/hexu/projects/computeU/examples/LaNiO3/wannier90_hr.dat"
    !call tb_ham%read_from_wannier90_hr(fname)
    !call tb_ham%fill_supercell_with_matrix([ncell,0,0,0,1,0,0,0,1], sc_tb_ham)
    kmesh=[8,8,8]
    !ne=1.0d-6
    ne=25*ncell
    nspin=1
    beta=0.1
    call tb_model%initialize( tb_ham, kmesh, ne,nspin, beta)
    call tb_model%run()
    print *, "Fermi energy:", tb_model%efermi
    print*, sum(tb_model%den)
    print*, tb_model%den
    !call  assertion(all_close_real(tb_model%den, [(ne/ncell, i=1, ncell)], 1d-4 ), "wrong density " )
  end subroutine test_tbmodel_w90_hr



end module m_tb_model

