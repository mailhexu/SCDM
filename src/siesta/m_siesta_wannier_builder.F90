
#include "../common/abi_common.h"

module m_siesta_wannier_builder
  use defs_basis, only: dp, fnlen
  use m_scdm_math, only:  build_Rgrid
  use m_wannier_builder, only: WannierBuilder_t
  use m_siesta_wf_netcdf, only: siesta_wf_t
  use m_load_json, only: wannier_parameters_t
  use m_wann_netcdf, only: IOWannNC
  use m_parameters, only: atomic_masses
  implicit none
  private

  type, extends(WannierBuilder_t), public :: siesta_wannier_builder_t
     type(wannier_parameters_t) :: params
     character(len=500) :: wfk_fname
     type(siesta_wf_t) :: wf
     integer:: current_ik = -1
     real(dp), allocatable:: evals_cache(:)
     complex(dp), pointer:: psi_k_cache(:, :)
     complex(dp), pointer :: Sk_cache(:, :)
   contains
     procedure :: load
     procedure :: read_config
     procedure :: set_wfk_file
     procedure :: get_psi_k
     procedure :: get_Sk
     procedure :: get_evals_k
     procedure :: run_all
     procedure :: finalize
  end type siesta_wannier_builder_t

contains

  subroutine load(self, config_fname, wfk_fname)
    class(siesta_wannier_builder_t), intent(inout) :: self
    character(len=*), intent(in) :: config_fname
    character(len=*), intent(in) :: wfk_fname
    integer, allocatable ::  Rlist(:,:)
    integer :: disentangle_func_type

    call self%read_config(config_fname)
    call self%set_wfk_file(wfk_fname)

    call build_Rgrid(self%params%kmesh, Rlist)

    !if (self%params%disentangle_func_type(1)=='G' .or. self%params%disentangle_func_type(1)=='g'):
    select case(self%params%disentangle_func_type(1:1))
    case ('g', 'G')
       disentangle_func_type=3
    case('f', 'F')
       disentangle_func_type=2
    case('u', 'U')
       disentangle_func_type=1
    case default
       disentangle_func_type=1
    end select

    call self%WannierBuilder_t%initialize(self%wf%kpts, self%wf%kweights, Rlist, self%params%nwann, &
         & self%wf%nbasis, self%wf%nband, &
         &  disentangle_func_type, self%params%mu, self%params%sigma, self%params%exclude_bands, &
         & self%params%project_to_anchor, self%params%method)

    ABI_MALLOC(self%psi_k_cache, (self%nbasis, self%nband))
    ABI_MALLOC(self%Sk_cache, (self%nbasis, self%nbasis))

 if (self%params%method ==1 ) then
    !if(allocated(self%params%anchor_ibands)) then
       if(self%params%anchor_ibands(1) > 0) then
          call self%set_anchor( self%params%anchor_kpt, self%params%anchor_ibands)
       !endif
    else
       call self%set_anchor(anchor_kpt = self%params%anchor_kpt)
    end if

    print *, "anchor_ibands: ", self%anchor_ibands

    ! output information
 else if(self%params%method == 2) then
    ! TODO: projected lattice wannier function
     print *,  ' Constructing LWF with projected wannier function method.'
     stop "Projected Wannier function not implemented yet. "
    !call self%set_disp_projector(params%lwf_projector)
    !write(msg, '(2a)')  ' The projectors: ', trim(ltoa(dtset%lwf_projector))
  end if
  end subroutine load

  subroutine finalize(self)
    class(siesta_wannier_builder_t), intent(inout) :: self
    call self%params%finalize()
    call self%wf%finalize()
    call self%WannierBuilder_t%finalize()
    ABI_FREE(self%Sk_cache)
    ABI_FREE(self%psi_k_cache)
    Nullify(self%psi_k_cache)
  end subroutine finalize

  subroutine read_config(self, fname)
    class(siesta_wannier_builder_t), intent(inout) :: self
    character(len=*), intent(in) :: fname
    call self%params%read_config(fname)
  end subroutine read_config

  subroutine set_wfk_file(self, fname)
    class(siesta_wannier_builder_t), intent(inout) :: self
    character(len=*), intent(in) :: fname
    self%wfk_fname=trim(fname)
    call self%wf%read_from_netcdf(fname, read_evecs=.False.)
  end subroutine set_wfk_file

  function get_psi_k(self, ikpt) result (psik)
    class(siesta_wannier_builder_t),   intent(inout):: self
    integer, intent(in):: ikpt
    complex(dp), pointer:: psik(:, :)
    if (self%current_ik /= ikpt) then
       call self%wf%get_evecs_for_one_kpoint(ikpt, self%psi_k_cache)
       self%current_ik = ikpt
    end if
    psik => self%psi_k_cache
  end function get_psi_k

  function get_evals_k(self, ikpt) result (ek)
    class(siesta_wannier_builder_t),  intent(inout):: self
    integer, intent(in):: ikpt
    real(dp), pointer:: ek(:)
    ek=> self%wf%eigenvalues(:, ikpt)
  end function get_evals_k

  function get_Sk(self, ikpt) result (Sk)
    class(siesta_wannier_builder_t),  intent(inout):: self
    integer, intent(in):: ikpt
    complex(dp), pointer :: Sk(:, :)
    call self%wf%get_Sk(ikpt, self%Sk_cache)
    Sk => self%Sk_cache
  end function get_Sk

  
  subroutine run_all(self, ncfilename, Amnkfilename)
    class(siesta_Wannier_Builder_t), intent(inout):: self
    character(*), intent(in):: ncfilename
    character(*), intent(in):: Amnkfilename
    type(IOWannNC):: ncfile
    real(dp) :: masses(self%wf%natom)
    integer :: i
    call self%construct_wannier()
    call self%create_ncfile(ncfilename, ncfile)
    call self%write_wann_netcdf( ncfile,   &
         &wannR_unit='dimensionless', HwannR_unit='eV')
    do i=1, self%wf%natom
       masses(i) = atomic_masses(self%wf%atomic_numbers(i))
    end do
    print *, "writting masses"
    call ncfile%write_atoms(natom=self%wf%natom, cell=self%wf%cell, &
         & numbers=self%wf%atomic_numbers, masses=masses, xred=self%wf%xred)
    call self%close_ncfile(ncfile)
    call self%write_Amnk_w90(trim(Amnkfilename))
  end subroutine run_all


  subroutine test_siesta_wannier_builder_t(self)
    type(siesta_wannier_builder_t) :: self
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
  end subroutine test_siesta_wannier_builder_t

end module m_siesta_wannier_builder
