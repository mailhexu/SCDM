
#include "../common/abi_common.h"

module m_siesta_wannier_builder
    use defs_basis, only: dp, fnlen
    use m_scdm_math, only:  build_Rgrid
    use m_wannier_builder, only: WannierBuilder_t
    use m_siesta_wf_netcdf, only: siesta_wf_t
    use m_load_json, only: wannier_parameters_t
    implicit none
    private

    type, extends(WannierBuilder_t), public :: siesta_wannier_builder_t
        type(wannier_parameters_t) :: params
        character(len=500) :: wfk_fname
        type(siesta_wf_t) :: wf

        integer:: current_ik = -1
        real(dp), allocatable:: evals_cache(:)
        complex(dp), allocatable :: psi_k_cache(:, :)
    contains
        procedure :: read_config
        procedure :: set_wfk_file
        procedure :: run
        procedure :: finalize
     end type siesta_wannier_builder_t

contains


    subroutine finalize(self)
        class(siesta_wannier_builder_t), intent(inout) :: self
        call self%params%finalize()
    end subroutine finalize

    subroutine read_config(self, fname)
      class(siesta_wannier_builder_t), intent(inout) :: self
      character(len=*), intent(in) :: fname
      call self%params%read_config(fname)
    end subroutine read_config

    subroutine set_wfk_file(self, fname)
      class(siesta_wannier_builder_t), intent(inout) :: self
      character(:), intent(in) :: fname
      self%wfk_fname=trim(fname)
      call self%wf%loa
      !call self%WannierBuilder_t%intialize()
    end subroutine set_wfk_file

    function get_psi_k(self, ikpt)
      class(WannierBuilder_t), intent(inout):: self
      integer, intent(in):: ikpt
      complex(dp),  pointer:: psik(:, :)
      if (self%current_ik /= ikpt) then
         call self%wf%get_evecs_for_one_kpoint(self%psi_k_cache)
         self%current_ik = ikpt
      end if
      psik => self%psi_k_cache
    end function get_psi_k

    function get_evals_k(self, ikpt)
      class(WannierBuilder_t), intent(inout):: self
      integer, intent(in) :: ikpt
      complex(dp), pointer :: evals(:)
    end function get_evals_k

    subroutine run(self)
      class(siesta_wannier_builder_t), intent(inout) :: self
    end subroutine run

    subroutine test_siesta_wannier_builder_t(self)
        implicit none
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
