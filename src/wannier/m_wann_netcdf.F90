#include "abi_common.h"

module m_wann_netcdf
  use defs_basis
  use m_nctk, only: ab_define_var, netcdf_check
  use netcdf
  implicit none
  private

  type, public :: IOWannNC
     integer :: ncid
     integer :: d_nR, d_ndim, d_nwann, d_nbasis
     integer :: i_Rlist, i_wannR_real, i_wannR_imag, i_HwannR_real, i_HwannR_imag
     integer :: i_cell, i_numbers, i_masses, i_xred, i_xcart
   contains
     procedure :: write_wann
     procedure :: write_atoms
     procedure :: close
  end type IOWannNC

  contains

  subroutine write_wann( self,  filename, nR, ndim, nwann, nbasis, Rlist, WannR, HwannR)
    class(IOwannNC), intent(inout) :: self
    character(len=*),intent(in) :: filename
    integer, intent(in) :: nR, ndim, nwann, nbasis, Rlist(:,:)
    complex(dp) , intent(in) :: WannR(:,:,:), HwannR(:,:,:)
    ! id of variables
    integer :: ncerr

    ncerr = nf90_create(path=trim(filename), cmode=NF90_CLOBBER, ncid=self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when creating wannier netcdf  file")

    ! define dimensions
    ncerr=nf90_def_dim(self%ncid, "ndim", ndim , self%d_ndim)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension ndim in wannier nc file.")

    ncerr=nf90_def_dim(self%ncid, "wann_nwann", nwann,  self%d_nwann)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nwann in wannier nc file.")

    ncerr=nf90_def_dim(self%ncid, "wann_nbasis", nbasis,  self%d_nbasis)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nbasis in wannier nc file.")


    ncerr=nf90_def_dim(self%ncid, "nR", nR,  self%d_nR)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nR in wannier nc file.")

    ! define variables
    call ab_define_var(self%ncid, [self%d_ndim, self%d_nR], self%i_Rlist, &
         & NF90_INT, "Rlist", "List of R vectors", "dimensionless")

    call ab_define_var(self%ncid, [self%d_nbasis, self%d_nwann, self%d_nR], &
         & self%i_wannR_real, NF90_DOUBLE, "wann_wannier_function_real", "The real part of the Wannier function", "Angstrom")

    call ab_define_var(self%ncid, [self%d_nbasis, self%d_nwann, self%d_nR], &
         & self%i_wannR_imag, NF90_DOUBLE, "wann_wannier_function_imag", "The imaginary part of the Wannier function", "Angstrom")

    call ab_define_var(self%ncid, [self%d_nwann, self%d_nwann, self%d_nR], &
         & self%i_HwannR_real, NF90_DOUBLE, "wann_HamR_real", "The real part of the Wannier Hamiltonian", "eV")

    call ab_define_var(self%ncid, [self%d_nwann, self%d_nwann, self%d_nR], &
         & self%i_HwannR_imag, NF90_DOUBLE, "wann_HamR_imag", "The imaginary part of the Wannier Hamiltonian", "eV")
    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when ending def mode in wannier netcdf file")

    ! set variables
    
    ncerr = nf90_put_var(self%ncid, self%i_Rlist, Rlist, start=[1,1], count=[3, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting Rlist in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_wannR_real, real(real(wannR)), &
         & start=[1,1], count=[nbasis, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting WannR_real in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_wannR_imag, real(aimag(wannR)), &
         & start=[1,1], count=[nbasis, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting WannR_imag in wannier netcdf file.")


    ncerr = nf90_put_var(self%ncid, self%i_HwannR_real, real(real(HwannR)), &
         & start=[1,1], count=[nwann, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting HwannR_real in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_HwannR_imag, real(aimag(HwannR)), &
         & start=[1,1], count=[nwann, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting HwannR_imag in wannier netcdf file.")

  end subroutine write_wann

  subroutine close(self)
    class(IOwannNC), intent(inout) :: self
    integer :: ncerr
    ncerr = nf90_close(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error close wannier netcdf file.")
  end subroutine close

  subroutine write_atoms(self, natom, cell, numbers, masses, xred, xcart)
    class(IOwannNC), intent(inout) :: self
    character(len=*),intent(in) :: filename
    integer, intent(in) :: n_three, n_natom
    integer , intent(in) :: natom
    real(dp), intent(in) :: cell(:, :)
    integer, intent(in) :: numbers(:)
    real(dp), intent(in), masses(:), xred(:, :), xcart(:,:)

    ncerr = nf90_redef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error starting redef in wannier netcdf file")
    ncerr=nf90_def_dim(self%ncid, "three", 3, self%d_three)
    NCF_CHECK_MSG(ncerr, "Error defining three in wannier netcdf file")

    ncerr=nf90_def_dim(self%ncid, "natom", natom, self%d_natom)
    NCF_CHECK_MSG(ncerr, "Error defining natom in wannier netcdf file")

    !integer :: i_cell, i_numbers, i_masses, i_xred, i_xcart
    call ab_define_var(self%ncid, [self%d_three, self%d_natom], self%i_xred, NF90_DOUBLE, "Reduced positions for ")








    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error ending redef in wannier netcdf file")


  end subroutine write_atoms

end module m_wann_netcdf
