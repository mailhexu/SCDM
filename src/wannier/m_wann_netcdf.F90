module m_wann_netcdf
  use defs_basis
  use netcdf
  implicit none
  private

  type, public :: IOWannNC
     integer :: ncid
     integer :: d_nR, d_three, d_nwann, d_basis
     integer :: i_Rlist, i_wannR_real, i_wann_imag, i_HR_real, i_HR_imag
   contains
     procedure :: initialize
     procedure :: finalize
  end type IOWannNC

  contains

  subroutine write_Hwann( filename, nR, nwann, nbasis, Rlist, wannR, HR)
    character(len=*),intent(in) :: filename
    integer, intent(in) :: nR, nwann, Rlist(:,:)
    ! id of variables
    integer :: ncerr
    ncerr = nf90_create(path=trim(filename), cmode=NF90_CLOBBER, ncid=self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when creating wannier netcdf  file")

    ! define dimensions
    ncerr=nf90_def_dim(self%ncid, "three", 3, self%d_three)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension three in wannier nc file.")

    ncerr=nf90_def_dim(self%ncid, "nwann", nwann,  self.d_nwann)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nwann in wannier nc file.")

    ncerr=nf90_def_dim(self%ncid, "nbasis", nbasis,  self.d_basis)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nbasis in wannier nc file.")


    ncerr=nf90_def_dim(self%ncid, "nR", nR,  self%d_nR)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nR in wannier nc file.")

    ! define variables
    call ab_define_var(self%ncid, [self%d_three, self%d_nR], self%i_Rlist, NF90_INT, "Rlist", "List of R vectors", "dimensionless")

    call ab_define_var(self%ncid, [self%d_nbasis, self%d_nwann, self%d_nR], self%i_wannR_real, NF90_DOUBLE, "WannR_real", "The real part of the Wannier function", "Angstrom")

    call ab_define_var(self%ncid, [self%d_nbasis, self%d_nwann, self%d_nR], self%i_wannR_imag, NF90_DOUBLE, "WannR_imag", "The imaginary part of the Wannier function", "Angstrom")

    call ab_define_var(self%ncid, [self%d_nwann, self%d_nwann, self%d_nR], self%i_HR_real, NF90_DOUBLE, "HR_real", "The real part of the Wannier Hamiltonian", "eV")

    call ab_define_var(self%ncid, [self%d_nwann, self%d_nwann, self%d_nR], self%i_HR_imag, NF90_DOUBLE, "HR_imag", "The imaginary part of the Wannier Hamiltonian", "eV")

    ncerr = nf90_put_var(self%ncid, self%id_Rlist, Rlist, start=[1,1], count=[three, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting Rlist in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%id_wannR_real, real(real(wannR)), start=[1,1], count=[nbasis, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting WannR_real in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%id_wannR_imag, real(aimag(wannR)), start=[1,1], count=[nbasis, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting WannR_imag in wannier netcdf file.")



    ! set variables

    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when ending def mode in wannier netcdf file")
  end subroutine write_Hwann

  
end module m_wann_netcdf
