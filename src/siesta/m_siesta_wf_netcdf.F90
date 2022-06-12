#include "abi_common.h"

module m_siesta_wf_netcdf
  use defs_basis
  use netcdf
  implicit none
  private

  type, public:: ncfile_t
     integer :: ncid
   contains
     procedure :: open_ncfile
     procedure :: close
  end type ncfile_t

  type, public :: siesta_wf_t
     integer :: nkpt
     integer :: nband
     integer :: nbasis
     integer :: natom
     integer :: norb
     integer :: spintype
     integer :: three=3

     integer, allocatable :: atomic_numbers(:)
     real(dp) :: cell(3,3)
     real(dp), allocatable :: xred(:, :)
     real(dp), allocatable :: kpts(:, :)
     real(dp), allocatable :: kweights(:)
     real(dp), allocatable :: eigenvalues(:, :)
     complex(dp), allocatable :: eigenvectors(:, :, :)
     type(ncfile_t) :: wffile
   contains
     procedure :: read_from_netcdf
     procedure :: get_evecs_for_one_kpoint
     procedure :: finalize
  end type siesta_wf_t


contains

  subroutine check(status, message)
    integer, intent(in) :: status
    character (*), intent(in):: message
    if(status/=nf90_noerr) then
       ABI_ERROR(message)
    end if
  end subroutine check

  subroutine get_dim(ncid, name, var)
    integer, intent(inout) :: ncid, var
    character(len=*) :: name
    integer :: dim_id
    call check(nf90_inq_dimid(ncid, name, dim_id), &
         & "getting" // name// "dim id")
    call check(nf90_inquire_dimension(ncid, dim_id, len=var), &
         & "getting "//name)
  end subroutine get_dim


  subroutine open_ncfile(ncfile, fname)
    class(ncfile_t), intent(inout) :: self
    character(len=*), intent(in) :: fname
    call check(nf90_open(trim(fname), NF90_NOWRITE, self%ncid), &
         & "openning nc file"//fname)

  end subroutine open_ncfile

  subroutine close(ncfile)
    class(ncfile_t), intent(inout) :: self
    call check(nf90_close(self%ncid), "closing nc file "//fname)
  end subroutine close


  subroutine read_from_netcdf(wf, fname, read_evecs)
    ! Note that eigenvectors are not read.
    class(siesta_wf_t), intent(inout) :: wf
    character(len=*), intent(in) :: fname
    logical :: read_evecs
    real(dp), allocatable :: evecs_real(:,:,:), evecs_imag(:,:,:)
    integer :: ncid,   var_id
    ! open nc file
    call check(nf90_open(trim(fname), NF90_NOWRITE, ncid), &
         & "openning nc file"//fname)
    ! read dimensions

    call get_dim(ncid, "nkpt", wf%nkpt )
    call get_dim(ncid, "nband", wf%nband)
    call get_dim(ncid, "nbasis", wf%nbasis)
    call get_dim(ncid, "natom", wf%natom)
    call get_dim(ncid, "norb", wf%norb)
    call get_dim(ncid, "spintype", wf%spintype)
    call get_dim(ncid, "three", wf%three)

    ABI_MALLOC(wf%atomic_numbers, (wf%natom))
    ABI_MALLOC(wf%xred, (3, wf%natom))
    ABI_MALLOC(wf%kpts, (3, wf%nkpt))
    ABI_MALLOC(wf%kweights, (wf%nkpt))
    ABI_MALLOC(wf%eigenvalues, (wf%nband, wf%nkpt))
    
    if(read_evecs) then
        ABI_MALLOC(evecs_real,(wf%nbasis, wf%nband, wf%nkpt))
        ABI_MALLOC(evecs_imag,(wf%nbasis, wf%nband, wf%nkpt))
        ABI_MALLOC(wf%eigenvectors,(wf%nbasis, wf%nband, wf%nkpt))
    end if

    ! read vars
    ! atomic_numbers
    call check(nf90_inq_varid(ncid, "atomic_numbers", var_id), &
         & "inquire atomic_numbers id")
    call check(nf90_get_var(ncid, var_id, wf%atomic_numbers), &
         & "reading var atomic_numbers")

    ! cell
    call check(nf90_inq_varid(ncid, "cell", var_id), &
         & "inquire cell id")
    call check(nf90_get_var(ncid, var_id, wf%cell), &
         & "reading var cell")

    ! xred
    call check(nf90_inq_varid(ncid, "xred", var_id), &
         & "inquire xred id")
    call check(nf90_get_var(ncid, var_id, wf%xred), &
         & "reading var xred")

    ! kweights
    call check(nf90_inq_varid(ncid, "kweights", var_id), &
         & "inquire kweights id")
    call check(nf90_get_var(ncid, var_id, wf%kweights), &
         & "reading var kweights")

    ! eigenvalues
    call check(nf90_inq_varid(ncid, "eigenvalues", var_id), &
         & "inquire eigenvalues id")
    call check(nf90_get_var(ncid, var_id, wf%eigenvalues), &
         & "reading var eigenvalues")

    ! eigenvectors_real
    call check(nf90_inq_varid(ncid, "eigenvectors_real", var_id), &
         & "inquire eigenvectors_real id")
    call check(nf90_get_var(ncid, var_id, evecs_real), &
         & "reading var eigenvectors_real")

    ! eigenvectors_imag
    call check(nf90_inq_varid(ncid, "eigenvectors_imag", var_id), &
         & "inquire eigenvectors_imag id")
    call check(nf90_get_var(ncid, var_id, evecs_imag), &
         & "reading var eigenvectors_imag")

    if(read_evecs) then
       wf%eigenvectors = CMPLX(evecs_real, evecs_imag, kind=dp)
       ABI_SFREE(evecs_real)
       ABI_SFREE(evecs_imag)
    end if

    call check(nf90_close(ncid), "closing nc file "//fname)

  end subroutine read_from_netcdf

  subroutine get_evecs_for_one_kpoint(wf,  ik,  evecs)
    class(siesta_wf_t), intent(inout) :: wf
    integer, intent(in) ::  ik
    complex(dp), intent(inout) :: evecs(:,:)
    integer :: nc_id,var_id

    real(dp) :: evecs_real(wf%nbasis, wf%nband)
    real(dp) :: evecs_imag(wf%nbasis, wf%nband)

    nc_id = wf%wffile%ncid
    call check(nf90_inq_varid(ncid, "eigenvectors_real", var_id), &
         & "inquire eigenvectors_real id")
    call check(nf90_get_var(ncid, var_id, evecs_real, start=[1,1, ik], count=[wf%nbasis, wf%nband, 1]), &
         & "reading var eigenvectors_real")

    call check(nf90_inq_varid(ncid, "eigenvectors_imag", var_id), &
         & "inquire eigenvectors_imag id")
    call check(nf90_get_var(ncid, var_id, evecs_imag, start=[1,1, ik], count=[wf%nbasis, wf%nband, 1]), &
         & "reading var eigenvectors_imag")
    evecs(:,:) = cmplx(evecs_real, evecs_imag, kind=dp)
  end subroutine get_evecs_for_one_kpoint


  subroutine finalize(wf)
    class(siesta_wf_t), intent(inout) :: wf
    ABI_SFREE(wf%atomic_numbers)
    ABI_SFREE(wf%xred)
    ABI_SFREE(wf%kpts)
    ABI_SFREE(wf%kweights)
    ABI_SFREE(wf%eigenvalues)
    ABI_SFREE(wf%eigenvectors)
  end subroutine finalize


end module m_siesta_wf_netcdf
