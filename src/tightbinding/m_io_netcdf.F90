#include "abi_common.h"

module m_io_netcdf
  use defs_basis
  use netcdf
  implicit none
  private

  public :: read_tb_netcdf
  public :: save_band

  public :: test_read_tb_netcdf
  public :: test_save_band


contains

  subroutine check(status, message)
    integer, intent(in) :: status
    character (*), intent(in):: message
    if(status/=nf90_noerr) then
       ABI_ERROR(message)
    end if
  end subroutine check

  subroutine read_tb_netcdf(fname, Rlist, HamR_real, HamR_imag)
    character(len=*), intent(in) :: fname
    integer, allocatable, intent(inout) :: Rlist(:, :)
    real(dp), allocatable, intent(inout) :: HamR_real(:,:,:), HamR_imag(:,:,:)
    integer :: ncid, nwann, nR, nwann_id, nR_id, hamr_real_id, hamr_imag_id, Rlist_id
    ! open nc file
    call check(nf90_open(trim(fname), NF90_NOWRITE, ncid), &
         & "openning nc file"//fname)
    ! read dimensions
    call check(nf90_inq_dimid(ncid, 'wann_nwann', nwann_id), &
         & "getting nwann dim id")
    call check(nf90_inquire_dimension(ncid, nwann_id, len=nwann), &
         & "getting nwann")

    call check(nf90_inq_dimid(ncid, 'nR', nR_id), &
         & "getting nR dim id")
    call check(nf90_inquire_dimension(ncid, nR_id, len=nR), &
         & "getting nR")

    ! intialize
    ABI_MALLOC(Rlist, (3, nR))
    ABI_MALLOC(HamR_real,(nwann, nwann, nR))
    ABI_MALLOC(HamR_imag,(nwann, nwann, nR))

    ! read Rlist and HamR
    call check(nf90_inq_varid(ncid, "wann_Rlist", Rlist_id), &
         & "inquire Rlist id")
    call check(nf90_get_var(ncid, Rlist_id, Rlist), &
         & "reading var Rlist")

    call check(nf90_inq_varid(ncid, "wann_HamR_real", hamr_real_id), &
         & "inquire hamR_real id")
    call check(nf90_get_var(ncid, hamr_real_id, HamR_real), &
         & "reading var HamR_real")

    call check(nf90_inq_varid(ncid, "wann_HamR_imag", hamr_imag_id), &
         & "inquire wann_HamR_imag_id")
    call check(nf90_get_var(ncid, hamr_imag_id, HamR_imag), &
         & "reading var wann_HamR_imag")

    call check(nf90_close(ncid), "closing nc file "//fname)

  end subroutine read_tb_netcdf

  subroutine save_band(fname, kpoints, evals, evecs)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: kpoints(:,:)
    real(dp), intent(in) :: evals(:,:)
    complex(dp) , intent(in):: evecs(:,:,:)
    integer :: nk, norb, nband
    integer :: ncid,  norb_dimid, nband_dimid, nk_dimid, ndim_dimid
    integer :: kpoints_id, evals_id, evecs_real_id, evecs_imag_id

    norb=size(evecs, 1)
    nband=size(evecs, 2)
    nk=size(evecs, 3)

    ! open file
    call check(nf90_create(trim(fname), NF90_CLOBBER, ncid), "openning nc file"//fname)
    ! define dimensions
    call check(nf90_def_dim(ncid, "ndim", 3, ndim_dimid), "defining nc dim ndim")
    call check(nf90_def_dim(ncid, "norb", norb, norb_dimid), "defining nc dim norb")
    call check(nf90_def_dim(ncid, "nband", nband, nband_dimid), "defining nc dim nband")
    call check(nf90_def_dim(ncid, "nk", nk, nk_dimid), "defining nc dim nk")
    ! define vars
    call check(nf90_def_var(ncid, "kpoints", NF90_DOUBLE, [ndim_dimid, nk_dimid], kpoints_id), &
         & "defining nc var kpoints")
    call check(nf90_def_var(ncid, "eigen_values", NF90_DOUBLE, [nband_dimid, nk_dimid], evals_id), &
         & "defining nc var kpoints")
    call check(nf90_def_var(ncid, "eigen_vectors_real", NF90_DOUBLE, [norb_dimid, nband_dimid, nk_dimid], evecs_real_id), &
         &"defining nc var evecs_real")
    call check(nf90_def_var(ncid, "eigen_vectors_imag", NF90_DOUBLE, [norb_dimid, nband_dimid, nk_dimid], evecs_imag_id), &
         &"defining nc var evecs_imag")

    ! writting vars
    call check(nf90_put_var(ncid, kpoints_id, kpoints), &
         & "writing nc variable kpoints")
    call check(nf90_put_var(ncid, evals_id, evals), &
         & "writing nc variable eigen_values")
    call check(nf90_put_var(ncid, evecs_real_id, real(evecs)), &
         & "writing nc variable real part of eigen_vectors")
    call check(nf90_put_var(ncid, evecs_imag_id, aimag(evecs)), &
         & "writing nc variable imag part of eigen_vectors")

    ! close

    call check(nf90_close(ncid), "closing nc file "//fname)
  end subroutine save_band

  subroutine save_densiy(fname, den)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: den(:)
    integer :: ncid, norb_dimid, den_id, norb
    norb=size(den, 1)
    call check(nf90_create(trim(fname), NF90_CLOBBER, ncid), "openning nc file"//fname)
    call check(nf90_def_dim(ncid, "norb", norb, norb_dimid), "defining nc dim norb")
    call check(nf90_def_var(ncid, "density", NF90_DOUBLE, [ norb_dimid], den_id), &
         & "defining nc var kpoints")
    call check(nf90_put_var(ncid, den_id, den), &
         & "writing nc variable den")
    call check(nf90_close(ncid), "closing nc file "//fname)
  end subroutine save_densiy

!####################### Unit tests #############################

  subroutine test_read_tb_netcdf()
  end subroutine test_read_tb_netcdf

  subroutine test_save_band()
    integer,parameter :: nk=3, norb=2
    real(dp) :: kpoints(3, nk), evals(norb, nk)
    complex(dp) :: evecs(norb,norb,nk)

    kpoints(:,1) = [0.0, 0.0,0.0]
    kpoints(:,1) = [0.5, 0.0,0.0]

    evals(:, 1) =[-2.0, 1.0]
    evals(:, 2) =[-3.0, 4.0]

    evecs=cmplx(0.0,0.0, kind=dp)
    evecs(1, 2, 1) =cmplx(1.2, 3.5, kind=dp)

    call save_band(fname='test_band.nc', kpoints=kpoints, evals=evals, evecs=evecs )
  end subroutine test_save_band


end module m_io_netcdf


