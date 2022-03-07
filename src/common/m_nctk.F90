#include "abi_common.h"

module m_nctk
use netcdf
implicit none
private

public :: netcdf_check
public:: ab_define_var
contains

subroutine netcdf_check(ierr, msg)
      integer, intent ( in) :: ierr
      character(*), intent(in) :: msg

      if(ierr/= nf90_noerr) then
      print *, trim(nf90_strerror(ierr))//msg
      stop "Stopped"
      end if

end subroutine netcdf_check

subroutine ab_define_var(ncid, var_dim_id, var_id, var_type, var_name, var_mnemo, var_units)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ncid
 integer, intent(out) :: var_id
 character(len=*), intent(in) :: var_mnemo,var_units,var_name
 integer,intent(in) :: var_type
!arrays
 integer,intent(in) :: var_dim_id(:)

!Local variables-------------------------------
!scalars
 integer :: ncerr

! *************************************************************************

!#ifdef HAVE_NETCDF
 ncerr = nf90_def_var(ncid, trim(var_name), var_type, var_dim_id, var_id)
 NCF_CHECK_MSG(ncerr," define variable "//trim(var_name))

 ncerr = nf90_put_att(ncid, var_id,  "units",trim(var_units))
 NCF_CHECK_MSG(ncerr," define attribute for "//trim(var_name))

 ncerr = nf90_put_att(ncid, var_id,  "mnemonics", trim(var_mnemo))
 NCF_CHECK_MSG(ncerr," define attribute for "//trim(var_name))
!#endif

end subroutine ab_define_var
!!***

end module m_nctk
