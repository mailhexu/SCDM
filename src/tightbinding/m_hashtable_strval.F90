!!****m* ABINIT/m_hastable
!! NAME
!! m_hashtable
!!
!! FUNCTION
!! This module provide a string: value pair hash table
!! COPYRIGHT
!! Taken from http://fortranwiki.org/fortran/show/hash+table+example
!! The code is originally written by Izaak Beekman under the LGPL license.
!! Adapted for usage in Abinit by hexu
!!
!! Note: the behavior is different from the origial version
!! The value will be overwritten in this version, whereas it is ignored in the
!! original version if the key is already in the table (why??!!). 
!!
!! Note2:!!!!!!!!!!!!!!!!! FIXME
!! It does not handle white space at the end of string correctly. It does not affect
!! the usage in Multibinit but BE CAREFUL. 
!!
!! Below is the original Copyright.
!!=======================================
!! Copyright (C) Izaak Beekman 2010
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_hashtable_strval

!!***
  use defs_basis
  use m_test_utils, only: assertion

  IMPLICIT NONE ! Use strong typing
  ! a 
  integer, private, parameter :: nan=transfer(-225179999,1)
  INTEGER, PARAMETER :: tbl_size = 50

  TYPE sllist
     TYPE(sllist), POINTER :: child => NULL()
     CHARACTER(len=:), ALLOCATABLE :: key
     integer :: val
   CONTAINS
     PROCEDURE :: put  => put_sll
     PROCEDURE :: get  => get_sll
     PROCEDURE :: free => free_sll
     procedure :: print_all => print_all_sll
  END TYPE sllist

  TYPE hash_table_t
     TYPE(sllist), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                 :: vec_len = 0
     LOGICAL                                 :: is_init = .FALSE.
   CONTAINS
     PROCEDURE :: init => init_hash_table_t
     PROCEDURE :: put  => put_hash_table_t
     PROCEDURE :: get  => get_hash_table_t
     PROCEDURE :: free => free_hash_table_t
     PROCEDURE :: print_all => print_all_hash_table_t
     procedure :: has_key
     procedure :: put_int3
     procedure :: get_int3
     procedure :: has_key_int3
  END TYPE hash_table_t

  PUBLIC :: hash_table_t
  public :: test_hash_table
  public :: test_hash_table_int3

CONTAINS

  RECURSIVE SUBROUTINE put_sll(list,key,val)
    CLASS(sllist),    INTENT(inout) :: list
    CHARACTER(len=*), INTENT(in)    :: key
    integer, intent(in)  :: val
    INTEGER                         :: keylen

    keylen = LEN(key)
    IF (ALLOCATED(list%key)) THEN
       IF (list%key /= key) THEN
          IF ( .NOT. ASSOCIATED(list%child)) then
             ABI_MALLOC_SCALAR(list%child)
          end IF

          CALL put_sll(list%child,key,val)
       else
          list%val=val
       END IF
    ELSE
       IF (.NOT. ALLOCATED(list%key)) &
            ABI_MALLOC_TYPE_SCALAR(CHARACTER(len=keylen), list%key)
       list%key = key
       list%val = val
    END IF
  END SUBROUTINE put_sll


  RECURSIVE SUBROUTINE get_sll(list,key,val)
    CLASS(sllist),                 INTENT(in)    :: list
    CHARACTER(len=*),              INTENT(in)    :: key
    integer,                      INTENT(out)   :: val
    INTEGER                                      :: vallen

    vallen = 0
    val=nan
    IF (ALLOCATED(list%key) .AND. (list%key == key)) THEN
       val = list%val
    ELSE IF(ASSOCIATED(list%child)) THEN ! keep going
       CALL get_sll(list%child,key,val)
    ELSE ! At the end of the list, no key found
       return
    END IF
  END SUBROUTINE get_sll


  RECURSIVE SUBROUTINE free_sll(list)
    CLASS(sllist), INTENT(inout) :: list
    IF (ASSOCIATED(list%child)) THEN
       CALL free_sll(list%child)
       ABI_FREE_SCALAR(list%child)
    END IF
    list%child => NULL()
    ABI_SFREE(list%key)
  END SUBROUTINE free_sll


  recursive subroutine print_all_sll(self)
    class(sllist), intent(in) :: self

    if (allocated(self%key)) then
      write(*, "(A40, 1X, ES13.5)") self%key, self%val 
      if(associated(self%child)) then
        call self%child%print_all()
      endif
    end if
  end subroutine print_all_sll


  SUBROUTINE init_hash_table_t(tbl,tbl_len)
    CLASS(hash_table_t),   INTENT(inout) :: tbl
    INTEGER,     OPTIONAL, INTENT(in)    :: tbl_len

    ABI_SFREE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
       ABI_MALLOC(tbl%vec, (tbl_len))
       tbl%vec_len = tbl_len
    ELSE
       ABI_MALLOC(tbl%vec, (tbl_size))
       tbl%vec_len = tbl_size
    END IF
    tbl%is_init = .TRUE.
  END SUBROUTINE init_hash_table_t

  ! The first part of the hashing procedure using the string
  ! collating sequence
  ELEMENTAL FUNCTION sum_string(str) RESULT(sig)
    CHARACTER(len=*), INTENT(in)   :: str
    INTEGER                        :: sig
    CHARACTER, DIMENSION(LEN(str)) :: tmp
    INTEGER :: i

    FORALL (i=1:LEN(str))
       tmp(i) = str(i:i)
    END FORALL
    sig = SUM(ICHAR(tmp))
  END FUNCTION sum_string


  SUBROUTINE put_hash_table_t(tbl,key,val)
    CLASS(hash_table_t), INTENT(inout) :: tbl
    CHARACTER(len=*),    INTENT(in)    :: key
    integer,            INTENT(in)    :: val
    INTEGER                            :: hash

    hash = MOD(sum_string(key),tbl%vec_len) +1
    CALL tbl%vec(hash)%put(key=key,val=val)
  END SUBROUTINE put_hash_table_t


  function get_hash_table_t(tbl,key) result(val)
    CLASS(hash_table_t),           INTENT(in)    :: tbl
    CHARACTER(len=*),              INTENT(in)    :: key
    integer                        :: val
    INTEGER                                      :: hash
    hash = MOD(sum_string(key),tbl%vec_len) + 1
    CALL tbl%vec(hash)%get(key=key,val=val)
  END function get_hash_table_t


  SUBROUTINE free_hash_table_t(tbl)
    CLASS(hash_table_t), INTENT(inout) :: tbl    
    INTEGER     :: i, low, high

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    IF (ALLOCATED(tbl%vec)) THEN
       DO i=low,high
          CALL tbl%vec(i)%free()
       END DO
       ABI_FREE(tbl%vec)
    END IF
    tbl%is_init = .FALSE.
  END SUBROUTINE free_hash_table_t

  function has_key(self, key)
    class(hash_table_t), intent(in) :: self
    character(*), intent(in) :: key
    logical :: has_key
    has_key=(self%get(key)/=nan)
  end function has_key
  

  subroutine print_all_hash_table_t(self)
    class(hash_table_t), intent(in) :: self
    integer :: i, low, high
    low  = LBOUND(self%vec,dim=1)
    high = UBOUND(self%vec,dim=1) 

    if (allocated(self%vec)) then
       do i =low, high
          call self%vec(i)%print_all()
       end do
    end if
  end subroutine print_all_hash_table_t

  subroutine put_int3(self, key, val)
    class(hash_table_t), intent(inout) :: self
    integer, intent(in) :: key(3), val
    character(len=12) :: tmp
    call self%put(transfer(key, tmp), val)
  end subroutine put_int3

  function get_int3(self, key) result(val)
    class(hash_table_t), intent(inout) :: self
    integer, intent(in) :: key(3)
    integer :: val
    character(len=12) :: tmp
    val = self%get(transfer(key, tmp))
  end function get_int3


  function has_key_int3(self, key) result(val)
    class(hash_table_t), intent(inout) :: self
    integer, intent(in) :: key(3)
    logical :: val
    character(len=12) :: tmp
    val = self%has_key(transfer(key, tmp))
  end function has_key_int3


  !===============================================================
  !
  !> @
  !===============================================================
  subroutine test_hash_table()
    type(hash_table_t) :: dict
    integer :: val
    call dict%init(3)
    call dict%put("one", 1)
    call dict%put("two", 2)
    call dict%put("three", 3)
    call dict%put("four", 4)
    call dict%put("four", 40)

    val=dict%get("one")
    call assertion(val==1, " wrong result from hashtable")

    val=dict%get("two")
    call assertion(val==2, " wrong result from hashtable")

    val=dict%get("four")
    call assertion(val==40, " wrong result from hashtable")

    call assertion(dict%has_key("four"), "wrong has_key")
    call assertion(.not. dict%has_key(" four"), "wrong has_key")
    call assertion(.not. dict%has_key("four "), "wrong has_key")
    call assertion(.not. dict%has_key(" four"), "wrong has_key")
    call dict%free()
  end subroutine test_hash_table

  subroutine test_hash_table_int3()
    type(hash_table_t) :: dict
    call dict%init(3)
    call dict%put_int3([1,1,1], 111)
    call dict%put_int3([1,1,2], 112)
    call dict%put_int3([1,1,3], 113)
    call dict%put_int3([-1,1,3], -113)
    call dict%put_int3([-1,0,0], -113)
    call dict%put_int3([0,0,0], -113)

    call assertion(dict%get_int3([1,1,1])==111, "wrong dict int3")
    call assertion(dict%get_int3([1,1,2])==112, "wrong dict int3")
    call assertion(dict%get_int3([1,1,3])==113, "wrong dict int3")
    call assertion(dict%get_int3([-1,1,3])==-113, "wrong dict int3")
    call assertion(dict%has_key_int3([1,1,1]), "wrong dict int3 has_key")
  end subroutine test_hash_table_int3


end module m_hashtable_strval
