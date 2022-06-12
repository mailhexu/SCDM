
#include "abi_common.h"


module m_load_json

  use defs_basis
  use json_module

  implicit none

  private

  public :: read_config
  public :: test_read_json
  public :: wannier_parameters_t

  type wannier_parameters_t
     integer :: nwann = 0
     integer :: method = 0
     integer, allocatable :: kmesh(:)
     CHARACTER(len=:), allocatable :: disentangle_func_type
     integer, allocatable :: exclude_bands(:)
     integer :: n_exclude_bands = 0
     real(dp):: mu = 0.0_dp
     real(dp) :: sigma = 0.0_dp
     logical :: project_to_anchor = .False.
     real(dp), allocatable :: anchor_kpt(:)
     integer, allocatable :: anchor_ibands(:)
   contains
     procedure :: read_config
     procedure :: print
     procedure :: finalize
  end type wannier_parameters_t

contains

  subroutine print(self)
    class(wannier_parameters_t), intent(inout) :: self
    print *, "nwann:", self%nwann
    print *, "method:", self%method
    print *, "kmesh:", self%kmesh
    print *, "disentangle_func_type:", self%disentangle_func_type
    print *, "mu:", self%mu
    print *, "sigma", self%sigma
    print *, "n_exclude_bands:", self%n_exclude_bands
    if (self%n_exclude_bands>0) then
       print *, "exclude_bands", self%exclude_bands
    end if
    print *, "project_to_anchor", self%project_to_anchor
    print *, "anchor_kpt:", self%anchor_kpt
    print *, "anchor_ibands:", self%anchor_ibands
  end subroutine print


  subroutine read_config(params, filename)
    class(wannier_parameters_t), intent(inout) :: params
    character(len=*), intent(in) :: filename
    type(json_file) :: json
    logical :: found

    call json%initialize()
    call json%load(trim(filename))
    print *, trim(filename)
    call json%print()

    call json%get('nwann', params%nwann, found)
    if ( .not. found ) stop "nwann not found in json config file."

    call json%get('method', params%method, found)
    if ( .not. found ) stop "method not found in json config file."
    if (params%method < 1 .or. params%method>2) then
       stop "method should be 1 (scdm-k) or 2 (projected-wannier-function)."
    end if

    call json%get('kmesh', params%kmesh, found)
    if ( .not. found ) then
       print *, "kmesh not found in json config file, use default value 5 5 5."
       ABI_MALLOC(params%kmesh, (3))
       params%kmesh(:) =5
    end if


    call json%get('disentangle_func_type', params%disentangle_func_type, found)
    if ( .not. found ) then
       stop "disentangle_func_type not found in json config file. Allowed values: unity, gaussian, fermi"
    end if


    call json%get('exclude_bands', params%exclude_bands, found)
    if ( .not. found ) then
       params%n_exclude_bands=0
    else
       params%n_exclude_bands=size(params%exclude_bands)
    end if

    call json%get('mu', params%mu, found)
    if ( .not. found ) then
       print *, "mu not found in json config file, use default value 0.0 ."
    end if

    call json%get('sigma', params%sigma, found)
    if ( .not. found ) then
       print *, "sigma not found in json config file, use default value 0.0 ."
    end if


    call json%get('project_to_anchor', params%project_to_anchor, found)

    call json%get('anchor_kpt', params%anchor_kpt, found)
    if ( .not. found ) then
      allocate(params%anchor_kpt(3))
      params%anchor_kpt(:) = 0.0_dp
      print *, "anchor_kpt not found in json config file, use default Gamma point [0.0,0.0,0.0]"
    end if


    call json%get('anchor_ibands', params%anchor_ibands, found)
    if ( .not. found ) then
       if (params%project_to_anchor) then
          print *, "anchor_ibands not found in json config file, will decide automatically."
          allocate(params%anchor_ibands(params%nwann))
          params%anchor_ibands(:) = -1
       end if
    end if
    call json%destroy()
    if (json%failed()) stop "Parse Json file failed."
    call params%print()
  end subroutine read_config

  subroutine finalize(params)
    class(wannier_parameters_t), intent(inout) :: params
    ABI_SFREE(params%disentangle_func_type)
    ABI_SFREE(params%exclude_bands)
    ABI_SFREE(params%kmesh)
    ABI_SFREE(params%anchor_kpt)
    ABI_SFREE(params%anchor_ibands)
  end subroutine finalize




  subroutine test_read_json()
    type(json_file) :: json
    logical :: found
    integer :: int_var, missing_int_var
    real(8) :: float_var
    integer, allocatable :: int_list_var(:)
    real(8), allocatable :: float_list_var(:)
    CHARACTER(len=:), allocatable :: string_var

    ! initialize the class
    call json%initialize()

    ! read the file
    call json%load(filename = 'input.json')

    ! print the file to the console
    call json%print()

    ! extract data from the file
    ! [found can be used to check if the data was really there]
    call json%get('int_var', int_var, found)
    if ( .not. found ) stop 1
    call json%get('float_var', float_var, found)
    if ( .not. found ) stop 1

    call json%get('int_list_var', int_list_var, found)
    if ( .not. found ) stop 1

    call json%get('float_list_var', float_list_var, found)
    if ( .not. found ) stop 1

    call json%get('string_var', string_var, found)
    if ( .not. found ) stop 1


    ! missing value: use default

    call json%get('missing_int_var', missing_int_var, found)
    if (.not. found) missing_int_var=-99

    !call json%get('data(1).number', k, found)
    !if ( .not. found ) stop 1

    print *, "int_var: ", int_var
    print *, "missing_int_var: ", int_var
    print *, "float_var: ", float_var
    print *, "int_list_var: ", int_list_var
    print *, "float_list_var: ", float_list_var
    print *, "string_var: ", string_var

    ! clean up

    deallocate(int_list_var)
    deallocate(float_list_var)
    deallocate(string_var)

    call json%destroy()
    if (json%failed()) stop 1
  end subroutine test_read_json

end module m_load_json
