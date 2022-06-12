#include "abi_common.h"

module m_siesta_wannier_main
  use defs_basis
  use m_load_json, only: wannier_parameters_t
  use m_siesta_wannier_builder, only: siesta_wannier_builder_t
  implicit none
  private

  public :: siesta_wannier_main_t

  type siesta_wannier_main_t
  end type siesta_wannier_main_t
contains



end module m_siesta_wannier_main
