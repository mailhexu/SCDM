
#include "abi_common.h"
program main
  use defs_basis
  use m_siesta_wannier_builder, only: siesta_wannier_builder_t
  implicit none

  type(siesta_wannier_builder_t) :: swann
  call swann%read_config(fname="input.json")
  !call swann%run()
  call swann%finalize()
end program main
