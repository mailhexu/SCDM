
#include "abi_common.h"
program main
  use defs_basis
  use m_siesta_wannier_builder, only: siesta_wannier_builder_t
  implicit none

  type(siesta_wannier_builder_t) :: swann
  call swann%load(config_fname="input.json", wfk_fname="wfk.nc")
  call swann%run_all(ncfilename="wannier.nc", Amnkfilename="Amnk.dat")
  call swann%finalize()
end program main
