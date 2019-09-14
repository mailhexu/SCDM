cmake_minimum_required(VERSION 3.1)
project(scdm)
enable_language(Fortran)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

if(NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR "Could not find Fortran MPI.  Please set MPI_Fortran_COMPILER to point to the mpifort wrapper.")
endif()


find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

#include_directories(${MPI_Fortran_INCLUDE_PATH})
#add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

#add_definitions(${MPI_Fortran_COMPILE_FLAGS})
#link_directories(${MPI_Fortran_LIBRARIES})

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mods)

add_library(scdm SHARED
    src/m_scdm.F90 
    )

target_link_libraries(scdm ${lapackblas_libraries} ${MPI_Fortran_LIBRARIES})


target_link_libraries(scdm )

if(MPI_COMPILE_FLAGS)
  set_target_properties(scdm PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(scdm PROPERTIES
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
endif()


add_executable(testscdm
  src/m_scdm.F90
  src/test.F90
  )


target_link_libraries(testscdm ${lapackblas_libraries} ${MPI_Fortran_LIBRARIES})
