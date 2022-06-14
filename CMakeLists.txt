cmake_minimum_required(VERSION 3.1)
enable_language(Fortran)
set(cmake_dir "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

list(INSERT CMAKE_MODULE_PATH 0 "${cmake_dir}")

project(scdm)


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

include(FindNetCDF)
set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDES})


include_directories(${MPI_Fortran_INCLUDE_PATH})
add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

add_definitions(${MPI_Fortran_COMPILE_FLAGS})
link_directories(${MPI_Fortran_LIBRARIES})


# libjson
set (Json_Fortran_INLCUDE_PATH /home/hexu/.local/jsonfortran-gnu-8.3.0/lib)
set (Json_Fortran_LIBRARY_PATH /home/hexu/.local/jsonfortran-gnu-8.3.0/lib)
set (Json_Fortran_LIBRARIES jsonfortran)
include_directories(${Json_Fortran_INLCUDE_PATH})
link_directories(${Json_Fortran_LIBRARY_PATH})

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mods)

set(HEADERS
    "src/common/*.h"
   )

file(GLOB sourcefile
    "src/*.F90"
    "src/common/*.F90"
    "src/wannier/*.F90"
    "src/tightbinding/*.F90"
    "src/phonon/*.F90"
    "src/unittests/*.F90"
    "src/siesta/*.F90"
    "external/*.F90"
	"src/json/*.F90"
    )

list(REMOVE_ITEM sourcefile
    "${CMAKE_CURRENT_SOURCE_DIR}/src/unittests/unittest.F90"
)


add_library(scdm SHARED
    ${sourcefile}
    )

target_include_directories(scdm PUBLIC
    "src/common"
    )
target_link_libraries(scdm ${lapackblas_libraries} ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARIES_F90} ${Json_Fortran_LIBRARIES})
#target_link_libraries(scdm ${lapackblas_libraries}  ${NETCDF_LIBRARIES_F90})


target_link_libraries(scdm )

if(MPI_COMPILE_FLAGS)
  set_target_properties(scdm PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(scdm PROPERTIES
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
endif()

add_executable(unittest
    src/unittests/unittest.F90
  )


#target_link_libraries(unittest scdm ${lapackblas_libraries} ${NETCDF_LIBRARIES_F90})
target_link_libraries(unittest scdm)


add_executable(sdown
	src/main/main.F90
	)	

target_link_libraries(sdown scdm)

#add_executable(testscdm
#  src/m_math.F90
#  src/unittest.F90
#  )

#target_link_libraries(testscdm ${lapackblas_libraries} ${MPI_Fortran_LIBRARIES})
