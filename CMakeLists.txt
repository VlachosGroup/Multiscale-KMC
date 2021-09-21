cmake_minimum_required(VERSION 3.10.2)
project(msa_kmc)

include_directories(include)
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_executable(msa-kmc ${SOURCES})

# Optional compilation of parallel MPI version
option(MPI "MPI" OFF)
IF(MPI)
  MESSAGE("Building parallel version with MPI")
  target_compile_definitions(msa-kmc PUBLIC -DMPI)
ENDIF()
