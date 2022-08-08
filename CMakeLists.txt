cmake_minimum_required(VERSION 3.10.2)
project(msa_kmc)


add_library(msa-kmcl STATIC
	src/KMC_traj.cpp
	src/Traj_stats.cpp
	src/KMC_traj_TTS.cpp 
	src/file_reader.cpp
)
set_target_properties(msa-kmcl PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(msa-kmcl PRIVATE -pedantic-errors -Wall -Wextra -Wconversion -Wsign-conversion)

# Executable
add_executable(msa-kmc src/main.cpp)
target_compile_options(msa-kmc PRIVATE -pedantic-errors -Wall -Wextra -Wconversion -Wsign-conversion)
target_link_libraries(msa-kmc PRIVATE msa-kmcl) 

# Unit tests
#add_subdirectory(googletest)

# Python bindings
#add_subdirectory(pybind11)
#pybind11_add_module(msa_kmc binder.cpp)

