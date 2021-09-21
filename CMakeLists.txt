cmake_minimum_required(VERSION 3.10.2)
project(msa_kmc)

include_directories(include)
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

add_executable(msa-kmc ${SOURCES})
