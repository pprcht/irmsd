cmake_minimum_required(VERSION 3.22)
project(irmsd_all LANGUAGES C Fortran)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        add_compile_options(
            -O0 -g
            -Wall -Wextra
            -fcheck=all
            -fbacktrace
            -finit-real=snan
        )
    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        add_compile_options(
            -O0 -g
            -warn all
            -check all
            -traceback
        )
    endif()
endif()

add_subdirectory(src/irmsd/_fortran)

