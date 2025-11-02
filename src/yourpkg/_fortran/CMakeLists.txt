cmake_minimum_required(VERSION 3.22)
project(yourpkg_fortran LANGUAGES C Fortran)

# Build shared library with Fortran kernels
add_library(yourpkg_fortran SHARED
    fmods/linalg.f90
)

set_target_properties(yourpkg_fortran PROPERTIES
    OUTPUT_NAME "yourpkg_fortran"
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)

# Install the shared library into the Python package folder inside the wheel
install(TARGETS yourpkg_fortran
    LIBRARY DESTINATION yourpkg
    RUNTIME DESTINATION yourpkg
    ARCHIVE DESTINATION yourpkg
)

