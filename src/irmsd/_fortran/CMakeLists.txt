cmake_minimum_required(VERSION 3.22)
project(irmsd_fortran LANGUAGES C Fortran)

# Add srcs
set(srcs)
add_subdirectory("fmods")
add_subdirectory("crest_mods")

# Build shared library with Fortran kernels
add_library(irmsd_fortran SHARED
    ${srcs}
)

set_target_properties(irmsd_fortran PROPERTIES
    OUTPUT_NAME "irmsd_fortran"
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)

# Install the shared library into the Python package folder inside the wheel
install(TARGETS irmsd_fortran
    LIBRARY DESTINATION irmsd
    RUNTIME DESTINATION irmsd
    ARCHIVE DESTINATION irmsd
)

