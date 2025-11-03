set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/xyz_bridge.f90"
  "${dir}/strucrd.f90"
  "${dir}/utilities.f90"
)

set(srcs ${srcs} PARENT_SCOPE)

