set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/linalg.f90"
  "${dir}/xyz_bridge.f90"
  "${dir}/strucrd.f90"
)

set(srcs ${srcs} PARENT_SCOPE)

