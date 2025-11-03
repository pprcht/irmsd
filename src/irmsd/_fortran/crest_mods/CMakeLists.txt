set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/atmasses.f90"
  "${dir}/axis_module.f90" 
)

set(srcs ${srcs} PARENT_SCOPE)

