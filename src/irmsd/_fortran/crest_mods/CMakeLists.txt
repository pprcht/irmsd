set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/atmasses.f90"
  "${dir}/axis_module.f90" 
  "${dir}/cn.f90"
  "${dir}/crest_pars.f90"
  "${dir}/geo.f90"
  "${dir}/hungarian.f90"
  "${dir}/ls_rmsd.f90"
  "${dir}/miscdata.f90"
)

set(srcs ${srcs} PARENT_SCOPE)

