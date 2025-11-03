set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/adjacency.f90"
  "${dir}/atmasses.f90"
  "${dir}/axis_module.f90" 
  "${dir}/canonical.f90"
  "${dir}/cn.f90"
  "${dir}/crest_pars.f90"
  "${dir}/dijkstra.f90"
  "${dir}/floydwarshall.f90"
  "${dir}/geo.f90"
  "${dir}/hungarian.f90"
  "${dir}/irmsd_module.f90"
  "${dir}/ls_rmsd.f90"
  "${dir}/miscdata.f90"
)

set(srcs ${srcs} PARENT_SCOPE)

