set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/strucrd.f90"
  "${dir}/utilities.f90"
  "${dir}/cn_exposed.f90"
  "${dir}/rotanisotrop.f90"
  "${dir}/term_ui.f90"
  "${dir}/axis_exposed.f90"
  "${dir}/canonical_exposed.f90"
  "${dir}/rmsd_exposed.f90"
  "${dir}/irmsd_exposed.f90"
  "${dir}/sorter_exposed.f90"
)

set(srcs ${srcs} PARENT_SCOPE)

