foreach(ITEM IN LISTS FILES_NEED_OPENACC_FLAGS)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -g -WF,-D_OPENMP -qsmp=omp -qoffload -qtbtable=full ")
endforeach()
