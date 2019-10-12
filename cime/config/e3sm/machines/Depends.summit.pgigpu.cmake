
foreach(ITEM IN LISTS FILES_NEED_OPENACC_FLAGS)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -Minline -ta=nvidia,cc70,fastmath,loadcache:L1,unroll,fma,managed,ptxinfo -Mcuda -Minfo=accel ")
endforeach()

foreach(ITEM IN LISTS FILES_NEED_CUDA_C_FLAGS)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -O3 -std=c++11 ")
endforeach()

