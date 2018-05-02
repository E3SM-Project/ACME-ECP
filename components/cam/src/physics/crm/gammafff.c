/*
 gamma-function for Fortran
 (C) Marat Khairoutdinov */

#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#pragma acc routine seq
float gammafff(float *x) {return (float)exp(lgamma(*x));}

#pragma acc routine seq
float gammafff_(float *x) {return (float)exp(lgamma(*x));}

#ifdef __cplusplus
}
#endif
