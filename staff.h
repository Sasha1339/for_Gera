#ifndef _STAFF_H_
#define _STAFF_H_
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "z_fft.h" 
#include "Config.h" 

void saveToExcel(const char *filename, COMPLEX_NUMBER *x, const size_t N, const float rate);

void chaaf(COMPLEX_NUMBER *x, size_t N, float rate);
void mltp(COMPLEX_NUMBER *x, size_t N, float q, size_t index1, size_t index2);
void s_mltp(COMPLEX_NUMBER *x, size_t N, float q, size_t index1, size_t index2);

#endif /* _STAFF_H_ */