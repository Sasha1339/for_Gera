#ifndef _Z_FFT_H
#define _Z_FFT_H

#include "Config.h"

#define TYPE_FFT_E     float    /* Type is the same with COMPLEX member */     

#ifndef PI
#define PI             (3.14159265f)
#endif

typedef COMPLEX COMPLEX_NUMBER;  /* Define COMPLEX in Config.h */

extern int fourier_transform(COMPLEX_NUMBER *buffer, uint32_t size);
extern int fourier_transform_real(COMPLEX_NUMBER *buffer, uint32_t size);
extern int inverted_fourier_transform(COMPLEX_NUMBER *buffer, uint32_t size);
extern int inverted_fourier_transform_real(COMPLEX_NUMBER *buffer, uint32_t size);
extern float get_sin(uint16_t factor);
extern float get_cos(uint16_t factor);

#endif /* _Z_FFT_H_ */