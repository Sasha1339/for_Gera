#include "z_fft.h"
#include "z_math.h"

int fourier_transform(COMPLEX_NUMBER *buffer, uint32_t size)
{
    uint32_t length = size >> 1;
    uint32_t index_next = length;
    uint32_t index_last = size - 2;

    int index_previous,shift;
	static uint32_t M = 0;
	static int length_changed_step,length_changed_step_two;
	static float sR,sI;

    for (index_previous=1; index_previous <= index_last; index_previous++) {
        if (index_previous < index_next) {
            float tR = buffer[index_next].real;
            float tI = buffer[index_next].imag;
            buffer[index_next].real = buffer[index_previous].real;
            buffer[index_next].imag = buffer[index_previous].imag;
            buffer[index_previous].real = tR;
            buffer[index_previous].imag = tI;
		}
        shift = length;
		while (shift <= index_next) {
            index_next = index_next - shift;
            shift = shift >> 1;
		}
        index_next = index_next + shift;
	}

	for (int length_changed=1; length_changed <= shift_log(size); length_changed++) {
        float tR, tI;
        length_changed_step  = (int)(1 << length_changed);
        length_changed_step_two = (int)(length_changed_step >> 1);
        float uR = 1;
        float uI = 0;

        shift = shift_log(length_changed_step_two);
        sR = get_cos(shift);
        sI = -get_sin(shift);
		for (index_next=1; index_next <= length_changed_step_two; index_next++) {
			for (index_previous = index_next - 1; index_previous < size; index_previous += length_changed_step) {
				index_last = index_previous + length_changed_step_two;
                tR = buffer[index_last].real * uR - buffer[index_last].imag * uI;
                tI = buffer[index_last].real * uI + buffer[index_last].imag * uR;
                buffer[index_last].real = buffer[index_previous].real - tR;
                buffer[index_last].imag = buffer[index_previous].imag - tI;
                buffer[index_previous].real += tR;
                buffer[index_previous].imag += tI;
			}
			tR = uR;
			uR = tR * sR - uI * sI;
			uI = tR * sI + uI *sR;
		}
	}

	return 0;
}

int fourier_transform_real(COMPLEX_NUMBER *buffer, uint32_t size)
{
	int index_prev,index_next,len,shift;
	static uint32_t M = 0;
    static uint32_t ND4 = 0;
	static float sR,sI,tR,tI,uR,uI;

    M = size >> 1;
    for (index_prev=0; index_prev < M; index_prev++) {
        buffer[index_prev].real = buffer[index_prev << 1].real;
        buffer[index_prev].imag = buffer[(index_prev << 1) + 1].real;
    }

    fourier_transform(buffer, M);

    ND4 = size >> 2;
    for (index_prev=1; index_prev < ND4; index_prev++) {
        index_next = M - index_prev;
        shift = index_prev + M;
        len = index_next + M;
        buffer[shift].real = (buffer[index_prev].imag + buffer[index_next].imag) / 2;
        buffer[len].real = buffer[shift].real;
        buffer[shift].imag = -(buffer[index_prev].real - buffer[index_next].real) / 2;
        buffer[len].imag = -buffer[shift].imag;
        buffer[index_prev].real = (buffer[index_prev].real + buffer[index_next].real) / 2;
        buffer[index_next].real = buffer[index_prev].real;
        buffer[index_prev].imag = (buffer[index_prev].imag - buffer[index_next].imag) / 2;
        buffer[index_next].imag = -buffer[index_prev].imag;
    }
    buffer[size - ND4].real = buffer[ND4].imag;
    buffer[M].real = buffer[0].imag;
    buffer[size - ND4].imag = 0;
    buffer[M].imag = 0;
    buffer[ND4].imag = 0;
    buffer[0].imag = 0;

    uR = 1;
    uI = 0;
    shift = shift_log(M);
    sR = get_cos(shift);
    sI = -get_sin(shift);

    for (index_prev=0; index_prev < M; index_prev++) {
        shift = index_prev + M;
        tR = buffer[shift].real * uR - buffer[shift].imag * uI;
        tI = buffer[shift].real * uI + buffer[shift].imag * uR;
        buffer[shift].real = buffer[index_prev].real - tR;
        buffer[shift].imag = buffer[index_prev].imag - tI;
        buffer[index_prev].real += tR;
        buffer[index_prev].imag += tI;

        tR = uR;
        uR = tR * sR - uI * sI;
        uI = tR * sI + uI *sR;
    }
	return 0;
}

int inverted_fourier_transform(COMPLEX_NUMBER *buffer, uint32_t size)
{
	for (int k=0; k <= size - 1; k++) {
        buffer[k].imag = -buffer[k].imag;
	}

    fourier_transform(buffer, size);    /* using FFT */

	for (int k=0; k <= size - 1; k++) {
        buffer[k].real = buffer[k].real / size;
        buffer[k].imag = -buffer[k].imag / size;
	}

	return 0;
}


int inverted_fourier_transform_real(COMPLEX_NUMBER *buffer, uint32_t size)
{
	for (int k= (size >> 1) + 1; k < size; k++) {
        buffer[k].real = buffer[size - k].real;
        buffer[k].imag = -buffer[size - k].imag;
	}

    for (int k=0; k < size; k++) {
        buffer[k].real += buffer[k].imag;
    }

    fourier_transform_real(buffer, size);

	for (int k=0; k < size; k++) {
        buffer[k].real = (buffer[k].real + buffer[k].imag) / size;
        buffer[k].imag = 0;
	}

	return 0;
}

float get_sin(uint16_t factor){
    return sinf(PI / powf(2, factor));
}

float get_cos(uint16_t factor){
    return cosf(PI / powf(2, factor));
}