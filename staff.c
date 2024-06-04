#include "staff.h"
#include "main.h"

void saveToExcel(const char *filename, COMPLEX_NUMBER *x, const size_t N, const float rate) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return;
    }

    // Заголовки
    fprintf(file, "Frequency, Amplitude\n");

    // Вычисление частоты и амплитуды
    for (size_t i = 0; i < N/2; i++) {
        float frequency = i * rate / N; 
        float amplitude = 2.0 * sqrt(x[i].real * x[i].real + x[i].imag * x[i].imag) / N; 
        fprintf(file, "%f  %f\n", frequency, amplitude);
        printf("%f\n",amplitude);
    }

    fclose(file);
}

// глобалы на индексы частоты
size_t index_0 = 0;
size_t index_3600 = 0;
size_t index_7200 = 0;
size_t index_10800 = 0;
size_t index_14400 = 0;
size_t index_18000 = 0;
size_t index_22050 = 0;

// вычисление индексов если нужно
void chaaf(COMPLEX_NUMBER *x, size_t N, float rate)
{
    index_0 = (size_t)(0.0f / rate * N);
    index_3600 = (size_t)(3600.0f / rate * N);
    index_7200 = (size_t)(7200.0f / rate * N);
    index_10800 = (size_t)(10800.0f / rate * N);
    index_14400 = (size_t)(14400.0f / rate * N);
    index_18000 = (size_t)(18000.0f / rate * N);
    index_22050 = (size_t)(22050.0f / rate * N);

    if (sliders[0].value > 0) mltp(x, N, sliders[0].value, index_0, index_3600);
    if (sliders[0].value < 0) s_mltp(x, N, sliders[0].value, index_0, index_3600);

    if (sliders[1].value > 0) mltp(x, N, sliders[1].value, index_3600, index_7200);
    if (sliders[1].value < 0) s_mltp(x, N, sliders[1].value, index_3600, index_7200);

    if (sliders[2].value > 0) mltp(x, N, sliders[2].value, index_7200, index_10800);
    if (sliders[2].value < 0) s_mltp(x, N, sliders[2].value, index_7200, index_10800);

    if (sliders[3].value > 0) mltp(x, N, sliders[3].value, index_10800, index_14400);
    if (sliders[3].value < 0) s_mltp(x, N, sliders[3].value, index_10800, index_14400);

    if (sliders[4].value > 0) mltp(x, N, sliders[4].value, index_14400, index_18000);
    if (sliders[4].value < 0) s_mltp(x, N, sliders[4].value, index_14400, index_18000);

    if (sliders[5].value > 0) mltp(x, N, sliders[5].value, index_18000, index_22050);
    if (sliders[5].value < 0) s_mltp(x, N, sliders[5].value, index_18000, index_22050);
}

void mltp(COMPLEX_NUMBER *x, size_t N, float q, size_t index1, size_t index2)
{   
    // Компоненты спектра в диапазоне до index Гц
    for (size_t i = index1; i <= index2; i++) {
        x[i].real *= (float)q;
        x[i].imag *= (float)q;
    }
    // Симметричное для отрицательных частот
    size_t mirror_low_index = N - index1;
    size_t mirror_high_index = N - index2;

    for (size_t i = mirror_low_index; i <= mirror_high_index; i++) {
        x[i].real *= (float)q;
        x[i].imag *= (float)q;
    }
}

void s_mltp(COMPLEX_NUMBER *x, size_t N, float q, size_t index1, size_t index2)
{
    // Компоненты спектра в диапазоне до index Гц
    for (size_t i = index1; i <= index2; i++) {
        x[i].real /= (float)q;
        x[i].imag /= (float)q;
    }
    // Симметричное для отрицательных частот
    size_t mirror_low_index = N - index1;
    size_t mirror_high_index = N - index2;

    for (size_t i = mirror_low_index; i <= mirror_high_index; i++) {
        x[i].real /= (float)q;
        x[i].imag /= (float)q;
    }
}
