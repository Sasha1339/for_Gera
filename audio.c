#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "audio.h"
#include "staff.h"
#include "main.h"

int play(const char *file_name) {
    g_print("Processing file: %s\n", file_name);
    char filename[256];
    sprintf(filename, "%s", file_name);

    mpg123_handle *mh;
    unsigned char *buffer;
    size_t buffer_size;
    size_t done;
    int err;

    int driver;
    ao_device *dev;

    ao_sample_format format;
    int channels, encoding;
    long rate;

    // иниты
    ao_initialize();
    driver = ao_default_driver_id();
    mpg123_init();
    mh = mpg123_new(NULL, &err);
    buffer = (unsigned char*) malloc(SAMPLE_NODES * sizeof(unsigned char));

    // читаем параметры 
    mpg123_open(mh, filename);
    mpg123_getformat(mh, &rate, &channels, &encoding);

    // формат вывод и девайс
    format.bits = mpg123_encsize(encoding) * BITS;
    format.rate = rate;
    format.channels = channels;
    format.byte_format = AO_FMT_NATIVE;
    format.matrix = 0;
    dev = ao_open_live(driver, &format, NULL);

    // открыть файл на запись
    #if OUTPUT == 1
        SF_INFO sfinfo;
        memset(&sfinfo, 0, sizeof(sfinfo));
        sfinfo.samplerate = rate;
        sfinfo.channels = channels;
        sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16; // Output file format (WAV)
        SNDFILE *outfile = sf_open("output.wav", SFM_WRITE, &sfinfo);
    #endif

    // дэкод и воспроизведение

    while ((mpg123_read(mh, buffer, SAMPLE_NODES, &done) == MPG123_OK) && done > 0) {

        // Нормализация и преобразование в действительные числа для FFT
        short *samples = (short*) buffer;
        size_t sample_count = done / sizeof(short);
        
        // Нормализация 16-битных данных
        for (size_t i = 0; i < done; i++) {
            x[i].real = samples[i] / 32768.0f;  
            x[i].imag = 0.0f;
        }

        fourier_transform_real(x, done);
        
        #if DEBUG == 1
            saveToExcel("test_csv/fft_out.csv", x, done, rate);
        #endif

        chaaf(x, done, rate);
        
        #if DEBUG == 1
            saveToExcel("test_csv/chaaf_out.csv", x, done, rate);
        #endif

        inverted_fourier_transform_real(x, done);
     
        for (size_t i = 0; i < done; i++) {
            float normalized_value = x[i].real;
            if (normalized_value > 1.0f) {
                normalized_value = 1.0f;
            } else if (normalized_value < -1.0f) {
                normalized_value = -1.0f;
            }
            // Обратная нормализация в 16-битные данные
            samples[i] = (short)(normalized_value * 32767.0f ); 
        }

        #if OUTPUT == 1
            sf_write_raw(outfile, (char*)samples, done);
        #endif
        
        #if PLAY == 1
            ao_play(dev, (char*)samples, done);
        #endif
    }

    #if OUTPUT == 1
        sf_close(outfile);
    #endif

    // чистим 
    free(buffer);
    ao_close(dev);
    mpg123_close(mh);
    mpg123_delete(mh);
    mpg123_exit();
    ao_shutdown();

    return 0;
}
