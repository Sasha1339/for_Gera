#ifndef AUDIO_H
#define AUDIO_H

#include <gtk/gtk.h>
#include <mpg123.h>
#include <ao/ao.h>
#include <sndfile.h>
#include "z_fft.h" 
#include "Config.h" 


#define  SAMPLE_NODES              (32768*2)  
COMPLEX x[SAMPLE_NODES];
#define BITS 8

int play(const char *file_name);

#endif /* AUDIO_H */