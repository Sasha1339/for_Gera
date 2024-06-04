#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "buttons.h"
#include "audio.h"

#define DEBUG 0
#define OUTPUT 1
#define PLAY 1

GtkWidget *window;
GtkWidget *button;
GtkWidget *hbox;
GtkWidget *vbox;
GtkWidget *slider;
GtkWidget *progress_bar;

GtkWidget *vbox_main;

// Структура для хранения имени слайдера и его значения 
typedef struct {
    int value_name;
    int value;
} Slider;
Slider sliders[6];

void extractNumber(const char *str, int *num);
Slider* add_slider_value(const char* slider_name, int slider_value);

void w_slider(GtkWidget *widget, gpointer data);
void main_window(void);
void slider_array(void);

void default_value_sliders();
