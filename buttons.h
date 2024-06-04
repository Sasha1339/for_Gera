#ifndef BUTTONS_H
#define BUTTONS_H

#include <gtk/gtk.h>

void button_array(GtkWidget *hbox);

void load_button(GtkWidget *widget, gpointer data);
void play_button(GtkWidget *widget, gpointer data);

#endif /* BUTTONS_H */