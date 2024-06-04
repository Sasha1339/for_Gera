#include "buttons.h"
#include "audio.h"

GtkWidget *b_load;
GtkWidget *b_play;

// Глобал на путь
char *file_name = NULL;

// Загрузка
void load_button(GtkWidget *widget, gpointer data) {
    g_print("Load!\n");

    GtkWidget *dialog;
    GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
    gint res;

    // Диалог средатсвами ГТК
    dialog = gtk_file_chooser_dialog_new("Open File",
                                         GTK_WINDOW(data),
                                         action,
                                         "_Cancel",
                                         GTK_RESPONSE_CANCEL,
                                         "_Open",
                                         GTK_RESPONSE_ACCEPT,
                                         NULL);

    res = gtk_dialog_run(GTK_DIALOG(dialog));
    if (res == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
        file_name = gtk_file_chooser_get_filename(chooser);

        g_print("File selected: %s\n", file_name);
    }
    gtk_widget_destroy(dialog);
}

void play_button(GtkWidget *widget, gpointer data) {
    g_print("Play!\n");
    play(file_name);
}

void button_array(GtkWidget *hbox){

    GtkWidget *right_section = gtk_vbox_new(FALSE, 5);  
    gtk_container_set_border_width(GTK_CONTAINER(right_section), 20);

    gtk_box_pack_start(GTK_BOX(hbox), right_section, FALSE, FALSE, 0);

    b_load = gtk_button_new_with_label("Load");
    g_signal_connect(b_load, "clicked", G_CALLBACK(load_button), NULL);
    gtk_box_pack_start(GTK_BOX(right_section), b_load, TRUE, TRUE, 0);

    b_play = gtk_button_new_with_label("Play");
    g_signal_connect(b_play, "clicked", G_CALLBACK(play_button), NULL);
    gtk_box_pack_start(GTK_BOX(right_section), b_play, TRUE, TRUE, 0);
}