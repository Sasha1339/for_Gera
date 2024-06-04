/*
Разработать программное обеспечение с графическим пользовательским интерфейсом на любом языке программирования. 
Входные тестовые данные для каждого варианта ДЗ должны загружаться из текстового, 
либо звукового файла (в зависимости от условия задания). 
Генерирование входных тестовых файлов допускается осуществлять с помощью скриптов Matlab. 
Использование библиотечных функций для реализации алгоритмов обработки 
сигналов не допускается (если в задании не указано иное).
*/

/*
Разработать шестиканальный цифровой эквалайзер звукового сигнала. 
Должна быть возможность регулировки параметров ослабления или усиления каждого из каналов. 
Ширина полосы пропускания всех каналов должна быть одинаковой. 
Для чтения и записи данных аудиофайлов допускается использовать сторонние библиотеки. 
Обработка звука в реальном времени и его воспроизведение не требуются. 
Для демонстрации работы эквалайзера допускается использовать готовую библиотеку БПФ.
*/

#include "main.h"

void default_value_sliders(){
    for (int i = 0; i < 6; i++)
        sliders[i].value_name = i;
}

void extractNumber(const char *str, int *num) {
    while (*str) 
    { 
        if (*str >= '0' && *str <= '9') 
        {
            *num = strtol(str, NULL, 10);
            break;
        }
        str++;
    }
}

// Обработчик значений слайдеров
Slider* add_slider_value(const char* slider_name, int slider_value) {
    int value_name;
    extractNumber(slider_name, &value_name);
    if (value_name < 7 && value_name > 0) 
    {   
        switch (value_name)
        {
        case 1:
            sliders[0].value_name = value_name;
            sliders[0].value = slider_value;
            break;
        case 2:
            sliders[1].value_name = value_name;
            sliders[1].value = slider_value;
            break;
        case 3:
            sliders[2].value_name = value_name;
            sliders[2].value = slider_value;
            break;
        case 4:
            sliders[3].value_name = value_name;
            sliders[3].value = slider_value;
            break;
        case 5:
            sliders[4].value_name = value_name;
            sliders[4].value = slider_value;
            break;
        case 6:
            sliders[5].value_name = value_name;
            sliders[5].value = slider_value;
            break;
        default:
            break;
        }
        
        #if DEBUG == 1
            for(int i = 0; i < 6; i++)
            printf("Слайдер %d - значение %d \n",sliders[i].value_name, sliders[i].value);
        #endif
    } 

    return sliders;
}


// Обработчик слайдеров
void w_slider(GtkWidget *widget, gpointer data) {
    gint value = gtk_range_get_value(GTK_RANGE(widget));
    gchar *slider_name = g_object_get_data(G_OBJECT(widget), "name");
    
    add_slider_value(slider_name, value);
    
    #if DEBUG == 2 
        g_print("%s slider value: %d\n", slider_name, value);
    #endif
}


// Окно
void main_window(void)
{
    window =  gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "Эквалайзер");
    gtk_window_set_default_size (GTK_WINDOW (window), 1100, 700);
    gtk_widget_show(window);
    g_signal_connect(window, "destroy", G_CALLBACK(gtk_main_quit), NULL);
}

void slider_array(void)
{
    const gchar *slider_names[] = {"Slider 1", "Slider 2", "Slider 3", "Slider 4", "Slider 5", "Slider 6"};
    const gchar *freq_array[] = {"0-3600", "3600-7200", "7200-10800", "10800-14400", "14400-18000", "18000-22050"};

    for (int i = 0; i < 6; i++) {
        
        vbox = gtk_vbox_new(FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);

        // Лэйблы
        gchar *label_text = g_strdup_printf("Slider %d: \n%s\n", i+1, freq_array[i]);
        GtkWidget *label = gtk_label_new(label_text);
        // printf("%s",label_text );
        g_free(label_text);
        gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);

        // Слайдеры
        slider = gtk_vscale_new_with_range(-5, 5, 1);
        gtk_range_set_inverted(slider, TRUE);
        gtk_range_set_value(slider, 0);
        gtk_box_pack_start(GTK_BOX(vbox), slider, TRUE, TRUE, 0);
        
        g_object_set_data(G_OBJECT(slider), "name", slider_names[i]);
        g_signal_connect(slider, "value-changed", G_CALLBACK(w_slider), label);

        default_value_sliders();
    }
}

int main(int argc, char **argv) {

    gtk_init(&argc, &argv);
    main_window();
    
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER(hbox), 50);
    gtk_container_add(GTK_CONTAINER(window), hbox);

    slider_array();
    button_array(hbox);
    gtk_widget_show_all(window);

    gtk_main();
    return 0;
}