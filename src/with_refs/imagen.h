#ifndef IMAGEN_H
#define IMAGEN_H
#include <stdio.h>

/**
 * Codigos de error para biblioteca de entrada/salida imagenes
 */
typedef enum codigo_error {
    PNM_OK=0,
    PNM_ERROR_LECTURA=1,
    PNM_ARCHIVO_INEXISTENTE=2,
    PNM_ENCABEZADO_INVALIDO=3,
    PNM_DATOS_INVALIDOS=4,
    PNM_ERROR_ESCRITURA=5,
    PNM_FORMATO_INVALIDO=6
} CodigoError;


const char* mensaje_error(CodigoError c);

/**
 * Tipos de imagen representadas
 */
typedef enum tipo_imagen {
    BIN_ASC=1,
    GRISES_ASC=2,
    BIN_BIN=4,
    GRISES_BIN=5,
    COLOR_ASC=3,
    COLOR_BIN=6
} TipoImagen;

/**
 * Canales de color
 */
typedef enum canal {
    ROJO=0,
    VERDE=1,
    AZUL=2
} Canal;

/**
 * pixel es un entero sin signo de al menos 32 bits
 */
typedef unsigned int Pixel;

typedef struct imagen {
    TipoImagen tipo;
    int ancho;
    int alto;
    int valor_maximo;
    Pixel* pixels;
    int npixels;
} Imagen;

/**
 * Dados un ancho y un alto, y una estructura de Imagen, reserva espacio en memoria para ancho*alto pixels en la imagen.
 */
int inicializar_imagen(int ancho, int alto, TipoImagen tipo, Imagen* pimg);

/**
 * Crea una imagen de mismos atributos que la original
 */
int duplicar_imagen(const Imagen* origen, Imagen* destino);

/**
 * Libera la memoria asociada a los pixels de la imagen dada.
 */
void destruir_imagen(Imagen* pimg);

CodigoError leer_encabezado(FILE* fimg, Imagen* pimg);
CodigoError leer_datos(FILE* fimg,Imagen* pimg);

/**
 * Carga una imagen PNM desde  un archivo dado.
 */
CodigoError leer_imagen(const char* ruta, Imagen* pimg);

/**
 * Escribe una imagen PNM desde en un archivo en la ruta especificada.
 */
CodigoError escribir_imagen(const Imagen* pimg, const char* ruta);

#endif
