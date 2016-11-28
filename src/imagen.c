#include "imagen.h"
#include <stdlib.h>
#include <stdio.h>

const char* tabla_mensaje_error[] = { "Todo bien","Error de lectura", "Archivo inexistente", "Encabezado invalido", "Datos invalidos", "Error de escritura" };

const char* mensaje_error(CodigoError c) {
    return tabla_mensaje_error[c];
}
/*
 * TIEMPOS:
 *   IMPLEMENTAR imagen.h + imagen.c : 1:00
 *   DEPURAR imagen.h + imagen.c     : 0:40
 */
int inicializar_imagen(int ancho, int alto, TipoImagen tipo, Imagen* pimg) {
    if (alto*ancho > 0) {
        pimg->pixels = (Pixel*) malloc(ancho*alto*sizeof(Pixel));
        pimg->ancho = ancho;
        pimg->alto = alto;
        pimg->npixels = ancho*alto;
        pimg->tipo = tipo;
	pimg->valor_maximo = tipo == 4 ? 1 : 255;
        return (!pimg)? 1: 0;
    }
    return 1;
}

int duplicar_imagen(const Imagen* po, Imagen* pd) {
    pd->valor_maximo = po->valor_maximo;
    return inicializar_imagen(po->ancho,po->alto,po->tipo,pd);
}

void destruir_imagen(Imagen* pimg) {
    free(pimg->pixels);
    pimg->ancho = 0;
    pimg->alto = 0;
    pimg->pixels = NULL;
}

static CodigoError leer_datos_p2(FILE* fimg, Imagen* pimg);
static CodigoError leer_datos_p4(FILE* fimg, Imagen* pimg);
static CodigoError leer_datos_p5(FILE* fimg, Imagen* pimg);
static CodigoError leer_datos_p6(FILE* fimg, Imagen* pimg);

CodigoError leer_datos(FILE* fimg,Imagen* pimg) {
    switch (pimg->tipo) {
    case GRISES_ASC:
        return leer_datos_p2(fimg,pimg);
        break;
    case GRISES_BIN:
        return leer_datos_p5(fimg,pimg);
        break;
    case COLOR_BIN:
        return leer_datos_p6(fimg,pimg);
        break;
    case BIN_BIN:
        return leer_datos_p4(fimg,pimg);
        break;
    default:
        printf("Tipo de imagen %d no soportado.\n",pimg->tipo);
        return PNM_FORMATO_INVALIDO;
    }
}

CodigoError leer_imagen(const char* ruta, Imagen* pimg) {
    CodigoError ret;
    /* abrimos archivo */
    FILE* fimg = fopen(ruta, "r");
    if (fimg == NULL) {
        return PNM_ARCHIVO_INEXISTENTE;
    }
    /* leemos encabezado */
    ret  = leer_encabezado(fimg,pimg);
    if (ret != PNM_OK) {
        fclose(fimg);
        return ret;
    }
    /* leemos datos según tipo */
    inicializar_imagen(pimg->ancho,pimg->alto,pimg->tipo,pimg);
    leer_datos(fimg,pimg);
    if (ret != PNM_OK)
        destruir_imagen(pimg);
    fclose(fimg);
    return ret;
}

CodigoError leer_datos_p2(FILE* fimg, Imagen* pimg) {
    register int i,aux;
    const int n = pimg->ancho*pimg->alto;
    aux = 1;
    for (i = 0; (i < n) && (aux > 0); ++i) {
        aux = fscanf(fimg," %d",(pimg->pixels + i));
    }
    if (i < n)
        return PNM_DATOS_INVALIDOS;
    else
        return PNM_OK;
}

CodigoError leer_datos_p4(FILE* fimg, Imagen* pimg) {
  register int i,li;
    const int m = pimg->alto;
    const int n = pimg->ancho;
    for (i = 0, li=0; i < m; ++i) {
      unsigned short mask = 0;
      int j;
      short r = 0;
      for (j = 0; j < n; ++j, ++li) {
	if (!mask) {
	  r = fgetc(fimg); 
	  mask = 0x80;
	  if (r < 0) return PNM_DATOS_INVALIDOS;
	}
        pimg->pixels[li] = (r & mask) ? 1 : 0;
	mask >>= 1;
      }
    }
    if (li < n)
        return PNM_DATOS_INVALIDOS;
    else
        return PNM_OK;
}

CodigoError leer_datos_p6(FILE* fimg, Imagen* pimg) {
    register int i, aux;
    unsigned char r,g,b;
    const int n = pimg->ancho*pimg->alto;
    aux = 3;
    for (i = 0; (i < n) && (aux == 3); ++i) {
        aux = fread(&r,sizeof(char),1,fimg);    /* maximo asumido de 8 bits por canal! */
        aux += fread(&g,sizeof(char),1,fimg);
        aux += fread(&b,sizeof(char),1,fimg);
        pimg->pixels[i] = (r << 16) + (g << 8) + b;
    }
    if (i < n)
        return PNM_DATOS_INVALIDOS;
    else
        return PNM_OK;
}

CodigoError leer_datos_p5(FILE* fimg, Imagen* pimg) {
    register int i;
    unsigned char p;
    const int n = pimg->ancho*pimg->alto;
    for (i = 0; i < n; ++i) {
        int r = fread(&p,sizeof(char),1,fimg);    /* maximo asumido de 8 bits por canal! */
        if (!r) return PNM_DATOS_INVALIDOS;
        pimg->pixels[i] = p;
    }
    if (i < n)
        return PNM_DATOS_INVALIDOS;
    else
        return PNM_OK;
}

CodigoError leer_encabezado(FILE* fimg, Imagen* pimg) {
    char magic;
    int ancho,alto,maxval = 0, aux = 0;
    magic = fgetc(fimg);
    /* 1. numero magico */
    if (magic != 'P')
        return PNM_ENCABEZADO_INVALIDO;
    magic = fgetc(fimg);
    pimg->tipo = magic - '0';
    if (pimg->tipo >6)
        return PNM_ENCABEZADO_INVALIDO;

    /* ancho */
    aux = fscanf(fimg," %d",&ancho);
    if ((ancho == 0) || (aux == 0))
        return PNM_ENCABEZADO_INVALIDO;
    /* alto */
    aux = fscanf(fimg," %d",&alto);
    if ((alto == 0) || (aux == 0))
        return PNM_ENCABEZADO_INVALIDO;
    if (pimg->tipo != 4) {
    /* maxval */
    aux = fscanf(fimg," %d",&maxval);
    if ((maxval == 0) || (aux == 0))
        return PNM_ENCABEZADO_INVALIDO;
    if ((maxval > 255) && (magic == '6'))
        return PNM_ENCABEZADO_INVALIDO;
    } else {
      maxval = 1;
    }
    /* un espacio al final */
    aux = fgetc(fimg);
    pimg->ancho = ancho;
    pimg->alto = alto;
    pimg->valor_maximo = maxval;
    return PNM_OK;
}

CodigoError escribir_datos_p2(const Imagen* pimg, FILE* fimg);
CodigoError escribir_datos_p4(const Imagen* pimg, FILE* fimg);
CodigoError escribir_datos_p5(const Imagen* pimg, FILE* fimg);
CodigoError escribir_datos_p6(const Imagen* pimg, FILE* fimg);

/**
 * Escribe una imagen PNM desde en un archivo en la ruta especificada.
 */
CodigoError escribir_imagen(const Imagen* pimg, const char* ruta) {
    FILE* fimg;
    int ret;
    fimg = fopen(ruta,"w");
    if (fimg == NULL)
        return PNM_ARCHIVO_INEXISTENTE;

    /* encabezado */
    if (pimg->tipo != 4)
      fprintf(fimg,"P%c %d %d\n%d\n",pimg->tipo+'0', pimg->ancho, pimg->alto, pimg->valor_maximo);
    else 
      fprintf(fimg,"P%c %d %d\n",pimg->tipo+'0', pimg->ancho, pimg->alto);
    switch (pimg->tipo) {
    case GRISES_ASC:
        ret = escribir_datos_p2(pimg,fimg);
        break;
    case BIN_BIN:
        ret = escribir_datos_p4(pimg,fimg);
        break;
    case GRISES_BIN:
        ret = escribir_datos_p5(pimg,fimg);
        break;
    case COLOR_BIN:
        ret = escribir_datos_p6(pimg,fimg);
        break;
    default:
        printf("Tipo de imagen %d no soportado.\n",pimg->tipo);
        ret = PNM_FORMATO_INVALIDO;
    }
    fclose(fimg);
    return ret;
}


CodigoError escribir_datos_p2(const Imagen* pimg, FILE* fimg) {
    /* datos */
    const int ancho = pimg->ancho;
    const int alto = pimg->alto;
    register int i,j;
    Pixel* pix = pimg->pixels;
    for (i = 0; i < alto; ++i) {
        for (j = 0; j < ancho; ++j, ++pix) {
            fprintf(fimg," %d", *pix);
        }
        fputc('\n',fimg);
    }
    /* fin */
    return PNM_OK;
}

CodigoError escribir_datos_p4(const Imagen* pimg, FILE* fimg) {
    const int ancho = pimg->ancho;
    const int alto = pimg->alto;
    const Pixel* pixels = pimg->pixels;
    register int i,j;
    for (i = 0; i < alto; ++i) {
      unsigned short mask = 0x80;
      unsigned short word = 0x00;
      for (j = 0; j < ancho; ++j,++pixels) {
        unsigned short p = *pixels;
	if (p) 
	  word |= mask;
	mask >>=1;
	if (!mask) {
	  fputc(word, fimg);
	  mask = 0x80;
	  word = 0;
	}
      }
      if ((mask > 0) && (mask != 0x80))
          fputc(word, fimg);      
    }
    return PNM_OK;
}

CodigoError escribir_datos_p5(const Imagen* pimg, FILE* fimg) {
    const int ancho = pimg->ancho;
    const int alto = pimg->alto;
    const int n = ancho*alto;
    const Pixel* pixels = pimg->pixels;
    register int i;
    for (i = 0; i < n; ++i,++pixels) {
        int p = *pixels;
        fputc(p & 0xff, fimg);
    }
    return PNM_OK;
}

CodigoError escribir_datos_p6(const Imagen* pimg, FILE* fimg) {
    const int ancho = pimg->ancho;
    const int alto = pimg->alto;
    const int n = ancho*alto;
    const Pixel* pixels = pimg->pixels;
    register int i;
    for (i = 0; i < n; ++i,++pixels) {
        int p = *pixels;
        fputc((p >> 16) & 0xff, fimg); /* r */
        fputc((p  >> 8) & 0xff, fimg); /* g */
        fputc((p      ) & 0xff, fimg); /* b */
    }
    return PNM_OK;
}


