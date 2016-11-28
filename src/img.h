#ifndef IMG_H
#define IMG_H
#include <cstdlib>

typedef unsigned int pixel_t;  // NEEDED FOR COLOR EXAMPLES TO WORK!!!
//typedef unsigned char pixel_t; // THIS IS NEEDED IN ORDER TO USE IQA LIBRARY, COLOR BECOMES BROKEN!!
/** very very simple 8 bit grayscale  image type */
class Image {
 public:
  int m;
  int n;
  int npixels;
  int maxval;
  pixel_t* pixels;

  inline Image(): m(0), n(0), npixels(0), maxval(0), pixels(NULL) {}

  inline Image(int _m, int _n, int _mv): m(_m), n(_n), npixels(_m*_n), 
    maxval(_mv) {
    pixels = new pixel_t[npixels];
  }

  inline ~Image() { 
    if (pixels)
      delete[] pixels; 
    m = n = npixels = 0; 
  }

  inline pixel_t& operator[](int i) {
    return pixels[i];
  }

  inline const pixel_t& operator[](int i) const {
    return pixels[i];
  }
  
  inline pixel_t& operator()(int i, int j) {
    i = i < 0 ? 0 : ( i >= m ? m-1 : i );
    j = j < 0 ? 0 : ( j >= n ? n-1 : j );
    return pixels[i*n+j];
  }

  inline const pixel_t& operator()(int i, int j) const {
    i = i < 0 ? 0 : ( i >= m ? m-1 : i );
    j = j < 0 ? 0 : ( j >= n ? n-1 : j );
    return pixels[i*n+j];
  }

};

#endif
