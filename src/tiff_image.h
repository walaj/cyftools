#ifndef TIFF_IMAGE_H
#define TIFF_IMAGE_H

#ifdef HAVE_TIFFLIB
#include <string>
#include <cassert>
#include <vector>
#include <tiffio.h>
#include <sstream>
#include <memory>

#include "tiff_cp.h"
#include "tiff_reader.h"
#include "tiff_utils.h"

#define PIXEL_GRAY 999
#define PIXEL_RED 0
#define PIXEL_GREEN 1
#define PIXEL_BLUE 2
#define PIXEL_ALPHA 3

/*
#include <random>
int getRandomInt() {
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    // Set the range of the generator to [0, 255]
    std:
    :uniform_int_distribution<> dis(0, 255);

    // Generate and return a random integer
    return dis(gen);
    }*/

// copy tags from one TIFF to another
static int TiffTagCopy(TIFF* o, TIFF* d);

class TiffImage {

 public:

  friend class TiffWriter;
  
  // create a new image from a TIFF file (libtiff)
  TiffImage(const TiffReader& tr);

  // create an empty iamge
  TiffImage() {}

  ~TiffImage() {
    clear_raster();
  }

  TiffImage(uint32_t width, uint32_t height, uint8_t mode);
  
  // access the data pointer
  void* data() const { return m_data; }

  // read the image tif to a 2D raster
  int ReadToRaster();

  // access a pixel as type T (uint8_t or uint32_t)
  //template <typename T>  
  uint8_t pixel(uint64_t x, uint64_t y, int p) const;

  void SetPixel(uint64_t x, uint64_t y, int p, uint8_t v);
  
  // access a pixel as type T (uint8_t or uint32_t) from the flattened index
  template <typename T>  
  T element(uint64_t e) const;

  // write the raster to the image using an already opened TIFF
  int write(TIFF* otif) const;

  // clear the raster and free the memory
  void clear_raster();

  // mean pixel value
  double mean() const;

  /// return the total number of pixels
  uint64_t numPixels() const { return m_pixels; }

  // set the verbose output flag
  void setverbose(bool v) { verbose = v; }

  // scale the image downx
  int Scale(double scale, bool ismode, int threads);

  // get the mode (gray 8-bit, RBG etc)
  uint8_t GetMode() const;

  uint32_t width() const { return m_width; }
  uint32_t height() const { return m_height; }   

  void setthreads(size_t threads) { m_threads = threads; }

  // template functions
  template <typename T>
  void SetPixel(uint64_t x, uint64_t y, int p, T v) {

  assert(p == PIXEL_GRAY || p == PIXEL_RED || p == PIXEL_GREEN || p == PIXEL_BLUE || p == PIXEL_ALPHA);

  // check that the pixel is in bounds
  uint64_t dpos = y * m_width + x; 
  if (dpos > m_pixels || x > m_width || y > m_height) {
    fprintf(stderr, "ERROR: Accesing out of bound pixel at (%llu,%llu)\n",x,y);
    assert(false);
  }

  int mode = GetMode();
  
  if (p == PIXEL_GRAY) {
    assert(mode == 8 || mode == 16);
    static_cast<T*>(m_data)[dpos] = v;
  } else {
    assert(mode == 3);
    // else its rgb and so add the pixel offset
    static_cast<uint8_t*>(m_data)[dpos*3 + p] = v;
  }
  
  return;
}

  
  
 private:

  bool verbose= false;

  size_t m_threads = 1;
  
  TiffReader m_tr;

  TiffIFD m_ifd;
  
  // raster of the image
  void* m_data = NULL; 

  // flag for whether data is stored as 1 byte or 4 bytes
  uint32_t m_width = 0, m_height = 0;

  // total number of pixels in the image
  uint64_t m_pixels = 0;

  // various TIFF flags
  // useful to know that options for photometric are:
  // - PHOTOMETRIC_MINISBLACK
  // - PHOTOMETRIC_MINISWHITE
  // - PHOTOMETRIC_RGB
  // - PHOTOMETRIC_MASK
  // - PHOTOMETRIC_DEPTH  
  // - PHOTOMETRIC_PALATTE
  // useful to know that options for planar config are:
  // PLANARCONFIG_CONTIG 1 -- RGBA RGBA RGBA etc (interleaved)
  // PLANARCONFIG_SEPARATE 2 -- RRRRRR GGGGGG BBBBB etc (separate planes)
  uint16_t m_photometric, m_planar, m_bits_per_sample, m_samples_per_pixel;

  // allocate the memory for the image raster
  // used for starting with a blank image
  int __alloc();
  
  // has the image had its raster stored to memory
  bool __is_rasterized() const;

};

#endif
#endif
