#include "tiff_image.h"
#include "tiff_utils.h"
#include <cassert>
#include <iostream>
#include "cell_utils.h"

#include <string.h>
#define DP(x) fprintf(stderr,"DEBUG %d\n", x)
#define DEBUG(var) printf(#var " = %d\n", var)

/*TiffImage::TiffImage(TIFF* tif) {

  // clear the raster data
  // note that TiffImage class is NOT responsible for managing TIFF memory
  clear_raster();

  // add the tiff
  __give_tiff(tif);
  
  }*/

TiffImage::TiffImage(const TiffReader& tr) {

  // make a local copy of this tiff reader
  m_tr = tr;

  m_ifd = tr.CurrentIFD();
  
  // Initial image properties come from the file
  m_height = m_ifd.height;
  m_width = m_ifd.width;
  m_photometric = m_ifd.photometric;
  m_planar = m_ifd.planar;
  m_samples_per_pixel = m_ifd.samples_per_pixel;
  m_bits_per_sample = m_ifd.bits_per_sample;

  m_pixels = static_cast<uint64_t>(m_width) * m_height;
  
}


//template <typename T>
//T TiffImage::pixel(uint64_t x, uint64_t y, int p) const {
uint8_t TiffImage::pixel(uint64_t x, uint64_t y, int p) const {

  assert(p == PIXEL_GRAY || p == PIXEL_RED || p == PIXEL_GREEN || p == PIXEL_BLUE || p == PIXEL_ALPHA);
  
  // right now only two types allowed
  //static_assert(std::is_same<T, uint32_t>::value || std::is_same<T, uint8_t>::value,
  //              "T must be either uint32_t or uint8_t");

  // check that the data has been read
  assert(__is_rasterized());
  
  // check that the pixel is in bounds
  uint64_t dpos = y * m_width + x; 
  if (dpos > m_pixels || x > m_width || y > m_height) {
    fprintf(stderr, "ERROR: Accesing out of bound pixel at (%llu,%llu)\n",x,y);
    assert(false);
  }
  
  if (p == PIXEL_GRAY) 
    return static_cast<uint8_t*>(m_data)[dpos];

  // else its rgb and so add the pixel offset
  return static_cast<uint8_t*>(m_data)[dpos*3 + p];
}

template <typename T>
T TiffImage::element(uint64_t e) const {

  // right now only two types allowed
  static_assert(std::is_same<T, uint32_t>::value || std::is_same<T, uint8_t>::value,
                "T must be either uint32_t or uint8_t");

  // check that the data has been read
  assert(__is_rasterized());
  
  // check that the element is in bounds
  if (e > m_pixels) {
    fprintf(stderr, "ERROR: Accesing out of bound pixel (flattened) at (%ul)\n",e);
    assert(false);
  }
    
  return static_cast<T*>(m_data)[e];
}

double TiffImage::mean() const {

  if (!__is_rasterized()) {
    std::cerr << "ERROR: needs to be read into raster first" << std::endl;
    return -1;
  }

  size_t mode = GetMode();
  
  double sum = 0;
  
  switch (mode) {
  case 8:
    for (size_t i = 0; i < m_pixels; ++i) {
      sum += static_cast<uint8_t*>(m_data)[i];
    }
    break;
  case 32:
    for (size_t i = 0; i < m_pixels; ++i)
      sum += static_cast<uint32_t*>(m_data)[i];
    break;
  }

  return (sum / m_pixels);

}

int TiffImage::Scale(double scale, bool ismode, int threads) {

  if (scale > 1) {
    fprintf(stderr, "Scale must be < 1\n");
    return 1;
  }

  //debug
  //  for (int x = 0; x < 100; x++)
  //  std::cerr << (unsigned int)static_cast<uint8_t*>(m_data)[x] << std::endl;
  
  size_t mode = GetMode();
  
  // New image as a row-major vector
  void* new_image;

  // New image dimensions
  uint64_t new_width = std::ceil(scale * m_width);
  uint64_t new_height = std::ceil(scale * m_height);
  uint64_t w_chunk = std::ceil(1/scale);
  uint64_t h_chunk = std::ceil(1/scale);

  // store the values you see
  std::vector<uint8_t*> sum3(3);
  uint8_t* sum8;
  uint64_t chunk_pixels = h_chunk * w_chunk;
  
  switch (mode) {
  case 8:
    new_image = calloc(new_width * new_height, sizeof(uint8_t));
    sum8 = (uint8_t*)calloc(chunk_pixels, sizeof(int));    
    break;
  case 3:
    new_image = calloc(new_width * new_height, sizeof(uint8_t) * 3);
    for (int x = 0; x < 3; x++)
      sum3[x] = (uint8_t*)calloc(chunk_pixels, sizeof(int));
    break;
  default:
    fprintf(stderr, "not able to understand mode %zu\n", mode);
    assert(false);
  }
  
  std::cerr << "new dim " << new_width << " x " << new_height <<
    " chunk " << w_chunk << " x " << h_chunk << std::endl;
  
  size_t debugint = 0;
  
  // Iterate over the new image
  #pragma omp parallel for num_threads(m_threads)
  for (uint64_t i = 0; i < new_height; i++) {
    for (uint64_t j = 0; j < new_width; j++) {
      ++debugint;
      if (debugint % 10000 == 0)
	std::cerr << " debugint " << AddCommas(debugint) << " of " << AddCommas(new_height * new_width)<< std::endl;
      size_t count = 0;
      
      // loop and do the averaging
      uint64_t ii_s = std::round(static_cast<double>(i) / scale);
      uint64_t jj_s = std::round(static_cast<double>(j) / scale);
      for (uint64_t ii = ii_s; ii < ii_s + h_chunk; ii++) {
	for (uint64_t jj = jj_s; jj < jj_s + w_chunk; jj++) {
	  //std::cout << "OLD: " << PAIRSTRING(ii,jj) << std::endl;
	  if (jj < m_width && ii < m_height) {
	    switch(mode) {
	    case 8:
	      sum8[count] = static_cast<uint8_t*>(m_data)[ii * m_width + jj];
	      break;
	    case 3:
	      for (int x = 0; x < 3; ++x) 
		sum3[x][count] = static_cast<uint8_t*>(m_data)[(ii * m_width + jj)*3 + x];
	      break;
	    }
	    ++count;
	  }
	}
      }
      
      funcmm_t mmfunc = ismode ? getMode : getMean;
      
      switch(mode) {
      case 8:
	static_cast<uint8_t*>(new_image)[i * new_width + j] = mmfunc(sum8, count);
	memset(sum8, 0, chunk_pixels * sizeof(uint8_t));
	break;
      case 3:
	for (int x = 0; x < 3; ++x) {
	  static_cast<uint8_t*>(new_image)[(i * new_width + j)*3 + x] = mmfunc(sum3[x], count);
	  memset(sum3[x], 0, chunk_pixels * sizeof(uint8_t));
	}
	break;
      }
      
    }// end new image j
  }// end new image i

  switch(mode) {
    case 8:
      free(sum8);
      break;
  case 3:
    for (int x = 0; x < 3; ++x)
      free(sum3[x]);
    break;
  }
  
  // transfer the data
  free(m_data);
  m_data = new_image;

  m_width = new_width;
  m_height = new_height;
  m_pixels = static_cast<uint64_t>(m_width) * m_height;
  std::cerr << "New " << PAIRSTRING(new_width, new_height) << " pixels " << m_pixels << std::endl;
  
  // set the new width and height
  /*if (!TIFFSetField(m_outtif, TIFFTAG_IMAGEWIDTH, new_width)) {
    fprintf(stderr, "ERROR: image (shrunk) failed to set width\n");
    m_width = 0;
    return 1;
  }
  
  if (!TIFFSetField(m_outtif, TIFFTAG_IMAGELENGTH, new_height)) {
    fprintf(stderr, "ERROR: image (shrunk) failed to set height\n");
    m_height = 0;
    return 1;
    }*/
  
  return 0;
  
}
int TiffImage::__alloc() {

  assert(m_pixels);
  
  size_t mode = GetMode(); 

  if (m_data) {
    std::cerr << "m_data not empty, freeing. Are you sure you want this?" << std::endl;
    free(m_data);
  }

  //debug
  //std::cerr << "m_pixels " << m_pixels << " mode " << mode << std::endl;
  
  switch (mode) {
      case 8:
	m_data = calloc(m_pixels, sizeof(uint8_t));
	//fprintf(stderr, "alloc size %f GB\n", m_pixels / 1e9);
	break;
       case 16:
         m_data = calloc(m_pixels, sizeof(uint16_t));
      case 32:
        m_data = calloc(m_pixels, sizeof(uint32_t));
	//fprintf(stderr, "alloc size %f GB\n", m_pixels * 4 / 1e9);
	break;
  case 3:
    m_data = calloc(m_pixels * 3, sizeof(uint8_t));
    break;
  default:
    fprintf(stderr, "not able to understand mode %zu\n", mode);
    assert(false);

    }
  
  if (m_data == NULL) {
    fprintf(stderr, "ERROR: unable to allocate image raster\n");
    return 1;
  }
  
  return 0;
}

int TiffImage::ReadToRaster() {

  m_data = m_ifd.ReadRaster();
  
  if (m_data == NULL)
    return 1;
  
  return 0;
}

void TiffImage::clear_raster() {

  if (m_data == NULL) {
    //fprintf(stderr, "Warning: no raster to clear\n");
    return;
  }

  _TIFFfree(m_data);

  m_data = NULL;

}

bool TiffImage::__is_rasterized() const {

  if (m_data == NULL) {
    fprintf(stderr, "Warning: Trying to write empty image. Try running ReadToRaster first\n");    
    return false;
  }
  return true;
}


uint8_t TiffImage::GetMode() const {

  // RGB
  if (m_samples_per_pixel == 3) 
    return 3;
  else if (m_bits_per_sample == 8) 
    return 8;
  else if (m_bits_per_sample == 16) 
    return 16;
  else if (m_bits_per_sample == 32)
    return 32;
  else if (m_samples_per_pixel == 4)
    return 4;
  return 0;
  
}

TiffImage::TiffImage(uint32_t width, uint32_t height, uint8_t mode) {

  // defaults
  m_photometric = PHOTOMETRIC_MINISBLACK;
  m_samples_per_pixel = 1;
  m_planar = PLANARCONFIG_CONTIG;

  switch(mode) {
  case 3:
    m_samples_per_pixel = 3;
    m_bits_per_sample = 8;
    m_photometric = PHOTOMETRIC_RGB;
    break;
  case 8:
    m_bits_per_sample = 8;
    break;
  case 16:
    m_bits_per_sample = 16;
    break;
  case 32:
    m_bits_per_sample = 32;
    break;
  default:
    fprintf(stderr, "Error: unknown mode: %d\n", mode);
    assert(false);
  }

  m_width = width;
  m_height = height;
  m_pixels = static_cast<uint64_t>(width) * height;

  __alloc();
}
