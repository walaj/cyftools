#ifndef TIFF_IFD_H
#define TIFF_IFD_H

#ifdef HAVE_TIFFLIB
#include <vector>
#include <iostream>
#include <tiffio.h>

// this always belongs as a member of the m_ifds vector
// in TiffReader or as a member of another TiffIFD
class TiffIFD {

 public:

  TiffIFD() {}

  TiffIFD(TIFF* tif);

  // this directoy id
  uint16_t dir = 0;

  // directory number of the SUB ifd
  uint16_t curr_ifd = 0;
  
  // required fields
  uint64_t height = 0;
  uint64_t width  = 0;  
  uint64_t samples_per_pixel = 0;
  uint64_t bits_per_sample  = 0;  

  uint16_t photometric = 0;
  uint16_t planar = 0;
  
  // optional fields
  uint64_t tile_height = 0;
  uint64_t tile_width  = 0;  
  uint64_t sample_format = 0;

  //friend std::ostream& operator<<(std::ostream& out, const TiffIFD& o);

  // print libtiff style info to stdout
  void print(); 

  bool isTiled();

  uint8_t GetMode() const;

  std::vector<double> mean();

  void* ReadRaster();
  
 private:

  TIFF* m_tif = NULL;

  void* __tiled_read();

  void* __lined_read();

  // allocate the memory for the raster
  // this is passed to the reader method, which then
  // passes it out. This class does not store raster info
  void* __alloc(); 
  
  std::vector<TiffIFD> m_subifds;

  // get the tag and if not found, print an error
  template <typename T>
  void __get_sure_tag(int tag, T& value);

  // get the tag and if not found, print an error
  template <typename T>
  void __get_tag(int tag, T& value);
  
};

#endif
#endif
