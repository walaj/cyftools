#include "tiff_ifd.h"
#include "tiff_utils.h"
#include "cell_utils.h"
#include <cstring>
#include <cassert>

template <typename T>  
void TiffIFD::__get_sure_tag(int tag, T& value) {

  if (!TIFFGetField(m_tif, tag, &value)) {
    fprintf(stderr, "ERROR: image does not have tag specified: %d\n", tag);
    value = 0;
  }
  
}

template <typename T>  
void TiffIFD::__get_tag(int tag, T& value) {
  
  if (!TIFFGetField(m_tif, tag, &value)) {
    value = 0;
  }
  
}

/*std::ostream& operator<<(std::ostream& out, const TiffIFD& o) {
  out << "----- IFD " << o.m_dir << "-----" << std::endl;
  out << "Dim: " << o.width << " x " << o.height << std::endl;
  
  TIFFPrintDirectory(o.m_tif, stdout);
  
  return out;
  }*/

void TiffIFD::print() {

  // store it and then move it back. Kludgy
  uint16_t cur_dir = TIFFCurrentDirectory(m_tif);

  TIFFSetDirectory(m_tif, dir);

  TIFFPrintDirectory(m_tif, stdout);

  // and put it back
  TIFFSetDirectory(m_tif, cur_dir); 
  
}

TiffIFD::TiffIFD(TIFF* tif) {

  // what we are doing here is creating a copy of the TIFF
  // filepointer that refers to this particular IFD
  // It's basically the same as every other TIFF except that it
  // will contain the tif_diroff set to THIS directory
  // this does NOT copy the image itself (aka the pixels)
  //std::memcpy(&m_tif, tif, sizeof(*tif));

  m_tif = tif;

  dir = TIFFCurrentDirectory(tif);

  // store basic image properties that must be ther
  __get_sure_tag(TIFFTAG_IMAGEWIDTH, width);
  __get_sure_tag(TIFFTAG_IMAGELENGTH, height);
  //__get_sure_tag(TIFFTAG_BITSPERSAMPLE, bits_per_sample);
  //__get_sure_tag(TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);

  if (!TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel)) {
    fprintf(stderr, "WARNING: image does not have a SAMPLESPERPIXEL tag. Assuming 1\n");
    samples_per_pixel = 1;
  }
  
  if (!TIFFGetField(m_tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample)) {
    fprintf(stderr, "WARNING: image does not have a BITSPERSAMPLE tag. Assuming 8\n");
    bits_per_sample = 8;
  }
  
  if (!TIFFGetField(m_tif, TIFFTAG_PHOTOMETRIC, &photometric)) {
    fprintf(stderr, "WARNING: image does not have a PHOTOMETRIC tag. Assuming MINISBLACK\n");
    photometric = PHOTOMETRIC_MINISBLACK;
  }
  
  // note that optinos for planar config are:
  // PLANARCONFIG_CONTIG 1 -- RGBA RGBA RGBA etc (interleaved)
  // PLANARCONFIG_SEPARATE 2 -- RRRRRR GGGGGG BBBBB etc (separate planes)
  if (!TIFFGetField(m_tif, TIFFTAG_PLANARCONFIG, &planar)) {
    fprintf(stderr, "WARNING: image does not have a PLANARCONFIG tag. Assuming PLANARCONFIG_CONTIG\n");
    planar = PLANARCONFIG_CONTIG;
  }

  
  // store other image properties to keep them handye
  __get_tag(TIFFTAG_SAMPLEFORMAT, sample_format);
  __get_tag(TIFFTAG_TILEWIDTH, tile_width);
  __get_tag(TIFFTAG_TILELENGTH, tile_height);      

  /*
   * Check the image to see if TIFFReadRGBAImage can deal with it.
   * 1/0 is returned according to whether or not the image can
   * be handled.  If 0 is returned, emsg contains the reason
   * why it is being rejected.
   */
  char emsg[1024];
  int rgba_ok = TIFFRGBAImageOK(m_tif, emsg);
  if (!rgba_ok) {
    fprintf(stderr, "%s", emsg);
  }

  
}

bool TiffIFD::isTiled() {

  // store it and then move it back. Kludgy
  uint16_t tmp_dir = TIFFCurrentDirectory(m_tif);
  TIFFSetDirectory(m_tif, dir);
  
  bool out = TIFFIsTiled(m_tif);
  
  // and put it back
  TIFFSetDirectory(m_tif, tmp_dir); 

  return out;
  
}

std::vector<double> TiffIFD::mean() {

  // store it and then move it back. Kludgy
  uint16_t tmp_dir = TIFFCurrentDirectory(m_tif);
  TIFFSetDirectory(m_tif, dir);
  
  
  uint8_t mode = GetMode();
    
  void* buf;
  
  if (isTiled()) {
    // allocate memory for a single tile
    buf  = _TIFFmalloc(TIFFTileSize(m_tif));
  } else {
    buf  = _TIFFmalloc(TIFFScanlineSize(m_tif));
  }

  // loop through the tiles
  int x, y;
  double np = static_cast<uint64_t>(width) * height;

  std::vector<double> out = {0,0,0,0}; 
    
  /*  if (verbose)
      std::cerr << "...starting image " <<
      ifd << " with " << AddCommas(static_cast<uint64_t>(np)) <<
      " pixels" << std::endl;
      
      uint64_t m_pix = 0; // running count for verbose command
  */
  
  
  if (isTiled()) {
    for (y = 0; y < height; y += tile_height) {
      for (x = 0; x < width; x += tile_width) {
	
	// Read the tile
	if (TIFFReadTile(m_tif, buf, x, y, 0, 0) < 0) {
	  fprintf(stderr, "Error reading tile at (%d, %d)\n", x, y);
	  return {-1};
	}
	
	// Get the mean
	for (int ty = 0; ty < tile_height; ty++) {
	  for (int tx = 0; tx < tile_width; tx++) {
	    if (x + tx < width && y + ty < height) {
	      
	      /*		  ++m_pix;
		  if (verbose && m_pix % static_cast<uint64_t>(1e9) == 0)
		    std::cerr << "...working on pixel: " <<
		      AddCommas(static_cast<uint64_t>(m_pix)) << std::endl;
	      */  
	      switch (mode) {
	      case 8:
		out[3] += static_cast<uint8_t*>(buf)[ty * tile_width + tx];
		break;
	      case 3:
		out[0] += static_cast<uint8_t*>(buf)[(ty * tile_width + tx) * 3    ];
		out[1] += static_cast<uint8_t*>(buf)[(ty * tile_width + tx) * 3 + 1];
		out[2] += static_cast<uint8_t*>(buf)[(ty * tile_width + tx) * 3 + 2];
		break;
	      case 32:
		out[3] += static_cast<uint32_t*>(buf)[ty * tile_width + tx];	      
		break;
	      }
		}
	      } // tile x loop
	    } // tile y loop
	  } // image x loop
	} // image y loop
      }
      
      // lined image
      else {
	
	uint64_t ls = TIFFScanlineSize(m_tif);
	for (uint64_t y = 0; y < height; y++) {
	  for (uint64_t x = 0; x < ls; x++) {

	    // Read the line
	    if (TIFFReadScanline(m_tif, buf, y) < 0) {
	      fprintf(stderr, "Error reading line at row %llu\n", y);
	      return {-1};
	    }
	    
	    /*	    ++m_pix;
	    if (verbose && m_pix % static_cast<uint64_t>(1e9) == 0)
	      std::cerr << "...working on pixel: " <<
		AddCommas(static_cast<uint64_t>(m_pix)) << std::endl;
	    */
	    
	    switch (mode) {
	    case 8:
	      out[3] += static_cast<uint8_t*>(buf)[x];
	    break;
	    case 3:
	      out[0] += static_cast<uint8_t*>(buf)[x * 3    ];
	      out[1] += static_cast<uint8_t*>(buf)[x * 3 + 1];
	      out[2] += static_cast<uint8_t*>(buf)[x * 3 + 2];
	      break;
	    case 4:
	      out[3] += static_cast<uint32_t*>(buf)[x];	      
	      break;
	    }
	  } // x loop
	} // y loop
      } //scanline else

  // get the mean
  for (auto& a : out)
    a /= np;

  // and put it back
  TIFFSetDirectory(m_tif, tmp_dir); 
  
  return out;
}

uint8_t TiffIFD::GetMode() const {

  // RGB
  if (samples_per_pixel == 3) 
    return 3;
  else if (bits_per_sample == 8) 
    return 8;
  else if (bits_per_sample == 16) 
    return 16;
  else if (bits_per_sample == 32)
    return 32;
  else if (samples_per_pixel == 4)
    return 4;
  return 0;
  
}

void* TiffIFD::ReadRaster() {

  if (isTiled())
    return(__tiled_read());
  else
    return(__lined_read());
  
}


void* TiffIFD::__alloc() {
  
  size_t mode = GetMode(); 

  uint64_t pixels = static_cast<uint64_t>(width) * height;
  
  void* data; 
  
  switch (mode) {
      case 8:
	data = calloc(pixels, sizeof(uint8_t));
	//fprintf(stderr, "alloc size %f GB\n", m_pixels / 1e9);
	break;
      case 4:
        data = calloc(pixels, sizeof(uint32_t));
	//fprintf(stderr, "alloc size %f GB\n", m_pixels * 4 / 1e9);
	break;
      case 3:
        data = calloc(pixels*3, sizeof(uint8_t));
	//fprintf(stderr, "alloc size %f GB\n", m_pixels * 4 / 1e9);
	break;
  default:
    fprintf(stderr, "mode not recognized: %zu\n", mode);
    assert(false);
    }
  
  if (data == NULL) {
    fprintf(stderr, "ERROR: unable to allocate image raster\n");
    return NULL;
  }
  
  return data;
}


void* TiffIFD::__tiled_read() {

  // store it and then move it back. Kludgy
  uint16_t tmp_dir = TIFFCurrentDirectory(m_tif);
  TIFFSetDirectory(m_tif, dir);

  uint8_t mode = GetMode();
  
  // allocate memory for a single tile
  void* tile = (void*)calloc(TIFFTileSize(m_tif), sizeof(uint8_t));

  // allocate the memory for this buffer
  // IFD object is NOT in charge of storing this
  void* data = __alloc();

  uint64_t tile_size = static_cast<uint64_t>(tile_height) * tile_width;						 

  // loop through the tiles
  uint64_t x, y;
  uint64_t m_pix = 0;
  for (y = 0; y < height; y += tile_height) {
    for (x = 0; x < width; x += tile_width) {
      
      // Read the tile
      if (TIFFReadTile(m_tif, tile, x, y, 0, 0) < 0) {
	fprintf(stderr, "Error reading tile at (%llu, %llu)\n", x, y);
	return NULL;
      }

      bool verbose = true;
      
      // Copy the tile data into the image buffer
      for (uint64_t ty = 0; ty < tile_height; ty++) {
	for (uint64_t tx = 0; tx < tile_width; tx++) {
	  uint64_t ind = (y + ty) * width + x + tx;
	  if (x + tx < width && y + ty < height) {

	    ++m_pix;
	    if (verbose && m_pix % static_cast<uint64_t>(1e9) == 0)
	      std::cerr << "...working on pixel: " <<
		AddCommas(static_cast<uint64_t>(m_pix)) << std::endl;

	    uint64_t tile_ind = ty * tile_width + tx;
	   
	    switch (mode) {
	    case 3:
	      ind *= 3;
	      tile_ind *= 3;
	      assert(static_cast<uint8_t*>(data)[ind] == 0);	      
	      memcpy(static_cast<uint8_t*>(data)+ind, static_cast<uint8_t*>(tile) + tile_ind, 3);
	      break;
	    case 8:
	      assert(static_cast<uint8_t*>(data)[ind] == 0);	      
	      memcpy(static_cast<uint8_t*>(data)+ind, static_cast<uint8_t*>(tile) + tile_ind, 1);
	      break;
	      //static_cast<uint8_t*>(data)[ind] =
	      //  static_cast<uint8_t*>(tile)[ty * tile_width + tx];
	    case 32:
	      static_cast<uint32_t*>(data)[ind] =
		static_cast<uint32_t*>(tile)[ty * tile_width + tx];
	    break;
          default:
	    fprintf(stderr, "not able to understand mode %d\n", mode);
	    assert(false);
	  }
	  }
	} // tile x loop
      } // tile y loop

    } // image x loop
  } // image y loop
  
  _TIFFfree(tile);

  // and put it back
  TIFFSetDirectory(m_tif, tmp_dir); 
  
  return data;
  
}

void* TiffIFD::__lined_read() {

  // store it and then move it back. Kludgy
  uint16_t tmp_dir = TIFFCurrentDirectory(m_tif);
  TIFFSetDirectory(m_tif, dir);
  
  size_t mode = GetMode();

  uint8_t* data = static_cast<uint8_t*>(__alloc());
  
  uint64_t ls = TIFFScanlineSize(m_tif);

  // allocate memory for a single line
  uint8_t* buf = static_cast<uint8_t*>(_TIFFmalloc(TIFFScanlineSize(m_tif)));

  size_t offset = 0;
  uint64_t m_pix = 0;
  for (uint64_t y = 0; y < height; y++) {

    // Read the line
    if (TIFFReadScanline(m_tif, buf, y) < 0) {
      fprintf(stderr, "Error reading line at row %llu\n", y);
      return NULL;
    }

    // data is already row-majored in TIFF, so just copy it right over
    memcpy(data + offset, buf, ls);

    offset += ls;

  }
  
  _TIFFfree(buf);

  // and put it back
  TIFFSetDirectory(m_tif, tmp_dir); 
  
  return data;
}




