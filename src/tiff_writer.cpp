#ifdef HAVE_TIFFLIB
#include "tiff_writer.h"

int TiffWriter::SetTag(uint32_t tag, ...) {

  va_list ap;
  int status;
  
  va_start(ap, tag);
  status = TIFFVSetField(m_tif.get(), tag, ap);
  va_end(ap);
  return (status);
  
}

TiffWriter::TiffWriter(const char* c) {

  m_tif = std::shared_ptr<TIFF>(TIFFOpen(c, "w8"), TIFFClose);
  
  // Open the output TIFF file
  if (!m_tif) {
    fprintf(stderr, "Error opening %s for writing\n", c);
    return;
  }

  m_filename = std::string(c);

  // only no compression currently supported
  TIFFSetField(m_tif.get(), TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  
}

void TiffWriter::CopyFromReader(const TiffReader& tr) {

  // copy the tags from reader tif to writer tif
  tiffcp(tr.m_tif.get(), m_tif.get(), false);

  // only no compression currently supported
  TIFFSetField(m_tif.get(), TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  
}

void TiffWriter::SetTile(int h, int w) {
  
  if (!TIFFSetField(m_tif.get(), TIFFTAG_TILEWIDTH, w)) {
    fprintf(stderr, "ERROR: unable to set tile width %ul on writer\n", w);
    return;
  }
  if (!TIFFSetField(m_tif.get(), TIFFTAG_TILELENGTH, h)) {
    fprintf(stderr, "ERROR: unable to set tile height %ul on writer\n", h);
    return;
  }

  // turn off the tiled flag
  if (h == 0 || w == 0) {
    TIFFSetTiledOff(m_tif.get());
    TIFFUnsetField(m_tif.get(), TIFFTAG_TILEWIDTH);
    TIFFUnsetField(m_tif.get(), TIFFTAG_TILELENGTH);    
  }
  
}

bool TiffWriter::isTiled() {

  assert(m_tif);
  return(TIFFIsTiled(m_tif.get()));
}

int TiffWriter::Write(const TiffImage& ti) {

  uint32_t m_width, m_height;
  assert(TIFFGetField(m_tif.get(), TIFFTAG_IMAGEWIDTH, &m_width));
  assert(TIFFGetField(m_tif.get(), TIFFTAG_IMAGELENGTH, &m_height));  

  if (ti.m_width != m_width) {
    std::cerr << "Warning: Image writer has width " << m_width <<
      " and image width " << ti.m_width << ". You may want to run writer.UpdateDims(image)" << std::endl;
  }
  if (ti.m_height != m_height) {
    std::cerr << "Warning: Image writer has height " << m_height <<
      " and image height " << ti.m_height << ". You may want to run writer.UpdateDims(image)" << std::endl;
  }

  if (isTiled())
    return(__tiled_write(ti));
  else
    return(__lined_write(ti));
      
}

void TiffWriter::UpdateDims(const TiffImage& ti) {

  assert(TIFFSetField(m_tif.get(), TIFFTAG_IMAGEWIDTH, ti.m_width));
  assert(TIFFSetField(m_tif.get(), TIFFTAG_IMAGELENGTH, ti.m_height));  
  
}

int TiffWriter::__tiled_write(const TiffImage& ti) const {

  // sanity check
  assert(TIFFIsTiled(m_tif.get()));
  assert(TIFFTileSize(m_tif.get()));
  
  // pull tile and image dims directly from libtiff writer
  uint32_t o_tile_width, o_tile_height;
  assert(TIFFGetField(m_tif.get(), TIFFTAG_TILEWIDTH, &o_tile_width));
  assert(TIFFGetField(m_tif.get(), TIFFTAG_TILELENGTH, &o_tile_height));  
  uint32_t m_width, m_height;
  assert(TIFFGetField(m_tif.get(), TIFFTAG_IMAGEWIDTH, &m_width));
  assert(TIFFGetField(m_tif.get(), TIFFTAG_IMAGELENGTH, &m_height));  

  uint8_t mode = GetMode(); 

  if (mode != ti.GetMode()) {
    fprintf(stderr, "Mode of writer: %d, mode of image: %d\n", mode, ti.GetMode());
    fprintf(stderr, "Will attempt to continue, but not sure this is what you want\n");
  }
  
  // error check that tile dimensions
  if (o_tile_width % 16) {
    fprintf(stderr, "ERROR: output tile dims must be multiple of 16, it is: %d\n", o_tile_width);
    return 1;
  }
  if (o_tile_height % 16) {
    fprintf(stderr, "ERROR: output tile dims must be multiple of 16, it is: %d\n", o_tile_height);
    return 1;
  }

  // Allocate a buffer for the output tiles
  void* buf = _TIFFmalloc(TIFFTileSize(m_tif.get()));
  if (buf == NULL) {
    fprintf(stderr, "Error allocating memory for tiles\n");
    return 1;
  }

  // loop through the output image tiles and write the pixels
  // from the raster
  uint64_t x, y;

  for (y = 0; y < m_height; y += o_tile_height) {
    for (x = 0; x < m_width; x += o_tile_width) {
      // Fill the tile buffer with data
      for (uint64_t ty = 0; ty < o_tile_height; ty++) {
	for (uint64_t tx = 0; tx < o_tile_width; tx++) {
	  uint64_t xp = x + tx; // want 64 bits to be able to store really big d_pos, since this is the nth pixel
	  uint64_t yp = y + ty;
	  uint64_t d_pos = yp * m_width + xp; // index of pixel in the entire image
	  uint64_t t_pos = ty * o_tile_width + tx; // index of pixel in this particular tile
	  if (xp < m_width && yp < m_height) {
	    switch(mode) {
	    case 8:
	      static_cast<uint8_t*>(buf)[t_pos] = ti.pixel(xp, yp, PIXEL_GRAY);
	      break;
	      //case 32:
	      //static_cast<uint32_t*>(buf)[t_pos] = ti.pixel<uint32_t>(xp, yp, PIXEL_GRAY);
	      //break;
	    case 3:
	      t_pos = t_pos * 3;
	      static_cast<uint8_t*>(buf)[t_pos]     = ti.pixel(xp, yp, PIXEL_RED);
	      static_cast<uint8_t*>(buf)[t_pos + 1] = ti.pixel(xp, yp, PIXEL_GREEN);
	      static_cast<uint8_t*>(buf)[t_pos + 2] = ti.pixel(xp, yp, PIXEL_BLUE);
	      break;
	    }
	  }
	}
      }

      // Write the tile to the TIFF file
      // this function will automatically calculate memory size from TIFF tags
      if (TIFFWriteTile(m_tif.get(), buf, x, y, 0, 0) < 0) { 
	fprintf(stderr, "Error writing tile at (%llu, %llu)\n", x, y);
	return 1;
      }
    }
  }

  // Clean up
  _TIFFfree(buf);
  
  return 0;
}

void TiffWriter::MatchTagsToRaster(const TiffImage& ti) {

  TIFFSetField(m_tif.get(), TIFFTAG_IMAGEWIDTH, ti.m_width);
  TIFFSetField(m_tif.get(), TIFFTAG_IMAGELENGTH, ti.m_height);
  
}

int TiffWriter::__lined_write(const TiffImage& ti) const {

  uint64_t ls = TIFFScanlineSize(m_tif.get());
  assert(ls);

  uint32_t m_width, m_height;
  assert(TIFFGetField(m_tif.get(), TIFFTAG_IMAGEWIDTH, &m_width));
  assert(TIFFGetField(m_tif.get(), TIFFTAG_IMAGELENGTH, &m_height));  

  if (m_width != ti.m_width) {
    fprintf(stderr, "Warning: Image dim %d x %d does not match writer tag %d x %d\n", ti.m_width, ti.m_height, m_width, m_height);
    fprintf(stderr, "         Attempting to write anyway\n");
  }
  
  uint8_t m_mode = GetMode(); 
  
  if (m_mode != ti.GetMode()) {
    fprintf(stderr, "Mode of writer: %d, mode of image: %d\n", m_mode, ti.GetMode());
    fprintf(stderr, "Will attempt to continue, but not sure this is what you want\n");
  }

  // this should take into account the mode, since TIFFScanlineSize should
  // factor in samples per pixel. But we'll always work with it as a uint8_t array
  // for ease
  //uint8_t* buf = (uint8_t*)_TIFFmalloc(TIFFScanlineSize(otif));

  for (uint64_t y = 0; y < m_height; y++) {

    // copy directly from m_data since the order of writing
    // in a scanline file is same as order that m_data is stored
    // this should be really clean because it doesn't make me guess the
    // number of bytes per pixel to deal with
    if (TIFFWriteScanline(m_tif.get(), static_cast<uint8_t*>(ti.m_data) + y * ls, y, 0) < 0) {
      fprintf(stderr, "Error writing line row %llu\n", y);
      return 1;
    }
    
    //          memcpy(buf, m_data + y * ls, ls);
  }
  return 0;
}



/*void TiffWriter::SetField(ttag+t tag, ...) {

  TIFFSetField(m_tif.get(), tag, ...);

  }*/

/*TiffWriter::~TiffWriter() {

  if (m_tif != NULL) 
    TIFFClose(m_tif);
  m_tif = NULL;
  }*/


uint8_t TiffWriter::GetMode() const {

  uint16_t samples_per_pixel, bits_per_sample;

  assert(TIFFGetField(m_tif.get(), TIFFTAG_BITSPERSAMPLE, &bits_per_sample));
  assert(TIFFGetField(m_tif.get(), TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel));
  
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
#endif
