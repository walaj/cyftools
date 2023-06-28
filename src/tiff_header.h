#ifndef TIFF_HEADER_H
#define TIFF_HEADER_H

#include <iostream>
#include <cassert>
#include <string>
#include <memory>
#include <fstream>
#include <cstring>

#define HEADER_BUFF 1024

/* 
 // TIFF header
typedef	struct {
	uint16	tiff_magic;	// magic number (defines byte order) 
#define TIFF_MAGIC_SIZE		2
	uint16	tiff_version;	// TIFF version number 
#define TIFF_VERSION_SIZE	2
	uint32	tiff_diroff;	// byte offset to first directory 
#define TIFF_DIROFFSET_SIZE	4
} TIFFHeader;
*/


/*
// Default Read/Seek/Write definitions (tiffiop.h)
#ifndef ReadOK
#define	ReadOK(tif, buf, size) \
	(TIFFReadFile(tif, (tdata_t) buf, (tsize_t)(size)) == (tsize_t)(size))
#endif
#ifndef SeekOK
#define	SeekOK(tif, off) \
	(TIFFSeekFile(tif, (toff_t) off, SEEK_SET) == (toff_t) off)
#endif
#ifndef WriteOK
#define	WriteOK(tif, buf, size) \
	(TIFFWriteFile(tif, (tdata_t) buf, (tsize_t) size) == (tsize_t) size)
#endif
*/


 // show whether the image file is big or little endian
static const char _little[2] = {'I','I'};
static const char _big[2] = {'M','M'};

class TiffHeader {

 public:

  TiffHeader(const char* c);
  
  TiffHeader(const std::string& fn);
  
  char* data() { return m_data.get(); }

  // send to human readable stdout
  void view_stdout() const;
    
 private:

  // TiffHeader class is directly from libtiff
  // see https://github.com/LuaDist/libtiff/blob/43d5bd6d2da90e9bf254cd42c377e4d99008f00b/libtiff/tiff.h#L95
  //TiffHeader m_tif_header;

  std::string m_filename;
  
  // smart pointer to the 
  std::unique_ptr<char[]> m_data;

  // the tiff id (e.g. 42 for 32-bit, 43 for 64-bit)
  uint16_t m_tid = 0;

  // offset parameters, namely where does the location of the first offset value
  // start and how long is the offset value (e.g. 32 or 64 bits)
  size_t m_offset_start, m_offset_len;

  // location of the first IFD;
  uint64_t m_first_offset;

  // call internal routines to make the header
  int __construct_header();
  
  // get the location of the first IFD
  int __get_offsets();

  // check that the endianness is specified
  int __endianess_error_check() const;

  // check that the 3rd word of a BigTIFF is 0
  int __big_tiff_error_check() const;
  
};

#endif
