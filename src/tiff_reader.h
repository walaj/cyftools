#ifndef TIFF_READER_H
#define TIFF_READER_H

#include <string>
#include <cassert>
#include <vector>
#include <tiffio.h>
#include <iostream>
#include <memory>

#include "tiff_ifd.h"

// this class does not store pixel data, but
// does contain the TIFF pointer to the original image.
// it will *never* change the properties of the TIFF pointer
class TiffReader {

 public:

  friend class TiffWriter;
  
  // create an empty TiffReader
  TiffReader() {}

  // destroy and close
  ~TiffReader() {}
  
  // create a new reader from a filename
  TiffReader(const char* c);

  // return the number of IFDs
  size_t NumDirs() const;

  // print it for debugging
  friend std::ostream& operator<<(std::ostream& out, const TiffReader& o);

  TiffIFD CurrentIFD() const {
    assert(curr_ifd < m_ifds.size());
    return m_ifds.at(curr_ifd);
  }

  void light_mean();

  void print();

  void print_means();

  uint32_t width() const;
  uint32_t height() const;

  TIFF* get() const { return m_tif.get(); }
  
 private:
  
  std::string m_filename;

  size_t m_num_dirs;

  std::shared_ptr<TIFF> m_tif;

  std::vector<TiffIFD> m_ifds;

  size_t curr_ifd = 0;
  
};

#endif
