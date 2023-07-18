#ifdef HAVE_TIFFLIB
#include "tiff_reader.h"

void TiffReader::print_means() {

  for (int i = 0; i < m_num_dirs; i++) {

    TiffIFD& ifd = m_ifds[i];
    
    uint8_t mode = ifd.GetMode();

    switch(mode) {
    case 8:
      std::cout << "Mean: " << ifd.mean()[3] << std::endl;
      break;
    case 3:
      std::vector<double> d = ifd.mean();
      std::cout << "Mean: (R) " << d[0] << " (G) " << d[1] <<
	" (B) " << d[2] << std::endl;
      break;
    }
    
  }
    

}

TiffReader::TiffReader(const char* c) {
  
  m_tif = std::shared_ptr<TIFF>(TIFFOpen(c, "rm"), TIFFClose);
  
  // Open the input TIFF file
  if (!m_tif) {
    fprintf(stderr, "Error opening %s for reading\n", c);
    return;
  }
  
  // set the number of directories
  m_num_dirs = TIFFNumberOfDirectories(m_tif.get());

   // set the filename
   m_filename = std::string(c);

   // set the IFDs
   for (size_t i = 0; i < m_num_dirs; i++) {
     TIFFSetDirectory(m_tif.get(), i);
     m_ifds.push_back(TiffIFD(m_tif.get()));
   }
   TIFFSetDirectory(m_tif.get(), 0);
   
 }

void TiffReader::print() {

  std::cout << "-- Image file: " << m_filename << std::endl;
  std::cout << "-- Num dirs: " << m_num_dirs << std::endl;
  for (auto& i : m_ifds) {
    std::cout << " ----- IFD " << i.dir << " -----" << std::endl;
    i.print();
  }
  
}

size_t TiffReader::NumDirs() const {
  return m_num_dirs;
}


uint32_t TiffReader::width() const {
  uint32_t width;
  TIFFGetField(m_tif.get(), TIFFTAG_IMAGEWIDTH, &width);
  return width;
}

uint32_t TiffReader::height() const {
  uint32_t height;
  TIFFGetField(m_tif.get(), TIFFTAG_IMAGELENGTH, &height);
  return height;

}

#endif
