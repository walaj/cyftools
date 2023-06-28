#include "tiff_header.h"

TiffHeader::TiffHeader(const char* c) {

  m_filename = std::string(c);

  __construct_header();
  
}

TiffHeader::TiffHeader(const std::string& fn) {
  m_filename = fn;
  __construct_header();
}

int TiffHeader::__construct_header() {
  
  // messing around with the tiff header
  std::ifstream inputFile;
  inputFile.open(m_filename.c_str());
  
  // check that file is open
  if (!inputFile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }

  // allocate the buffe
  m_data.reset(new char[HEADER_BUFF]);
  //  m_data = std::make_unique<char[]>(HEADER_BUFF);
  
  //m_data = std::unique_ptr<char[]> = new char[1024];
  
  // stream in the header
  inputFile.read(m_data.get(), HEADER_BUFF);

  // get the tiff id
  memcpy(&m_tid, m_data.get() + 2, sizeof(uint16_t));

  // set the offsets
  __get_offsets();

  // make sure endianness is specified
  __endianess_error_check();

  // check that the bigTiff main 8 byte header is properly formatted
  if (m_tid == 43)
    __big_tiff_error_check();
  
  // close the file
  inputFile.close();

  return 0;
}

// get the offsets using a header
int TiffHeader::__get_offsets() {

  assert(m_tid);
  
  size_t ostart;
  size_t offset_len = 0;
  
  if (m_tid == 42) {
    m_offset_start = 4;
    m_offset_len = 4;
  } else if (m_tid == 43) {
    m_offset_start = 8;
    m_offset_len = 8;
  } else {
    std::cerr << "ERROR: TIF FORMAT IS NOT PROPERLY SPECIFIED (should be 42 or 43), value is " << m_tid << std::endl;
    return 1;
  }
  
  memcpy(&m_first_offset, m_data.get() + m_offset_start, m_offset_len);
  return 0;
}

// make sure that the bigTIFF header is properly formatted
int TiffHeader::__big_tiff_error_check() const {
  
  // ensure that for bigtiff that bytes with offset 6-7 are 0
  uint16_t _zero;
  memcpy(&_zero, m_data.get() + 6, sizeof(uint16_t));
  if (_zero) {
    std::cerr << "ERROR: Improperly formatted BigTIFF: Offsets 6-7 should be zero. They are: " << _zero << std::endl;
    return 1;
  }
  return 0;
  
}

// check that the TIFF header specifies an endianness
int TiffHeader::__endianess_error_check() const {

  
  // check that endianness is either big or little
  if (!std::equal(_big, _big + 2, m_data.get()) && !std::equal(_little, _little + 2, m_data.get())){
    std::cerr << "ERROR: FIRST TWO BYTES MUST BE \"II\" FOR LITTLE-ENDIAN OR \"MM\" FOR BIG-ENDIAN" << std::endl;
    return 1;
  }
  return 0;
}


// print the header to stdout
void TiffHeader::view_stdout() const {
  
  // disply endian-ness
  std::cout << "Endian ID: " << std::string(m_data.get(), 2) << " -- " <<
    (std::equal(_little, _little + 2, m_data.get()) ? "little" : "big") << std::endl;
  
  // tiff id
  std::cout << "TIFF iD: " << m_tid << " -- " << (m_tid == 42 ? "32" : "64") <<
    "-bit TIFF" << std::endl;
  
  // offsets
  uint64_t first_offset; // offset to first IFD (in bytes)
  std::cout << "First IFD offset: 0x" << std::hex <<
    (unsigned int)m_first_offset << std::dec << std::endl;
  
  return;
}
