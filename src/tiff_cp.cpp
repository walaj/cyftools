#ifdef HAVE_TIFFLIB
#include "tiff_cp.h"
#include <string.h>
#include <sstream>
#include <iostream>

// copied from libtiff tiffcp.c
#define	CopyField(tag, v) \
    if (TIFFGetField(in, tag, &v)) TIFFSetField(out, tag, v)
#define	CopyField2(tag, v1, v2) \
    if (TIFFGetField(in, tag, &v1, &v2)) TIFFSetField(out, tag, v1, v2)
#define	CopyField3(tag, v1, v2, v3) \
    if (TIFFGetField(in, tag, &v1, &v2, &v3)) TIFFSetField(out, tag, v1, v2, v3)
#define	CopyField4(tag, v1, v2, v3, v4) \
    if (TIFFGetField(in, tag, &v1, &v2, &v3, &v4)) TIFFSetField(out, tag, v1, v2, v3, v4)

static void
cpTag(TIFF* in, TIFF* out, uint16_t tag, uint16_t count, TIFFDataType type)
{
  switch (type) {
  case TIFF_SHORT:
    if (count == 1) {
      uint16_t shortv;
      CopyField(tag, shortv);
    } else if (count == 2) {
      uint16_t shortv1, shortv2;
      CopyField2(tag, shortv1, shortv2);
    } else if (count == 4) {
      uint16_t *tr, *tg, *tb, *ta;
      CopyField4(tag, tr, tg, tb, ta);
    } else if (count == static_cast<uint16_t>(-1)) {
      uint16_t shortv1;
      uint16_t* shortav;
      CopyField2(tag, shortv1, shortav);
    }
    break;
  case TIFF_LONG:
    { uint32_t longv;
      CopyField(tag, longv);
    }
    break;
  case TIFF_RATIONAL:
    if (count == 1) {
      float floatv;
      CopyField(tag, floatv);
    } else if (count == static_cast<uint16_t>(-1)) {
      float* floatav;
      CopyField(tag, floatav);
    }
		break;
  case TIFF_ASCII:
		{ char* stringv;
		  CopyField(tag, stringv);
		}
		break;
  case TIFF_DOUBLE:
    if (count == 1) {
      double doublev;
      CopyField(tag, doublev);
    } else if (count == static_cast<uint16_t>(-1)) {
      double* doubleav;
      CopyField(tag, doubleav);
    }
    break;
  default:
    TIFFError(TIFFFileName(in),
	      "Data type %d is not supported, tag %d skipped.",
	      tag, type);
  }
}

int tiffcp(TIFF* in, TIFF* out, bool verbose) {
  // one of:
  // COMPRESSION_NONE, COMPRESSION_PACKBITS, COMPRESSION_JPEG
  // COMPRESSION_CCITTFAX3, COMPRESSION_CCITTFAX4, COMPRESSION_LZW
  // COMPRESSION_ADOBE_DEFLATE
  // see libtiff function processCompressOptions(char* opt)
  // https://github.com/LuaDist/libtiff/blob/43d5bd6d2da90e9bf254cd42c377e4d99008f00b/tools/tiffcp.c#L324
  // default is just to copy from previous image
  uint16_t compression = static_cast<uint16_t>(-1);

  // force fill order to be either FILLORDER_LSB2MSB, FILLORDER_MSB2LSB
  // default is just 0 which is don't change it in the new image
  uint16_t fillorder= 0;  

  // image orientation
  uint16_t orientation;

  // if we are dealing with jpegs....
  int jpegcolormode = JPEGCOLORMODE_RGB;

  // tile width and length
  // default is copy from input
  uint32_t tilewidth = 0; 
  uint32_t tilelength = 0; 
  
  // output should be tiled
  // default is just to copy the input
  int outtiled = -1;

  // planar config
  // default is copy from input
  uint16_t config = static_cast<uint16_t>(-1);

  // number of rows per strip
  // one of: PLANARCONFIG_SEPARATE, PLANARCONFIG_CONTIG
  // default is copy from input
  uint32_t rowsperstrip = 0;

  // jpeg quality
  int quality = 75;

  // predictor used if COMPRESSION_DEFLATE
  uint16_t predictor = static_cast<uint16_t>(-1);

  // option used for COMPRESSION_CCITTFAX3
  uint32_t g3opts = static_cast<uint32_t>(-1);
  
  uint16_t bitspersample, samplesperpixel;
  copyFunc cf;
  uint32_t width, length;
  struct cpTag* p;
  
  CopyField(TIFFTAG_IMAGEWIDTH, width);
  CopyField(TIFFTAG_IMAGELENGTH, length);
  CopyField(TIFFTAG_BITSPERSAMPLE, bitspersample);
  CopyField(TIFFTAG_SAMPLESPERPIXEL, samplesperpixel);

  //
  // compression
  //
  if (compression != static_cast<uint16_t>(-1))
    TIFFSetField(out, TIFFTAG_COMPRESSION, compression);
  else
    CopyField(TIFFTAG_COMPRESSION, compression);
  
  if (compression == COMPRESSION_JPEG) {
    uint16_t input_compression, input_photometric;
    
    if (TIFFGetField(in, TIFFTAG_COMPRESSION, &input_compression)
	&& input_compression == COMPRESSION_JPEG) {
      TIFFSetField(in, TIFFTAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB);
    }
    if (TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &input_photometric)) {
      if(input_photometric == PHOTOMETRIC_RGB) {
	if (jpegcolormode == JPEGCOLORMODE_RGB)
	  TIFFSetField(out, TIFFTAG_PHOTOMETRIC,
		       PHOTOMETRIC_YCBCR);
	else
	  TIFFSetField(out, TIFFTAG_PHOTOMETRIC,
		       PHOTOMETRIC_RGB);
      } else
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC,
		     input_photometric);
    }
  }
  else if (compression == COMPRESSION_SGILOG
	   || compression == COMPRESSION_SGILOG24)
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC,
		 samplesperpixel == 1 ?
		 PHOTOMETRIC_LOGL : PHOTOMETRIC_LOGLUV);
  else
    CopyTag(TIFFTAG_PHOTOMETRIC, 1, TIFF_SHORT);

  // copy the fill order
  if (fillorder != 0)
    TIFFSetField(out, TIFFTAG_FILLORDER, fillorder);
  else
    CopyTag(TIFFTAG_FILLORDER, 1, TIFF_SHORT);
  
  //
  // Will copy `Orientation' tag from input image
  //
  TIFFGetFieldDefaulted(in, TIFFTAG_ORIENTATION, &orientation);
  switch (orientation) {
  case ORIENTATION_BOTRIGHT:
  case ORIENTATION_RIGHTBOT:
    TIFFWarning(TIFFFileName(in), "using bottom-left orientation");
    orientation = ORIENTATION_BOTLEFT;
    // fall thru...
  case ORIENTATION_LEFTBOT: 
  case ORIENTATION_BOTLEFT:
    break;
  case ORIENTATION_TOPRIGHT:
  case ORIENTATION_RIGHTTOP:
  default:
    TIFFWarning(TIFFFileName(in), "using top-left orientation");
    orientation = ORIENTATION_TOPLEFT;
    //
  case ORIENTATION_LEFTTOP: 
  case ORIENTATION_TOPLEFT:
    break;
  }
  TIFFSetField(out, TIFFTAG_ORIENTATION, orientation);
  

  /*
   * Choose tiles/strip for the output image according to
   * the command line arguments (-tiles, -strips) and the
   * structure of the input image.
   */
  if (outtiled == -1)
    outtiled = TIFFIsTiled(in);
  if (outtiled) {
    /*
     * Setup output file's tile width&height.  If either
     * is not specified, use either the value from the
     * input image or, if nothing is defined, use the
     * library default.
     */
    if (!tilewidth) {
      TIFFGetField(in, TIFFTAG_TILEWIDTH, &tilewidth);
    }
    if (!tilelength)
      TIFFGetField(in, TIFFTAG_TILELENGTH, &tilelength);
    TIFFDefaultTileSize(out, &tilewidth, &tilelength);
    TIFFSetField(out, TIFFTAG_TILEWIDTH, tilewidth);
    TIFFSetField(out, TIFFTAG_TILELENGTH, tilelength);
  } else {
    /*
     * RowsPerStrip is left unspecified: use either the
     * value from the input image or, if nothing is defined,
     * use the library default.
     */
    if (rowsperstrip == 0) {
      if (!TIFFGetField(in, TIFFTAG_ROWSPERSTRIP,
			&rowsperstrip)) {
	rowsperstrip =
	  TIFFDefaultStripSize(out, rowsperstrip);
      }
      if (rowsperstrip > length)
	rowsperstrip = length;
    }
    else if (rowsperstrip == static_cast<uint32_t>(-1))
      rowsperstrip = length;
    TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
  }
  if (config != static_cast<uint16_t>(-1))
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, config);
  else
    CopyField(TIFFTAG_PLANARCONFIG, config);
  if (samplesperpixel <= 4)
    CopyTag(TIFFTAG_TRANSFERFUNCTION, 4, TIFF_SHORT);
  CopyTag(TIFFTAG_COLORMAP, 4, TIFF_SHORT);


  //
  // compression
  //
  switch (compression) {
  case COMPRESSION_JPEG:
    TIFFSetField(out, TIFFTAG_JPEGQUALITY, quality);
    TIFFSetField(out, TIFFTAG_JPEGCOLORMODE, jpegcolormode);
    break;
  case COMPRESSION_LZW:
  case COMPRESSION_ADOBE_DEFLATE:
  case COMPRESSION_DEFLATE:
    if (predictor != static_cast<uint16_t>(-1))
      TIFFSetField(out, TIFFTAG_PREDICTOR, predictor);
    else
      CopyField(TIFFTAG_PREDICTOR, predictor);
    break;
  case COMPRESSION_CCITTFAX3:
  case COMPRESSION_CCITTFAX4:
    if (compression == COMPRESSION_CCITTFAX3) {
      if (g3opts != static_cast<uint32_t>(-1))
	TIFFSetField(out, TIFFTAG_GROUP3OPTIONS,
		     g3opts);
      else
	CopyField(TIFFTAG_GROUP3OPTIONS, g3opts);
    } else
      CopyTag(TIFFTAG_GROUP4OPTIONS, 1, TIFF_LONG);
    CopyTag(TIFFTAG_BADFAXLINES, 1, TIFF_LONG);
    CopyTag(TIFFTAG_CLEANFAXDATA, 1, TIFF_LONG);
    CopyTag(TIFFTAG_CONSECUTIVEBADFAXLINES, 1, TIFF_LONG);
    CopyTag(TIFFTAG_FAXRECVPARAMS, 1, TIFF_LONG);
    CopyTag(TIFFTAG_FAXRECVTIME, 1, TIFF_LONG);
    CopyTag(TIFFTAG_FAXSUBADDRESS, 1, TIFF_ASCII);
    break;
  }

  // not sure what this does
  uint32_t len32;
  void** data;
  if (TIFFGetField(in, TIFFTAG_ICCPROFILE, &len32, &data))
    TIFFSetField(out, TIFFTAG_ICCPROFILE, len32, data);

  // not sure what this does
  uint16_t ninks;
  const char* inknames;
  if (TIFFGetField(in, TIFFTAG_NUMBEROFINKS, &ninks)) {
    TIFFSetField(out, TIFFTAG_NUMBEROFINKS, ninks);
    if (TIFFGetField(in, TIFFTAG_INKNAMES, &inknames)) {
      int inknameslen = strlen(inknames) + 1;
      const char* cp = inknames;
      while (ninks > 1) {
	cp = strchr(cp, '\0');
	if (cp) {
	  cp++;
	  inknameslen += (strlen(cp) + 1);
	}
	ninks--;
      }
      TIFFSetField(out, TIFFTAG_INKNAMES, inknameslen, inknames);
    }
  }

  // not sure what this does
  unsigned short pg0, pg1;
  int pageNum = 0;
  if (TIFFGetField(in, TIFFTAG_PAGENUMBER, &pg0, &pg1)) {
    if (pageNum < 0) /* only one input file */
      TIFFSetField(out, TIFFTAG_PAGENUMBER, pg0, pg1);
    else 
      TIFFSetField(out, TIFFTAG_PAGENUMBER, pageNum++, 0);
  }

  size_t ntags = sizeof (tags) / sizeof (tags[0]);
  for (p = tags; p < &tags[ntags]; p++) {
    if (verbose)
      fprintf(stderr, "copying tag: %-21s %-5d %-20s\n",
	      tiffTagNames.at(p->tag).c_str(),
	      p->count, dataTypeNames.at(p->type).c_str());
    CopyTag(p->tag, p->count, p->type);
  }

  // this starts a rabbit hole of copying the actual image data
  // cf = pickCopyFunc(in, out, bitspersample, samplesperpixel);

  return 0;
}

int GetIntTag(TIFF* tif, int tag) {
  uint32_t t;
  if (!TIFFGetField(tif, tag, &t))
    return -1;
  return t;
}

std::string tiffprint(TIFF* tif) {

  std::stringstream ss;
  
  /*
  size_t ntags = sizeof (tags) / sizeof (tags[0]);
  struct cpTag* p;
  for (p = tags; p < &tags[ntags]; p++) {
    fprintf(stderr, "tag: %-21s %-5d %-20s\n",
	    tiffTagNames.at(p->tag).c_str(),
	    p->count, dataTypeNames.at(p->type).c_str());

    short s_s = 0;
    long ll = 0;
    char buf[100024]; 

    switch (p->type) {
    case TIFF_LONG:
      if (TIFFGetField(tif, p->tag, ll))
	fprintf(stderr, "LONG tag: %-21s %-20ul\n",
		tiffTagNames.at(p->tag).c_str(), ll);
      break;
    case TIFF_SHORT:
      if (TIFFGetField(tif, p->tag, ll))
	fprintf(stderr, "SHORT tag: %-21s %-20d\n",
		tiffTagNames.at(p->tag).c_str(), ll);
      break;
    case TIFF_ASCII:
      if (TIFFGetField(tif, p->tag, buf))
	fprintf(stderr, "ASCII tag: %-21s %-20s\n",
		tiffTagNames.at(p->tag).c_str(), buf);
      break;
    default:
      ;
      //fprintf(stderr, "UNKNOWN TAG\n");
       
    }
  }  
  */
  
  ss << "TIFF name: " << TIFFFileName(tif) << std::endl;

  // *** uint32		tif_flags;
  // sets a number of options for this TIFF
  // see https://github.com/LuaDist/libtiff/blob/43d5bd6d2da90e9bf254cd42c377e4d99008f00b/libtiff/tiffiop.h#L100

  ss << "------ Image flags ------" << std::endl;  
  ss << "  TIFF Tiled: " << TIFFIsTiled(tif) << std::endl;
  ss << "  TIFF Byte swapped: " << TIFFIsByteSwapped(tif) << std::endl;
  ss << "  TIFF Upsampled: " << TIFFIsUpSampled(tif) << std::endl;
  ss << "  TIFF BigEndian: " << TIFFIsBigEndian(tif) << std::endl;
  ss << "  TIFF MSB-2-LSB: " << TIFFIsMSB2LSB(tif) << std::endl;


  ss << "------ Image tags ------" << std::endl;
  ss << "  Image size: " <<  GetIntTag(tif, TIFFTAG_IMAGEWIDTH) << "x" <<
    GetIntTag(tif, TIFFTAG_IMAGELENGTH) << std::endl;

  if (GetIntTag(tif, TIFFTAG_EXTRASAMPLES) > 0)
    ss << "  Extra samples: " << GetIntTag(tif, TIFFTAG_EXTRASAMPLES) << std::endl;
  ss << "  Samples per pixel: " << GetIntTag(tif, TIFFTAG_SAMPLESPERPIXEL) << std::endl;
  ss << "  Bits per sample: " << GetIntTag(tif, TIFFTAG_BITSPERSAMPLE) << std::endl;
  ss << "  Planar config: " << GetIntTag(tif, TIFFTAG_PLANARCONFIG) << std::endl;    
  ss << "------ Internals -------" << std::endl;
  ss << "  Current offset: 0x" << std::hex << TIFFCurrentDirOffset(tif) << std::dec << std::endl;

  if (!TIFFIsTiled(tif)) {
    // *** uint32		tif_row;	/* current scanline */
    ss << "  Current row: " << TIFFCurrentRow(tif) << std::endl;
    // *** tstrip_t	tif_curstrip;	/* current strip for read/write */
    ss << "  Current strip: " << TIFFCurrentStrip(tif) << std::endl;
  }
  
  // *** tdir_t tif_curdir;	/* current directory (index) */
  ss << "  Current dir: " << TIFFCurrentDirectory(tif) << std::endl;
  
  ss << "  BigTiff: " << (TIFFIsBigTIFF(tif) ? "True" : "False") << std::endl;

  // *** int tif_mode;	/* open mode (O_*) */
  // https://github.com/LuaDist/libtiff/blob/43d5bd6d2da90e9bf254cd42c377e4d99008f00b/contrib/acorn/convert#L145
  // O_RDONLY 0
  // O_WRONLY 1
  // O_RDWR 2
  // O_APPEND 8
  // O_CREAT		0x200
  // O_TRUNC		0x40
  ss << "  TIFF read/write mode: " << TIFFGetMode(tif) << std::endl;
  
  // *** int tif_fd -- open file descriptor
  ss << "  TIFF I/O descriptor: " << TIFFFileno(tif) << std::endl;
  
  if (TIFFIsTiled(tif)) {
    ss << "------ Tile Info -----" << std::endl;    
    ///
    /// tile support
    ///
    ss << "  Tiles: " << GetIntTag(tif, TIFFTAG_TILEWIDTH) << "x" <<
      GetIntTag(tif, TIFFTAG_TILELENGTH) << std::endl;
    if (TIFFCurrentTile(tif) != static_cast<uint32_t>(-1))
      ss << "  Current tile: " << TIFFCurrentTile(tif) << std::endl;
    
    // return number of bytes of row-aligned tile
    ss << "  Tile size (bytes): " << std::dec << TIFFTileSize(tif) << std::endl;
    
    // return the number of bytes in each row of a tile
    ss << "  Tile row size (bytes): " << TIFFTileRowSize(tif) << std::endl;
    
    // compute the number of tiles in an image
    ss << "  Tile count: " << TIFFNumberOfTiles(tif) << std::endl;
    
    }
    
    // *** toff_t	 tif_nextdiroff;	/* file offset of following directory */
    
    // *** toff_t*  tif_dirlist;	/* list of offsets to already seen */
  
    // *** uint16	 tif_dirnumber;  /* number of already seen directories */
    
    // *** TIFFDirectory	tif_dir;	/* internal rep of current directory */

  // *** TIFFHeader	tif_header;	/* file's header block */
  /* https://github.com/LuaDist/libtiff/blob/43d5bd6d2da90e9bf254cd42c377e4d99008f00b/libtiff/tiff.h#L95
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

  return ss.str();

}
#endif
