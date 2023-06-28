#ifndef TIFF_COPY_H
#define TIFF_COPY_H


/*
 * Copyright (c) 1988-1997 Sam Leffler
 * Copyright (c) 1991-1997 Silicon Graphics, Inc.
 *
 *  Revised:  2/18/01 BAR -- added syntax for extracting single images from
 *                          multi-image TIFF files.
 *
 *    New syntax is:  sourceFileName,image#
 *
 * image# ranges from 0..<n-1> where n is the # of images in the file.
 * There may be no white space between the comma and the filename or
 * image number.
 *
 *    Example:   tiffcp source.tif,1 destination.tif
 *
 * Copies the 2nd image in source.tif to the destination.
 *
 *****
 * Permission to use, copy, modify, distribute, and sell this software and 
 * its documentation for any purpose is hereby granted without fee, provided
 * that (i) the above copyright notices and this permission notice appear in
 * all copies of the software and related documentation, and (ii) the names of
 * Sam Leffler and Silicon Graphics may not be used in any advertising or
 * publicity relating to the software without the specific, prior written
 * permission of Sam Leffler and Silicon Graphics.
 * 
 * THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY 
 * WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 * IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
 * ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
 * OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF 
 * LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE 
 * OF THIS SOFTWARE.
 */

#include <tiffio.h>
#include <string>
#include <map>

#include <map>
#include <string>

// Lookup table mapping TIFF tag values to TIFFTAG names
const std::map<int, std::string> tiffTagNames = {
    { TIFFTAG_SUBFILETYPE,	"SubFileType" },
    { TIFFTAG_OSUBFILETYPE,	"OldSubFileType" },
    { TIFFTAG_IMAGEWIDTH,	"ImageWidth" },
    { TIFFTAG_IMAGELENGTH,	"ImageLength" },
    { TIFFTAG_BITSPERSAMPLE,	"BitsPerSample" },
    { TIFFTAG_COMPRESSION,	"Compression" },
    { TIFFTAG_PHOTOMETRIC,	"Photometric" },
    { TIFFTAG_THRESHHOLDING,	"Threshholding" },
    { TIFFTAG_CELLWIDTH,	"CellWidth" },
    { TIFFTAG_CELLLENGTH,	"CellLength" },
    { TIFFTAG_FILLORDER,	"FillOrder" },
    { TIFFTAG_DOCUMENTNAME,	"DocumentName" },
    { TIFFTAG_IMAGEDESCRIPTION,	"ImageDescription" },
    { TIFFTAG_MAKE,		"Make" },
    { TIFFTAG_MODEL,		"Model" },
    { TIFFTAG_STRIPOFFSETS,	"StripOffsets" },
    { TIFFTAG_ORIENTATION,	"Orientation" },
    { TIFFTAG_SAMPLESPERPIXEL,	"SamplesPerPixel" },
    { TIFFTAG_ROWSPERSTRIP,	"RowsPerStrip" },
    { TIFFTAG_STRIPBYTECOUNTS,	"StripByteCounts" },
    { TIFFTAG_MINSAMPLEVALUE,	"MinSampleValue" },
    { TIFFTAG_MAXSAMPLEVALUE,	"MaxSampleValue" },
    { TIFFTAG_XRESOLUTION,	"XResolution" },
    { TIFFTAG_YRESOLUTION,	"YResolution" },
    { TIFFTAG_PLANARCONFIG,	"PlanarConfig" },
    { TIFFTAG_PAGENAME,		"PageName" },
    { TIFFTAG_XPOSITION,	"XPosition" },
    { TIFFTAG_YPOSITION,	"YPosition" },
    { TIFFTAG_FREEOFFSETS,	"FreeOffsets" },
    { TIFFTAG_FREEBYTECOUNTS,	"FreeByteCounts" },
    { TIFFTAG_GRAYRESPONSEUNIT,	"GrayResponseUnit" },
    { TIFFTAG_GRAYRESPONSECURVE,"GrayResponseCurve" },
    { TIFFTAG_GROUP3OPTIONS,	"Group3Options" },
    { TIFFTAG_GROUP4OPTIONS,	"Group4Options" },
    { TIFFTAG_RESOLUTIONUNIT,	"ResolutionUnit" },
    { TIFFTAG_PAGENUMBER,	"PageNumber" },
    { TIFFTAG_COLORRESPONSEUNIT,"ColorResponseUnit" },
    { TIFFTAG_TRANSFERFUNCTION,	"TransferFunction" },
    { TIFFTAG_SOFTWARE,		"Software" },
    { TIFFTAG_DATETIME,		"DateTime" },
    { TIFFTAG_ARTIST,		"Artist" },
    { TIFFTAG_HOSTCOMPUTER,	"HostComputer" },
    { TIFFTAG_PREDICTOR,	"Predictor" },
    { TIFFTAG_WHITEPOINT,	"Whitepoint" },
    { TIFFTAG_PRIMARYCHROMATICITIES,"PrimaryChromaticities" },
    { TIFFTAG_COLORMAP,		"Colormap" },
    { TIFFTAG_HALFTONEHINTS,	"HalftoneHints" },
    { TIFFTAG_TILEWIDTH,	"TileWidth" },
    { TIFFTAG_TILELENGTH,	"TileLength" },
    { TIFFTAG_TILEOFFSETS,	"TileOffsets" },
    { TIFFTAG_TILEBYTECOUNTS,	"TileByteCounts" },
    { TIFFTAG_BADFAXLINES,	"BadFaxLines" },
    { TIFFTAG_CLEANFAXDATA,	"CleanFaxData" },
    { TIFFTAG_CONSECUTIVEBADFAXLINES, "ConsecutiveBadFaxLines" },
    { TIFFTAG_SUBIFD,		"SubIFD" },
    { TIFFTAG_INKSET,		"InkSet" },
    { TIFFTAG_INKNAMES,		"InkNames" },
    { TIFFTAG_NUMBEROFINKS,	"NumberOfInks" },
    { TIFFTAG_DOTRANGE,		"DotRange" },
    { TIFFTAG_TARGETPRINTER,	"TargetPrinter" },
    { TIFFTAG_EXTRASAMPLES,	"ExtraSamples" },
    { TIFFTAG_SAMPLEFORMAT,	"SampleFormat" },
    { TIFFTAG_SMINSAMPLEVALUE,	"SMinSampleValue" },
    { TIFFTAG_SMAXSAMPLEVALUE,	"SMaxSampleValue" },
    { TIFFTAG_JPEGPROC,		"JPEGProcessingMode" },
    { TIFFTAG_JPEGIFOFFSET,	"JPEGInterchangeFormat" },
    { TIFFTAG_JPEGIFBYTECOUNT,	"JPEGInterchangeFormatLength" },
    { TIFFTAG_JPEGRESTARTINTERVAL,"JPEGRestartInterval" },
    { TIFFTAG_JPEGLOSSLESSPREDICTORS,"JPEGLosslessPredictors" },
    { TIFFTAG_JPEGPOINTTRANSFORM,"JPEGPointTransform" },
    { TIFFTAG_JPEGTABLES,       "JPEGTables" },
    { TIFFTAG_JPEGQTABLES,	"JPEGQTables" },
    { TIFFTAG_JPEGDCTABLES,	"JPEGDCTables" },
    { TIFFTAG_JPEGACTABLES,	"JPEGACTables" },
    { TIFFTAG_YCBCRCOEFFICIENTS,"YCbCrCoefficients" },
    { TIFFTAG_YCBCRSUBSAMPLING,	"YCbCrSubsampling" },
    { TIFFTAG_YCBCRPOSITIONING,	"YCbCrPositioning" },
    { TIFFTAG_REFERENCEBLACKWHITE, "ReferenceBlackWhite" },
    { TIFFTAG_REFPTS,		"IgReferencePoints (Island Graphics)" },
    { TIFFTAG_REGIONTACKPOINT,	"IgRegionTackPoint (Island Graphics)" },
    { TIFFTAG_REGIONWARPCORNERS,"IgRegionWarpCorners (Island Graphics)" },
    { TIFFTAG_REGIONAFFINE,	"IgRegionAffine (Island Graphics)" },
    { TIFFTAG_MATTEING,		"OBSOLETE Matteing (Silicon Graphics)" },
    { TIFFTAG_DATATYPE,		"OBSOLETE DataType (Silicon Graphics)" },
    { TIFFTAG_IMAGEDEPTH,	"ImageDepth (Silicon Graphics)" },
    { TIFFTAG_TILEDEPTH,	"TileDepth (Silicon Graphics)" },
    { 32768,			"OLD BOGUS Matteing tag" },
    { TIFFTAG_COPYRIGHT,	"Copyright" },
    { TIFFTAG_ICCPROFILE,	"ICC Profile" },
    { TIFFTAG_JBIGOPTIONS,	"JBIG Options" },
    { TIFFTAG_STONITS,		"StoNits" }
};

// Lookup table mapping integer values to TIFFDataType names
static const std::map<int, std::string> dataTypeNames = {
    {1,  "TIFF_BYTE"},
    {2,  "TIFF_ASCII"},
    {3,  "TIFF_SHORT"},
    {4,  "TIFF_LONG"},
    {5,  "TIFF_RATIONAL"},
    {6,  "TIFF_SBYTE"},
    {7,  "TIFF_UNDEFINED"},
    {8,  "TIFF_SSHORT"},
    {9,  "TIFF_SLONG"},
    {10, "TIFF_SRATIONAL"},
    {11, "TIFF_FLOAT"},
    {12, "TIFF_DOUBLE"},
    {13, "TIFF_IFD"},
    {16, "TIFF_LONG8"},
    {17, "TIFF_SLONG8"},
    {18, "TIFF_IFD8"}
};


static struct cpTag {
  uint16_t	tag;
  uint16_t	count;
  TIFFDataType type;
} tags[] = {
  { TIFFTAG_SUBFILETYPE,		1, TIFF_LONG },
  { TIFFTAG_THRESHHOLDING,	1, TIFF_SHORT },
  { TIFFTAG_DOCUMENTNAME,		1, TIFF_ASCII },
  { TIFFTAG_IMAGEDESCRIPTION,	1, TIFF_ASCII },
  { TIFFTAG_MAKE,			1, TIFF_ASCII },
  { TIFFTAG_MODEL,		1, TIFF_ASCII },
  { TIFFTAG_MINSAMPLEVALUE,	1, TIFF_SHORT },
  { TIFFTAG_MAXSAMPLEVALUE,	1, TIFF_SHORT },
  { TIFFTAG_XRESOLUTION,		1, TIFF_RATIONAL },
  { TIFFTAG_YRESOLUTION,		1, TIFF_RATIONAL },
  { TIFFTAG_PAGENAME,		1, TIFF_ASCII },
  { TIFFTAG_XPOSITION,		1, TIFF_RATIONAL },
  { TIFFTAG_YPOSITION,		1, TIFF_RATIONAL },
  { TIFFTAG_RESOLUTIONUNIT,	1, TIFF_SHORT }, 
  { TIFFTAG_SOFTWARE,		1, TIFF_ASCII },
  { TIFFTAG_DATETIME,		1, TIFF_ASCII },
  { TIFFTAG_ARTIST,		1, TIFF_ASCII },
  { TIFFTAG_HOSTCOMPUTER,		1, TIFF_ASCII },
  { TIFFTAG_WHITEPOINT,		(uint16_t) -1, TIFF_RATIONAL },
  { TIFFTAG_PRIMARYCHROMATICITIES,(uint16_t) -1,TIFF_RATIONAL },
  { TIFFTAG_HALFTONEHINTS,	2, TIFF_SHORT },
  { TIFFTAG_INKSET,		1, TIFF_SHORT },
  { TIFFTAG_DOTRANGE,		2, TIFF_SHORT },
  { TIFFTAG_TARGETPRINTER,	1, TIFF_ASCII },
  { TIFFTAG_SAMPLEFORMAT,		1, TIFF_SHORT },
  { TIFFTAG_YCBCRCOEFFICIENTS,	(uint16_t) -1,TIFF_RATIONAL },
  { TIFFTAG_YCBCRSUBSAMPLING,	2, TIFF_SHORT },
  { TIFFTAG_YCBCRPOSITIONING,	1, TIFF_SHORT },
  { TIFFTAG_REFERENCEBLACKWHITE,	(uint16_t) -1,TIFF_RATIONAL },
  { TIFFTAG_EXTRASAMPLES,		(uint16_t) -1, TIFF_SHORT },
  { TIFFTAG_SMINSAMPLEVALUE,	1, TIFF_DOUBLE },
  { TIFFTAG_SMAXSAMPLEVALUE,	1, TIFF_DOUBLE },
  { TIFFTAG_STONITS,		1, TIFF_DOUBLE },
};

#define	CopyTag(tag, count, type) cpTag(in, out, tag, count, type)

static void
cpTag(TIFF* in, TIFF* out, uint16_t tag, uint16_t count, TIFFDataType type);

#define	TRUE	1
#define	FALSE	0

// copied from libtiff tiffcp.c
typedef int (*copyFunc)
    (TIFF* in, TIFF* out, uint32_t l, uint32_t w, uint16_t samplesperpixel);

// copied from libtiff tiffcp.c
// this will help choose how to copy pixel data (e.g. rows->tiles etc)
//static	copyFunc pickCopyFunc(TIFF*, TIFF*, uint16_t, uint16_t);

int tiffcp(TIFF* in, TIFF* out, bool verbose);

std::string tiffprint(TIFF* tif);

int GetIntTag(TIFF* tif, int tag);

#endif
