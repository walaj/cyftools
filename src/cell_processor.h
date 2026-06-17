#pragma once

#include "cell_header.h"
#include "cell_row.h"
#include "cysift.h"
#include "cell_selector.h"
#include "polygon.h"

#include <set>
#include <map>
#include <string>
#include <cstdint>
#include <cassert>
#include <random>

#include <cereal/types/vector.hpp>
#include <cereal/archives/portable_binary.hpp>

#include "cell_archive.h"  // OutArchive: writes cereal or CYF behind one interface

// operation type for SelectProcessor
enum optype {
  GREATER_THAN,
  LESS_THAN,
  EQUAL_TO,
  GREATER_THAN_OR_EQUAL,
  LESS_THAN_OR_EQUAL
};

typedef std::pair<optype, float> SelectOp;
typedef std::vector<SelectOp> SelectOpVec;
typedef std::unordered_map<std::string, SelectOpVec> SelectOpMap;
typedef std::unordered_map<int, SelectOpVec> SelectOpNumMap;

const float dummy_float = 23.28323130082;

class CellProcessor {
 public:

  static const int NO_WRITE_CELL = 0;
  static const int WRITE_CELL = 1;
  static const int SAVE_CELL = 2;
  static const int SAVE_NODATA_CELL = 3;
  static const int SAVE_NODATA_NOGRAPH_CELL = 4;  
  static const int HEADER_NO_ACTION = 0;
  static const int ONLY_WRITE_HEADER = 1;
  static const int SAVE_HEADER = 2;
  
  virtual ~CellProcessor() = default;
  
  virtual int ProcessHeader(CellHeader& header) = 0;

  virtual int ProcessLine(Cell& cell) = 0;

  void SetupOutputStream() {

    // Format follows the output extension (.cyf text, .byf binary); ".ocyf"
    // throws (cereal is read-only); stdout "-" defaults to binary.
    const cyf::OutFormat fmt = cyf::formatForPath(m_output_file);
    if (m_output_file == "-") {
      m_archive = std::make_unique<OutArchive>(std::cout, fmt);
    } else {
      m_os = std::make_unique<std::ofstream>(m_output_file, std::ios::binary);
      m_archive = std::make_unique<OutArchive>(*m_os, fmt);
    }
   assert(m_archive);
  }
  
  void OutputLine(const Cell& cell) const { 
    assert(m_archive);
    (*m_archive)(cell);
  }

  void SetCommonParams(const std::string& output_file,
		       const std::string& cmd,
		       bool verbose) {
    m_output_file = output_file;
    m_cmd = cmd;
    m_verbose = verbose;
  }
  
protected:

  // copy of the header to manipiulate and output
  CellHeader m_header;

  // name of the output file. empty if none
  std::string m_output_file;

  // string representation of the input command,
  // to track in PG tag in header
  std::string m_cmd;
  
  // common archiver objects
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<OutArchive> m_archive;

  // increase verbosity
  bool m_verbose = false;

  // count the lines as they come, for verbosity
  size_t m_count = 0; 
};

class LineProcessor {
 public:
  virtual ~LineProcessor() = default;
  
  virtual int ProcessHeader(CellHeader& header) = 0;
  
  virtual int ProcessLine(const std::string& line) = 0;
};

class ScatterProcessor : public CellProcessor {

public:
  
  void SetParams(int width, int height, int seed) {
    m_width = width;
    m_height = height;
    m_seed = seed;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

private:

  int m_width = 0;
  int m_height = 0;
  int m_seed = 42;
  
  std::mt19937 m_gen; // random number generator
  std::uniform_int_distribution<> m_wdistrib;
  std::uniform_int_distribution<> m_hdistrib;   
  
};

class DebugProcessor : public CellProcessor {

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

private:

  size_t m_cell_id = 0;
  
};

class RescaleProcessor : public CellProcessor {

public:
  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

private:
  
};


class MagnifyProcessor : public CellProcessor {

public:
  void SetParams(float factor) {
    m_factor = factor;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

private:

  float m_factor = 1;
  
};

class HallucinateProcessor : public CellProcessor {

public:

  void SetParams(size_t n, int seed) {
    m_nphenotypes = n;
    m_seed = seed;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

private:

  size_t m_nphenotypes = 0;

  int m_seed = 42;
  
  std::unordered_set<size_t> m_to_remove;

  std::mt19937 m_gen; // random number generator
  std::uniform_int_distribution<> m_distrib;

};


class CellCountProcessor : public CellProcessor {

public:

  void SetParams(const std::vector<uint64_t>& af) {
    m_additional_flags = af;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void EmitCell() const;
  
 private:

  std::vector<size_t> m_counts;

  std::vector<uint64_t> m_additional_flags;

  size_t m_num_marker_tags = 0;

};

class HeadProcessor : public CellProcessor {

public:

  void SetParams(size_t n) {
    m_n = n;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  
 private:
  
  size_t m_n = 10; // head limit
  size_t m_current_n = 0; // current count
  
};

// Cut processor
class CutProcessor : public CellProcessor {
  
 public:
  
  void SetParams(const std::unordered_set<std::string>& include) {
    m_include = include;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:
  
  std::unordered_set<std::string> m_include;
  
  std::unordered_set<size_t> m_data_to_keep;

};

class CleanProcessor : public CellProcessor {

public:

  void SetParams(bool clean_programs, bool clean_meta, bool clean_marker,
		 bool clean_cflags, bool clean_pflags,
		 cy_uint cflag, cy_uint pflag) {
    
    m_clean_programs = clean_programs;
    m_clean_meta = clean_meta;
    m_clean_marker = clean_marker;
    m_clean_cflags = clean_cflags;
    m_clean_pflags = clean_pflags;
    m_cflag = cflag;
    m_pflag = pflag;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;


private:

  cy_uint m_pflag = 0;
  cy_uint m_cflag = 0;
  bool m_clean_programs  = false; 
  bool m_clean_meta   = false; 
  bool m_clean_marker = false;
  bool m_clean_pflags = false;
  bool m_clean_cflags = false;

  std::unordered_set<size_t> m_to_remove;

};

// Adds (or merges into) a single header tag, then streams every cell through
// unchanged. Backs `cyftools addtag` — e.g. an @IM image/acquisition record.
class AddTagProcessor : public CellProcessor {

public:

  // Build one tag from a class + ID field + "KEY:VALUE" fields (used by `addtag`).
  // tag_id is the tag's ID field (sample id for @IM/@SA).
  void SetParams(uint8_t tag_type, const std::string& tag_id,
                 const std::vector<std::string>& fields) {
    std::string data;
    for (size_t i = 0; i < fields.size(); ++i) { if (i) data += "\t"; data += fields[i]; }
    m_tags.emplace_back(tag_type, tag_id, data);
  }

  // Add several pre-built tags at once (used by `addroi`).
  void SetTags(const std::vector<Tag>& tags) { m_tags = tags; }

  // ROI add policy for `addroi`: how to handle a file that already has @RO tags.
  // ROI_REQUIRE aborts (with a message) unless ADD/OVERWRITE is chosen explicitly.
  enum RoiMode : uint8_t { ROI_NONE = 0, ROI_REQUIRE, ROI_ADD, ROI_OVERWRITE };
  void SetRoiMode(RoiMode m) { m_roi_mode = m; }
  bool aborted() const { return m_aborted; }

  // Set a sample group/category label, stored as the @SA tag's GP field. It is
  // merged into the file's existing @SA tag; if there is none, a fresh @SA is
  // created keyed by default_sa_id (e.g. the input filename stem).
  void SetGroup(const std::string& label, const std::string& default_sa_id) {
    m_group = label; m_group_default_id = default_sa_id;
  }

  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

private:

  std::vector<Tag> m_tags;
  std::string m_group;             // sample group/category (@SA GP field), if set
  std::string m_group_default_id;  // @SA id to use when the file has no @SA tag
  RoiMode m_roi_mode = ROI_NONE;   // addroi existing-@RO policy
  bool m_aborted = false;          // set when ROI_REQUIRE hits existing ROIs
};

// Reads polygons from the @RO header tags and sets one cflag/pflag bit on every
// cell that falls inside a matching region. Backs `cyftools flagroi` — the
// "lasso in the viewer -> a flag bit on disk" write-back primitive.
// Multiply every header @RO polygon coordinate by a factor (e.g. microns-per-
// pixel), then stream cells through unchanged. Header-only transform.
class ScaleRoiProcessor : public CellProcessor {

public:

  // Per @RO vertex transform, applied as scale -> invert -> offset.
  void SetParams(double factor, double xoff, double yoff, bool flipx, bool flipy) {
    m_factor = factor; m_xoff = xoff; m_yoff = yoff; m_flipx = flipx; m_flipy = flipy;
  }

  int ProcessHeader(CellHeader& header) override;
  int ProcessLine(Cell& cell) override;

private:

  double m_factor = 1.0, m_xoff = 0.0, m_yoff = 0.0;
  bool m_flipx = false, m_flipy = false;
};

// Remove @RO polygon tags from the header (all, or filtered by name/sample),
// then stream cells through unchanged. Header-only transform.
class ClearRoiProcessor : public CellProcessor {

public:

  void SetParams(const std::string& name_filter, long sample_filter) {
    m_name = name_filter; m_sample = sample_filter;
  }

  int ProcessHeader(CellHeader& header) override;
  int ProcessLine(Cell& cell) override;

private:

  std::string m_name;       // only @RO whose NM contains this ("" = all)
  long        m_sample = -1; // only @RO with this SA sample id (< 0 = all)
};

// Read-only header check for `cyftools validate`: pulls the @HD format fields
// (VN/SO and the required MP microns-per-pixel + UN coordinate units) so the
// caller can report and instruct. Stops after the header (no cell reading).
class ValidateProcessor : public CellProcessor {

public:

  int ProcessHeader(CellHeader& header) override;
  int ProcessLine(Cell& cell) override { return NO_WRITE_CELL; }

  std::string version, sort_order, mpp, units;
};

class FlagRoiProcessor : public CellProcessor {

public:

  // reg: 'c' or 'p' (which flag register); bit: 0-based bit index to set;
  // name_filter: only @RO whose NM contains this substring ("" = all);
  // sample_filter: only @RO with this SA sample id (< 0 = all).
  void SetParams(char reg, int bit, const std::string& name_filter, long sample_filter) {
    m_reg = reg;
    m_bit = bit;
    m_name_filter = name_filter;
    m_sample_filter = sample_filter;
  }

  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

private:

  char m_reg = 'p';
  int m_bit = 0;
  std::string m_name_filter;
  long m_sample_filter = -1;
  std::vector<Polygon> m_polys;
};

// Streams a table once, accumulating columns, and writes a dependency-free
// column-major binary "viewer pack" (magic "CYFV") for a GPU front-end to load
// with typed arrays. Backs `cyftools export`. Call finalize() after streaming.
class ExportProcessor : public CellProcessor {

public:

  void SetParams(const std::string& outpath) { m_outpath = outpath; }

  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

  void finalize();   // write the pack; call after StreamTable returns

private:

  std::string m_outpath;
  std::vector<std::string> m_markers;        // data-column names (color/gate attrs)
  std::vector<char> m_kinds;                 // per column: 'M'=marker(MA), 'C'=calculated(CA)
  std::vector<uint64_t> m_id;
  std::vector<float> m_x, m_y;
  std::vector<uint64_t> m_cflag, m_pflag;
  std::vector<std::vector<float>> m_cols;    // one f32 column per marker
};

// Cohort meta-data builder, used by `cyftools cohort`. Streams one CYF table,
// collecting per-cell x/y/cflag/pflag plus the header's @FL bit names and @SA
// sample. After StreamTable, WriteSampleJSON computes — for each cflag-defined
// region (and a synthetic "All" over every cell) — the painted tissue area (the
// union of fixed-radius discs over that region's cells, overlaps counted once)
// and the per-pflag cell density in cells/mm^2, and emits one JSON object per
// sample for a cohort-level summary file.
class CohortProcessor : public CellProcessor {

public:

  void SetParams(double radius_um, double pixel_um, size_t threads) {
    m_radius   = radius_um;
    m_pixel    = pixel_um;
    m_nthreads = threads ? threads : 1;
  }

  int ProcessHeader(CellHeader& header) override;
  int ProcessLine(Cell& cell) override;

  // Compute regions/areas/densities and write this sample's JSON object to `os`,
  // each line prefixed with `indent`. `file` is the source path (recorded as-is,
  // and used to derive the sample label if the header carries no @SA). Call after
  // StreamTable returns; progress is logged to stderr.
  void WriteSampleJSON(std::ostream& os, const std::string& file,
                       const std::string& indent);

  size_t CellCount() const { return m_x.size(); }

private:

  // params
  double m_radius   = 20.0;   // paint disc radius, microns
  double m_pixel    = 1.0;    // raster pixel size, microns
  size_t m_nthreads = 1;

  // collected per-cell columns
  std::vector<float>    m_x, m_y;
  std::vector<uint64_t> m_cflag, m_pflag;

  // sample identity: id high 32 bits are the sample id ((sample_id<<32)|cell_id)
  uint64_t    m_sample_id = 0;
  std::string m_sample_name;        // from @SA, else derived from the filename
  std::string m_sample_group;       // from @SA GP field (set via `addtag --group`)

  // bit -> symbolic name, from @FL tags (RG:cflag / RG:pflag), when present
  std::map<int,std::string> m_cflag_names;
  std::map<int,std::string> m_pflag_names;

  // pflag bit -> @MA marker name by column order (the convention the viewers use:
  // 1st marker = bit 0). Fallback for files with no @FL pflag declarations.
  std::map<int,std::string> m_pflag_auto;
};

// a single collection of value to calculate a mean on
// so it's one data type (e.g. dist_tls) and one group by
// e.g. tls_id
struct SumElement {

  std::vector<double> elems;
  
  // add a new
  void Add(float val) {
    elems.push_back(val);
  }

  float Mean() const {

    double n = static_cast<double>(elems.size());
    
    if (n == 0) {
      std::cerr << "Error - Unable to get mean, no elements" << std::endl;
      assert(false);
    }
    
    // calculate the means
    double mean = 0;
    for (size_t i = 0; i < elems.size(); i++) {
      mean += elems.at(i) / n;
    }
    
    return mean;
  }
};

// a collection of SumElements for one group (e.g. tls_id)
// but for all data
struct SumGroup {
  
  // the key is the unique group-by value (e.g. tls id)
  std::vector<SumElement> data;

  // add a single data point
  void Add(size_t cols_index, float value) {

    // already have this group, just add to it
    if (cols_index < data.size()) {
      ; 
    // a new one
    } else if (cols_index == data.size()) {
      data.push_back(SumElement());
    } else {
      assert(false);
    }
    
    data[cols_index].Add(value);

  }

  // construct a cell with means
  Cell EmitCell() const {

    assert(data.size());
    
    Cell cell;

    // calculate the means and push
    for (size_t i = 0; i < data.size(); i++) {
      cell.cols.push_back(data.at(i).Mean());
    }

    // zero the hard data since they are meaningless
    cell.x = 0;
    cell.y = 0;  
    cell.pflag = 0;
    cell.cflag = 0;
    cell.id = 0;

    return cell;
  }
  
};

class AverageProcessor : public CellProcessor {

public:

  void SetParams(const std::string& by) {
    m_group_by = by;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void EmitCells() const;
  
private:

  // here, the key is the group and each SumGroup
  // holds all of the numerical data for all columns with that group-by
  std::unordered_map<float, SumGroup> allgroups;

  // group means by this column
  std::string m_group_by;
  int m_group_by_id = -1;
  
  // make a set to get uniqe tls ids
  //int m_tls_column = -1;
  //std::set<int> m_tls_set;
  
};


// Pheno processor
class PhenoProcessor : public CellProcessor {
  
 public:
  
  void SetParams(PhenoMap p, float scale,
		 float random_scale) {
    m_p = p;
    m_scale = scale;
    m_random_scale = random_scale;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:

  // map connecting marker names to indicies
  // marker id -> { column index (all data cols), pflag bit (@MA order) }
  std::unordered_map<std::string, std::pair<size_t, size_t>> m_marker_map;

  // map of all of the gates (string - pair<float,float>)
  PhenoMap m_p;

  float m_scale = 1;

  float m_random_scale = 0;

};

// Tumor processor
/*class TumorProcessor : public CellProcessor {
  
 public:
  
  void SetParams(int n, cy_uint flag, float frac) {
    m_n = n;
    m_flag = flag;
    m_frac = frac;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:

  int m_n;
  cy_uint m_flag;
  float m_frac;
 
};
*/

class SummaryProcessor : public CellProcessor {

 public:
  
  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void Print();
  
 private:
  std::vector<float> m_min;
  std::vector<float> m_max;
  std::vector<float> m_mean;  
  size_t m_count = 0;
};
  
// Count processor
class CountProcessor : public CellProcessor {
  
 public:
  
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void PrintCount();
  
 private:
  size_t m_count = 0;
};


// filter processor
class TrimProcessor : public CellProcessor {

public:
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

};

// Filter processor
class FilterProcessor : public CellProcessor {
  
 public:
  
  void SetFlagParams(CellSelector select, bool mark1, bool mark2) {
    m_select = select;
    m_mark1 = mark1;
    m_mark2 = mark2;
  }

  void SetFieldParams(const SelectOpMap criteria, bool or_toggle) {
    m_criteria = criteria;
    m_or_toggle = or_toggle;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:

  CellSelector m_select;

  bool m_mark1 = false;
  bool m_mark2 = false;  
  
  // field selectors
  SelectOpMap m_criteria;
  SelectOpNumMap m_criteria_int;  // created by ProcessHeader
  bool m_or_toggle = false;
};

// Hidden processor just to confirm flags are cleared
class MarkCheckProcessor : public CellProcessor {
   public:
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
};

// Process flag logic ops
class FlagsetProcessor : public CellProcessor {
  
  public:

  void SetParams(cy_uint pflag_set, cy_uint cflag_set,
		 cy_uint pflag_clear, cy_uint cflag_clear) {
    m_pflag_to_set = pflag_set;
    m_cflag_to_set = cflag_set;
    m_pflag_to_clear = pflag_clear;
    m_cflag_to_clear = cflag_clear;    
  }

  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

private:

  cy_uint m_pflag_to_set = 0;
  cy_uint m_cflag_to_set = 0;  
  cy_uint m_pflag_to_clear = 0;
  cy_uint m_cflag_to_clear = 0;
  
};

// Divide processor
class DivideProcessor : public CellProcessor {

 public:

  void SetParams(const std::string& num,
		 const std::string& den,
		 float div)  {
    m_numer_string = num;
    m_denom_string = den;
    m_div_zero = div;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  std::string m_numer_string; // name numerator
  std::string m_denom_string; // name denominator
  size_t m_numer = -1; // index of numerator in data column
  size_t m_denom = -1; // index of denominator in data column  
  float m_div_zero; // value to be given to divde-by-zero
  size_t m_count = 0; // for verbose
  size_t m_existing_column = -1; // if the column for the answer already exists, overwrite
};

// Log10 processor
class LogProcessor : public CellProcessor {
  
 public:

  void SetParams()  { }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  std::unordered_set<size_t> m_to_log;

  bool m_bool_warning_emitted = false;
  
};


class ROIProcessor : public CellProcessor {

 public:
  
  void SetParams(bool label,
		 const std::vector<Polygon>& rois,
		 bool bl) {
    m_label = label;
    m_rois = rois;
    m_blacklist_remove = bl;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

 private:

  std::vector<Polygon> m_rois;

  bool m_label;

  bool m_blacklist_remove = false;

  bool m_roi_region = false; // set true if any polygon/roi is a region
};

class ViewProcessor : public CellProcessor { 

 public:
  
  void SetParams(bool print_header,
		 bool header_only,
		 bool rheader,
		 bool adjacent,
		 bool crevasse,
		 int round,
		 const std::unordered_set<std::string>& include,
		 bool tabprint,
		 bool strict_cut,
		 bool list_markers
		 ) {
    
    m_print_header = print_header;
    m_header_only = header_only;
    m_round = round;
    m_crevasse = crevasse;
    m_to_view = include;
    m_csv_header = rheader;
    m_adjacent = adjacent;
    m_tabprint = tabprint;
    m_strict_cut = strict_cut;
    m_list_markers = list_markers;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  bool m_tabprint = false;
  
  bool m_csv_header = false;

  bool m_adjacent = false;

  bool m_crevasse = false;
  
  bool m_header_only;
  
  bool m_print_header;

  bool m_strict_cut = false; // skip the CellID etc

  bool m_list_markers = false;
  
  // number of integers to round output to
  int m_round;

  std::unordered_set<std::string> m_to_view;

  std::unordered_set<size_t> m_to_view_indicies;
};

class BuildProcessor : public CellProcessor { 

 public:
  
  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

};


class SubsampleProcessor : public CellProcessor {

public:

  void SetParams(float rate, int seed, std::unordered_set<uint32_t> samplenum) {
    m_rate = rate;
    m_seed = seed;
    m_samplenum = samplenum;
  }

  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

private:

  std::unordered_set<uint32_t> m_samplenum;
  float m_rate = 0;

  size_t m_count = 0;
  size_t m_kept = 0;
  
  int m_seed = 42;
  std::mt19937 m_gen; // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> m_dis; 
  
  
};


class OffsetProcessor : public CellProcessor {

public:

  void SetParams(float x, float y) {
    m_x = x;
    m_y = y;
  }

  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

private:

  float m_x = 0;
  float m_y = 0;
  
};

class CheckProcessor : public CellProcessor {
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
private:
  
};

class FlipProcessor : public CellProcessor {
  
public:
  
  void SetParams(float x, float y,
		 float xmax, float ymax) {
    m_x = x;
    m_y = y;
    m_xmax = xmax;
    m_ymax = ymax;
    
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
private:
  
  float m_x = 0;
  float m_y = 0;
  float m_xmax = 0;
  float m_ymax = 0;
  
};

class CatProcessor : public CellProcessor { 

 public:
  
  void SetParams(int offset,
		 int sample) {
    
    m_offset = offset;
    m_sample = sample;
    
  }

  void SetOffset(size_t offset) { m_offset = offset; }

  void SetSample(size_t sample) { m_sample = sample; }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

  size_t GetMaxCellID() const { return m_max_cellid; }
  
private:
  
  size_t m_cellid_index = static_cast<size_t>(-1);
  
  size_t m_offset;
  size_t m_sample;

  bool m_error_emitted = false; // print error just once per file
  
  // the header to compare with (and print if needed)
  CellHeader m_master_header;
  bool m_master_set = false;

  size_t m_max_cellid = 0;

  std::vector<size_t> m_graph_indicies;
  
};

class CerealProcessor : public LineProcessor { 

 public:
  
  void SetParams(const std::string& filename,
		 const std::string& cmd,
		 uint32_t sampleid) {
    m_filename = filename;
    m_cmd = cmd;
    m_sampleid = sampleid;
  }

  // which are the columns for X and Y
  void SetXInd(int x) { m_x_index = x; }
  void SetYInd(int y) { m_y_index = y; }
  void SetIDInd(int i) { m_id_index = i; }

  // Per-CSV-column conversion plan, built by StreamTableCSV. Each column maps to
  // one action; gate (*p) columns carry the pflag bit to flip in m_gate_bit.
  enum ColKind : uint8_t {
    COL_IGNORE = 0, COL_ID, COL_X, COL_Y,
    COL_MARKER,   // marker intensity   -> cols (markers, in column order)
    COL_CA,       // non-marker numeric -> cols (metas,   in column order)
    COL_IC,       // -> cols (meta) AND set cflag bit 5 when non-zero
    COL_GATE,     // -> set pflag bit m_gate_bit[col] when non-zero
    COL_REGION,   // -> set cflag bit 3 when value == 1
  };
  void SetColumnPlan(std::vector<uint8_t> kinds, std::vector<int> gate_bits) {
    m_col_kind = std::move(kinds);
    m_gate_bit = std::move(gate_bits);
  }

  // Required @HD scale tags stamped on convert: MP (microns per pixel) and UN
  // (coordinate units, "micron" or "pixel").
  void SetScale(const std::string& mpp, const std::string& units) {
    m_mpp = mpp; m_units = units;
  }

  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(const std::string& line) override;
  //int ProcessLine(const std::string& line) override;

private:

  uint32_t m_cellid = 0;
  uint32_t m_sampleid = 0;  
  std::vector<float> vec1;

  std::string m_cmd;
  std::string m_filename;
  std::string m_mpp;      // @HD MP (microns per pixel); empty if not set
  std::string m_units;    // @HD UN ("micron"/"pixel"); empty if not set

  // index tracking
  int m_x_index = -1;
  int m_y_index = -1; 
  int m_id_index = -1; 
  
  CellHeader m_header;
  
  // per-CSV-column conversion plan (see ColKind); empty -> legacy positional parse
  std::vector<uint8_t> m_col_kind;
  std::vector<int>     m_gate_bit;

  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<OutArchive> m_archive;

};

/*class RadialProcessor : public CellProcessor { 

 public:

  void SetParams(const std::vector<cy_uint>& inner, const std::vector<cy_uint>& outer,
		 const std::vector<cy_uint>& logor, const std::vector<cy_uint>& logand,
		 const std::vector<std::string>& label, bool normalize) {

    // check the radial geometry parameters
    assert(inner.size());
    assert(inner.size() == outer.size());
    assert(inner.size() == logor.size());
    assert(inner.size() == logand.size());
    assert(inner.size() == label.size());    

    m_inner = inner;
    m_outer = outer;
    m_logor = logor;
    m_logand = logand;
    m_label = label;
    m_normalize = normalize;
  }
  
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:
  
  std::vector<cy_uint> m_inner, m_outer, m_logor, m_logand;
  std::vector<std::string> m_label;
  bool m_normalize;
};
*/
