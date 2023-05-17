#include "cell_column.h"
#include "cell_table.h"
#include "cell_processor.h"

#include <unistd.h> // or #include <getopt.h> on Windows systems
#include <getopt.h>

#include "cell_row.h"

namespace opt {
  static bool verbose = false;
  static std::string infile;
  static std::string quantfile;
  static std::string outfile;
  static std::string module;
  
  static std::vector<std::string> infile_vec;
  
  static bool header = false;
  static bool header_only = false;
  
  static bool csv = false; // should we print as csv instead of screen readable
  
  static std::string roifile;

  static std::string cropstring;

  static int width = 50;
  static int height = 50;
  
  static int threads = 1;

  static int seed = 1337;
  static int n = 0;

  static bool sort = false;
}

#define DEBUG(x) std::cerr << #x << " = " << (x) << std::endl

#define TVERB(msg) \
  if (opt::verbose) {		   \
    std::cerr << msg << std::endl; \
 }

static CellTable table;

static void build_table(bool convert);

static const char* shortopts = "jhHNyvr:t:a:d:g:b:q:c:s:k:n:r:w:l:x:X:o:R:f:";
static const struct option longopts[] = {
  { "verbose",                    no_argument, NULL, 'v' },
  { "threads",                    required_argument, NULL, 't' },
  { "quant-file",                 required_argument, NULL, 'q' },
  { "seed",                       required_argument, NULL, 's' },
  { "roi",                        required_argument, NULL, 'r' },  
  { "numrows",                    required_argument, NULL, 'n' },
  { "w",                          required_argument, NULL, 'w' },
  { "l",                          required_argument, NULL, 'l' },
  { "crop",                       required_argument, NULL, 'c' },
  { "cut",                        required_argument, NULL, 'x' },
  { "strict-cut",                 required_argument, NULL, 'X' },    
  { "sort",                       no_argument, NULL, 'y' },  
  { "csv",                        no_argument, NULL, 'j'},
  { NULL, 0, NULL, 0 }
};

static const char *RUN_USAGE_MESSAGE =
"Usage: cysift [module] <options> \n"
"  view       - View the cell table\n"
"  cut        - Select only given markers and metas\n"
"  cat        - Concatenate multiple samples\n"  
"  subsample  - Subsample cells randomly\n"
"  plot       - Generate an ASCII style plot\n"
"  roi        - Trim cells to a region of interest within a given polygon\n"
"  histogram  - Create a histogram of the data\n"
"  log10      - Apply a base-10 logarithm transformation to the data\n"
"  correlate  - Calculate the correlation between variables\n"
"  info       - Display information about the dataset"
"  knn        - Construct the marker space KNN graph\n"
"  spatial    - Construct the spatial KNN graph\n"
"  select     - Select by cell phenotype flags\n"
"  pheno      - Phenotype cells to set the flag\n"
"  radialdens - Calculate density of cells within a radius\n"
"\n";

static int cerealfunc(int argc, char** argv);
static int catfunc(int argc, char** argv);
static int radialdensfunc(int argc, char** argv);
static int subsamplefunc(int argc, char** argv);
static int viewfunc(int argc, char** argv);
static int infofunc(int argc, char** argv);
static int roifunc(int argc, char** argv);
static int correlatefunc(int argc, char** argv);
static int histogramfunc(int argc, char** argv);
static int plotfunc(int argc, char** argv);
static int debugfunc(int argc, char** argv);
static int cropfunc(int argc, char** argv);
static int subsamplefunc(int argc, char** argv);
static int log10func(int argc, char** argv);
static int cutfunc(int argc, char** argv);
static int knnfunc(int argc, char** argv);
static int radiusfunc(int argc, char** argv);
static int selectfunc(int argc, char** argv);
static int spatialfunc(int argc, char** argv); 
static int phenofunc(int argc, char** argv);

static void parseRunOptions(int argc, char** argv);

int main(int argc, char **argv) {
  
  // Check if a command line argument was provided
  if (argc < 2) {
    std::cerr << RUN_USAGE_MESSAGE;
    return 1;
  }

  parseRunOptions(argc, argv);

  int val = 0;
  
/*
  switch(hash(opt::module)) {
  case hash("debug"):     val = debugfunc(argc, argv); break;
  case hash("subsample"): val = subsamplefunc(argc, argv); break;
  case hash("plot"):     return plotfunc(argc, argv); 
  case hash("roi"):       val = roifunc(argc, argv); break;
  case hash("crop"):      val = cropfunc(argc, argv); break;
  case hash("correlate"): return correlatefunc(argc, argv);
  case hash("histogram"); return histogramfunc(argc, argv);
  case hash("log10"):     val = log10func(argc, argv); break;
  case hash("cut"):       val = cutfunc(argc, argv); break;
  case hash("info"):      return(infofunc(argc, argv));
  case hash("view"):      val = viewfunc(argc, argv); break;
  case hash("knn"):       val = knnfunc(argc, argv); break;
  default: assert(false);
  }
*/

  // get the module
  if (opt::module == "debug") {
    val = debugfunc(argc, argv);
  } else if (opt::module == "radialdens") {
    val = radialdensfunc(argc, argv);
  } else if (opt::module == "subsample") {
    val = subsamplefunc(argc, argv);
  } else if (opt::module == "plot") {
    return(plotfunc(argc, argv));
  } else if (opt::module == "roi") {
    val = roifunc(argc, argv);
  } else if (opt::module == "crop") {
    val = cropfunc(argc, argv);
  } else if (opt::module == "correlate") {
    return(correlatefunc(argc, argv));
  } else if (opt::module == "histogram") {
    return(histogramfunc(argc, argv));
  } else if (opt::module == "log10") {
    return(log10func(argc, argv));
  } else if (opt::module == "cut") {
    val = cutfunc(argc, argv);
  } else if (opt::module == "cereal") {
    val = cerealfunc(argc, argv);
  } else if (opt::module == "info") {
    return(infofunc(argc, argv));
  } else if (opt::module == "view") {
    return(viewfunc(argc, argv));
  } else if (opt::module == "cat") {
    return (catfunc(argc, argv));
  } else if (opt::module == "knn") {
    val = knnfunc(argc, argv);
  } else if (opt::module == "spatial") {
    val = spatialfunc(argc, argv);    
  } else if (opt::module == "select") {
    val = selectfunc(argc, argv);
  } else if (opt::module == "pheno") {
    val = phenofunc(argc, argv);
  } else {
    assert(false);
  }

  return 0;
}

// build the table into memory
static void build_table(bool convert) {

  // set table params
  if (opt::verbose)
    table.SetVerbose();
  
  if (opt::threads > 1)
    table.SetThreads(opt::threads);
    
  // read the table into memory
  table.BuildTable(opt::infile);

  if (convert) {
    if (opt::verbose)
      std::cerr << "...read table: converting graph and flag columns" << std::endl;
    table.ConvertColumns();
  }
}


/* static int radiusfunc(int argc, char** argv) {

  int radius = 20;
  bool die = false;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'r' : arg >> radius; break;
    case 'v' : opt::verbose = true; break;
    case 'h' : opt::header = true; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift radius [csvfile] -r <dist>\n"
      "  Select cells with centroids within spatial distance of r\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -r [10]                   Radius to select out\n"
      "    -v, --verbose             Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  read_table();

  return 0;
}
*/

static int knnfunc(int argc, char** argv) {

  int n = 10;
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'h' : opt::header = true; break;
    case 't' : arg >> opt::threads; break;
    case 'k' : arg >> n; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift knn [csvfile]\n"
      "  Construct the KNN graph (in marker space)\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -k [10]                   Number of neighbors\n"
      "    -t [1]                   Number of threads\n"      
      "    -v, --verbose             Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // build the table
  // but don't have to convert columns
  // since we don't use pre-existing Graph or Flags for this
  build_table(false);

  // build the KNN graph in marker-space
  table.KNN_marker(n);

  // print it
  table.PrintTable(opt::header);
  
  return 0;
}

static int catfunc(int argc, char** argv) {

  std::string samples;
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 's' : arg >> samples; break;
    case 'h' : opt::header = true; break;      
    default: die = true;
    }
  }

  optind++;
  
  // Process any remaining no-flag options
  while (optind < argc) {
    opt::infile_vec.push_back(argv[optind]);
    optind++;
  }

  //
  if (opt::verbose)
    for (const auto& v : opt::infile_vec)
      std::cerr << "...set to read " << v << std::endl;
  
  // display help if no input
  if (opt::infile_vec.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift cat [csvfile]\n"
      "  Concatenate together multiple cell tables\n"
      "    csvfile: filepaths of cell tables\n"
      "    -v, --verbose             Increase output to stderr\n"
      "    -s                        Sample numbers, to be in same number as inputs and comma-sep\n"      
      "    -h                        Output with the header\n"            
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // parse the sample numbers
  std::vector<int> sample_nums;
  std::unordered_set<std::string> tokens;
  std::stringstream ss(samples);
  std::string token;
  while (std::getline(ss, token, ',')) {
    sample_nums.push_back(std::stoi(token));
  }

  if (sample_nums.size() != opt::infile_vec.size() && sample_nums.size()) {
    throw std::runtime_error("Sample number csv line should have same number of tokens as number of input files");
  }
  
  size_t offset = 0;

  CatProcessor catp;
  catp.SetParams(opt::header, opt::header_only, offset, 0);
  
  for (size_t sample_num = 0; sample_num < opt::infile_vec.size(); sample_num++) {

    // can make a new table for each iteration, since we are dumping right to stdout
    CellTable this_table;
    
    // set table params
    if (opt::verbose)
      this_table.SetVerbose();

    if (opt::verbose)
      std::cerr << "...reading and concatenating " << opt::infile_vec.at(sample_num) << std::endl;
    
    // update the offset and sample num
    catp.SetOffset(offset);

    if (sample_nums.size())
      catp.SetSample(sample_nums.at(sample_num));
    else
      catp.SetSample(sample_num);

    // stream in the lines
    this_table.StreamTable(catp, opt::infile_vec.at(sample_num));
    
    offset = catp.GetMaxCellID() + 1; // + 1 to avoid dupes if new cellid starts at 0
  }

  return 0;
}

static int infofunc(int argc, char** argv) {

  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift info [csvfile]\n"
      "  Display basic information on the cell table\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -v, --verbose             Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }
  
  // build it into memory and then provide information
  build_table(false);

  // provide information to stdout
  std::cout << table;
  
  return 0;
}

static int cutfunc(int argc, char** argv) {

  bool strict_cut = false;
  std::string cut; // list of markers, csv separated, to cut on
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'x' : arg >> cut; break;
    case 'X' : strict_cut = true; arg >> cut; break;
    case 'n' : arg >> opt::n; break;
    case 'h' : opt::header = true; break;
    case 'H' : opt::header_only = true; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  // display help if no input
  if (opt::infile.empty() || die || cut.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift cut [csvfile] <marker1,markers2>\n"
      "  Cut the file to only certain markers\n"
      "    <file>: filepath or a '-' to stream to stdin\n"
      "    -x, --cut                 Comma-separated list of markers to cut to\n"
      "    -X, --strict-cut          Comma-separated list of markers to cut to\n"      
      "    -v, --verbose             Increase output to stderr\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }
  
  // parse the markers
  std::unordered_set<std::string> tokens;
  std::stringstream ss(cut);
  std::string token;
  while (std::getline(ss, token, ',')) {
    tokens.insert(token);
  }

  // set table params
  if (opt::verbose)
    table.SetVerbose();
  
  // setup the cut processor
  CutProcessor cutp;
  cutp.SetParams(tokens, opt::header, strict_cut);

  // process
  table.StreamTable(cutp, opt::infile);

  return 0;
}

static int log10func(int argc, char** argv)  {
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> opt::n; break;
    case 'h' : opt::header = true; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift log10 [csvfile] <options>\n"
      "  Calculate the log10 of marker intensities\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -v, --verbose             Increase output to stderr"
      "  -h                        Output with the header\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // set table params
  if (opt::verbose)
    table.SetVerbose();

  LogProcessor log;
  log.SetParams(opt::header);

  table.StreamTable(log, opt::infile);

  return 0;
}

// parse the command line options
static void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 1) 
    die = true;
  else {
    // get the module
    opt::module = argv[1];
  }

  // make sure module is implemented
  if (! (opt::module == "debug" || opt::module == "subsample" ||
	 opt::module == "plot"  || opt::module == "roi" ||
	 opt::module == "histogram" || opt::module == "log10" ||
	 opt::module == "crop"  || opt::module == "knn" ||
	 opt::module == "cat" || opt::module == "cereal" || 
	 opt::module == "correlate" || opt::module == "info" ||
	 opt::module == "cut" || opt::module == "view" ||
	 opt::module == "spatial" || opt::module == "radialdens" || 
	 opt::module == "select" || opt::module == "pheno")) {
    std::cerr << "Module " << opt::module << " not implemented" << std::endl;
    die = true;
  }

  /*
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> opt::threads; break;
    case 'c' : arg >> opt::cropstring; break;
    case 'j' : opt::csv = true; break;
    case 'x' : arg >> opt::cut; break;
    case 'X' : opt::strict_cut = true; arg >> opt::cut; break;
    case 'q' : arg >> opt::quantfile; break;
    case 's' : arg >> opt::seed; break;
    case 'n' : arg >> opt::n; break;
    case 'r' : arg >> opt::roifile; break;
    case 'h' : opt::header = true; break;
    case 'H' : opt::header_only = true; break;
    case 'w' : arg >> opt::width; break;
    case 'l' : arg >> opt::height; break;
    case 'y' : opt::sort = true; break;
    default: die = true;
    }
  }
*/
  /*
  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  */
  
  if (die) {
      std::cerr << "\n" << RUN_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

static int roifunc(int argc, char** argv) {

  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> opt::n; break;
    case 'h' : opt::header = true; break;
    case 'r' : arg >> opt::roifile;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }


  // display help if no input
  if (opt::infile.empty() || opt::roifile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift roi [csvfile] <options>\n"
      "  Subset or label the cells to only those contained in the rois\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -r                        ROI file\n"
      "  -l                        Output all cells and add \"roi\" column with ROI label\n"      
      "  -h                        Output with the header\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // read in the roi file
  std::vector<Polygon> rois = read_polygons_from_file(opt::roifile);

  if (opt::verbose)
    for (const auto& c : rois)
      std::cerr << c << std::endl;

  ROIProcessor roip;
  roip.SetParams(opt::header, false, rois);// false is placeholder for label function, that i need to implement

  table.StreamTable(roip, opt::infile);

  return 0;
  
}

static int viewfunc(int argc, char** argv) {

  int precision = -1;

  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> precision; break;
    case 't' : arg >> opt::threads; break;      
    case 'h' : opt::header = true; break;
    case 'H' : opt::header_only = true; break;
    default: die = true;
    }
  }
  
  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift view [csvfile] <options>\n"
      "  View the contents of a cell table\n" 
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -n  [-1]                  Number of decimals to keep (-1 is no change)\n"
      "  -H                        View only the header\n"      
      "  -h                        Output with the header\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // set table params
  if (opt::verbose)
    table.SetVerbose();
  
  ViewProcessor viewp;
  viewp.SetParams(opt::header, opt::header_only, precision);
  
  table.StreamTable(viewp, opt::infile);
  
  return 0;  
}

static int histogramfunc(int argc, char** argv) {

  int n_bins = 50;
  int w_bins = 50;
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> n_bins; break;
    case 'w' : arg >> w_bins; break;      
    case 'r' : arg >> opt::roifile;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift histogram [csvfile] <options>\n"
      "  Calculate the histogram of a set of markers\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -n  [50]                  Number of bins\n"
      "  -w  [50]                  Binwidth\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  //read_table();
  //table.histogram(opt::n, opt::width);

  return 0;
  
}

static int plotfunc(int argc, char** argv) {

  int length = 50;
  int width = 50;
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'l' : arg >> length; break;
    case 'w' : arg >> width; break;      
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift plot [csvfile] <options>\n"
      "  Outputs an ASCII-style plot of cell locations\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -l, --length        [50]  Height (length) of output plot, in characters\n"   
      "    -w, --width         [50]  Width of output plot, in characters\n"
      "    -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table(false);
  
  // make an ASCII plot of this
  table.PlotASCII(width, length);

  return 0;
}

static int subsamplefunc(int argc, char** argv) {

  int n = 0;
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> n; break;
    case 's' : arg >> opt::seed; break;      
    case 'h' : opt::header = true; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *SUBSAMPLE_USAGE_MESSAGE =
      "Usage: cysift subsample [csvfile] <options>\n"
      "  Subsamples a cell quantification table, randomly.\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -n, --numrows             Number of rows to subsample\n"
      "    -s, --seed         [1337] Seed for random subsampling\n"
      "    -h                        Output with the header\n"            
      "    -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << SUBSAMPLE_USAGE_MESSAGE;
    return 1;
  }

  build_table(false);
  
  // subsample
  table.Subsample(n, opt::seed);

  table.PrintTable(opt::header);
  
  return 0;
  
}

static int correlatefunc(int argc, char** argv) {

  bool sorted = false;
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 's' : sorted = true; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift correlate [csvfile] <options>\n"
      "  Outputs an ASCII-style plot of marker intensity pearson correlations\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -j                        Output as a csv file\n"
      "    -s                        Sort the output by Pearson correlation\n"
      "    -v, --verbose             Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table(false);
  
  //
  table.PrintPearson(opt::csv, opt::sort);
  
  return 0;
}

static int cropfunc(int argc, char** argv) {

  std::string cropstring;
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'c' : arg >> cropstring; break;
    case 'h' : opt::header = true; break;      
     default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  
  // display help if no input
  if (opt::infile.empty() || cropstring.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift crop [csvfile] <options>\n"
      "  Crop the table to a given rectangle (in pixels)\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    --crop                    String of form xlo,xhi,ylo,yhi\n"
      "    -h                        Output with the header\n"
      "    -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  int xlo, xhi, ylo, yhi;
  std::vector<int*> coordinates = {&xlo, &xhi, &ylo, &yhi};

  std::istringstream iss(opt::cropstring);
  std::string token;
  
  size_t idx = 0;
  while (std::getline(iss, token, ',')) {
    if (idx >= coordinates.size()) {
      throw std::runtime_error("Error: More than 4 numbers provided");
    }

    try {
      *coordinates[idx] = std::stoi(token);
    } catch (const std::invalid_argument&) {
      throw std::runtime_error("Error: Non-numeric value encountered");
    } catch (const std::out_of_range&) {
      throw std::runtime_error("Error: Numeric value out of range");
    }

    ++idx;
  }

  if (idx < coordinates.size()) {
    throw std::runtime_error("Error: Fewer than 4 numbers provided");
  }

  build_table(false);
  
  //
  table.Crop(xlo, xhi, ylo, yhi);

  return 0;
  
}

static int spatialfunc(int argc, char** argv) {

  int n = 10;
  int d = -1;
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'h' : opt::header = true; break;
    case 't' : arg >> opt::threads; break;
    case 'k' : arg >> n; break;
    case 'd' : arg >> d; break;            
    default: die = true;
    }
  }
  
  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift spatial [csvfile]\n"
      "  Construct the Euclidean KNN spatial graph\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -k [10]               Number of neighbors\n"
      "    -d [-1]               Max distance to include as neighbor (-1 = none)\n"
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table(false);

  table.KNN_spatial(n, d);

  table.PrintTable(opt::header);

  return 0;
}

static int selectfunc(int argc, char** argv) {

  uint64_t logor = 0;
  uint64_t logand = 0;
  bool lognot = false;
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'h' : opt::header = true; break;
    case 'o' : arg >> logor; break;
    case 'a' : arg >> logand; break;
    case 'N' : lognot = true; break;      
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift select [csvfile]\n"
      "  Select cells by phenotype flag\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -o                    Logical OR flags\n"
      "    -a                    Logical AND flags\n"
      "    -N                    Not flag\n"
      "    -h                    Output with the header\n"      
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // read it in 
  if (!opt::infile.empty()) {
    table = CellTable(); //opt::quantfile.c_str(), opt::verbose, opt::header_only, false);
  }

  // setup the selector processor
  SelectProcessor select;
  select.SetParams(logor, logand, lognot, opt::header);

  // process
  table.StreamTable(select, opt::infile);
  
  return 0;
}

static int phenofunc(int argc, char** argv) {

  std::string file;
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'h' : opt::header = true; break;
    case 't' : arg >> file; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift pheno [csvfile]\n"
      "  Phenotype cells (set the flags) with threshold file\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -t                    Phenotype file\n"
      "    -h                    Output with the header\n"      
      "    -v, --verbose         Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table(false);

  std::unordered_map<std::string, std::pair<float,float>> tt =
    table.phenoread(file);
  
  if (opt::verbose)
    for (const auto& c : tt)
      std::cerr << c.first << " -- " << c.second.first << "," << c.second.second << std::endl;

  table.phenotype(tt);

  table.PrintTable(opt::header);
  
  return 0;
}

static int radialdensfunc(int argc, char** argv) {

  uint64_t inner = 0;
  uint64_t outer = 20;
  uint64_t logor = 0;
  uint64_t logand = 0;
  std::string label;

  std::string file;
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'h' : opt::header = true; break;
    case 't' : arg >> opt::threads; break;
    case 'R' : arg >> inner; break;
    case 'r' : arg >> outer; break;
    case 'o' : arg >> logor; break;
    case 'a' : arg >> logand; break;
    case 'l' : arg >> label; break;
    case 'f' : arg >> file; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  if (label.empty() && file.empty())
    die = true;
  
  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift radialdens [csvfile]\n"
      "  Calculate the density of cells away from individual cells\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -r [20]               Outer radius\n"
      "    -R [0]                Inner radius\n"
      "    -o                    Logical OR flags\n"
      "    -a                    Logical AND flags\n"
      "    -l                    Label the column\n"
      "    -f                    File for multiple labels [r,R,o,a,l]\n"
      "    -h                    Output with the header\n"      
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  if (inner >= outer) {
    std::cerr << "Inner radius should be smaller than outer, or else no cells are included" << std::endl;
    return 1;
  }

  std::vector<RadialSelector> rsv;
  if (!file.empty()) {
    std::ifstream input_file(file);
   
    if (!input_file.is_open()) {
      throw std::runtime_error("Failed to open file: " + file);
    }
    
    std::string line;
    while (std::getline(input_file, line)) {
      rsv.push_back(RadialSelector(line));
    }
    
    if (opt::verbose) {
      for (const auto& rr : rsv)
	std::cerr << rr << std::endl;
    }
  }
  
  build_table(true);

  if (rsv.size()) {
    for (const auto& rr : rsv) {
      if (opt::verbose)
	std::cerr << "...working on radial density group: " << rr << std::endl;
      
      table.RadialDensity(rr.int_data[0], rr.int_data[1],
			  rr.int_data[2], rr.int_data[3],
			  rr.label);
    }
  } else {
    table.RadialDensity(inner, outer, logor, logand, label);
  }
  
  table.PrintTable(opt::header);

  return 0;
}

int debugfunc(int argc, char** argv) {
  return 1;

}

static int cerealfunc(int argc, char** argv) {
  
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'h' : opt::header = true; break;
    default: die = true;
    }
  }

  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } 
    optind++;
  }

  // display help if no input
  if (opt::infile.empty() || die) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift cereal [csvfile]\n"
      "  ***\n" 
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -h                    Output with the header\n"      
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  CerealProcessor cerp;

  table.StreamTable(cerp, opt::infile);

  return 0;
  
  std::ifstream file("cereal.bin", std::ios::binary);
  cereal::BinaryInputArchive inputArchive(file);
  
  // Keep loading objects from the archive until you reach the end of the file:
  std::vector<Cell> vec;
  size_t i = 0;
  while(file.peek() != EOF) {
    Cell cell;
    inputArchive(cell);
    vec.push_back(cell);
    if (i % 100000 == 0)
      std::cerr << i << std::endl;
    i++;
    //objects.push_back(obj);
  }
  
  return 0;
}
