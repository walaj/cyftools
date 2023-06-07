#include "cell_column.h"
#include "cell_table.h"
#include "cell_processor.h"
#include "cell_lda.h"

#include <unistd.h> // or #include <getopt.h> on Windows systems
#include <getopt.h>
#include <ctime>
#include <regex>

#include "cell_row.h"

static std::string cmd_input;
static bool die = false;

namespace opt {
  static bool verbose = false;

  // file i/o
  static std::string infile;
  static std::string outfile;
  static std::string module;
  static std::vector<std::string> infile_vec;
 
  static bool header = false;
  static bool header_only = false;

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

// process in and outfile cmd arguments
static bool in_out_process(int argc, char** argv);

// process in cmd arguments
static bool in_only_process(int argc, char** argv);

static CellTable table;

static void build_table();

static const char* shortopts = "jhHNyvmMGAr:t:a:d:g:b:c:s:k:n:r:w:l:x:X:o:R:f:";
static const struct option longopts[] = {
  { "verbose",                    no_argument, NULL, 'v' },
  { "threads",                    required_argument, NULL, 't' },
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
"  count      - Output number of cells in table\n"  
"  cut        - Select only given markers and metas\n"
"  clean      - Removes data to decrease disk size\n"
"  average    - Average all of the data columns\n"  
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
"  cereal     - Create a .cys format file from a CSV\n"
"\n";

static int ldafunc(int argc, char** argv);
static int averagefunc(int argc, char** argv);
static int cleanfunc(int argc, char** argv);
static int countfunc(int argc, char** argv);
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

  // store the cmd input
  for (int i = 0; i < argc; ++i) {
    cmd_input += argv[i];
    cmd_input += " ";
  }

  // Get the current time
  std::time_t now = std::time(nullptr);
  char timestamp[100];
  
  // Format the time as a string: YYYY-MM-DD HH:MM:SS
  std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
  
  // Append the timestamp to the string
  cmd_input = cmd_input + "\t -- " + timestamp;
  
  // get the module
  if (opt::module == "debug") {
    val = debugfunc(argc, argv);
  } else if (opt::module == "lda") {
    val = ldafunc(argc, argv);
  } else if (opt::module == "clean") {
    val = cleanfunc(argc, argv);
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
  } else if (opt::module == "average") {
    val = averagefunc(argc, argv);
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
  } else if (opt::module == "count") {
    countfunc(argc, argv);
  } else {
    assert(false);
  }

  return 0;
}

// build the table into memory
static void build_table() {

  // set table params
  if (opt::verbose)
    table.SetVerbose();
  
  if (opt::threads > 1)
    table.SetThreads(opt::threads);

  // stream into memory
  BuildProcessor buildp;
  buildp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  table.StreamTable(buildp, opt::infile);
  
  table.SetCmd(cmd_input);
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
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> opt::threads; break;
    case 'k' : arg >> n; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
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
  build_table();

  // build the KNN graph in marker-space
  table.KNN_marker(n);

  // print it
  table.PrintTable(opt::header);
  
  return 0;
}

static int ldafunc(int argc, char** argv) {

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift lda [csvfile]\n"
      "  Perform topic-model learning using Latent Dirichlet Allocation\n"
      "    <file>: filepath or a '-' to stream to stdin\n"
      "    -v, --verbose             Increase output to stderr\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // stream into memory
  BuildProcessor buildp;
  buildp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  table.StreamTable(buildp, opt::infile);

  //
  std::cerr << "...writing" << std::endl;
  table.HDF5Write("output.h5");
  
  //CellLDA lda_topic;
  //lda_topic.run();
  
  return 0;
  
}

static int cleanfunc(int argc, char** argv) {

  bool clean_graph = false;
  bool clean_markers = false;
  bool clean_meta = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'm' : clean_markers = true; break;
    case 'M' : clean_meta = true; break;
    case 'G' : clean_graph = true; break;
    case 'A' : clean_graph = true; clean_markers = true; clean_meta = true; break;      
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift clean [csvfile]\n"
      "  Clean up the data to reduce disk space\n"
      "    <file>: filepath or a '-' to stream to stdin\n"
      "    -m,          Remove all marker data\n"
      "    -M,          Remove all meta data\n"
      "    -G,          Remove all graph data\n"      
      "    -A,          Remove all data\n"      
      "    -v, --verbose             Increase output to stderr\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }
  
  CleanProcessor cleanp;
  cleanp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  cleanp.SetParams(clean_graph, clean_meta, clean_markers);  

  // process 
  if (!table.StreamTable(cleanp, opt::infile)) {
    return 1;
  }

  return 0;
  
}

static int catfunc(int argc, char** argv) {

  std::string samples;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 's' : arg >> samples; break;
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
  catp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  catp.SetParams(offset, 0);
  
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

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  // display help if no input
  if (die || in_only_process(argc, argv)) {
    
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
  build_table();
  
  // provide information to stdout
  std::cout << table;
  
  return 0;
}

static int cutfunc(int argc, char** argv) {

  bool strict_cut = false;
  std::string cut; // list of markers, csv separated, to cut on
  
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

  if (die || in_out_process(argc, argv)) {
    
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
  cutp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  cutp.SetParams(tokens); 

  // process 
  if (!table.StreamTable(cutp, opt::infile)) {
    return 1;
  }

  return 0;
}

static int averagefunc(int argc, char** argv) {
 
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift average [csvfile] <options>\n"
      "  Calculate the average of each data column\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // set table params
  if (opt::verbose)
    table.SetVerbose();

  AverageProcessor avgp;
  avgp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);  

  if (table.StreamTable(avgp, opt::infile))
    return 1; // non-zero status on StreamTable

  // write the one line with the averages
  avgp.EmitCell();
  
  return 0;

  
}

static int log10func(int argc, char** argv)  {
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> opt::n; break;
    default: die = true;
    }
  }
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift log10 [csvfile] <options>\n"
      "  Calculate the log10 of marker intensities\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // set table params
  if (opt::verbose)
    table.SetVerbose();

  LogProcessor logp;
  logp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);  

  if (table.StreamTable(logp, opt::infile))
    return 1; // non-zero status on StreamTable

  return 0;
}



// parse the command line options
static void parseRunOptions(int argc, char** argv) {

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
	 opt::module == "count" || opt::module == "clean" ||
	 opt::module == "cat" || opt::module == "cereal" || 
	 opt::module == "correlate" || opt::module == "info" ||
	 opt::module == "cut" || opt::module == "view" ||
	 opt::module == "average" || opt::module == "lda" || 
	 opt::module == "spatial" || opt::module == "radialdens" || 
	 opt::module == "select" || opt::module == "pheno")) {
    std::cerr << "Module " << opt::module << " not implemented" << std::endl;
    die = true;
  }
  
  if (die) {
    std::cerr << "\n" << RUN_USAGE_MESSAGE;
    if (die)
      exit(EXIT_FAILURE);
    else 
      exit(EXIT_SUCCESS);	
  }
}

static int roifunc(int argc, char** argv) {

  std::string roifile;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> opt::n; break;
    case 'r' : arg >> roifile;
    default: die = true;
    }
  }
  
  if (die || roifile.empty() || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift roi [csvfile] <options>\n"
      "  Subset or label the cells to only those contained in the rois\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -r                        ROI file\n"
      "  -l                        Output all cells and add \"roi\" column with ROI label\n"      
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // read in the roi file
  std::vector<Polygon> rois = read_polygons_from_file(roifile);

  if (opt::verbose)
    for (const auto& c : rois)
      std::cerr << c << std::endl;

  ROIProcessor roip;
  roip.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  roip.SetParams(false, rois);// false is placeholder for label function, that i need to implement

  if (table.StreamTable(roip, opt::infile))
    return 1; // non-zero status in StreamTable

  return 0;
  
}

static int viewfunc(int argc, char** argv) {

  int precision = -1;

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
  
  if (die || in_only_process(argc, argv)) {
    
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
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> n_bins; break;
    case 'w' : arg >> w_bins; break;      
    default: die = true;
    }
  }

  if (die || in_only_process(argc, argv)) {
    
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
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'l' : arg >> length; break;
    case 'w' : arg >> width; break;      
    default: die = true;
    }
  }

  if (die || in_only_process(argc, argv)) {
    
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

  build_table();
  
  // make an ASCII plot of this
  table.PlotASCII(width, length);

  return 0;
}

static int subsamplefunc(int argc, char** argv) {

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> opt::n; break;
    case 's' : arg >> opt::seed; break;      
    default: die = true;
    }
  }

  // Process any remaining no-flag options
  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE =
      "Usage: cysift subsample [csvfile] <options>\n"
      "  Subsamples a cell quantification table, randomly.\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -n, --numrows             Number of rows to subsample\n"
      "    -s, --seed         [1337] Seed for random subsampling\n"
      "    -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();
  
  // subsample
  table.Subsample(opt::n, opt::seed);

  table.PrintTable(opt::header);
  
  return 0;
  
}

static int correlatefunc(int argc, char** argv) {

  bool sorted = false;
  bool csv_print = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'j' : csv_print = true; break;
    case 's' : sorted = true; break;
    default: die = true;
    }
  }

  if (die || in_only_process(argc, argv)) {
    
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

  build_table();
  
  //
  table.PrintPearson(csv_print, opt::sort);
  
  return 0;
}

static int cropfunc(int argc, char** argv) {

  std::string cropstring;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'c' : arg >> cropstring; break;
     default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift crop [csvfile] <options>\n"
      "  Crop the table to a given rectangle (in pixels)\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    --crop                    String of form xlo,xhi,ylo,yhi\n"
      "    -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  int xlo, xhi, ylo, yhi;
  std::vector<int*> coordinates = {&xlo, &xhi, &ylo, &yhi};

  std::istringstream iss(cropstring);
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

  build_table();
  
  //
  table.Crop(xlo, xhi, ylo, yhi);

  return 0;
  
}

static int spatialfunc(int argc, char** argv) {

  int n = 10;
  int d = -1;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> opt::threads; break;
    case 'k' : arg >> n; break;
    case 'd' : arg >> d; break;            
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
  
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

  build_table();

  // no table to work with
  if (!table.size())
    return 1;
  
  table.SetupOutputWriter(opt::outfile);
  
  table.KNN_spatial(n, d);

  //table.PrintTable(opt::header);

  return 0;
}

static int selectfunc(int argc, char** argv) {

  cy_uint logor = 0;
  cy_uint logand = 0;
  bool lognot = false;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'o' : arg >> logor; break;
    case 'a' : arg >> logand; break;
    case 'N' : lognot = true; break;      
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
  
    const char *USAGE_MESSAGE =
      "Usage: cysift select [csvfile]\n"
      "  Select cells by phenotype flag\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -o                    Logical OR flags\n"
      "    -a                    Logical AND flags\n"
      "    -N                    Not flag\n"
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // setup the selector processor
  SelectProcessor select;
  select.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  select.SetParams(logor, logand, lognot);

  // process
  table.StreamTable(select, opt::infile);
  
  return 0;
}

static int phenofunc(int argc, char** argv) {

  std::string file;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> file; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE =
      "Usage: cysift pheno [csvfile]\n"
      "  Phenotype cells (set the flags) with threshold file\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -t               File that holds gates: marker(string), low(float), high(float)\n"
      "    -v, --verbose    Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // read the phenotype filoe
  PhenoMap pheno = phenoread(file);

  if (!pheno.size()) {
    std::cerr << "Unable to read phenotype file or its empty: " << file << std::endl;
    return 2;
  }
  
  if (opt::verbose)
    for (const auto& c : pheno)
      std::cerr << c.first << " -- " << c.second.first << "," << c.second.second << std::endl;

  PhenoProcessor phenop;
  phenop.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  phenop.SetParams(pheno);

  if (table.StreamTable(phenop, opt::infile))
    return 1; // non-zero StreamTable status
  
  return 0;
}

static int radialdensfunc(int argc, char** argv) {

  cy_uint inner = 0;
  cy_uint outer = 20;
  cy_uint logor = 0;
  cy_uint logand = 0;
  std::string label;

  std::string file;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
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

  if (die || in_out_process(argc, argv)) {
    
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
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  if (inner >= outer) {
    std::cerr << "Inner radius should be smaller than outer, or else no cells are included" << std::endl;
    return 1;
  }

  // read in the multiple selection file
  std::vector<RadialSelector> rsv;
  if (!file.empty()) {
    std::ifstream input_file(file);
   
    if (!input_file.is_open()) {
      throw std::runtime_error("Failed to open file: " + file);
    }
    
    std::string line;
    std::regex pattern("^[-+]?[0-9]*\\.?[0-9]+,");
    while (std::getline(input_file, line)) {

      // remove white space
      line.erase(std::remove_if(line.begin(), line.end(), ::isspace),line.end());
	  
      // make sure starts with number
      if (!std::regex_search(line, pattern))
	continue;
     
      rsv.push_back(RadialSelector(line));
    }
    
    if (opt::verbose) {
      for (const auto& rr : rsv)
	std::cerr << rr << std::endl;
    }
  }

  /////// NEW WAY
  RadialProcessor radp;
  radp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  
  if (rsv.empty()) {
    radp.SetParams({inner},{outer},{logor},{logand},{label});
  } else {
    std::vector<cy_uint> innerV(rsv.size());
    std::vector<cy_uint> outerV(rsv.size());  
    std::vector<cy_uint> logorV(rsv.size());
    std::vector<cy_uint> logandV(rsv.size());  
    std::vector<std::string> labelV(rsv.size());
    for (size_t i = 0; i < innerV.size(); i++) {
      innerV[i] = rsv.at(i).int_data.at(0);
      outerV[i] = rsv.at(i).int_data.at(1);
      logorV[i] = rsv.at(i).int_data.at(2);
      logandV[i] = rsv.at(i).int_data.at(3);     
      labelV[i] = rsv.at(i).label;
    }
    radp.SetParams(innerV, outerV, logorV, logandV, labelV);
  }

  if (table.StreamTable(radp, opt::infile))
    return 1; // non-zero exit from StreamTable

  return 0;
}

int debugfunc(int argc, char** argv) {
  return 1;

}

static int cerealfunc(int argc, char** argv) {
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift cys [csvfile]\n"
      "  Create a .cys formatted file from a csv file\n" 
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  CerealProcessor cerp;
  cerp.SetParams(opt::outfile, cmd_input);

  table.StreamTableCSV(cerp, opt::infile);

  return 0;
}

static int countfunc(int argc, char** argv) {
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if (die || in_only_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift count [cysfile]\n"
      "  Output the number of cells in a file\n"
      "    cysfile: filepath or a '-' to stream to stdin\n"
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  CountProcessor countp;

  if (table.StreamTable(countp, opt::infile)) 
    return 1; // non-zero status on StreamTable
  
  countp.PrintCount();

  return 0;
}

// return TRUE if you want the process to die and print message
static bool in_out_process(int argc, char** argv) {
  
  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } else if (opt::outfile.empty()) {
      opt::outfile = argv[optind];      
    }
    optind++;
  }
  
  return opt::infile.empty() || opt::outfile.empty();

}

// return TRUE if you want the process to die and print message
static bool in_only_process(int argc, char** argv) {
  
  optind++;
  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    }
    optind++;
  }
  
  return opt::infile.empty();

}
