#include "cell_column.h"
#include "cell_table.h"

#include <unistd.h> // or #include <getopt.h> on Windows systems
#include <getopt.h>

namespace opt {
  static bool verbose = false;
  static std::string infile;
  static std::string quantfile;
  static std::string outfile;
  static std::string module;

  static bool header = false;
  
  static bool csv = false; // should we print as csv instead of screen readable
  
  static std::string roifile;

  static std::string cropstring;

  static std::string cut; // list of markers, csv separated, to cut on
  static bool strict_cut = false; // if true, don't carry ID or dims
  
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

static const char* shortopts = "jhyvr:t:g:b:q:c:s:n:r:w:l:x:X:";
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
"Modules:\n"
"  cut - Select only given markers and metas\n"
"  subsample - Subsample cells randomly\n"
"  plot - Generate an ASCII style plot\n"
"  roi - Trim cells to a region of interest within a given polygon\n"
"  histogram - Create a histogram of the data\n"
"  log10 - Apply a base-10 logarithm transformation to the data\n"
"  correlate - Calculate the correlation between variables\n"
"  info - Display information about the dataset\n"
"\n";

static int subsamplefunc();
static int infofunc();
static int roifunc();
static int correlatefunc();
static int histogramfunc();
static int plotfunc();
static int debugfunc();
static int cropfunc();
static int roundfunc();
static int subsamplefunc();
static int log10func();
static int cutfunc();

static void parseRunOptions(int argc, char** argv);

int main(int argc, char **argv) {
  
  // Check if a command line argument was provided
  if (argc < 2) {
    std::cerr << "Error: missing command line argument" << std::endl;
    std::cerr << RUN_USAGE_MESSAGE;
    return 1;
  }

  parseRunOptions(argc, argv);
  
  // get the module
  if (opt::module == "debug") {
    return(debugfunc());
  } else if (opt::module == "subsample") {
    return(subsamplefunc());
  } else if (opt::module == "plot") {
    return(plotfunc());
  } else if (opt::module == "round") {
    return(roundfunc());
  } else if (opt::module == "roi") {
    return(roifunc());
  } else if (opt::module == "crop") {
    return(cropfunc());
  } else if (opt::module == "correlate") {
    return(correlatefunc());
  } else if (opt::module == "histogram") {
    return(histogramfunc());
  } else if (opt::module == "log10") {
    return (log10func());
  } else if (opt::module == "cut") {
    cutfunc();
  } else if (opt::module == "info") {
    infofunc();
  } else {
    assert(false);
  }
  
  return 1;
}

static int infofunc() {
  
  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift info [csvfile]\n"
      "  Display basic information on the cell table\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -v, --verbose             Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  CellTable table(opt::quantfile.c_str(), false);

  std::cout << table;
  
  return 0;
}

static int cutfunc() {
  
  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift cut [csvfile] <marker1,markers2>\n"
      "  Cut the file to only certain markers\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -x, --cut                 Comma-separated list of markers to cut to\n"
      "  -v, --verbose             Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  CellTable table(opt::quantfile.c_str(), false);
  
  // parse the markers
  std::set<std::string> tokens;
  std::stringstream ss(opt::cut);
  std::string token;
  while (std::getline(ss, token, ',')) {
    tokens.insert(token);
  }

  // if not a strict cut, keep dims and id
  if (!opt::strict_cut) {
    tokens.insert(table.GetHeader().GetX());
    tokens.insert(table.GetHeader().GetY());
    tokens.insert(table.GetHeader().GetZ());
    tokens.insert(table.GetHeader().GetID());
  }
  tokens.erase(std::string());
  
  // cut from the table and header
  table.Cut(tokens);

  table.PrintTable(opt::header);
  
  return 0;
}

static int log10func()  {

  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift log10 [csvfile] <options>\n"
      "  Calculate the log10 of marker intensities\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  CellTable table(opt::quantfile.c_str(), false);

  table.Log10();

  table.PrintTable(opt::header);

  return 0;
}

// parse the command line options
static void parseRunOptions(int argc, char** argv) {
    bool die = false;

  if (argc <= 1) 
    die = true;

  bool help = false;
  std::stringstream ss;

  std::string tmp;
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
    case 'w' : arg >> opt::width; break;
    case 'l' : arg >> opt::height; break;
    case 'y' : opt::sort = true; break;
    default: die = true;
    }
  }

  // Process any remaining no-flag options
  while (optind < argc) {
    if (opt::module.empty()) {
      opt::module = argv[optind];      
    } else if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } else if (opt::outfile.empty()) {
      opt::outfile = argv[optind];
    } 
    optind++;
  }

  if (! (opt::module == "debug" || opt::module == "subsample" ||
	 opt::module == "plot"  || opt::module == "roi" ||
	 opt::module == "histogram" || opt::module == "log10" ||
	 opt::module == "crop"  || opt::module == "round" || 
	 opt::module == "correlate" || opt::module == "info" ||
	 opt::module == "cut")) {
    std::cerr << "Module " << opt::module << " not implemented" << std::endl;
    die = true;
  }
  
  if (die || help) 
    {
      std::cerr << "\n" << RUN_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

static int roifunc() {

  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift roi [csvfile] <options>\n"
      "  Subset the cells to only those contained in the rois\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -r                        ROI file\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  // read in the roi file
  std::vector<Polygon> rois = read_polygons_from_file(opt::roifile);

  // read in the file
  CellTable table(opt::quantfile.c_str(), opt::verbose);  

  // subset the vertices
  table.SubsetROI(rois);

  // write to stdout
  table.PrintTable(opt::header);
  
  return 0;
  
}

static int histogramfunc() {

  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift histogram [csvfile] <options>\n"
      "  Calculate the histogram of a set of markers\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -n  [50]                  Number of bins\n"
      "  -w  [50]                  Binwidth\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  // default is to use numbins. If binwidth is specified,
  // then that takes precedence
  if (opt::n == 0)
    opt::n = 50;

  if (opt::width != 50)
    opt::n = 0;

  // read in the file
  CellTable table(opt::quantfile.c_str(), opt::verbose);    
  
  //table.histogram(opt::n, opt::width);

  return 0;
  
}

static int plotfunc() {

  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift plot [csvfile] <options>\n"
      "  Outputs an ASCII-style plot of cell locations\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -l, --length        [50]  Height (length) of output plot, in characters\n"   
      "  -w, --width         [50]  Width of output plot, in characters\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  // read in the file
  CellTable table(opt::quantfile.c_str(), opt::verbose);    
  
  // make an ASCII plot of this
  table.PlotASCII(opt::width, opt::height);

  return 0;
}

static int roundfunc() {

  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift round [csvfile] <options>\n"
      "  Rounds all of the numerics to <-n> decimals\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -n  [0]                   Number of decimals to keep\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  // read in the file
  CellTable table(opt::quantfile.c_str(), opt::verbose);      

  table.SetPrecision(opt::n);

  std::cerr << " priting" << std::endl;
  
  table.PrintTable(opt::header);
  
  return 0;
  
}

static int subsamplefunc() {

  // display help if no input
  if (opt::infile.empty()) {
    
    const char *SUBSAMPLE_USAGE_MESSAGE =
      "Usage: cysift subsample [csvfile] <options>\n"
      "  Subsamples a cell quantification table, randomly.\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -n, --numrows             Number of rows to subsample\n"
      "  -s, --seed         [1337] Seed for random subsampling\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << SUBSAMPLE_USAGE_MESSAGE;
    return 0;
  }

  // read in the file
  CellTable table(opt::quantfile.c_str(), opt::verbose);      

  // subsample
  table.Subsample(opt::n, opt::seed);

  // output
  table.PrintTable(opt::header);
  
  return 0;
  
}

static int correlatefunc() {

  // display help if no input
  if (opt::infile.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift correlate [csvfile] <options>\n"
      "  Outputs an ASCII-style plot of marker intensity pearson correlations\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  -v, --verbose             Increase output to stderr\n"
      "  -j                        Output as a csv file\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
  }

  // read in the file
  CellTable table(opt::quantfile.c_str(), opt::verbose);      

  //
  table.PrintPearson(opt::csv, opt::sort);
  
  return 0;
}

static int cropfunc() {

  // display help if no input
  if (opt::infile.empty() || opt::cropstring.empty()) {
    
    const char *USAGE_MESSAGE =
      "Usage: cysift crop [csvfile] <options>\n"
      "  Crop the table to a given rectangle (in pixels)\n"
      "  csvfile: filepath or a '-' to stream to stdin\n"
      "  --crop                    String of form xlo,xhi,ylo,yhi\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 0;
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

  // read the table
  CellTable table(opt::quantfile.c_str(), opt::verbose);        

  //
  table.Crop(xlo, xhi, ylo, yhi);

  // write it out
  table.PrintTable(opt::header);
  
  return 0;
  
}

static int debugfunc() {

  CellTable table(opt::quantfile.c_str(), false);

  std::cerr << table << std::endl;
  
  table.PrintTable(opt::header);
  
  return 1;
    
  std::cerr << "...subset roi" << std::endl;
  std::vector<std::pair<float,float>> roi(
{{21910.26,41337.85},
{20878.64,41712.79},
{20631.33,41524.77},
{20500.91,41412.35},
{20245.42,41305.65},
{20084.74,41347.66},
{19981.91,41051.02},
{19792.77,40859.37},
{19675.93,40768.16},
{19532.80,40742.42},
{19394.83,40805.17},
{19483.95,40630.09},
{19309.55,40565.94},
{19619.63,39964.34},
{19688.64,38782.54},
{19835.26,38234.93},
{20322.71,37796.40},
{19803.00,37498.01},
{19718.38,36439.55},
{20790.33,37171.06},
{20306.44,36209.11},
{29731.14,39150.87},
{29008.49,40268.42},
{28977.95,42623.04},
{28238.65,41010.24},
{27152.96,41107.28},
{27827.88,43090.00},
{26482.79,43447.49},
{26122.26,41982.23},
{25130.31,42846.87},
{23918.49,42554.43}});

  //table.SubsetROI(roi);
  //std::cerr << "subsetted ROI has size " << table.num_cells() << std::endl;

  // make an ASCII plot of this
  std::cerr << "...plotting" << std::endl;
  table.PlotASCII(50, 25);
  
  return 0;
  
}
