#include "cell_column.h"
#include "cell_table.h"
#include "cell_processor.h"
#include "color_map.h"

#include <unistd.h> // or #include <getopt.h> on Windows systems
#include <getopt.h>
#include <ctime>
#include <regex>

#include "cell_row.h"
#include "cell_selector.h"

#ifdef HAVE_TIFFLIB
#include "tiff_reader.h"
#include "tiff_writer.h"
#endif

static std::string cmd_input;
static bool die = false;

namespace opt {
  static bool verbose = false;

  // file i/o
  static std::string infile;
  static std::string outfile;
  static std::string module;
  static StringVec infile_vec;
 
  static bool header = false;
  static bool header_only = false;

  static int width = 50;
  static int height = 50;
  
  static int threads = 1;

  static int seed = 1337;
  static int n = 0;

  static bool sort = false;
}

// process in and outfile cmd arguments
static bool in_out_process(int argc, char** argv);
static bool out_only_process(int argc, char** argv);
static bool in_only_process(int argc, char** argv);

static CellTable table;

static void build_table();

//static const char* ushortopts = "CJjhHNyvmMPr:e:g:G:p:t:a:i:A:O:d:b:c:s:k:n:r:w:l:L:x:X:o:R:f:D:V:z:S:";
static const struct option longopts[] = {
  { NULL, 0, NULL, 0 }
};


static const char *RUN_USAGE_MESSAGE =
"Usage: cyftools [module] <options> \n"
"VERSION: 0.1.0\n"
" --- Information --- \n"
"  view        - View the cell table as character data\n"
"  info        - Display detailed information\n"
"  summary     - Display brief information\n"    
"  count       - Count cells\n"
  //"  head        - Returns first lines of a file\n"
" --- Low-level processing ---\n"
"  convert     - Create a .cyf format file from a CSV\n"    
"  cut         - Select only given markers and metas\n"
"  reheader    - Change the header\n"
"  clean       - Remove classes of data (e.g. all meta)\n"
"  cat         - Concatenate multiple files\n"
"  sort        - Sort the cells\n"
"  subsample   - Subsample cells randomly\n"
"  sampleselect- Select a particular sample from a multi-sample file\n"
"  check       - Simply read and stream without any other operations, format check and start pipelines\n"
" --- Spatial transformations ---\n"  
"  offset      - Offset the x-y positions of the cells\n"
"  flip        - Flip x and/or y positions\n"
"  magnify     - Rescale the x and y coordinates\n"    
" --- Marker ops ---\n"  
"  pheno       - Phenotype cells (set the phenotype flags)\n"
"  filter      - Mark cellsfor later processing (filtering etc)\n"  
"  cellcount   - Provide total slide cell count for each cell type\n"
"  roi         - Trim cells to a region of interest within a given polygon\n"  
"  rescale     - Rescale the marker intensities similar to scimap\n"
" --- Numeric ---\n"
"  mean        - Collapse to mean of each column\n"  
"  divide      - Divide two columns\n"
"  log10       - Apply a base-10 logarithm transformation to the data\n"
"  jaccard     - Calculate the Jaccard similarty coeffiecient for cell flags\n"
"  pearson     - Calculate the Pearson correlation between marker intensities\n"
"  dist        - Calculate cell-cell distances\n"  
" --- Clustering ---\n"
"  dbscan      - Cluster cells based on the DBSCAN method\n"
"  tls         - Identify TLS structures\n"
" --- PDF/PNG ---\n"
"  png         - Plot PNG\n"
  //"  plot       - Generate an ASCII style plot\n"
  //"  histogram  - Create a histogram of the data\n"
" --- Graph ops ---\n"
"  delaunay    - Calculate the Delaunay triangulation\n"
"  umap        - Construct the marker space UMAP\n"
  //"  spatial     - Construct the spatial KNN graph\n"
"  tumor       - Set the tumor flag using KNN approach\n"
"  margin      - Set the tumor margin using KNN approach\n"
"  island      - Remove islands of stroma within tumor\n"  
"  radialdens  - Calculate density of cells within a radius\n"  
" --- Convolution ---\n"
"  convolve    - Density convolution to produce TIFF\n"
" --- Neighborhood ---\n"
"  ldacreate   - Create a Latent Dirichlet model\n"
"  ldarun      - Run a Latent Dirichlet model\n"
" --- Null models ---\n"  
"  scramble    - Randomly permute the phenotype flags\n"
"  scatter     - Randomly assign x and y positions throughout the slide\n"
"  hallucinate - Randomly assign cell phenotypes\n"
"  synth       - Various approaches for creating synthetic data\n"  
"\n";

static int checkfunc(int argc, char** argv);
static int flipfunc(int argc, char** argv);
static int offsetfunc(int argc, char** argv);
static int forprintfunc(int argc, char** argv);
static int distfunc(int argc, char** argv);
static int tlsfunc(int argc, char** argv);
static int markcheckfunc(int argc, char** argv);
static int magnifyfunc(int argc, char** argv);
static int rescalefunc(int argc, char** argv);
static int islandfunc(int argc, char** argv);
static int marginfunc(int argc, char** argv);
static int trimfunc(int argc, char** argv);
static int syntheticfunc(int argc, char** argv);
static int dbscanfunc(int argc, char** argv);
static int sampleselectfunc(int argc, char** argv);
static int cellcountfunc(int argc, char** argv);
static int jaccardfunc(int argc, char** argv);
static int summaryfunc(int argc, char** argv);
static int hallucinatefunc(int argc, char** argv);
static int scramblefunc(int argc, char** argv);
static int scatterfunc(int argc, char** argv);
static int plotpngfunc(int argc, char** argv);
static int reheaderfunc(int argc, char** argv);
static int dividefunc(int argc, char** argv);
static int sortfunc(int argc, char** argv);
static int headfunc(int argc, char** argv);
static int convolvefunc(int argc, char** argv);
static int delaunayfunc(int argc, char** argv);
static int tumorfunc(int argc, char** argv);
static int ldacreatefunc(int argc, char** argv);
static int ldarunfunc(int argc, char** argv);
static int meanfunc(int argc, char** argv);
static int cleanfunc(int argc, char** argv);
static int countfunc(int argc, char** argv);
static int convertfunc(int argc, char** argv);
static int catfunc(int argc, char** argv);
static int radialdensfunc(int argc, char** argv);
static int subsamplefunc(int argc, char** argv);
static int viewfunc(int argc, char** argv);
static int infofunc(int argc, char** argv);
static int roifunc(int argc, char** argv);
static int pearsonfunc(int argc, char** argv);
static int histogramfunc(int argc, char** argv);
static int plotfunc(int argc, char** argv);
static int debugfunc(int argc, char** argv);
static int cropfunc(int argc, char** argv);
static int subsamplefunc(int argc, char** argv);
static int log10func(int argc, char** argv);
static int cutfunc(int argc, char** argv);
static int umapfunc(int argc, char** argv);
static int radiusfunc(int argc, char** argv);
static int filterfunc(int argc, char** argv);
static int flagsetfunc(int argc, char** argv); 
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
  } else if (opt::module == "tls") {
    val = tlsfunc(argc, argv);
  } else if (opt::module == "flagset") {
    val = flagsetfunc(argc, argv);
  } else if (opt::module == "ldacreate") {
    val = ldacreatefunc(argc, argv);
  } else if (opt::module == "check") {
    val = checkfunc(argc, argv);
  } else if (opt::module == "markcheck") {
    val = markcheckfunc(argc, argv);
  } else if (opt::module == "island") {
    val = islandfunc(argc, argv);
  } else if (opt::module == "ldarun") {
    val = ldarunfunc(argc, argv);
  } else if (opt::module == "clean") {
    val = cleanfunc(argc, argv);
  } else if (opt::module == "png") {
    val = plotpngfunc(argc, argv);
  } else if (opt::module == "radialdens") {
    val = radialdensfunc(argc, argv);
  } else if (opt::module == "cellcount") {
    val = cellcountfunc(argc, argv);
  } else if (opt::module == "margin") {
    val = marginfunc(argc, argv);
  } else if (opt::module == "for") {
    val = forprintfunc(argc, argv);
  } else if (opt::module == "dist") {
    val = distfunc(argc, argv);
  } else if (opt::module == "subsample") {
    val = subsamplefunc(argc, argv);
  } else if (opt::module == "plot") {
    return(plotfunc(argc, argv));
  } else if (opt::module == "divide") {
    return (dividefunc(argc, argv));
  } else if (opt::module == "offset") {
    return (offsetfunc(argc, argv));
  } else if (opt::module == "roi") {
    val = roifunc(argc, argv);
  } else if (opt::module == "crop") {
    val = cropfunc(argc, argv);
  } else if (opt::module == "pearson") {
    return(pearsonfunc(argc, argv));
  } else if (opt::module == "jaccard") {
    return(jaccardfunc(argc, argv));    
  } else if (opt::module == "histogram") {
    return(histogramfunc(argc, argv));
  } else if (opt::module == "log10") {
    return(log10func(argc, argv));
  } else if (opt::module == "annotate") {
    return(tumorfunc(argc, argv));
  } else if (opt::module == "cut") {
    val = cutfunc(argc, argv);
  } else if (opt::module == "rescale") {
    val = rescalefunc(argc, argv);
  } else if (opt::module == "magnify") {
    val = magnifyfunc(argc, argv);
  } else if (opt::module == "convert") {
    val = convertfunc(argc, argv);
  } else if (opt::module == "mean") {
    val = meanfunc(argc, argv);
  } else if (opt::module == "info") {
    return(infofunc(argc, argv));
  } else if (opt::module == "view") {
    return(viewfunc(argc, argv));
  } else if (opt::module == "cat") {
    return (catfunc(argc, argv));
  } else if (opt::module == "umap") {
    val = umapfunc(argc, argv);
  } else if (opt::module == "filter") {
    val = filterfunc(argc, argv);
  } else if (opt::module == "flip") {
    val = flipfunc(argc, argv);
  } else if (opt::module == "sampleselect") {
    val = sampleselectfunc(argc, argv);
  } else if (opt::module == "convolve") {
    val = convolvefunc(argc, argv);
    //  } else if (opt::module == "head") {
    //  val = headfunc(argc, argv);
  } else if (opt::module == "sort") {
    val = sortfunc(argc, argv);
  } else if (opt::module == "delaunay") {
    val = delaunayfunc(argc, argv);
  } else if (opt::module == "reheader") {
    val = reheaderfunc(argc, argv);
  } else if (opt::module == "synth") {
    val = syntheticfunc(argc, argv);
  } else if (opt::module == "dbscan") {
    val = dbscanfunc(argc, argv);
  } else if (opt::module == "pheno") {
    val = phenofunc(argc, argv);
  } else if (opt::module == "summary") {
    val = summaryfunc(argc, argv);
  } else if (opt::module == "count") {
    countfunc(argc, argv);
  } else if (opt::module == "scramble") {
    scramblefunc(argc, argv);
  } else if (opt::module == "scatter") {
    scatterfunc(argc, argv);
  // } else if (opt::module == "trim") {
  //   trimfunc(argc, argv);
  } else if (opt::module == "hallucinate") {
    hallucinatefunc(argc, argv);
  } else {
    assert(false);
  }

  return 0;
}


// build the table into memory
static void build_table() {

  // set table params
  table.setVerbose(opt::verbose);
  table.setThreads(opt::threads);
  
  // stream into memory
  BuildProcessor buildp;
  buildp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  table.StreamTable(buildp, opt::infile);

  // have to set the PG tag here, because BuildProcessor is a reader only
  table.setCmd(cmd_input);

  // check we were able to read the table
  if (table.CellCount() == 0) {
    std::cerr << "Warning: Table with no cells? Error in upstream operation?" << std::endl;
  }
  
}

static int magnifyfunc(int argc, char** argv) {
  const char* shortopts = "vf:";

  float factor = 1;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;      
    case 'f' : arg >> factor; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools magnify [.cyf file]\n"
      "  Rescale the x and y coordinates\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -f <float>     Scale factor (default 1)\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  MagnifyProcessor res;
  res.SetCommonParams(opt::outfile, cmd_input, opt::verbose); 
  res.SetParams(factor);
  
  if (table.StreamTable(res, opt::infile)) 
    return 1; // non-zero status on StreamTable
  return 0;
}

static int checkfunc(int argc, char** argv) {

  const char* shortopts = "v"; 
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
      std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true;
    default: die = true;
    }
  }
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools check <input.cyf> <output.cyf file>\n"
      "  Simply read and stream cyf. Useful to validate and to start pipelines that have variable first \"real\" module\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input file path or '-' to stream from stdin.\n"
      "  <output.cyf file>     Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Options:\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools check input.cyf - | cyftools othercmd\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // anything that is as "stream" function wher ethere is no dependency between cells, gets a "Processor" objects
  CheckProcessor checkp;
  checkp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);

  // run the processor and return if it fails
  if (!table.StreamTable(checkp, opt::infile)) {
    return 1;
  }

  return 0;

}

static int flipfunc(int argc, char** argv) {

  // parse the command line options
  int DEFAULT=-100000;
  int xflip = DEFAULT;
  int yflip = DEFAULT;
  int xmax = DEFAULT;
  int ymax = DEFAULT;
  const char* shortopts = "vx:y:X:Y:"; // : means theres and argument. No colon = flag
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
      std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true;
    case 'x' : arg >> xflip; break;
    case 'y' : arg >> yflip; break;
    case 'X' : arg >> xmax; break;
    case 'Y' : arg >> ymax; break;      
    default: die = true;
    }
  }

  if (xflip == DEFAULT && yflip == DEFAULT) {
    std::cerr << "********************************" << std::endl;
    std::cerr << " cyftools flip - need to specify a flip axis (x and or y)" << std::endl;
    std::cerr << "********************************" << std::endl;
    die = true;
  }

  if ( (xflip != DEFAULT && xmax == DEFAULT) || (yflip != DEFAULT && ymax == DEFAULT)) {
    std::cerr << "********************************" << std::endl;
    std::cerr << " cyftools flip - need to specify an X or Y size if flipping on that axis" << std::endl;
    std::cerr << "********************************" << std::endl;
    die = true;
  }
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools flip <input.cyf> <output.cyf file> [options]\n"
      "  Flip the x or y position (reflection operation)\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input file path or '-' to stream from stdin.\n"
      "  <output.cyf file>          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Options:\n"
      "  -x                        Position of a X-flip line\n"
      "  -y                        Position of a Y-flip line\n"
      "  -X                        X dimension\n"
      "  -Y                        Y dimension\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools flip input.cyf cleaned_output.cyf -x 100 -X 14032\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // anything that is as "stream" function wher ethere is no dependency between cells, gets a "Processor" objects
  FlipProcessor flipp;
  flipp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  flipp.SetParams(xflip, yflip, xmax, ymax); 
  
  // run the processor and return if it fails
  if (!table.StreamTable(flipp, opt::infile)) {
    return 1;
  }

  return 0;

}

// hidden function for easy looping in bash
static int forprintfunc(int argc, char** argv) {

  bool xargs = false;
  const char* shortopts = "xp";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
      std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'x' : xargs=true; break;
    default: die = true;
    }
  }

  if (xargs) {
    std::cout << "find . -name \"*.cyf\" " <<
      "| xargs -I{} -P 4 bash -c " <<
      "'./script.sh {}'" << std::endl;
    std::cout << std::endl;
    std::cout <<
      "#!/bin/bash" << std::endl <<
      "### place this as a script.sh" << std::endl <<
      "### remember you will have to chmod+x script.sh" << std::endl << 
      "filename=$(basename $1)" << std::endl <<
      "base_name=$(echo $filename | cut -d \".\" -f 1)" << std::endl <<
      "echo \"...working on $base_name\"" << std::endl <<
      "cyftools cmd $1 ${base_name}.cyf" << std::endl;
  } else {
    std::cout << "for file in *.cyf; do" <<
      " filename=$(basename $file); " << 
      " base_name=$(echo $filename | cut -d \".\" -f 1); " <<
      " echo \"...working on $base_name\"; " <<
      "cyftools cmd $file ${base_name}.out.cyf; done" << std::endl;
  }
  
  return 0;
}

// hidden function for debugging
static int markcheckfunc(int argc, char** argv) {

  if (in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = "hidden function to check if MARK_FLAG is set\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  
  MarkCheckProcessor res;
  res.SetCommonParams(opt::outfile, cmd_input, opt::verbose); 
  
  if (table.StreamTable(res, opt::infile)) 
    return 1; // non-zero status on StreamTable
  return 0;

  
}

static int distfunc(int argc, char** argv) {

  std::string id;
  const char* shortopts = "vi:";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'i' : arg >> id; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools dist <in.cyf> <out.cyf>\n"
      "  Calculate the cell-to-cell distance, based on marked cells from cytools filter\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "Requred Input:\n"
      "     -i <string>               ID to mark the output column (will be dist_<id>)\n"
      "Optional Options:\n"
      "     -v, --verbose             Increase output to stderr.\n"
      "Example:\n"
      " cyftools filter -a 1024 -M <in> - | cyftools dist - <out>\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // build the table into memory
  build_table();

  table.SetupOutputWriter(opt::outfile);

  // get the distance
  table.Distances(id);

  // write it out
  table.OutputTable();
  
  return 0;
}

static int tlsfunc(int argc, char** argv) {

  cy_uint bcell_marker = 0;
  cy_uint immune_marker = 0;  
  int min_cluster_size = 100;
  int dist_max = 100;
  
  const char* shortopts = "vm:d:b:i:";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'b' : arg >> bcell_marker; break;
    case 'i' : arg >> immune_marker; break;      
    case 'd' : arg >> dist_max; break;
    case 'm' : arg >> min_cluster_size; break;
    default: die = true;
    }
  }

  if (bcell_marker == 0 || immune_marker == 0) {
    std::cerr << "Error: cyftools tls - need to specify -b (Bcell marker, as base-10 number) and immune-cell marker(s) (as base-10-number)" << std::endl;
    die = true;
  }
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools tls <in.cyf> <out.cyf>\n"
      "  Identify tertiary lymphoid structures and enumarte in tls_id column\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "Requred Input:\n"
      "     -b [int]                  B-cell bit marker (e.g. 64)\n"
      "     -i [int]                  Immune-cells bits marker (e.g. 382032)\n"
      "Optional Options:\n"
      "     -v, --verbose             Increase output to stderr.\n"
      "     -m [100]                  Minimum number of cells to call a cluster\n"
      "     -d [100]                  Distance threshold (microns) to consider in KNN calcs\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // build the table into memory
  build_table();

  table.SetupOutputWriter(opt::outfile);
  
  table.CallTLS(bcell_marker, immune_marker, min_cluster_size,
		dist_max); 

  table.OutputTable();
  
  return 0;
}

static int rescalefunc(int argc, char** argv) {
  
  const char* shortopts = "v";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools rescale [.cyf file]\n"
      "  Rescale the marker intensities similar to scimap\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  RescaleProcessor res;
  res.SetCommonParams(opt::outfile, cmd_input, opt::verbose); 
  //res.SetParams(factor);
  
  if (table.StreamTable(res, opt::infile)) 
    return 1; // non-zero status on StreamTable
  return 0;
}

static int flagsetfunc(int argc, char** argv) {
  
  const char* shortopts = "vp:P:c:C:";
  cy_uint pflag_set = 0;
  cy_uint cflag_set = 0;
  cy_uint pflag_clear = 0;
  cy_uint cflag_clear = 0;  
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'p' : arg >> pflag_set; break;
    case 'P' : arg >> pflag_clear; break;
    case 'c' : arg >> cflag_set; break;
    case 'C' : arg >> cflag_clear; break;      
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv) || !(pflag_set + cflag_set + pflag_clear + cflag_clear)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools flagset <in> <out>\n"
      "  Set flags on or off based on filtering from cyftools filter\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "  -p <int>        P-flag to set if mark is on\n"
      "  -P <int>        P-flag to clear if mark is on\n"
      "  -c <int>        C-flag to set if mark is on\n"
      "  -C <int>        C-flag to clear if mark is on\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  FlagsetProcessor fs;
  fs.SetCommonParams(opt::outfile, cmd_input, opt::verbose); 
  fs.SetParams(pflag_set, cflag_set, pflag_clear, cflag_clear);
  
  if (table.StreamTable(fs, opt::infile)) 
    return 1; // non-zero status on StreamTable
  return 0;
}


static int reheaderfunc(int argc, char** argv) {

  const char* shortopts = "vr:s:";
  
  StringVec rename;
  StringVec sample;
  std::string str;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'r' : arg >> str; rename.push_back(str); break;
    case 's' : arg >> str; sample.push_back(str); break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools reheader [.cyf file]\n"
      "  Change or reheader the header only \n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -r       Rename a marker of meta column (old:new) \n"
      "    -s       Add a sample tag (<id>:<field> e.g. microns_to_pixel:0.325) \n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  /////////// read the header
  // std::istream *inputStream = nullptr;
  // std::unique_ptr<std::ifstream> fileStream;
  // // set input from file or stdin
  // if (opt::infile == "-") {
  //   inputStream = &std::cin;
  // } else {
  //   fileStream = std::make_unique<std::ifstream>(opt::infile, std::ios::binary);
  //   if (!fileStream->good()) {
  //     std::cerr << "Error opening: " << opt::infile << " - file may not exist" << std::endl;
  //     return 1;
  //   }
  //   inputStream = fileStream.get();
  // }

  // cereal::PortableBinaryInputArchive inputArchive(*inputStream);

  // CellHeader header;
  
  // // First read the CellHeader
  // try {
  //   inputArchive(header);
  // } catch (const std::bad_alloc& e) {
  //   // Handle bad_alloc exception
  //   std::cerr << "Memory allocation failed during deserialization: " << e.what() << std::endl;
  //   return 1;  // or handle the error appropriately for your program
  // } catch (const cereal::Exception& e) {
  //   // Handle exception if any error occurs while deserializing header
  //   std::cerr << "Error while deserializing header: " << e.what() << std::endl;
  //   return 1;  // or handle the error appropriately for your program
  // }
  
  /*  is.seekg(0, is.end);
  int length = is.tellg() - is.tellg(); // calculate the remaining bytes in the file
  is.seekg(is.tellg(), is.beg); // go to the current position
  */

  return 0;
}

static int islandfunc(int argc, char** argv) {

  const char* shortopts = "vd:t:n:TS";
  float dist = 100;
  bool tumor_fill = false;
  bool stroma_fill = false;
  int n = 100;

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> n; break;
    case 'd' : arg >> dist; break;
    case 'T' : tumor_fill = true; break;
    case 'S' : stroma_fill = true; break;      
    case 't' : arg >> opt::threads; break;      
    default: die = true;
    }
  }

    if (die || in_out_process(argc, argv) || (tumor_fill == stroma_fill)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools island[.cyf file]\n"
      "  Fill-in islands with <= N cells\n"
      "\n"
      "Arguments:\n"
      "  [.cyf file]                 Input .cyf file path or '-' to stream from stdin.\n"
      "  -n [100]                  Minimum size for an island to remain a stroma island\n"
      "\n"
      "Optional Options:\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "  -T                        Fill in islands of stroma within tumors\n"
      "  -S                        Fill in islands of tumors within stroma\n"      
      "\n"
      "Example:\n"
      "  cyftools island -n 100 -T input.cyf output.cyf\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();

  table.SetupOutputWriter(opt::outfile);

  table.ClearCFlag(MARGIN_FLAG);
  
  // fills NON-MARK with MARK
  if (tumor_fill)
    table.IslandFill(n, TUMOR_FLAG, false, TUMOR_FLAG, false);
  else if (stroma_fill)
    table.IslandFill(n, TUMOR_FLAG, true, TUMOR_FLAG, true);
  
  table.OutputTable();
  
  return 0;

}

static int marginfunc(int argc, char** argv) {

  const char* shortopts = "vd:t:T:M:";
  float dist = 100;
  cy_uint tumor_flag = TUMOR_FLAG;
  cy_uint margin_flag = MARGIN_FLAG;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'd' : arg >> dist; break;
    case 't' : arg >> opt::threads; break;
    case 'T' : arg >> tumor_flag; break;
    case 'M' : arg >> margin_flag; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools margin [.cyf file]\n"
      "  Identify cells at the tumor / stroma interface margin\n"
      "\n"
      "Arguments:\n"
      "  [.cyf file]               Input .cyf file path or '-' to stream from stdin.\n"
      "  -d [100]                  Radius to look at the margin for\n"
      "  -T [1]                    Flag to label tumor as\n"
      "  -M [4]                    Flag to label margin as\n"      
      "\n"
      "Optional Options:\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools margin -d 200 input.cyf output.cyf\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();

  table.SetupOutputWriter(opt::outfile);
  
  table.TumorMargin(dist, tumor_flag, margin_flag); 
  
  table.OutputTable();
  
  return 0;
  
}

static int dbscanfunc(int argc, char** argv) {
  
  // dbscan params
  float epsilon = 100.0f;
  int min_size = 25;
  int min_cluster_size = 100;

  const char* shortopts = "ve:s:c:";
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'e' : arg >> epsilon; break;
    case 's' : arg >> min_size; break;
    case 'c' : arg >> min_cluster_size; break;            
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv) || min_size < 0) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools filter <in> - -M <filter flags> | cyftools dbscan <in> <out> \n"
      "  Cluster MARKED cells based on the DBSCAN algorithm (need to run filter -M first)\n"
      "  NB: Cluster assignment of 0 means not in a cluster. First cluster is integer 1\n"
      "\n"
      "Arguments:\n"
      "  [.cyf file]                 Input .cyf file path or '-' to stream from stdin.\n"
      "\n"
      "Optional Options:\n"
      "     -v, --verbose             Increase output to stderr.\n"
      "  -e [100]                Epsilon parameter for DBSCAN\n"
      "  -s [25]                 Min size parameter for DBSCAN\n"
      "  -c [100]                Min cluster size (otherwise mark as zero)\n"
      "\n"
      "Example:\n"
      "  cyftools filter -A 32 <in> - -M | dbscan - output.cyf\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();

  table.SetupOutputWriter(opt::outfile);
  
  table.clusterDBSCAN(epsilon, min_size, min_cluster_size);
  
  table.OutputTable();
  
  return 0;
  
}

static int syntheticfunc(int argc, char** argv) {

  string module;
  int width = 1000;
  int height = 1000;
  int num_clusters = 10;
  int num_points = 100;
  double sigma = 100;

  const char* shortopts = "vs:m:w:l:p:n:e:";
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 's' : arg >> opt::seed; break;
    case 'm' : arg >> module; break;
    case 'w' : arg >> width; break;
    case 'l' : arg >> height; break;
    case 'p' : arg >> num_points; break;
    case 'n' : arg >> num_clusters; break;
    case 'e' : arg >> sigma; break;
    default: die = true;
    }
  }

  if (die || out_only_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools synthetic [.cyf file]\n"
      "  Various modules for creating synthetic data\n"
      "\n"
      "Arguments:\n"
      "  [.cyf file]                 Output .cyf file path or '-' to stream to stdout.\n"
      "\n"
      "Optional Options:\n"
      "  -s <int>                  Random seed.\n"
      "  -w <int>                  Width of synthetic .cyf file\n"
      "  -l <int>                  Length of synthetic .cyf file\n"
      "  -p <int>                  Number of points per cluster\n"
      "  -n <int>                  Number of clusters\n"
      "  -e <float>                Sigma parameter on cluster distribution\n"            
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools synthetic -m cluster\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }


  CellSynth synth(width, height);

  synth.SetSeed(opt::seed);
  
  synth.Clusters(num_clusters, num_points, sigma, 1);

  synth.WriteTable(opt::outfile);

  return 0;
  
}

static int jaccardfunc(int argc, char** argv) {

  bool sorted = false;
  bool csv_print = false;
  bool subset_score = false;

  const char* shortopts = "vjsS";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'j' : csv_print = true; break;
    case 's' : sorted = true; break;
    case 'S' : subset_score = true; break;      
    default: die = true;
    }
  }

  if (die || in_only_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools jaccard [.cyf file]\n"
      "  Compute the Jaccard similarity score for the cell phenotype flags.\n"
      "\n"
      "Arguments:\n"
      "  [.cyf file]                 Input .cyf file path or '-' to stream from stdin.\n"
      "\n"
      "Optional Options:\n"
      "  -j                        Output as a CSV file.\n"
      "  -s                        Sort the output by Jaccard score.\n"
      "  -S                        Modified Jaccard such that divisand is min of size of A and B.\n"
      "                            (This effectively gives score of 1 if A is subset of B or vice versa).\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools jaccard input.cyf -j -s\n"
      "  cyftools jaccard - -v\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }
  
  build_table();

  table.PrintJaccardSimilarity(csv_print, sorted, subset_score);

  return 0;
}

static int cellcountfunc(int argc, char** argv) {

  const char* shortopts = "va:";
  std::vector<uint64_t> additional_flags;
  std::string tmp_string;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'a' : arg >> tmp_string; additional_flags.push_back(std::stoull(tmp_string)); break;
    default: die = true;
    }
  }
  
  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools cellcount [.cyf file]\n"
      "  Compute the number of cells per marker and output as one \"cell\".\n"
      "\n"
      "Arguments:\n"
      "  [.cyf file]                 Input .cyf file path or '-' to stream from stdin.\n"
      "  -a                        Additional AND flag combos to test\n"
      "\n"
      "Optional Options:\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools cellcount input.cyf\n"
      "  cyftools cellcount - -v\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }
  
  
  CellCountProcessor cellp;
  cellp.SetCommonParams(opt::outfile, cmd_input, opt::verbose); // really shouldn't need any of these
  cellp.SetParams(additional_flags);

  if (table.StreamTable(cellp, opt::infile)) 
    return 1; // non-zero status on StreamTable
  
  cellp.EmitCell();

  return 0;
}

static int scatterfunc(int argc, char** argv) {

  int width = 1000;
  int height = 1000;

  const char* shortopts = "vw:l:s:";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'w' : arg >> width; break;
    case 'l' : arg >> height; break;
    case 's' : arg >> opt::seed; break;      
    default: die = true;
    }
  }

  if (width <= 0 || height <= 0) {
    std::cerr << "Error: Width and height must be positive integer" << std::endl;
  }

  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools scatter <input.cyf> [output.cyf file]\n"
      "  Randomly assign cell x,y position and output to a .cyf file.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input .cyf file path or '-' to stream from stdin.\n"
      "  [output.cyf file]          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Required Options:\n"
      "  -w <int>                  Width of the \"slide\" to scatter on.\n"
      "  -l <int>                  Length of the \"slide\" to scatter on.\n"
      "\n"
      "Optional Options:\n"
      "  -s <int>                  Random seed.\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools scatter input.cyf output_scattered.cyf -w 100 -l 200\n"
      "  cyftools scatter input.cyf - -w 100 -l 200 -s 12345 -v\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  ScatterProcessor scatterp;
  scatterp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  scatterp.SetParams(width, height, opt::seed); 

  // process 
  if (!table.StreamTable(scatterp, opt::infile)) {
    return 1;
  }
  
  return 0;

  
}

static int hallucinatefunc(int argc, char** argv) {

  const char* shortopts = "vn:s:";
  int n_phenotypes = 10;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> n_phenotypes; break;
    case 's' : arg >> opt::seed; break;
    default: die = true;
    }
  }

#ifdef USE_64_BIT
  if (n_phenotypes <= 0 || n_phenotypes > 64) {
    const std::string plim = "64";
#else
  if (n_phenotypes < 0 || n_phenotypes > 32) {
    const std::string plim = "32";    
#endif
    std::cerr << "Error: Number of phenotypes must be > 0 and < " << plim << std::endl;
  }
  
  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools hallucinate <input.cyf> [output.cyf file]\n"
      "  Randomly assign cell phenotypes and output to a .cyf file.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input .cyf file path or '-' to stream from stdin.\n"
      "  [output.cyf file]          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Required Options:\n"
      "  -n <int>                  Number of cell types possible.\n"
      "\n"
      "Optional Options:\n"
      "  -s <int>                  Random seed.\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools hallucinate input.cyf output.cyf -n 5\n"
      "  cyftools hallucinate input.cyf - -n 5 -s 12345 -v\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  HallucinateProcessor hallp;
  hallp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  hallp.SetParams(n_phenotypes, opt::seed); 

  // process 
  if (!table.StreamTable(hallp, opt::infile)) {
    return 1;
  }
  
  return 0;
}


static int scramblefunc(int argc, char** argv) {

  const char* shortopts = "vPs:p";
  bool lock_flags = false;
  bool phenotype_only = false; // scramble only the phenotype flags, leaving cell flags intact
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'P' : lock_flags = true; break;
    case 'p' : phenotype_only = true; break;
    case 'v' : opt::verbose = true; break;
    case 's' : arg >> opt::seed; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools scramble <input.cyf> [output.cyf file]\n"
      "  Scrambles the phenotype and cell flags among cells. Essentially, cell frequencies stay the same\n"
      "  and slide morphology, but the labels on each cell are scrambled.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input .cyf file path or '-' to stream from stdin.\n"
      "  [output.cyf file]          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Optional Options:\n"
      "  -P                        Flag to lock phenotype flags (just permute which cells they go to).\n"
      "  -p                        Flag to scramble only the phenotype flags (keeping cell flags intact)\n"
      "  -s <int>                  Random seed.\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools scramble input.cyf output_scrambled.cyf\n"
      "  cyftools scramble input.cyf - -P -s 12345 -v\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();

  table.SetupOutputWriter(opt::outfile);

  table.ScramblePflag(opt::seed, lock_flags, phenotype_only);

  // print it
  table.OutputTable();
  
  return 0;
}

static int convolvefunc(int argc, char** argv) {

  const char* shortopts = "vi:d:t:w:";
  
  int width = 200;
  std::string intiff;
  float microns_per_pixel = 0;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'i' : arg >> intiff; break;
    case 'd' : arg >> microns_per_pixel; break;
    case 't' : arg >> opt::threads; break;      
    case 'v' : opt::verbose = true; break;
    case 'w' : arg >> width; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv) || microns_per_pixel <= 0) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools convolve <input.cyf> <output_tiff> [options]\n"
      "  Perform a convolution to produce a TIFF.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input .cyf file path or '-' to stream from stdin.\n"
      "  <output_tiff>             Path to the output TIFF file.\n"
      "\n"
      "Required Options:\n"
      "  -i <tiff_file>            Input TIFF file to set params for output.\n"
      "  -d <float>                Number of microns per pixel (e.g. 0.325).\n"
      "  -w <int>                  Width of the convolution box (in pixels). Default: 200.\n"
      "\n"
      "Optional Options:\n"
      "  -t <int>                  Number of threads. Default: 1.\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools convolve input.cyf output.tiff -i input_params.tiff -d 0.325 -w 200\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // build the table
  // but don't have to convert columns
  // since we don't use pre-existing Graph or Flags for this
  build_table();

  // check we were able to read the table
  if (table.CellCount() == 0) {
    std::cerr << "Ending with no cells? Error in upstream operation?" << std::endl;
    std::cerr << "no tiff created" << std::endl;
    return 0;
  }

#ifdef HAVE_TIFFLIB  
  // open the TIFFs
  TiffReader itif(intiff.c_str());
  int inwidth, inheight;
  TIFFGetField(itif.get(), TIFFTAG_IMAGEWIDTH,  &inwidth);
  TIFFGetField(itif.get(), TIFFTAG_IMAGELENGTH, &inheight);

  // set the output tiff
  TiffWriter otif(opt::outfile.c_str());

  TIFFSetField(otif.get(), TIFFTAG_IMAGEWIDTH, inwidth); 
  TIFFSetField(otif.get(), TIFFTAG_IMAGELENGTH, inheight);

  TIFFSetField(otif.get(), TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(otif.get(), TIFFTAG_BITSPERSAMPLE, 16);
  TIFFSetField(otif.get(), TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(otif.get(), TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(otif.get(), TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  
  // build the umap in marker-space
  table.Convolve(&otif, width, microns_per_pixel);
#else
  std::cerr << "Unable to proceed, need to include and link libtiff and have preprocessor directory -DHAVE_TIFFLIB" << std::endl;
#endif
  
  return 0;
  
}

static int umapfunc(int argc, char** argv) {
  const char* shortopts = "vt:D:k:w:l:";
  int n = 15;
  int width = 1000;
  int height = 1000;
  std::string pdf_file;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> opt::threads; break;
    case 'D' : arg >> pdf_file; break;      
    case 'k' : arg >> n; break;
    case 'w' : arg >> width; break;
    case 'l' : arg >> height; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools umap <input.cyf> <output.cyf file> [options]\n"
      "  Construct the UMAP (in marker space) and output to a .cyf file with added umap1 and umap2 columns.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input .cyf file path or '-' to stream from stdin.\n"
      "  <output.cyf file>          Output .cyf file path with umap1 and umap2 columns added.\n"
      "\n"
      "Options:\n"
      "  -k <int>                  Number of neighbors. Default: 15.\n"
      "  -t <int>                  Number of threads. Default: 1.\n"
      "  -D <file>                 Optional output to PDF.\n"
      "  -w <int>                  Width of the PDF. Default: 1000.\n"
      "  -l <int>                  Length of the PDF. Default: 1000.\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools umap input.cyf output_umap.cyf -k 15 -t 2 -D output.pdf -w 1200 -l 800\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // build the table
  // but don't have to convert columns
  // since we don't use pre-existing Graph or Flags for this
  build_table();

  table.SetupOutputWriter(opt::outfile);
  
  // build the umap in marker-space
  if (!pdf_file.empty() && table.HasColumn("umap1")) {
    std::cerr << "...umap already exists. Skipping build and outputing PDF." << std::endl;
    std::cerr << "...if you want to rebuild UMAP, call without -D flag and then try again" << std::endl;
  } else if (!table.HasColumn("umap1")) {
    table.UMAP(n);
  } else {
    std::cerr << "Warning: Overwriting prior umap" << std::endl;
    table.UMAP(n);    
  }

  // pdf
  if (!pdf_file.empty()) {
    if (opt::verbose)
      std::cerr << "...outputting UMAP pdf" << std::endl;
    table.UMAPPlot(pdf_file, width, height); 
  } else {
    
    // print it
    table.OutputTable();
  }
  
  return 0;
}

static int ldacreatefunc(int argc, char** argv) {

  const char* shortopts = "vo:s:n:r:i:t:";
#ifdef HAVE_LDAPLUSPLUS
  std::string model_out, fields;
  int n_topics = 10;
  int n_iterations = 10;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'o' : arg >> model_out; break;
    case 's' : arg >> opt::seed; break;
      //    case 'm' : arg >> model_in; break;
    case 'n' : arg >> n_topics; break;
    case 'r' : arg >> fields; break;
    case 'i' : arg >> n_iterations; break;
    case 't' : arg >> opt::threads; break;      
    default: die = true;
    }
  }

  if (model_out.empty()) {
    die = true;
    std::cerr << "Error: Must specify model output file with -o" << std::endl;
  }

  if (fields.empty()) {
    die = true;
    std::cerr << "Error: Must specify columns to score" << std::endl;
  }
  
  if (die || in_only_process(argc, argv)) {
    const char *USAGE_MESSAGE = 
      "Usage: cyftools ldacreate <input_file> [options]\n"
      "  Create a topic-model learning using Latent Dirichlet Allocation.\n"
      "\n"
      "Arguments:\n"
      "  <input_file>              Input file path or '-' to stream from stdin.\n"
      "\n"
      "Required Options:\n"
      "  -r <list>                 Comma-separated list of input data columms to LDA model.\n"
      "  -o <file>                 Output model file.\n"
      "\n"
      "Optional Options:\n"
      "  -n <int>                  Number of topics. Default: 10.\n"
      "  -i <int>                  Max number of iterations. Default: 10.\n"
      "  -t <int>                  Number of threads. Default: 1.\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools ldacreate input.cyf -r input1,input2,input3 -o model_output -n 15 -i 20\n"
      "  cyftools ldacreate - -r input1,input2 -o model_output -t 4 -v\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // parse the markers
  std::set<std::string> tokens;
  std::stringstream ss(fields);
  std::string token;
  while (std::getline(ss, token, ',')) {
    tokens.insert(token);
  }
  
  // stream into memory
  build_table();

  // error check the markers
  StringVec markers;
  for (const auto& s : tokens) {

    if (!table.HasColumn(s)) {
      std::cerr << "Error: Requested column " << s << " not in table" << std::endl;
      return 1;
    }
    
    if (opt::verbose)
      std::cerr << "...lda - column " << s << std::endl;
    markers.push_back(s);
  }

  // create the model 
  table.LDA_create_model(markers,
			 n_topics,
			 n_iterations,
			 opt::seed);

  // write the model
  table.LDA_write_model(model_out);
  
  return 0;

#else
  std::cerr << "Error: Unable to run LDA without including ldaplusplus header library to build." << std::endl;
  return 1;
#endif

}

static int ldarunfunc(int argc, char** argv) {

  const char* shortopts = "vi:l:c:D:";
#ifdef HAVE_LDAPLUSPLUS
  
  std::string model_in, pdf;
  int topic_highlight = 0;
  float cont_cutoff = 0.10f;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'i' : arg >> model_in; break;
    case 'l' : arg >> topic_highlight; break;
    case 'c' : arg >> cont_cutoff; break;
    case 'D' : arg >> pdf; break;
    default: die = true;
    }
  }

  if (model_in.empty()) {
    die = true;
    std::cerr << "Error: Must specify model input file with -i" << std::endl;
  }

  if (!check_readable(model_in)) {
    die = true;
    std::cerr << "Error: Model file " << model_in << " not readable/exists" << std::endl;
  }

  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools ldarun <input.cyf> <output.cyf file> [options]\n"
      "  Score cells using a pre-computed topic-model learning using Latent Dirichlet Allocation.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input file path or '-' to stream from stdin.\n"
      "  <output.cyf file>          Output .cyf file path with scored cells.\n"
      "\n"
      "Required Options:\n"
      "  -i <model_file>           Input model file (won't re-run LDA).\n"
      "\n"
      "Optional Options:\n"
      "  -D <pngfile>              Output a PNG plot of the topics.\n"
      "  -c <float>                Display cutoff (between 0 and 1) to only plot topics for a given cell when that topic exceeds this threshold.\n"
      "  -l <topic_num>            Highlight topic number in the PDF output (only used with -D).\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools ldarun input.cyf output_scored.cyf -i model_input\n"
      "  cyftools ldarun input.cyf output_scored.cyf -i model_input -D -l 3 -c 0.1 -v\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // stream into memory
  build_table();
  
  // loading the model
  table.LDA_load_model(model_in);

  // score the scells
  table.LDA_score_cells(pdf, topic_highlight, cont_cutoff);

  // print the table
  table.SetupOutputWriter(opt::outfile);
  table.OutputTable();
  
  return 0;

#else
  std::cerr << "Error: Unable to run LDA without including ldaplusplus header library to build." << std::endl;
  return 1;
#endif
  
}

static int plotpngfunc(int argc, char** argv) {
  const char* shortopts = "vf:m:r:t:p:";
  
#ifdef HAVE_CAIRO
  
  float scale_factor = 0.25f;
  std::string roifile;
  std::string module;
  std::string title;
  std::string palette;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'f' : arg >> scale_factor; break;
    case 'r' : arg >> roifile; break;
    case 't' : arg >> title; break;
      //case 'm' : arg >> module; break;
    case 'p' : arg >> palette; break;
    default: die = true;
    }
  }

  if (module.empty() && palette.empty()) {
    std::cerr << "ERROR: cyftools png -- need to input module with -m or palette with -p" << std::endl;
    die = true;
  }

  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools png <input.cyf> <output_png> -p palette\n"
      "  Plot the input as a PNG file.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>               Input file path or '-' to stream from stdin.\n"
      "  <output_png>              Output PNG file path.\n"
      "  -p <file>                 File containing a palette which is: cflag value, pflag value, R, G, B, alpha, label\n."
      "                            Lines are in order, so that cell that meets two criteria gets the one with lower line number.\n"
      "                            Example line: 0,1024,255,0,0,1,Tcell\n"
      "                            NB: A label of \"nolabel\" will exclude that label from the legend\n"
      "                            NB: Should include a pflag=0,cflag=0 line at end to label cells not belonging to above (e.g. gray)\n"
      "\n"
      "Options:\n"
      "  -r <roifile>              Plot rois\n"
      "  -f <float>                Fractional scale factor. Default: 0.25. (1 means each pixel is 1 x-unit; smaller values result in a smaller image)\n"
      "  -t <string>               Title to display\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools png input.cyf output.png -f 0.5 -m tls\n"
      "  cyftools png - output.png -f 0.75 -v -m tumor\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  ColorLabelVec colors;
  
  // read the palette file
  if (!palette.empty()) {
    std::ifstream input_file(palette);
    if (!input_file.is_open()) {
      throw std::runtime_error("Failed to open file: " + palette);
    }
    std::string line;
    std::regex pattern("^[-+]?[0-9]*\\.?[0-9]+,");

    // loop the lines
    while (std::getline(input_file, line)) {

      // remove white space
      //line.erase(std::remove_if(line.begin(), line.end(), ::isspace),line.end());

      // skip empty and comment
      if (line.empty()) continue;
      if (line.at(0) == '#') continue;

      // tokenize
      StringVec tokens = tokenize_comma_delimited<StringVec>(line);
      if (tokens.size() != 7) {
	std::cerr << "Error: cytools png -- palette line " << line << " does not have 7 elements (pflag, cflag,r,g,b,a,label)" << std::endl;
	return 1;
      }

      // build the color
      CellColor col;
      col.cflag   = std::stoi(tokens[0]);
      col.pflag   = std::stoi(tokens[1]);
      col.c.red   = std::stoi(tokens[2]);
      col.c.green = std::stoi(tokens[3]);
      col.c.blue  = std::stoi(tokens[4]);
      col.c.alpha = std::stof(tokens[5]);
      col.label   = tokens[6];

      colors.push_back(col);

      if (opt::verbose)
	std::cerr << ".... read " << col << std::endl;
    }
  }
  
  // stream into memory
  build_table();
  
  table.PlotPNG(opt::outfile, scale_factor, module, roifile, title, colors);
  
  return 0;

#else
  std::cerr << "Error: Unable to run PNG without including cairo library in build." << std::endl;
  return 1;
#endif
  
}

 static int offsetfunc(int argc, char** argv) {

   const char* shortopts = "vx:y:";
   float x = 0;
   float y = 0;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'x' : arg >> x; break;
    case 'y' : arg >> y; break;
    default: die = true;
    }
  }

  if (!x & !y && opt::verbose) {
    std::cerr << "******************************" << std::endl;
    std::cerr << "No offsets provided, see usage" << std::endl;
    std::cerr << "******************************" << std::endl;    
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE = 
      "Usage: cyftools offset <input.cyf> <output.cyf file> [options]\n"
      "  Offset (shift) the cells in x and y (to align with e.g. ROIs)\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input file path or '-' to stream from stdin.\n"
      "  <output.cyf file>          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Options:\n"
      "  -x                        X offset [0]\n"
      "  -y                        Y offset [0]\n"      
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools offset input.cyf cleaned_output.cyf -x 100 -y 100\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  OffsetProcessor offp;
  offp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  offp.SetParams(x, y); 
  
  // process 
  if (!table.StreamTable(offp, opt::infile)) {
    return 1;
  }

  return 0;

  
 }
 
static int cleanfunc(int argc, char** argv) {

  const char* shortopts = "vmMGRcpC:P:";
  bool clean_markers = false;
  bool clean_meta = false;
  bool clean_cflags = false;
  bool clean_pflags = false;
  bool clean_programs = false;
  cy_uint cflag_reset = -1;
  cy_uint pflag_reset = -1; 
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'c' : clean_cflags= true; break;
    case 'p' : clean_pflags= true; break;      
    case 'm' : clean_markers = true; break;
    case 'M' : clean_meta = true; break;
    case 'P' : arg >> pflag_reset; break;
    case 'C' : arg >> cflag_reset; break;      
    case 'G' : clean_programs = true; break;
    case 'R' : clean_programs = true; clean_meta = true; clean_pflags = true; clean_cflags = true; break;
    default: die = true;
    }
  }

  
  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools clean <input.cyf> <output.cyf file> [options]\n"
      "  Clean up the data to reduce disk space.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input file path or '-' to stream from stdin.\n"
      "  <output.cyf file>          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Options:\n"
      "  -c                        Clear all cell flags\n"
      "  -p                        Clear all phenotype flags\n"
      "  -C <flags>                Clear specific cell flags\n"
      "  -P <flags>                Clear the phenotype flags\n"
      "  -m                        Remove all marker data.\n"
      "  -M                        Remove all meta data.\n"
      "  -G                        Remove all @PG tags from header.\n"
      "  -R                        Remove all flag, meta and program @PG data (resets to just original cell table)\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools clean input.cyf cleaned_output.cyf -m -M\n"
      "  cyftools clean input.cyf - -R\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }
  
  CleanProcessor cleanp;
  cleanp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  cleanp.SetParams(clean_programs, clean_meta, clean_markers, clean_cflags, clean_pflags,
		   cflag_reset, pflag_reset);  

  // process 
  if (!table.StreamTable(cleanp, opt::infile)) {
    return 1;
  }

  return 0;
  
}


 
static int catfunc(int argc, char** argv) {

  const char* shortopts = "v";
  std::string samples;
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
      "Usage: cyftools cat <.cyf file1> <.cyf file2> ... [options]\n"
      "  Concatenate together multiple cell tables and stream to stdout in .cyf format.\n"
      "\n"
      "Arguments:\n"
      "  <.cyf file1> <.cyf file2> ...  Filepaths of cell tables to concatenate.\n"
      "\n"
      "Options:\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools cat table1.cyf table2.cyf > tablecat.cyf\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // parse the sample numbers
  /*std::vector<int> sample_nums;
  std::unordered_set<std::string> tokens;
  std::stringstream ss(samples);
  std::string token;
  while (std::getline(ss, token, ',')) {
    sample_nums.push_back(std::stoi(token));
  }

  if (sample_nums.size() != opt::infile_vec.size() && sample_nums.size()) {
    throw std::runtime_error("Sample number csv line should have same number of tokens as number of input files");
    }*/
  
  size_t offset = 0;

  // force output to stdout
  opt::outfile = "-";

  // check files are readable at least
  for (const auto& v : opt::infile_vec) {
    if (!check_readable(v)) {
      std::cerr << "Error: " << v << " not readable/exists" << std::endl;
      return 1;
    }
  }
  
  CatProcessor catp;
  catp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  catp.SetParams(offset, 0);
  
  for (size_t sample_num = 0; sample_num < opt::infile_vec.size(); sample_num++) {

    // can make a new table for each iteration, since we are dumping right to stdout
    CellTable this_table;
    
    // set table params
    this_table.setVerbose(opt::verbose);

    if (opt::verbose)
      std::cerr << "...reading and concatenating " << opt::infile_vec.at(sample_num) <<
	" of " << opt::infile_vec.size() << " files " << std::endl;
    
    // update the offset and sample num
    catp.SetOffset(offset);
    catp.SetSample(sample_num);

    // stream in the lines
    this_table.StreamTable(catp, opt::infile_vec.at(sample_num));
    
    offset = catp.GetMaxCellID() + 1; // + 1 to avoid dupes if new cellid starts at 0
  }

  return 0;
}

static int summaryfunc(int argc, char** argv) {

  // display help if no input
  if (die || in_only_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools summary <.cyf file> [options]\n"
      "  Display a brief summary of the cell table and print to stdout.\n"
      "\n"
      "Arguments:\n"
      "  <.cyf file>                 Input .cyf file path or '-' to stream from stdin.\n"
      "\n"
      "Example:\n"
      "  cyftools summary table1.cyf\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  SummaryProcessor summp;
  summp.SetCommonParams(opt::outfile, cmd_input, opt::verbose); // really shouldn't need any of these

  if (table.StreamTable(summp, opt::infile)) 
    return 1; // non-zero status on StreamTable
  
  summp.Print();
  
  return 0;
}


static int infofunc(int argc, char** argv) {

  // display help if no input
  if (die || in_only_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools info <.cyf file> [options]\n"
      "  Display detailed information on the cell table and print to stdout.\n"
      "\n"
      "Arguments:\n"
      "  <.cyf file>                 Input .cyf file path or '-' to stream from stdin.\n"
      "\n"
      "\n"
      "Example:\n"
      "  cyftools info table1.cyf\n";
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

  const char* shortopts = "vx:";
  std::string cut; // list of markers, csv separated, to cut on
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'x' : arg >> cut; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools cut <input.cyf> <output.cyf file> -x <marker1,marker2,...> [options]\n"
      "  Cut the file to only certain meta or markers and output to a .cyf file or stream to stdout.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input .cyf file path or '-' to stream from stdin.\n"
      "  <output.cyf file>          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Required Options:\n"
      "  -x <markers>              Comma-separated list of meta or markers to cut to.\n"
      "\n"
      "Optional Options:\n"
      "  -v                        Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools cut input.cyf output_cut.cyf -x marker1,meta1\n"
      "  cyftools cut input.cyf - -x marker1,marker2,meta1 -v\n";
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
  table.setVerbose(opt::verbose);

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

static int meanfunc(int argc, char** argv) {

  std::string group_by;
  const char* shortopts = "vb:";  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'b' : arg >> group_by; break;
    default: die = true;
    }
  }
  
  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools mean <input.cyf> <output.cyf> [options]\n"
      "  Calculate the mean of each data column and output to a .cyf file or stream to stdout.\n"
      "  The output will contain a single 'cell' with the means for each column.\n"
      "\n"
      "Arguments:\n"
      "  <input.cyf>           Input .cyf file path or '-' to stream from stdin.\n"
      "  <output.cyf file>          Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Options:\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "  -b <string>               Column name to group means by, so each unique elem of groupby column is own \"meta\" cell\n"
      "\n"
      "Example:\n"
      "  cyftools mean input.cyf output_mean.cyf -b tls_id\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // set table params
  table.setVerbose(opt::verbose);

  AverageProcessor avgp;
  avgp.SetParams(group_by);
  avgp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);

  if (table.StreamTable(avgp, opt::infile))
    return 1; // non-zero status on StreamTable

  // write the one line with the averages
  avgp.EmitCells();
  
  return 0;

}

static int tumorfunc(int argc, char** argv) {

  const char* shortopts = "vt:k:f:d:F:B";

  float dist = 100000;
  int n = 25;
  cy_uint flag_to_set = 0;
  float frac = 0.50;

  // when set, will only build KNN graph with BUILD_MARK + cells
  // not really used yet, so not exposed explicitly
  //bool build_tree_with_marked_only = false; 
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> opt::threads; break;
    case 'k' : arg >> n; break;
    case 'f' : arg >> frac; break;
    case 'F' : arg >> flag_to_set; break;
    case 'd' : arg >> dist; break;
      //case 'B' : build_tree_with_marked_only = true; break;
    default: die = true;
    }
  }

  if (flag_to_set != 1 && flag_to_set != 16 && flag_to_set != 32) {
    std::cerr << "-F (c flag to set) must be 1,16,32" << std::endl;
    die = true;
  }

  if (die || in_out_process(argc, argv)) {
  
    const char *USAGE_MESSAGE =
      "Usage: cyftools annotate [.cyf file]\n"
      "  Using KNN, set the flag on whether a cell is in the region, using marked cells\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -k [25]               Number of neighbors\n"
      "    -f [0.50]             Fraction of neighbors\n"
      "    -d [100000]           Max distance to consider\n"
      "    -F [1]                Cflag to set (1=tumor, 16=Tcell)\n"
      "    -v, --verbose         Increase output to stderr\n"
      "  Example:\n"
      "           # Check if a cell has >50% of 25 nearest neighbors with pflag 4 set\n"
      "           cyftools filter -a 4 <in.cyf> - -M | cyftools annotate -f 0.5 -k 25  - <out.cyf>\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();

  // no table to work with
  if (!table.size())
    return 1;
  
  table.SetupOutputWriter(opt::outfile);

  table.AnnotateCall(n, frac, dist, flag_to_set);

  table.OutputTable();

  return 0;
  
}

static int log10func(int argc, char** argv)  {

  const char* shortopts = "vn:";
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
      "Usage: cyftools log10 [.cyf file] <options>\n"
      "  Calculate the log10 of marker intensities\n"
      "  .cyf file: filepath or a '-' to stream to stdin\n"
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // set table params
  table.setVerbose(opt::verbose);

  LogProcessor logp;
  logp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);  

  if (table.StreamTable(logp, opt::infile))
    return 1; // non-zero status on StreamTable

  return 0;
}

static int dividefunc(int argc, char** argv)  {

  const char* shortopts = "vn:d:z:";
  std::string numerator, denominator;
  float div_zero = -1;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> numerator; break;
    case 'd' : arg >> denominator; break;
    case 'z' : arg >> div_zero; break;            
    default: die = true;
    }
  }

  if (numerator.empty() || denominator.empty())
    die = true;
  
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools divide [.cyf file] <options>\n"
      "  Divide two columns by each other\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -v, --verbose             Increase output to stderr\n"
      "    -n                        Name of numerator column\n"
      "    -d                        Name of denominator column\n"
      "    -z [-1]                   Value given to divide-by-zero\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }
  
  // set table params
  table.setVerbose(opt::verbose);

  DivideProcessor divp;
  divp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  divp.SetParams(numerator, denominator, div_zero);

  if (table.StreamTable(divp, opt::infile))
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
	 opt::module == "crop"  || opt::module == "umap" ||
	 opt::module == "count" || opt::module == "clean" ||
	 opt::module == "annotate" || opt::module == "convolve" ||
	 opt::module == "flagset" || opt::module == "markcheck" || 
	 opt::module == "cat" || opt::module == "convert" ||
	 opt::module == "sort" || opt::module == "divide" || 
	 opt::module == "pearson" || opt::module == "info" ||
	 opt::module == "cut" || opt::module == "view" ||
	 opt::module == "delaunay" ||  
	 opt::module == "mean" || opt::module == "ldacreate" ||
	 opt::module == "ldarun" || opt::module == "png" ||
	 opt::module == "tls" || opt::module == "dist" ||
	 opt::module == "for" || opt::module == "offset" || 
	 opt::module == "radialdens" || opt::module == "margin" ||
	 opt::module == "scramble" || opt::module == "scatter" ||
	 opt::module == "hallucinate" || opt::module == "summary" ||
	 opt::module == "magnify" || opt::module == "flip" ||
	 opt::module == "check" || 
	 opt::module == "island" || opt::module == "rescale" || 
	 opt::module == "jaccard" || opt::module == "cellcount" ||
	 opt::module == "synth" || 
	 opt::module == "sampleselect" || opt::module == "dbscan" || 
	 opt::module == "filter" || opt::module == "pheno")) {
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
  const char* shortopts = "vr:m:";
  bool blacklist = true;
  float micron_per_pixel = 0;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
      //case 'n' : arg >> opt::n; break;
    case 'r' : arg >> roifile; break;
    case 'm' : arg >> micron_per_pixel; break;
      //case 'b' : blacklist = true; break;
    default: die = true;
    }
  }

  if (micron_per_pixel == 0) {
    std::cerr << "WARNING: Must input a micron-per-pixel scale factor" << std::endl;
    die = true;
  }
    
  // make sure ROI file exists/readable
  if (!check_readable(roifile)) {
    std::cerr << "Error: " << roifile << " not readable/exists" << std::endl;
    die = true;
  }
  
  if (die || roifile.empty() || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools roi [.cyf file] <options>\n"
      "  Subset or label the cells to only those contained in the rois\n"
      "  .cyf file: filepath or a '-' to stream to stdin\n"
      "  -r <roifile>              ROI file\n"
      "  -m <float>                Microns per pixel scale factor\n"
      //      "  -l                        Output all cells and add \"roi\" column with ROI label\n"      
      "  -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // read in the roi file
  std::vector<Polygon> rois = read_polygons_from_file(roifile);

  // scale the polygons
  float maxr = 01;
  for (auto& p : rois) {
    for (auto& v : p.vertices) {
      if (v.x > maxr)
	maxr = v.x; 
      v.x *= micron_per_pixel;
      v.y *= micron_per_pixel;
    }
  }

  if (opt::verbose)
    for (const auto& c : rois)
      std::cerr << c << std::endl;

  ROIProcessor roip;
  roip.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  roip.SetParams(false, rois, blacklist);// false is placeholder for label function, that i need to implement

  if (table.StreamTable(roip, opt::infile))
    return 1; // non-zero status in StreamTable

  return 0;
  
}

static int viewfunc(int argc, char** argv) {
  
  const char* shortopts = "vn:hHRAtCx:X:";
  int precision = 2;
  bool rheader = false;  // view as csv with csv header
  bool adjacent = false; // view as name:value
  bool crevasse = false; // view as output for crevasse
  bool strict_cut = false; // only print the cut field
  bool tabprint = false; // view as tab-delimited and columns justified
  std::string cut; // list of markers, csv separated, to cut on
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> precision; break;
    case 'C' : crevasse = true; break;
    case 'X' : arg >> cut; strict_cut = true; break;
    case 'x' : arg >> cut; break;
    case 't' : tabprint = true; break;
    case 'A' : adjacent = true; break;
    case 'R' : rheader = true; break;
    case 'h' : opt::header = true; break;
    case 'H' : opt::header_only = true; break;
    default: die = true;
    }
  }

  if (rheader + adjacent + crevasse + opt::header_only + opt::header > 1) {
    std::cerr << "Warning: Can only chose one of A, R, h, H, C" << std::endl;
    die = true;
  }
  
  if (die || in_only_process(argc, argv)) {

    const char *USAGE_MESSAGE = 
      "Usage: cyftools view <.cyf file> [options]\n"
      "  View the contents of a cell table and optionally modify the output.\n"
      "\n"
      "Arguments:\n"
      "  <.cyf file>                  File path or '-' to stream from stdin.\n"
      "\n"
      "Optional Options:\n"
      "  -R                         Print (with header) in csv format\n"
      "  -A                         Print as name:value format for viewing\n"
      "  -C                         Print as CellID,x,y,...markers... for Crevasse\n"
      "  -t                         Print tab-delimited and const space\n"
      "  -x <fields>                Comma-separated list of fields to trim output to.\n"
      "  -X <fields>                Like -x, but ONLY output those columns, no CellID, etc\n"
      "  -n <decimals>              Number of decimals to keep. Default is -1 (no change).\n"
      "  -H                         View only the header.\n"
      "  -h                         Output with the header.\n"
      "  -v, --verbose              Increase output to stderr.\n"
      "\n"
      "Example:\n"
      "  cyftools view input.cyf -n 2\n"
      "  cyftools view - -H\n";
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
  table.setVerbose(opt::verbose);
  
  ViewProcessor viewp;
  viewp.SetParams(opt::header, opt::header_only, rheader, adjacent, crevasse,
		  precision, tokens, tabprint, strict_cut);

  table.StreamTable(viewp, opt::infile);
  
   return 0;  
}

static int histogramfunc(int argc, char** argv) {
  const char* shortopts = "vn:w:";
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
      "Usage: cyftools histogram [.cyf file] <options>\n"
      "  Calculate the histogram of a set of markers\n"
      "  .cyf file: filepath or a '-' to stream to stdin\n"
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
  const char* shortopts = "vl:w:";      
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
      "Usage: cyftools plot [.cyf file] <options>\n"
      "  Outputs an ASCII-style plot of cell locations\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
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

static int headfunc(int argc, char** argv) {
  const char* shortopts = "vn:";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> opt::n; break;
    default: die = true;
    }
  }

  // Process any remaining no-flag options
  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE =
      "Usage: cyftools head [.cyf file] <options>\n"
      "  Keep only the first n cells\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -n, --numrows             Number of rows to keep\n"
      "    -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // set table params
  table.setVerbose(opt::verbose);

  HeadProcessor headp;
  headp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);  
  headp.SetParams(opt::n);
    
  if (table.StreamTable(headp, opt::infile))
    return 1; // non-zero status on StreamTable

  
  //build_table();
  
  // subsample
  //table.Subsample(opt::n, opt::seed);

  //table.SetupOutputWriter(opt::outfile);
  
  // print it
  //table.OutputTable();
  
  return 0;
  
}

 static int sampleselectfunc(int argc, char** argv) {
   const char* shortopts = "vs:";
   int samplenum = -1;
   std::string tmpstring;
   std::string samplestring;
   for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
     std::istringstream arg(optarg != NULL ? optarg : "");
     switch (c) {
     case 'v' : opt::verbose = true; break;
     case 's' :
       arg >> samplenum;
       // add to existing
       if (samplestring.length()) {
	 tmpstring.clear();
	 arg >> tmpstring;
	 samplestring = tmpstring + "," + samplestring;
       } else {
	 arg >> samplestring;
       }
       break;
     default: die = true;
    }
   }

   std::unordered_set<uint32_t> samples;
   std::stringstream ss(samplestring);
   std::string token;
   while (std::getline(ss, token, ',')) {
     int sam = std::stoi(token);
     assert(sam >= 0);
     samples.insert(sam);
   }
   
   // Process any remaining no-flag options
   if (die || in_out_process(argc, argv) || samples.size() == 0) {
     
     const char *USAGE_MESSAGE =
       "Usage: cyftools sampleselect [.cyf file] <options>\n"
       "  Selects a particular sample from the cell quantification table\n"
       "    .cyf file: filepath or a '-' to stream to stdin\n"
       "    -s                        Sample number to select on\n"
       "    -v, --verbose             Increase output to stderr"
       "\n";
     std::cerr << USAGE_MESSAGE;
     return 1;
   }

  SubsampleProcessor subs;
  subs.SetCommonParams(opt::outfile, cmd_input, opt::verbose);      
  subs.SetParams(1, opt::seed, samples); // will just ignore seed. Rate is 1
  table.StreamTable(subs, opt::infile);

  return 0;
  
 }

 
 
 static int subsamplefunc(int argc, char** argv) {

   const char* shortopts = "vn:s:r:";
  float rate = 0;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'n' : arg >> opt::n; break;
    case 's' : arg >> opt::seed; break;
    case 'r' : arg >> rate; break;
    default: die = true;
    }
  }

  // Process any remaining no-flag options
  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE =
      "Usage: cyftools subsample [.cyf file] <options>\n"
      "  Subsamples a cell quantification table, randomly.\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -n, --numrows             Number of rows to subsample. Reads full file into memory!\n"
      "    -r                        Subsample at given rate (0,1]. Overrides -n. Streams without memory load, but doesn't guarentee given n outcome\n"
      "    -s, --seed           [42] Seed for random subsampling\n"
      "    -v, --verbose             Increase output to stderr"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // build the table
  if (rate <= 0) {

    std::cerr << "...reading table into memory. If not desired, use -r instead (no guarentee on N but no memory overhead)" << std::endl;
    
    build_table();
  
    // subsample
    table.Subsample(opt::n, opt::seed);
    
    table.SetupOutputWriter(opt::outfile);
    
    // print it
    table.OutputTable();
    
  } else {

    SubsampleProcessor subs;
    subs.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
    subs.SetParams(rate, opt::seed, std::unordered_set<uint32_t>()); // the -1 is so that all samples are selected
    table.StreamTable(subs, opt::infile);
    
  }
  return 0;
  
}

static int pearsonfunc(int argc, char** argv) {
  const char* shortopts = "vjs";
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
      "Usage: cyftools pearson [.cyf file] <options>\n"
      "  Outputs an ASCII-style plot of marker intensity Pearson correlations\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
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
  const char* shortopts = "vc:";
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
      "Usage: cyftools crop [.cyf file] <options>\n"
      "  Crop the table to a given rectangle (in pixels)\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
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
/*
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
      "Usage: cyftools spatial [.cyf file]\n"
      "  Construct the Euclidean KNN spatial graph\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
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

  return 0;
}
*/
 
 // static int trimfunc(int argc, char** argv) {
 //   const char* shortopts = "v";
   
 //   for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
 //     std::istringstream arg(optarg != NULL ? optarg : "");
 //     switch (c) {
 //     case 'v': opt::verbose = true; break;
 //     }   
 //   }

 //  if (die || in_out_process(argc, argv)) {
  
 //    const char *USAGE_MESSAGE =
 //      "Usage: cyftools trim [.cyf file]\n"
 //      "  Keep only marked cells\n"
 //      "    .cyf file: filepath or a '-' to stream to stdin\n"
 //      "  Options\n"
 //      "    -v, --verbose         Increase output to stderr\n"
 //      "  Example\n"
 //      "    cyftools trim <in> <out>"
 //      "\n";
 //    std::cerr << USAGE_MESSAGE;
 //    return 1;
 //  }

 //  TrimProcessor trim;
 //  trim.SetCommonParams(opt::outfile, cmd_input, opt::verbose);

 //  // process
 //  if (table.StreamTable(trim, opt::infile))
 //    return 1;

 //  return 0;
  
 // }
 
 static int filterfunc(int argc, char** argv) {
   
   const char* shortopts = "vr:a:n:A:N:of:g:l:G:L:e:jr:s:S:MB";

   bool mark1_cells = false;
   bool mark2_cells = false;   
   
   // field selection
   bool or_toggle = false;
   std::string sholder;
   float holder;
   
   // radius for mask select
   float radius = 0.0f;

   // SelectOp = pair<enum of comparators, float value>
   SelectOpVec num_vec;
   SelectOpMap criteria;
   StringVec field_vec;

   // short hand for selecting multiple flags as OR
   // (-s or -S), better than enumerating every flag with -a -o construct
   cy_uint single_flag_por = 0;
   cy_uint single_flag_cor = 0;   
   
   // a SelectionUnit is a logical set of flag conditions
   // and filter will OR all of the different selection units
   std::vector<SelectionUnit> selections;
   selections.push_back(SelectionUnit());
   
   for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
     std::istringstream arg(optarg != NULL ? optarg : "");
     switch (c) {
     case 'v': opt::verbose = true; break;
     case 'M' : mark1_cells = true; break;
     case 'B' : mark2_cells = true; break;       
     case 'a':
       if (selections.back().pand) // already have an and condition, so start new selection unit
	 selections.push_back(SelectionUnit());
       arg >> selections.back().pand; break;
     case 'n':
       if (selections.back().pnot)
	 selections.push_back(SelectionUnit());
       arg >> selections.back().pnot; break;
     case 'A':
       if (selections.back().cand)
	 selections.push_back(SelectionUnit());	
       arg >> selections.back().cand; break;
     case 'N':
       if (selections.back().cnot) 
	 selections.push_back(SelectionUnit());		
       arg >> selections.back().cnot; break;
     case 'o': selections.push_back(SelectionUnit()); break;
     case 'f' : arg >> sholder; field_vec.push_back(sholder); break;
     case 's':  arg >> selections.back().por; break;
     case 'S':  arg >> selections.back().cor; break;
     case 'g' : arg >> holder; num_vec.push_back({optype::GREATER_THAN, holder}); break;
     case 'l' : arg >> holder; num_vec.push_back({optype::LESS_THAN, holder}); break;
     case 'G' : arg >> holder; num_vec.push_back({optype::GREATER_THAN_OR_EQUAL, holder}); break;
     case 'L' : arg >> holder; num_vec.push_back({optype::LESS_THAN_OR_EQUAL, holder}); break;
     case 'e' : arg >> holder; num_vec.push_back({optype::EQUAL_TO, holder}); break;
     case 'j' : or_toggle = true; break;
     case 'r' : arg >> radius; break;
     default: die = true; break;  
     }
   }

   // mark2 not really exposted, so shouldn't get here
   if (mark1_cells && mark2_cells) {
     std::cerr << "Error: cyftools filter - select only -m or -M, can't mark both\n" << std::endl;
     die = true;
   }

   // send the selector unit to the CellSelector object
   CellSelector cellselect;
   for (const auto& a : selections) {
     if (!a.isEmpty())
       cellselect.AddSelectionUnit(a);
   }

   if (opt::verbose)
     std::cerr << cellselect;
   
   if (field_vec.size() > 1) {
     
     // make sure fields and criteria line up
     if (num_vec.size() != field_vec.size()) {
       std::cerr << "Must specify same number of fields (or just 1 field) as criteria" << std::endl;
       die = true;
     } else {
       for (size_t i = 0; i < num_vec.size(); i++) {
	 criteria[field_vec.at(i)].push_back(num_vec.at(i));
       }
     }
     
   } else if (field_vec.size() == 1) {
     
     if (num_vec.size() == 0) {
       std::cerr << "Must specify a numeric criteria" << std::endl;
       die = true;
     } else {
       for (const auto& c : num_vec) {
	 criteria[field_vec.at(0)].push_back(c);
       }
     }
     
   }

   // print the numeric selectors
   // criteria is string keyed map of vector(pair<op enum, float value>)
   if (opt::verbose) {
     for (const auto& m : criteria) {
       std::cerr <<   "--- Marker/Meta field: " << m.first << std::endl;
       for (const auto& v : m.second) {
	 std::cerr << "           - " << static_cast<int>(v.first) << " - value " << v.second << std::endl;     
       }
     }
   }
   
  if (die || in_out_process(argc, argv)) {
  
    const char *USAGE_MESSAGE =
      "Usage: cyftools filter [.cyf file]\n"
      "  Filter cells and mark (-M) for downstream processing or trim output (default)\n"
      "    DEFAULT behavior is to remove cells not in filter criteria and clear any marks\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "  Flag selection\n"
      "    -a <int>              Phenotype logical AND\n"
      "    -n <int>              Phenotype logical NOT\n"
      "    -s <int>              Phenotype logical OR\n"
      "    -A <int>              Cell logical AND\n"  
      "    -N <int>              Cell logical NOT\n"
      "    -S <int>              Cell logical OR\n"
      "    -o                    Flag to mark OR separator for units above\n"
      "  Marker/meta selection\n"
      "    -f <string>           Marker / meta field to select on\n"
      "    -g <int>              > - Greater than\n"
      "    -G <int>              >= - Greater than or equal to \n"      
      "    -l <int>              < - Less than\n"
      "    -L <int>              <= - Less than or equal to\n"      
      "    -e <int>              Equal to (can use with -g or -m for >= or <=)\n"
      "    -j                    Flag to OR the operations (default is and)\n"
      "  Mask selection\n"
      "    -r                    Radius (in x,y coords) of region around selected cells to also include\n"
      "  Options\n"
      "    -M                    Set the mark (keeps all cells, marks some for downstream analysis)\n"
      //"    -B                    Set the build mark ON/OFF (used to limit downstream tree building to only build-marked cells)\n"
      "    -v, --verbose         Increase output to stderr\n"
      "  Example\n"
      "    cyftools filter -a 4160 -o -a 4352 <in> <out>  # select cells with bits 4096+256 OR 4096+64\n"
      "    cyftools filter -s 4160 <in> <out>             # select cells with bits markers 4096 OR 256 on\n"
      "    cyftools filter -a 4096 -o 256                 # equivalent to above\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // setup the selector processor for zero radius
  if (radius <= 0) {

    FilterProcessor select;
    select.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
    select.SetFlagParams(cellselect, mark1_cells, mark2_cells); 
    select.SetFieldParams(criteria, or_toggle);
    
    // process
    table.StreamTable(select, opt::infile);

    return 0;
  }
  
  // or else we read table and then select
  build_table();

  table.SetupOutputWriter(opt::outfile);

  table.Select(cellselect, 
	       criteria,
	       or_toggle,
	       radius);

  table.OutputTable();

  return 0;
}
  
static int phenofunc(int argc, char** argv) {
  
  const char* shortopts = "vt:s:r:";
  std::string file;
  float scale = 1;
  float random_scale = 0;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> file; break;
    case 's' : arg >> scale; break;
    case 'r' : arg >> random_scale; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {

    const char *USAGE_MESSAGE =
      "Usage: cyftools pheno [.cyf file]\n"
      "  Phenotype cells (set the flags) with threshold file\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -t <file>        File that holds gates: marker(string), low(float), high(float)\n"
      "    -s <float>       How much to scale the gate, for purposes of robustness testing (s=1 is no scaling)\n"
      "    -r <float>       Randomly vary the gates by +/- factor of r (e.g. r=0.2 -> s=[0.8-1.2])\n"
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
  phenop.SetParams(pheno, scale, random_scale);

  if (table.StreamTable(phenop, opt::infile))
    return 1; // non-zero StreamTable status
  
  return 0;
}

static void assign_uint64_t_from_int(std::istream& stream, uint64_t& variable) {
    int temp;
    if (stream >> temp) {
        if (temp < 0) {
            variable = std::numeric_limits<uint64_t>::max();
        } else {
            variable = static_cast<uint64_t>(temp);
        }
    } else {
        std::cerr << "Invalid input for uint64_t assignment" << std::endl;
    }
}
 
static int radialdensfunc(int argc, char** argv) {
  const char* shortopts = "vR:r:o:a:l:f:jJt:";
  cy_uint inner = 0;
  cy_uint outer = 20;
  cy_uint logor = 0;
  cy_uint logand = 0;
  std::string label;

  bool normalize_local = false; // normalize to cells count in radius
  bool normalize_global = false; // normalize to cells count in slide
  std::string file;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'f' : arg >> file; break;      
    case 't' : arg >> opt::threads; break;
    case 'r' : arg >> inner; break;
    case 'R' : arg >> outer; break;
    case 'o' : assign_uint64_t_from_int(arg, logor); break;
    case 'a' : assign_uint64_t_from_int(arg, logand); break;      
    case 'j' : normalize_local = true; break;
    case 'J' : normalize_global = true; break;
    case 'l' : arg >> label; break;
    default: die = true;
    }
  }

  if (label.empty() && file.empty()) {
    std::cerr << "Error: Must specify file with -f or individual params with -l etc" << std::endl;
    die = true;
  }
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools radialdens [.cyf file]\n"
      "  Calculate the density of cells away from individual cells\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -f                    File for multiple labels [r,R,o,a,j,J,l]\n"
      "    -r [0]                Inner radius\n"
      "    -R [20]               Outer radius\n"
      "    -o                    Logical OR flags\n"
      "    -a                    Logical AND flags\n"
      "    -j                    Flag for normalizing to cell count in radius\n"
      "    -J                    Flag for normalizing to cell count in slide\n"
      "    -l                    Label the column\n"
      "    -t                    Number of threads (default 1)\n"
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
    
  }
 
  
  // streaming way
  //RadialProcessor radp;
  //radp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);

  std::vector<cy_uint> innerV(rsv.size());
  std::vector<cy_uint> outerV(rsv.size());  
  std::vector<cy_uint> logorV(rsv.size());
  std::vector<cy_uint> logandV(rsv.size());
  std::vector<int>     normalize_localV(rsv.size());
  std::vector<int>     normalize_globalV(rsv.size());  
  StringVec labelV(rsv.size());
  if (rsv.empty()) {
    innerV = {inner};
    outerV = {outer};
    logorV = {logor};
    logandV= {logand};
    labelV = {label};
    normalize_localV={normalize_local};
    normalize_globalV={normalize_global};    
  } else {
    for (size_t i = 0; i < innerV.size(); i++) {
      innerV[i] = rsv.at(i).int_data.at(0);
      outerV[i] = rsv.at(i).int_data.at(1);
      logorV[i] = rsv.at(i).int_data.at(2);
      logandV[i] = rsv.at(i).int_data.at(3);
      normalize_localV[i] = rsv.at(i).int_data.at(4);
      normalize_globalV[i] = rsv.at(i).int_data.at(5);
      labelV[i] = rsv.at(i).label;
    }
  }

  // building way
  build_table();

  table.SetupOutputWriter(opt::outfile);
  
  table.RadialDensityKD(innerV, outerV, logorV, logandV, labelV,
			normalize_localV, normalize_globalV);

  table.OutputTable();
  
  return 0;
  
  /*if (table.StreamTable(radp, opt::infile))
    return 1; // non-zero exit from StreamTable

    return 0;*/
}

static void cyftools_cat(const StringVec& inputFiles, const std::string& outputFile) {
    std::ofstream os(outputFile, std::ios::binary);
    cereal::PortableBinaryOutputArchive oarchive(os);

    for (size_t i = 0; i < inputFiles.size(); ++i) {
        std::ifstream is(inputFiles[i], std::ios::binary);
        cereal::PortableBinaryInputArchive iarchive(is);

	std::cerr << " WORKING ON " << i << " " << inputFiles[i] << std::endl;
	
        CellHeader header;
        if (i == 0) {  // For the first file, serialize the header to the output
	  std::cerr << " deserializing header " << std::endl;
            iarchive(header);
	  std::cerr << " erializing header " << std::endl;	    
            oarchive(header);
        } else {  // For subsequent files, skip the header
            iarchive(header);
        }

        // Read, then immediately write each cell to the output
        Cell cell;
        while (is) {
	  //std::cerr << "reading " << std::endl;
            iarchive(cell);
            oarchive(cell);
        }
    }
}

 int moranfunc(int argc, char** argv) {

  const char* shortopts = "vR:r:o:a:l:f:t:";
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
      "Usage: cyftools moran [.cyf file]\n"
      "  Calculate the Moran's I cells in neighborhoods\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -r [20]               Outer radius\n"
      "    -R [0]                Inner radius\n"
      "    -o                    Logical OR flags\n"
      "    -a                    Logical AND flags\n"
      "    -l                    Label the column\n"
      "    -f                    File for multiple labels [r,R,o,a,l]\n"
      "    -t                    Number of threads (default 1)\n"
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();

  table.SetupOutputWriter(opt::outfile);

  //table.MoranI();

  table.OutputTable();

  return 0;
 }
 
int debugfunc(int argc, char** argv) {
  const char* shortopts = "v";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    std::cerr << "invalid input" << std::endl;
    return 1;
  }

  DebugProcessor debug;
  debug.SetCommonParams(opt::outfile, cmd_input, opt::verbose); // really shouldn't need any of these

  if (table.StreamTable(debug, opt::infile)) 
    return 1; // non-zero status on StreamTable

  return 0;
}

static int convertfunc(int argc, char** argv) {
  const char* shortopts = "s:vm:u:c:";
  uint32_t sampleid = 0; //static_cast<uint32_t>(-1);
  std::string metacols;
  std::string units;
  std::string microns_to_pixels;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 's' : arg >> sampleid; break;
    case 'v' : opt::verbose = true; break;
    case 'm' : arg >> metacols; break;
    case 'u' : arg >> units; break;
    case 'c' : arg >> microns_to_pixels; break;      
    default: die = true;
    }
  }

  /*if (sampleid == static_cast<uint32_t>(-1)) {
    die = true;
    }*/
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools convert <csvfile> [.cyf file]\n"
      "  Convert a CSV file to a .cyf formatted file or stream.\n"
      "\n"
      "Arguments:\n"
      "  <csvfile>                 Input CSV file path.\n"
      "  [.cyf file]                 Output .cyf file path or '-' to stream as a cyf-formatted stream to stdout.\n"
      "\n"
      "Required Options:\n"
      "  -s, --sampleid <int>      Provide a unique sample id (>= 0).\n"
      "\n"
      "Optional Options:\n"
      "  -v, --verbose             Increase output to stderr.\n"
      "  -h,                       Optionally, supply the header text file here, rather than as appended to csv\n"
      "  -m,                       List metadata headers as a comma-delimited string (i.e. 'Area,Size,LDA')"
      "\n"
      "Example:\n"
      "  cyftools convert input.csv output.cyf\n"
      "  cyftools convert input.csv - -s 12345\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  CerealProcessor cerp;
  cerp.SetParams(opt::outfile, cmd_input, sampleid);

  table.StreamTableCSV(cerp, opt::infile, metacols);

  return 0;
}

static int delaunayfunc(int argc, char** argv) {
  const char* shortopts = "vt:D:V:l:";
  std::string delaunay;
  std::string voronoi;
  int limit = -1;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 't' : arg >> opt::threads; break;
    case 'D' : arg >> delaunay; break;
    case 'V' : arg >> voronoi; break;
    case 'l' : arg >> limit; break;
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools delaunay [.cyf file]\n"
      "  Perform a the Delaunay triangulation of a cell table\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -t [1]                    Number of threads\n"
      "    -D                        Filename of PDF to output of Delaunay triangulation\n"
      "    -V                        Filename of PDF to output of Voronoi diagram\n"
      "    -l                        Size limit of an edge in the Delaunay triangulation\n"
      "    -v, --verbose             Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // build the table
  // but don't have to convert columns
  // since we don't use pre-existing Graph or Flags for this
  build_table();

  // check we were able to read the table
  if (table.CellCount() == 0) {
    std::cerr << "Ending with no cells? Error in upstream operation?" << std::endl;
    return 0;
  }

  table.SetupOutputWriter(opt::outfile);

  table.Delaunay(delaunay, voronoi, limit);
  
  table.OutputTable();
  
  return 0;
  
}

static int sortfunc(int argc, char** argv) {
  const char* shortopts = "yx:jv";
  bool xy = false;
  std::string field;
  bool reverse = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'y' : xy = true; break;
    case 'x' : arg >> field; break;
    case 'j' : reverse = true; break;
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if ( !xy && field.empty() ) {
    die = true;
    std::cerr << "Must select only one of -y flag or -x <arg>\n" << std::endl;
  }
  
  if (die || in_out_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools sort [.cyf file]\n"
      "  Sort cells\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -y                    Flag to have cells sort by (x,y), in increasing distance from 0\n"
      "    -x                    Field to sort on\n"
      "    -j                    Reverse sort order\n"
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  build_table();

  // check we were able to read the table
  if (table.CellCount() == 0) {
    std::cerr << "Ending with no cells? Error in upstream operation?" << std::endl;
    return 0;
  }
  
  // setup output writer
  table.SetupOutputWriter(opt::outfile);

  if (xy)
    table.sortxy(reverse);

  if (!field.empty())
    table.sort(field, reverse);
  
  // print it
  table.OutputTable();

  return 0;
}


static int countfunc(int argc, char** argv) {
  
  const char* shortopts = "v";
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    default: die = true;
    }
  }

  if (die || in_only_process(argc, argv)) {
    
    const char *USAGE_MESSAGE =
      "Usage: cyftools count [.cyf file]\n"
      "  Output the number of cells in a file\n"
      "    .cyf file: filepath or a '-' to stream to stdin\n"
      "    -v, --verbose         Increase output to stderr\n"
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  CountProcessor countp;
  countp.SetCommonParams(opt::outfile, cmd_input, opt::verbose); // really shouldn't need any of these

  if (table.StreamTable(countp, opt::infile)) 
    return 1; // non-zero status on StreamTable
  
  countp.PrintCount();

  return 0;
}


// return TRUE if you want the process to die and print message
static bool in_out_process(int argc, char** argv) {
  
  optind++;
  // Process any remaining no-flag options
  size_t count = 0;
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    } else if (opt::outfile.empty()) {
      opt::outfile = argv[optind];
    }
    count++;
    optind++;
  }

  // there should be only 2 non-flag input
  if (count > 2) {
    std::cerr << "Error: >2 in/out detected. Need 2 (input and output)" << std::endl;
    return true;
  }
  // die if no inputs provided
  if (count == 0) {
    std::cerr << "Error:NO in/out detected. Need 2 (input and output)" << std::endl;    
    return true;
  }

  // check if the input file is readable
  if (!check_readable(opt::infile) && opt::infile != "-") {
    std::cerr << "Error: Input file " << opt::infile << " not readable/exists" << std::endl;
    return true;
  }

  // check if the input file is readable
  if (!check_writeable_folder(opt::outfile) && opt::outfile != "-") {
    std::cerr << "Error: Output file " << opt::outfile << " not readable/exists" << std::endl;
    return true;
  }

  return opt::infile.empty() || opt::outfile.empty();


}
// return TRUE if you want the process to die and print message
static bool in_only_process(int argc, char** argv) {
  
  optind++;
  // Process any remaining no-flag options
  size_t count = 0;
  while (optind < argc) {
    if (opt::infile.empty()) {
      opt::infile = argv[optind];
    }
    count++;
    optind++;
  }

  // there should be only 1 non-flag input
  if (count > 1)
    return true;
  
  // die if no inputs provided
  if (count == 0)
    return true;
  
  if (!check_readable(opt::infile) && opt::infile != "-") {
    std::cerr << "Error: File " << opt::infile << " not readable/exists" << std::endl;
    return true;
  }
  
  return opt::infile.empty();

}

static bool out_only_process(int argc, char** argv) {
  
  optind++;
  // Process any remaining no-flag options
  size_t count = 0;
  while (optind < argc) {
    if (opt::outfile.empty()) {
      opt::outfile = argv[optind];
    }
    count++;
    optind++;
  }
  
  // there should be only 1 non-flag input
  if (count > 1)
    return true;
  
  // die if no inputs provided
  if (count == 0)
    return true;
  
  return opt::outfile.empty();

}

 Polygon::Polygon(const std::vector<JPoint>& vec) {
   Id = -1;
   Name = "na";
   type = "na";
   vertices = vec;
 }
