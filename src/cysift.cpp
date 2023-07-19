#include "cell_column.h"
#include "cell_table.h"
#include "cell_processor.h"

#include <unistd.h> // or #include <getopt.h> on Windows systems
#include <getopt.h>
#include <ctime>
#include <regex>

#include "cell_row.h"

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

// process in and outfile cmd arguments
static bool in_out_process(int argc, char** argv);

// process in cmd arguments
static bool in_only_process(int argc, char** argv);

static CellTable table;

static void build_table();

static const char* shortopts = "jhHNyvmMPr:e:g:G:t:a:i:A:O:d:b:c:s:k:n:r:w:l:L:x:X:o:R:f:D:V:z:";
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
  //{ "sort",                       no_argument, NULL, 'y' },  
  { "csv",                        no_argument, NULL, 'j'},
  { NULL, 0, NULL, 0 }
};

static const char *RUN_USAGE_MESSAGE =
"Usage: cysift [module] <options> \n"
"  view       - View the cell table\n"
"  count      - Output number of cells in table\n"  
"  cut        - Select only given markers and metas\n"
"  divide     - Divide two columns\n"
"  head       - Keep the first lines of a file\n"
"  clean      - Removes data to decrease disk size\n"
"  delaunay   - Calculate the Delaunay triangulation\n"
"  average    - Average all of the data columns\n"  
"  cat        - Concatenate multiple samples\n"
"  sort       - Sort the cells\n"
"  subsample  - Subsample cells randomly\n"
  //"  plot       - Generate an ASCII style plot\n"
"  roi        - Trim cells to a region of interest within a given polygon\n"
  //"  histogram  - Create a histogram of the data\n"
"  log10      - Apply a base-10 logarithm transformation to the data\n"
  //"  correlate  - Calculate the correlation between variables\n"
"  info       - Display information about the dataset\n"
"  umap       - Construct the marker space UMAP\n"
"  spatial    - Construct the spatial KNN graph\n"
"  tumor      - Set the tumor flag\n"
"  select     - Select by cell phenotype flags\n"
"  pheno      - Phenotype cells to set the flag\n"
"  convolve   - Density convolution to produce TIFF\n"
"  radialdens - Calculate density of cells within a radius\n"
"  cereal     - Create a .cys format file from a CSV\n"
"  ldacreate  - Create a Latent Dirichlet model\n"
"  ldarun     - Run a Latent Dirichlet model\n"    
"\n";

static int dividefunc(int argc, char** argv);
static int sortfunc(int argc, char** argv);
static int headfunc(int argc, char** argv);
static int convolvefunc(int argc, char** argv);
static int delaunayfunc(int argc, char** argv);
static int tumorfunc(int argc, char** argv);
static int ldacreatefunc(int argc, char** argv);
static int ldarunfunc(int argc, char** argv);
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
static int umapfunc(int argc, char** argv);
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
  } else if (opt::module == "ldacreate") {
    val = ldacreatefunc(argc, argv);
  } else if (opt::module == "ldarun") {
    val = ldarunfunc(argc, argv);    
  } else if (opt::module == "clean") {
    val = cleanfunc(argc, argv);
  } else if (opt::module == "radialdens") {
    val = radialdensfunc(argc, argv);
  } else if (opt::module == "subsample") {
    val = subsamplefunc(argc, argv);
  } else if (opt::module == "plot") {
    return(plotfunc(argc, argv));
  } else if (opt::module == "divide") {
    return (dividefunc(argc, argv));
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
  } else if (opt::module == "tumor") {
    return(tumorfunc(argc, argv));
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
  } else if (opt::module == "umap") {
    val = umapfunc(argc, argv);
  } else if (opt::module == "spatial") {
    val = spatialfunc(argc, argv);    
  } else if (opt::module == "select") {
    val = selectfunc(argc, argv);
  } else if (opt::module == "convolve") {
    val = convolvefunc(argc, argv);
  } else if (opt::module == "head") {
    val = headfunc(argc, argv);
  } else if (opt::module == "sort") {
    val = sortfunc(argc, argv);
  } else if (opt::module == "delaunay") {
    val = delaunayfunc(argc, argv);
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
  table.setVerbose(opt::verbose);
  table.setThreads(opt::threads);

  // stream into memory
  BuildProcessor buildp;
  buildp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  table.StreamTable(buildp, opt::infile);
  
  table.setCmd(cmd_input);
}

static int convolvefunc(int argc, char** argv) {
 
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
      "Usage: cysift convolve [csvfile]\n"
      "  Perform a convolution to produce a TIFF\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -i                        Input TIFF file to set params for output\n"
      "    -d                        Number of microns per pixel (e.g. 0.325). Required\n"
      "    -w [200]                  Width of the convolution box (in pixels)\n"
      "    -t [1]                    Number of threads\n"      
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
      "Usage: cysift umap [csvfile]\n"
      "  Construct the UMAP (in marker space)\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -k [15]                   Number of neighbors\n"
      "    -t [1]                    Number of threads\n"
      "    -D <file>                 Optional output to PDF\n"
      "    -w [1000]                 Width of the PDF\n"
      "    -l [1000]                 Length of the PDF\n"
      "    -v, --verbose             Increase output to stderr\n"
      "\n";
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
  
#ifdef HAVE_LDAPLUSPLUS
  std::string model_out, fields;
  int n_topics = 10;
  int n_iterations = 10;
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'o' : arg >> model_out; break;

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
      "Usage: cysift ldacreate [csvfile]\n"
      "  Create a topic-model learning using Latent Dirichlet Allocation\n"
      "    <file>: filepath or a '-' to stream to stdin\n"
      "    -r                        Comma-separated list of inputs to LDA model\n"
      "    -o                        Output model file\n"
      //      "    -m                        Input model (won't re-run LDA)\n"            
      "    -n                        Number of topics [10]\n"
      "    -i                        Max number of iterations [10]\n"
      "    -t [1]                    Number of threads\n"      
      "    -v, --verbose             Increase output to stderr\n";
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
  std::vector<std::string> markers;
  for (const auto& s : tokens) {

    if (!table.ContainsColumn(s)) {
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
			 n_iterations);

  // write the model
  table.LDA_write_model(model_out);
  
  return 0;

#else
  std::cerr << "Error: Unable to run LDA without including ldaplusplus header library to build." << std::endl;
  return 1;
#endif

}

static int ldarunfunc(int argc, char** argv) {

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
      "Usage: cysift ldarun [csvfile]\n"
      "  Score cells using a pre-computed topic-model learning using Latent Dirichlet Allocation\n"
      "    <file>: filepath or a '-' to stream to stdin\n"
      "    -i                        Input model (won't re-run LDA)\n"
      "    -D                        Output a PDF plot of the topics\n"
      "    -c                        Cont cutoff\n"
      "    -l                        Highlight topic\n"      
      "    -v, --verbose             Increase output to stderr\n";
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
    case 'P' : clean_graph = true; clean_markers = true; clean_meta = true; break;      
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
      "    -P,          Remove all graph data\n"      
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
      "\n";
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
      std::cerr << "...reading and concatenating " << opt::infile_vec.at(sample_num) << std::endl;
    
    // update the offset and sample num
    catp.SetOffset(offset);
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
  table.setVerbose(opt::verbose);

  AverageProcessor avgp;
  avgp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);  

  if (table.StreamTable(avgp, opt::infile))
    return 1; // non-zero status on StreamTable

  // write the one line with the averages
  avgp.EmitCell();
  
  return 0;

}

static int tumorfunc(int argc, char** argv) {

  int n = 20;
  float frac = 0.75;
  cy_uint orflag = 0;
  cy_uint andflag = 0;  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 't' : arg >> opt::threads; break;
    case 'k' : arg >> n; break;
    case 'f' : arg >> frac; break;
    case 'o' : arg >> orflag; break;
    case 'a' : arg >> andflag; break;      
    default: die = true;
    }
  }

  if (die || in_out_process(argc, argv)) {
  
    const char *USAGE_MESSAGE =
      "Usage: cysift tumor [csvfile]\n"
      "  Set the flag on whether a cell is in the tumor region\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "    -k [20]               Number of neighbors\n"
      "    -f [0.75]             Fraction of neighbors\n"
      "    -o                    Flag OR for tumor\n"
      "    -a                    Flag AND for tumor\n"      
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

  table.TumorCall(n, frac, orflag, andflag, 600);

  table.OutputTable();

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
  table.setVerbose(opt::verbose);

  LogProcessor logp;
  logp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);  

  if (table.StreamTable(logp, opt::infile))
    return 1; // non-zero status on StreamTable

  return 0;
}

static int dividefunc(int argc, char** argv)  {

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
      "Usage: cysift divide [cysfile] <options>\n"
      "  Divide two columns by each other\n"
      "    cysfile: filepath or a '-' to stream to stdin\n"
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
	 opt::module == "tumor" || opt::module == "convolve" || 
	 opt::module == "cat" || opt::module == "cereal" ||
	 opt::module == "sort" || opt::module == "divide" || 
	 opt::module == "correlate" || opt::module == "info" ||
	 opt::module == "cut" || opt::module == "view" ||
	 opt::module == "delaunay" || opt::module == "head" || 
	 opt::module == "average" || opt::module == "ldacreate" ||
	 opt::module == "ldarun" || 
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
  table.setVerbose(opt::verbose);
  
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

static int headfunc(int argc, char** argv) {

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
      "Usage: cysift head [csvfile] <options>\n"
      "  Keep only the first n cells\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
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


static int subsamplefunc(int argc, char** argv) {

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
      "Usage: cysift subsample [csvfile] <options>\n"
      "  Subsamples a cell quantification table, randomly.\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
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
    subs.SetParams(rate, opt::seed);
    table.StreamTable(subs, opt::infile);
    
  }
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

  return 0;
}

static int selectfunc(int argc, char** argv) {


  // flag selection
  cy_uint plogor = 0;
  cy_uint plogand = 0;
  bool plognot = false;

  cy_uint clogor = 0;
  cy_uint clogand = 0;
  bool clognot = false;


  // field selection
  bool or_toggle = false;
  std::string sholder;
  float holder;
  
  SelectOpVec num_vec;
  SelectOpMap criteria;
  std::vector<std::string> field_vec;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v' : opt::verbose = true; break;
    case 'o' : arg >> plogor; break;
    case 'a' : arg >> plogand; break;
    case 'N' : plognot = true; break;
    case 'O' : arg >> clogor; break;
    case 'A' : arg >> clogand; break;
    case 'M' : clognot = true; break;
    case 'f' : arg >> sholder; field_vec.push_back(sholder); break;
    case 'g' : arg >> holder; num_vec.push_back({optype::GREATER_THAN, holder}); break;
    case 'l' : arg >> holder; num_vec.push_back({optype::LESS_THAN, holder}); break;
    case 'G' : arg >> holder; num_vec.push_back({optype::GREATER_THAN_OR_EQUAL, holder}); break;
    case 'L' : arg >> holder; num_vec.push_back({optype::LESS_THAN_OR_EQUAL, holder}); break;
    case 'e' : arg >> holder; num_vec.push_back({optype::EQUAL_TO, holder}); break;
    case 'j' : or_toggle = true; break;
    default: die = true;
    }
  }


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
  
  if (die || in_out_process(argc, argv)) {
  
    const char *USAGE_MESSAGE =
      "Usage: cysift select [csvfile]\n"
      "  Select cells by phenotype flag\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
      "  Flag selection\n"
      "    -o                    Cell phenotype: Logical OR flags\n"
      "    -a                    Cell phenotype: Logical AND flags\n"
      "    -N                    Cell phenotype: Not flag\n"
      "    -O                    Cell flag: Logical OR flags\n"
      "    -A                    Cell flag: Logical AND flags\n"
      "    -M                    Cell flag: Not flag\n"
      "  Marker selection\n"
      "    -f                    Marker / meta field to select on\n"
      "    -g                    > - Greater than\n"
      "    -G                    >= - Greater than or equal to \n"      
      "    -l                    < - Less than\n"
      "    -L                    <= - Less than or equal to\n"      
      "    -e                    Equal to (can use with -g or -m for >= or <=)\n"
      "    -j                    Or the operations (default is and)\n"
      "  Options\n"
      "    -v, --verbose         Increase output to stderr\n"      
      "\n";
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  // setup the selector processor
  SelectProcessor select;
  select.SetCommonParams(opt::outfile, cmd_input, opt::verbose);
  select.SetFlagParams(plogor, plogand, plognot, clogor, clogand, clognot); 
  select.SetFieldParams(criteria, or_toggle);
			
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

  bool normalize = false; // normalize to cells count
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
    case 'j' : normalize = true; break;
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
      "    -j                    Flag for normalizing to cell count\n"
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

  
  // streaming way
  //RadialProcessor radp;
  //radp.SetCommonParams(opt::outfile, cmd_input, opt::verbose);

  std::vector<cy_uint> innerV(rsv.size());
  std::vector<cy_uint> outerV(rsv.size());  
  std::vector<cy_uint> logorV(rsv.size());
  std::vector<cy_uint> logandV(rsv.size());
  std::vector<int>     normalizeV(rsv.size());    
  std::vector<std::string> labelV(rsv.size());
  if (rsv.empty()) {
    innerV = {inner};
    outerV = {outer};
    logorV = {logor};
    logandV= {logand};
    labelV = {label};
    normalizeV={normalize};
  } else {
    for (size_t i = 0; i < innerV.size(); i++) {
      innerV[i] = rsv.at(i).int_data.at(0);
      outerV[i] = rsv.at(i).int_data.at(1);
      logorV[i] = rsv.at(i).int_data.at(2);
      logandV[i] = rsv.at(i).int_data.at(3);
      normalizeV[i] = rsv.at(i).int_data.at(4);
      labelV[i] = rsv.at(i).label;
    }
  }

  // building way
  build_table();

  table.SetupOutputWriter(opt::outfile);
  
  table.RadialDensityKD(innerV, outerV, logorV, logandV, labelV, normalizeV);

  table.OutputTable();
  
  return 0;
  
  /*if (table.StreamTable(radp, opt::infile))
    return 1; // non-zero exit from StreamTable

    return 0;*/
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

static int delaunayfunc(int argc, char** argv) {

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
      "Usage: cysift delaunay [csvfile]\n"
      "  Perform a the Delaunay triangulation of a cell table\n"
      "    csvfile: filepath or a '-' to stream to stdin\n"
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
      "Usage: cysift sort [cysfile]\n"
      "  Sort cells\n"
      "    cysfile: filepath or a '-' to stream to stdin\n"
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
  if (count > 2)
    return true;
  
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
  
  return opt::infile.empty();

}

