#include <functional>

// Define the function wrapper type
//typedef std::function<bool(const std::string& in, std::string& out)> LineStreamerWrapper;

// define phenotype map
typedef std::pair<float,float> Pheno;
typedef std::unordered_map<std::string, Pheno> PhenoMap;
