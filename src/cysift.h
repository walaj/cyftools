#include <unordered_map>

// define phenotype map
typedef std::pair<float,float> Pheno;
typedef std::unordered_map<std::string, Pheno> PhenoMap;

#ifdef USE_64_BIT
typedef uint64_t cy_uint;
#else
typedef uint32_t cy_uint;
#endif
