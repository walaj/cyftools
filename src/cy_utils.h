#ifndef SPADF_UTILS_H
#define SPADF_UTILS_H

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <numeric>
#include <iostream>
#include <random>
#include <set>


/** Format an integer to include commas
 * @param data Number to format
 * @return String with formatted number containing commas
 */
template <typename T> inline std::string AddCommas(T data) {
  std::stringstream ss; 
  ss << std::fixed << data; 
  std::string s = ss.str();
  if (s.length() > 3)
    for (int i = s.length()-3; i > 0; i -= 3)
      s.insert(i,",");
  return s;
}

/*std::set<int> sampleWithoutReplacement(int n, int N, int seed) {

  assert(n <= N);
  std::set<int> s;
  std::vector<int> v(N);
  for (int i = 0; i < N; i++) {
    v[i] = i; // initialize the vector with 0-based indices
  }
  
  // use a random device to generate random numbers
  std::random_device rd;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> dis(0, N-1);
  
  // randomly sample n indices without replacement
    while (s.size() < n) {
      int index = dis(gen); // generate a random index
      if (s.count(index) == 0) { // check if the index has already been sampled
	s.insert(index); // if not, insert it into the set
      }
    }
    
    // convert the sampled indices back to 0-based indexing and return the set
    std::set<int> sampledIndices;
    for (int index : s) {
      sampledIndices.insert(v[index]);
    }
    return sampledIndices;
}
*/

#endif
