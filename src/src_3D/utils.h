#include "data_structures.h"
#include "config.h"

#include <vector>
#include <unordered_map>

using namespace std;

void readImage(string const &filename, file_format const &format, vector<double> &image, vector<uint32_t> &shape);
void tokenize(string const &str, const char delim, vector<string> &out);
void print_result(CubicalGridComplex* const cgc0, const CubicalGridComplex* const cgc1, const CubicalGridComplex* const cgcComp, 
                    const vector<vector<Pair>>& pairs0, const vector<vector<Pair>>& pairs1, const vector<vector<Match>>& matched, 
                    unordered_map<index_t, bool>& isMatched0, unordered_map<index_t, bool>& isMatched1);