#include "data_structures.h"
#include "../data_structures.h"
#include <unordered_map>


namespace dim1 {
    class BettiMatching {
        public:
        BettiMatching(vector<value_t> input0, vector<value_t> input1, vector<value_t> comparison, vector<index_t> shape,
                        Config& config);
        BettiMatching(BettiMatching &&other);
        void computeMatching();
        void computeVoxels();
        void printResult();

        // Read-only properties
        const vector<VoxelMatch> &matched = _matched;
        const vector<VoxelPair> &unmatched0 = _unmatched0;
        const vector<VoxelPair> &unmatched1 = _unmatched1;

        private:
        CubicalGridComplex cgc0;
        CubicalGridComplex cgc1;
        CubicalGridComplex cgcComp;
        vector<Pair> pairs0;
        vector<Pair> pairs1;
        vector<Pair> pairsComp;
        vector<Match> matches;
        unordered_map<uint64_t, bool> isMatched0;
        unordered_map<uint64_t, bool> isMatched1;
        vector<VoxelMatch> _matched;
        vector<VoxelPair> _unmatched0;
        vector<VoxelPair> _unmatched1;
        Config& config;
    };
}