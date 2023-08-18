#include "data_structures.h"


namespace dim1 {
    class BettiMatching {
        public:
        BettiMatching(vector<value_t> input0, vector<value_t> input1, vector<value_t> comparison, vector<index_t> shape,
                        Config& config);
        void computeMatching();
        void computeVoxels();
        void printResult();

        private:
        const CubicalGridComplex cgc0;
        const CubicalGridComplex cgc1;
        const CubicalGridComplex cgcComp;
        vector<Pair> pairs0;
        vector<Pair> pairs1;
        vector<Pair> pairsComp;
        vector<Match> matches;
        unordered_map<uint64_t, bool> isMatched0;
        unordered_map<uint64_t, bool> isMatched1;
        vector<VoxelMatch> matched;
        vector<VoxelPair> unmatched0;
        vector<VoxelPair> unmatched1;
        Config& config;
    };
}