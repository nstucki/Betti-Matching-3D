#include "data_structures.h"


namespace dimN {
    class BettiMatching {
        public:
        BettiMatching(vector<value_t> input0, vector<value_t> input1, vector<value_t> comparison, vector<index_t> shape,
                        Config& config);
        void computeMatching();
        void computeVoxels();
        void printResult();

        private:
        size_t dim;
        const CubicalGridComplex cgc0;
        const CubicalGridComplex cgc1;
        const CubicalGridComplex cgcComp;
        vector<vector<Pair>> pairs0;
        vector<vector<Pair>> pairs1;
        vector<vector<Pair>> pairsComp;
        vector<vector<Match>> matches;
        unordered_map<index_t, bool> isMatched0;
        unordered_map<index_t, bool> isMatched1;
        vector<vector<VoxelMatch>> matched;
        vector<vector<VoxelPair>> unmatched0;
        vector<vector<VoxelPair>> unmatched1;
        Config& config;
    };
}