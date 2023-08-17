#include "data_structures.h"


namespace D2 {
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
        vector<vector<Pair>> pairs0;
        vector<vector<Pair>> pairs1;
        vector<vector<Pair>> pairsComp;
        vector<vector<Match>> matches;
        vector<unordered_map<uint64_t, bool>> isMatched0;
        vector<unordered_map<uint64_t, bool>> isMatched1;
        vector<vector<VoxelMatch>> matched;
        vector<vector<VoxelPair>> unmatched0;
        vector<vector<VoxelPair>> unmatched1;
        vector<Cube> ctr0;
        vector<Cube> ctr1;
        vector<Cube> ctrComp;
        Config& config;
    };
}