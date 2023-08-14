#include "data_structures.h"
#include "config.h"
#include "betti_matching_nD.h"
#include "top_dimension.h"
#include "inter_dimensions.h"
#include "dimension_0.h"

#include <unordered_map>
#include <cassert>

BettiMatchingND::BettiMatchingND(
    const CubicalGridComplex &cgc0,
    const CubicalGridComplex &cgc1,
    const CubicalGridComplex &cgcComp,
    const Config &config
    // TODO: Depending on config is not optimal, since this includes fields not related to this class
    )
    : cgc0(cgc0), cgc1(cgc1), cgcComp(cgcComp),
      _pairs0(cgc0.dim), _pairs1(cgc0.dim), _pairsComp(cgc0.dim), _matches(cgc0.dim),
      _unmatchedPairs0(cgc0.dim), _unmatchedPairs1(cgc0.dim),
      _isMatched0(), _isMatched1(),
      dim(cgc0.dim), config(config)
{
    assert(cgc0.shape == cgc1.shape && cgc0.shape == cgcComp.shape);
};

void BettiMatchingND::computeMatching()
{
    // Reset vectors to be safe if this method is run multiple times
    _pairs0 = vector<vector<Pair>>(dim);
    _pairs1 = vector<vector<Pair>>(dim);
    _pairsComp = vector<vector<Pair>>(dim);
    _matches = vector<vector<Match>>(dim);

    if (dim > 1)
    {
        TopDimension topDim(cgc0, cgc1, cgcComp, config, _pairs0[dim - 1], _pairs1[dim - 1], _pairsComp[dim - 1], _matches[dim - 1],
                            _isMatched0, _isMatched1);
        topDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
    }

    if (dim > 2)
    {
        InterDimensions interDim(cgc0, cgc1, cgcComp, config, _pairs0, _pairs1, _pairsComp, _matches, _isMatched0, _isMatched1);
        interDim.computePairsAndMatch(ctr0, ctr1, ctrComp);
    }

    Dimension0 dim0(cgc0, cgc1, cgcComp, config, _pairs0[0], _pairs1[0], _pairsComp[0], _matches[0], _isMatched0, _isMatched1);
    dim0.computePairsAndMatch(ctr0, ctr1, ctrComp);

    computeUnmatched();
}

void BettiMatchingND::computeUnmatched()
{
    _unmatchedPairs0 = vector<vector<Pair>>(dim);
    _unmatchedPairs1 = vector<vector<Pair>>(dim);

    for (int d = 0; d < dim; d++)
    {
        for (auto &pair : _pairs0[d])
        {
            if (!_isMatched0[cgc0.getCubeIndex(pair.birth)])
            {
                _unmatchedPairs0[d].push_back(pair);
            }
        }
        for (auto &pair : _pairs1[d])
        {
            if (!_isMatched1[cgc1.getCubeIndex(pair.birth)])
            {
                _unmatchedPairs1[d].push_back(pair);
            }
        }
    }
}
