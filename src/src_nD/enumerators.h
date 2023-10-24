#pragma once

#include "data_structures.h"
#include "../utils.h"


namespace dimN {
class CubeEnumerator {
    public:
    CubeEnumerator(const CubicalGridComplex& cgc, const index_t dim);
    bool hasNextCube();
    Cube getNextCube() const;

    private:
    const CubicalGridComplex& cgc;
    Cube nextCube;
    vector<bool> nonDegen;
    index_t nonDegenMax;
    index_t degenMax;
    
    bool canIncreaseByOne(const index_t& axis) const;
    bool canIncreaseByTwo(const index_t& axis) const;
    void IncreaseByOne(const index_t& axis);
    void IncreaseByTwo(const index_t& axis);
};


class BoundaryEnumerator {
    public:
    Cube nextFace;

    BoundaryEnumerator(const CubicalGridComplex& cgc);
    void setBoundaryEnumerator(const Cube& cube, const size_t& dim);
    bool hasNextFace();
    
    private:
    const CubicalGridComplex& cgc;
    Cube cube;
    index_t dim;
    vector<index_t> nonDegenAxes;
    index_t position;
    int8_t shift;
};
}