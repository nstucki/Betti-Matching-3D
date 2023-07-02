#pragma once

#include "data_structures.h"
#include "utils.h"


class CubeEnumerator{
    private:
    const CubicalGridComplex& cgc;
    const uint64_t dim;
    Cube nextCube;
    vector<bool> nonDegen;
    uint64_t nonDegenMax;
    uint64_t degenMax;
    
    bool canIncreaseByOne(const uint64_t& axis) const;
    bool canIncreaseByTwo(const uint64_t& axis) const;
    void IncreaseByOne(const uint64_t& axis);
    void IncreaseByTwo(const uint64_t& axis);
    

    public:
    CubeEnumerator(const CubicalGridComplex& cgc, const uint64_t dim);
    bool hasNextCube();
    Cube getNextCube() const;
};


class BoundaryEnumerator {
private:
	const CubicalGridComplex& cgc;
	Cube cube;
    uint64_t dim;
    Cube nextFace;
    vector<uint64_t> nonDegenAxes;
    uint64_t position;
    int8_t shift;

public:
	BoundaryEnumerator(const CubicalGridComplex& cgc);
	void setBoundaryEnumerator(const Cube& cube, uint64_t dim);
	bool hasNextFace();
    Cube getNextFace() const;
};


class DualEdgeEnumerator{
    private:
    const CubicalGridComplex& cgc;
    Cube nextCube;
    uint64_t degenAxis;

    public:
    DualEdgeEnumerator(const CubicalGridComplex& cgc);
    bool hasNextCube();
    Cube getNextCube() const;
};