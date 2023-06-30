#pragma once

#include "data_structures.h"


class DualEdgeEnumerator{

    private:
    const CubicalGridComplex& cgc;
    Cube cube;
    uint64_t degen;

    public:
    DualEdgeEnumerator(const CubicalGridComplex& cgc);
    bool hasNextCube();
    Cube getCube() const;
};


class BoundaryEnumerator {
private:
	const CubicalGridComplex& cgc;
	Cube cube;
    uint64_t dim;
    Cube nextFace;
    vector<uint64_t> nonDegen;
    uint64_t position;
    int8_t shift;

public:
	BoundaryEnumerator(const CubicalGridComplex& cgc);
	void setBoundaryEnumerator(const Cube& cube, uint64_t dim);
	bool hasNextFace();
    Cube getNextFace() const;
};