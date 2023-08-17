#pragma once

#include "data_structures.h"

class BoundaryEnumerator {
	public:
	Cube nextFace;
	
	BoundaryEnumerator(const CubicalGridComplex& cgc);
	void setBoundaryEnumerator(const Cube& cube);
	bool hasNextFace();

	private:
	const CubicalGridComplex& cgc;
	Cube cube;
	uint8_t position;
};
