#pragma once

#include "data_structures.h"

class BoundaryEnumerator {
	public:
	Cube nextFace;
	
	BoundaryEnumerator(const CubicalGridComplex* const _cgc);
	void setBoundaryEnumerator(const Cube& cube);
	bool hasNextFace();

	private:
	const CubicalGridComplex* const cgc;
	Cube cube;
	uint8_t position;
};
