#pragma once

#include "data_structures.h"

class BoundaryEnumerator {
	public:
	BoundaryEnumerator(const CubicalGridComplex* const _cgc);
	void setBoundaryEnumerator(const Cube& cube);
	bool hasNextFace();
	Cube getNextFace() const;

	private:
	const CubicalGridComplex* const cgc;
	Cube cube;
	Cube nextFace;
	uint8_t position;
};
