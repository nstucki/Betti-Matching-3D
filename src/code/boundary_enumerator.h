#pragma once

#include "cubical_grid_complex.h"

class BoundaryEnumerator {
private:
	CubicalGridComplex* cgc;
	Cube cube;
	uint8_t position;
public:
	Cube nextFace;

	BoundaryEnumerator(CubicalGridComplex* _cgc);
	void setBoundaryEnumerator(Cube& _cube);
	bool hasNextFace();
};
