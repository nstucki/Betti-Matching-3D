#pragma once

#include "data_structures.h"


namespace dim3 {
class BoundaryEnumerator {
	public:
	Cube nextFace;
	
	BoundaryEnumerator(const CubicalGridComplex& cgc);
	void setBoundaryEnumerator(const Cube& cube);
	bool hasPreviousFace();
	bool hasNextFace();

	private:
	const CubicalGridComplex& cgc;
	Cube cube;
	uint8_t position;
};


class CoboundaryEnumerator {
	public:
	Cube nextCoface;
	
	CoboundaryEnumerator(const CubicalGridComplex& cgc);
	void setCoboundaryEnumerator(const Cube& cube);
	bool hasNextCoface();

	private:
	const CubicalGridComplex& cgc;
	Cube cube;
	uint8_t position;
};
}