#include "enumerators.h"

#include <iostream>
using namespace std;

using namespace dim3;


BoundaryEnumerator::BoundaryEnumerator(const CubicalGridComplex& _cgc) : cgc(_cgc) { nextFace = Cube(); }

void BoundaryEnumerator::setBoundaryEnumerator(const Cube& _cube) {
	cube = _cube;
	position = 0; 
}

bool BoundaryEnumerator::hasPreviousFace() {
	if (position == 4) { return false; } 
	else {
		index_t x = cube.x();
		index_t y = cube.y();
		index_t z = cube.z();
		value_t birth;
		switch(cube.type()) {
			case 0:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x, y+1, z, 2, 1);
				nextFace = Cube(birth, x, y+1, z, 2);
				break;

				case 1:
				birth = cgc.getBirth(x, y, z+1, 1, 1);
				nextFace = Cube(birth, x, y, z+1, 1);
				break;

				case 2:
				birth = cgc.getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 3:
				birth = cgc.getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;
			}
			break;

			case 1:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x+1, y, z, 2, 1);
				nextFace = Cube(birth, x+1, y, z, 2);
				break;

				case 1:
				birth = cgc.getBirth(x, y, z+1, 0, 1);
				nextFace = Cube(birth, x, y, z+1, 0);
				break;

				case 2:
				birth = cgc.getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 3:
				birth = cgc.getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;
			}
			break;
			
			case 2:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x+1, y, z, 1, 1);
				nextFace = Cube(birth, x+1, y, z, 1);
				break;

				case 1:
				birth = cgc.getBirth(x, y+1, z, 0, 1);
				nextFace = Cube(birth, x, y+1, z, 0);
				break;

				case 2:
				birth = cgc.getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;

				case 3:
				birth = cgc.getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;
			}
			break;
		}
		++position;
		return true;
	}	
}

bool BoundaryEnumerator::hasNextFace() {
	if (position == 4) { return false; } 
	else {
		index_t x = cube.x();
		index_t y = cube.y();
		index_t z = cube.z();
		value_t birth;
		switch(cube.type()) {
			case 0:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;

				case 1:
				birth = cgc.getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 2:
				birth = cgc.getBirth(x, y, z+1, 1, 1);
				nextFace = Cube(birth, x, y, z+1, 1);
				break;

				case 3:
				birth = cgc.getBirth(x, y+1, z, 2, 1);
				nextFace = Cube(birth, x, y+1, z, 2);
				break;
			}
			break;

			case 1:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;

				case 1:
				birth = cgc.getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 2:
				birth = cgc.getBirth(x, y, z+1, 0, 1);
				nextFace = Cube(birth, x, y, z+1, 0);
				break;

				case 3:
				birth = cgc.getBirth(x+1, y, z, 2, 1);
				nextFace = Cube(birth, x+1, y, z, 2);
				break;
			}
			break;
			
			case 2:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;

				case 1:
				birth = cgc.getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;

				case 2:
				birth = cgc.getBirth(x, y+1, z, 0, 1);
				nextFace = Cube(birth, x, y+1, z, 0);
				break;

				case 3:
				birth = cgc.getBirth(x+1, y, z, 1, 1);
				nextFace = Cube(birth, x+1, y, z, 1);
				break;
			}
			break;
		}
		++position;
		return true;
	}	
}


CoboundaryEnumerator::CoboundaryEnumerator(const CubicalGridComplex& _cgc) : cgc(_cgc) { nextCoface = Cube(); }

void CoboundaryEnumerator::setCoboundaryEnumerator(const Cube& _cube) {
	cube = _cube;
	position = 0; 
}

bool CoboundaryEnumerator::hasNextCoface() {
	if (position == 4) { return false; } 
	else {
		index_t x = cube.x();
		index_t y = cube.y();
		index_t z = cube.z();
		value_t birth;
		switch(cube.type()) {
			case 0:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x, y-1, z, 2, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y-1, z, 2);
					break;
				} else { ++position; }

				case 1:
				birth = cgc.getBirth(x, y, z-1, 1, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z-1, 1);
					break;
				} else { ++position; }
				
				case 2:
				birth = cgc.getBirth(x, y, z, 1, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z, 1);
					break;
				} else { ++position; }

				case 3:
				birth = cgc.getBirth(x, y, z, 2, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z, 2);
					break;
				} else { ++position; }

				case 4:
				return false;
			}
			break;

			case 1:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x-1, y, z, 2, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x-1, y, z, 2);
					break;
				} else { ++position; }

				case 1:
				birth = cgc.getBirth(x, y, z-1, 0, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z-1, 0);
					break;
				} else { ++position; }

				case 2:
				birth = cgc.getBirth(x, y, z, 0, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z, 0);
					break;
				} else { ++position; }

				case 3:
				birth = cgc.getBirth(x, y, z, 2, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z, 2);
					break;
				} else { ++position; }

				case 4:
				return false;
			}
			break;
			
			case 2:
			switch (position) {
				case 0:
				birth = cgc.getBirth(x-1, y, z, 1, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x-1, y, z, 1);
					break;
				} else { ++position; }

				case 1:
				birth = cgc.getBirth(x, y-1, z, 0, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y-1, z, 0);
					break;
				} else { ++position; }

				case 2:
				birth = cgc.getBirth(x, y, z, 0, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z, 0);
					break;
				} else { ++position; }

				case 3:
				birth = cgc.getBirth(x, y, z, 1, 2);
				if (birth != INFTY) {
					nextCoface = Cube(birth, x, y, z, 1);
					break;
				} else { ++position; }

				case 4:
				return false;
			}
			break;
		}
		++position;
		return true;
	}	
}