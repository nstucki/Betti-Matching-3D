#include "enumerators.h"

BoundaryEnumerator::BoundaryEnumerator(const CubicalGridComplex* const _cgc) : cgc(_cgc) { nextFace = Cube(); }

void BoundaryEnumerator::setBoundaryEnumerator(const Cube& _cube) {
	cube = _cube;
	position = 0; 
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
				#ifdef USE_EMERGENT_PAIRS
				case 0:
				birth = cgc->getBirth(x, y+1, z, 2, 1);
				nextFace = Cube(birth, x, y+1, z, 2);
				break;

				case 1:
				birth = cgc->getBirth(x, y, z+1, 1, 1);
				nextFace = Cube(birth, x, y, z+1, 1);
				break;

				case 2:
				birth = cgc->getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 3:
				birth = cgc->getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;
				#else
				case 0:
				birth = cgc->getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;

				case 1:
				birth = cgc->getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 2:
				birth = cgc->getBirth(x, y, z+1, 1, 1);
				nextFace = Cube(birth, x, y, z+1, 1);
				break;

				case 3:
				birth = cgc->getBirth(x, y+1, z, 2, 1);
				nextFace = Cube(birth, x, y+1, z, 2);
				break;
				#endif
			}
			break;

			case 1:
			switch (position) {
				#ifdef USE_EMERGENT_PAIRS
				case 0:
				birth = cgc->getBirth(x+1, y, z, 2, 1);
				nextFace = Cube(birth, x+1, y, z, 2);
				break;

				case 1:
				birth = cgc->getBirth(x, y, z+1, 0, 1);
				nextFace = Cube(birth, x, y, z+1, 0);
				break;

				case 2:
				birth = cgc->getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 3:
				birth = cgc->getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;
				#else
				case 0:
				birth = cgc->getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;

				case 1:
				birth = cgc->getBirth(x, y, z, 2, 1);
				nextFace = Cube(birth, x, y, z, 2);
				break;

				case 2:
				birth = cgc->getBirth(x, y, z+1, 0, 1);
				nextFace = Cube(birth, x, y, z+1, 0);
				break;

				case 3:
				birth = cgc->getBirth(x+1, y, z, 2, 1);
				nextFace = Cube(birth, x+1, y, z, 2);
				break;
				#endif
			}
			break;
			
			case 2:
			switch (position) {
				#ifdef USE_EMERGENT_PAIRS
				case 0:
				birth = cgc->getBirth(x+1, y, z, 1, 1);
				nextFace = Cube(birth, x+1, y, z, 1);
				break;

				case 1:
				birth = cgc->getBirth(x, y+1, z, 0, 1);
				nextFace = Cube(birth, x, y+1, z, 0);
				break;

				case 2:
				birth = cgc->getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;

				case 3:
				birth = cgc->getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;
				#else
				case 0:
				birth = cgc->getBirth(x, y, z, 0, 1);
				nextFace = Cube(birth, x, y, z, 0);
				break;

				case 1:
				birth = cgc->getBirth(x, y, z, 1, 1);
				nextFace = Cube(birth, x, y, z, 1);
				break;

				case 2:
				birth = cgc->getBirth(x, y+1, z, 0, 1);
				nextFace = Cube(birth, x, y+1, z, 0);
				break;

				case 3:
				birth = cgc->getBirth(x+1, y, z, 1, 1);
				nextFace = Cube(birth, x+1, y, z, 1);
				break;
				#endif
			}
			break;
		}
		position ++;
		return true;
	}	
}