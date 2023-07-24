#include "enumerators.h"

BoundaryEnumerator::BoundaryEnumerator(CubicalGridComplex* _cgc) {
	nextFace = Cube();
	cgc = _cgc;
}

void BoundaryEnumerator::setBoundaryEnumerator(Cube& _cube) {
	cube = _cube;
	position = 0; 
}

bool BoundaryEnumerator::hasNextFace() {
	if (position == 4) {
		return false;
	} else {
		auto x = cube.x;
		auto y = cube.y;
		auto z = cube.z;
		double birth;
		switch (cube.type) {
			case 0:
			switch (position) {
				case 0:
				birth = cgc->getBirth(x,y,z,1,1);
				nextFace = Cube(birth,x,y,z,1,1);
				break;
				case 1:
				birth = cgc->getBirth(x,y,z,2,1);
				nextFace = Cube(birth,x,y,z,2,1);
				break;
				case 2:
				birth = cgc->getBirth(x,y+1,z,2,1);
				nextFace = Cube(birth,x,y+1,z,2,1);
				break;
				case 3:
				birth = cgc->getBirth(x,y,z+1,1,1);
				nextFace = Cube(birth,x,y,z+1,1,1);
				break;
			}
			break;
			case 1:
			switch (position) {
				case 0:
				birth = cgc->getBirth(x,y,z,0,1);
				nextFace = Cube(birth,x,y,z,0,1);
				break;
				case 1:
				birth = cgc->getBirth(x,y,z,2,1);
				nextFace = Cube(birth,x,y,z,2,1);
				break;
				case 2:
				birth = cgc->getBirth(x+1,y,z,2,1);
				nextFace = Cube(birth,x+1,y,z,2,1);
				break;
				case 3:
				birth = cgc->getBirth(x,y,z+1,0,1);
				nextFace = Cube(birth,x,y,z+1,0,1);
				break;
			}
			break;
			case 2:
			switch (position) {
				case 0:
				birth = cgc->getBirth(x,y,z,0,1);
				nextFace = Cube(birth,x,y,z,0,1);
				break;
				case 1:
				birth = cgc->getBirth(x,y,z,1,1);
				nextFace = Cube(birth,x,y,z,1,1);
				break;
				case 2:
				birth = cgc->getBirth(x+1,y,z,1,1);
				nextFace = Cube(birth,x+1,y,z,1,1);
				break;
				case 3:
				birth = cgc->getBirth(x,y+1,z,0,1);
				nextFace = Cube(birth,x,y+1,z,0,1);
				break;
			}
			break;
		}
		position ++;
		return true;
	}	
}
