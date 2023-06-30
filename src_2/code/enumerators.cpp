#include "enumerators.h"

#include <iostream>


DualEdgeEnumerator::DualEdgeEnumerator(const CubicalGridComplex& _cgc) : cgc(_cgc) {
    vector<uint64_t> coordinates(cgc.shape.size(), 1);
    coordinates[0] = 0;
    float birth = cgc.getBirth(coordinates);
    cube = Cube(birth, coordinates);
    degen = 0;
}

bool DualEdgeEnumerator::hasNextCube() {
    uint64_t axis = cgc.shape.size()-1;
    while (axis <= cgc.shape.size()-1) {
        if ((axis > degen) && (cube.coordinates[axis]+2 < 2*cgc.shape[axis]-1)) {
            cube.coordinates[axis] += 2;
            for (uint64_t a = axis+1; a < cgc.shape.size(); a++) {
                cube.coordinates[a] = 1;
            }
            break;
        }
        if ((axis <= degen) && (cube.coordinates[axis]+1 < 2*cgc.shape[axis]-1)) {     
            if (axis == degen) {
                if (axis == cgc.shape.size()-1) {
                    cube.coordinates[cgc.shape.size()-1] += 2;
                } else {
                    cube.coordinates[axis] += 1;
                    cube.coordinates[axis+1] = 0;
                    degen = axis+1;
                    for (uint64_t a = axis+2; a < cgc.shape.size(); a++) {
                        cube.coordinates[a] = 1;
                    }
                }           
            } else {
                cube.coordinates[axis] += 1;
                degen = axis;
                for (uint64_t a = axis+1; a < cgc.shape.size(); a++) {
                    cube.coordinates[a] = 1;
                }
            }
            break;
        }
        --axis;
    }
    if (axis <= cgc.shape.size()-1) {
        cube.birth = cgc.getBirth(cube.coordinates);
        return true;
    } else {
        return false;
    }
}

Cube DualEdgeEnumerator::getCube() const { return cube; }



BoundaryEnumerator::BoundaryEnumerator(const CubicalGridComplex& _cgc) : cgc(_cgc) {}


void BoundaryEnumerator::setBoundaryEnumerator(const Cube& _cube, uint64_t _dim) {
	cube = _cube;
    dim = _dim;
    nonDegen.clear();
    nonDegen.reserve(dim);
    uint64_t counter = 0;
    for (uint64_t i = 0; i < cgc.dim; i++) {
        if (cube.coordinates[i]%2 != 0) {
            nonDegen.push_back(i);
            counter++;
            if (counter == dim) { break; }
        }
    }
    position = 0;
    shift = -1;
}

bool BoundaryEnumerator::hasNextFace() {
	if ( position == -1 ) { return false; } 
    else {
        vector<uint64_t> coordinates = cube.coordinates;
        coordinates[nonDegen[position]] += shift;
        float birth = cgc.getBirth(coordinates);
        nextFace = Cube(birth, coordinates);
        if (position == dim-1 && shift == -1) {
            shift = 1;
        } else if (shift == -1) {
            position++;
        } else {
            position--;
        }
        return true;
    }
}


Cube BoundaryEnumerator::getNextFace() const { return nextFace; }