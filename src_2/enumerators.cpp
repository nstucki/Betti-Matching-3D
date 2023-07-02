#include "enumerators.h"

#include <iostream>


CubeEnumerator::CubeEnumerator(const CubicalGridComplex& _cgc, const uint64_t _dim) : cgc(_cgc), dim(_dim) {
    degenMax = cgc.dim-dim-1;
    nonDegenMax = cgc.dim-1;
    nonDegen = vector<bool>(cgc.dim, false);

    vector<uint64_t> coordinates(cgc.dim, 0);
    for (uint64_t i = cgc.dim; i-- > cgc.dim-dim;) {
        coordinates[i] = 1;
        nonDegen[i] = true;
    }
    float birth = cgc.getBirth(coordinates);
    nextCube = Cube(birth, coordinates);
}

bool CubeEnumerator::hasNextCube() {
    uint64_t axis = cgc.dim;

    while (axis > 0) {
        --axis;
        if (canIncreaseByOne(axis)) {
            IncreaseByOne(axis);
            return true;
        } else if (nextCube.coordinates[axis]+2 < 2*cgc.shape[axis]-1) {
            IncreaseByTwo(axis);
            return true;
        } else { continue; }
    }

    return false;
}

bool CubeEnumerator::canIncreaseByOne(const uint64_t& axis) const {
    if (nextCube.coordinates[axis]+1 == 2*cgc.shape[axis]-1) { return false; }
    else if (nonDegen[axis] &&  axis < degenMax) { return true; }
    else if ( !nonDegen[axis] && axis < nonDegenMax) { return true; }
    else { return false; }
}

bool CubeEnumerator::canIncreaseByTwo(const uint64_t& axis) const {
    if (nextCube.coordinates[axis]+2 >= 2*cgc.shape[axis]-1) { return false; }
    else { return true; }
}

void CubeEnumerator::IncreaseByOne(const uint64_t& axis) {
    nextCube.coordinates[axis] += 1;
    nonDegen[axis] = !nonDegen[axis];
    if (nonDegen[axis]) {
        nonDegen[degenMax+1] = false;
        degenMax += 1;
        if (degenMax == nonDegenMax) { nonDegenMax = axis; }
    } else {
        nonDegen[degenMax] = true;
         if (axis == nonDegenMax) { nonDegenMax = degenMax; }
        degenMax -= 1;
    }
    for (uint64_t a = axis+1; a < cgc.dim; a++) {
        nextCube.coordinates[a] = nonDegen[a];
    }
    nextCube.birth = cgc.getBirth(nextCube.coordinates);
}

void CubeEnumerator::IncreaseByTwo(const uint64_t& axis) {
    nextCube.coordinates[axis] += 2;
    for (uint64_t a = axis+1; a < cgc.dim; a++) {
        nextCube.coordinates[a] = nonDegen[a];
    }
    nextCube.birth = cgc.getBirth(nextCube.coordinates);
}

Cube CubeEnumerator::getNextCube() const { return nextCube; }

uint64_t CubeEnumerator::getNumberOfCubes() const {
    vector<vector<bool>> subsets = getSubsets(cgc.dim, dim);
    uint64_t numCubes = 0;
    uint64_t numCubesForSubset;
    for (auto& subset : subsets) {
        numCubesForSubset = 1;
        for (uint64_t i = 0; i < cgc.dim; i++) {
            numCubesForSubset *= (cgc.shape[i] - subset[i]);
        }
        numCubes += numCubesForSubset;
    }
    return numCubes;
}


BoundaryEnumerator::BoundaryEnumerator(const CubicalGridComplex& _cgc) : cgc(_cgc) {}

void BoundaryEnumerator::setBoundaryEnumerator(const Cube& _cube, uint64_t _dim) {
	cube = _cube;
    dim = _dim;
    nonDegenAxes.clear();
    nonDegenAxes.reserve(dim);
    uint64_t counter = 0;
    for (uint64_t i = 0; i < cgc.dim; i++) {
        if (cube.coordinates[i]%2 != 0) {
            nonDegenAxes.push_back(i);
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
        coordinates[nonDegenAxes[position]] += shift;
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


DualEdgeEnumerator::DualEdgeEnumerator(const CubicalGridComplex& _cgc) : cgc(_cgc) {
    vector<uint64_t> coordinates(cgc.dim, 1);
    coordinates[0] = 0;
    degenAxis = 0;
    float birth = cgc.getBirth(coordinates);
    nextCube = Cube(birth, coordinates);
}

bool DualEdgeEnumerator::hasNextCube() {
    uint64_t axis = cgc.dim-1;

    while (axis <= cgc.dim-1) {
        if ((axis > degenAxis) && (nextCube.coordinates[axis]+2 < 2*cgc.shape[axis]-1)) {
            nextCube.coordinates[axis] += 2;
            for (uint64_t a = axis+1; a < cgc.dim; a++) {
                nextCube.coordinates[a] = 1;
            }
            break;
        }
        if ((axis <= degenAxis) && (nextCube.coordinates[axis]+1 < 2*cgc.shape[axis]-1)) {     
            if (axis == degenAxis) {
                if (axis == cgc.dim-1) {
                    nextCube.coordinates[cgc.dim-1] += 2;
                } else {
                    nextCube.coordinates[axis] += 1;
                    nextCube.coordinates[axis+1] = 0;
                    degenAxis = axis+1;
                    for (uint64_t a = axis+2; a < cgc.dim; a++) {
                        nextCube.coordinates[a] = 1;
                    }
                }           
            } else {
                nextCube.coordinates[axis] += 1;
                degenAxis = axis;
                for (uint64_t a = axis+1; a < cgc.dim; a++) {
                    nextCube.coordinates[a] = 1;
                }
            }
            break;
        }
        --axis;
    }
    if (axis <= cgc.dim-1) {
        nextCube.birth = cgc.getBirth(nextCube.coordinates);
        return true;
    } else {
        return false;
    }
}

Cube DualEdgeEnumerator::getNextCube() const { return nextCube; }