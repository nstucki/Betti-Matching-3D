#include "enumerators.h"

#include <iostream>

using namespace dimN;

CubeEnumerator::CubeEnumerator(const CubicalGridComplex &_cgc,
                               const index_t _dim)
    : cgc(_cgc) {
    degenMax = cgc.dim - _dim - 1;
    nonDegenMax = cgc.dim - 1;
    nonDegen = vector<bool>(cgc.dim, false);
    vector<index_t> coordinates(cgc.dim, 0);
    for (index_t i = cgc.dim; i-- > cgc.dim - _dim;) {
        coordinates[i] = 1;
        nonDegen[i] = true;
    }
    value_t birth = cgc.getBirth(coordinates);
    nextCube = Cube(birth, coordinates);
}

Cube CubeEnumerator::getNextCube() const { return nextCube; }

bool CubeEnumerator::hasNextCube() {
    index_t axis = cgc.dim;
    while (axis > 0) {
        --axis;
        if (canIncreaseByOne(axis)) {
            IncreaseByOne(axis);
            return true;
        } else if (nextCube.coordinates[axis] + 2 < 2 * cgc.shape[axis] - 1) {
            IncreaseByTwo(axis);
            return true;
        } else {
            continue;
        }
    }
    return false;
}

bool CubeEnumerator::canIncreaseByOne(const index_t &axis) const {
    if (nextCube.coordinates[axis] + 1 == 2 * cgc.shape[axis] - 1) {
        return false;
    } else if (nonDegen[axis] && axis < degenMax) {
        return true;
    } else if (!nonDegen[axis] && axis < nonDegenMax) {
        return true;
    } else {
        return false;
    }
}

bool CubeEnumerator::canIncreaseByTwo(const index_t &axis) const {
    if (nextCube.coordinates[axis] + 2 >= 2 * cgc.shape[axis] - 1) {
        return false;
    } else {
        return true;
    }
}

void CubeEnumerator::IncreaseByOne(const index_t &axis) {
    nextCube.coordinates[axis] += 1;
    nonDegen[axis] = !nonDegen[axis];
    if (nonDegen[axis]) {
        nonDegen[degenMax + 1] = false;
        degenMax += 1;
        if (degenMax == nonDegenMax) {
            nonDegenMax = axis;
        }
    } else {
        nonDegen[degenMax] = true;
        if (axis == nonDegenMax) {
            nonDegenMax = degenMax;
        }
        degenMax -= 1;
    }
    for (size_t a = axis + 1; a < cgc.dim; ++a) {
        nextCube.coordinates[a] = nonDegen[a];
    }
    nextCube.birth = cgc.getBirth(nextCube.coordinates);
}

void CubeEnumerator::IncreaseByTwo(const index_t &axis) {
    nextCube.coordinates[axis] += 2;
    for (size_t a = axis + 1; a < cgc.dim; ++a) {
        nextCube.coordinates[a] = nonDegen[a];
    }
    nextCube.birth = cgc.getBirth(nextCube.coordinates);
}

BoundaryEnumerator::BoundaryEnumerator(const CubicalGridComplex &_cgc)
    : cgc(_cgc) {}

void BoundaryEnumerator::setBoundaryEnumerator(const Cube &_cube,
                                               const size_t &_dim) {
    cube = _cube;
    dim = _dim;
    nonDegenAxes.clear();
    nonDegenAxes.reserve(dim);
    size_t counter = 0;
    for (size_t axis = 0; axis < cgc.dim; ++axis) {
        if (cube.coordinates[axis] % 2 != 0) {
            nonDegenAxes.push_back(axis);
            if (++counter == dim) {
                break;
            }
        }
    }
    position = 0;
#ifdef USE_EMERGENT_PAIRS
    shift = 1;
#else
    shift = -1;
#endif
}

bool BoundaryEnumerator::hasNextFace() {
    if (position == -1) {
        return false;
    } else {
        vector<index_t> coordinates = cube.coordinates;
        coordinates[nonDegenAxes[position]] += shift;
        nextFace = Cube(cgc.getBirth(coordinates), coordinates);
#ifdef USE_EMERGENT_PAIRS
        if (position == dim - 1 && shift == 1) {
            shift = -1;
        } else if (shift == 1) {
            ++position;
        }
#else
        if (position == dim - 1 && shift == -1) {
            shift = 1;
        } else if (shift == -1) {
            ++position;
        }
#endif
        else {
            --position;
        }
        return true;
    }
}
