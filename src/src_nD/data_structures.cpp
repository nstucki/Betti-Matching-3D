#include "data_structures.h"
#include "utils.h"
#include "template_functions.h"

#include <iostream>
#include <functional>
#include <algorithm>

using namespace std;


Cube::Cube() {
	birth = 0; 
	coordinates.push_back(NONE);
}

Cube::Cube(value_t _birth, vector<index_t>_coordinates) : birth(_birth), coordinates(_coordinates) {}

Cube::Cube(const Cube &cube) {
    birth = cube.birth;
    coordinates = cube.coordinates;
}

bool Cube::operator==(const Cube& rhs) const { return (birth == rhs.birth && coordinates == rhs.coordinates); }

bool Cube::operator!=(const Cube& rhs) const { return (birth != rhs.birth || coordinates != rhs.coordinates); }

void Cube::print() const {
    cout << "(" << birth << ";";
    for (auto c = coordinates.cbegin(); c != coordinates.cend(); ++c) {
        cout << *c;
        if (c != coordinates.cend() - 1) {
            cout << ",";
        }
    }
    cout << ")";
}


bool CubeComparator::operator()(const Cube& cube1, const Cube& cube2) const {
    if (cube1.birth == cube2.birth) {
        for (index_t i = 0; i < cube1.coordinates.size(); i++) {
            if (cube1.coordinates[i] == cube2.coordinates[i]) {
                continue;
            } else { return cube1.coordinates[i] < cube2.coordinates[i]; } 
        }
        return false;
    } else { return cube1.birth < cube2.birth; }
}

Pair::Pair() {}

Pair::Pair(const Cube& _birth, const Cube& _death) : birth(_birth), death(_death) {}

Pair::Pair(const Pair& pair) : birth(pair.birth), death(pair.death) {}

bool Pair::operator==(const Pair& rhs) const { return (birth == rhs.birth && death == rhs.death); }

void Pair::print() const { cout << "("; birth.print(); cout << ","; death.print(); cout << ")"; }


Match::Match(const Pair &_pair0, const Pair &_pair1) : pair0(_pair0), pair1(_pair1) {}

void Match::print() const {
    cout << "matched "; pair0.print(); cout << " with "; pair1.print(); cout << endl;
}


CubicalGridComplex::CubicalGridComplex(vector<value_t> _image, const vector<index_t>& _shape) : 
image(_image), shape(_shape), dim(_shape.size()) {
	pixelCoordFactor = multiplyFromRightExclusively(_shape);
	vector<index_t> cubeGridShape = addToContainerElementwise(multiplyContainerElementwise(_shape, 2), -1);
	cubeCoordFactor = multiplyFromRightExclusively(cubeGridShape);
}

index_t CubicalGridComplex::getNumberOfCubes(index_t _dim) const {
    vector<vector<bool>> degenCoords = getSubsets(dim, _dim);
    index_t numCubes = 0;
    index_t numCubesForDegenCoord;
    for (auto& degenCoord : degenCoords) {
        numCubesForDegenCoord = 1;
        for (index_t i = 0; i < dim; i++) {
            numCubesForDegenCoord *= (shape[i] - degenCoord[i]);
        }
        numCubes += numCubesForDegenCoord;
    }
    return numCubes;
}

index_t CubicalGridComplex::getCubeIndex(const Cube& cube) const {
	index_t idx = 0;
	for (index_t i = dim; i-- > 0;) {
		idx += cube.coordinates[i]*cubeCoordFactor[i];
	}
	return idx;
}

index_t CubicalGridComplex::getCubeIndex(const vector<index_t>& coordinates) const {
	index_t idx = 0;
	for (index_t i = dim; i-- > 0;) {
		idx += coordinates[i]*cubeCoordFactor[i];
	}
	return idx;
}

vector<index_t> CubicalGridComplex::getCubeCoordinates(index_t idx) const {
	vector<index_t> coordinates;
	coordinates.reserve(dim);
	index_t c;
	for (index_t i = dim-1; i > 0; i--) {
		c = idx % (2*shape[i-1]-1);
		coordinates.push_back(c);
		idx = (idx-c)/(2*shape[i-1]-1);
	}
	coordinates.push_back(idx);
	std::reverse(coordinates.begin(), coordinates.end());
	return coordinates;
}

value_t CubicalGridComplex::getBirth(const vector<index_t>& coordinatesCube) const {
	index_t idx = getIndex(divideContainerElementwise(coordinatesCube, 2));
	value_t birth = image[idx];
	vector<index_t> indices{idx};
	index_t newIdx;
	for (index_t axis = dim; axis-- > 0;) {
		if (coordinatesCube[axis]%2 == 1) {
			vector<index_t>::size_type size = indices.size();
			for (std::vector<index_t>::size_type i = 0; i < size; ++i) {
				newIdx = indices[i]+pixelCoordFactor[axis];
	 			birth = max(birth, image[newIdx]);
	 			indices.push_back(newIdx);
			}
		}
	}
	return birth;
}

vector<index_t> CubicalGridComplex::getParentVoxel(const Cube &c) const
{
	auto workingCoordinates(c.coordinates);
	auto birth = this->getBirth(c.coordinates);
	std::function<bool(int)> getParentVoxelRecursively;
	getParentVoxelRecursively = [&workingCoordinates, &getParentVoxelRecursively, this, birth](int dim)
	{
		index_t coord = workingCoordinates[dim];
		if (dim < workingCoordinates.size())
		{
			if (coord % 2 == 0)
			{
				return getParentVoxelRecursively(dim + 1);
			}
			else
			{
				workingCoordinates[dim] = coord - 1;
				if (getParentVoxelRecursively(dim + 1))
				{
					return true;
				}
				workingCoordinates[dim] = coord + 1;
				return getParentVoxelRecursively(dim + 1);
			}
		}
		return this->getBirth(workingCoordinates) == birth;
	};
	if (!getParentVoxelRecursively(0))
	{
		throw std::runtime_error("no adjacent voxel with same birth found");
	}
	vector<index_t> pixelCoordinates;
	std::transform(std::begin(workingCoordinates), std::end(workingCoordinates), back_inserter(pixelCoordinates), [](auto coordinate)
				   { return coordinate / 2; });
	return pixelCoordinates;
}

void CubicalGridComplex::printImage() const {
    value_t birth;
	if (dim == 2) {
		for (index_t i = 0; i < shape[0]; i++) {
			for (index_t j = 0; j < shape[1]; j++) {
				birth = getValue(vector<index_t> {i,j});
				if (birth < 10){
					cout << ' ' << birth << ' ';
				}else{
					cout << birth << ' ';
				}            
			}
			cout << endl;
		}
	} else if (dim == 3) {
		for (index_t i = 0; i < shape[0]; i++) {
			for (index_t j = 0; j < shape[1]; j++) {
				for (index_t k = 0; k < shape[2]; k++) {
					birth = getValue(vector<index_t> {i,j,k});
					if (birth < 10){
						cout << ' ' << birth << ' ';
					}else{
						cout << birth << ' ';
					}            
				}
				cout << endl;
			}
			cout << endl;
		}
	}
}

void CubicalGridComplex::printCubes() const {
    value_t birth;
	if (dim == 2) {
		for (index_t i = 0; i < 2*shape[0]-1; i++) {
			for (index_t j = 0; j < 2*shape[1]-1; j++) {
				birth = getBirth({i,j});
				if (birth < 10){
					cout << ' ' << birth << ' ';
				}else{
					cout << birth << ' ';
				}            
			}
			cout << endl;
		}
	} else if (dim == 3) {
		for (index_t i = 0; i < 2*shape[0]-1; i++) {
			for (index_t j = 0; j < 2*shape[1]-1; j++) {
				for (index_t k = 0; k < 2*shape[2]-1; k++) {
					birth = getBirth({i,j,k});
					if (birth < 10){
						cout << ' ' << birth << ' ';
					}else{
						cout << birth << ' ';
					}            
				}
				cout << endl;
			}
			cout << endl;
		}
	}
}

index_t CubicalGridComplex::getIndex(const vector<index_t>& pixelCoordinates) const {
	index_t idx = 0;
	for (index_t i = dim; i-- >0;) {
		idx += pixelCoordinates[i]*pixelCoordFactor[i];
	}
	return idx;
}

value_t CubicalGridComplex::getValue(const vector<index_t>& pixelCoordinates) const { return image[getIndex(pixelCoordinates)]; }


UnionFind::UnionFind(const CubicalGridComplex& _cgc) : cgc(_cgc) {
	index_t n = cgc.getNumberOfCubes(0);
	parent.reserve(n);
	birthtime.reserve(n);
	for (index_t i = 0; i < n; i++) {
		parent.push_back(i);
		birthtime.push_back(cgc.getBirth(getCoordinates(i)));
	}
}

index_t UnionFind::find(index_t x) {
	index_t y = x, z = parent[y];
	while (z != y) {
		y = z;
		z = parent[y];
	}
	y = parent[x];
	while (z != y) {
		parent[x] = z;
		x = y;
		y = parent[x];
	}
	return z;
}

index_t UnionFind::link(const index_t& x, const index_t& y){
	if (birthtime[x] > birthtime[y]){
		parent[x] = y; 
		return x;
	} else if (birthtime[x] < birthtime[y]){
		parent[y] = x;
		return y;
	} else {
		if (x > y){
			parent[x] = y;
			return x;
		} else {
			parent[y] = x;
			return y;
		}
	}
}

value_t UnionFind::getBirth(const index_t& idx) const { return birthtime[idx]; }

index_t UnionFind::getIndex(const vector<index_t>& coordinates) const {
	index_t idx = 0;
	index_t factor = 1;
	for (index_t d = cgc.dim; d-- > 0;) {
		idx += coordinates[d]/2 * factor;
		factor *= cgc.shape[d];
	}
	return idx;
}

vector<index_t> UnionFind::getCoordinates(index_t idx) const {
	vector<index_t> coordinates(cgc.dim, 0);
	index_t remainder;
	for (index_t d = cgc.dim; d-- > 0;) {
		remainder = idx%cgc.shape[d];
		idx -= remainder;
		idx /= cgc.shape[d];
		coordinates[d] = 2*remainder;
	}
	return coordinates;
}

void UnionFind::reset() { for (index_t i = 0; i < parent.size(); i++) { parent[i] = i; } }


UnionFindDual::UnionFindDual(const CubicalGridComplex &_cgc) : cgc(_cgc) {
	star = cgc.getNumberOfCubes(cgc.dim);
	parent.reserve(star+1);
	birthtime.reserve(star+1);
	for (index_t i = 0; i < star; i++) {
		parent.push_back(i);
		birthtime.push_back(cgc.getBirth(getCoordinates(i)));
	}
	parent.push_back(star);
	birthtime.push_back(numeric_limits<value_t>::infinity());
}

index_t UnionFindDual::find(index_t x) {
	index_t y = x, z = parent[y];
	while (z != y) {
		y = z;
		z = parent[y];
	}
	y = parent[x];
	while (z != y) {
		parent[x] = z;
		x = y;
		y = parent[x];
	}
	return z;
}

index_t UnionFindDual::link(const index_t& x, const index_t& y){
	if (birthtime[x] < birthtime[y]){
		parent[x] = y; 
		return x;
	} else if (birthtime[x] > birthtime[y]){
		parent[y] = x;
		return y;
	} else {
		if (x < y){
			parent[x] = y;
			return x;
		} else {
			parent[y] = x;
			return y;
		}
	}
}

value_t UnionFindDual::getBirth(const index_t& idx) const { return birthtime[idx]; }

index_t UnionFindDual::getIndex(const vector<index_t>& coordinates) const {
	index_t idx = 0;
	index_t factor = 1;
	for (index_t d = cgc.dim; d-- > 0;) {
		idx += ((coordinates[d]-1)/2 * factor);
		factor *= (cgc.shape[d]-1);
	}
	return idx;
}

vector<index_t> UnionFindDual::getCoordinates(index_t idx) const {
	vector<index_t> coordinates(cgc.dim, 0);
	index_t remainder;
	for (index_t d = cgc.dim; d-- > 0;) {
		remainder = idx%(cgc.shape[d]-1);
		idx -= remainder;
		idx /= (cgc.shape[d]-1);
		coordinates[d] = 2*remainder + 1;
	}
	return coordinates;
}

void UnionFindDual::reset() { for (index_t i = 0; i < parent.size(); i++) { parent[i] = i; } }