#include "data_structures.h"
#include <iostream>
#include <algorithm>

using namespace dim2;
using namespace std;



Cube::Cube() : birth(0), index(NONE_INDEX) {}


Cube::Cube(const Cube& cube) : birth(cube.birth), index(cube.index) {}


Cube::Cube(value_t _birth, index_t _x, index_t _y, uint8_t _type) : birth(_birth) {
	index = ((uint64_t)_x << 34) | ((uint64_t)_y<<4) | (uint64_t)_type;
}


index_t Cube::x() const { return ((index >> 34) & 0xfffff); }


index_t Cube::y() const { return ((index >> 4) & 0xfffff); }


uint8_t Cube::type() const { return( index & 0xf); }


bool Cube::operator==(const Cube& rhs) const{ return (index == rhs.index); }


void Cube::print() const { cout << "(" << birth << "," << x() << "," << y() << "," << unsigned(type()) << ")"; }



bool CubeComparator::operator()(const Cube& cube1, const Cube& cube2) const{
	if (cube1.birth == cube2.birth) { return (cube1.index < cube2.index); } 
	else { return (cube1.birth < cube2.birth); }
}



Pair::Pair() {}


Pair::Pair(const Cube& _birth, const Cube& _death) : birth(_birth), death(_death) {}


Pair::Pair(const Pair& pair) : birth(pair.birth), death(pair.death) {}


bool Pair::operator==(const Pair &rhs) const { return (birth == rhs.birth && death == rhs.death); }


void Pair::print() const { cout << "("; birth.print(); cout << ";"; death.print(); cout << ")"; }



Match::Match(Pair _pair0, Pair _pair1) : pair0(_pair0), pair1(_pair1) {}


void Match::print() const { pair0.print(); cout << " <-> "; pair1.print(); cout << endl; }



CubicalGridComplex::CubicalGridComplex(const vector<value_t>& image, const vector<index_t>& _shape) : 
	shape(_shape), m_x(shape[0]-1), m_y(shape[1]-1), m_xy(m_x*m_y), n_xy(shape[0]*shape[1]) { getGridFromVector(image); }


CubicalGridComplex::CubicalGridComplex(CubicalGridComplex &&other) : 
	m_x(other.m_x), m_y(other.m_y), m_xy(other.m_xy), n_xy(other.n_xy), shape(std::move(other.shape)) {
	grid = other.grid;
	other.grid = nullptr;
}


CubicalGridComplex::~CubicalGridComplex() {
	if (grid != nullptr) {
		for (index_t i = 0; i < shape[0]+2; ++i) { delete[] grid[i]; }
		delete[] grid;
	}
}


size_t CubicalGridComplex::getNumberOfCubes(const uint8_t& dim) const {
	switch (dim) {
		case 0:
			return n_xy;

		case 1:
			return m_x*shape[1] + shape[0]*m_y;

		case 2:
			return m_xy;
	}
	cerr << "no cubes in dim " << unsigned(dim) << endl;
	return -1;
}


value_t CubicalGridComplex::getBirth(const index_t& x, const index_t&  y) const { return grid[x+1][y+1]; }


value_t CubicalGridComplex::getBirth(const index_t& x, const index_t& y, const uint8_t& type, const uint8_t& dim) const {
	switch (dim) {
		case 0:
			return getBirth(x, y);

		case 1:
			switch (type) {
				case 0:
					return max(getBirth(x, y), getBirth(x+1, y));

				case 1:
					return max(getBirth(x, y), getBirth(x, y+1));
			}
		
		case 2:
			return max({getBirth(x, y),getBirth(x, y+1),getBirth(x+1, y),getBirth(x+1, y+1)});
	}
	cerr << "birth not found!" << endl;
	return INFTY;
}


vector<index_t> CubicalGridComplex::getParentVoxel(const Cube& cube, const uint8_t& dim) const {
	index_t x = cube.x();
	index_t y = cube.y();
	switch (dim) {
		case 0:
			return {x,y};

		case 1:
			switch (cube.type()) {
				case 0:
					if (cube.birth == getBirth(x+1, y)) { return {x+1,y}; }
					else { return {x,y}; }

				case 1:
					if (cube.birth == getBirth(x, y+1)) { return {x,y+1}; }
					else { return {x,y}; }
			}
			
		case 2:
			if (cube.birth == getBirth(x+1, y+1)) { return {x+1,y+1}; }
			else if (cube.birth == getBirth(x+1, y)) { return {x+1,y}; }
			else if (cube.birth == getBirth(x, y+1)) { return {x,y+1}; }
			else { return {x,y}; }
	}
	cerr << "parent voxel not found!" << endl;
	return {0,0};
}


void CubicalGridComplex::printImage() const {
    value_t birth;
    for (index_t x = 0; x < shape[0]; ++x) {
        for (index_t y = 0; y < shape[1]; ++y) {
			birth = getBirth(x, y);
			if (birth < 10) { cout << ' ' << birth << ' '; }
			else { cout << birth << ' '; }            
        }
        cout << endl;
    }
	cout << endl;
}


void CubicalGridComplex::printRepresentativeCycle(const vector<vector<index_t>>& reprCycle) const {
	for (index_t x = 0; x < shape[0]; ++x) {
		for (index_t y = 0; y < shape[1]; ++y) {
			auto it = find(reprCycle.begin(), reprCycle.end(), vector<index_t>{x,y});
			if (it == reprCycle.begin()) { cout << "2  "; }
			else if (it == reprCycle.end()-1) { cout << "-1 "; }
			else if (it != reprCycle.end()) { cout << "1  "; }
			else { cout << "0  "; }   
		}
		cout << endl;
	}
}


value_t** CubicalGridComplex::allocateMemory() const {
	value_t** g = new value_t*[shape[0]+2];
    for (index_t i = 0; i < shape[0]+2; ++i) { g[i] = new value_t[shape[1]+2]; }
	if (g == NULL) { cerr << "out of memory!" << endl; }
	return g;
}


void CubicalGridComplex::getGridFromVector(const vector<value_t>& vec) {
	size_t counter = 0;
	grid = allocateMemory();
	for (index_t x = 0; x < shape[0]+2; ++x) {
		for (index_t y = 0; y < shape[1]+2; ++y) {
			if (x == 0 || x == shape[0]+1 || y == 0 || y == shape[1]+1) { grid[x][y] = INFTY; }
			else { grid[x][y] = vec[counter++]; }
		}
	}
}



UnionFind::UnionFind(const CubicalGridComplex& _cgc) : cgc(_cgc) {
	size_t n = cgc.getNumberOfCubes(0);
	parent.reserve(n);
	birthtime.reserve(n);
	index_t counter = 0;
	for (index_t x = 0; x < _cgc.shape[0]; ++x) {
		for (index_t y = 0; y < _cgc.shape[1]; ++y) {
			parent.push_back(counter++);
			birthtime.push_back(cgc.getBirth(x, y));
		}
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


index_t UnionFind::link(index_t x, index_t y) {
	if (birthtime[x] > birthtime[y]) {
		parent[x] = y; 
		return x;
	} else if (birthtime[x] < birthtime[y]) {
		parent[y] = x;
		return y;
	} else {
		if (x > y) {
			parent[x] = y;
			return x;
		} else {
			parent[y] = x;
			return y;
		}
	}
}


value_t UnionFind::getBirth(const index_t& idx) const { return birthtime[idx]; }


vector<index_t> UnionFind::getCoordinates(index_t idx) const { return {idx/cgc.shape[1],idx % cgc.shape[1]}; }


vector<index_t> UnionFind::getBoundaryIndices(const Cube& edge) const {
	vector<index_t> boundaryIndices(2);
	switch (edge.type()) {
		case 0:
			boundaryIndices[0] = edge.x()*cgc.shape[1] + edge.y();
			boundaryIndices[1] = (edge.x()+1)*cgc.shape[1] + edge.y();
			return boundaryIndices;

		case 1:
			boundaryIndices[0] = edge.x()*cgc.shape[1] + edge.y();
			boundaryIndices[1] = edge.x()*cgc.shape[1] + edge.y()+1;
			return boundaryIndices;
	}
	return boundaryIndices;
}


void UnionFind::reset() { for (size_t i = 0; i < parent.size(); ++i) { parent[i] = i; } }



UnionFindDual::UnionFindDual(const CubicalGridComplex& _cgc) : cgc(_cgc) {
	index_t n = cgc.getNumberOfCubes(2) + 1;
	parent.reserve(n);
	birthtime.reserve(n);
	index_t counter = 0;
	for (index_t x = 0; x < _cgc.m_x; ++x) {
		for (index_t y = 0; y < _cgc.m_y; ++y) {
			parent.push_back(counter++);
			birthtime.push_back(cgc.getBirth(x, y, 0, 2));
		}
	}
	parent.push_back(counter);
	birthtime.push_back(INFTY);
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


index_t UnionFindDual::link(index_t x, index_t y) {
	if (birthtime[x] < birthtime[y]) {
		parent[x] = y; 
		return x;
	} else if (birthtime[x] > birthtime[y]) {
		parent[y] = x;
		return y;
	} else {
		if (x < y) {
			parent[x] = y;
			return x;
		} else {
			parent[y] = x;
			return y;
		}
	}
}


value_t UnionFindDual::getBirth(const index_t& idx) const { return birthtime[idx]; }


vector<index_t> UnionFindDual::getCoordinates(index_t idx) const { return {idx/cgc.m_y,idx % cgc.m_y}; }


vector<index_t> UnionFindDual::getBoundaryIndices(const Cube& edge) const {
	vector<index_t> boundaryIndices(2);
	switch (edge.type()) {
		case 0:
			if (edge.y() == 0) { boundaryIndices[0] = cgc.m_xy; }
			else { boundaryIndices[0] = edge.x()*cgc.m_y + edge.y()-1; }
			if (edge.y() == cgc.m_y) { boundaryIndices[1] = cgc.m_xy; }
			else { boundaryIndices[1] = edge.x()*cgc.m_y + edge.y(); }
			return boundaryIndices;

		case 1:
			if (edge.x() == 0) { boundaryIndices[0] = cgc.m_xy; }
			else { boundaryIndices[0] = (edge.x()-1)*cgc.m_y + edge.y(); }
			if (edge.x() == cgc.m_x) { boundaryIndices[1] = cgc.m_xy; }
			else { boundaryIndices[1] = edge.x()*cgc.m_y + edge.y(); }	
			return boundaryIndices;
	}
	return boundaryIndices;
}


void UnionFindDual::reset() { for (size_t i = 0; i < parent.size(); ++i) { parent[i] = i; } }