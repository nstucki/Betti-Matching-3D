#include "data_structures.h"
#include <iostream>

using namespace std;


Cube::Cube() {
    birth = 0;
    index = NONE;
}

Cube::Cube(const Cube& cube) : birth(cube.birth), index(cube.index) {}

Cube::Cube(value_t _birth, index_t _x, index_t _y, index_t _z, uint8_t _type) : birth(_birth) {
    index = ((uint64_t)_x << 44) | ((uint64_t)_y<<24) | ((uint64_t)_z<<4) | (uint64_t)_type;
}

index_t Cube::x() const { return ((index >> 44) & 0xfffff); }

index_t Cube::y() const { return ((index >> 24) & 0xfffff); }

index_t Cube::z() const { return ((index >> 4) & 0xfffff); }

uint8_t Cube::type() const { return( index & 0xf); }

bool Cube::operator==(const Cube& rhs) const{ return (index == rhs.index); }

void Cube::print() const {
    cout << "(" << birth << "," << x() << "," << y() << "," << z() << "," << unsigned(type()) << ")";
}


bool CubeComparator::operator()(const Cube& cube1, const Cube& cube2) const{
	if (cube1.birth == cube2.birth) { return (cube1.index < cube2.index); } 
	else { return (cube1.birth < cube2.birth); }
}


Pair::Pair() {}

Pair::Pair(const Cube& _birth, const Cube& _death) : birth(_birth), death(_death) {}

Pair::Pair(const Pair& pair) : birth(pair.birth), death(pair.death) {}

bool Pair::operator==(const Pair &rhs) const { return (birth == rhs.birth && death == rhs.death); }

void Pair::print() const { cout << "("; birth.print(); cout << ","; death.print(); cout << ")"; }


Match::Match(Pair* _pair0, Pair* _pair1) : pair0(_pair0), pair1(_pair1) {}

Match::Match(Pair& _pair_0, Pair& _pair_1) : pair0(&_pair_0), pair1(&_pair_1) {}

void Match::print() const { cout << "matched "; pair0->print(); cout << " with "; pair1->print(); cout << endl; }


CubicalGridComplex::CubicalGridComplex(const vector<value_t>& image, const vector<index_t>& _shape) : 
	shape(_shape), n_x(shape[0]-1), n_y(shape[1]-1), n_z(shape[2]-1), n_yz(n_y*n_z), n_xyz(n_x*n_yz), 
	m_yz(shape[1]*shape[2]), m_xyz(shape[0]*m_yz) {
	getGridFromVector(image);
}

CubicalGridComplex::~CubicalGridComplex(){
	for (int i = 0; i < shape[0]+2; i++){
        for (int j = 0; j < shape[1]+2; j++){
            delete[] grid[i][j];
        }
        delete[] grid[i];
    }
    delete[] grid;
}

index_t CubicalGridComplex::getNumberOfCubes(uint8_t dim) const {
	switch (dim) {
		case 0:
		return shape[0]*shape[1]*shape[2];

		case 1:
		return n_x*shape[1]*shape[2] + shape[0]*n_y*shape[2] + shape[0]*shape[1]*n_z;

		case 2:
		return n_x*n_y*shape[2] + n_x*shape[1]*n_z + shape[0]*n_y*n_z;

		case 3:
		return n_x*n_y*n_z;
	}
	return -1;
}

value_t CubicalGridComplex::getBirth(index_t x, index_t y, index_t z) const {
	return grid[x+1][y+1][z+1];
}

value_t CubicalGridComplex::getBirth(index_t x, index_t y, index_t z, uint8_t type, uint8_t dim) const {
	switch (dim) {
		case 0:
			return grid[x+1][y+1][z+1];
		case 1:
			switch (type) {
			case 0:
				return max(grid[x+1][y+1][z+1],grid[x+2][y+1][z+1]);
			case 1:
				return max(grid[x+1][y+1][z+1],grid[x+1][y+2][z+1]);
			case 2:
				return max(grid[x+1][y+1][z+1],grid[x+1][y+1][z+2]);
			}
		case 2:
			switch (type) {
			case 0:
				return max({grid[x+1][y+1][z+1],grid[x+1][y+2][z+1],grid[x+1][y+1][z+2],grid[x+1][y+2][z+2]});
			case 1:
				return max({grid[x+1][y+1][z+1],grid[x+2][y+1][z+1],grid[x+1][y+1][z+2],grid[x+2][y+1][z+2]});
			case 2:
				return max({grid[x+1][y+1][z+1],grid[x+2][y+1][z+1],grid[x+1][y+2][z+1],grid[x+2][y+2][z+1]});
			}
		case 3:
			return max({grid[x+1][y+1][z+1], grid[x+2][y+1][z+1], grid[x+1][y+2][z+1], grid[x+1][y+1][z+2],
						grid[x+2][y+2][z+1], grid[x+2][y+1][z+2], grid[x+1][y+2][z+2], grid[x+2][y+2][z+2] });
	}
	return numeric_limits<value_t>::infinity();
}

vector<index_t> CubicalGridComplex::getParentVoxel(const Cube& cube, uint8_t dim) const {
	index_t x = cube.x();
	index_t y = cube.y();
	index_t z = cube.z();
	switch (dim) {
		case 0:
			return {x,y,z};

		case 1:
		switch (cube.type()) {
			case 0:
			if (cube.birth == grid[x+2][y+1][z+1]) { return {x+1,y,z}; }
			else { return {x,y,z}; }

			case 1:
			if (cube.birth == grid[x+1][y+2][z+1]) { return {x,y+1,z}; }
			else { return {x,y,z}; }

			case 2:
			if (cube.birth == grid[x+1][y+1][z+2]) { return {x,y,z+1}; }
			else { return {x,y,z}; }
		}

		case 2:
		switch(cube.type()) {
			case 0:
			if (cube.birth == grid[x+1][y+2][z+2]) { return {x,y+1,z+1}; }
			else if (cube.birth == grid[x+1][y+2][z+1]) { return {x,y+1,z}; }
			else if(cube.birth == grid[x+1][y+1][z+2]) { return {x,y,z+1}; }
			else { return {x,y,z}; }

			case 1:
			if (cube.birth == grid[x+2][y+1][z+2]) { return {x+1,y,z+1}; }
			else if(cube.birth == grid[x+2][y+1][z+1]) { return {x+1,y,z}; }
			else if(cube.birth == grid[x+1][y+1][z+2]) { return {x,y,z+1}; }
			else { return {x,y,z}; }

			case 2:
			if (cube.birth == grid[x+2][y+2][z+1]) { return {x+1,y+1,z}; } 
			else if (cube.birth == grid[x+2][y+1][z+1]) { return {x+1,y,z}; }
			else if (cube.birth == grid[x+1][y+2][z+1]) { return {x,y+1,z}; } 
			else { return {x,y,z}; }
		}

		case 3:
		if (cube.birth == grid[x+2][y+2][z+2]){ return {x+1,y+1,z+1}; }
		else if (cube.birth == grid[x+2][y+2][z+1]) { return {x+1,y+1,z}; }
		else if (cube.birth == grid[x+2][y+1][z+2]) { return {x+1,y,z+1}; }
		else if (cube.birth == grid[x+2][y+1][z+1]) { return {x+1,y,z}; }
		else if (cube.birth == grid[x+1][y+2][z+2]) { return {x,y+1,z+1}; }
		else if (cube.birth == grid[x+1][y+2][z+1]) { return {x,y+1,z}; }
		else if(cube.birth == grid[x+1][y+1][z+2]) { return {x,y,z+1}; } 
		else { return {x,y,z}; }	
	}
	cerr << "parent voxel not found!" << endl;
	return {0,0,0};
}

void CubicalGridComplex::printImage() const {
    value_t birth;
    for (index_t x = 0; x < shape[0]; x++) {
        for (index_t y = 0; y < shape[1]; y++) {
            for (index_t z = 0; z < shape[2]; z++) {
                birth = getBirth(x, y, z);
                if (birth < 10) { cout << ' ' << birth << ' '; }
				else { cout << birth << ' '; }            
            }
            cout << endl;
        }
        cout << endl;
    }
}

value_t*** CubicalGridComplex::allocateMemory() const {
	value_t*** g = new value_t**[shape[0]+2];
    for (index_t i = 0; i < shape[0]+2; i++) {
        g[i] = new value_t*[shape[1]+2];
        for (index_t j = 0; j < shape[1]+2; j++) {
            g[i][j] = new value_t[shape[2]+2];
        }
    }
	if (g == NULL) {
		cerr << "Out of memory!" << endl;
	}
	return g;
}

void CubicalGridComplex::getGridFromVector(const vector<value_t> vec) {
	index_t i = 0;
	grid = allocateMemory();
	for (index_t x = 0; x < shape[0]+2; x++) {
		for (index_t y = 0; y < shape[1]+2; y++) {
			for (index_t z = 0; z < shape[2]+2; z++) {
				if (x == 0 || x == shape[0]+1 || y == 0 || y == shape[1]+1 || z == 0 || z == shape[2]+1){
					grid[x][y][z] = numeric_limits<value_t>::infinity();
				}
				else{
					grid[x][y][z] = vec[i++];
				}
			}
		}
	}
}


UnionFind::UnionFind(const CubicalGridComplex* const _cgc) : cgc(_cgc) {
	index_t n = _cgc->getNumberOfCubes(0);
	parent.reserve(n);
	birthtime.reserve(n);
	index_t i = 0;
	for (index_t x = 0; x < _cgc->shape[0]; x++) {
		for (index_t y = 0; y < _cgc->shape[1]; y++) {
			for(index_t z = 0; z < _cgc->shape[2]; z++){
				parent.push_back(i++);
				birthtime.push_back(_cgc->getBirth(x, y, z));
			}
		}
	}
}

index_t UnionFind::find(index_t x){
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

value_t UnionFind::getBirth(index_t x) const { return birthtime[x]; }

vector<index_t> UnionFind::getCoordinates(index_t x) const { 
	return {x/(cgc->m_yz), x/(cgc->shape[2]) % (cgc->shape[1]), x % (cgc->shape[2])};
}

vector<index_t> UnionFind::getBoundaryIndices(const Cube& edge) const {
	vector<index_t> boundaryIndices(2);
	switch (edge.type()) {
		case 0:
			boundaryIndices[0] = (edge.x())*(cgc->m_yz) + (edge.y())*(cgc->shape[2]) + (edge.z());
			boundaryIndices[1] = (edge.x()+1)*(cgc->m_yz) + (edge.y())*(cgc->shape[2]) + (edge.z());
			return boundaryIndices;

		case 1:
			boundaryIndices[0] = (edge.x())*(cgc->m_yz) + (edge.y())*(cgc->shape[2]) + (edge.z());
			boundaryIndices[1] = (edge.x())*(cgc->m_yz) + (edge.y()+1)*(cgc->shape[2]) + (edge.z());
			return boundaryIndices;

		case 2:
			boundaryIndices[0] = (edge.x())*(cgc->m_yz) + (edge.y())*(cgc->shape[2]) + (edge.z());
			boundaryIndices[1] = (edge.x())*(cgc->m_yz) + (edge.y())*(cgc->shape[2]) + (edge.z()+1);
			return boundaryIndices;
	}
	return boundaryIndices;
}

void UnionFind::reset() { for (index_t i = 0; i < parent.size(); i++) { parent[i] = i; } }


UnionFindDual::UnionFindDual(const CubicalGridComplex* const _cgc) : cgc(_cgc) {
	index_t n = _cgc->getNumberOfCubes(3) + 1;
	parent.reserve(n);
	birthtime.reserve(n);
	index_t i = 0;
	for (index_t x = 0; x < _cgc->n_x; x++) {
		for (index_t y = 0; y < _cgc->n_y; y++) {
			for(index_t z = 0; z < _cgc->n_z; z++){
				parent.push_back(i++);
				birthtime.push_back(_cgc->getBirth(x, y, z, 0, 3));
			}
		}
	}
	parent.push_back(i);
	birthtime.push_back(numeric_limits<value_t>::infinity());
}

index_t UnionFindDual::find(index_t x){
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

index_t UnionFindDual::link(index_t x, index_t y){
	if (birthtime[x] < birthtime[y]){
		parent[x] = y; 
		return x;
	} else if (birthtime[x] > birthtime[y]) {
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

value_t UnionFindDual::getBirth(index_t x) const { return birthtime[x]; }

vector<index_t> UnionFindDual::getCoordinates(index_t x) const { 
	return {x/(cgc->n_yz), x/(cgc->n_z) % (cgc->n_y), x % (cgc->n_z)};
}

vector<index_t> UnionFindDual::getBoundaryIndices(const Cube& edge) const {
	vector<index_t> boundaryIndices(2);
	switch (edge.type()) {
		case 0:
			if (edge.x() == 0) { boundaryIndices[0] = (cgc->n_xyz); }
			else { boundaryIndices[0] = (edge.x()-1)*(cgc->n_yz) + (edge.y())*(cgc->n_z) + (edge.z()); }
			if (edge.x() == cgc->n_x) { boundaryIndices[1] = (cgc->n_xyz); }
			else { boundaryIndices[1] = (edge.x())*(cgc->n_yz) + (edge.y())*(cgc->n_z) + (edge.z()); }
			return boundaryIndices;

		case 1:
			if (edge.y() == 0) { boundaryIndices[0] = (cgc->n_xyz); }
			else { boundaryIndices[0] = (edge.x())*(cgc->n_yz) + (edge.y()-1)*(cgc->n_z) + (edge.z()); }
			if (edge.y() == cgc->n_y) { boundaryIndices[1] = (cgc->n_xyz); }
			else { boundaryIndices[1] = (edge.x())*(cgc->n_yz) + (edge.y())*(cgc->n_z) + (edge.z()); }	
			return boundaryIndices;

		case 2:
			if (edge.z() == 0) { boundaryIndices[0] = (cgc->n_xyz); }
			else { boundaryIndices[0] = (edge.x())*(cgc->n_yz) + (edge.y())*(cgc->n_z) + (edge.z()-1); }
			if (edge.z() == cgc->n_z) { boundaryIndices[1] = (cgc->n_xyz); }
			else { boundaryIndices[1] = (edge.x())*(cgc->n_yz) + (edge.y())*(cgc->n_z) + (edge.z()); }
			return boundaryIndices;
	}
	return boundaryIndices;
}

void UnionFindDual::reset() { for (index_t i = 0; i < parent.size(); i++) { parent[i] = i; } }