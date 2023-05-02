#include "cubical_grid_complex.h"
#include <iostream>

CubicalGridComplex::CubicalGridComplex(const vector<double>& image, const vector<uint32_t> _shape){
	dim = 0;
	shape = _shape;
	for (uint8_t i : shape){
		if(shape[i] != 1){
			dim += 1;
		}
	}
	n_x = shape[0]-1;
	n_y = shape[1]-1;
	n_z = shape[2]-1;
	n_yz = n_y*n_z;
	n_xyz = n_x*n_yz;
	m_yz = shape[1]*shape[2];
	m_xyz = shape[0]*m_yz;
	gridFromVector(image);
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

double*** CubicalGridComplex::allocate_memory() {
	double*** a = new double**[shape[0]+2];
 
    for (int i = 0; i < shape[0]+2; i++) {
        a[i] = new double*[shape[1]+2];
        for (int j = 0; j < shape[1]+2; j++) {
            a[i][j] = new double[shape[2]+2];
        }
    }
	if (a == NULL) {
		cerr << "Out of memory!" << endl;
	}
	return a;
}

void CubicalGridComplex::gridFromVector(const vector<double> vec) {
	uint64_t i = 0;
	grid = allocate_memory();
	for (uint32_t x = 0; x < shape[0]+2; x++) {
		for (uint32_t y = 0; y < shape[1]+2; y++) {
			for (uint32_t z = 0; z < shape[2]+2; z++) {
				if (x == 0 || x == shape[0]+1 || y == 0 || y == shape[1]+1 || z == 0 || z == shape[2]+1){
					grid[x][y][z] = numeric_limits<double>::infinity();
				}
				else{
					grid[x][y][z] = vec[i++];
				}
			}
		}
	}
}

double CubicalGridComplex::getBirth(uint32_t x, uint32_t y, uint32_t z){
	return grid[x+1][y+1][z+1];
}

double CubicalGridComplex::getBirth(uint32_t x, uint32_t y, uint32_t z, uint8_t t, uint8_t dim){
	switch (dim) {
		case 0:
			return grid[x+1][y+1][z+1];
		case 1:
			switch (t) {
			case 0:
				return max(grid[x+1][y+1][z+1],grid[x+2][y+1][z+1]);
			case 1:
				return max(grid[x+1][y+1][z+1],grid[x+1][y+2][z+1]);
			case 2:
				return max(grid[x+1][y+1][z+1],grid[x+1][y+1][z+2]);
			}
		case 2:
			switch (t) {
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
	return numeric_limits<double>::infinity();
}

vector<uint32_t> CubicalGridComplex::ParentVoxel(Cube &c){
	uint32_t x = c.x;
	uint32_t y = c.y;
	uint32_t z = c.z;

	switch(c.dim){
		case 0:
			return {x,y,z};
		case 1:
		switch(c.type){
			case 0:
			if (c.birth == grid[x+2][y+1][z+1]){
				return {x+1,y,z};
			}else{
				return {x,y,z};
			}
			case 1:
			if (c.birth == grid[x+1][y+2][z+1]){
				return {x,y+1,z};
			}else{
				return {x,y,z};
			}
			case 2:
			if (c.birth == grid[x+1][y+1][z+2]){
				return {x,y,z+1};
			}else{
				return {x,y,z};
			}
		}
		case 2:
		switch(c.type){
			case 0:
			if (c.birth == grid[x+1][y+2][z+2]){
				return {x,y+1,z+1};
			}else if(c.birth == grid[x+1][y+2][z+1]){
				return {x,y+1,z};
			}else if(c.birth == grid[x+1][y+1][z+2]){
				return {x,y,z+1};
			}else{
				return {x,y,z};
			}
			case 1:
			if (c.birth == grid[x+2][y+1][z+2]){
				return {x+1,y,z+1};
			}else if(c.birth == grid[x+2][y+1][z+1]){
				return {x+1,y,z};
			}else if(c.birth == grid[x+1][y+1][z+2]){
				return {x,y,z+1};
			}else{
				return {x,y,z};
			}
			case 2:
			if (c.birth == grid[x+2][y+2][z+1]){
				return {x+1,y+1,z};
			}else if(c.birth == grid[x+2][y+1][z+1]){
				return {x+1,y,z};
			}else if(c.birth == grid[x+1][y+2][z+1]){
				return {x,y+1,z};
			}else{
				return {x,y,z};
			}
		}
		case 3:
		if (c.birth == grid[x+2][y+2][z+2]){
			return {x+1,y+1,z+1};
		}else if(c.birth == grid[x+2][y+2][z+1]){
			return {x+1,y+1,z};
		}else if(c.birth == grid[x+2][y+1][z+2]){
			return {x+1,y,z+1};
		}else if(c.birth == grid[x+2][y+1][z+1]){
			return {x+1,y,z};
		}else if(c.birth == grid[x+1][y+2][z+2]){
			return {x,y+1,z+1};
		}else if(c.birth == grid[x+1][y+2][z+1]){
			return {x,y+1,z};
		}else if(c.birth == grid[x+1][y+1][z+2]){
			return {x,y,z+1};
		}else{
			return {x,y,z};
		}	
	}
	cerr << "parent voxel not found!" << endl;
	return {0,0,0};
}