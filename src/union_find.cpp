#include "union_find.h"

UnionFind::UnionFind(CubicalGridComplex* const _cgc){
	uint64_t n = _cgc->shape[0] * _cgc->shape[1] * _cgc->shape[2];
	parent.resize(n);
	birthtime.resize(n);
	time_max.resize(n);

	uint32_t i=0;
	for (uint32_t x = 0; x < _cgc->shape[0]; x++) {
		for (uint32_t y = 0; y < _cgc->shape[1]; y++) {
			for(uint32_t z = 0; z < _cgc->shape[2]; z++){
				parent[i] = i;
				birthtime[i] = _cgc->getBirth(x,y,z);
				time_max[i] = birthtime[i];
				i++;
			}
		}
	}
}

uint32_t UnionFind::find(uint32_t x){
	uint32_t y = x, z = parent[y];
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

uint32_t UnionFind::link(uint32_t x, uint32_t y){
	if (x == y) return -1;
	if (birthtime[x] > birthtime[y]){
		parent[x] = y; 
		time_max[y] = max(time_max[x], time_max[y]);
		return x;
	} else if (birthtime[x] < birthtime[y]){
		parent[y] = x;
		time_max[x] = max(time_max[x], time_max[y]);
		return y;
	} else {
		if (x > y){
			parent[x] = y;
			time_max[y] = max(time_max[x], time_max[y]);
			return x;
		} else {
			parent[y] = x;
			time_max[x] = max(time_max[x], time_max[y]);
			return y;
		}
	}
}

UnionFindDual::UnionFindDual(CubicalGridComplex* const _cgc){
	uint64_t n = (_cgc->n_x) * (_cgc->n_y) * (_cgc->n_z) + 1;
	parent.resize(n);
	birthtime.resize(n);
	time_max.resize(n);

	uint32_t i=0;
	for (uint32_t x = 0; x < _cgc->n_x; x++) {
		for (uint32_t y = 0; y < _cgc->n_y; y++) {
			for(uint32_t z = 0; z < _cgc->n_z; z++){
				parent[i] = i;
				birthtime[i] = _cgc->getBirth(x,y,z,0,3);
				time_max[i] = birthtime[i];
				i++;
			}
		}
	}
	parent[i] = i;
	birthtime[i] = numeric_limits<double>::infinity();;
	time_max[i] = birthtime[i];
}

uint32_t UnionFindDual::find(uint32_t x){
	uint32_t y = x, z = parent[y];
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

uint32_t UnionFindDual::link(uint32_t x, uint32_t y){
	if (x == y) return -1;
	if (birthtime[x] < birthtime[y]){
		parent[x] = y; 
		time_max[y] = max(time_max[x], time_max[y]);
		return x;
	} else if (birthtime[x] > birthtime[y]){
		parent[y] = x;
		time_max[x] = max(time_max[x], time_max[y]);
		return y;
	} else {
		if (x < y){
			parent[x] = y;
			time_max[y] = max(time_max[x], time_max[y]);
			return x;
		} else {
			parent[y] = x;
			time_max[x] = max(time_max[x], time_max[y]);
			return y;
		}
	}
}