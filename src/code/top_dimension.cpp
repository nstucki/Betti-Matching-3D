#include "top_dimension.h"
#include "union_find.h"

TopDimension::TopDimension(CubicalGridComplex* _cgc, vector<WritePair>& _pairs){
	cgc = _cgc;
	pairs = &_pairs;
}

void TopDimension::enum_edges(vector<Cube>& ctr){
	ctr.clear();
	for (uint32_t x = 0; x < cgc->shape[0]; x++){
		for (uint32_t y = 0; y < cgc->shape[1]; y++){
			for(uint32_t z = 0; z < cgc->shape[2]; z++){
				for (uint8_t type = 0; type < 3; type++){
					double birth = cgc->getBirth(x,y,z,type,2);
					// add threshold
					if (birth != numeric_limits<double>::infinity()){
						ctr.push_back(Cube(birth,x,y,z,type,2));
					}
				}
			}
		}
	}
	sort(ctr.begin(), ctr.end(), CubeComparator());
}

void TopDimension::compute_pairs(vector<Cube>& ctr){
	UnionFindDual uf = UnionFindDual(cgc);
	enum_edges(ctr);
	uint32_t uind;
	uint32_t vind;
	uint32_t u;
	uint32_t v;
	uint32_t death_ind;
	double death;
	for (auto c = ctr.rbegin(), last = ctr.rend(); c != last; ++c){
		switch (c->type){
			case 0:
				if (c->x == 0){
					uind = (cgc->n_xyz);
				}
				else{
					uind = (c->x-1)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}
				if (c->x == cgc->n_x){
					vind = (cgc->n_xyz);
				}
				else{
					vind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}
				break;
			case 1:
				if (c->y == 0){
					uind = (cgc->n_xyz);
				}
				else{
					uind = (c->x)*(cgc->n_yz) + (c->y-1)*(cgc->n_z) + (c->z);
				}
				if (c->y == cgc->n_y){
					vind = (cgc->n_xyz);
				}
				else{
					vind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}	
				break;
			case 2:
				if (c->z == 0){
					uind = (cgc->n_xyz);
				}
				else{
					uind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z-1);
				}
				if (c->z == cgc->n_z){
					vind = (cgc->n_xyz);
				}
				else{
					vind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}
				break;
		}
		u = uf.find(uind);
		v = uf.find(vind);
		if (u!=v){
			death_ind = uf.link(u,v);
			death = uf.birthtime[death_ind];
			if (c->birth != death && death_ind != cgc->n_xyz){
				pairs->push_back(WritePair(*c, Cube(death, death_ind/(cgc->n_yz), death_ind/(cgc->n_z) % (cgc->n_y), death_ind % (cgc->n_z), 0, 3)));
			}
			c->index = NONE;
		}
	}
	// new list instead of deletions
	auto new_end = std::remove_if(ctr.begin(), ctr.end(),
								[](const Cube& c){ return c.index == NONE; });
	ctr.erase(new_end, ctr.end());
}


TopDimensionImage::TopDimensionImage(CubicalGridComplex* _cgc, CubicalGridComplex* _cgc_comp, vector<WritePair>& _pairs, vector<WritePair>& _pairs_im){
	cgc = _cgc;
	cgc_comp = _cgc_comp;
	pairs = &_pairs;
	pairs_im = &_pairs_im;
}

void TopDimensionImage::enum_edges(vector<Cube>& ctr){
	ctr.clear();
	for (uint32_t x = 0; x < cgc->shape[0]; x++){
		for (uint32_t y = 0; y < cgc->shape[1]; y++){
			for(uint32_t z = 0; z < cgc->shape[2]; z++){
				for (uint8_t type = 0; type < 3; type++){
					double birth = cgc->getBirth(x,y,z,type,2);
					if (birth != numeric_limits<double>::infinity()){
						ctr.push_back(Cube(birth,x,y,z,type,2));
					}
				}
			}
		}
	}
	sort(ctr.begin(), ctr.end(), CubeComparator());
}

void TopDimensionImage::compute_pairs(vector<Cube>& ctr){
	UnionFindDual uf = UnionFindDual(cgc);
	UnionFindDual uf_im = UnionFindDual(cgc_comp);
	enum_edges(ctr);
	uint32_t uind;
	uint32_t vind;
	uint32_t u;
	uint32_t v;
	uint32_t birth_ind;
	uint32_t death_ind;
	double birth;
	double death;
	for (auto c = ctr.rbegin(), last = ctr.rend(); c != last; ++c){
		switch (c->type){
			case 0:
				if (c->x == 0){
					uind = (cgc->n_xyz);
				}
				else{
					uind = (c->x-1)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}
				if (c->x == cgc->n_x){
					vind = (cgc->n_xyz);
				}
				else{
					vind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}
				break;
			case 1:
				if (c->y == 0){
					uind = (cgc->n_xyz);
				}
				else{
					uind = (c->x)*(cgc->n_yz) + (c->y-1)*(cgc->n_z) + (c->z);
				}
				if (c->y == cgc->n_y){
					vind = (cgc->n_xyz);
				}
				else{
					vind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}	
				break;
			case 2:
				if (c->z == 0){
					uind = (cgc->n_xyz);
				}
				else{
					uind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z-1);
				}
				if (c->z == cgc->n_z){
					vind = (cgc->n_xyz);
				}
				else{
					vind = (c->x)*(cgc->n_yz) + (c->y)*(cgc->n_z) + (c->z);
				}
				break;
		}
		u = uf.find(uind);
		v = uf.find(vind);
		if (u!=v){
			death_ind = uf.link(u,v);
			death = uf.birthtime[death_ind];
			if (c->birth != death && death_ind != cgc->n_xyz){
				pairs->push_back(WritePair(*c, Cube(death, death_ind/(cgc->n_yz), death_ind/(cgc->n_z) % (cgc->n_y), death_ind % (cgc->n_z), 0, 3)));
			}
			u = uf_im.find(uind);
			v = uf_im.find(vind);
			death_ind = uf_im.link(u,v);
			death = uf_im.birthtime[death_ind];
			if (death_ind != cgc_comp->n_xyz){
				pairs_im->push_back(WritePair(*c, Cube(death, death_ind/(cgc_comp->n_yz), death_ind/(cgc_comp->n_z) % (cgc_comp->n_y), death_ind % (cgc_comp->n_z), 0, 3)));
			}
			c->index = NONE;
		}
	}
	auto new_end = std::remove_if(ctr.begin(), ctr.end(),
								[](const Cube& c){ return c.index == NONE; });
	ctr.erase(new_end, ctr.end());
}
