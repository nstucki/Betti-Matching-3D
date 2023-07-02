#include "dimension_0.h"
#include "union_find.h"

Dimension0::Dimension0(CubicalGridComplex* _cgc, vector<WritePair>& _pairs){
	cgc = _cgc;
	pairs = &_pairs;
}

void Dimension0::enum_edges(vector<Cube>& ctr){
	ctr.clear();
	for (uint32_t x = 0; x < cgc->shape[0]; x++){
		for (uint32_t y = 0; y < cgc->shape[1]; y++){
			for (uint32_t z = 0; z < cgc->shape[2]; z++){
				for (uint8_t type = 0; type < 3; type++){
					double birth = cgc->getBirth(x,y,z,type,1);
					if (birth != numeric_limits<double>::infinity()){
						ctr.push_back(Cube(birth,x,y,z,type,1));
					}	
				}				
			}
		}
	}
	sort(ctr.begin(), ctr.end(), CubeComparator());
}

void Dimension0::compute_pairs(vector<Cube>& ctr){
	UnionFind uf = UnionFind(cgc);
	enum_edges(ctr);
	uint32_t uind;
	uint32_t vind;
	uint32_t u;
	uint32_t v;
	uint32_t birth_ind;
	double birth;
	for (auto c = ctr.begin(), last = ctr.end(); c != last; ++c){
		switch (c->type){
			case 0:
				uind = (c->x)*(cgc->m_yz) + (c->y)*(cgc->shape[2]) + (c->z);
				vind = (c->x+1)*(cgc->m_yz) + (c->y)*(cgc->shape[2]) + (c->z);
				break;
			case 1:
				uind = (c->x)*(cgc->m_yz) + (c->y)*(cgc->shape[2]) + (c->z);
				vind = (c->x)*(cgc->m_yz) + (c->y+1)*(cgc->shape[2]) + (c->z);
				break;
			case 2:
				uind = (c->x)*(cgc->m_yz) + (c->y)*(cgc->shape[2]) + (c->z);
				vind = (c->x)*(cgc->m_yz) + (c->y)*(cgc->shape[2]) + (c->z+1);
				break;
		}
		u = uf.find(uind);
		v = uf.find(vind);
		if (u!=v){
			birth_ind = uf.link(u,v);
			birth = uf.birthtime[birth_ind];
			if (birth != c->birth){
				pairs->push_back(WritePair(Cube(birth, birth_ind/(cgc->m_yz), birth_ind/(cgc->shape[2]) % (cgc->shape[1]), birth_ind % (cgc->shape[2]), 0, 0), *c));
			}
		}
	}
}


Dimension0Image::Dimension0Image(CubicalGridComplex* _cgc_0, CubicalGridComplex* _cgc_1, CubicalGridComplex* _cgc_comp, vector<WritePair>& _pairs_im_0, vector<WritePair>& _pairs_im_1, vector<WritePair>& _pairs_comp){
	cgc_0 = _cgc_0;
	cgc_1 = _cgc_1;
	cgc_comp = _cgc_comp;
	pairs_im_0 = &_pairs_im_0;
	pairs_im_1 = &_pairs_im_1;
	pairs_comp = &_pairs_comp;
}

void Dimension0Image::enum_edges(vector<Cube>& ctr){
	ctr.clear();
	for (uint32_t x = 0; x < cgc_comp->shape[0]; x++){
		for (uint32_t y = 0; y < cgc_comp->shape[1]; y++){
			for (uint32_t z = 0; z < cgc_comp->shape[2]; z++){
				for (uint8_t type = 0; type < 3; type++){
					double birth = cgc_comp->getBirth(x,y,z,type,1);
					if (birth != numeric_limits<double>::infinity()){
						ctr.push_back(Cube(birth,x,y,z,type,1));
					}	
				}				
			}
		}
	}
	sort(ctr.begin(), ctr.end(), CubeComparator());
}

void Dimension0Image::compute_pairs(vector<Cube>& ctr){
	UnionFind uf_0 = UnionFind(cgc_0);
	UnionFind uf_1 = UnionFind(cgc_1);
	UnionFind uf_comp = UnionFind(cgc_comp);
	enum_edges(ctr);
	uint32_t uind;
	uint32_t vind;
	uint32_t u;
	uint32_t v;
	uint32_t birth_ind;
	double birth;
	for (auto c = ctr.begin(), last = ctr.end(); c != last; ++c){
		switch (c->type){
			case 0:
				uind = (c->x)*(cgc_comp->m_yz) + (c->y)*(cgc_comp->shape[2]) + (c->z);
				vind = (c->x+1)*(cgc_comp->m_yz) + (c->y)*(cgc_comp->shape[2]) + (c->z);
				break;
			case 1:
				uind = (c->x)*(cgc_comp->m_yz) + (c->y)*(cgc_comp->shape[2]) + (c->z);
				vind = (c->x)*(cgc_comp->m_yz) + (c->y+1)*(cgc_comp->shape[2]) + (c->z);
				break;
			case 2:
				uind = (c->x)*(cgc_comp->m_yz) + (c->y)*(cgc_comp->shape[2]) + (c->z);
				vind = (c->x)*(cgc_comp->m_yz) + (c->y)*(cgc_comp->shape[2]) + (c->z+1);
				break;
		}
		u = uf_comp.find(uind);
		v = uf_comp.find(vind);
		if (u!=v){
			birth_ind = uf_comp.link(u,v);
			birth = uf_comp.birthtime[birth_ind];
			if (birth != c->birth){
				pairs_comp->push_back(WritePair(Cube(birth, birth_ind/(cgc_comp->m_yz), birth_ind/(cgc_comp->shape[2]) % (cgc_comp->shape[1]), birth_ind % (cgc_comp->shape[2]), 0, 0), *c));
			}
			u = uf_0.find(uind);
			v = uf_0.find(vind);
			birth_ind = uf_0.link(u,v);
			birth = uf_0.birthtime[birth_ind];
			pairs_im_0->push_back(WritePair(Cube(birth, birth_ind/(cgc_comp->m_yz), birth_ind/(cgc_comp->shape[2]) % (cgc_comp->shape[1]), birth_ind % (cgc_comp->shape[2]), 0, 0), *c));
			u = uf_1.find(uind);
			v = uf_1.find(vind);
			birth_ind = uf_1.link(u,v);
			birth = uf_1.birthtime[birth_ind];
			pairs_im_1->push_back(WritePair(Cube(birth, birth_ind/(cgc_comp->m_yz), birth_ind/(cgc_comp->shape[2]) % (cgc_comp->shape[1]), birth_ind % (cgc_comp->shape[2]), 0, 0), *c));
		}
	}
}