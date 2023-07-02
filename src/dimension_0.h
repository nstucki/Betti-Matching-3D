#pragma once

#include "cubical_grid_complex.h"
#include "write_pair.h"

class Dimension0{
	private:
	CubicalGridComplex* cgc;
	vector<WritePair>* pairs;

	public:
	Dimension0(CubicalGridComplex* _cgc, vector<WritePair>& _pairs);
	void enum_edges(vector<Cube>& ctr);
	void compute_pairs(vector<Cube>& ctr);
};

class Dimension0Image{
	private:
	CubicalGridComplex* cgc_0;
	CubicalGridComplex* cgc_1;
	CubicalGridComplex* cgc_comp;
	vector<WritePair>* pairs_im_0;
	vector<WritePair>* pairs_im_1;
	vector<WritePair>* pairs_comp;
	
	public:
	Dimension0Image(CubicalGridComplex* _cgc_0, CubicalGridComplex* _cgc_1, CubicalGridComplex* _cgc_comp, vector<WritePair>& pairs_im_0, vector<WritePair>& pairs_im_1, vector<WritePair>& pairs_comp);
	void enum_edges(vector<Cube>& ctr);
	void compute_pairs(vector<Cube>& ctr);
};
