#pragma once

#include "cubical_grid_complex.h"
#include "write_pair.h"

class TopDimension{
	private:
	vector<WritePair>* pairs;

	public:
	CubicalGridComplex* cgc;

	TopDimension(CubicalGridComplex* _cgc, vector<WritePair>& _pairs);
	void enum_edges(vector<Cube>& ctr);
	void compute_pairs(vector<Cube>& ctr);
};

class TopDimensionImage{
private:
	vector<WritePair>* pairs;
	vector<WritePair>* pairs_im;

public:
	CubicalGridComplex* cgc;
	CubicalGridComplex* cgc_comp;

	TopDimensionImage(CubicalGridComplex* _cgc, CubicalGridComplex* _cgc_comp, vector<WritePair>& _pairs, vector<WritePair>& _pairs_im);
	void enum_edges(vector<Cube>& ctr);
	void compute_pairs(vector<Cube>& ctr);
};