#pragma once

#include "data_structures.h"
#include "enumerators.h"

#include <cstdint>
#include <optional>
#include <queue>
#include <unordered_map>


namespace dim3 {
	typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQueue;

	template<typename _Tp>
    class Cube1Map {
		public:
		Cube1Map(vector<index_t> shape);
		void emplace(uint64_t cube_index, _Tp element);
		std::optional<_Tp> find(uint64_t cube_index) const;
		void clear();
		optional<_Tp>& operator[](int index);
		
		private:
		vector<std::optional<_Tp>> elements;
		uint64_t computeCoordinateIndex(uint64_t cube_index) const;
		vector<index_t> shape;
		const int size_x_direction;
		const int size_y_direction;
		const int size_z_direction;
	};
	

	class Dimension1 {
		public:
		Dimension1(const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
					const Config& config, vector<Pair>& pairs0, vector<Pair>& pairs1, vector<Pair>& pairsComp,
					vector<Match>& matches,
					Cube1Map<bool>& _isMatched0, Cube1Map<bool>& _isMatched1
					);
		void computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp);


		private:
		const CubicalGridComplex& cgc0;
		const CubicalGridComplex& cgc1;
		const CubicalGridComplex& cgcComp;
		const Config& config;
		vector<Pair>& pairs0;
		vector<Pair>& pairs1;
		vector<Pair>& pairsComp;
		vector<Match>& matches;
		Cube1Map<bool> &isMatched0;
		Cube1Map<bool> &isMatched1;
		Cube1Map<bool> isMatchedComp;

		Cube1Map<Pair> matchMap0;
		Cube1Map<Pair> matchMap1;
		Cube1Map<uint64_t> matchMapIm0;
		Cube1Map<uint64_t> matchMapIm1;
		Cube1Map<uint64_t> pivotColumnIndex;
		Cube1Map<CubeQueue> cache;

		void computePairs(const vector<Cube>& ctr, uint8_t k);
		void computePairsComp(vector<Cube>& ctr);
		void computeImagePairs(vector<Cube>& ctr, uint8_t k);
		void computeMatching();
		void enumerateEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const;
		Cube popPivot(CubeQueue& column) const;
		Cube getPivot(CubeQueue& column) const;
		void addCache(const index_t& i, CubeQueue& working_boundary);
	};
}
