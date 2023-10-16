#pragma once

#include "../config.h"

#include <cstdint>
#include <vector>
#include <set>

using namespace std;


namespace dim3 {
		class Cube {
		public:
		Cube();
		Cube(value_t birth, index_t x, index_t y, index_t z, uint8_t type);
		Cube(const Cube& cube);
		bool operator==(const Cube& rhs) const;
		index_t x() const;
		index_t y() const;
		index_t z() const;
		uint8_t type() const;
		void print() const;
		value_t birth;
		uint64_t index;
	};



	struct CubeComparator{ bool operator()(const Cube& Cube1, const Cube& Cube2) const; };



	class Pair {
		public:
		Pair();
		Pair(const Cube& birth, const Cube& death);
		Pair(const Pair& pair);
		bool operator==(const Pair &rhs) const;
		void print() const;
		const Cube birth;
		const Cube death;
	};



	class Match {
		public:
		Match(Pair pair0, Pair pair1);
		void print() const;
		Pair pair0;
		Pair pair1;
	};



	class CubicalGridComplex {
		public:
		CubicalGridComplex(const vector<value_t>& image, const vector<index_t>& shape);
		CubicalGridComplex(CubicalGridComplex &&other);
		~CubicalGridComplex();
		size_t getNumberOfCubes(const uint8_t& dim) const;
		value_t getBirth(const index_t& x, const index_t& y, const index_t& z) const;
		value_t getBirth(const index_t& x, const index_t& y, const index_t& z, const uint8_t& type, const uint8_t& dim) const;
		vector<index_t> getParentVoxel(const Cube& c, const uint8_t& dim) const;
		void printImage() const;
		void printRepresentativeCycle(const set<vector<index_t>>& reprCycle) const;
		const vector<index_t> shape;
		const index_t m_x;
		const index_t m_y;
		const index_t m_z;
		const index_t m_yz;
		const index_t m_xyz;
		const index_t n_yz;
		const index_t n_xyz;

		private:
		value_t*** allocateMemory() const;
		void getGridFromVector(const vector<value_t>& vector);
		value_t*** grid;
	};


	class UnionFind {
		public:
		UnionFind(const CubicalGridComplex& cgc);
		index_t find(index_t x);
		index_t link(index_t x, index_t y);
		value_t getBirth(const index_t& idx) const;
		vector<index_t> getCoordinates(index_t idx) const;
		vector<index_t> getBoundaryIndices(const Cube& edge) const;
		void reset();

		private:
		vector<index_t> parent;
		vector<value_t> birthtime;
		const CubicalGridComplex& cgc;
	};

	class UnionFindDual {
		public:
		UnionFindDual(const CubicalGridComplex& cgc);
		index_t find(index_t x);
		index_t link(index_t x, index_t y);
		value_t getBirth(const index_t& idx) const;
		vector<index_t> getCoordinates(index_t idx) const;
		vector<index_t> getBoundaryIndices(const Cube& edge) const;
		void reset();

		private:
		vector<index_t> parent;
		vector<value_t> birthtime;
		const CubicalGridComplex& cgc;
	};
}
