#pragma once

#include "../config.h"

#include <cstdint>
#include <vector>

using namespace std;


namespace dim3 {
	class CubicalGridComplex;

	class Cube {
		public:
		value_t birth;
		uint64_t index;
#ifdef CUBE_DEBUG_INFO
		int dimension = -1;
		// const std::optional<std::reference_wrapper<CubicalGridComplex>> parentCgc;
		const CubicalGridComplex *parentCgc = nullptr;
#endif

		Cube();
		Cube(value_t birth, index_t x, index_t y, index_t z, uint8_t type
#ifdef CUBE_DEBUG_INFO
	, const CubicalGridComplex *parentCgc = nullptr
#endif
		);
		Cube(const Cube& cube);
		index_t x() const;
		index_t y() const;
		index_t z() const;
		uint8_t type() const;
		bool operator==(const Cube& rhs) const;
		void print() const;
	};

	class Cube0 : public Cube {
		public:
		Cube0(value_t birth, index_t x, index_t y, index_t z, uint8_t type, const CubicalGridComplex *parentCgc = nullptr) : Cube(birth, x, y, z, type, parentCgc) {
#ifdef CUBE_DEBUG_INFO
			dimension = 0;
#endif
		}
		Cube0() : Cube() {dimension = 0;}
	};

	class Cube1 : public Cube {
		public:
		Cube1(value_t birth, index_t x, index_t y, index_t z, uint8_t type, const CubicalGridComplex *parentCgc = nullptr) : Cube(birth, x, y, z, type, parentCgc) {
#ifdef CUBE_DEBUG_INFO
			dimension = 1;
#endif
		}
		Cube1() : Cube() {dimension = 1;}
	};

	class Cube2 : public Cube {
		public:
		Cube2(value_t birth, index_t x, index_t y, index_t z, uint8_t type, const CubicalGridComplex *parentCgc = nullptr) : Cube(birth, x, y, z, type, parentCgc) {
#ifdef CUBE_DEBUG_INFO
			dimension = 2;
#endif
		}
		Cube2() : Cube() {dimension = 2;}
	};

	struct CubeComparator{ bool operator()(const Cube& Cube1, const Cube& Cube2) const; };


	class Pair {
		public:
		const Cube birth;
		const Cube death;

		Pair();
		Pair(const Cube& birth, const Cube& death);
		Pair(const Pair& pair);
		bool operator==(const Pair &rhs) const;
		void print() const;
	};


	class Match {
		public:
		Pair pair0;
		Pair pair1;

		Match(Pair pair0, Pair pair1);
		void print() const;
	};

	class CubicalGridComplex {
		public:
		const vector<index_t> shape;
		const index_t m_x;
		const index_t m_y;
		const index_t m_z;
		const index_t m_yz;
		const index_t m_xyz;
		const index_t n_yz;
		const index_t n_xyz;
		CubicalGridComplex(const vector<value_t>& image, const vector<index_t>& shape);
		CubicalGridComplex(CubicalGridComplex &&other);
		~CubicalGridComplex();
		size_t getNumberOfCubes(const uint8_t& dim) const;
		value_t getBirth(const index_t& x, const index_t& y, const index_t& z) const;
		value_t getBirth(const index_t& x, const index_t& y, const index_t& z, const uint8_t& type, const uint8_t& dim) const;
		vector<index_t> getParentVoxel(const Cube& c, const uint8_t& dim) const;
		void printImage() const;

		private:
		value_t*** grid;
		value_t*** allocateMemory() const;
		void getGridFromVector(const vector<value_t>& vector);
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
