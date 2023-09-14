#pragma once

#include "../config.h"

#include <cstdint>
#include <vector>
#include <optional>

using namespace std;


namespace dim3 {
		class Cube {
		public:
		value_t birth;
		uint64_t index;

		Cube();
		Cube(value_t birth, index_t x, index_t y, index_t z, uint8_t type);
		Cube(const Cube& cube);
		index_t x() const;
		index_t y() const;
		index_t z() const;
		uint8_t type() const;
		bool operator==(const Cube& rhs) const;
		void print() const;
	};


	struct CubeComparator{ bool operator()(const Cube& Cube1, const Cube& Cube2) const; };


	class Pair {
		public:
		Cube birth;
		Cube death;

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



	template<typename _Tp>
    class Cube1Map {
		public:
		Cube1Map(vector<index_t> shape);
		void emplace(uint64_t cube_index, _Tp element);
		const std::optional<_Tp>& find(uint64_t cube_index) const;
		void clear();
		optional<_Tp>& operator[](int index);

		private:
		vector<std::optional<_Tp>> elements;
		uint64_t computeCoordinateIndex(uint64_t cube_index) const;
		vector<index_t> shape;
		std::optional<_Tp> none = {};
		const int size_x_direction;
		const int size_y_direction;
		const int size_z_direction;
	};

	template<class _Tp>
	Cube1Map<_Tp>::Cube1Map(vector<index_t> shape) : shape(shape), size_x_direction(shape[1] * shape[2] * 3), size_y_direction(shape[2] * 3), size_z_direction(3), elements(shape[0] * shape[1] * shape[2] * 3) {}

	template<class _Tp>
	void Cube1Map<_Tp>::emplace(uint64_t cube_index, _Tp element) {
		if (cube_index != NONE) {
			elements[computeCoordinateIndex(cube_index)] = element;
		} else {
			throw runtime_error("Cube1Map::emplace may not be called with NONE magic number");
		}
	}

	template<class _Tp>
	const std::optional<_Tp>& Cube1Map<_Tp>::find(uint64_t cube_index) const {
		if (cube_index != NONE) {
			return elements[computeCoordinateIndex(cube_index)];
		}
		return none;
	}

	template<class _Tp>
	optional<_Tp>& Cube1Map<_Tp>::operator[](int cube_index) {
		if (cube_index != NONE) {
			return elements[computeCoordinateIndex(cube_index)];
		}
		throw runtime_error("Cube1Map subscript operator may not be called with NONE magic number");
	}

	template<class _Tp>
	uint64_t Cube1Map<_Tp>::computeCoordinateIndex(uint64_t cube_index) const {
		int x = (cube_index >> 44) & 0xfffff;
		int y = (cube_index >> 24) & 0xfffff;
		int z = (cube_index >> 4) & 0xfffff;
		int type = cube_index & 0xf;

		return x * size_x_direction + y * size_y_direction + z * size_z_direction + type;
	}

	template<class _Tp>
	void Cube1Map<_Tp>::clear() {
		elements.clear();
		elements.resize(shape[0] * shape[1] * shape[2] * 3);
	}

}
