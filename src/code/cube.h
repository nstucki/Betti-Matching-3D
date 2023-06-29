#pragma once

#include <cstdint>

#define NONE 0xffffffffffffffff

class Cube
{
public:
	double birth;
    uint64_t index;
	uint32_t x;
    uint32_t y;
    uint32_t z;
    uint8_t type;
    uint8_t dim;

	Cube();
    Cube(const Cube& v);
    Cube(double _b, uint32_t _x, uint32_t _y, uint32_t _z, uint8_t _t, uint8_t _dim);
    void copyCube(const Cube& c);
    void print();
    bool operator==(const Cube& rhs) const;
};

struct CubeComparator{
    bool operator()(const Cube& Cube1, const Cube& Cube2) const;
};