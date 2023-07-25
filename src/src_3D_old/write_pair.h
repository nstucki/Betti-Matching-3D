#pragma once

#include "cube.h"

class WritePair{
public:
	Cube birth;
	Cube death;

    WritePair(Cube _birth, Cube _death);
	bool operator==(const WritePair &rhs) const;
	void print();
};
