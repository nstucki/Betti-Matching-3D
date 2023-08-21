#pragma once

#include "config.h"
#include <vector>

using namespace std;

class VoxelPair
{
public:
    const vector<index_t> birth;
    const vector<index_t> death;

    VoxelPair(const vector<index_t> &birth, const vector<index_t> &death);
    void print() const;
};

class VoxelMatch
{
public:
    const VoxelPair pair0;
    const VoxelPair pair1;

    VoxelMatch(const VoxelPair &pair0, const VoxelPair &pair1);
    void print() const;
};