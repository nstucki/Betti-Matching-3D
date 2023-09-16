#include "data_structures.h"

#include <iostream>

VoxelPair::VoxelPair(const vector<index_t> &_birth, const vector<index_t> &_death) : birth(_birth), death(_death) {}

void VoxelPair::print() const
{
    cout << "((";
    for (int i = 0; i < birth.size()-1; i++) {
        cout << birth[i] << " ";
    };
    cout << birth[birth.size()-1] << ");(";
    for (int i = 0; i < death.size()-1; i++) {
        cout << death[i] << " ";
    }
    cout << death[death.size()-1] << "))";
}

VoxelMatch::VoxelMatch(const VoxelPair &_pair0, const VoxelPair &_pair1) : pair0(_pair0), pair1(_pair1) {}

void VoxelMatch::print() const
{
    pair0.print();
    cout << " <-> ";
    pair1.print();
    cout << endl;
}
