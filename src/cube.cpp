#include "cube.h"
#include <iostream>

using namespace std;

Cube::Cube() {
    birth = 0;
    index = NONE;
    x = 0;
    y = 0;
    z = 0;
    type = 0;
    dim = 0;
}

Cube::Cube(const Cube& c) {
    birth = c.birth;
    index = c.index;
    x = c.x;
    y = c.y;
    z = c.z;
    type = c.type;
    dim = c.dim;
}


Cube::Cube(double _b, uint32_t _x, uint32_t _y, uint32_t _z, uint8_t _type, uint8_t _dim) {
    birth = _b;
    x = _x;
    y = _y;
    z = _z;
    type = _type;
    dim = _dim;
    index = (uint64_t)_x | ((uint64_t)_y<<20) | ((uint64_t)_z<<40) | ((uint64_t)_type<<60);
}

    
void Cube::copyCube(const Cube& c) {
    birth = c.birth;
    index = c.index;
    x = c.x;
    y = c.y;
    z = c.z;
    type = c.type;
    dim = c.dim;
}

void Cube::print() {
    cout << "(" << birth << "," << x << "," << y << "," << z << "," << unsigned(type) << "," << unsigned(dim) << ")";
}


bool Cube::operator==(const Cube& rhs) const{
    return (x == rhs.x && y == rhs.y && z == rhs.z && type == rhs.type && dim == rhs.dim);
}


bool CubeComparator::operator()(const Cube& Cube1, const Cube& Cube2) const{
    if(Cube1.birth == Cube2.birth){
        if(Cube1.x == Cube2.x){
            if(Cube1.y == Cube2.y){
                if(Cube1.z == Cube2.z){
                    return Cube1.type < Cube2.type;
                }else{
                    return Cube1.z < Cube2.z;
                }
            }else{
                return Cube1.y < Cube2.y;
            }
        }else{
            return Cube1.x < Cube2.x;
        }
    }else{
        return Cube1.birth < Cube2.birth;
    }
}