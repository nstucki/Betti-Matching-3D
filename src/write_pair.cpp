#include "write_pair.h"
#include <iostream>

using namespace std;

WritePair::WritePair(Cube _birth, Cube _death) {
    birth = _birth;
    death = _death;
}

void WritePair::print(){
    cout << "("; birth.print(); cout << ","; death.print(); cout << ")";
}