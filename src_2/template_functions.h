#pragma once

#include <iostream>

using namespace std;

template <typename Container, typename Scalar>
Container DivideContainerElementwise(const Container &c, Scalar k){
	Container r = c;
	transform(c.begin(), c.end(), r.begin(), [k](auto &c){ return c/k; });
	return r;
}

template <typename Container>
void printContainer(const Container& c) {
    for (auto e : c) {
        cout << e << ' ';
    }
    cout << endl;
}