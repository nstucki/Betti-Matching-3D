#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>

using namespace std;

namespace dimN {
template <typename Container>
Container multiplyFromRightExclusively(Container c) {
    std::reverse(c.begin(), c.end());
    std::exclusive_scan(c.begin(), c.end(), c.begin(), 1, std::multiplies<>{});
    std::reverse(c.begin(), c.end());
    return c;
}

template <typename Container, typename Scalar>
Container multiplyContainerElementwise(const Container &c, Scalar k) {
    Container r(c.size());
    transform(c.begin(), c.end(), r.begin(), [k](auto &c) { return c * k; });
    return r;
}

template <typename Container, typename Scalar>
Container addToContainerElementwise(const Container &c, Scalar k) {
    Container r(c.size());
    transform(c.begin(), c.end(), r.begin(), [k](auto &c) { return c + k; });
    return r;
}

template <typename Container, typename Scalar>
Container divideContainerElementwise(const Container &c, Scalar k) {
    Container r(c.size());
    transform(c.begin(), c.end(), r.begin(), [k](auto &c) { return c / k; });
    return r;
}

template <typename Container> void printContainer(const Container &c) {
    for (auto e : c) {
        cout << e << ' ';
    }
    cout << endl;
}
} // namespace dimN
