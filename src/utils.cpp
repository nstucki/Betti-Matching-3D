#include "utils.h"
#include "npy.hpp"

#include <cassert>
#include <cfloat>

void readImage(const string &filename, const fileFormat &format,
               vector<double> &image, vector<index_t> &shape) {
    switch (format) {
    case DIPHA: {
        ifstream fin(filename, ios::in | ios::binary);
        int64_t d;
        fin.read((char *)&d, sizeof(int64_t));
        assert(d == 8067171840);
        fin.read((char *)&d, sizeof(int64_t));
        assert(d == 1);
        fin.read((char *)&d, sizeof(int64_t));
        fin.read((char *)&d, sizeof(int64_t));
        uint8_t dim = d;
        assert(dim < 4);
        uint64_t n;
        fin.read((char *)&d, sizeof(int64_t));
        shape.push_back(d);
        n = d;
        if (dim > 1) {
            fin.read((char *)&d, sizeof(int64_t));
            shape.push_back(d);
            n *= d;
        }
        if (dim > 2) {
            fin.read((char *)&d, sizeof(int64_t));
            shape.push_back(d);
            n *= d;
        }
        double value;
        image.reserve(n);
        while (!fin.eof()) {
            fin.read((char *)&value, sizeof(double));
            image.push_back(value);
        }
        fin.close();
        return;
    }

    case PERSEUS: {
        ifstream reading_file;
        reading_file.open(filename.c_str(), ios::in);
        string reading_line_buffer;
        getline(reading_file, reading_line_buffer);
        uint8_t dim = atoi(reading_line_buffer.c_str());
        uint64_t n;
        getline(reading_file, reading_line_buffer);
        shape.push_back(atoi(reading_line_buffer.c_str()));
        n = shape[0];
        if (dim > 1) {
            getline(reading_file, reading_line_buffer);
            shape.push_back(atoi(reading_line_buffer.c_str()));
            n *= shape[1];
        }
        if (dim > 2) {
            getline(reading_file, reading_line_buffer);
            shape.push_back(atoi(reading_line_buffer.c_str()));
            n *= shape[2];
        }
        image.reserve(n);
        double value;
        while (!reading_file.eof()) {
            getline(reading_file, reading_line_buffer);
            value = atof(reading_line_buffer.c_str());
            if (value != -1) {
                image.push_back(value);
            } else {
                image.push_back(DBL_MAX);
            }
        }
        reading_file.close();
        return;
    }

    case NUMPY: {
        vector<unsigned long> _shape;
        bool fortran_order;
        try {
            npy::LoadArrayFromNumpy(filename.c_str(), _shape, fortran_order,
                                    image);
        } catch (...) {
            cerr << "The data type of an numpy array should be numpy.float64."
                 << endl;
            exit(-2);
        }
        if (fortran_order) {
            cerr << "The array at " << filename
                 << " is not C-contiguous. Convert using "
                    "numpy.ascontiguousarray() before saving."
                 << endl;
            exit(-2);
        }
        for (uint32_t i : _shape) {
            shape.push_back(i);
        }
        return;
    }
    }
}

vector<vector<bool>> getSubsets(index_t n, index_t k) {
    if (k == 0) {
        vector<vector<bool>> subsets(1, vector<bool>(n, false));
        return subsets;
    } else {
        vector<vector<bool>> subsets;
        vector<vector<bool>> prevSubsets = getSubsets(n, k - 1);
        vector<bool> subset;
        for (auto &prevSubset : prevSubsets) {
            for (index_t i = 0; i < n; i++) {
                if (prevSubset[i]) {
                    break;
                } else {
                    subset = prevSubset;
                    subset[i] = true;
                    subsets.push_back(subset);
                }
            }
        }
        return subsets;
    }
}
