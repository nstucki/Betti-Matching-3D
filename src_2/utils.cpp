#include "utils.h"

#include "npy.hpp"

#include <vector>
#include <cfloat>

using namespace std;


vector<vector<bool>> getSubsets(index_t n, index_t k) {
    if (k == 0) {
        vector<vector<bool>> subsets(1, vector<bool>(n, false));
        return subsets;
    } else {
        vector<vector<bool>> subsets;
        vector<vector<bool>> prevSubsets = getSubsets(n, k-1);
        vector<bool> subset;
        for (auto& prevSubset : prevSubsets) {
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


void readImage(string const &filename, fileFormat const &format, vector<double> &image, vector<index_t> &shape) {
    switch(format){
        case DIPHA: {
            ifstream fin(filename, ios::in | ios::binary );
            int64_t d;
            fin.read((char *) &d, sizeof(int64_t));
            assert(d == 8067171840);
            fin.read((char *) &d, sizeof(int64_t));
            assert(d == 1);
            fin.read((char *) &d, sizeof(int64_t));
            fin.read((char *) &d, sizeof(int64_t));
            uint8_t dim = d;
            index_t n;
            fin.read((char *) &d, sizeof(int64_t));
            shape.push_back(d);
            n = d;
            if (dim>1) {
                fin.read((char *) &d, sizeof(int64_t));
                shape.push_back(d);
                n *= d;
            }
            if (dim>2) {
                fin.read((char *)&d, sizeof(int64_t));
                shape.push_back(d);
                n *= d;
            }
            double value;
            image.reserve(n);
            while (!fin.eof()){
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
            if (dim>1) {
                getline(reading_file, reading_line_buffer); 
                shape.push_back(atoi(reading_line_buffer.c_str()));
                n *= shape[1];
            }
            if (dim>2) {
                getline(reading_file, reading_line_buffer);
                shape.push_back(atoi(reading_line_buffer.c_str()));
                n *= shape[2];
            }
            image.reserve(n);
            double value;
            while(!reading_file.eof()){
                getline(reading_file, reading_line_buffer);
                value = atof(reading_line_buffer.c_str());
                if (value != -1) {
                    image.push_back(value);
                }else{
                    image.push_back(DBL_MAX);
                }
            }
            reading_file.close();
            return;
		}
        case NUMPY: {
            vector<unsigned long> _shape;
            try{
                npy::LoadArrayFromNumpy(filename.c_str(), _shape, image);
            } catch (...) {
                cerr << "The data type of an numpy array should be numpy.float64." << endl;
                exit(-2);
            }
            uint8_t dim = shape.size();
            for (uint32_t i : _shape) {
                shape.push_back(i);
            }
            return;
        }
    }
}

void printResult(index_t dim, const CubicalGridComplex& cgc0, const CubicalGridComplex& cgc1, const CubicalGridComplex& cgcComp, 
    const vector<vector<Pair>>& pairs0, const vector<vector<Pair>>& pairs1, const vector<vector<Pair>>& pairsComp, 
    const vector<vector<Match>>& matches, unordered_map<index_t, bool>& isMatched0, unordered_map<index_t, bool>& isMatched1) {
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 0:" << endl; cout << endl; cgc0.printImage();
    cout << "pairs:" << endl;
    for (index_t d = 0; d < dim; d++) {
        cout << "in dimension " << d << ":" << endl;
        for (auto &pair : pairs0[d]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 1" << endl; cout << endl; cgc1.printImage();
    cout << "pairs:" << endl;
    for (index_t d = 0; d < dim; d++) {
        cout << "in dimension " << d << ":" << endl;
        for (auto &pair : pairs1[d]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Comparison" << endl; cout << endl; cgcComp.printImage();
    cout << "pairs:" << endl;
    for (index_t d = 0; d < dim; d++) {
        cout << "in dimension " << d << ":" << endl;
        for (auto &pair : pairsComp[d]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;

    cout << "Betti Matching:" << endl;
    for (index_t d = 0; d < dim; d++) {
        cout << "in dimension " << d << ":" << endl;
        for (auto &match : matches[d]) {
            match.print();
        }
        cout << endl;
    }
    cout << "unmatched in Input 0:" << endl;
    for (index_t d = 0; d < dim; d++) {
        cout << "dim " << d << ":" << endl;
        for (auto &pair : pairs0[d]) {
            if (!isMatched0[cgc0.getCubeIndex(pair.birth)]) {
                pair.print(); cout << endl;
            }
        }
        cout << endl;
    }
    cout << "unmatched in Input 2:" << endl;
    for (index_t d = 0; d < dim; d++) {
        cout << "dim " << d << ":" << endl;
        for (auto &pair : pairs1[d]) {
            if (!isMatched1[cgc1.getCubeIndex(pair.birth)]) {
                pair.print(); cout << endl;
            }
        }
        cout << endl;
    }
}