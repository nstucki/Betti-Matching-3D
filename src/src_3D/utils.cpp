#include "utils.h"
#include "npy.hpp"

#include <cfloat>
#include <cassert>

void readImage(const string& filename, const fileFormat& format, vector<double>& image, vector<index_t>& shape) {
    switch (format) {
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
            assert(dim < 4);
            uint64_t n;
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
            assert(dim < 4);
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
            while(!reading_file.eof()) {
                getline(reading_file, reading_line_buffer);
                value = atof(reading_line_buffer.c_str());
                if (value != -1) { image.push_back(value); }
                else { image.push_back(DBL_MAX); }
            }
            reading_file.close();
            return;
		}

        case NUMPY: {
            vector<unsigned long> _shape;
            try { npy::LoadArrayFromNumpy(filename.c_str(), _shape, image); } 
            catch (...) {
                cerr << "The data type of an numpy array should be numpy.float64." << endl;
                exit(-2);
            }
            uint8_t dim = shape.size();
            assert (dim < 4);
            for (uint32_t i : _shape) { shape.push_back(i); }
            return;
        }
    }
}

void tokenize(const string& str, const char delim, vector<string>& out) { 
    stringstream ss(str); 
    string s; 
    while (getline(ss, s, delim)) { out.push_back(s); } 
} 

void printResult(const CubicalGridComplex* const cgc0, const CubicalGridComplex* const cgc1, const CubicalGridComplex* const cgcComp, 
                    const vector<vector<Pair>>& pairs0, const vector<vector<Pair>>& pairs1, const vector<vector<Pair>>& pairsComp,
                    const vector<vector<Match>>& matches, 
                    vector<unordered_map<uint64_t, bool>>& isMatched0, vector<unordered_map<uint64_t, bool>>& isMatched1) {
    index_t count;
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 0:" << endl;
    if (cgc0->shape[0] < 10 && cgc0->shape[1] < 10 && cgc0->shape[2] < 10) { cgc0->printImage(); }
    for (uint8_t d = 0; d < 3; d++) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs0[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs0[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Input 1" << endl; 
    if (cgc1->shape[0] < 10 && cgc1->shape[1] < 10 && cgc1->shape[2] < 10) { cgc1->printImage(); }
    for (uint8_t d = 0; d < 3; d++) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairs1[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs1[d]) { pair.print(); cout << endl; }
        } else { cout << count << endl; }
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Comparison" << endl; 
    if (cgcComp->shape[0] < 10 && cgcComp->shape[1] < 10 && cgcComp->shape[2] < 10) { cgcComp->printImage(); }
    for (uint8_t d = 0; d < 3; d++) {
        cout << "dim " << unsigned(d) << ": ";
        count = pairsComp[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairsComp[d]) { pair.print(); cout << endl; }
        } else { cout << "number of pairs: " << count << endl; }
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;

    cout << "Betti Matching:" << endl;
    cout << "matched: " << endl;
    for (uint8_t d = 0; d < 3; d++) {
        cout << "dim " << unsigned(d) << ": ";
        count = matches[d].size();
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &match : matches[d]) { match.print(); }
        } else { cout << count << endl; }
    }
    cout << "unmatched in Input 0" << endl;
    for (uint8_t d = 0; d < 3; d++) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs0[d]) { if (!isMatched0[d][pair.birth.index]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs0[d]) { if (!isMatched0[d][pair.birth.index]) { pair.print(); cout << endl; } }
        } else { cout << count << endl; }
    }
    cout << "unmatched in Input 1" << endl;
    for (uint8_t d = 0; d < 3; d++) {
        cout << "dim " << unsigned(d) << ": ";
        count = 0;
        for (auto &pair : pairs1[d]) {if (!isMatched1[d][pair.birth.index]) { ++count; } }
        if (0 < count && count < 10) {
            cout << endl;
            for (auto &pair : pairs1[d]) { if (!isMatched1[d][pair.birth.index]) { pair.print(); cout << endl; } }
        } else { cout << count << endl; }
    }
}