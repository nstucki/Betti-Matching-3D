#include "top_dimension.h"
#include "dimension_1.h"
#include "dimension_0.h"
#include "match.h"
#include "npy.hpp"
#include <iostream>
#include <sstream> 
#include <fstream>
#include <cfloat>
#include <chrono>

using namespace std::chrono;

void print_image(CubicalGridComplex* cgc) {
    double birth;
    for (int i = 0; i < cgc->shape[0]; i++) {
        for (int j = 0; j < cgc->shape[1]; j++) {
            for (int k = 0; k < cgc->shape[2]; k++) {
                birth = cgc->getBirth(i,j,k);
                if (birth < 10){
                    cout << ' ' << birth << ' ';
                }else{
                    cout << birth << ' ';
                }            
            }
            cout << '\n';
        }
        cout << '\n';
    }
}

void tokenize(string const &str, const char delim, vector<string> &out) { 
    stringstream ss(str); 
    string s; 
    while (getline(ss, s, delim)) { 
        out.push_back(s); 
    } 
} 

void read_image(string const &filename, file_format const &format, vector<double> &image, vector<uint32_t> &shape) {
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
            assert (dim < 4);
            for (uint32_t i : _shape) {
                shape.push_back(i);
            }
            return;
        }
    }
}

void write_matched_pairs(const string &filename, vector<WriteMatch> (&matched)[3], CubicalGridComplex* cgc_0, CubicalGridComplex* cgc_1) {
    ofstream writing_file;
    vector<uint32_t> coordinates;
	if(filename.find(".csv")!= string::npos) {
		writing_file.open(filename, ios::out);
		if(!writing_file.is_open()){
			cerr << " error: open "<< filename << " failed! " << endl;
			exit(-3);
		}
        for (uint8_t dim = 0; dim < 3; dim++) {
            writing_file << unsigned(dim) << endl;
            for (WriteMatch &match : matched[dim]) {
                coordinates = cgc_0->ParentVoxel(match.pair_0.birth);
                writing_file << "(" << match.pair_0.birth.birth << " " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " , ";
                coordinates = cgc_0->ParentVoxel(match.pair_0.death);
                writing_file << match.pair_0.death.birth << " " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << ") - (";
                coordinates = cgc_1->ParentVoxel(match.pair_1.birth);
                writing_file << match.pair_1.birth.birth << " " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " , ";
                coordinates = cgc_1->ParentVoxel(match.pair_1.death);
                writing_file << match.pair_1.death.birth << " " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << ")" << endl;
            }
        }
		writing_file.close();
	}
}

void write_unmatched_pairs(const string &filename, vector<WritePair> unmatched[3], CubicalGridComplex* cgc) {
    ofstream writing_file;
    vector<uint32_t> coordinates;
    if(filename.find(".csv")!= string::npos) {
		writing_file.open(filename, ios::out);
		if(!writing_file.is_open()){
			cerr << " error: open " << filename << " failed! " << endl;
			exit(-3);
		}
        for (uint8_t dim = 0; dim < 3; dim++) {
            writing_file << unsigned(dim) << endl;
            for (WritePair &pair : unmatched[dim]) {
                coordinates = cgc->ParentVoxel(pair.birth);
                writing_file << "(" << pair.birth.birth << " " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " , ";
                coordinates = cgc->ParentVoxel(pair.death);
                writing_file << pair.death.birth << " " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << ")" << endl;
            }
        }
		writing_file.close();
	}
}

void print_result(CubicalGridComplex* cgc_0, CubicalGridComplex* cgc_1, CubicalGridComplex* cgc_comp, 
    vector<WritePair> (&pairs_0)[3], vector<WritePair> (&pairs_1)[3], vector<WritePair> (&pairs_comp)[3], 
    vector<WriteMatch> (&matched)[3], vector<WritePair> (&unmatched_0)[3], vector<WritePair> (&unmatched_1)[3]) {
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Image 0:" << endl; cout << endl; print_image(cgc_0);
    cout << "pairs:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (WritePair &pair : pairs_0[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Image 1" << endl; cout << endl; print_image(cgc_1);
    cout << "pairs:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (WritePair &pair : pairs_1[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Image comp" << endl; cout << endl; print_image(cgc_comp);
    cout << "pairs:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (WritePair &pair : pairs_comp[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;

    cout << "Betti Matching:" << endl;
    cout << "matched:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (WriteMatch &match : matched[dim]) {
            match.print();
        }
        cout << endl;
    }
    cout << "unmatched in 1st image:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (WritePair &pair : unmatched_0[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "unmatched in 2nd image:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (WritePair &pair : unmatched_1[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
}

void print_usage_and_exit(int exit_code) {
    cout << endl;
    cerr << "Usage: "
         << "Betti Matching "
         << "[options] [input_0_filename] [input_1_filename]" << endl
         << endl
         << "Options:" << endl
         << endl
         << "  --help, -h                      print this screen" << endl
         << "  --verbose, -v                   display detailed processing information" << endl
         << "  --min_recursion_to_cache, -mc   minimum number of recursion for a reduced column to be cached" << endl
         << "  --cache_size, -c	          maximum number of reduced columns to be cached" << endl
         << "  --matched, -m                   name of the file containing matched pairs" << endl
         << "  --unmatchedmatched_0, -u0       name of the file containing unmatched pairs of input 0" << endl
         << "  --unmatchedmatched_1, -u1       name of the file containing unmatched pairs of input 1" << endl
         << "  --print, -p                     print result in console" << endl
         << endl;
	exit(exit_code);
}

int main(int argc, char** argv){
    Config config;
    for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help" || arg == "-h") {
			print_usage_and_exit(0);
		} else if (arg == "--verbose" || arg == "-v") {
			config.verbose = true;
		} else if (arg == "--matched" || arg == "-m") {
			config.matched_filename = string(argv[++i]);
        } else if (arg == "--unmatched_0" || arg == "-u0") {
			config.unmatched_0_filename = string(argv[++i]);
        } else if (arg == "--unmatched_1" || arg == "-u1") {
			config.unmatched_1_filename = string(argv[++i]);
		} else if (arg == "--min_recursion_to_cache" || arg == "-mc"){
            config.min_recursion_to_cache = stoi(argv[++i]);
		} else if (arg == "--cache_size" || arg == "-c"){
            config.cache_size = stoi(argv[++i]);
		} else if (arg == "--print" || arg == "-p"){
			config.print = true;
		} else {
            if (config.filename_0.empty()) {
                config.filename_0 = argv[i];
            } else {
                config.filename_1 = argv[i];
            }
		}
	}
    if (config.filename_0.empty()) {print_usage_and_exit(-1);} 
    if(config.filename_0.find(".txt")!= std::string::npos) {
		config.format_0 = PERSEUS;
	} else if(config.filename_0.find(".npy")!= std::string::npos) {
		config.format_0 = NUMPY;
	} else if(config.filename_0.find(".complex")!= std::string::npos) {
		config.format_0 = DIPHA;
	} else {
		cerr << "unknown input file format! (the filename extension should be .txt/.npy/.complex): " << config.filename_0 << endl;
		exit(-1);
	}
    ifstream file_stream_0(config.filename_0);
	if (!config.filename_0.empty() && file_stream_0.fail()) {
		cerr << "couldn't open file " << config.filename_0 << endl;
		exit(-1);
	}
    if (config.filename_1.empty()) {print_usage_and_exit(-1);}
    if(config.filename_1.find(".txt")!= std::string::npos) {
		config.format_1 = PERSEUS;
	} else if(config.filename_1.find(".npy")!= std::string::npos) {
		config.format_1 = NUMPY;
	} else if(config.filename_1.find(".complex")!= std::string::npos) {
		config.format_1 = DIPHA;
	} else {
		cerr << "unknown input file format! (the filename extension should be .txt/.npy/.complex): " << config.filename_1 << endl;
		exit(-1);
	}
    ifstream file_stream_1(config.filename_1);
	if (!config.filename_1.empty() && file_stream_1.fail()) {
		cerr << "couldn't open file " << config.filename_1 << endl;
		exit(-1);
	}

    if (config.verbose) {cout << "reading images ... ";}
    auto start = high_resolution_clock::now();
    vector<double> image_0;
    vector<uint32_t> shape_0;
    vector<double> image_1;
    vector<uint32_t> shape_1;
    read_image(config.filename_0, config.format_0, image_0, shape_0);
    read_image(config.filename_1, config.format_1, image_1, shape_1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    auto duration_total = duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}
    
    assert (shape_0 == shape_1);
    vector<double> image_comp;
    transform(image_0.begin(), image_0.end(), image_1.begin(), back_inserter(image_comp), [](double a, double b){return min(a,b);});
    CubicalGridComplex* cgc_0 = new CubicalGridComplex(move(image_0), shape_0);
    CubicalGridComplex* cgc_1 = new CubicalGridComplex(move(image_1), shape_1);
    CubicalGridComplex* cgc_comp = new CubicalGridComplex(move(image_comp), shape_0);
    vector<WritePair> pairs_0[3];
    vector<WritePair> pairs_1[3];
    vector<WritePair> pairs_comp[3];
    vector<WritePair> pairs_im_0[3];
    vector<WritePair> pairs_im_1[3];  
    vector<Cube> ctr;

    if (config.verbose) {cout << "comoputing dim 2 of comparison ... ";}
    start = high_resolution_clock::now();
    TopDimension* TD_comp = new TopDimension(cgc_comp, pairs_comp[2]);
    TD_comp->compute_pairs(ctr);
    delete TD_comp;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose){ cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 1 of comparison ... ";}
    start = high_resolution_clock::now();
    Dimension1* D1_comp = new Dimension1(cgc_comp, pairs_comp[1], config);
    D1_comp->compute_pairs(ctr);
    delete D1_comp;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 1 of im_0 ... ";}
    start = high_resolution_clock::now();
    Dimension1* D1_im_0 = new Dimension1(cgc_0, pairs_im_0[1], config);
    D1_im_0->compute_pairs(ctr, true);
    delete D1_im_0;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 1 of im_1 ... ";}
    start = high_resolution_clock::now();
    Dimension1* D1_im_1 = new Dimension1(cgc_1, pairs_im_1[1], config);
    D1_im_1->compute_pairs(ctr, true);
    delete D1_im_1;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 2 of input_0 and im_0 ... ";}
    start = high_resolution_clock::now();
    TopDimensionImage* TDI_0 = new TopDimensionImage(cgc_0, cgc_comp, pairs_0[2], pairs_im_0[2]);
    TDI_0->compute_pairs(ctr);
    delete TDI_0;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 1 of input_0 ... ";}
    start = high_resolution_clock::now();
    Dimension1* D1_0 = new Dimension1(cgc_0, pairs_0[1], config);
    D1_0->compute_pairs(ctr);
    delete D1_0;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 2 of input_1 and im_1 ... ";}
    start = high_resolution_clock::now();
    TopDimensionImage* TDI_1 = new TopDimensionImage(cgc_1, cgc_comp, pairs_1[2], pairs_im_1[2]);
    TDI_1->compute_pairs(ctr);
    delete TDI_1;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 1 of input_1 ... ";}
    start = high_resolution_clock::now();
    Dimension1* D1_1 = new Dimension1(cgc_1, pairs_1[1], config);
    D1_1->compute_pairs(ctr);
    delete D1_1;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 0 of input_0 ... ";}
    start = high_resolution_clock::now();
    Dimension0* D0_0 = new Dimension0(cgc_0, pairs_0[0]);
    D0_0->compute_pairs(ctr);
    delete D0_0;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 0 of input_1 ... ";}
    start = high_resolution_clock::now();
    Dimension0* D0_1 = new Dimension0(cgc_1, pairs_1[0]);
    D0_1->compute_pairs(ctr);
    delete D0_1;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing dim 0 of comparison, im_0 and im_1 ... ";}
    start = high_resolution_clock::now();
    Dimension0Image* D0I = new Dimension0Image(cgc_0, cgc_1, cgc_comp, pairs_im_0[0], pairs_im_1[0], pairs_comp[0]);
    D0I->compute_pairs(ctr);
    delete D0I;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    if (config.verbose) {cout << "comoputing the matching ... ";}
    start = high_resolution_clock::now();
    vector<WriteMatch> matched[3];
    vector<WritePair> unmatched_0[3];
    vector<WritePair> unmatched_1[3];
    Match* match = new Match(pairs_0, pairs_1, pairs_comp, pairs_im_0, pairs_im_1, matched, unmatched_0, unmatched_1);
    match->compute_matching();
    delete match;
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    duration_total += duration;
    if (config.verbose) {cout << "took " << duration.count() << " ms" << endl;}

    write_matched_pairs(config.matched_filename, matched, cgc_0, cgc_1);
    write_unmatched_pairs(config.unmatched_0_filename, unmatched_0, cgc_0);
    write_unmatched_pairs(config.unmatched_1_filename, unmatched_1, cgc_1);

    if (config.verbose) {cout << "total computation time: " << duration_total.count() << " ms" << endl;}

    if (config.print) {
        print_result(cgc_0, cgc_1, cgc_comp, pairs_0, pairs_1, pairs_comp, matched, unmatched_0, unmatched_1);
    }
}