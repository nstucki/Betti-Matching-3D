#include "top_dimension.h"
#include "dimension_1.h"
#include "dimension_0.h"
#include "match.h"
#include <iostream>

void print_image(CubicalGridComplex* cgc){
    double birth;
    for (int i = 0; i < cgc->shape[0]; i++){
        for (int j = 0; j < cgc->shape[1]; j++){
            for (int k = 0; k < cgc->shape[2]; k++){
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


int main(){
    vector<double> image_0 {1, 2, 3, 8, 9, 4, 7, 6, 5, 21, 22, 23, 28, 29, 24, 27, 26, 25, 11, 12, 13, 18, 19, 14, 17, 16, 15};
    //vector<double> image_0 {1,1,1,1,2,1,1,1,1,4,4,4,4,5,4,4,4,4,2,2,2,2,3,2,2,2,2};
    vector<uint32_t> shape_0 {3, 3, 3};

    vector<double> image_1 {11, 12, 13, 18, 19, 14, 17, 16, 15, 21, 22, 23, 28, 29, 24, 27, 26, 25, 1, 2, 3, 8, 9, 4, 7, 6, 5};
    //vector<double> image_1 {3,3,3,3,4,3,3,3,3,2,2,2,2,3,2,2,2,2,1,1,1,1,2,1,1,1,1};
    vector<uint32_t> shape_1 {3, 3, 3};

    assert (shape_0==shape_1);

    vector<double> image_comp;
    transform(image_0.begin(), image_0.end(), image_1.begin(), back_inserter(image_comp), [](int a, int b){return std::min(a,b);});

    CubicalGridComplex* cgc_0 = new CubicalGridComplex(image_0, shape_0);
    CubicalGridComplex* cgc_1 = new CubicalGridComplex(image_1, shape_1);
    CubicalGridComplex* cgc_comp = new CubicalGridComplex(image_comp, shape_0);

    vector<WritePair> pairs_0[3];
    vector<WritePair> pairs_1[3];
    vector<WritePair> pairs_comp[3];
    vector<WritePair> pairs_im_0[3];
    vector<WritePair> pairs_im_1[3];
    vector<Cube> ctr;

    TopDimension* TD_comp = new TopDimension(cgc_comp, pairs_comp[2]);
    TD_comp->compute_pairs(ctr);
    delete TD_comp;

    Dimension1* D1_comp = new Dimension1(cgc_comp, pairs_comp[1]);
    D1_comp->compute_pairs(ctr);
    delete D1_comp;

    Dimension1* D1_im_0 = new Dimension1(cgc_0, pairs_im_0[1]);
    D1_im_0->compute_pairs(ctr, true);
    delete D1_im_0;

    Dimension1* D1_im_1 = new Dimension1(cgc_1, pairs_im_1[1]);
    D1_im_1->compute_pairs(ctr, true);
    delete D1_im_1;

    TopDimensionImage* TDI_0 = new TopDimensionImage(cgc_0, cgc_comp, pairs_0[2], pairs_im_0[2]);
    TDI_0->compute_pairs(ctr);
    delete TDI_0;

    Dimension1* D1_0 = new Dimension1(cgc_0, pairs_0[1]);
    D1_0->compute_pairs(ctr);
    delete D1_0;

    TopDimensionImage* TDI_1 = new TopDimensionImage(cgc_1, cgc_comp, pairs_1[2], pairs_im_1[2]);
    TDI_1->compute_pairs(ctr);
    delete TDI_1;

    Dimension1* D1_1 = new Dimension1(cgc_1, pairs_1[1]);
    D1_1->compute_pairs(ctr);
    delete D1_1;

    Dimension0* D0_0 = new Dimension0(cgc_0, pairs_0[0]);
    D0_0->compute_pairs(ctr);
    delete D0_0;

    Dimension0* D0_1 = new Dimension0(cgc_1, pairs_1[0]);
    D0_1->compute_pairs(ctr);
    delete D0_1;

    Dimension0Image* D0I = new Dimension0Image(cgc_0, cgc_1, cgc_comp, pairs_im_0[0], pairs_im_1[0], pairs_comp[0]);
    D0I->compute_pairs(ctr);
    delete D0I;

    vector<WriteMatch> matches[3];

    Match* match = new Match(pairs_0, pairs_1, pairs_comp, pairs_im_0, pairs_im_1, matches);
    match->compute_matching();
    
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Image 0:" << endl; cout << endl; print_image(cgc_0);
    cout << "pairs:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (auto pair : pairs_0[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Image 1" << endl; cout << endl; print_image(cgc_1);
    cout << "pairs:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (auto pair : pairs_1[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Image comp" << endl; cout << endl; print_image(cgc_comp);
    cout << "pairs:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (auto pair : pairs_comp[dim]) {
            pair.print(); cout << endl;
        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------------" << endl;

    //cout << "im 0" << endl;
    //for (auto wp : wp_im_0) {
    //    wp.print();
    //}
    //cout << "------------------------------" << endl;

    //cout << "im 1" << endl;
    //for (auto wp : wp_im_1) {
    //    wp.print();
    //}
    //cout << "------------------------------" << endl;

    cout << "Betti Matching:" << endl;
    for (uint8_t dim = 0; dim < 3; dim++) {
        cout << "dim " << unsigned(dim) << ":" << endl;
        for (WriteMatch match : matches[dim]) {
            match.print();
        }
        cout << endl;
    }
}