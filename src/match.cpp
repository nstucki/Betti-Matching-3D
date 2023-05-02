#include "match.h"
#include <iostream>

WriteMatch::WriteMatch(WritePair& _pair0, WritePair& _pair1) : pair0(_pair0), pair1(_pair1) {}

void WriteMatch::print() {
    cout << "paired "; pair0.print(); cout << " with "; pair1.print(); cout << endl;
}

Match::Match(vector<WritePair> (&_pairs_0)[3], vector<WritePair> (&_pairs_1)[3], vector<WritePair> (&_pairs_comp)[3], 
vector<WritePair> (&_pairs_im_0)[3], vector<WritePair> (&_pairs_im_1)[3], vector<WriteMatch> (&_matches)[3]) : 
pairs_0(_pairs_0), pairs_1(_pairs_1), pairs_comp(_pairs_comp), pairs_im_0(_pairs_im_0), pairs_im_1(_pairs_im_1), matches(_matches) {}

void Match::compute_matching() {
    for (uint8_t dim = 0; dim < 3; dim++) {
        for (WritePair& pair_comp : pairs_comp[dim]) {
            for (WritePair& pair_im_0 : pairs_im_0[dim]) {
                if (pair_im_0.death == pair_comp.death) {
                    for (WritePair& pair_0 : pairs_0[dim]) {
                        if (pair_0.birth == pair_im_0.birth) {
                            for (WritePair& pair_im_1 : pairs_im_1[dim]) {
                                if (pair_im_1.death == pair_comp.death) {
                                    for (WritePair& pair_1 : pairs_1[dim]) {
                                        if (pair_1.birth == pair_im_1.birth) {
                                            matches[dim].push_back(WriteMatch(pair_0, pair_1));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}