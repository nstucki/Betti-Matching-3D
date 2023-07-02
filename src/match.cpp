#include "match.h"
#include <iostream>

WriteMatch::WriteMatch(WritePair &_pair_0, WritePair &_pair_1) : pair_0(_pair_0), pair_1(_pair_1) {}

void WriteMatch::print() {
    cout << " matched "; pair_0.print(); cout << " with "; pair_1.print(); cout << endl;
}

Match::Match(vector<WritePair> (&_pairs_0)[3], vector<WritePair> (&_pairs_1)[3], vector<WritePair> (&_pairs_comp)[3], 
    vector<WritePair> (&_pairs_im_0)[3], vector<WritePair> (&_pairs_im_1)[3], vector<WriteMatch> (&_matches)[3], 
    vector<WritePair> (&_unmatched_0)[3], vector<WritePair> (&_unmatched_1)[3]) : 
    pairs_0(_pairs_0), pairs_1(_pairs_1), pairs_comp(_pairs_comp), pairs_im_0(_pairs_im_0), 
    pairs_im_1(_pairs_im_1), matches(_matches), unmatched_0(_unmatched_0), unmatched_1(_unmatched_1) {}

void Match::compute_matching() {
    for (uint8_t dim = 0; dim < 3; dim++) {
        for (WritePair &pair_comp : pairs_comp[dim]) {
            for (WritePair &pair_im_0 : pairs_im_0[dim]) {
                if (pair_im_0.death == pair_comp.death) {
                    for (WritePair &pair_0 : pairs_0[dim]) {
                        if (pair_0.birth == pair_im_0.birth) {
                            for (WritePair &pair_im_1 : pairs_im_1[dim]) {
                                if (pair_im_1.death == pair_comp.death) {
                                    for (WritePair &pair_1 : pairs_1[dim]) {
                                        if (pair_1.birth == pair_im_1.birth) {
                                            matches[dim].push_back(WriteMatch(pair_0, pair_1));
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        }
        bool matched;
        for (WritePair &pair : pairs_0[dim]) {
            matched = false;
            for (WriteMatch &match : matches[dim]) {
                if (pair == match.pair_0) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                unmatched_0[dim].push_back(pair);
            }
        }
        for (WritePair &pair : pairs_1[dim]) {
            matched = false;
            for (WriteMatch &match : matches[dim]) {
                if (pair == match.pair_1) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                unmatched_1[dim].push_back(pair);
            }
        }
    }
}