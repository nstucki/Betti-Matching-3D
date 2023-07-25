#pragma once

#include "data_structures.h"

using namespace std;


class ComputeMatching {
    public:
    vector<Pair> (&pairs_0)[3];
    vector<Pair> (&pairs_1)[3];
    vector<Pair> (&pairs_comp)[3];
    vector<Pair> (&pairs_im_0)[3];
    vector<Pair> (&pairs_im_1)[3];
    vector<Match> (&matches)[3];
    vector<Pair> (&unmatched_0)[3];
    vector<Pair> (&unmatched_1)[3];

    ComputeMatching(vector<Pair> (&_pairs_0)[3], vector<Pair> (&_pairs_1)[3], vector<Pair> (&_pairs_comp)[3], 
        vector<Pair> (&_pairs_im_0)[3], vector<Pair> (&_pairs_im_1)[3], vector<Match> (&_matches)[3], 
        vector<Pair> (&_unmatched_0)[3], vector<Pair> (&_unmatched_1)[3]);
    void compute_matching();
};