#pragma once

#include "write_pair.h"
#include <vector>

using namespace std;

class WriteMatch{
public:
	WritePair& pair_0;
	WritePair& pair_1;

    WriteMatch(WritePair& _pair_0, WritePair& _pair_1);
	void print();
};

class Match{
    public:
    vector<WritePair> (&pairs_0)[3];
    vector<WritePair> (&pairs_1)[3];
    vector<WritePair> (&pairs_comp)[3];
    vector<WritePair> (&pairs_im_0)[3];
    vector<WritePair> (&pairs_im_1)[3];
    vector<WriteMatch> (&matches)[3];
    vector<WritePair> (&unmatched_0)[3];
    vector<WritePair> (&unmatched_1)[3];

    Match(vector<WritePair> (&_pairs_0)[3], vector<WritePair> (&_pairs_1)[3], vector<WritePair> (&_pairs_comp)[3], 
        vector<WritePair> (&_pairs_im_0)[3], vector<WritePair> (&_pairs_im_1)[3], vector<WriteMatch> (&_matches)[3], 
        vector<WritePair> (&_unmatched_0)[3], vector<WritePair> (&_unmatched_1)[3]);
    void compute_matching();
};