#pragma once

#include "write_pair.h"
#include <vector>

using namespace std;

class WriteMatch{
public:
	WritePair& pair0;
	WritePair& pair1;

    WriteMatch(WritePair& _pair0, WritePair& _pair1);
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

    Match(vector<WritePair> (&_pairs_0)[3], vector<WritePair> (&_pairs_1)[3], vector<WritePair> (&_pairs_comp)[3], vector<WritePair> (&_pairs_im_0)[3], vector<WritePair> (&_pairs_im_1)[3], vector<WriteMatch> (&_matches)[3]);
    void compute_matching();
};