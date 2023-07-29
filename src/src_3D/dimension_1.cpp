#include "dimension_1.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>

using namespace std::chrono;


Dimension1::Dimension1(const CubicalGridComplex* const _cgc0, const CubicalGridComplex* const _cgc1, 
						const CubicalGridComplex* const _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), pairs0(_pairs0), pairs1(_pairs1),
						pairsComp(_pairsComp), matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1) {}

void Dimension1::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	if (config.verbose) { cout << "comoputing dimension 1 ... "; }
	auto start = high_resolution_clock::now();

	computePairsComp(ctrComp);

	#ifdef USE_CLEARING_DIM_0
	vector<Cube> ctrImage = ctrComp;
	ctrComp.clear();
	assembleColumnsToReduce(cgcComp, ctrComp);
	#endif
	
	pivotColumnIndex.clear();
	cache.clear();
	computePairs(ctr0, 0);

	#ifdef USE_CLEARING_DIM_0
	ctr0.clear();
	assembleColumnsToReduce(cgc0, ctr0);
	#endif
	
	pivotColumnIndex.clear();
	cache.clear();
	computePairs(ctr1, 1);

	#ifdef USE_CLEARING_DIM_0
	ctr1.clear();
	assembleColumnsToReduce(cgc1, ctr1);
	#endif
	
	#ifdef USE_CLEARING_DIM_0
	pivotColumnIndex.clear();
	cache.clear();
	computeImagePairs(ctrImage, 0);
	#else
	pivotColumnIndex.clear();
	cache.clear();
	computeImagePairs(ctrComp, 0); 
	#endif

	#ifdef USE_CLEARING_DIM_0
	pivotColumnIndex.clear();
	cache.clear();
	computeImagePairs(ctrImage, 1);
	#else
	pivotColumnIndex.clear();
	cache.clear();
	computeImagePairs(ctrComp, 1); 
	#endif
	
	computeMatching();
	
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }
}

void Dimension1::computePairsComp(vector<Cube>& ctr) {
	index_t ctrSize = ctr.size();
	pivotColumnIndex.reserve(ctrSize);	
	cache.reserve(min(config.cacheSize, ctrSize));
	BoundaryEnumerator faces = BoundaryEnumerator(cgcComp);
	queue<index_t> cachedColumnIdx;
	index_t numRecurse;
	index_t j;
	bool shouldClear = false;
	for(index_t i = 0; i < ctrSize; i++) {
		CubeQueue workingBoundary;
		j = i;
		numRecurse = 0;
		while (true) {
			bool cacheHit = false;
			if (i != j) {
                auto findCb = cache.find(j);
                if (findCb != cache.end()) {
                    cacheHit = true;
                    auto cachedBoundary = findCb -> second;
                    while (!cachedBoundary.empty()) {
                        workingBoundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				faces.setBoundaryEnumerator(ctr[j]);
				while (faces.hasNextFace()) { workingBoundary.push(faces.getNextFace()); }
			}
			Cube pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair -> second;
					numRecurse++;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, workingBoundary);
						cachedColumnIdx.push(i);
						if (cachedColumnIdx.size() > config.cacheSize) {
							cache.erase(cachedColumnIdx.front());
							cachedColumnIdx.pop();
						}
                    }
					pivotColumnIndex[pivot.index] = i;
					if (pivot.birth != ctr[i].birth) {
						pairsComp.push_back(Pair(pivot, ctr[i]));
						isMatchedComp.emplace(ctr[i].index, true);
					}
					break;
				}
			} else {
				ctr[i].index = NONE;
				shouldClear = true;
				break;
			}
		}
	}
	if (shouldClear) {
		auto newEnd = remove_if(ctr.begin(), ctr.end(), [](const Cube& cube) { return cube.index == NONE; });
		ctr.erase(newEnd, ctr.end());
	}	
}

void Dimension1::computePairs(const vector<Cube>& ctr, uint8_t k) {
	const CubicalGridComplex* const cgc = (k == 0) ? cgc0 : cgc1;
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<uint64_t, Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;

	index_t ctrSize = ctr.size();
	pivotColumnIndex.reserve(ctrSize);
	cache.reserve(min(config.cacheSize, ctrSize));
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	queue<index_t> cachedColumnIdx;
	index_t numRecurse;
	index_t j;	
	for(index_t i = 0; i < ctrSize; i++) {
		CubeQueue workingBoundary;
		j = i;
		numRecurse = 0;
		while (true) {
			bool cacheHit = false;
			if (i != j) {
                auto findCb = cache.find(j);
                if (findCb != cache.end()) {
                    cacheHit = true;
                    auto cachedBoundary = findCb -> second;
                    while (!cachedBoundary.empty()) {
                        workingBoundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				faces.setBoundaryEnumerator(ctr[j]);
				while (faces.hasNextFace()) { workingBoundary.push(faces.getNextFace()); }
			}
			Cube pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair -> second;
					numRecurse++;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, workingBoundary);
						cachedColumnIdx.push(i);
						if (cachedColumnIdx.size() > config.cacheSize) {
							cache.erase(cachedColumnIdx.front());
							cachedColumnIdx.pop();
						}
                    }
					pivotColumnIndex[pivot.index] = i;
					if (pivot.birth != ctr[i].birth) {
						pairs.push_back(Pair(pivot, ctr[i]));
						matchMap.emplace(pivot.index, pairs.back());
					}
					break;
				}
			} else { break; }
		}
	}
}

void Dimension1::computeImagePairs(const vector<Cube>& ctr, uint8_t k) {
	const CubicalGridComplex* const cgc = (k == 0) ? cgc0 : cgc1;
	unordered_map<uint64_t, uint64_t>& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;

	index_t ctrSize = ctr.size();
	pivotColumnIndex.reserve(ctrSize);
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMapIm.reserve(pairsComp.size());
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	queue<index_t> cachedColumnIdx;
	index_t numRecurse;
	index_t j;
	for (index_t i = 0; i < ctrSize; i++) {
		CubeQueue workingBoundary;
		j = i;
		numRecurse = 0;
		while (true) {
			bool cacheHit = false;
			if (i != j) {
                auto findCb = cache.find(j);
                if (findCb != cache.end()) {
                    cacheHit = true;
                    auto cachedBoundary = findCb -> second;
                    while (!cachedBoundary.empty()) {
                        workingBoundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				faces.setBoundaryEnumerator(ctr[j]);
				while (faces.hasNextFace()) { workingBoundary.push(faces.getNextFace()); }
			}
			Cube pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
					numRecurse++;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, workingBoundary);
						cachedColumnIdx.push(i);
						if (cachedColumnIdx.size() > config.cacheSize) {
							cache.erase(cachedColumnIdx.front());
							cachedColumnIdx.pop();
						}
                    }
					pivotColumnIndex[pivot.index] = i;
					if (isMatchedComp[ctr[i].index]) { matchMapIm.emplace(ctr[i].index, pivot.index); }
					break;
				}
			} else { break; }
		}
	}
}

void Dimension1::computeMatching() {
	uint64_t birthIndex0;
	uint64_t birthIndex1;
	for (Pair& pair : pairsComp) {
		auto find0 = matchMapIm0.find(pair.death.index);
		auto find1 = matchMapIm1.find(pair.death.index);
		if (find0 != matchMapIm0.end() && find1 != matchMapIm1.end()) {
			birthIndex0 = find0->second;
			birthIndex1 = find1->second;
			auto find0 = matchMap0.find(birthIndex0);
			auto find1 = matchMap1.find(birthIndex1);
			if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
				matches.push_back(Match(find0->second, find1->second));
				isMatched0.emplace(find0->second.birth.index, true);
				isMatched1.emplace(find1->second.birth.index, true);
			}
		}
	}
}

void Dimension1::assembleColumnsToReduce(const CubicalGridComplex* const cgc, vector<Cube>& ctr) const {
	ctr.reserve(cgc->getNumberOfCubes(1));
	value_t birth;
	Cube cube;
	for (index_t x = 0; x < cgc->shape[0]; x++) {
		for (index_t y = 0; y < cgc->shape[1]; y++) {
			for (index_t z = 0; z < cgc->shape[2]; z++) {
				for (uint8_t type = 0; type < 3; type++) {
					birth = cgc->getBirth(x, y, z, type, 1);
					if (birth < config.threshold) {
						cube = Cube(birth, x, y, z, type);
						auto find = pivotColumnIndex.find(cube.index);
						if (find != pivotColumnIndex.end()) { ctr.push_back(cube); }
					}	
				}				
			}
		}
	}
	sort(ctr.begin(), ctr.end(), CubeComparator());
}

Cube Dimension1::popPivot(CubeQueue& column) const {
    if (column.empty()) { return Cube(); } else {
        Cube pivot = column.top();
        column.pop();
        while (!column.empty() && column.top() == pivot) {
            column.pop();
            if (column.empty()) {return Cube(); }
            else {
                pivot = column.top();
                column.pop();
            }
        }
        return pivot;
    }
}

Cube Dimension1::getPivot(CubeQueue& column) const {
	Cube result = popPivot(column);
	if (result.index != NONE) { column.push(result); }
	return result;
}

void Dimension1::addCache(index_t i, CubeQueue& workingBoundary) {
	CubeQueue cleanWb;
	while (!workingBoundary.empty()) {
		Cube c = workingBoundary.top();
		workingBoundary.pop();
		if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
		else { cleanWb.push(c); }
	}
	cache.emplace(i, cleanWb);
}