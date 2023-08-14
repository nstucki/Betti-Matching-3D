#include "inter_dimensions.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>
#include <algorithm>

using namespace std::chrono;


InterDimensions::InterDimensions(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, const CubicalGridComplex& _cgcComp, 
									const Config& _config, vector<vector<Pair>>& _pairs0, vector<vector<Pair>>& _pairs1, 
									vector<vector<Pair>>& _pairsComp, vector<vector<Match>>& _matches, 
									unordered_map<index_t, bool>& _isMatched0, unordered_map<index_t, bool>& _isMatched1) : 
									cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), pairs0(_pairs0), pairs1(_pairs1),
									pairsComp(_pairsComp), matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1) {
	computeDim = cgc0.dim-2;
}

void InterDimensions::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	vector<Cube> ctrImage;

	while (computeDim > 0) {
		if (config.verbose) { cout << "comoputing dimension " << computeDim << " ... "; }
        auto start = high_resolution_clock::now();
		
		computePairsComp(ctrComp);
		if (computeDim > 1) { 
			ctrImage = ctrComp;
			assembleColumnsToReduce(cgcComp, ctrComp);
		}
		
		computePairs(ctr0, 0);
		if (computeDim > 1) { assembleColumnsToReduce(cgc0, ctr0); }
		
		computePairs(ctr1, 1);
		if (computeDim > 1) { assembleColumnsToReduce(cgc1, ctr1); }
		
		if (computeDim > 1) { 
			computeImagePairs(ctrImage, 0); 
			computeImagePairs(ctrImage, 1); 
		} else { 
			computeImagePairs(ctrComp, 0); 
			computeImagePairs(ctrComp, 1); 
		}
		
		computeMatching();
		
		auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }

		--computeDim;
	}
}

void InterDimensions::computePairsComp(vector<Cube>& ctr) {
	index_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);	
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMapComp.clear();
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
				faces.setBoundaryEnumerator(ctr[j], computeDim+1);
				while (faces.hasNextFace()) { workingBoundary.push(faces.getNextFace()); }
			}
			Cube pivot = getPivot(workingBoundary);
			if (pivot.coordinates[0] != NONE) {
				auto pair = pivotColumnIndex.find(cgcComp.getCubeIndex(pivot));
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
					pivotColumnIndex[cgcComp.getCubeIndex(pivot)] = i;
					if (pivot.birth != ctr[i].birth) {
						pairsComp[computeDim].push_back(Pair(pivot, ctr[i]));
						matchMapComp.emplace(cgcComp.getCubeIndex(ctr[i]), true);
					}
					break;
				}
			} else {
				ctr[i].coordinates[0] = NONE;
				shouldClear = true;
				break;
			}
		}
	}
	if (shouldClear) {
		auto newEnd = remove_if(ctr.begin(), ctr.end(), [](const Cube& cube) { return cube.coordinates[0] == NONE; });
		ctr.erase(newEnd, ctr.end());
	}	
}

void InterDimensions::computePairs(const vector<Cube>& ctr, uint8_t k) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	vector<vector<Pair>>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t, Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;

	index_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMap.clear();
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
				faces.setBoundaryEnumerator(ctr[j], computeDim+1);
				while (faces.hasNextFace()) { workingBoundary.push(faces.getNextFace()); }
			}
			Cube pivot = getPivot(workingBoundary);
			if (pivot.coordinates[0] != NONE) {
				auto pair = pivotColumnIndex.find(cgc.getCubeIndex(pivot));
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
					pivotColumnIndex[cgc.getCubeIndex(pivot)] = i;
					if (pivot.birth != ctr[i].birth) {
						pairs[computeDim].push_back(Pair(pivot, ctr[i]));
						matchMap.emplace(cgc.getCubeIndex(pivot), pairs[computeDim].back());
					}
					break;
				}
			} else { break; }
		}
	}
}

void InterDimensions::computeImagePairs(const vector<Cube>& ctr, uint8_t k) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	unordered_map<index_t, Cube>& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;

	index_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMapIm.clear();
	matchMapIm.reserve(pairsComp[computeDim].size());
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
				faces.setBoundaryEnumerator(ctr[j], computeDim+1);
				while (faces.hasNextFace()) { workingBoundary.push(faces.getNextFace()); }
			}
			Cube pivot = getPivot(workingBoundary);
			if (pivot.coordinates[0] != NONE) {
				auto pair = pivotColumnIndex.find(cgc.getCubeIndex(pivot));
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
					pivotColumnIndex[cgc.getCubeIndex(pivot)] = i;
					if (matchMapComp[cgcComp.getCubeIndex(ctr[i])]) { matchMapIm.emplace(cgcComp.getCubeIndex(ctr[i]), pivot); }
					break;
				}
			} else { break; }
		}
	}
}

void InterDimensions::computeMatching() {
	Cube birth0;
	Cube birth1;
	Pair pair0;
	Pair pair1;
	for (Pair& pair : pairsComp[computeDim]) {
		auto find0 = matchMapIm0.find(cgcComp.getCubeIndex(pair.death));
		auto find1 = matchMapIm1.find(cgcComp.getCubeIndex(pair.death));
		if (find0 != matchMapIm0.end() && find1 != matchMapIm1.end()) {
			birth0 = find0->second;
			birth1 = find1->second;
			auto find0 = matchMap0.find(cgc0.getCubeIndex(birth0));
			auto find1 = matchMap1.find(cgc1.getCubeIndex(birth1));
			if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
				pair0 = find0->second;
				pair1 = find1->second;
				matches[computeDim].push_back(Match(pair0, pair1));
				isMatched0.emplace(cgc0.getCubeIndex(pair0.birth), true);
				isMatched1.emplace(cgc1.getCubeIndex(pair1.birth), true);
			}
		}
	}
}

void InterDimensions::assembleColumnsToReduce(const CubicalGridComplex& cgc, vector<Cube>& ctr) const {
	ctr.clear();
	ctr.reserve(cgc.getNumberOfCubes(computeDim));	
	CubeEnumerator cubeEnum(cgc, computeDim);
	Cube cube = cubeEnum.getNextCube();
	if (cube.birth < config.threshold) {
		auto find = pivotColumnIndex.find(cgc.getCubeIndex(cube));
		if (find != pivotColumnIndex.end()) { ctr.push_back(cube); }
	}
	while (cubeEnum.hasNextCube()) {
		cube = cubeEnum.getNextCube();
		if (cube.birth < config.threshold) {
			auto find = pivotColumnIndex.find(cgc.getCubeIndex(cube));
			if (find != pivotColumnIndex.end()) { ctr.push_back(cube); }
		}
	}
	sort(ctr.begin(), ctr.end(), CubeComparator());
}

Cube InterDimensions::popPivot(CubeQueue& column) const {
    if (column.empty()) { return Cube(); } 
	else {
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

Cube InterDimensions::getPivot(CubeQueue& column) const {
	Cube result = popPivot(column);
	if (result.coordinates[0] != NONE) { column.push(result); }
	return result;
}

void InterDimensions::addCache(index_t i, CubeQueue& workingBoundary) {
	CubeQueue cleanWb;
	while (!workingBoundary.empty()) {
		Cube c = workingBoundary.top();
		workingBoundary.pop();
		if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
		else { cleanWb.push(c); }
	}
	cache.emplace(i, cleanWb);
}