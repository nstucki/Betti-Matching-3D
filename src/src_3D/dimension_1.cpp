#include "dimension_1.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>

using namespace std::chrono;


Dimension1::Dimension1(const CubicalGridComplex* const _cgc0, const CubicalGridComplex* const _cgc1, 
						const CubicalGridComplex* const _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<index_t, bool>& _isMatched0, unordered_map<index_t, bool>& _isMatched1) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), pairs0(_pairs0), pairs1(_pairs1),
						pairsComp(_pairsComp), matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1) {}

void Dimension1::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	if (config.verbose) { cout << "comoputing dimension 1 ... "; }
	auto start = high_resolution_clock::now();
	
	computePairsComp(ctrComp);
	//assembleColumnsToReduce(cgcComp, ctrComp);
	
	computePairs(ctr0, 0);
	//assembleColumnsToReduce(cgc0, ctr0);
	
	computePairs(ctr1, 1);
	//assembleColumnsToReduce(cgc1, ctr1);
	
	computeImagePairs(ctrComp, 0); 
	computeImagePairs(ctrComp, 1);
	
	computeMatching();
	
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }
}

void Dimension1::computePairsComp(vector<Cube>& ctr) {
	index_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);	
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	BoundaryEnumerator faces = BoundaryEnumerator(cgcComp);
	queue<index_t> cachedColumnIdx;
	index_t numRecurse;
	index_t j;
	bool shouldClear = false;
	for(index_t i = 0; i < ctrSize; i++) {
		CubeQue workingBoundary;
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
						matchMapComp.emplace(ctr[i].index, true);
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
		auto new_end = remove_if(ctr.begin(), ctr.end(), [](const Cube& edge) { return edge.index == NONE; });
		ctr.erase(new_end, ctr.end());
	}	
}

void Dimension1::computePairs(const vector<Cube>& ctr, uint8_t k) {
	const CubicalGridComplex* const cgc = (k == 0) ? cgc0 : cgc1;
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<uint64_t, Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;

	index_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	queue<index_t> cachedColumnIdx;
	index_t numRecurse;
	index_t j;	
	for(index_t i = 0; i < ctrSize; i++) {
		CubeQue workingBoundary;
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
			if (pivot.coordinates[0] != NONE) {
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
						matchMap.emplace(pivot.index, &pairs.back());
					}
					break;
				}
			} else { break; }
		}
	}
}

void Dimension1::computeImagePairs(const vector<Cube>& ctr, uint8_t k) {
	const CubicalGridComplex* const cgc = (k == 0) ? cgc0 : cgc1;
	unordered_map<uint64_t, Cube*>& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;

	index_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMapIm.reserve(pairsComp.size());
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	queue<index_t> cachedColumnIdx;
	index_t numRecurse;
	index_t j;
	for (index_t i = 0; i < ctrSize; i++) {
		CubeQue workingBoundary;
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
	for (auto& pair : pairsComp[computeDim]) {
		auto findIm0 = matchMapIm0.find(cgcComp.getCubeIndex(pair.death));
		auto findIm1 = matchMapIm1.find(cgcComp.getCubeIndex(pair.death));
		if (findIm0 != matchMapIm0.end() && findIm1 != matchMapIm1.end()) {
			birth0 = findIm0->second;
			birth1 = findIm1->second;
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
	ctr.reserve(cgc.getNumberOfCubes(computeDim-1));	
	CubeEnumerator cubeEnum(cgc, computeDim-1);
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

Cube Dimension1::popPivot(CubeQue& column) const {
    if (column.empty()) {
        return Cube();
    } else {
        Cube pivot = column.top();
        column.pop();
        while (!column.empty() && column.top() == pivot) {
            column.pop();
            if (column.empty())
                return Cube();
            else {
                pivot = column.top();
                column.pop();
            }
        }
        return pivot;
    }
}

Cube InterDimensions::getPivot(CubeQue& column) const {
	Cube result = popPivot(column);
	if (result.coordinates[0] != NONE) { column.push(result); }
	return result;
}

void InterDimensions::addCache(index_t i, CubeQue& workingBoundary) {
	CubeQue cleanWb;
	while (!workingBoundary.empty()) {
		auto c = workingBoundary.top();
		workingBoundary.pop();
		if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
		else { cleanWb.push(c); }
	}
	cache.emplace(i, cleanWb);
}