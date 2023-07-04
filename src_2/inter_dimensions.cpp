#include "inter_dimensions.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>

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
	vector<Cube> ctrIm;

	while (computeDim > 0) {
		if (config.verbose) { cout << "comoputing dimension " << computeDim << " ... "; }
        auto start = high_resolution_clock::now();
		
		computePairsComp(ctrComp);
		if (computeDim > 1) { 
			ctrIm = ctrComp;
			assembleColumnsToReduce(cgcComp, ctrComp);
		}
		
		computePairs(ctr0, 0);
		if (computeDim > 1) { assembleColumnsToReduce(cgc0, ctr0); }
		
		computePairs(ctr1, 1);
		if (computeDim > 1) { assembleColumnsToReduce(cgc1, ctr1); }
		
		if (computeDim > 1) { computeImagePairs(ctrIm, 0); computeImagePairs(ctrIm, 1); }
		else { computeImagePairs(ctrComp, 0); computeImagePairs(ctrComp, 1); }
		
		computeMatching();
		
		auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        if (config.verbose) { cout << "took " << duration.count() << " ms" << endl; }

		--computeDim;
	}
	
}

void InterDimensions::computePairsComp(vector<Cube>& ctr) {
	index_t ctrSize = ctr.size();
	pivot_column_index.clear();
	pivot_column_index.reserve(ctrSize);	
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMapComp.clear();
	BoundaryEnumerator faces = BoundaryEnumerator(cgcComp);
	queue<index_t> cached_column_idx;
	index_t numRecurse;
	index_t j;
	bool shouldClear = false;
	for(index_t i = 0; i < ctrSize; i++) {
		CubeQue working_boundary;
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
                        working_boundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				faces.setBoundaryEnumerator(ctr[j], computeDim+1);
				while (faces.hasNextFace()) {
					working_boundary.push(faces.getNextFace());
				}
			}
			Cube pivot = getPivot(working_boundary);
			if (pivot.coordinates[0] != NONE) {
				auto pair = pivot_column_index.find(cgcComp.getCubeIndex(pivot));
				if (pair != pivot_column_index.end()) {
					j = pair->second;
					numRecurse++;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, working_boundary);
						cached_column_idx.push(i);
						if (cached_column_idx.size() > config.cacheSize) {
							cache.erase(cached_column_idx.front());
							cached_column_idx.pop();
						}
                    }
					pivot_column_index[cgcComp.getCubeIndex(pivot)] = i;
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
	// create new vector instead?
	if (shouldClear) {
		auto new_end = remove_if(ctr.begin(), ctr.end(), [](const Cube& edge){ return edge.coordinates[0] == NONE; });
		ctr.erase(new_end, ctr.end());
	}
	
}

void InterDimensions::computePairs(const vector<Cube>& ctr, uint8_t k) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	vector<vector<Pair>>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<index_t,Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;

	index_t ctrSize = ctr.size();
	pivot_column_index.clear();
	pivot_column_index.reserve(ctrSize);
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMap.clear();
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	queue<index_t> cached_column_idx;
	index_t numRecurse;
	index_t j;	
	for(index_t i = 0; i < ctrSize; i++) {
		CubeQue working_boundary;
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
                        working_boundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				faces.setBoundaryEnumerator(ctr[j], computeDim+1);
				while (faces.hasNextFace()) {
					working_boundary.push(faces.getNextFace());
				}
			}
			Cube pivot = getPivot(working_boundary);
			if (pivot.coordinates[0] != NONE) {
				auto pair = pivot_column_index.find(cgc.getCubeIndex(pivot));
				if (pair != pivot_column_index.end()) {
					j = pair->second;
					numRecurse++;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, working_boundary);
						cached_column_idx.push(i);
						if (cached_column_idx.size() > config.cacheSize) {
							cache.erase(cached_column_idx.front());
							cached_column_idx.pop();
						}
                    }
					pivot_column_index[cgc.getCubeIndex(pivot)] = i;
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
	unordered_map<index_t,Cube>& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;

	index_t ctrSize = ctr.size();
	pivot_column_index.clear();
	pivot_column_index.reserve(ctrSize);
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMapIm.clear();
	matchMapIm.reserve(pairsComp[computeDim].size());
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	queue<index_t> cached_column_idx;
	index_t numRecurse;
	index_t j;
	for (index_t i = 0; i < ctrSize; i++) {
		CubeQue working_boundary;
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
                        working_boundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				faces.setBoundaryEnumerator(ctr[j], computeDim+1);
				while (faces.hasNextFace()) {
					working_boundary.push(faces.getNextFace());
				}
			}
			Cube pivot = getPivot(working_boundary);
			if (pivot.coordinates[0] != NONE) {
				auto pair = pivot_column_index.find(cgc.getCubeIndex(pivot));
				if (pair != pivot_column_index.end()) {
					j = pair->second;
					numRecurse++;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, working_boundary);
						cached_column_idx.push(i);
						if (cached_column_idx.size() > config.cacheSize) {
							cache.erase(cached_column_idx.front());
							cached_column_idx.pop();
						}
                    }
					pivot_column_index[cgc.getCubeIndex(pivot)] = i;

					if (matchMapComp[cgcComp.getCubeIndex(ctr[i])]) {
						matchMapIm.emplace(cgcComp.getCubeIndex(ctr[i]), pivot);
					}
					break;
				}
			} else { break; }
		}
	}
}

Cube InterDimensions::popPivot(CubeQue& column) const {
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
	if (cube.birth <= config.threshold) {
		auto find = pivot_column_index.find(cgc.getCubeIndex(cube));
		if (find != pivot_column_index.end()) {
			ctr.push_back(cube);
		}
	}
	while (cubeEnum.hasNextCube()) {
		cube = cubeEnum.getNextCube();
		if (cube.birth <= config.threshold) {
			auto find = pivot_column_index.find(cgc.getCubeIndex(cube));
			if (find != pivot_column_index.end()) {
				ctr.push_back(cube);
			}
		}
	}
	sort(ctr.begin(), ctr.end(), CubeComparator());
}

Cube InterDimensions::getPivot(CubeQue& column) const {
	Cube result = popPivot(column);
	if (result.coordinates[0] != NONE) {
		column.push(result);
	}
	return result;
}

void InterDimensions::addCache(index_t i, CubeQue& working_boundary) {
	CubeQue cleanWb;
	while (!working_boundary.empty()) {
		auto c = working_boundary.top();
		working_boundary.pop();
		if (!working_boundary.empty() && c == working_boundary.top()) {
			working_boundary.pop();
		} else {
			cleanWb.push(c);
		}
	}
	cache.emplace(i, cleanWb);
}