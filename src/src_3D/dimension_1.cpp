#include "dimension_1.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>
#include <algorithm>

using namespace dim3;
using namespace std::chrono;


Dimension1::Dimension1(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, 
						const CubicalGridComplex& _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), pairs0(_pairs0), pairs1(_pairs1),
						pairsComp(_pairsComp), matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1) {}

void Dimension1::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
#ifdef RUNTIME
	cout << endl << "input 0: ";
#endif
	computePairs(ctr0, 0);
#ifdef USE_CLEARING_DIM0
	ctr0.clear();
	enumerateEdges(cgc0, ctr0);
#endif

#ifdef RUNTIME
	cout << endl << "input 1: ";
#endif
	pivotColumnIndex.clear();
	cache.clear();
	computePairs(ctr1, 1);
#ifdef USE_CLEARING_DIM0
	ctr1.clear();
	enumerateEdges(cgc1, ctr1);
#endif

#ifdef RUNTIME
	cout << endl << "comparison: ";
#endif
	pivotColumnIndex.clear();
	cache.clear();
	computePairsComp(ctrComp);
#ifdef USE_CLEARING_DIM0
	vector<Cube> ctrImage = ctrComp;
	ctrComp.clear();
	enumerateEdges(cgcComp, ctrComp);
#endif

#ifdef RUNTIME
	cout << endl << "image 0: ";
#endif
	pivotColumnIndex.clear();
	cache.clear();
#ifdef USE_CLEARING_DIM0
	computeImagePairs(ctrImage, 0);
#else
	computeImagePairs(ctrComp, 0); 
#endif

#ifdef RUNTIME
	cout << endl << "image 1: ";
#endif
	pivotColumnIndex.clear();
	cache.clear();
#ifdef USE_CLEARING_DIM0
	computeImagePairs(ctrImage, 1);
#else
	computeImagePairs(ctrComp, 1); 
#endif
	
#ifdef RUNTIME
	cout << endl << "matching: ";
#endif
	computeMatching();
}

void Dimension1::computePairs(const vector<Cube>& ctr, uint8_t k) {
#ifdef RUNTIME
	cout << "barcode ";
	auto start = high_resolution_clock::now();
#endif
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	vector<Pair>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<uint64_t, Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;

	size_t ctrSize = ctr.size();
	pivotColumnIndex.reserve(ctrSize);
	cache.reserve(min(config.cacheSize, ctrSize));
	BoundaryEnumerator enumerator(cgc);
	Cube pivot;
	queue<index_t> cachedColumnIdx;
	size_t numRecurse;
	size_t j;
	bool cacheHit;
#ifdef USE_EMERGENT_PAIRS
	vector<Cube> faces;
	size_t numEmergentPairs = 0;
	bool checkEmergentPair;
	bool foundPair;
#endif
	for (size_t i = 0; i < ctrSize; ++i) {
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
		foundPair = false;
#endif
		CubeQueue workingBoundary;
		j = i;
		numRecurse = 0;
		while (true) {
			cacheHit = false;
			if (i != j) {
                auto findCb = cache.find(j);
                if (findCb != cache.end()) {
                    cacheHit = true;
                    auto cachedBoundary = findCb->second;
                    while (!cachedBoundary.empty()) {
                        workingBoundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				enumerator.setBoundaryEnumerator(ctr[j]);
#ifdef USE_EMERGENT_PAIRS
				faces.clear();
				while (enumerator.hasPreviousFace()) {
					if (checkEmergentPair && ctr[j].birth == enumerator.nextFace.birth) {
						if (pivotColumnIndex.find(enumerator.nextFace.index) == pivotColumnIndex.end()) {
							pivot = enumerator.nextFace;
                            foundPair = true;
							break;
						} else { checkEmergentPair = false; }
					}
					faces.push_back(enumerator.nextFace);
				}
				if (foundPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
                    ++numEmergentPairs;
                    break;
                } else { for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); } }
#else			
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			}
			pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
					++numRecurse;
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
					pivotColumnIndex.emplace(pivot.index, i);
					if (pivot.birth != ctr[i].birth) {
						pairs.push_back(Pair(pivot, ctr[i]));
						matchMap.emplace(pivot.index, pairs.back());
					}
					break;
				}
			} else { break; }
		}
	}
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#ifdef USE_EMERGENT_PAIRS
	cout <<" with " << numEmergentPairs << " emergent pairs";
#endif
#endif
}

void Dimension1::computePairsComp(vector<Cube>& ctr) {
#ifdef RUNTIME
	cout << "barcode ";
	auto start = high_resolution_clock::now();
#endif
	size_t ctrSize = ctr.size();
	pivotColumnIndex.reserve(ctrSize);	
	cache.reserve(min(config.cacheSize, ctrSize));
	BoundaryEnumerator enumerator = BoundaryEnumerator(cgcComp);
	Cube pivot;
	queue<index_t> cachedColumnIdx;
	size_t numRecurse;
	size_t j;
	bool cacheHit;
	bool shouldClear = false;
#ifdef USE_EMERGENT_PAIRS
	vector<Cube> faces;
	size_t numEmergentPairs = 0;
	bool checkEmergentPair;
	bool foundPair;
#endif
	for (size_t i = 0; i < ctrSize; ++i) {
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
		foundPair = false;
#endif
		CubeQueue workingBoundary;
		j = i;
		numRecurse = 0;
		while (true) {
			cacheHit = false;
			if (i != j) {
                auto findCb = cache.find(j);
                if (findCb != cache.end()) {
                    cacheHit = true;
                    auto cachedBoundary = findCb->second;
                    while (!cachedBoundary.empty()) {
                        workingBoundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				enumerator.setBoundaryEnumerator(ctr[j]);
#ifdef USE_EMERGENT_PAIRS
				faces.clear();
				while (enumerator.hasPreviousFace()) {
					if (checkEmergentPair && ctr[j].birth == enumerator.nextFace.birth) {
						if (pivotColumnIndex.find(enumerator.nextFace.index) == pivotColumnIndex.end()) {
							pivot = enumerator.nextFace;
                            foundPair = true;
							break;
						} else { checkEmergentPair = false; }
					}
					faces.push_back(enumerator.nextFace);
				}
				if (foundPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
                    ++numEmergentPairs;
                    break;
                } else { for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); } }
#else			
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			}
			pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
					++numRecurse;
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
					pivotColumnIndex.emplace(pivot.index, i);
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
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#ifdef USE_EMERGENT_PAIRS
	cout << " with " << numEmergentPairs << " emergent pairs";
#endif
#endif
}

void Dimension1::computeImagePairs(vector<Cube>& ctr, uint8_t k) {
#ifdef RUNTIME
	cout << "barcode ";
	auto start = high_resolution_clock::now();
#endif
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	unordered_map<uint64_t, uint64_t>& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;

	size_t ctrSize = ctr.size();
	pivotColumnIndex.reserve(ctrSize);
	cache.reserve(min(config.cacheSize, ctrSize));
	matchMapIm.reserve(pairsComp.size());
	BoundaryEnumerator enumerator = BoundaryEnumerator(cgc);
	Cube pivot;
	value_t birth;
	queue<index_t> cachedColumnIdx;
	size_t numRecurse;
	size_t j;
	bool cacheHit;
	bool shouldClear = false;
#ifdef USE_EMERGENT_PAIRS
	vector<Cube> faces;
	size_t numEmergentPairs = 0;
	bool checkEmergentPair;
	bool foundPair;
#endif
	for (size_t i = 0; i < ctrSize; ++i) {
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
		foundPair = false;
#endif
		CubeQueue workingBoundary;
		j = i;
		numRecurse = 0;
		while (true) {
			cacheHit = false;
			if (i != j) {
                auto findCb = cache.find(j);
                if (findCb != cache.end()) {
                    cacheHit = true;
                    auto cachedBoundary = findCb->second;
                    while (!cachedBoundary.empty()) {
                        workingBoundary.push(cachedBoundary.top());
                        cachedBoundary.pop();
                    }
                }
			}
			if (!cacheHit) {
				enumerator.setBoundaryEnumerator(ctr[j]);
#ifdef USE_EMERGENT_PAIRS
				faces.clear();
				while (enumerator.hasPreviousFace()) {
					if (checkEmergentPair) {
						birth = cgc.getBirth(ctr[j].x(), ctr[j].y(), ctr[j].z(), ctr[j].type(), 2);
						if (birth == enumerator.nextFace.birth) {
							if (pivotColumnIndex.find(enumerator.nextFace.index) == pivotColumnIndex.end()) {
								pivot = enumerator.nextFace;
                            	foundPair = true;
								break;
							} else { checkEmergentPair = false; }
						}
					}
					faces.push_back(enumerator.nextFace);
				}
				if (foundPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
					if (isMatchedComp[ctr[i].index]) { matchMapIm.emplace(ctr[i].index, pivot.index); }
                    ++numEmergentPairs;
                    break;
                } else { for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); } }
#else			
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			}
			pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
					++numRecurse;
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
					pivotColumnIndex.emplace(pivot.index, i);
					if (isMatchedComp[ctr[i].index]) { matchMapIm.emplace(ctr[i].index, pivot.index); }
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
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#ifdef USE_EMERGENT_PAIRS
	cout << " with " << numEmergentPairs << " emergent pairs";
#endif
#endif
}

void Dimension1::computeMatching() {
#ifdef RUNTIME
	auto start = high_resolution_clock::now();
#endif
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
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}

void Dimension1::enumerateEdges(const CubicalGridComplex& cgc, vector<Cube>& edges) const {
#ifdef RUNTIME
	cout << ", enumeration ";
	auto start = high_resolution_clock::now();
#endif
	edges.reserve(cgc.getNumberOfCubes(1));
	value_t birth;
	Cube cube;
	for (index_t x = 0; x < cgc.shape[0]; ++x) {
		for (index_t y = 0; y < cgc.shape[1]; ++y) {
			for (index_t z = 0; z < cgc.shape[2]; ++z) {
				for (uint8_t type = 0; type < 3; ++type) {
					birth = cgc.getBirth(x, y, z, type, 1);
					if (birth < config.threshold) {
						cube = Cube(birth, x, y, z, type);
						auto find = pivotColumnIndex.find(cube.index);
						if (find == pivotColumnIndex.end()) { edges.push_back(cube); }
					}	
				}				
			}
		}
	}
	sort(edges.begin(), edges.end(), CubeComparator());
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
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

void Dimension1::addCache(const index_t& i, CubeQueue& workingBoundary) {
	CubeQueue cleanWb;
	while (!workingBoundary.empty()) {
		Cube c = workingBoundary.top();
		workingBoundary.pop();
		if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
		else { cleanWb.push(c); }
	}
	cache.emplace(i, cleanWb);
}