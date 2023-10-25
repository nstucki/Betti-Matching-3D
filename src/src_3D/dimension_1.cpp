#include "dimension_1.h"
#include "data_structures.h"
#include "enumerators.h"

#include <iostream>
#include <chrono>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <future>

using namespace dim3;
using namespace std::chrono;


Dimension1::Dimension1(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, 
						const CubicalGridComplex& _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches,
						Cube1Map<bool>& _isMatched0, Cube1Map<bool>& _isMatched1
						) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), pairs0(_pairs0), pairs1(_pairs1),
						pairsComp(_pairsComp), matches(_matches),
						isMatched0(_isMatched0), isMatched1(_isMatched1),
						isMatchedComp(_cgc0.shape),
						matchMap0(_cgc0.shape),
						matchMap1(_cgc0.shape),
						matchMapIm0(_cgc0.shape),
						matchMapIm1(_cgc0.shape),
						pivotColumnIndex0(_cgc0.shape),
						pivotColumnIndex1(_cgc0.shape),
						pivotColumnIndexComp(_cgcComp.shape)
						// isMatched0(_isMatched0), isMatched1(_isMatched1), isMatchedComp(cgc0.shape[0] * cgc0.shape[1] * cgc0.shape[2] * 3, false)
						{}

void Dimension1::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {


	auto input0Future = std::async([this, &ctr0]() {
#ifdef RUNTIME
		cout << endl << "input 0: ";
#endif
		computePairs(ctr0, 0);
#ifdef USE_CLEARING_DIM0
		ctr0.clear();
		enumerateEdges(cgc0, ctr0, pivotColumnIndex0);
#endif
	});

	auto input1Future = std::async([this, &ctr1]() {
#ifdef RUNTIME
		cout << endl << "input 1: ";
#endif
		computePairs(ctr1, 1);
#ifdef USE_CLEARING_DIM0
		ctr1.clear();
		enumerateEdges(cgc1, ctr1, pivotColumnIndex1);
#endif
	});

	auto inputCompFuture = std::async([this, &ctrComp]() {
#ifdef RUNTIME
		cout << endl << "comparison: ";
#endif
		computePairsComp(ctrComp);
#ifdef USE_CLEARING_DIM0
		vector<Cube> ctrImage = ctrComp;
		ctrComp.clear();
		enumerateEdges(cgcComp, ctrComp, pivotColumnIndexComp);
		return ctrImage;
#else
	return std::optional<vector<Cube>>();
#endif

	});

	input0Future.wait();
	input1Future.wait();
	auto ctrImage = inputCompFuture.get();

	auto input0AndComparisonFuture = std::async([this, &ctrImage, &ctrComp](){
#ifdef RUNTIME
		cout << endl << "image 0: ";
#endif
		pivotColumnIndex0.clear();
#ifdef USE_CLEARING_DIM0
		computeImagePairs(ctrImage, 0);
#else
		computeImagePairs(ctrComp, 0);
#endif
	});

	auto input1AndComparisonFuture = std::async([this, &ctrImage, &ctrComp](){
#ifdef RUNTIME
		cout << endl << "image 1: ";
#endif
		pivotColumnIndex1.clear();
#ifdef USE_CLEARING_DIM0
		computeImagePairs(ctrImage, 1);
#else
		computeImagePairs(ctrComp, 1);
#endif
	});
	
	input0AndComparisonFuture.wait();
	input1AndComparisonFuture.wait();

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
	Cube1Map<Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;
	auto &pivotColumnIndex = (k == 0) ? pivotColumnIndex0 : pivotColumnIndex1;

	size_t ctrSize = ctr.size();
	// pivotColumnIndex.reserve(ctrSize);
	// cache.reserve(min(config.cacheSize, ctrSize));
	vector<std::optional<vector<Cube>>> cache = vector<std::optional<vector<Cube>>>(ctrSize);
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
                auto &cachedBoundary = cache[j];
                if (cachedBoundary.has_value()) {
                    cacheHit = true;
					for (auto &face : *cachedBoundary) {
						workingBoundary.push(face);
					}
                }
			}
			if (!cacheHit) {
				enumerator.setBoundaryEnumerator(ctr[j]);
#ifdef USE_EMERGENT_PAIRS
				faces.clear();
#endif
				while (enumerator.hasNextFace()) {
#ifdef USE_EMERGENT_PAIRS
					if (checkEmergentPair && ctr[j].birth == enumerator.nextFace.birth) {
						if (!pivotColumnIndex.find(enumerator.nextFace.index).has_value()) {
							pivot = enumerator.nextFace;
                            foundPair = true;
							break;
						} else { checkEmergentPair = false; }
					}
					faces.push_back(enumerator.nextFace);
#else
					workingBoundary.push(enumerator.nextFace);
#endif
				}
#ifdef USE_EMERGENT_PAIRS
				if (foundPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
                    ++numEmergentPairs;
                    break;
                }
				for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#endif
			}
			pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair.has_value()) {
					j = *pair;
					++numRecurse;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, workingBoundary, cache);
						cachedColumnIdx.push(i);
						// if (cachedColumnIdx.size() > config.cacheSize) {
						// 	cache.erase(cachedColumnIdx.front());
						// 	cachedColumnIdx.pop();
						// }
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
	// pivotColumnIndex.reserve(ctrSize);	
	// cache.reserve(min(config.cacheSize, ctrSize));
	vector<std::optional<vector<Cube>>> cache = vector<std::optional<vector<Cube>>>(ctrSize);
	BoundaryEnumerator enumerator = BoundaryEnumerator(cgcComp);
	auto &pivotColumnIndex = pivotColumnIndexComp;
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
                auto &cachedBoundary = cache[j];
                if (cachedBoundary.has_value()) {
                    cacheHit = true;
					for (auto &face : *cachedBoundary) {
						workingBoundary.push(face);
					}
                }
			}
			if (!cacheHit) {
#ifdef USE_EMERGENT_PAIRS
				faces.clear();
#endif
				enumerator.setBoundaryEnumerator(ctr[j]);
				while (enumerator.hasNextFace()) {
#ifdef USE_EMERGENT_PAIRS
					if (checkEmergentPair && ctr[j].birth == enumerator.nextFace.birth) {
						if (!pivotColumnIndex.find(enumerator.nextFace.index).has_value()) {
							pivot = enumerator.nextFace;
                            foundPair = true;
							break;
						} else { checkEmergentPair = false; }
					}
					faces.push_back(enumerator.nextFace);
#else
					workingBoundary.push(enumerator.nextFace);
#endif
				}
#ifdef USE_EMERGENT_PAIRS
				if (foundPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
                    ++numEmergentPairs;
                    break;
                }
				for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#endif
			}
			pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair.has_value()) {
					j = *pair;
					++numRecurse;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, workingBoundary, cache);
						cachedColumnIdx.push(i);
						// if (cachedColumnIdx.size() > config.cacheSize) {
						// 	cache.erase(cachedColumnIdx.front());
						// 	cachedColumnIdx.pop();
						// }
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
	auto& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;
	auto &pivotColumnIndex = (k == 0) ? pivotColumnIndex0 : pivotColumnIndex1;

	size_t ctrSize = ctr.size();
	// pivotColumnIndex.reserve(ctrSize);
	// cache.reserve(min(config.cacheSize, ctrSize));
	// matchMapIm.reserve(pairsComp.size());
	vector<std::optional<vector<Cube>>> cache = vector<std::optional<vector<Cube>>>(ctrSize);
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
                auto &cachedBoundary = cache[j];
                if (cachedBoundary.has_value()) {
                    cacheHit = true;
					for (auto &face : *cachedBoundary) {
						workingBoundary.push(face);
					}
                }
			}
			if (!cacheHit) {
				enumerator.setBoundaryEnumerator(ctr[j]);
#ifdef USE_EMERGENT_PAIRS
				faces.clear();
#endif
				while (enumerator.hasNextFace()) {
#ifdef USE_EMERGENT_PAIRS
					if (checkEmergentPair) {
						birth = cgc.getBirth(ctr[j].x(), ctr[j].y(), ctr[j].z(), ctr[j].type(), 2);
						if (birth == enumerator.nextFace.birth) {
							if (!pivotColumnIndex.find(enumerator.nextFace.index).has_value()) {
								pivot = enumerator.nextFace;
                            	foundPair = true;
								break;
							} else { checkEmergentPair = false; }
						}
					}
					faces.push_back(enumerator.nextFace);
#else
					workingBoundary.push(enumerator.nextFace);
#endif
				}
#ifdef USE_EMERGENT_PAIRS
				if (foundPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
					if (isMatchedComp.find(ctr[i].index).value_or(false)) {
						matchMapIm.emplace(ctr[i].index, pivot.index);
					}
                    ++numEmergentPairs;
                    break;
                }
				for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#endif
			}
			pivot = getPivot(workingBoundary);
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair.has_value()) {
					j = *pair;
					++numRecurse;
					continue;
				} else {
					if (numRecurse >= config.minRecursionToCache) {
                        addCache(i, workingBoundary, cache);
						cachedColumnIdx.push(i);
						// if (cachedColumnIdx.size() > config.cacheSize) {
						// 	cache.erase(cachedColumnIdx.front());
						// 	cachedColumnIdx.pop();
						// }
                    }
					pivotColumnIndex.emplace(pivot.index, i);
					if (isMatchedComp.find(ctr[i].index).value_or(false)) {
						matchMapIm.emplace(ctr[i].index, pivot.index);
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

void Dimension1::computeMatching() {
#ifdef RUNTIME
	auto start = high_resolution_clock::now();
#endif
	uint64_t birthIndex0;
	uint64_t birthIndex1;
	for (Pair& pair : pairsComp) {
		auto find0 = matchMapIm0.find(pair.death.index);
		auto find1 = matchMapIm1.find(pair.death.index);
		if (find0.has_value() && find1.has_value()) {
			birthIndex0 = *find0;
			birthIndex1 = *find1;
			auto find0 = matchMap0.find(birthIndex0);
			auto find1 = matchMap1.find(birthIndex1);
			if (find0.has_value() && find1.has_value()) {
				matches.push_back(Match(*find0, *find1));
				isMatched0.emplace(find0->birth.index, true);
				isMatched1.emplace(find1->birth.index, true);
			}
		}
	}
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}

void Dimension1::enumerateEdges(const CubicalGridComplex& cgc, vector<Cube>& edges, Cube1Map<uint64_t> &pivotColumnIndex) const {
#ifdef RUNTIME
	cout << ", enumeration ";
	auto start = high_resolution_clock::now();
#endif
	edges.reserve(cgc.getNumberOfCubes(1));
	value_t birth;
	Cube cube;
	bool binaryInputs = true;
	for (index_t x = 0; x < cgc.shape[0]; ++x) {
		for (index_t y = 0; y < cgc.shape[1]; ++y) {
			for (index_t z = 0; z < cgc.shape[2]; ++z) {
				for (uint8_t type = 0; type < 3; ++type) {
					birth = cgc.getBirth(x, y, z, type, 1);
					if (birth < config.threshold) {
						cube = Cube(birth, x, y, z, type);
						auto find = pivotColumnIndex.find(cube.index);
						if (!find.has_value()) {
							edges.push_back(cube);
							if (binaryInputs && cube.birth != 0 && cube.birth != 1) binaryInputs = false;
						}
					}	
				}				
			}
		}
	}
	if (binaryInputs) {
		std::stable_partition(edges.begin(), edges.end(), [](Cube &cube) { return cube.birth == 0; });
	} else {
		std::stable_sort(edges.begin(), edges.end(), [](const Cube &cube1, const Cube &cube2) { return cube1.birth < cube2.birth; }); //CubeComparator
	}
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms" << endl;
	cout << "Binary inputs: " << (binaryInputs ? "true" : "false") << endl;
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

void Dimension1::addCache(const index_t& i, CubeQueue& workingBoundary, vector<std::optional<vector<Cube>>> &cache) {
	std::vector<Cube> cleanWb;
	while (!workingBoundary.empty()) {
		Cube c = workingBoundary.top();
		workingBoundary.pop();
		if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
		else { cleanWb.emplace_back(c); }
	}
	cache[i] = std::move(cleanWb);
}
