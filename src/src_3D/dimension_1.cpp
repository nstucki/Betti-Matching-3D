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
	computePairs(ctr1, 1);
#ifdef USE_CLEARING_DIM0
	ctr1.clear();
	enumerateEdges(cgc1, ctr1);
#endif

#ifdef RUNTIME
	cout << endl << "comparison: ";
#endif
	computePairsComp(ctrComp);
#ifdef USE_CLEARING_DIM0
	vector<Cube> ctrImage = ctrComp;
	ctrComp.clear();
	enumerateEdges(cgcComp, ctrComp);
#endif

#ifdef RUNTIME
	cout << endl << "image 0: ";
#endif
#ifdef USE_CLEARING_DIM0
	computePairsImage(ctrImage, 0);
#else
	computePairsImage(ctrComp, 0); 
#endif

#ifdef RUNTIME
	cout << endl << "image 1: ";
#endif
#ifdef USE_CLEARING_DIM0
	computePairsImage(ctrImage, 1);
#else
	computePairsImage(ctrComp, 1); 
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
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);
	BoundaryEnumerator enumerator(cgc);
	CoboundaryEnumerator coEnumerator(cgc);
	Cube pivot;
	size_t j;
	size_t numRecurse;
	bool cacheHit;
#ifdef USE_CACHE
	queue<index_t> cachedColumnIdx;
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
#endif
#ifdef USE_APPARENT_PAIRS
	bool foundApparentPair;
#endif
#ifdef USE_EMERGENT_PAIRS
	size_t numEmergentPairs = 0;
	bool checkEmergentPair;
	bool foundEmergentPair;
#endif
#if defined(USE_APPARENT_PAIRS) || defined(USE_EMERGENT_PAIRS)
	vector<Cube> faces;
#endif
	for (size_t i = 0; i < ctrSize; ++i) {
		CubeQueue workingBoundary;
		j = i;
		cacheHit = false;
		numRecurse = 0;
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
#ifdef USE_CACHE
			if (i != j) { cacheHit = tryCache(j, workingBoundary); }
#endif
			if (!cacheHit) {
#ifdef USE_EMERGENT_PAIRS
				foundEmergentPair = isEmergentPair(ctr[j], pivot, faces, checkEmergentPair, enumerator, coEnumerator);
				if (foundEmergentPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
                    ++numEmergentPairs;
                    break;
                } else { for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); } }
#else			
				enumerator.setBoundaryEnumerator(ctr[j]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			}
#ifdef USE_APPARENT_PAIRS
			while (true) {
				pivot = getPivot(workingBoundary);
				foundApparentPair = isApparentPair(pivot, faces, enumerator, coEnumerator);	
				if (foundApparentPair) {
					for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
				} else { break; }
			}
#else
			pivot = getPivot(workingBoundary);
#endif
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
					++numRecurse;
					continue;
				} else {
					pivotColumnIndex.emplace(pivot.index, i);
					if (pivot.birth != ctr[i].birth) {
						pairs.push_back(Pair(pivot, ctr[i]));
						matchMap.emplace(pivot.index, pairs.back());
					}
#ifdef USE_CACHE
					if (numRecurse >= config.minRecursionToCache) { addCache(i, workingBoundary, cachedColumnIdx); }
#endif
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
	cout << " with " << numEmergentPairs << " emergent pairs";
#endif
#endif
}

void Dimension1::computePairsComp(vector<Cube>& ctr) {
#ifdef RUNTIME
	cout << "barcode ";
	auto start = high_resolution_clock::now();
#endif
	size_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);	
	BoundaryEnumerator enumerator(cgcComp);
	CoboundaryEnumerator coEnumerator(cgcComp);
	Cube pivot;
	size_t j;
	size_t numRecurse;
	bool cacheHit;
#ifdef USE_CACHE
	queue<index_t> cachedColumnIdx;
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
#endif
#ifdef USE_APPARENT_PAIRS
	bool foundApparentPair;
#endif
#ifdef USE_EMERGENT_PAIRS
	size_t numEmergentPairs = 0;
	bool checkEmergentPair;
	bool foundEmergentPair;
#endif
#if defined(USE_APPARENT_PAIRS) || defined(USE_EMERGENT_PAIRS)
	vector<Cube> faces;
#endif
#ifdef USE_CLEARING_IMAGE
	bool shouldClear = false;
#endif
	for (size_t i = 0; i < ctrSize; ++i) {
		CubeQueue workingBoundary;
		j = i;
		cacheHit = false;
		numRecurse = 0;
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
#ifdef USE_CACHE
			if (i != j) { cacheHit = tryCache(j, workingBoundary); }
#endif
			if (!cacheHit) {
#ifdef USE_EMERGENT_PAIRS
				foundEmergentPair = isEmergentPair(ctr[j], pivot, faces, checkEmergentPair, enumerator, coEnumerator);
				if (foundEmergentPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
                    ++numEmergentPairs;
                    break;
                } else { for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); } }					
#else			
				enumerator.setBoundaryEnumerator(ctr[j]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			}
#ifdef USE_APPARENT_PAIRS
			while (true) {
				pivot = getPivot(workingBoundary);
				foundApparentPair = isApparentPair(pivot, faces, enumerator, coEnumerator);	
				if (foundApparentPair) {
					for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
				} else { break; }
			}
#else
			pivot = getPivot(workingBoundary);
#endif
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
					++numRecurse;
					continue;
				} else {
					pivotColumnIndex.emplace(pivot.index, i);
					if (pivot.birth != ctr[i].birth) {
						pairsComp.push_back(Pair(pivot, ctr[i]));
						isPairedComp.emplace(ctr[i].index, true);
					}
#ifdef USE_CACHE
					if (numRecurse >= config.minRecursionToCache) { addCache(i, workingBoundary, cachedColumnIdx); }
#endif
					break;
				}
			} else {
#ifdef USE_CLEARING_IMAGE
				ctr[i].index = NONE;
				shouldClear = true;
#endif
				break;
			}
		}
	}
#ifdef USE_CLEARING_IMAGE
	if (shouldClear) {
		auto newEnd = remove_if(ctr.begin(), ctr.end(), [](const Cube& cube) { return cube.index == NONE; });
		ctr.erase(newEnd, ctr.end());
	}
#endif
#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#ifdef USE_EMERGENT_PAIRS
	cout << " with " << numEmergentPairs << " emergent pairs";
#endif
#endif
}

void Dimension1::computePairsImage(vector<Cube>& ctr, uint8_t k) {
#ifdef RUNTIME
	cout << "barcode ";
	auto start = high_resolution_clock::now();
#endif
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	unordered_map<uint64_t, uint64_t>& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;
	matchMapIm.reserve(pairsComp.size());
	size_t ctrSize = ctr.size();
	pivotColumnIndex.clear();
	pivotColumnIndex.reserve(ctrSize);
	BoundaryEnumerator enumerator = BoundaryEnumerator(cgc);
	CoboundaryEnumerator coEnumerator(cgcComp);
	Cube pivot;
	value_t birth;
	size_t j;
	size_t numRecurse;
	bool cacheHit;
#ifdef USE_CACHE
	queue<index_t> cachedColumnIdx;
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
#endif
#ifdef USE_APPARENT_PAIRS
	bool foundApparentPair;
#endif
#ifdef USE_EMERGENT_PAIRS
	size_t numEmergentPairs = 0;
	bool checkEmergentPair;
	bool foundEmergentPair;
#endif
#if defined(USE_APPARENT_PAIRS) || defined(USE_EMERGENT_PAIRS)
	vector<Cube> faces;
#endif
#ifdef USE_CLEARING_IMAGE
	bool shouldClear = false;
#endif
	for (size_t i = 0; i < ctrSize; ++i) {
		CubeQueue workingBoundary;
		j = i;
		numRecurse = 0;
		cacheHit = false;
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
#ifdef USE_CACHE
			if (i != j) { cacheHit = tryCache(j, workingBoundary); }
#endif
			if (!cacheHit) {
#ifdef USE_EMERGENT_PAIRS
				foundEmergentPair = isEmergentPairImage(ctr[j], pivot, faces, checkEmergentPair, enumerator, coEnumerator, cgc);
				if (foundEmergentPair) {
                    pivotColumnIndex.emplace(pivot.index, i);
					if (isPairedComp[ctr[i].index]) { matchMapIm.emplace(ctr[i].index, pivot.index); }
                    ++numEmergentPairs;
                    break;
                } else { for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); } }
#else			
				enumerator.setBoundaryEnumerator(ctr[j]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			}
#ifdef USE_APPARENT_PAIRS
			while (true) {
				pivot = getPivot(workingBoundary);
				foundApparentPair = isApparentPairImage(pivot, faces, enumerator, coEnumerator);	
				if (foundApparentPair) {
					for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
				} else { break; }
			}
#else
			pivot = getPivot(workingBoundary);
#endif
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
					++numRecurse;
					continue;
				} else {
					pivotColumnIndex.emplace(pivot.index, i);
					if (isPairedComp[ctr[i].index]) { matchMapIm.emplace(ctr[i].index, pivot.index); }
#ifdef USE_CACHE
					if (numRecurse >= config.minRecursionToCache) { addCache(i, workingBoundary, cachedColumnIdx); }
#endif
					break;
				}
			} else {
#ifdef USE_CLEARING_IMAGE
				ctr[i].index = NONE;
				shouldClear = true;
#endif
				break;
			}
		}
	}
#ifdef USE_CLEARING_IMAGE
	if (shouldClear) {
		auto newEnd = remove_if(ctr.begin(), ctr.end(), [](const Cube& cube) { return cube.index == NONE; });
		ctr.erase(newEnd, ctr.end());
	}
#endif
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

#ifdef USE_CACHE
bool Dimension1::tryCache(const size_t& j, CubeQueue& workingBoundary) const {
	auto findCb = cache.find(j);
	if (findCb != cache.end()) {
		auto cachedBoundary = findCb->second;
		while (!cachedBoundary.empty()) {
			workingBoundary.push(cachedBoundary.top());
			cachedBoundary.pop();
		}
		return true;
	} else { return false; }
}

void Dimension1::addCache(const index_t& i, CubeQueue& workingBoundary, queue<index_t>& cachedColumnIdx) {
	CubeQueue cleanWb;
	Cube c;
	while (!workingBoundary.empty()) {
		c = workingBoundary.top();
		workingBoundary.pop();
		if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
		else { cleanWb.push(c); }
	}
	cache.emplace(i, cleanWb);
	cachedColumnIdx.push(i);
	if (cachedColumnIdx.size() > config.cacheSize) {
		cache.erase(cachedColumnIdx.front());
		cachedColumnIdx.pop();
	}
}
#endif

#ifdef USE_APPARENT_PAIRS
bool Dimension1::isApparentPair(const Cube& face, vector<Cube>& faces, 
									BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const {
	bool foundApparentPair = false;
	coEnumerator.setCoboundaryEnumerator(face);
	while (coEnumerator.hasNextCoface()) {
		if (coEnumerator.nextCoface.birth == face.birth) {
			faces.clear();
			enumerator.setBoundaryEnumerator(coEnumerator.nextCoface);
			while (enumerator.hasPreviousFace()) {
				if (enumerator.nextFace == face) { foundApparentPair = true; } 
				else if (!foundApparentPair && enumerator.nextFace.birth == coEnumerator.nextCoface.birth) { return false; }
				faces.push_back(enumerator.nextFace);
			}
			break;
		}
	}
	return foundApparentPair;
}

bool Dimension1::isApparentPairImage(const Cube& face, vector<Cube>& faces, 
									BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const {
	bool foundApparentPair = false;
	value_t birth = cgcComp.getBirth(face.x(), face.y(), face.z(), face.type(), 1);
	coEnumerator.setCoboundaryEnumerator(face);
	while (coEnumerator.hasNextCoface()) {
		if (coEnumerator.nextCoface.birth == birth) {
			faces.clear();
			enumerator.setBoundaryEnumerator(coEnumerator.nextCoface);
			while (enumerator.hasPreviousFace()) {
				if (enumerator.nextFace == face) { foundApparentPair = true; } 
				else {
					birth = cgcComp.getBirth(enumerator.nextFace.x(), enumerator.nextFace.y(), enumerator.nextFace.z(), enumerator.nextFace.type(), 1);
					if (!foundApparentPair && birth == coEnumerator.nextCoface.birth) { return false; }
				}
				faces.push_back(enumerator.nextFace);
			}
			break;
		}
	}
	return foundApparentPair;
}

bool Dimension1::pivotIsApparentPair(const Cube& pivot, const Cube& column, CoboundaryEnumerator& coEnumerator) const {
	coEnumerator.setCoboundaryEnumerator(pivot);
	while (coEnumerator.hasNextCoface()) {
		if (coEnumerator.nextCoface.birth == pivot.birth) {
			if (coEnumerator.nextCoface == column) { return false; }
			else { return true; }
		}
	}
	return false;
}

bool Dimension1::pivotIsApparentPairImage(const Cube& pivot, const Cube& column, CoboundaryEnumerator& coEnumerator) const {
	value_t birth = cgcComp.getBirth(pivot.x(), pivot.y(), pivot.z(), pivot.type(), 1);
	coEnumerator.setCoboundaryEnumerator(pivot);
	while (coEnumerator.hasNextCoface()) {
		if (coEnumerator.nextCoface.birth == birth) {
			if (coEnumerator.nextCoface == column) { return false; }
			else { return true; }
		}
	}
	return false;
}
#endif

#ifdef USE_EMERGENT_PAIRS
bool Dimension1::isEmergentPair(const Cube& column, Cube& pivot, vector<Cube>& faces, bool& checkEmergentPair, BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const {
	faces.clear();
	enumerator.setBoundaryEnumerator(column);
	while (enumerator.hasPreviousFace()) {
		if (checkEmergentPair && enumerator.nextFace.birth == column.birth) {
#ifdef USE_APPARENT_PAIRS
			if (!pivotIsApparentPair(enumerator.nextFace, column, coEnumerator) && 
					pivotColumnIndex.find(enumerator.nextFace.index) == pivotColumnIndex.end()) {
				pivot = enumerator.nextFace;
                return true;
			} else { checkEmergentPair = false; }
#else
			if (pivotColumnIndex.find(enumerator.nextFace.index) == pivotColumnIndex.end()) {
				pivot = enumerator.nextFace;
                return true;
			} else { checkEmergentPair = false; }
#endif
		}
		faces.push_back(enumerator.nextFace);
	}
	return false;
}

bool Dimension1::isEmergentPairImage(const Cube& column, Cube& pivot, vector<Cube>& faces, bool& checkEmergentPair, BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator, const CubicalGridComplex& cgc) const {
	faces.clear();
	enumerator.setBoundaryEnumerator(column);
	while (enumerator.hasPreviousFace()) {
		if (checkEmergentPair) {
			value_t birth = cgc.getBirth(column.x(), column.y(), column.z(), column.type(), 2);
			if (enumerator.nextFace.birth == birth) {
#ifdef USE_APPARENT_PAIRS
				if (!pivotIsApparentPairImage(enumerator.nextFace, column, coEnumerator) && 
						pivotColumnIndex.find(enumerator.nextFace.index) == pivotColumnIndex.end()) {
					pivot = enumerator.nextFace;
					return true;
				} else { checkEmergentPair = false; }
#else
				if (pivotColumnIndex.find(enumerator.nextFace.index) == pivotColumnIndex.end()) {
					pivot = enumerator.nextFace;
					return true;
				} else { checkEmergentPair = false; }
#endif
			}
		}
		faces.push_back(enumerator.nextFace);
	}
	return false;
}
#endif