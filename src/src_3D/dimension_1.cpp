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

void Dimension1::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp, vector<Cube>& ctrImage) {
#ifdef RUNTIME
	cout << endl << "input 0: ";
#endif
	computePairs(ctr0, 0);
	enumerateEdges(ctr0, cgc0);

#ifdef RUNTIME
	cout << endl << "input 1: ";
#endif
	computePairs(ctr1, 1);
	enumerateEdges(ctr1, cgc1);

#ifdef RUNTIME
	cout << endl << "comparison: ";
#endif
	computePairsComp(ctrComp);
#if not defined(USE_APPARENT_PAIRS_COMP)
	ctrImage = ctrComp;
#endif
	enumerateEdges(ctrComp, cgcComp);

#ifdef RUNTIME
	cout << endl << "image 0: ";
#endif
#if defined(USE_APPARENT_PAIRS_COMP) or defined(USE_CLEARING_DIM0)
	computePairsImage(ctrImage, 0);
#else
	computePairsImage(ctrComp, 0);
#endif

#ifdef RUNTIME
	cout << endl << "image 1: ";
#endif
#if defined(USE_APPARENT_PAIRS_COMP) or defined(USE_CLEARING_DIM0)
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
	Cube pivot;
	size_t j;
	bool cacheHit;
#ifdef USE_CACHE
	queue<index_t> cachedColumnIdx;
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	size_t numRecurse;
#endif
#ifdef USE_EMERGENT_PAIRS
	bool checkEmergentPair;
	size_t numEmergentPairs = 0;
#endif
#if defined(USE_APPARENT_PAIRS) or defined(USE_EMERGENT_PAIRS)
	vector<Cube> faces;
	BoundaryEnumerator enumeratorAP(cgc);
	CoboundaryEnumerator coEnumeratorAP(cgc);
#endif

	for (size_t i = 0; i < ctrSize; ++i) {
		CubeQueue workingBoundary;
		j = i;
		cacheHit = false;
#ifdef USE_CACHE
		numRecurse = 0;
#endif
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
#ifdef USE_CACHE
			if (j != i) { cacheHit = tryCache(j, workingBoundary); }
#endif		
			if (cacheHit) { pivot = getPivot(workingBoundary); }
			else {
#ifdef USE_EMERGENT_PAIRS
				if (j == i) {
					if (isEmergentPair(ctr[i], pivot, j, faces, checkEmergentPair, enumerator, enumeratorAP, coEnumeratorAP)) {
						pivotColumnIndex.emplace(pivot.index, i);
						++numEmergentPairs;
						break;
					} else {
						for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_CACHE
						++numRecurse;
#endif
						if (j != i) { continue; }
						else { pivot = getPivot(workingBoundary); }
					}
				} else {
					enumerator.setBoundaryEnumerator(ctr[j]);
					while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
					pivot = getPivot(workingBoundary);
				}
#else			
				enumerator.setBoundaryEnumerator(ctr[j]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
				pivot = getPivot(workingBoundary);
#endif
			}
#ifdef USE_APPARENT_PAIRS
			while (true) {
				faces.clear();
				if (pivotIsApparentPair(pivot, faces, enumeratorAP, coEnumeratorAP)) {
					for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_CACHE
					++numRecurse;
#endif
					pivot = getPivot(workingBoundary);
				} else { break; }
			}
#endif		
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
#ifdef USE_CACHE
					++numRecurse;
#endif
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
	Cube pivot;
	size_t j;
	bool cacheHit;
#ifdef USE_CACHE
	queue<index_t> cachedColumnIdx;
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	size_t numRecurse;
#endif
#ifdef USE_EMERGENT_PAIRS
	bool checkEmergentPair;
	size_t numEmergentPairs = 0;
#endif
#if defined(USE_APPARENT_PAIRS_COMP) or defined(USE_EMERGENT_PAIRS)
	vector<Cube> faces;
	BoundaryEnumerator enumeratorAP(cgcComp);
	CoboundaryEnumerator coEnumeratorAP(cgcComp);
#endif
#if defined (USE_CLEARING_IMAGE) and not defined(USE_APPARENT_PAIRS_COMP)
	bool shouldClear = false;
#endif

	for (size_t i = 0; i < ctrSize; ++i) {
		CubeQueue workingBoundary;
		j = i;
		cacheHit = false;
#ifdef USE_CACHE
		numRecurse = 0;
#endif
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
#ifdef USE_CACHE
			if (j != i) { cacheHit = tryCache(j, workingBoundary); }
#endif
			if (cacheHit) { pivot = getPivot(workingBoundary); } 
			else {
#ifdef USE_EMERGENT_PAIRS
				if (j == i) {
					if (isEmergentPairComp(ctr[i], pivot, j, faces, checkEmergentPair, enumerator, enumeratorAP, coEnumeratorAP)) {
						pivotColumnIndex.emplace(pivot.index, i);
						++numEmergentPairs;
						break;
					} else {
						for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_CACHE
						++numRecurse;
#endif
						if (j != i) { continue; }
						else { pivot = getPivot(workingBoundary); }
					}
				} else {
					enumerator.setBoundaryEnumerator(ctr[j]);
					while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
					pivot = getPivot(workingBoundary);
				}
#else			
				enumerator.setBoundaryEnumerator(ctr[j]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
				pivot = getPivot(workingBoundary);
#endif
			}
#ifdef USE_APPARENT_PAIRS_COMP
			while (true) {
				faces.clear();
				if (pivotIsApparentPair(pivot, faces, enumeratorAP, coEnumeratorAP)) {
					for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_CACHE
					++numRecurse;
#endif
					pivot = getPivot(workingBoundary);
				} else { break; }
			}
#endif
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
#ifdef USE_CACHE
					++numRecurse;
#endif
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
#if defined(USE_CLEARING_IMAGE) and not defined(USE_APPARENT_PAIRS_COMP)
				ctr[i].index = NONE;
				shouldClear = true;
#endif
				break;
			}
		}
	}
#if defined(USE_CLEARING_IMAGE) and not defined(USE_APPARENT_PAIRS_COMP)
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
	Cube pivot;
	value_t birth;
	size_t j;
	bool cacheHit;
#ifdef USE_CACHE
	queue<index_t> cachedColumnIdx;
	cache.clear();
	cache.reserve(min(config.cacheSize, ctrSize));
	size_t numRecurse;
#endif
#ifdef USE_EMERGENT_PAIRS
	bool checkEmergentPair;
	size_t numEmergentPairs = 0;
#endif
#ifdef USE_EMERGENT_PAIRS
	vector<Cube> faces;
#endif
#ifdef USE_CLEARING_IMAGE
	bool shouldClear = false;
#endif

	for (size_t i = 0; i < ctrSize; ++i) {
		CubeQueue workingBoundary;
		j = i;
		cacheHit = false;
#ifdef USE_CACHE
		numRecurse = 0;
#endif
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
#ifdef USE_CACHE
			if (j != i) { cacheHit = tryCache(j, workingBoundary); }
#endif		
			if (cacheHit) { pivot = getPivot(workingBoundary); }
			else {
#ifdef USE_EMERGENT_PAIRS
				if (j == i) {
					if (isEmergentPairImage(ctr[i], pivot, j, faces, checkEmergentPair, cgc, enumerator)) {
						pivotColumnIndex.emplace(pivot.index, i);
						if (isPairedComp[ctr[i].index]) { matchMapIm.emplace(ctr[i].index, pivot.index); }
						++numEmergentPairs;
						break;
					} else {
						for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_CACHE
						++numRecurse;
#endif
						if (j != i) { continue; }
						else { pivot = getPivot(workingBoundary); }
					}
				} else {
					enumerator.setBoundaryEnumerator(ctr[j]);
					while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
					pivot = getPivot(workingBoundary);
				}
#else		
				enumerator.setBoundaryEnumerator(ctr[j]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
				pivot = getPivot(workingBoundary);
#endif
			}
			if (pivot.index != NONE) {
				auto pair = pivotColumnIndex.find(pivot.index);
				if (pair != pivotColumnIndex.end()) {
					j = pair->second;
#ifdef USE_CACHE
					++numRecurse;
#endif
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

void Dimension1::enumerateEdges(vector<Cube>& edges, const CubicalGridComplex& cgc) const {
#ifdef RUNTIME
	cout << ", enumeration ";
	auto start = high_resolution_clock::now();
#endif
	edges.clear();
	edges.reserve(cgc.getNumberOfCubes(1));
	value_t birth;
#ifdef USE_CLEARING_DIM0
	Cube cube;
#endif
	for (index_t x = 0; x < cgc.shape[0]; ++x) {
		for (index_t y = 0; y < cgc.shape[1]; ++y) {
			for (index_t z = 0; z < cgc.shape[2]; ++z) {
				for (uint8_t type = 0; type < 3; ++type) {
					birth = cgc.getBirth(x, y, z, type, 1);
					if (birth < config.threshold) {
#ifdef USE_CLEARING_DIM0
						cube = Cube(birth, x, y, z, type);
						auto find = pivotColumnIndex.find(cube.index);
						if (find == pivotColumnIndex.end()) { edges.push_back(cube); }
#else
						edges.push_back(Cube(birth, x, y, z, type));
#endif
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
	auto pair = cache.find(j);
	if (pair != cache.end()) {
		auto cachedBoundary = pair->second;
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

#ifdef USE_EMERGENT_PAIRS
bool Dimension1::isEmergentPair(const Cube& column, Cube& pivot, size_t& j, vector<Cube>& faces, bool& checkEmergentPair, 
								BoundaryEnumerator& enumerator, BoundaryEnumerator& enumeratorAP, 
								CoboundaryEnumerator& coEnumeratorAP) const {
	faces.clear();
	enumerator.setBoundaryEnumerator(column);
	while (enumerator.hasPreviousFace()) {
		if (checkEmergentPair && enumerator.nextFace.birth == column.birth) {
			auto pair = pivotColumnIndex.find(enumerator.nextFace.index);
#ifdef USE_APPARENT_PAIRS
			if (pair != pivotColumnIndex.end()) {
				checkEmergentPair = false;
				j = pair->second;
			} else if (pivotOfColumnIsApparentPair(enumerator.nextFace, column, faces, enumeratorAP, coEnumeratorAP)) {
				checkEmergentPair = false;
			} else {
				pivot = enumerator.nextFace;
				return true;
			}
#else
			if (pair != pivotColumnIndex.end()) {
				checkEmergentPair = false;
				j = pair->second;
			} else {
				pivot = enumerator.nextFace;
				return true;
			}
#endif
		}
		faces.push_back(enumerator.nextFace);
	}
	return false;
}

bool Dimension1::isEmergentPairComp(const Cube& column, Cube& pivot, size_t& j, vector<Cube>& faces, bool& checkEmergentPair, 
										BoundaryEnumerator& enumerator, BoundaryEnumerator& enumeratorAP, 
										CoboundaryEnumerator& coEnumeratorAP) const {
	faces.clear();
	enumerator.setBoundaryEnumerator(column);
	while (enumerator.hasPreviousFace()) {
		if (checkEmergentPair && enumerator.nextFace.birth == column.birth) {
			auto pair = pivotColumnIndex.find(enumerator.nextFace.index);
#ifdef USE_APPARENT_PAIRS_COMP
			if (pair != pivotColumnIndex.end()) {
				checkEmergentPair = false;
				j = pair->second;
			} else if (pivotOfColumnIsApparentPair(enumerator.nextFace, column, faces, enumeratorAP, coEnumeratorAP)) {
				checkEmergentPair = false;
			} else {
				pivot = enumerator.nextFace;
				return true;
			}
#else
			if (pair != pivotColumnIndex.end()) {
				checkEmergentPair = false;
				j = pair->second;
			} else {
				pivot = enumerator.nextFace;
				return true;
			}
#endif
		}
		faces.push_back(enumerator.nextFace);
	}
	return false;
}

bool Dimension1::isEmergentPairImage(const Cube& column, Cube& pivot, size_t& j, vector<Cube>& faces, bool& checkEmergentPair, 
										const CubicalGridComplex& cgc, BoundaryEnumerator& enumerator) const {
	value_t birth = cgc.getBirth(column.x(), column.y(), column.z(), column.type(), 2);
	faces.clear();
	enumerator.setBoundaryEnumerator(column);
	while (enumerator.hasPreviousFace()) {
		if (checkEmergentPair && enumerator.nextFace.birth == birth) {
			auto pair = pivotColumnIndex.find(enumerator.nextFace.index);
			if (pair != pivotColumnIndex.end()) {
				checkEmergentPair = false;
				j = pair->second;
			} else {
				pivot = enumerator.nextFace;
				return true;
			}
		}
		faces.push_back(enumerator.nextFace);
	}
	return false;
}
#endif

#if defined(USE_APPARENT_PAIRS) or defined(USE_APPARENT_PAIRS_COMP)
bool Dimension1::pivotIsApparentPair(const Cube& pivot, vector<Cube>& faces, 
												BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const {
	bool foundApparentPair = false;
	coEnumerator.setCoboundaryEnumerator(pivot);
	while (coEnumerator.hasNextCoface()) {
		if (coEnumerator.nextCoface.birth == pivot.birth) {
			vector<Cube> facesCopy = faces;
			enumerator.setBoundaryEnumerator(coEnumerator.nextCoface);
			while (enumerator.hasPreviousFace()) {
				if (enumerator.nextFace == pivot) { foundApparentPair = true; } 
				else if (!foundApparentPair && enumerator.nextFace.birth == coEnumerator.nextCoface.birth) { return false; }
				facesCopy.push_back(enumerator.nextFace);
			}
			if (foundApparentPair) { faces = facesCopy; }
			break;
		}
	}
	return foundApparentPair;
}

bool Dimension1::pivotOfColumnIsApparentPair(const Cube& pivot, const Cube& column, vector<Cube>& faces, 
												BoundaryEnumerator& enumerator, CoboundaryEnumerator& coEnumerator) const {
	bool foundApparentPair = false;
	coEnumerator.setCoboundaryEnumerator(pivot);
	while (coEnumerator.hasNextCoface()) {
		if (coEnumerator.nextCoface == column) { break; }
		if (coEnumerator.nextCoface.birth == pivot.birth) {
			vector<Cube> facesCopy = faces;
			enumerator.setBoundaryEnumerator(coEnumerator.nextCoface);
			while (enumerator.hasPreviousFace()) {
				if (enumerator.nextFace == pivot) { foundApparentPair = true; } 
				else if (!foundApparentPair && enumerator.nextFace.birth == coEnumerator.nextCoface.birth) { return false; }
				facesCopy.push_back(enumerator.nextFace);
			}
			if (foundApparentPair) { faces = facesCopy; }
			break;
		}
	}
	return foundApparentPair;
}
#endif