#include "dimension_1.h"
#include "enumerators.h"

#include <algorithm>
#include <chrono>
#include <future>
#include <iostream>
#include <stdexcept>

using namespace dim3;
using namespace std::chrono;


Dimension1::Dimension1(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, 
						const CubicalGridComplex& _cgcComp, const Config& _config, 
						vector<Pair>& _pairs0, vector<Pair>& _pairs1, vector<Pair>& _pairsComp, vector<Match>& _matches, 
						unordered_map<uint64_t, bool>& _isMatched0, unordered_map<uint64_t, bool>& _isMatched1) : 
						cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config), pairs0(_pairs0), pairs1(_pairs1),
						pairsComp(_pairsComp), matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1),
						matchMap0(_cgc0.shape), matchMap1(_cgc0.shape), matchMapIm0(_cgc0.shape), matchMapIm1(_cgc0.shape),
#ifdef USE_REDUCTION_MATRIX
						reductionMatrix(_cgc0.shape),
#endif
						pivotColumnIndexInput0(_cgc0.shape),
						pivotColumnIndexInput1(_cgc1.shape),
						pivotColumnIndexComp(_cgcComp.shape),
						pivotColumnIndexImage0(_cgc0.shape),
						pivotColumnIndexImage1(_cgc1.shape)
						{}


void Dimension1::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp, vector<Cube>& ctrImage) {
	auto processInput0 = [this, &ctr0]() {
#ifdef RUNTIME
		cout << endl << "input 0: ";
#endif
    	computePairs(ctr0, 0);
#ifdef USE_CLEARING_DIM0
    	enumerateEdges(ctr0, cgc0, pivotColumnIndexInput0);
#endif
	};

	auto processInput1 = [this, &ctr1]() {
#ifdef RUNTIME
    	cout << endl << "input 1: ";
#endif

	    computePairs(ctr1, 1);
#ifdef USE_CLEARING_DIM0
    	enumerateEdges(ctr1, cgc1, pivotColumnIndexInput1);
#endif
	};

	auto processComparison = [this, &ctrComp, &ctrImage]() {
#ifdef RUNTIME
    	cout << endl << "comparison: ";
#endif

	    computePairsComp(ctrComp);
#ifdef USE_CLEARING_DIM0
#ifndef USE_APPARENT_PAIRS_COMP
    	ctrImage = ctrComp;
#endif
	    enumerateEdges(ctrComp, cgcComp, pivotColumnIndexComp);
#endif
	};

	auto processImageInput0Comparison = [this, &ctrComp, &ctrImage]() {
#ifdef RUNTIME
    	cout << endl << "image 0: ";
#endif

#if defined(USE_APPARENT_PAIRS_COMP) or defined(USE_CLEARING_DIM0)
    	computePairsImage(ctrImage, 0);
#else
    	computePairsImage(ctrComp, 0);
#endif
	};

	auto processImageInput1Comparison = [this, &ctrComp, &ctrImage]() {
#ifdef RUNTIME
		cout << endl << "image 1: ";
#endif

#if defined(USE_APPARENT_PAIRS_COMP) or defined(USE_CLEARING_DIM0)
          computePairsImage(ctrImage, 1);
#else
          computePairsImage(ctrComp, 1);
#endif
	};

#ifdef PARALLELIZE_INDEPENDENT_BARCODES_DIM1
	auto comparisonFuture = std::async(processComparison);
	auto input0Future = std::async(processInput0);
	auto input1Future = std::async(processInput1);
#if defined(USE_APPARENT_PAIRS_COMP) or defined(USE_CLEARING_DIM0) or defined(USE_ISPAIRED)
	comparisonFuture.wait();
#endif
	auto imageInput0ComparisonFuture = std::async(processImageInput0Comparison);
	auto imageInput1ComparisonFuture = std::async(processImageInput1Comparison);
	input0Future.wait();
	input1Future.wait();
	imageInput0ComparisonFuture.wait();
	imageInput1ComparisonFuture.wait();
#if not(defined(USE_APPARENT_PAIRS_COMP) or defined(USE_CLEARING_DIM0) or defined(USE_ISPAIRED))
	comparisonFuture.wait();
#endif
#else
	processInput0();
	processInput1();
	processComparison();
	processImageInput0Comparison();
	processImageInput1Comparison();
#endif

#ifdef RUNTIME
	cout << endl << "matching: ";
#endif
	computeMatching();
}

void Dimension1::computeInput0Pairs(vector<Cube>& ctr0)  {
	computePairs(ctr0, 0);
#ifdef USE_CLEARING_DIM0
	enumerateEdges(ctr0, cgc0, pivotColumnIndexInput0);
#endif
}

vector<vector<index_t>> Dimension1::getRepresentativeCycle(const Pair& pair, const CubicalGridComplex& cgc) {
	vector<vector<index_t>> reprCycle;
	reprCycle.push_back(cgc.getParentVoxel(pair.birth, 1));
	
	vector<Cube> ctr;
	enumerateColumnsToReduce(ctr, cgc);
	size_t ctrSize = ctr.size();
	auto& pivotColumnIndex = pivotColumnIndexInput0; // Use input 0 pivot column index (doesn't really matter which one for this method)
	pivotColumnIndex.clear();
	BoundaryEnumerator enumerator(cgc);
	Cube pivot;
	size_t j;
	vector<index_t> vertex;
#ifdef USE_REDUCTION_MATRIX
	reductionMatrix.clear();
	vector<Cube> reductionColumn;
#endif
#ifdef USE_CACHE
	CubeMap<2, vector<Cube>> cache(cgc.shape);
	queue<uint64_t> cachedColumnIdx;
	size_t numRecurse;
#endif
#ifdef USE_EMERGENT_PAIRS
	bool checkEmergentPair;
#endif
#ifdef USE_EMERGENT_PAIRS
	vector<Cube> faces;
	BoundaryEnumerator enumeratorAP(cgc);
	CoboundaryEnumerator coEnumeratorAP(cgc);
#endif

	for (size_t i = 0; i < ctrSize; ++i) {
		if (pivot == pair.birth) { break; }
		CubeQueue workingBoundary;
		j = i;
#ifdef USE_CACHE
		numRecurse = 0;
#endif
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
			if (j == i) {
#ifdef USE_EMERGENT_PAIRS
				if (isEmergentPair<INPUT_PAIRS>(ctr[i], pivot, j, faces, checkEmergentPair, cgc, enumerator, enumeratorAP, coEnumeratorAP, pivotColumnIndex)) {
					pivotColumnIndex.emplace(pivot.index, i);
					break;
				} else {
					for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_CACHE
					++numRecurse;
#endif
					if (j != i) { continue; }
#ifdef USE_REDUCTION_MATRIX
					else { reductionColumn.push_back(coEnumeratorAP.nextCoface); }
#endif
				}
#else			
				enumerator.setBoundaryEnumerator(ctr[i]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			} else {
#ifdef USE_REDUCTION_MATRIX
				reductionColumn.push_back(ctr[j]);
#endif
#ifdef USE_CACHE
				if (!columnIsCached(ctr[j], workingBoundary, cache)) {
#endif
					enumerator.setBoundaryEnumerator(ctr[j]);
					while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#ifdef USE_REDUCTION_MATRIX
					useReductionMatrix(ctr[j], workingBoundary, enumerator
#ifdef USE_CACHE
										, cache
#endif
										);
#endif
#ifdef USE_CACHE
				}
#endif
			}
			pivot = getPivot(workingBoundary);
			if (pivot.index != NONE_INDEX) {
				auto it = pivotColumnIndex.find(pivot.index);
				if (it.has_value()) {
					j = *it;
#ifdef USE_CACHE
					++numRecurse;
#endif
					continue;
				} else {
					if (pivot == pair.birth) {
						Cube c;
						while (!workingBoundary.empty()) {
							c = workingBoundary.top();
							workingBoundary.pop();
							if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
							else {
								vertex = {c.x(), c.y(), c.z()};
								if(find(reprCycle.begin(), reprCycle.end(), vertex) == reprCycle.end()) { reprCycle.push_back(vertex); }
								switch(c.type()) {
									case 0:
										vertex = {c.x()+1, c.y(), c.z()};
										if(find(reprCycle.begin(), reprCycle.end(), vertex) == reprCycle.end()) { reprCycle.push_back(vertex); }
										break;

									case 1:
										vertex = {c.x(), c.y()+1, c.z()};
										if(find(reprCycle.begin(), reprCycle.end(), vertex) == reprCycle.end()) { reprCycle.push_back(vertex); }
										break;

									case 2:
										vertex = {c.x(), c.y(), c.z()+1};
										if(find(reprCycle.begin(), reprCycle.end(), vertex) == reprCycle.end()) { reprCycle.push_back(vertex); }
										break;
								}
								
							}
						}
						break;
					}
					pivotColumnIndex.emplace(pivot.index, i);
#ifdef USE_CACHE
					if (numRecurse >= config.minRecursionToCache) {
						addCache(ctr[i], workingBoundary, cachedColumnIdx, cache);
						break;
					}
#endif
#ifdef USE_REDUCTION_MATRIX
					if (reductionColumn.size() > 0) {
						reductionMatrix.emplace(ctr[i].index, reductionColumn);
						reductionColumn.clear();
					}
#endif
					break;
				}
			} else { break; }
		}
	}

	reprCycle.push_back(cgc.getParentVoxel(pair.death, 2));

	return reprCycle;
}


void Dimension1::computePairs(vector<Cube>& ctr, uint8_t k) {
	computePairsUnified<ComputePairsMode::INPUT_PAIRS>(ctr, k);
}

void Dimension1::computePairsComp(vector<Cube>& ctr) {
	computePairsUnified<ComputePairsMode::COMPARISON_PAIRS>(ctr, 0);
}

void Dimension1::computePairsImage(vector<Cube>& ctr, uint8_t k) {
	computePairsUnified<ComputePairsMode::IMAGE_PAIRS>(ctr, k);
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


void Dimension1::enumerateEdges(vector<Cube>& edges, const CubicalGridComplex& cgc, CubeMap<1, size_t>& pivotColumnIndex) const {
#ifdef RUNTIME
	cout << "; enumeration ";
	auto start = high_resolution_clock::now();
#endif

	edges.clear();
	edges.reserve(cgc.getNumberOfCubes(1));
	value_t birth;
#ifdef USE_CLEARING_DIM0
	Cube cube;
#endif
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
	bool binaryInputs = true;
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
						if (!find.has_value())
						{
							edges.push_back(cube);
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
							if (binaryInputs && birth != 0 && birth != 1) binaryInputs = false;
#endif
						}
#else
						edges.push_back(Cube(birth, x, y, z, type));
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
						if (binaryInputs && birth != 0 && birth != 1) binaryInputs = false;
#endif
#endif
					}	
				}				
			}
		}
	}
#ifdef USE_STABLE_SORT_OR_STABLE_PARTITION
	if (binaryInputs) {
		std::stable_partition(edges.begin(), edges.end(), [](Cube &cube) { return cube.birth == 0; });
	} else {
		std::stable_sort(edges.begin(), edges.end(), [](const Cube &cube1, const Cube &cube2) { return cube1.birth < cube2.birth; });
	}
#else
	std::sort(edges.begin(), edges.end(), CubeComparator());
#endif

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#endif
}


void Dimension1::enumerateColumnsToReduce(vector<Cube>& ctr, const CubicalGridComplex& cgc) const {
	ctr.reserve(cgc.getNumberOfCubes(2));
	value_t birth;
	for (index_t x = 0; x < cgc.shape[0]; ++x) {
		for (index_t y = 0; y < cgc.shape[1]; ++y) {
			for (index_t z = 0; z < cgc.shape[2]; ++z) {
				for (uint8_t type = 0; type < 3; ++type) {
					birth = cgc.getBirth(x, y, z, type, 2);
					if (birth < config.threshold) { ctr.push_back(Cube(birth, x, y, z, type)); }
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
	if (result.index != NONE_INDEX) { column.push(result); }
	return result;
}

#ifdef USE_REDUCTION_MATRIX
void Dimension1::useReductionMatrix(const Cube& column, CubeQueue& workingBoundary, BoundaryEnumerator& enumerator
#ifdef USE_CACHE
									, CubeMap<2, vector<Cube>>& cache
#endif
									) const {
	auto reductionColumn = reductionMatrix.find(column.index);
	if (reductionColumn.has_value()) {
		for (Cube& row : *reductionColumn) {
#ifdef USE_CACHE
			if (!columnIsCached(row, workingBoundary, cache)) {
#endif
				enumerator.setBoundaryEnumerator(row);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
				useReductionMatrix(row, workingBoundary, enumerator
#ifdef USE_CACHE
									, cache
#endif
									);
#ifdef USE_CACHE
			}
#endif
		}
	}
}
#endif

#ifdef USE_CACHE
bool Dimension1::columnIsCached(const Cube& column, CubeQueue& workingBoundary, CubeMap<2, vector<Cube>>& cache) const {
	auto& cachedBoundary = cache.find(column.index);
	if (cachedBoundary.has_value()) {
		for (auto &face : *cachedBoundary) {
			workingBoundary.push(face);
		}
		return true;
	} else { return false; }
}

void Dimension1::addCache(const Cube& column, CubeQueue& workingBoundary, queue<uint64_t>& cachedColumnIdx, CubeMap<2, vector<Cube>>& cache) {
	std::vector<Cube> cleanWb;
	while (!workingBoundary.empty()) {
		Cube c = workingBoundary.top();
		workingBoundary.pop();
		if (!workingBoundary.empty() && c == workingBoundary.top()) { workingBoundary.pop(); } 
		else { cleanWb.emplace_back(c); }
	}

	cache[column.index] = std::move(cleanWb);
	cachedColumnIdx.push(column.index);
	if (cachedColumnIdx.size() > config.cacheSize) {
		cache[cachedColumnIdx.front()] = {};
		cachedColumnIdx.pop();
	}
}
#endif

#ifdef USE_EMERGENT_PAIRS
template<Dimension1::ComputePairsMode computePairsMode>
bool Dimension1::isEmergentPair(const Cube&column, Cube& pivot, size_t& j, vector<Cube>& faces, bool& checkEmergentPair,
		const CubicalGridComplex& cgc, BoundaryEnumerator& enumerator, BoundaryEnumerator& enumeratorAP,
		CoboundaryEnumerator& coEnumeratorAP, CubeMap<1, size_t>& pivotColumnIndex) const {
    auto birth = (computePairsMode == IMAGE_PAIRS) ? cgc.getBirth(column.x(), column.y(), column.z(), column.type(), 2) : column.birth;

    const bool useApparentPairs =
#ifdef USE_APPARENT_PAIRS
        true;
#else
        false;
#endif
    const bool useApparentPairsComp =
#ifdef USE_APPARENT_PAIRS_COMP
        true;
#else
        false;
#endif

	faces.clear();
	enumerator.setBoundaryEnumerator(column);
	while (enumerator.hasPreviousFace()) {
		if (checkEmergentPair && enumerator.nextFace.birth == birth) {
			auto nextColumnIndex = pivotColumnIndex.find(enumerator.nextFace.index);
            if ((useApparentPairs && computePairsMode == INPUT_PAIRS)
                    || (useApparentPairsComp && computePairsMode == COMPARISON_PAIRS)
                    || computePairsMode == IMAGE_PAIRS) {
				if (nextColumnIndex.has_value()) {
                    checkEmergentPair = false;
                    j = *nextColumnIndex;
                }
#if defined(USE_APPARENT_PAIRS) or defined(USE_APPARENT_PAIRS_COMP)
				else if ((computePairsMode == INPUT_PAIRS || computePairsMode == COMPARISON_PAIRS)
                            && pivotOfColumnIsApparentPair(enumerator.nextFace, column, faces, enumeratorAP, coEnumeratorAP)) {
                    checkEmergentPair = false;
                }
#endif
				else {
                    pivot = enumerator.nextFace;
                    return true;
                }
            }
            else if (computePairsMode == INPUT_PAIRS || computePairsMode == COMPARISON_PAIRS) {
                if (nextColumnIndex.has_value()) {
                    checkEmergentPair = false;
                    j = *nextColumnIndex;
                } else {
                    pivot = enumerator.nextFace;
                    return true;
                }
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

template <Dimension1::ComputePairsMode computePairsMode>
void Dimension1::computePairsUnified(vector<Cube>& ctr, uint8_t k) {
#ifdef RUNTIME
	cout << "barcode ";
	auto start = high_resolution_clock::now();
#endif
    const CubicalGridComplex& cgc = (computePairsMode == COMPARISON_PAIRS) ? cgcComp : ((k == 0) ? cgc0 : cgc1);
    vector<Pair>& pairs = (computePairsMode == COMPARISON_PAIRS) ? pairsComp : ((k == 0) ? pairs0 : pairs1);
    CubeMap<1, Pair>& matchMap = (computePairsMode == INPUT_PAIRS) ? ((k == 0) ? matchMap0 : matchMap1) : matchMap0;
    CubeMap<1, uint64_t>& matchMapIm = (computePairsMode == IMAGE_PAIRS) ? ((k==0) ? matchMapIm0 : matchMapIm1) : matchMapIm0;
	CubeMap<1, size_t>& pivotColumnIndex =
		(computePairsMode == INPUT_PAIRS) ? ((k == 0) ? pivotColumnIndexInput0 : pivotColumnIndexInput1) :
		(computePairsMode == COMPARISON_PAIRS) ? pivotColumnIndexComp :
		((k == 0) ? pivotColumnIndexImage0 : pivotColumnIndexImage1);

    const bool useApparentPairs = 
#ifdef USE_APPARENT_PAIRS
        true;
#else
        false;
#endif
    const bool useApparentPairsComp = 
#ifdef USE_APPARENT_PAIRS_COMP
        true;
#else
        false;
#endif

	size_t ctrSize = ctr.size();
	BoundaryEnumerator enumerator(cgc);
	Cube pivot;
	size_t j;
#ifdef USE_REDUCTION_MATRIX
	reductionMatrix.clear();
	vector<Cube> reductionColumn;
#ifdef RUNTIME
	size_t numReductionColumns = 0;
#endif
#endif
#ifdef USE_CACHE
	CubeMap<2, vector<Cube>> cache(cgc.shape);
	queue<uint64_t> cachedColumnIdx;
	size_t numRecurse;
#ifdef RUNTIME
	size_t numCached = 0;
#endif
#endif
#ifdef USE_EMERGENT_PAIRS
	bool checkEmergentPair;
#ifdef RUNTIME
	size_t numEmergentPairs = 0;
#endif
#endif
	vector<Cube> faces;
	BoundaryEnumerator enumeratorAP(cgc);
	CoboundaryEnumerator coEnumeratorAP(cgc);
	bool shouldClear;
    if (computePairsMode == COMPARISON_PAIRS) {
#if defined (USE_CLEARING_IMAGE) and not defined(USE_APPARENT_PAIRS_COMP)
	    shouldClear = false;
#endif
    }
    if (computePairsMode == IMAGE_PAIRS) {
#ifdef USE_CLEARING_IMAGE
	    shouldClear = false;
#endif
    }

	for (size_t i = 0; i < ctrSize; ++i) {
		CubeQueue workingBoundary;
		j = i;
#ifdef USE_CACHE
		numRecurse = 0;
#endif
#ifdef USE_EMERGENT_PAIRS
		checkEmergentPair = true;
#endif
		while (true) {
			if (j == i) {
#ifdef USE_EMERGENT_PAIRS
				if (isEmergentPair<computePairsMode>(ctr[i], pivot, j, faces, checkEmergentPair, cgc, enumerator, enumeratorAP, coEnumeratorAP, pivotColumnIndex)) {
					pivotColumnIndex.emplace(pivot.index, i);

                    if (computePairsMode == IMAGE_PAIRS
#ifdef USE_ISPAIRED
						 && isPairedComp[ctr[i].index]
#endif
					) {
						matchMapIm.emplace(ctr[i].index, pivot.index);
					}


#ifdef RUNTIME
					++numEmergentPairs;
#endif
					break;
				} else {
					for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_CACHE
					++numRecurse;
#endif
					if (j != i) { continue; }
#ifdef USE_REDUCTION_MATRIX
					else if (computePairsMode == INPUT_PAIRS || computePairsMode == COMPARISON_PAIRS) { reductionColumn.push_back(coEnumeratorAP.nextCoface); }
#endif
				}
#else			
				enumerator.setBoundaryEnumerator(ctr[i]);
				while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#endif
			} else {
#ifdef USE_REDUCTION_MATRIX
				reductionColumn.push_back(ctr[j]);
#endif
#ifdef USE_CACHE
				if (!columnIsCached(ctr[j], workingBoundary, cache)) {
#endif
					enumerator.setBoundaryEnumerator(ctr[j]);
					while (enumerator.hasNextFace()) { workingBoundary.push(enumerator.nextFace); }
#ifdef USE_REDUCTION_MATRIX
					useReductionMatrix(ctr[j], workingBoundary, enumerator
#ifdef USE_CACHE
									, cache
#endif
									);
#endif
#ifdef USE_CACHE
				}
#endif
			}
			pivot = getPivot(workingBoundary);

#if defined(USE_APPARENT_PAIRS) or defined(USE_APPARENT_PAIRS_COMP)
            if ((useApparentPairs && computePairsMode == INPUT_PAIRS) || (useApparentPairsComp && computePairsMode == COMPARISON_PAIRS)) {
                while (true) {
                    faces.clear();
                    if (pivotIsApparentPair(pivot, faces, enumeratorAP, coEnumeratorAP)) {
                        for (auto face = faces.rbegin(), last = faces.rend(); face != last; ++face) { workingBoundary.push(*face); }
#ifdef USE_REDUCTION_MATRIX
                        reductionColumn.push_back(coEnumeratorAP.nextCoface);
#endif
#ifdef USE_CACHE
                        ++numRecurse;
#endif
                        pivot = getPivot(workingBoundary);
                    } else { break; }
                }
            }
#endif
			if (pivot.index != NONE_INDEX) {
				auto cachedIndex = pivotColumnIndex.find(pivot.index);
				if (cachedIndex.has_value()) {
					j = *cachedIndex;
#ifdef USE_CACHE
					++numRecurse;
#endif
					continue;
				} else {
					pivotColumnIndex.emplace(pivot.index, i);
					if (computePairsMode == INPUT_PAIRS || computePairsMode == COMPARISON_PAIRS) {
						if (pivot.birth != ctr[i].birth) {
							pairs.push_back(Pair(pivot, ctr[i]));
							if (computePairsMode == INPUT_PAIRS) {
								matchMap.emplace(pivot.index, pairs.back());
							}
#ifdef USE_ISPAIRED
							if (computePairsMode == COMPARISON_PAIRS) {
								isPairedComp.emplace(ctr[i].index, true);
							}
#endif
						}
					}
					if (computePairsMode == IMAGE_PAIRS
#ifdef USE_ISPAIRED
						&& isPairedComp[ctr[i].index]
#endif
					) {
						matchMapIm.emplace(ctr[i].index, pivot.index);
					}
#ifdef USE_CACHE
					if (numRecurse >= config.minRecursionToCache) {
						addCache(ctr[i], workingBoundary, cachedColumnIdx, cache);
#ifdef RUNTIME
						++numCached;
#endif
						break;
					}
#endif
#ifdef USE_REDUCTION_MATRIX
					if (reductionColumn.size() > 0) {
						reductionMatrix.emplace(ctr[i].index, reductionColumn);
						reductionColumn.clear();
#ifdef RUNTIME
						++numReductionColumns;
#endif
					}
#endif
					break;
				}
			} else {
#if defined(USE_CLEARING_IMAGE)
                if ((computePairsMode == COMPARISON_PAIRS && !useApparentPairsComp) || computePairsMode == IMAGE_PAIRS) {
                    ctr[i].index = NONE_INDEX;
                    shouldClear = true;
                }
#endif
                break;
            }
		}
	}

#if defined(USE_CLEARING_IMAGE)
	if ((computePairsMode == COMPARISON_PAIRS && !useApparentPairsComp) || computePairsMode == IMAGE_PAIRS) {
		if (shouldClear) {
			auto newEnd = remove_if(ctr.begin(), ctr.end(), [](const Cube& cube) { return cube.index == NONE_INDEX; });
			ctr.erase(newEnd, ctr.end());
		}
	}
#endif

#ifdef RUNTIME
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << duration.count() << " ms";
#ifdef USE_REDUCTION_MATRIX
	cout << ", " << numReductionColumns << " reduction columns";
#endif
#ifdef USE_CACHE
	cout << ", " << numCached << " cached columns";
#endif
#ifdef USE_EMERGENT_PAIRS
	cout << ", " << numEmergentPairs << " emergent pairs";
#endif
#endif
}
