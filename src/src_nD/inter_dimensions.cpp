#include "inter_dimensions.h"
#include "enumerators.h"

#include <algorithm>
#include <chrono>
#include <iostream>

using namespace dimN;
using namespace std::chrono;

InterDimensions::InterDimensions(
    const CubicalGridComplex &_cgc0, const CubicalGridComplex &_cgc1,
    const CubicalGridComplex &_cgcComp, const Config &_config,
    vector<vector<Pair>> &_pairs0, vector<vector<Pair>> &_pairs1,
    vector<vector<Pair>> &_pairsComp, vector<vector<Match>> &_matches,
    unordered_map<index_t, bool> &_isMatched0,
    unordered_map<index_t, bool> &_isMatched1)
    : cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config),
      pairs0(_pairs0), pairs1(_pairs1), pairsComp(_pairsComp),
      matches(_matches), isMatched0(_isMatched0), isMatched1(_isMatched1) {
    computeDim = cgc0.dim - 2;
}

void InterDimensions::computePairsAndMatch(vector<Cube> &ctr0,
                                           vector<Cube> &ctr1,
                                           vector<Cube> &ctrComp) {
    vector<Cube> ctrImage;
    while (computeDim > 0) {
#ifdef RUNTIME
        if (computeDim != cgc0.dim - 2) {
            cout << endl << endl;
        }
        cout << "dimension " << computeDim << ":";
        auto start = high_resolution_clock::now();
        cout << endl << "input 0: ";
#endif
        if (computeDim != cgc0.dim - 2) {
            pivotColumnIndex.clear();
            cache.clear();
            matchMap0.clear();
        }
        computePairs(ctr0, 0);
#ifdef USE_CLEARING_DIM0
        ctr0.clear();
        assembleColumnsToReduce(cgc0, ctr0);
#else
        if (computeDim != 1) {
            ctr0.clear();
            assembleColumnsToReduce(cgc0, ctr0);
        }
#endif

#ifdef RUNTIME
        cout << endl << "input 1: ";
#endif
        pivotColumnIndex.clear();
        cache.clear();
        if (computeDim != cgc0.dim - 2) {
            matchMap1.clear();
        }
        computePairs(ctr1, 1);
#ifdef USE_CLEARING_DIM0
        ctr1.clear();
        assembleColumnsToReduce(cgc1, ctr1);
#else
        if (computeDim != 1) {
            ctr1.clear();
            assembleColumnsToReduce(cgc1, ctr1);
        }
#endif

#ifdef RUNTIME
        cout << endl << "comparison: ";
#endif
        pivotColumnIndex.clear();
        cache.clear();
        if (computeDim != cgc0.dim - 2) {
            isMatchedComp.clear();
        }
        computePairsComp(ctrComp);
#ifdef USE_CLEARING_DIM0
        ctrImage = ctrComp;
        ctrComp.clear();
        assembleColumnsToReduce(cgcComp, ctrComp);
#else
        if (computeDim != 1) {
            ctrImage = ctrComp;
            ctrComp.clear();
            assembleColumnsToReduce(cgcComp, ctrComp);
        }
#endif

#ifdef RUNTIME
        cout << endl << "image 0: ";
#endif
        pivotColumnIndex.clear();
        cache.clear();
        if (computeDim != cgc0.dim - 2) {
            matchMapIm0.clear();
        }
#ifdef USE_CLEARING_DIM0
        computeImagePairs(ctrImage, 0);
#else
        if (computeDim != 1) {
            computeImagePairs(ctrImage, 0);
        } else {
            computeImagePairs(ctrComp, 0);
        }
#endif

#ifdef RUNTIME
        cout << endl << "image 1: ";
#endif
        pivotColumnIndex.clear();
        cache.clear();
        if (computeDim != cgc0.dim - 2) {
            matchMapIm1.clear();
        }
#ifdef USE_CLEARING_DIM0
        computeImagePairs(ctrImage, 1);
#else
        if (computeDim != 1) {
            computeImagePairs(ctrImage, 1);
        } else {
            computeImagePairs(ctrComp, 1);
        }
#endif

#ifdef RUNTIME
        cout << endl << "matching: ";
#endif
        computeMatching();

        --computeDim;
    }
}

void InterDimensions::computeInput0Pairs(vector<Cube> &ctr0) {
    vector<Cube> ctrImage;
    while (computeDim > 0) {
        if (computeDim != cgc0.dim - 2) {
            pivotColumnIndex.clear();
            cache.clear();
            matchMap0.clear();
        }
        computePairs(ctr0, 0);
#ifdef USE_CLEARING_DIM0
        ctr0.clear();
        assembleColumnsToReduce(cgc0, ctr0);
#else
        if (computeDim != 1) {
            ctr0.clear();
            assembleColumnsToReduce(cgc0, ctr0);
        }
#endif

        --computeDim;
    }
}

void InterDimensions::computePairs(const vector<Cube> &ctr, uint8_t k) {
#ifdef RUNTIME
    cout << "barcode ";
    auto start = high_resolution_clock::now();
#endif
    const CubicalGridComplex &cgc = (k == 0) ? cgc0 : cgc1;
    vector<vector<Pair>> &pairs = (k == 0) ? pairs0 : pairs1;
    unordered_map<index_t, Pair> &matchMap = (k == 0) ? matchMap0 : matchMap1;

    size_t ctrSize = ctr.size();
    pivotColumnIndex.reserve(ctrSize);
    cache.reserve(min(config.cacheSize, ctrSize));
    BoundaryEnumerator enumerator = BoundaryEnumerator(cgc);
    Cube pivot;
    queue<index_t> cachedColumnIdx;
    size_t numRecurse;
    index_t j;
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
#ifdef USE_EMERGENT_PAIRS
                faces.clear();
#endif
                enumerator.setBoundaryEnumerator(ctr[j], computeDim + 1);
                while (enumerator.hasNextFace()) {
#ifdef USE_EMERGENT_PAIRS
                    if (checkEmergentPair &&
                        ctr[j].birth == enumerator.nextFace.birth) {
                        if (pivotColumnIndex.find(
                                cgc.getCubeIndex(enumerator.nextFace)) ==
                            pivotColumnIndex.end()) {
                            pivot = enumerator.nextFace;
                            foundPair = true;
                            break;
                        } else {
                            checkEmergentPair = false;
                        }
                    }
                    faces.push_back(enumerator.nextFace);
#else
                    workingBoundary.push(enumerator.nextFace);
#endif
                }
#ifdef USE_EMERGENT_PAIRS
                if (foundPair) {
                    pivotColumnIndex.emplace(cgc.getCubeIndex(pivot), i);
                    ++numEmergentPairs;
                    break;
                }
                for (auto face = faces.rbegin(), last = faces.rend();
                     face != last; ++face) {
                    workingBoundary.push(*face);
                }
#endif
            }
            pivot = getPivot(workingBoundary);
            if (pivot.coordinates[0] != NONE) {
                auto pair = pivotColumnIndex.find(cgc.getCubeIndex(pivot));
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
                    pivotColumnIndex[cgc.getCubeIndex(pivot)] = i;
                    if (pivot.birth != ctr[i].birth) {
                        pairs[computeDim].push_back(Pair(pivot, ctr[i]));
                        matchMap.emplace(cgc.getCubeIndex(pivot),
                                         pairs[computeDim].back());
                    }
                    break;
                }
            } else {
                break;
            }
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

void InterDimensions::computePairsComp(vector<Cube> &ctr) {
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
    index_t j;
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
#ifdef USE_EMERGENT_PAIRS
                faces.clear();
#endif
                enumerator.setBoundaryEnumerator(ctr[j], computeDim + 1);
                while (enumerator.hasNextFace()) {
#ifdef USE_EMERGENT_PAIRS
                    if (checkEmergentPair &&
                        ctr[j].birth == enumerator.nextFace.birth) {
                        if (pivotColumnIndex.find(
                                cgcComp.getCubeIndex(enumerator.nextFace)) ==
                            pivotColumnIndex.end()) {
                            pivot = enumerator.nextFace;
                            foundPair = true;
                            break;
                        } else {
                            checkEmergentPair = false;
                        }
                    }
                    faces.push_back(enumerator.nextFace);
#else
                    workingBoundary.push(enumerator.nextFace);
#endif
                }
#ifdef USE_EMERGENT_PAIRS
                if (foundPair) {
                    pivotColumnIndex.emplace(cgcComp.getCubeIndex(pivot), i);
                    ++numEmergentPairs;
                    break;
                }
                for (auto face = faces.rbegin(), last = faces.rend();
                     face != last; ++face) {
                    workingBoundary.push(*face);
                }
#endif
            }
            Cube pivot = getPivot(workingBoundary);
            if (pivot.coordinates[0] != NONE) {
                auto pair = pivotColumnIndex.find(cgcComp.getCubeIndex(pivot));
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
                    pivotColumnIndex.emplace(cgcComp.getCubeIndex(pivot), i);
                    if (pivot.birth != ctr[i].birth) {
                        pairsComp[computeDim].push_back(Pair(pivot, ctr[i]));
                        isMatchedComp.emplace(cgcComp.getCubeIndex(ctr[i]),
                                              true);
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
        auto newEnd = remove_if(ctr.begin(), ctr.end(), [](const Cube &cube) {
            return cube.coordinates[0] == NONE;
        });
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

void InterDimensions::computeImagePairs(vector<Cube> &ctr, uint8_t k) {
#ifdef RUNTIME
    cout << "barcode ";
    auto start = high_resolution_clock::now();
#endif
    const CubicalGridComplex &cgc = (k == 0) ? cgc0 : cgc1;
    unordered_map<index_t, Cube> &matchMapIm =
        (k == 0) ? matchMapIm0 : matchMapIm1;

    size_t ctrSize = ctr.size();
    pivotColumnIndex.reserve(ctrSize);
    cache.reserve(min(config.cacheSize, ctrSize));
    matchMapIm.reserve(pairsComp[computeDim].size());
    BoundaryEnumerator enumerator = BoundaryEnumerator(cgc);
    Cube pivot;
    value_t birth;
    queue<index_t> cachedColumnIdx;
    size_t numRecurse;
    index_t j;
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
                enumerator.setBoundaryEnumerator(ctr[j], computeDim + 1);
#ifdef USE_EMERGENT_PAIRS
                faces.clear();
#endif
                while (enumerator.hasNextFace()) {
#ifdef USE_EMERGENT_PAIRS
                    if (checkEmergentPair) {
                        birth = cgc.getBirth(ctr[j].coordinates);
                        if (birth == enumerator.nextFace.birth) {
                            if (pivotColumnIndex.find(
                                    cgc.getCubeIndex(enumerator.nextFace)) ==
                                pivotColumnIndex.end()) {
                                pivot = enumerator.nextFace;
                                foundPair = true;
                                break;
                            } else {
                                checkEmergentPair = false;
                            }
                        }
                    }
                    faces.push_back(enumerator.nextFace);
#else
                    workingBoundary.push(enumerator.nextFace);
#endif
                }
#ifdef USE_EMERGENT_PAIRS
                if (foundPair) {
                    pivotColumnIndex.emplace(cgc.getCubeIndex(pivot), i);
                    if (isMatchedComp[cgcComp.getCubeIndex(ctr[i])]) {
                        matchMapIm.emplace(cgcComp.getCubeIndex(ctr[i]), pivot);
                    }
                    ++numEmergentPairs;
                    break;
                }
                for (auto face = faces.rbegin(), last = faces.rend();
                     face != last; ++face) {
                    workingBoundary.push(*face);
                }
#endif
            }
            pivot = getPivot(workingBoundary);
            if (pivot.coordinates[0] != NONE) {
                auto pair = pivotColumnIndex.find(cgc.getCubeIndex(pivot));
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
                    pivotColumnIndex[cgc.getCubeIndex(pivot)] = i;
                    if (isMatchedComp[cgcComp.getCubeIndex(ctr[i])]) {
                        matchMapIm.emplace(cgcComp.getCubeIndex(ctr[i]), pivot);
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
        auto newEnd = remove_if(ctr.begin(), ctr.end(), [](const Cube &cube) {
            return cube.coordinates[0] == NONE;
        });
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

void InterDimensions::computeMatching() {
#ifdef RUNTIME
    auto start = high_resolution_clock::now();
#endif
    Cube birth0;
    Cube birth1;
    Pair pair0;
    Pair pair1;
    for (Pair &pair : pairsComp[computeDim]) {
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
#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms";
#endif
}

void InterDimensions::assembleColumnsToReduce(const CubicalGridComplex &cgc,
                                              vector<Cube> &ctr) const {
#ifdef RUNTIME
    cout << ", enumeration ";
    auto start = high_resolution_clock::now();
#endif
    ctr.reserve(cgc.getNumberOfCubes(computeDim));
    CubeEnumerator cubeEnum(cgc, computeDim);
    Cube cube = cubeEnum.getNextCube();
    if (cube.birth < config.threshold) {
        auto find = pivotColumnIndex.find(cgc.getCubeIndex(cube));
        if (find == pivotColumnIndex.end()) {
            ctr.push_back(cube);
        }
    }
    while (cubeEnum.hasNextCube()) {
        cube = cubeEnum.getNextCube();
        if (cube.birth < config.threshold) {
            auto find = pivotColumnIndex.find(cgc.getCubeIndex(cube));
            if (find == pivotColumnIndex.end()) {
                ctr.push_back(cube);
            }
        }
    }
    sort(ctr.begin(), ctr.end(), CubeComparator());
#ifdef RUNTIME
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << " ms";
#endif
}

Cube InterDimensions::popPivot(CubeQueue &column) const {
    if (column.empty()) {
        return Cube();
    } else {
        Cube pivot = column.top();
        column.pop();
        while (!column.empty() && column.top() == pivot) {
            column.pop();
            if (column.empty()) {
                return Cube();
            } else {
                pivot = column.top();
                column.pop();
            }
        }
        return pivot;
    }
}

Cube InterDimensions::getPivot(CubeQueue &column) const {
    Cube result = popPivot(column);
    if (result.coordinates[0] != NONE) {
        column.push(result);
    }
    return result;
}

void InterDimensions::addCache(const index_t &i, CubeQueue &workingBoundary) {
    CubeQueue cleanWb;
    while (!workingBoundary.empty()) {
        Cube c = workingBoundary.top();
        workingBoundary.pop();
        if (!workingBoundary.empty() && c == workingBoundary.top()) {
            workingBoundary.pop();
        } else {
            cleanWb.push(c);
        }
    }
    cache.emplace(i, cleanWb);
}
