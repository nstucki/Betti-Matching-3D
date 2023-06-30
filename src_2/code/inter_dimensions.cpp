#include "inter_dimensions.h"
#include "enumerators.h"

#include <iostream>



InterDimensions::InterDimensions(const CubicalGridComplex& _cgc0, const CubicalGridComplex& _cgc1, const CubicalGridComplex& _cgcComp, 
							const Config& _config) : cgc0(_cgc0), cgc1(_cgc1), cgcComp(_cgcComp), config(_config) {
	computeDim = cgc0.dim-2;
	pairs0 = vector<vector<Pair>>(computeDim+1);
	pairs1 = vector<vector<Pair>>(computeDim+1);
	pairsComp = vector<vector<Pair>>(computeDim+1);
	matches = vector<vector<Match>>(computeDim+1);
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

Cube InterDimensions::getPivot(CubeQue& column) const {
	Cube result = popPivot(column);
	if (result.coordinates[0] != NONE) {
		column.push(result);
	}
	return result;
}

void InterDimensions::addCache(uint64_t i, CubeQue& working_boundary, unordered_map<uint64_t, CubeQue>& cache) const {
	CubeQue clean_wb;
	while (!working_boundary.empty()) {
		auto c = working_boundary.top();
		working_boundary.pop();
		if (!working_boundary.empty() && c == working_boundary.top()) {
			working_boundary.pop();
		} else {
			clean_wb.push(c);
		}
	}
	cache.emplace(i, clean_wb);
}


void InterDimensions::computePairsComp(vector<Cube>& ctr) {
	uint64_t num_reduction_steps = 0;
	auto ctrSize = ctr.size();
	unordered_map<uint64_t, uint64_t> pivot_column_index;
	pivot_column_index.reserve(ctrSize);
	unordered_map<uint64_t, CubeQue > cache;
	cache.reserve(ctrSize);
	queue<uint64_t> cached_column_idx;
	BoundaryEnumerator faces = BoundaryEnumerator(cgcComp);
	uint64_t j;

	for(uint64_t i = 0; i < ctrSize; i++) {
		CubeQue working_boundary;
		j = i;
		int num_recurse = 0;

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
					num_recurse++;
					num_reduction_steps++;
					continue;
				} else {
					if (num_recurse >= config.min_recursion_to_cache) {
                        addCache(i, working_boundary, cache);
						cached_column_idx.push(i);
						if (cached_column_idx.size() > config.cache_size) {
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
			} else { break; }
		}
	}
	
	if (config.verbose) { cout << " with " << num_reduction_steps << " reduction steps ... "; }
}


void InterDimensions::computePairs(uint8_t k, vector<Cube>& ctr) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	vector<vector<Pair>>& pairs = (k == 0) ? pairs0 : pairs1;
	unordered_map<uint64_t,Pair>& matchMap = (k == 0) ? matchMap0 : matchMap1;
	matchMap.clear();

	uint64_t num_reduction_steps = 0;
	auto ctrSize = ctr.size();
	unordered_map<uint64_t, uint64_t> pivot_column_index;
	pivot_column_index.reserve(ctrSize);
	unordered_map<uint64_t, CubeQue > cache;
	cache.reserve(ctrSize);
	queue<uint64_t> cached_column_idx;
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	uint64_t j;

	for(uint64_t i = 0; i < ctrSize; i++) {
		CubeQue working_boundary;
		j = i;
		int num_recurse = 0;

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
					num_recurse++;
					num_reduction_steps++;
					continue;
				} else {
					if (num_recurse >= config.min_recursion_to_cache) {
                        addCache(i, working_boundary, cache);
						cached_column_idx.push(i);
						if (cached_column_idx.size() > config.cache_size) {
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

	if (config.verbose) { cout << " with " << num_reduction_steps << " reduction steps ... "; }
}


void InterDimensions::computePairsImage(uint8_t k, vector<Cube>& ctr) {
	const CubicalGridComplex& cgc = (k == 0) ? cgc0 : cgc1;
	unordered_map<uint64_t,Cube>& matchMapIm = (k==0) ? matchMapIm0 : matchMapIm1;
	matchMapIm.reserve(pairsComp[computeDim].size());
	
	uint64_t num_reduction_steps = 0;
	auto ctrSize = ctr.size();
	unordered_map<uint64_t, uint64_t> pivot_column_index;
	pivot_column_index.reserve(ctrSize);
	unordered_map<uint64_t, CubeQue > cache;
	cache.reserve(ctrSize);
	queue<uint64_t> cached_column_idx;
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	uint64_t j;

	for(uint64_t i = 0; i < ctrSize; i++) {
		CubeQue working_boundary;
		j = i;
		int num_recurse = 0;

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
					num_recurse++;
					num_reduction_steps++;
					continue;
				} else {
					if (num_recurse >= config.min_recursion_to_cache) {
                        addCache(i, working_boundary, cache);
						cached_column_idx.push(i);
						if (cached_column_idx.size() > config.cache_size) {
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
	
	if (config.verbose) { cout << " with " << num_reduction_steps << " reduction steps ... "; }
}


void InterDimensions::computeMatching() {
	Cube birth0;
	Cube birth1;
	for (auto& pair : pairsComp[computeDim]) {
		auto findIm0 = matchMapIm0.find(cgcComp.getCubeIndex(pair.death));
		auto findIm1 = matchMapIm1.find(cgcComp.getCubeIndex(pair.death));
		if (findIm0 != matchMapIm0.end() && findIm1 != matchMapIm1.end()) {
			birth0 = findIm0->second;
			birth1 = findIm1->second;
			auto find0 = matchMap0.find(cgc0.getCubeIndex(birth0));
			auto find1 = matchMap1.find(cgc1.getCubeIndex(birth1));
			if (find0 != matchMap0.end() && find1 != matchMap1.end()) {
				matches[computeDim].push_back(Match(find0->second, find1->second));
			}
		}
	}
};


void InterDimensions::computePairsAndMatch(vector<Cube>& ctr0, vector<Cube>& ctr1, vector<Cube>& ctrComp) {
	computePairsComp(ctrComp);
	computePairs(0, ctr0);
	computePairs(1, ctr1);
	computePairsImage(0, ctrComp);
	computePairsImage(1, ctrComp);
	computeMatching();
}