#include "compute_pairs.h"
#include "enumerators.h"

#include <iostream>

ComputePairs::ComputePairs(Config& _config) : config(_config) {}

void ComputePairs::computePairs(const CubicalGridComplex& cgc, vector<Cube>& ctr, bool image) {
	uint64_t num_reduction_steps = 0;
	auto ctr_size = ctr.size();
	unordered_map<uint64_t, uint64_t> pivot_column_index;
	pivot_column_index.reserve(ctr_size);
	unordered_map<uint64_t, CubeQue > cache;
	cache.reserve(ctr_size);
	queue<uint64_t> cached_column_idx;
	BoundaryEnumerator faces = BoundaryEnumerator(cgc);
	for(uint64_t i = 0; i < ctr_size; i++) {
		CubeQue working_boundary;
		uint64_t j = i;
		int num_recurse = 0;
		while (true) {
			bool cache_hit = false;
			if (i!=j) {
                auto findCb = cache.find(j);
                if (findCb != cache.end()) {
                    cache_hit = true;
                    auto cached_boundary = findCb -> second;
                    while (!cached_boundary.empty()) {
                        working_boundary.push(cached_boundary.top());
                        cached_boundary.pop();
                    }
                }
			}
			if (!cache_hit) {
				faces.setBoundaryEnumerator(ctr[j]);
				while (faces.hasNextFace()) {
					working_boundary.push(faces.getFace());
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
					if (image) {
						pairs->push_back(WritePair(pivot, ctr[i]));
					} else if (pivot.birth != ctr[i].birth) {
						pairs->push_back(WritePair(pivot, ctr[i]));
					}
					break;
				}
			} else if (!image) { 
				clearing = true;
				ctr[i].index = NONE;
				break; 
			}
		}
	}
	if (config.verbose) {cout << " with " << num_reduction_steps << " reduction steps ... ";}
	if (clearing) {
		auto new_end = std::remove_if(ctr.begin(), ctr.end(),
								[](const Cube& c){ return c.index == NONE; });
		ctr.erase(new_end, ctr.end());
	}
}

Cube ComputePairs::popPivot(CubeQue& column) const {
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

Cube ComputePairs::getPivot(CubeQue& column) const {
	Cube result = popPivot(column);
	if (result.coordinates[0] != NONE) {
		column.push(result);
	}
	return result;
}

void ComputePairs::addCache(uint64_t i, CubeQue& working_boundary, unordered_map<uint64_t, CubeQue>& cache) const {
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