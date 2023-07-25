#include "dimension_1.h"
#include "enumerators.h"

#include <iostream>

#include <chrono>

using namespace std::chrono;

Dimension1::Dimension1(CubicalGridComplex* _cgc, vector<Pair> &_pairs, Config &_config) : cgc(_cgc), pairs(&_pairs), config(_config) {}

void Dimension1::compute_pairs(vector<Cube>& ctr, bool image) {
	uint32_t num_reduction_steps = 0;
	auto ctr_size = ctr.size();
	unordered_map<uint64_t, uint32_t> pivot_column_index;
	pivot_column_index.reserve(ctr_size);
	unordered_map<uint32_t, CubeQue > cache;
	cache.reserve(ctr_size);
	queue<uint32_t> cached_column_idx;
	BoundaryEnumerator faces(cgc);
	bool shouldClear = false;

	for(uint32_t i = 0; i < ctr_size; i++) {
		CubeQue working_boundary;
		uint32_t j = i;
		int num_recurse = 0;
		while (true) {
			bool cache_hit = false;
			if (i != j) {
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
					working_boundary.push(faces.nextFace);
				}
			}
			Cube pivot = get_pivot(working_boundary);
			if (pivot.index != NONE) {
				auto pair = pivot_column_index.find(pivot.index);
				if (pair != pivot_column_index.end()) {
					j = pair->second;
					num_recurse++;
					num_reduction_steps++;
					continue;
				} else {
					if (num_recurse >= config.min_recursion_to_cache) {
                        add_cache(i, working_boundary, cache);
						cached_column_idx.push(i);
						if (cached_column_idx.size() > config.cache_size) {
							cache.erase(cached_column_idx.front());
							cached_column_idx.pop();
						}
                    }
					pivot_column_index[pivot.index] = i;
					if (image) {
						pairs->push_back(Pair(pivot, ctr[i]));
					} else if (pivot.birth != ctr[i].birth) {
						pairs->push_back(Pair(pivot, ctr[i]));
					}
					break;
				}
			} else if (!image) { 
				shouldClear = true;
				ctr[i].index = NONE;
				break; 
			}
		}
	}
	if (config.verbose) {cout << " with " << num_reduction_steps << " reduction steps ... ";}
	if (shouldClear) {
		auto new_end = std::remove_if(ctr.begin(), ctr.end(),
								[](const Cube& c){ return c.index == NONE; });
		ctr.erase(new_end, ctr.end());
	}
}

Cube Dimension1::pop_pivot(CubeQue& column) {
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

Cube Dimension1::get_pivot(CubeQue& column) {
	Cube result = pop_pivot(column);
	if (result.index != NONE) {
		column.push(result);
	}
	return result;
}

void Dimension1::add_cache(uint32_t i, CubeQue &working_boundary, unordered_map<uint32_t, CubeQue> &cache) {
	CubeQue clean_wb;
	while (!working_boundary.empty()) {
		auto c = working_boundary.top();
		working_boundary.pop();
		if (!working_boundary.empty() && c.index == working_boundary.top().index) {
			working_boundary.pop();
		} else {
			clean_wb.push(c);
		}
	}
	cache.emplace(i, clean_wb);
}