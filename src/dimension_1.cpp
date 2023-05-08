#include "dimension_1.h"
#include "boundary_enumerator.h"
#include <unordered_map>

Dimension1::Dimension1(CubicalGridComplex* _cgc, vector<WritePair>& _pairs) {
	cgc = _cgc;
	pairs = &_pairs;
}

void Dimension1::compute_pairs(vector<Cube>& ctr, bool image) {
	auto ctr_size = ctr.size();
	unordered_map<uint64_t, uint32_t> pivot_column_index;
	pivot_column_index.reserve(ctr_size);
	vector<Cube> face_entries;
	BoundaryEnumerator faces(cgc);
	Cube pivot;
	uint32_t j;
	bool clearing = false;
	for(uint32_t i = 0; i < ctr_size; i++) {
		CubeQue working_boundary;
		j = i;
		while (true) {
			face_entries.clear();
			faces.setBoundaryEnumerator(ctr[j]);
			while (faces.hasNextFace()) {
				face_entries.push_back(faces.nextFace);
			}
			for (auto f : face_entries) {
				working_boundary.push(f);
			}
			pivot = get_pivot(working_boundary);
			if (pivot.index != NONE) {
				auto pair = pivot_column_index.find(pivot.index);
				if (pair != pivot_column_index.end()) {
					j = pair->second;
					continue;
				} else {
					pivot_column_index[pivot.index] = i;
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
	if (clearing) {
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