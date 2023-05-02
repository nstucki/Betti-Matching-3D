#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <unordered_map>

using namespace std;

typedef float value_t;
typedef int64_t index_t;
typedef uint16_t coefficient_t;


template <typename S>
void print_vector(const vector<S>& vector,
                 string sep = " "){
    for (int i = 0; i < vector.size(); i++) {
        cout << vector[i] << sep;
    }
    cout << endl;
}


template <typename A, typename B>
void zip(
    const vector<A> &a, 
    const vector<B> &b, 
    vector<pair<A,B>> &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i], b[i]));
    }
}

template <typename A, typename B>
void unzip(
    const vector<pair<A, B>> &zipped, 
    vector<A> &a, 
    vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}

template <typename T>
bool contains(vector<T> vec, const T & elem){
    bool result = false;
    if( find(vec.begin(), vec.end(), elem) != vec.end() )
    {
        result = true;
    }
    return result;
}


struct filtration{
    vector<value_t> values_0;
    vector<index_t> indices_0;
    vector<vector<index_t>> indices_1;
    vector<vector<index_t>> indices_2;
    vector<vector<index_t>> indices_3;
};


class CubicalComplex{

    public:
    vector<value_t> Image_flattened;
    vector<index_t> shape;
    index_t dimensions;
    struct filtration filtration;

    CubicalComplex(vector<value_t> Image_flattened, vector<index_t> shape)
     : Image_flattened(Image_flattened), shape(shape), dimensions(shape.size()){};

     vector<index_t> get_coordinates(index_t idx){
        vector<index_t> coordinates;
        index_t coordinate;
        index_t dimension = dimensions-1;
        while(dimension >= 0){
            coordinate = idx % shape[dimension];
            coordinates.push_back(coordinate);
            idx -= coordinate;
            idx /= shape[dimension];
            dimension -= 1;
        };
        reverse(coordinates.begin() , coordinates.end());
        return coordinates;
     };

     index_t get_index(vector<index_t> coordinates){
        index_t idx = 0;
        index_t factor = 1;
        for (index_t i = dimensions-1; i>=0; i--){
            idx += coordinates[i]*factor;
            factor *= shape[i];
        }
        return idx;
     };


    index_t get_neighbor(index_t idx, vector<index_t> direction){
        vector<index_t> coordinates = get_coordinates(idx);
        std::transform(direction.begin(), direction.end(), coordinates.begin(), coordinates.begin(), std::plus<>());
        vector<index_t> zero(coordinates.size(),0);
        if (equal(begin(coordinates), end(coordinates),
                  begin(shape), end(shape),
                  [](int a, int b)->bool {return a <= b-1; })
            && equal(begin(coordinates), end(coordinates),
                  begin(zero), end(zero),
                  [](int a, int b)->bool {return a >= b; })){
            return get_index(coordinates);
        }
        else 
        {
            return -1;
        }
    };


     void initialize_filtration(){
        // initialize dim-0
        filtration.values_0 = Image_flattened;
        vector<index_t> indices_(Image_flattened.size());
        filtration.indices_0 = indices_;
        iota(begin(filtration.indices_0), end(filtration.indices_0), 0);
        vector<pair<value_t,index_t>> zipped;
        zip(filtration.values_0, filtration.indices_0, zipped);
        std::sort(std::begin(zipped), std::end(zipped), 
            [&](const auto& a, const auto& b){
            return a.first < b.first;
        });
        unzip(zipped, filtration.values_0, filtration.indices_0);
        
        //initialize dim-1 - dim-3
        index_t orientations[2] = {-1,1};
        vector<bool> inserted(filtration.indices_0.size(), false);
        vector<index_t> direction;
        index_t neighbor_1;
        index_t neighbor_2;
        index_t neighbor_3;
        index_t neighbor_12;
        index_t neighbor_13;
        index_t neighbor_23;
        index_t neighbor_123;
        vector<index_t> cell;

        for (int i = filtration.indices_0.size()-1; i >= 0; i--){
            index_t idx = filtration.indices_0[i];
            inserted[idx] = true;
            for (int j : orientations){
                direction = {j,0,0};
                neighbor_1 = get_neighbor(idx, direction);
                if (neighbor_1 != -1 && !inserted[neighbor_1]){
                    cell = {idx, neighbor_1};
                    filtration.indices_1.push_back(cell);
                }
                direction = {0,j,0};
                neighbor_2 = get_neighbor(idx, direction);
                if (neighbor_2 != -1 && !inserted[neighbor_2]){
                    cell = {idx, neighbor_2};
                    filtration.indices_1.push_back(cell);
                }
                direction = {0,0,j};
                neighbor_3 = get_neighbor(idx, direction);
                if (neighbor_3 != -1 && !inserted[neighbor_3]){
                    cell = {idx, neighbor_3};
                    filtration.indices_1.push_back(cell);
                }
            }
            for (int j : orientations){
                for (int k : orientations){    
                    direction = {j,0,0};
                    neighbor_1 = {get_neighbor(idx, direction)};
                    direction = {0,k,0};
                    neighbor_2 = {get_neighbor(idx, direction)};
                    direction = {j,k,0};
                    neighbor_12 = get_neighbor(idx, direction);
                    if (neighbor_12 != -1 && !inserted[neighbor_1] && !inserted[neighbor_2] && !inserted[neighbor_12]){
                        cell = {idx, neighbor_1, neighbor_2, neighbor_12};
                        filtration.indices_2.push_back(cell);
                    }
                    direction = {0,0,k};
                    neighbor_3 = {get_neighbor(idx, direction)};
                    direction = {j,0,k};
                    neighbor_13 = get_neighbor(idx, direction);
                    if (neighbor_13 != -1 && !inserted[neighbor_1] && !inserted[neighbor_3] && !inserted[neighbor_13]){
                        cell = {idx, neighbor_1, neighbor_3, neighbor_13};
                        filtration.indices_2.push_back(cell);
                    }
                    direction = {0,j,0};
                    neighbor_2 = {get_neighbor(idx, direction)};
                    direction = {0,j,k};
                    neighbor_23 = get_neighbor(idx, direction);
                    if (neighbor_23 != -1 && !inserted[neighbor_2] && !inserted[neighbor_3] && !inserted[neighbor_23]){
                        cell = {idx, neighbor_2, neighbor_3, neighbor_23};
                        filtration.indices_2.push_back(cell);
                    }
                }
            }
            for (int j : orientations){
                for (int k : orientations){
                    for (int l : orientations){  
                        direction = {j,0,0};
                        neighbor_1 = get_neighbor(idx, direction);
                        direction = {0,k,0};
                        neighbor_2 = get_neighbor(idx, direction);
                        direction = {0,0,l};
                        neighbor_3 = get_neighbor(idx, direction);
                        direction = {j,k,0};
                        neighbor_12 = get_neighbor(idx, direction);
                        direction = {j,0,l};
                        neighbor_13 = get_neighbor(idx, direction);
                        direction = {0,k,l};
                        neighbor_23 = get_neighbor(idx, direction);     
                        direction = {j,k,l};
                        neighbor_123 = get_neighbor(idx, direction);
                        if (neighbor_123 != -1 && !inserted[neighbor_1] && !inserted[neighbor_2] && !inserted[neighbor_3] && !inserted[neighbor_12] && !inserted[neighbor_13] && !inserted[neighbor_23] && !inserted[neighbor_123]){
                            cell = {idx, neighbor_1, neighbor_2, neighbor_3, neighbor_12, neighbor_13, neighbor_23, neighbor_123};
                            filtration.indices_3.push_back(cell);
                        }
                    }
                }
            }
        }
        reverse(filtration.indices_1.begin(), filtration.indices_1.end());
        reverse(filtration.indices_2.begin(), filtration.indices_2.end());
        reverse(filtration.indices_3.begin(), filtration.indices_3.end());
     }
};


int main()
{
    vector<vector<vector<value_t>>>  Image { {{{ 1, 2, 3},
                                               { 8, 9, 4},
                                               { 7, 6, 5}},
                                              {{21,22,23},
                                               {28,29,24},
                                               {27,26,25}},
                                              {{11,12,13},
                                               {18,19,14},
                                               {17,16,15}}} };

    vector<value_t> Image_flattened { {1,2,3,8,9,4,7,6,5,21,22,23,28,29,24,27,26,25,11,12,13,18,19,14,17,16,15} };
    vector<index_t> shape { {3,3,3} };

    //CubicalComplex CC(Image_flattened,shape);

    //CC.initialize_filtration();
    //print_vector(CC.filtration.values_0);
    //print_vector(CC.filtration.indices_0);
    //for (index_t i=0; i<CC.filtration.indices_3.size(); i++){
    //    for (index_t j=0; j<CC.filtration.indices_3[i].size(); j++){
    //        cout << CC.filtration.indices_3[i][j] << ' ';    
    //    }
    //    cout << '\n';
    //}
    return 0;
};