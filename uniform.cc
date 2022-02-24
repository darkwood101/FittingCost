#include <vector>
#include <cstddef>

// Make the distribution fully uniform:
// [0, 1, 2, 3, 4, ...]
void fill_vector(std::vector<int>& vec) {
    for (size_t i = 0; i != vec.size(); ++i) {
        vec[i] = i;
    }
}
