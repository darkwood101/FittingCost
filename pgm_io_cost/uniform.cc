#include <vector>
#include <cstddef>

void init_data(std::vector<int>& data) {
    for (size_t i = 0; i != data.size(); ++i) {
        data[i] = i;
    }
}

