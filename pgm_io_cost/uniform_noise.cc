#include <vector>
#include <cstddef>
#include <random>

// Make the distribution almost fully uniform:
// [0, 1, 2, 3, ...] but add a bit of Gaussian noise to
// each entry. The noise has mean 0 and standard deviation 100
void init_data(std::vector<int>& data) {
    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution<> d{0, 100};

    for (size_t i = 0; i != data.size(); ++i) {
        data[i] = i + std::round(d(gen));
    }
}

