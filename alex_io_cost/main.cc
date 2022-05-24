#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <errno.h>
#include <cmath>
#include <unistd.h>

#include "alex.h"

#ifndef EPS
#define EPS 128
#endif

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s input_datafile output_datafile\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    int in_fd = open(argv[1], O_RDONLY);
    if (in_fd < 0) {
        fprintf(stderr, "Error opening file %s: %s\n",
                argv[1], strerror(errno));
        return EXIT_FAILURE;
    }
    int out_fd = open(argv[2], O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
    if (out_fd < 0) {
        close(in_fd);
        fprintf(stderr, "Error opening file %s: %s\n",
                argv[2], strerror(errno));
        return EXIT_FAILURE;
    }

    size_t data_len = lseek(in_fd, 0, SEEK_END) / sizeof(int);
    assert(data_len % 4 == 0);
    lseek(in_fd, 0, SEEK_SET);

    int* data = static_cast<int*>(malloc(data_len * sizeof(*data)));
    if (!data) {
        close(in_fd);
        close(out_fd);
        fprintf(stderr, "Error allocating memory: %s\n",
                strerror(errno));
        return EXIT_FAILURE;
    }

    size_t n_read = read(in_fd, data, data_len * sizeof(int));
    assert(n_read == data_len * sizeof(int));

    std::pair<int, int>* pairs = static_cast<std::pair<int, int>*>(malloc(data_len * sizeof(*pairs)));
    if (!pairs) {
        close(in_fd);
        close(out_fd);
        free(data);
        fprintf(stderr, "Error allocating memory: %s\n",
                strerror(errno));
        return EXIT_FAILURE;
    }

    for (size_t i = 0; i != data_len; ++i) {
        pairs[i].first = i;
        pairs[i].second = data[i];
    }

    printf("loading %lu\n", data_len);

    alex::Alex<int, int> index;
    index.bulk_load(pairs, data_len);

    alex::Alex<int, int>::NodeIterator iter(&index);

    int* level_nums = static_cast<int*>(malloc(10 * sizeof(int)));
    memset(level_nums, 0, 10 * sizeof(int));
    if (!level_nums) {
        close(in_fd);
        close(out_fd);
        free(data);
        free(pairs);
        fprintf(stderr, "Error allocating memory: %s\n",
                strerror(errno));
        return EXIT_FAILURE;
    }

    for (alex::AlexNode<int, int>* node = iter.current(); node; node = iter.next()) {
        if (node->is_leaf_) {
            assert(node->level_ < 10);
            ++level_nums[node->level_];
        }
    }

    for (int i = 0; i != 10; ++i) {
        printf("Num of nodes with level %d: %d\n", i, level_nums[i]);
    }

    free(data);
    free(pairs);
    close(in_fd);
    close(out_fd);

    return EXIT_SUCCESS;
}

