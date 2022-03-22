#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <errno.h>
#include <cmath>

#include "pgm_index.h"

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

    // Construct the PGM-index
    constexpr int epsilon = EPS;
    // try 32
    pgm::PGMIndex<int, epsilon> index(in_fd, out_fd);

    close(in_fd);
    close(out_fd);

    size_t num_write_IO = std::ceil((float) index.size_in_bytes() / SECTORSZ);

    printf("%lu,%lu,%lu,%lu,%lu\n",
           index.n * sizeof(int),
           index.these_buffers.num_read_IO,
           num_write_IO,
           index.segments_count(),
           index.height());

    return EXIT_SUCCESS;
}

