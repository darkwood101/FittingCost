#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <errno.h>

#include "pgm_index.h"

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
    constexpr int epsilon = 128;
    pgm::PGMIndex<int, epsilon> index(in_fd, out_fd);

    close(in_fd);
    close(out_fd);

    printf("%lu\n", index.these_buffers.num_read_IO);

    return EXIT_SUCCESS;
}

