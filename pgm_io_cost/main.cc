#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <errno.h>

#include "pgm_index.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s datafile\n",
                argv[0]);
        return 1;
    }

    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "Error opening file %s: %s\n",
                argv[1], strerror(errno));
        return 1;
    }

    // Construct the PGM-index
    constexpr int epsilon = 128;
    pgm::PGMIndex<int, epsilon> index(fd);

    close(fd);

    printf("%lu\n", index.these_buffers.num_read_IO);

    return 0;
}

