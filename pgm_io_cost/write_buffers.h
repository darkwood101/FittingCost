#ifndef _WRITE_BUFFERS_H
#define _WRITE_BUFFERS_H

#ifdef __linux__
    #include <sys/types.h>
    #include <unistd.h>
#else
    #error "Not running on Linux"
#endif

#include "common.h"

template <typename K>
struct write_buffer {
    K buf[SECTORSZ / sizeof(K)];
    size_t file_pos = 0;
    size_t buf_sz = 0;
    uint64_t timestamp = 0;

    static_assert(SECTORSZ % sizeof(K) == 0);

    bool contains(size_t idx) {

    }
};

#endif

