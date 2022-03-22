#ifndef _WRITE_BUFFERS_H
#define _WRITE_BUFFERS_H

#ifdef __linux__
    #include <sys/types.h>
    #include <unistd.h>
#else
    #error "Not running on Linux"
#endif

#include "common.h"

#include <cassert>

template <typename K>
struct write_buffer {
    K buf[SECTORSZ / sizeof(K)];
    size_t buf_sz = 0;

    static_assert(SECTORSZ % sizeof(K) == 0);
    static_assert(sizeof(buf) == SECTORSZ);

    bool full() {
        return buf_sz == sizeof(buf) / sizeof(K);
    }
    void push_back(int fd, K& el) {
        if (full()) {
            evict(fd);
        }
        buf[buf_sz++] = el;
    }
    void evict(int fd) {
        size_t file_pos = lseek(fd, 0, SEEK_END) / sizeof(K);
        off_t write_pos = file_pos * sizeof(K);
        assert(lseek(fd, write_pos, SEEK_SET) == write_pos);
        write(fd, buf, buf_sz * sizeof(K));
        buf_sz = 0;
    }
};

#endif

