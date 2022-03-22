#ifndef _READ_BUFFERS_H_
#define _READ_BUFFERS_H_

#ifdef __linux__
    #include <sys/types.h>
    #include <unistd.h>
#else
    #error "Not running on Linux"
#endif

#include "common.h"

#include <cassert>
#include <cstdio>

#include <sys/types.h>
#include <unistd.h>

template <typename K>
struct read_buffer {
    K buf[SECTORSZ / sizeof(K)];
    size_t file_pos = 0;
    size_t buf_sz = 0;
    uint64_t timestamp = 0;

    static_assert(SECTORSZ % sizeof(K) == 0);
    static_assert(sizeof(buf) == SECTORSZ);

    bool contains(size_t idx) {
        return idx >= file_pos && idx < file_pos + buf_sz;
    }
    K get(size_t idx) {
        assert(contains(idx));
        return buf[idx - file_pos];
    }
    void evict(int fd, size_t idx) {
        assert(!contains(idx));
        off_t read_pos = round_down(idx * sizeof(K), sizeof(buf));
        assert(lseek(fd, read_pos, SEEK_SET) == read_pos);
        buf_sz = read(fd, buf, sizeof(buf)) / sizeof(K);
        assert(buf_sz <= sizeof(buf) / sizeof(K));
        file_pos = read_pos / sizeof(K);
    }
};

template <typename K, size_t N>
struct read_buffers {
    uint64_t current_time = 0;
    size_t num_read_IO = 0;
    int fd;
    read_buffer<K> bufs[N];

    read_buffers(int file_desc) : fd(file_desc) {}

    K get(size_t idx) {
        do {
            for (size_t i = 0; i != N; ++i) {
                if (bufs[i].contains(idx)) {
                    bufs[i].timestamp = ++current_time;
                    return bufs[i].get(idx);
                }
            }
            evict_lru(idx);
            ++num_read_IO;
        } while (true);
    }

    void evict_lru(size_t idx) {
        uint64_t min_timestamp = bufs[0].timestamp;
        size_t min_idx = 0;
        for (size_t i = 1; i != N; ++i) {
            if (bufs[i].timestamp < min_timestamp) {
                min_timestamp = bufs[i].timestamp;
                min_idx = i;
            }
        }
        bufs[min_idx].evict(fd, idx);
        bufs[min_idx].timestamp = ++current_time;
    }
};

#endif

