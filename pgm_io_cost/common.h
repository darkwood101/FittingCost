#ifndef _COMMON_H_
#define _COMMON_H_

#include <cstdint>

#ifdef __linux__
    #include <sys/types.h>
    #include <unistd.h>
#else
    #error "Not running on Linux"
#endif

#define SECTORSZ 512

// Round down `x` to the nearest multiple of `mul`
static constexpr uint64_t round_down(uint64_t x, uint64_t mul) {
    return x - x % mul;
}

static_assert(round_down(512, 512) == 512);
static_assert(round_down(13, 512) == 0);
static_assert(round_down(8, 7) == 7);
static_assert(round_down(237, 23) == 230);

#endif

