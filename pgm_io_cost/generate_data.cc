#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <vector>
#include <random>
#include <algorithm>

enum dist_type {
    UNIFORM = 0,
    UNIFORM_NOISE = 1,
    ERROR = -1
};

struct prog_args {
    const char* filename;
    dist_type   dist;
    size_t      sz;
};

void usage(const char* prog_name);
int parse_args(int argc, char** argv, prog_args* args);
dist_type get_dist_type(const char* dist);

int generate_data(prog_args* args);
void generate_uniform(std::vector<int>& data);
void generate_unoise(std::vector<int>& data);


int main(int argc, char** argv) {
    prog_args args;
    if (parse_args(argc, argv, &args) < 0) {
        usage(argv[0]);
        return 1;
    }

    if (generate_data(&args) < 0) {
        return 1;
    }

    return 0;
}

int generate_data(prog_args* args) {
    std::vector<int> data;
    data.resize(args->sz);

    FILE* fout = fopen(args->filename, "wb");
    if (!fout) {
        fprintf(stderr, "Error opening file %s: %s\n",
                args->filename, strerror(errno));
        return -1;
    }

    switch (args->dist) {
        case UNIFORM:
            generate_uniform(data);
            break;
        case UNIFORM_NOISE:
            generate_unoise(data);
            break;
        default:
            break; // won't ever happen
    }
    std::sort(data.begin(), data.end());

    int* data_ptr = data.data();
    if (fwrite(data_ptr, sizeof(*data_ptr), args->sz, fout) < args->sz) {
        fprintf(stderr, "Error writing to file %s: %s\n",
                args->filename, strerror(errno));
        fclose(fout);
        return -1;
    }
    fflush(fout);
    fclose(fout);
    return 0;
}

void generate_uniform(std::vector<int>& data) {
    for (size_t i = 0; i != data.size(); ++i) {
        data[i] = i;
    }
}

void generate_unoise(std::vector<int>& data) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, 100};

    for (size_t i = 0; i != data.size(); ++i) {
        data[i] = i + std::round(d(gen));
    }
}

void usage(const char* prog_name) {
    fprintf(stderr, "Usage: %s -f filename -d distribution -s size\n\n",
            prog_name);
    fprintf(stderr, "filename - The name of the output data file.\n");
    fprintf(stderr, "distribution - The output data distribution. "
                    "Possible values:\n"
                    "  u - uniform distribution\n"
                    "  unoise - uniform distribution with normal noise\n");
    fprintf(stderr, "size - Number of data elements.\n");
    fflush(stderr);
}

int parse_args(int argc, char** argv, prog_args* args) {
    if (argc != 7) {
        return -1;
    }

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-f") == 0) {
            args->filename = argv[++i];
        } else if (strcmp(argv[i], "-d") == 0) {
            args->dist = get_dist_type(argv[++i]);
            if (args->dist < 0) {
                return -1;
            }
        } else if (strcmp(argv[i], "-s") == 0) {
            args->sz = atol(argv[++i]);
        } else {
            return -1;
        }
    }
    return 0;
}

dist_type get_dist_type(const char* dist) {
    if (strcmp(dist, "u") == 0) {
        return UNIFORM;
    } else if (strcmp(dist, "unoise") == 0) {
        return UNIFORM_NOISE;
    } else {
        return ERROR;
    }
}

