#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <vector>
#include <random>
#include <algorithm>

// Distribution type
enum dist_type {
    UNIFORM = 0,
    NORMAL = 1,
    ZIPF = 2,
    ERROR = -1
};


// Parameters for uniform distribution:
// min and max value (inclusive)
struct uniform_params {
    int64_t min_;
    int64_t max_;
};

#define DEFAULT_UNIFORM_MIN     0
#define DEFAULT_UNIFORM_MAX     100


// Parameters for normal distribution:
// mean and standard deviation
struct normal_params {
    double mean_;
    double sd_;
};

#define DEFAULT_NORMAL_MEAN     0
#define DEFAULT_NORMAL_SD       10000


// Parameters for zipf distribution:
// alpha and N
struct zipf_params {
    double  alpha_;
    int64_t N_;
};

#define DEFAULT_ZIPF_ALPHA      0.5
#define DEFAULT_ZIPF_N          100

// Unified type `dist` for all distributions
union dist_params {
    struct uniform_params u_params_;
    struct normal_params  n_params_;
    struct zipf_params    z_params_;
};

struct dist {
    enum dist_type    type_;
    union dist_params params_;
};


// Program arguments
//   `filename_` - name of output file
//   `dist_`     - output distribution
//   `len_`      - number of sampled elements
struct prog_args {
    const char* filename_;
    struct dist dist_;
    size_t      len_;
};

// Parse arguments function
void usage(const char* prog_name);
int parse_args(int argc, char** argv, struct prog_args* args);
dist_type get_dist_type(const char* dist);
void set_default_params(struct dist* dist_);

int generate_data(prog_args* args);
void generate_uniform(std::vector<int>& data, struct uniform_params* u_params);
void generate_normal(std::vector<int>& data, struct normal_params* n_params);
void generate_zipf(std::vector<int>& data, struct zipf_params* z_params);


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
    data.resize(args->len_);

    FILE* fout = fopen(args->filename_, "wb");
    if (!fout) {
        fprintf(stderr, "Error opening file %s: %s\n",
                args->filename_, strerror(errno));
        return -1;
    }

    switch (args->dist_.type_) {
        case UNIFORM:
            generate_uniform(data, &args->dist_.params_.u_params_);
            break;
        case NORMAL:
            generate_normal(data, &args->dist_.params_.n_params_);
            break;
        case ZIPF:
            generate_zipf(data, &args->dist_.params_.z_params_);
            break;
        default:
            break; // won't ever happen
    }
    std::sort(data.begin(), data.end());

    int* data_ptr = data.data();
    if (fwrite(data_ptr, sizeof(*data_ptr), args->len_, fout) < args->len_) {
        fprintf(stderr, "Error writing to file %s: %s\n",
                args->filename_, strerror(errno));
        fclose(fout);
        return -1;
    }
    fflush(fout);
    fclose(fout);
    return 0;
}

void generate_uniform(std::vector<int>& data, struct uniform_params* u_params) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<int64_t> d{u_params->min_,
                                             u_params->max_};

    for (size_t i = 0; i != data.size(); ++i) {
        data[i] = d(gen);
    }
}

void generate_normal(std::vector<int>& data, struct normal_params* n_params) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{n_params->mean_,
                                 n_params->sd_};

    for (size_t i = 0; i != data.size(); ++i) {
        data[i] = std::round(d(gen));
    }
}

void generate_zipf(std::vector<int>& data, struct zipf_params* z_params) {
    // We'll need a uniform number generator
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<> d{0, 1};

    // Compute normalization constant
    double c = 0;
    for (int i = 1; i <= z_params->N_; ++i) {
        c += 1.0 / pow((double) i, z_params->alpha_);
    }
    c = 1.0 / c;

    double* sum_probs = (double*) malloc((z_params->N_ + 1) * sizeof(*sum_probs));
    sum_probs[0] = 0;
    for (int i = 1; i != z_params->N_; ++i) {
        sum_probs[i] = sum_probs[i - 1] + c / pow(i, z_params->alpha_);
    }

    for (size_t i = 0; i != data.size(); ++i) {
        double z = d(gen);

        int low = 1;
        int high = z_params->N_;
        int mid;
        int zipf_val = 0;

        do {
            mid = (low + high) / 2;
            if (sum_probs[mid] >= z && sum_probs[mid - 1] < z) {
                zipf_val = mid;
                break;
            } else if (sum_probs[mid] >= z) {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        } while (low <= high);

        data[i] = zipf_val;
    }
    free(sum_probs);
}

void usage(const char* prog_name) {
    fprintf(stderr, "Usage: %s -f filename -d distribution -s size [-p parameters]\n\n",
            prog_name);
    fprintf(stderr, "filename - The name of the output data file.\n");
    fprintf(stderr, "distribution - The output data distribution. "
                    "Possible values:\n"
                    "  uniform - uniform distribution\n"
                    "  normal - normal distribution\n"
                    "  zipf - zipf distribution");
    fprintf(stderr, "size - Number of data elements.\n");
    fprintf(stderr, "parameters - Optional distribution parameters.\n"
                    "These have to be set after distribution.\n"
                    "Possible values: \n"
                    "  MIN MAX - for uniform distribution\n"
                    "  MEAN SD - for normal distribution\n"
                    "  ALPHA N - for zipf distribution");
    fflush(stderr);
}


// Helper parsing macros
// There are 4 things that the parser needs to set:
//   (1) distribution
//   (2) file name
//   (3) length of sample
//   (4) distribution parameters
// These macros use bitwise arithmetic to identify which of these things
// have been set so far, and to mark newly set things

#define DIST_SET_BIT        1U
#define FILE_SET_BIT        2U
#define LEN_SET_BIT         3U
#define PARAM_SET_BIT       4U

#define IS_DIST_SET(x)         ((x) & (1 << DIST_SET_BIT))
#define IS_FILE_SET(x)         ((x) & (1 << FILE_SET_BIT))
#define IS_LEN_SET(x)          ((x) & (1 << LEN_SET_BIT))
#define IS_PARAM_SET(x)        ((x) & (1 << PARAM_SET_BIT))
#define IS_ALL_SET(x)          (IS_DIST_SET(x) && \
                                IS_FILE_SET(x) && \
                                IS_LEN_SET(x) &&  \
                                IS_PARAM_SET(x))

// Everything is mandatory besides parameters
#define IS_MANDATORY_SET(x)    (IS_DIST_SET(x) && \
                                IS_FILE_SET(x) && \
                                IS_LEN_SET(x))

#define SET_DIST(x)            ((x) |= (1 << DIST_SET_BIT))
#define SET_FILE(x)            ((x) |= (1 << FILE_SET_BIT))
#define SET_LEN(x)             ((x) |= (1 << LEN_SET_BIT))
#define SET_PARAM(x)           ((x) |= (1 << PARAM_SET_BIT))

int parse_args(int argc, char** argv, struct prog_args* args) {
    uint32_t set_bits = 0;

    for (int i = 1; i < argc; ++i) {
        // Set file name
        if (strcmp(argv[i], "-f") == 0) {
            if (i == argc - 1) {
                return -1;
            }
            args->filename_ = argv[++i];
            SET_FILE(set_bits);
        }
        // Set sample length
        else if (strcmp(argv[i], "-s") == 0) {
            if (i == argc - 1) {
                return -1;
            }
            args->len_ = atol(argv[++i]);
            SET_LEN(set_bits);
        }
        // Set distribution
        else if (strcmp(argv[i], "-d") == 0) {
            if (i == argc - 1) {
                return -1;
            }
            args->dist_.type_ = get_dist_type(argv[++i]);
            if (args->dist_.type_ == ERROR) {
                return -1;
            }
            SET_DIST(set_bits);
        }
        // Set distribution parameters
        else if (strcmp(argv[i], "-p") == 0) {
            // Can't set distribution parameters before the
            // distribution has been set!!!
            if (!IS_DIST_SET(set_bits)) {
                return -1;
            }
            switch (args->dist_.type_) {
                case UNIFORM:
                    if (i >= argc - 2) {
                        return -1;
                    }
                    args->dist_.params_.u_params_.min_ = atol(argv[++i]);
                    args->dist_.params_.u_params_.max_ = atol(argv[++i]);
                    break;
                case NORMAL:
                    if (i >= argc - 2) {
                        return -1;
                    }
                    args->dist_.params_.n_params_.mean_ = atol(argv[++i]);
                    args->dist_.params_.n_params_.sd_ = atol(argv[++i]);
                    break;
                case ZIPF:
                    if (i >= argc - 2) {
                        return -1;
                    }
                    args->dist_.params_.z_params_.alpha_ = atof(argv[++i]);
                    args->dist_.params_.z_params_.N_ = atol(argv[++i]);
                    break;
                default:
                    return -1;
            }
            SET_PARAM(set_bits);
        }
    }
    // If mandatory bits are not set, we can't proceed
    if (!IS_MANDATORY_SET(set_bits)) {
        printf("%u\n", set_bits);
        return -1;
    }
    // If mandatory bits are set, but not distribution parameter
    // bits are not set, we set the default values
    if (!IS_PARAM_SET(set_bits)) {
        set_default_params(&args->dist_);
    }

    return 0;
}

dist_type get_dist_type(const char* dist) {
    if (strcmp(dist, "uniform") == 0) {
        return UNIFORM;
    } else if (strcmp(dist, "normal") == 0) {
        return NORMAL;  
    } else if (strcmp(dist, "zipf") == 0) {
        return ZIPF;  
    } else {
        return ERROR;
    }
}

void set_default_params(struct dist* dist_) {
    switch (dist_->type_) {
        case UNIFORM:
            dist_->params_.u_params_.min_ = DEFAULT_UNIFORM_MIN;
            dist_->params_.u_params_.max_ = DEFAULT_UNIFORM_MAX;
            break;
        case NORMAL:
            dist_->params_.n_params_.mean_ = DEFAULT_NORMAL_MEAN;
            dist_->params_.n_params_.sd_ = DEFAULT_NORMAL_SD;
            break;
        case ZIPF:
            dist_->params_.z_params_.alpha_ = 1.0;
            dist_->params_.z_params_.N_ = 1000;
        default:
            break;
    }
}

