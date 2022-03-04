#!/bin/bash

buildmodel () {
    make clean >/dev/null && make model_uniform_noise >/dev/null
    if [ $? -ne 0 ]; then
        exit 1
    fi
}

# Prepare output files
mkdir -p "uniform_noise/"
mkdir -p "uniform_noise/data/"
outcsv="uniform_noise/data/uniform_noise.csv"
echo "data_size,total_upper_accesses,total_lower_accesses,total_data_accesses" > "$outcsv"

buildmodel

for size in $(seq 10000 50000 1000000); do
    for t in $(seq 1 1 10); do
        ./model_uniform_noise $size >> "$outcsv"
    done
done
