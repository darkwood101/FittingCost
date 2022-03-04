#!/bin/bash

buildmodel () {
    make clean >/dev/null && make model_uniform >/dev/null
    if [ $? -ne 0 ]; then
        exit 1
    fi
}

# Prepare output files
mkdir -p "uniform/"
mkdir -p "uniform/data/"
outcsv="uniform/data/uniform.csv"
echo "data_size,total_upper_accesses,total_lower_accesses,total_data_accesses" > "$outcsv"

buildmodel

for size in $(seq 1000 5000 100000); do
    ./model_uniform $size >> "$outcsv"
done
