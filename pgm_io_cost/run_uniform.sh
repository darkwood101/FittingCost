#!/bin/bash

build () {
    make clean >/dev/null && make EPS=$1 >/dev/null
    if [ $? -ne 0 ]; then
        exit 1
    fi
}

outdir="uniform/"
headerline="data_size,read_IO,write_IO,ll_segmentcount,height"

mkdir -p "$outdir"

for eps in 32 64 128; do
    build $eps

    outcsv="$outdir/eps_$eps.csv"
    echo "$headerline" > "$outcsv"
    for size in 1000000 10000000 100000000; do
        ./generate_data -d uniform -f uniform.in -s $size -p 0 10000000
        ./main uniform.in uniform.out >> "$outcsv"
        rm -f *.in
        rm -f *.out
    done
done
