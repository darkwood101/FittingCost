#!/bin/bash

build () {
    make clean >/dev/null && make EPS=$1 >/dev/null
    if [ $? -ne 0 ]; then
        exit 1
    fi
}

outdir="zipf/"
headerline="data_size,read_IO,write_IO,ll_segmentcount,height"

mkdir -p "$outdir"

for eps in 32 64 128; do
    build $eps

    outcsv="$outdir/eps_$eps.csv"
    echo "$headerline" > "$outcsv"
    for size in 1000000 10000000 100000000; do
        ./generate_data -d zipf -f zipf.in -s $size -p 3.0 10000
        ./main zipf.in zipf.out >> "$outcsv"
        rm -f *.in
        rm -f *.out
    done
done
