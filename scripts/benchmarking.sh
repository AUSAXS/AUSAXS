#!/bin/bash

warmup=3
iterations=10
pdb="data/A2M_native/A2M_native.pdb"
dat="data/A2M_native/A2M_native.RSR"
hyperfine_cmd="hyperfine --warmup $warmup --runs $iterations"

crysol_bench() {
    eval "$hyperfine_cmd --command-name crysol --export-json test/benchmarks/crysol.json 'crysol $dat $pdb --prefix=temp/crysol/out --constant --implicit-hydrogen 1'"
}

pepsi_bench() {
    eval "$hyperfine_cmd --command-name pepsi --export-json test/benchmarks/pepsi.json '~/tools/Pepsi-SAXS/Pepsi-SAXS $pdb $dat -o temp/pepsi/pepsi.fit'"
}

foxs_bench() {
    cp "$dat" temp/foxs
    cp "$pdb" temp/foxs
    pdb_filename=$(basename "$pdb")
    dat_filename=$(basename "$dat")
    cd "temp/foxs"
    eval "$hyperfine_cmd --command-name foxs --export-json ../../test/benchmarks/foxs.json 'foxs $pdb_filename $dat_filename'"
}

ausaxs_bench() {
    eval "$hyperfine_cmd --command-name ausaxs --export-json test/benchmarks/ausaxs.json 'build/bin/intensity_fitter $pdb $dat --output temp/ausaxs/'"
}

ausaxs_bench_st() {
    eval "$hyperfine_cmd --command-name ausaxs-st --export-json test/benchmarks/ausaxs_st.json 'build/bin/intensity_fitter $pdb $dat -t 1 --output temp/ausaxs/'"
}

hyperfine_serial_cmd="hyperfine --warmup 0 --runs 5"
export map="data/emd_24889/emd_24889.map"
export mapdat="data/SASDJG5/SASDJG5.dat"
export min=17
export max=18
export steps=100
export step_size=$(bc -l <<< "($max - $min) / $steps")
crysol_serial_cmd() {
    for ((i = 0; i <= $steps; i++)); do
        val=$(bc -l <<< "$min + $i * $step_size")
        eval "em2dam $map $mapdat --threshold=$val --prefix temp/crysol/em2dam"
    done
}

crysol_serial_bench() {
    export -f crysol_serial_cmd
    eval "$hyperfine_serial_cmd --command-name crysol-serial --show-output --export-json test/benchmarks/crysol_serial.json --shell=bash 'crysol_serial_cmd'"
}

ausaxs_serial_bench() {
    eval "$hyperfine_serial_cmd --command-name ausaxs --export-json test/benchmarks/ausaxs_serial.json --show-output 'build/bin/em_bench $map $mapdat --levelmin $min --levelmax $max --max-iterations $steps'"
}

"$@"
