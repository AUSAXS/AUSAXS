#!/bin/bash

pdb="data/consensus/SASDPP4/SASDPP4.pdb"
dat="data/consensus/SASDPP4/SASDPP4.dat"
iterations=100
warmup=3

hyperfine_cmd="hyperfine --warmup $warmup --runs $iterations"

crysol_bench() {
    eval "$hyperfine_cmd --command-name crysol --export-json temp/benchmarks/crysol.json 'crysol $dat $pdb --prefix=temp/crysol/out --constant'"
}

pepsi_bench() {
    eval "$hyperfine_cmd --command-name pepsi --export-json temp/benchmarks/pepsi.json '~/tools/Pepsi-SAXS/Pepsi-SAXS $pdb $dat -o temp/pepsi/pepsi.fit'"
}

foxs_bench() {
    cp "$dat" temp/foxs
    cp "$pdb" temp/foxs
    pdb_filename=$(basename "$pdb")
    dat_filename=$(basename "$dat")
    eval "$hyperfine_cmd --command-name foxs --export-json temp/benchmarks/foxs.json 'foxs $pdb_filename $dat_filename'"
}

ausaxs_bench() {
    eval "$hyperfine_cmd --command-name ausaxs --export-json temp/benchmarks/ausaxs.json 'build/bin/intensity_fitter $pdb $dat --output temp/ausaxs'"
}

"$@"
