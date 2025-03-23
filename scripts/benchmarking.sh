#!/bin/bash

index=14;
iterations=10;
warmup=3
#	Index	Name		Count
#	0		SASDPT4		1960
#	1		SASDPP4		2000
#	2		SASDPS4		2835
#	3		SASDE35		2955
#	4		SASDQ59		3540
#	5		SASDJG5		4734
#	6		SASDME4		5458
#	7		SASDPB9		7715
#	8		urateox		9436
#	9		SASDDD3		12144
#	10		SASDPR4		12376
#	11		SASDEL9		16640
#	12		SASDPQ4		18792
#	13		SASDA92		?
#	14		A2M_native	43652

pdbs=(
	data/consensus/SASDPT4/SASDPT4_stripped.pdb
	data/consensus/SASDPP4/SASDPP4_stripped.pdb
	data/consensus/SASDPS4/SASDPS4_stripped.pdb
	data/SASDE35/SASDE35_stripped.pdb
	data/rigidbody/SASDQ59/SASDQ59_stripped.pdb
	data/SASDJG5/SASDJG5_stripped.pdb
	data/SASDME4/SASDME4_stripped.pdb
	data/rigidbody/SASDPB9/SASDPB9_stripped.pdb
	data/rigidbody/urateox/urateox_stripped.pdb
	data/SASDDD3/SASDDD3_stripped.pdb
	data/consensus/SASDPR4/SASDPR4_stripped.pdb
	data/SASDEL9/SASDEL9_stripped.pdb
	data/consensus/SASDPQ4/SASDPQ4_stripped.pdb
	data/SASDA92/SASDA92_stripped.pdb
	data/A2M_native/A2M_native_stripped.pdb
)

dats=(
	data/consensus/SASDPT4/SASDPT4.dat
	data/consensus/SASDPP4/SASDPP4.dat
	data/consensus/SASDPS4/SASDPS4.dat
	data/SASDE35/SASDE35.dat
	data/rigidbody/SASDQ59/SASDQ59.dat
	data/SASDJG5/SASDJG5.dat
	data/SASDME4/SASDME4.dat
	data/rigidbody/SASDPB9/SASDPB9.dat
	data/rigidbody/urateox/urateox.dat
	data/SASDDD3/SASDDD3.dat
	data/consensus/SASDPR4/SASDPR4.dat
	data/SASDEL9/SASDEL9.dat
	data/consensus/SASDPQ4/SASDPQ4.dat
	data/SASDA92/SASDA92.dat
	data/A2M_native/A2M_native.dat
)

pathnames=(
	SASDPT4
	SASDPP4
	SASDPS4
	SASDE35
	SASDQ59
	SASDJG5
	SASDME4
	SASDPB9
	urateox
	SASDDD3
	SASDPR4
	SASDEL9
	SASDPQ4
	SASDA92
	A2M_native
)

pdb=${pdbs[index]}
dat=${dats[index]}
pathname=${pathnames[index]}

generate_stripped_all() {
	for ((i=0; i<14; i++)); do
		eval "grep -v ' H ' $(dirname "${pdbs[i]}")/${pathnames[i]}.pdb | grep -v 'H$' | grep -v 'OW' | grep 'ATOM  ' > $(dirname "${pdbs[i]}")/${pathnames[i]}_stripped.pdb"
	done
}

generate_stripped() {
	eval "grep -v ' H ' $(dirname "${pdbs[i]}")/${pathname}.pdb | grep -v 'H$' | grep -v 'OW' | grep 'ATOM  ' > $(dirname "${pdb}")/${pathname}_stripped.pdb"
}

hyperfine_cmd="hyperfine --warmup $warmup --time-unit second"

crysol_bench_all() {
	for ((i=0;i<14;i++)) do
		mkdir -p "tests/benchmarks/${pathnames[i]}"
		mkdir -p "temp/crysol/${pathnames[i]}"
		echo "temp/crysol/${pathnames[i]}"
		cp "${dats[i]}" "temp/crysol/${pathnames[i]}"
		cp "${pdbs[i]}" "temp/crysol/${pathnames[i]}"
		pdb_filename=$(basename "${pdbs[i]}")
		dat_filename=$(basename "${dats[i]}")
		cd "temp/crysol/${pathnames[i]}"
	    eval "$hyperfine_cmd --command-name crysol --export-json ../../../tests/benchmarks/${pathnames[i]}/crysol.json 'crysol $dat_filename $pdb_filename --constant --implicit-hydrogen 1'"
	    cd "../../.."
	done
}

crysol_bench() {
	mkdir -p "tests/benchmarks/$pathname"
	mkdir -p "temp/crysol/$pathname"
	cp "$dat" "temp/crysol/$pathname"
	cp "$pdb" "temp/crysol/$pathname"
	pdb_filename=$(basename "$pdb")
	dat_filename=$(basename "$dat")
	cd "temp/crysol/$pathname"
	eval "$hyperfine_cmd --command-name crysol --export-json ../../../tests/benchmarks/$pathname/crysol.json 'crysol $dat_filename $pdb_filename --constant --implicit-hydrogen 1'"
	cd "../../.."
}

pepsi_bench_all() {
	for ((i=0;i<14;i++))
	do
		mkdir -p "tests/benchmarks/${pathnames[i]}"
	    eval "$hyperfine_cmd --command-name pepsi-${pathnames[i]} --export-json tests/benchmarks/${pathnames[i]}/pepsi.json '~/tools/Pepsi-SAXS/Pepsi-SAXS ${pdbs[i]} ${dats[i]} -o temp/pepsi/${pathnames[i]}.fit'"
	done
}

pepsi_bench() {
    eval "$hyperfine_cmd --runs $iterations --command-name pepsi --export-json tests/benchmarks/$pathname/pepsi.json '~/tools/Pepsi-SAXS/Pepsi-SAXS $pdb $dat -o temp/pepsi/$pathname.fit'"
}

foxs_bench_all() {
    mkdir -p "temp/foxs"
    for ((i=0;i<14;i++)) do
        mkdir -p "tests/benchmarks/${pathnames[i]}"
        mkdir -p "temp/foxs/${pathnames[i]}"
        cp "${dats[i]}" "temp/foxs/${pathnames[i]}"
        cp "${pdbs[i]}" "temp/foxs/${pathnames[i]}"
        pdb_filename=$(basename "${pdbs[i]}")
        dat_filename=$(basename "${dats[i]}")
        cd "temp/foxs/${pathnames[i]}"
        eval "$hyperfine_cmd --command-name foxs-${pathnames[i]} --export-json ../../../tests/benchmarks/${pathnames[i]}/foxs.json 'foxs $pdb_filename $dat_filename'"
        cd "../../.."
    done
}

foxs_bench() {
    mkdir -p "temp/foxs"
    cp "$dat" "temp/foxs"
    cp "$pdb" "temp/foxs"
    pdb_filename=$(basename "$pdb")
    dat_filename=$(basename "$dat")
    cd "temp/foxs"
    eval "$hyperfine_cmd --show-output --runs $iterations --command-name foxs --export-json ../../tests/benchmarks/$pathname/foxs.json 'foxs $pdb_filename $dat_filename'"
}

ausaxs_bench_simple_all() {
	for ((i=0;i<14;i++)) do
		mkdir -p tests/benchmarks/${pathnames[i]}
	    eval "$hyperfine_cmd --command-name ausaxs-${pathnames[i]} --export-json tests/benchmarks/${pathnames[i]}/ausaxs_simple.json 'build/bin/saxs_fitter ${pdbs[i]} ${dats[i]} --output temp/ausaxs/ --ignore-unknown-atom'"
	done
}

ausaxs_bench_simple() {
    mkdir -p tests/benchmarks/$pathname
    eval "$hyperfine_cmd --runs $iterations --show-output --command-name ausaxs-simple --export-json tests/benchmarks/$pathname/ausaxs_simple.json 'build/bin/saxs_fitter $pdb $dat --output temp/ausaxs/ --ignore-unknown-atom'"
}

ausaxs_bench_fraser_all() {
	for ((i=0;i<14;i++)) do
		mkdir -p tests/benchmarks/${pathnames[i]}
	    eval "$hyperfine_cmd --command-name ausaxs-${pathnames[i]} --export-json tests/benchmarks/${pathnames[i]}/ausaxs_fraser.json 'build/bin/saxs_fitter ${pdbs[i]} ${dats[i]} --output temp/ausaxs/ exv --model fraser --fit --ignore-unknown-atom'"
	done
}

ausaxs_bench_fraser() {
    mkdir -p tests/benchmarks/$pathname
    eval "$hyperfine_cmd --runs $iterations --show-output --command-name ausaxs-fraser --export-json tests/benchmarks/$pathname/ausaxs_fraser.json 'build/bin/saxs_fitter $pdb $dat --output temp/ausaxs/ exv --model fraser --fit --ignore-unknown-atom'"
}

ausaxs_bench_grid_all() {
	for ((i=0;i<14;i++)) do
		mkdir -p tests/benchmarks/${pathnames[i]}
	    eval "$hyperfine_cmd --command-name ausaxs-${pathnames[i]} --export-json tests/benchmarks/${pathnames[i]}/ausaxs_grid.json 'build/bin/saxs_fitter ${pdbs[i]} ${dats[i]} --output temp/ausaxs/ exv --model grid --fit --ignore-unknown-atom'"
	done
}

ausaxs_bench_grid() {
    mkdir -p tests/benchmarks/$pathname
    eval "$hyperfine_cmd --runs $iterations --show-output --command-name ausaxs-grid --export-json tests/benchmarks/$pathname/ausaxs_grid.json 'build/bin/saxs_fitter $pdb $dat --output temp/ausaxs/ exv --model grid --fit --ignore-unknown-atom'"
}

ausaxs_bench_st_all() {
	for ((i=0;i<14;i++)) do
		mkdir -p tests/benchmarks/${pathnames[i]}
	    eval "$hyperfine_cmd --command-name ausaxs-${pathnames[i]} --export-json tests/benchmarks/${pathnames[i]}/ausaxs_st.json 'build/bin/saxs_fitter ${pdbs[i]} ${dats[i]} -t 1 --output temp/ausaxs/ --ignore-unknown-atom'"
	done
}

ausaxs_bench_st() {
    eval "$hyperfine_cmd --runs $iterations --command-name ausaxs-st --export-json tests/benchmarks/$pathname/ausaxs_st.json 'build/bin/saxs_fitter $pdb $dat -t 1 --output temp/ausaxs/ --ignore-unknown-atom'"
}

hyperfine_serial_cmd="hyperfine --warmup 0 --time-unit second"
export map="data/emd_24889/emd_24889.map"
export mapdat="data/SASDJG5/SASDJG5.dat"
#pathname="SASDJG5_serial"
export min=4
export max=5
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
    eval "$hyperfine_serial_cmd --command-name crysol-serial --export-json tests/benchmarks/$pathname/crysol_serial.json --shell=bash 'crysol_serial_cmd'"
}

ausaxs_serial_bench_all() {
	mkdir -p "tests/benchmarks/$pathname"
	for ((i=17;i>=2;i--)); do
		min=$i
		max=$((i+1))
	    eval "$hyperfine_serial_cmd --command-name ausaxs-serial --export-json tests/benchmarks/$pathname/ausaxs_serial_${min}_${max}.json --show-output 'build/bin/em_bench $map $mapdat --levelmin $min --levelmax $max --max-iterations $steps'"
	done
}

ausaxs_serial_bench() {
    eval "$hyperfine_serial_cmd --command-name ausaxs-serial --export-json tests/benchmarks/$pathname/ausaxs_serial.json --show-output 'build/bin/em_bench $map $mapdat --levelmin $min --levelmax $max --max-iterations $steps'"
}

"$@"
