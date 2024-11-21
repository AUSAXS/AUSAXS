#!/bin/bash

pdbprefix="output/waxsis/fitted/"
datprefix="data"
files=(
	SASDA45
	SASDA92
	SASDAW3
	SASDCQ8
	SASDDD3
	SASDE35
	SASDE45
	SASDE65
	SASDEL8
	SASDEL9
	SASDEM9
	SASDEY4
	SASDF86
	SASDFZ3
	SASDGD2
	SASDHP7
	SASDJF5
	SASDJG4
	SASDJG5
	SASDJP5
	SASDJQ4
	SASDJQ7
	SASDJU5
	SASDKG4
	SASDKH2
	SASDKP2
	SASDL82
	SASDLT3
	SASDME4
	SASDNQ3
	SASDNV5
	SASDNY2
	SASDPM2
	SASDT75
	SASDT85
	SASDT95
	SASDTT4
	SASDLQ6
	SASDMB5
	SASDMZ9
	SASDP39
	SASDPB9
	SASDQ59
	SASDQN4
)

size=${#files[@]}
pdbs=()
dats=()
all_exist=true
for ((i=0; i<${size}; i++)); do
    folder=$(find $pdbprefix -name *${files[i]})
    if [ -z "$folder" ]; then
        echo "Missing folder for ${files[i]}"
        all_exist=false
        continue
    fi

	pdb=$(find $folder -name prot+solvlayer_0.pdb)
	dat=$(find $datprefix -name ${files[i]}.dat)
	if [ -z "$pdb" ] || [ -z "$dat" ]; then
		echo "Missing file for ${files[i]}"
		all_exist=false
        continue
	fi
    pdbs+=("$pdb")
    dats+=("$dat")
done
if ! $all_exist; then
	exit 1
fi

for ((i=0; i<${size}; i++)); do
	folder_pdb=$(dirname ${pdbs[i]})
	stripped_pdb=${folder_pdb}/${files[i]}_stripped.pdb
    grep -E -v ' H |H$|[123]H|OW|HOH| H[ABCDEFGHZ][123]?|NA|CL| D[AGCT] |CYX' "${pdbs[i]}" | grep 'ATOM  ' > "${stripped_pdb}"
	build/bin/remove_fourth_column_data ${dats[i]}
    folder_dat=$(dirname ${dats[i]})
    stripped_dat=${folder_dat}/${files[i]}_stripped.dat
    cp ${stripped_dat} ${folder_pdb}

	ausaxs_output=$(build/bin/saxs_fitter ${stripped_pdb} ${stripped_dat} exv --model grid fit --atomic-debye-waller --exv-debye-waller --ignore-unknown-atom --output output/debye_waller_fitting/)
	echo "${ausaxs_output}"
	Ba=$(echo "$ausaxs_output" | grep "Ba" | awk -F '|' '{gsub (" ", "", $3); print $3}')
	Bx=$(echo "$ausaxs_output" | grep "Bx" | awk -F '|' '{gsub (" ", "", $3); print $3}')
	echo "${Ba} ${Bx}" >> output/debye_waller_fitting/fit_vals.dat
done
