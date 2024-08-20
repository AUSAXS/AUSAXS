#!/bin/bash

files=(
	data/SASDA45
	data/SASDA92
	data/SASDAW3
	data/SASDCQ8
	data/SASDDD3
	data/SASDE35
	data/SASDE45
	data/SASDE65
	data/SASDEL8
	data/SASDEL9
	data/SASDEM9
	data/SASDEY4
	data/SASDF86
	data/SASDFZ3
	data/SASDGD2
	data/SASDGX9
	data/SASDHP7
	data/SASDJ84
	data/SASDJF5
	data/SASDJG4
	data/SASDJG5
	data/SASDJP5
	data/SASDJQ4
	data/SASDJQ7
	data/SASDJU5
	data/SASDJY3
	data/SASDKG4
	data/SASDKH2
	data/SASDKP2
	data/SASDL82
	data/SASDLT3
	data/SASDME4
	data/SASDNQ3
	data/SASDNV5
	data/SASDNY2
	data/SASDPM2
	data/SASDT75
	data/SASDT85
	data/SASDT95
	data/SASDTT4
	data/SASDLQ6
	data/SASDMB5
	data/SASDMZ9
	data/SASDP39
	data/SASDPB9
	data/SASDQ59
	data/SASDQN4
)

size=${#files[@]}
all_exist=true
for ((i=0; i<${size}; i++)); do
	if [[ ! -d ${files[i]} ]]; then
		echo "Folder ${files[i]} does not exist!"
		all_exist=false
	fi
done
if ! $all_exist; then
	exit
fi

for ((i=0; i<${size}; i++)); do
	stem=$(basename ${files[i]})
	base=${files[i]}/${stem}
	stripped=${base}_stripped
	eval "grep -v ' H ' ${base}.pdb | grep -v 'H$' | grep -v 'OW' | grep 'ATOM  ' > ${stripped}.pdb"
	build/bin/remove_fourth_column_data ${base}.dat

	# ausaxs
	build/bin/fit_all_exv ${stripped}.pdb ${stripped}.dat --no-exit-on-unknown-atom
	mkdir -p output/saxs_fitter/${stem}

	# foxs
	rm -rf temp/foxs
	mkdir -p temp/foxs
	cp ${stripped}.pdb ${stripped}.dat temp/foxs
	cd temp/foxs
	foxs $(basename ${stripped}).pdb $(basename ${stripped}).dat --max_q=1
	cd ../..
	mv temp/foxs/*.fit output/fit_all_exv/${stem}/foxs.fit

	# pepsi
	rm -rf temp/pepsi
	mkdir -p temp/pepsi
	~/tools/Pepsi-SAXS/Pepsi-SAXS ${stripped}.pdb ${stripped}.dat -o "temp/pepsi/pepsi.fit" -ms 1
	mv temp/pepsi/pepsi.fit output/fit_all_exv/${stem}/pepsi.fit

	# crysol
	mkdir -p temp/crysol
	crysol_output=$(crysol ${stripped}.dat ${stripped}.pdb --prefix="temp/crysol/out" --constant --implicit-hydrogen=1 --smax=1)
    fit_file="output/fit_all_exv/${stem}/crysol.fit"
	mv temp/crysol/out.fit ${fit_file}

    # check if crysol is being stupid with the chi2 val
    if grep -q "Chi\^2: ******" "$fit_file"; then
        chi_square=$(echo "$crysol_output" | grep "Chi-square of fit" | awk -F': ' '{print $2}')
        sed -i "s/\*\*\*\*\*\*/${chi_square}/" "$fit_file"
    fi

	mkdir -p output/saxs_fitter/${stem}
	cp output/fit_all_exv/${stem}/*.fit output/saxs_fitter/$*
done
