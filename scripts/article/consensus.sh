#!/bin/bash

files=(
	data/consensus/SASDPS4
)

size=${#files[@]}
for ((i=0; i<${size}; i++)); do
	stem=$(basename ${files[i]})
	base=${files[i]}/${stem}
	eval "grep -v ' H ' ${base}.pdb | grep -v 'H$' | grep -v 'OW' | grep 'ATOM  ' > ${base}_stripped.pdb"
	eval "grep -v ' H ' ${base}_WAXSiS.pdb | grep -v 'H$' | grep -v 'OW' | grep 'ATOM  ' > ${base}_WAXSiS_stripped.pdb"
	build/bin/remove_fourth_column_data ${base}.dat

	# ausaxs
	build/bin/fit_all_exv ${base}_stripped.pdb ${base}_stripped.dat --no-exit-on-unknown-atom -o output/consensus/
	build/bin/fit_all_exv ${base}_WAXSiS_stripped.pdb ${base}_stripped.dat --no-exit-on-unknown-atom -o output/consensus/
	mkdir -p output/saxs_fitter/${stem}
	mkdir -p output/saxs_fitter/${stem}_WAXSiS

	# foxs
	rm -rf temp/foxs
	mkdir -p temp/foxs
	cp ${base}_stripped.pdb ${base}_WAXSiS_stripped.pdb ${base}_stripped.dat temp/foxs
	cd temp/foxs
	foxs $(basename ${base})_stripped.pdb $(basename ${base})_stripped.dat --max_q=1
	mv temp/foxs/*.fit output/consensus/${stem}/foxs.fit
	foxs $(basename ${base})_WAXSiS_stripped.pdb $(basename ${base})_WAXSiS_stripped.dat --max_q=1
	mv temp/foxs/*.fit output/consensus/${stem}_WAXSiS/foxs.fit
	cd ../..

	# pepsi
	rm -rf temp/pepsi
	mkdir -p temp/pepsi
	~/tools/Pepsi-SAXS/Pepsi-SAXS ${base}_stripped.pdb ${base}_stripped.dat -o "temp/pepsi/pepsi.fit" -ms 1
	mv temp/pepsi/pepsi.fit output/consensus/${stem}/pepsi.fit
	~/tools/Pepsi-SAXS/Pepsi-SAXS ${base}_WAXSiS_stripped.pdb ${base}_stripped.dat -o "temp/pepsi/pepsi.fit" -ms 1
	mv temp/pepsi/pepsi.fit output/consensus/${stem}_WAXSiS/pepsi.fit

	# crysol
	mkdir -p temp/crysol
	crysol_output=$(crysol ${base}_stripped.dat ${base}_stripped.pdb --prefix="temp/crysol/out" --constant --implicit-hydrogen=1 --smax=1)
    fit_file="output/consensus/${stem}/crysol.fit"
	mv temp/crysol/out.fit ${fit_file}
    if grep -q "Chi\^2: ******" "$fit_file"; then
        chi_square=$(echo "$crysol_output" | grep "Chi-square of fit" | awk -F': ' '{print $2}')
        sed -i "s/\*\*\*\*\*\*/${chi_square}/" "$fit_file"
    fi
	crysol_output=$(crysol ${base}_stripped.dat ${base}_WAXSiS_stripped.pdb --prefix="temp/crysol/out" --constant --implicit-hydrogen=1 --smax=1)
    fit_file="output/consensus/${stem}_WAXSiS/crysol.fit"
	mv temp/crysol/out.fit ${fit_file}
    if grep -q "Chi\^2: ******" "$fit_file"; then
        chi_square=$(echo "$crysol_output" | grep "Chi-square of fit" | awk -F': ' '{print $2}')
        sed -i "s/\*\*\*\*\*\*/${chi_square}/" "$fit_file"
    fi

	cp output/consensus/${stem}/*.fit output/saxs_fitter/$*
done