# generate lists of files for easy use as dependencies
pymol := pymol
simprog := pdb2mrc
pdbfixer := ~/tools/conda/bin/pdbfixer

cmake_threads := 6

source := $(addprefix source/, $(shell find source/ -printf "%P "))
include := $(addprefix include/, $(shell find include/ -printf "%P "))

# all targets passes the options string to the executable
options :=


.SECONDARY:
.SECONDEXPANSION:
#################################################################################
###				PYTHON PLOTS				      ###
#################################################################################
# Plot a SAXS dataset along with any accompanying fits
plot_fits/%: scripts/compare_fit.py
	@ measurement=$$(find output/intensity_fitter/ -name "$*.dat"); \
	for f in $${measurement}; do\
		python3 $< $${f}; \
	done

plot_em/%: scripts/plot_fit.py
	@ measurement=$$(find output/em_fitter/ -name "$*.RSR" -or -name "$*.dat"); \
	echo $${measurement}; \
	for f in $${measurement}; do\
		folder=$$(dirname $${f}); \
		fit=$${folder}/fit.fit; \
		python3 $< $${f} $${fit} 5; \
	done

# run the plotting script on all files in a folder
plot/%: scripts/plot.py
	python3 $< $*

# run the plotting script on all files in a folder, using larger text fonts than usual
bigplot/%: scripts/plot.py
	python3 $< $* --big


#################################################################################
###				UTILITY					      ###
#################################################################################
# compile & show the docs
docs: 
	@ make -C build doc
	firefox build/docs/html/index.html 

.PHONY:
gui: build/source/gui/gui
	build/gui

# show the coverage of a single test. requires code to be compiled with the --coverage option.
coverage/%: test/%
	@ rm -r temp/coverage
	@ rm -r build/CMakeFiles/test.dir/test/*.gcda
	@ mkdir -p temp/coverage/
	gcovr --filter source/ --filter include/ --exclude-throw-branches --html-details temp/coverage/coverage.html
	firefox temp/coverage/coverage.html

# show the coverage of all tests. requires code to be compiled with the --coverage option.
coverage: tests
	@ rm -r temp/coverage
	@ rm -r build/CMakeFiles/test.dir/test/*.gcda
	@ mkdir -p temp/coverage/
	gcovr --filter source/ --filter include/ --exclude-throw-branches --html-details temp/coverage/coverage.html
	firefox temp/coverage/coverage.html


###################################################################################
###				EXECUTABLES					###
###################################################################################

gmx/%: build/bin/gmx
	$< $*

casein/%: build/bin/casein
	@ tomo=$$(find data/ -name "$*.rec");\
	$< $${tomo} ${options}

silica/%: build/bin/silica
	@ tomo=$$(find data/ -name "$*.uc");\
	$< $${tomo} ${options}
	make plot/output/silica/$*
#	valgrind --track-origins=yes --log-file="valgrind.txt" $< $${tomo} ${options}

crystal/%: build/bin/crystal_scattering
	@ grid=$$(find data/ -name "$*.uc" -or -name "$*.grid" -or -name "$*.pdb");\
	$< $${grid} ${options}
	make plot/output/crystal/$*

crystal_compare/%: build/bin/crystal_comparison
	@ grid=$$(find data/ -name "$*.uc" -or -name "$*.grid" -or -name "$*.pdb");\
	$< $${grid} ${options}
	make plot/output/crystal_compare/$*

# fix a pdb file
data/%_fixed.pdb: data/%.pdb
	@ $(pdbfixer) $< --replace-nonstandard --add-atoms=all --add-residues --output=$@ --verbose

# calculate the scattering from a pdb structure
scatter/%: build/bin/scattering
	@ structure=$$(find data/ -name "$*.pdb"); \
	$< $${structure} ${options}
	make plot/output/scattering/$*

# hydrate a structure and show it in pymol
hydrate/%: build/bin/new_hydration
	@ structure=$$(find data/ -name "$*.pdb"); \
	$< $${structure} output/$*.pdb ${options}
	$(pymol) output/$*.pdb -d "hide all; show spheres, hetatm; color orange, hetatm"

# show a structure in pymol
view/%: 
	@ file=$$(find data -name "$*.pdb"); \
	$(pymol) $${file} -r scripts/pymol.py

viewmap/%: 
	@ file=$$(find data -name "$*.map" -or -name "$*.ccp4" -or -name "$*.mrc"); \
	if [ $${file} ]; then \
		folder=$$(dirname $${file}); \
		pdb=$$(find $${folder} -name "*.pdb" -or -name "*.ent"); \
		$(pymol) $${file} -r scripts/pymol.py -d "isomesh mesh, $*, 3; color black, mesh; bg_color white"; \
	else \
		echo "File \"$*\" not found."; \
	fi
#		$(pymol) $${file} $${pdb} -r scripts/pymol.py -d "isomesh mesh, $*, 3; color black, mesh; bg_color white"; \
#	fi

res := 10
# show a simulated structure in pymol
simview/%:
	@ structure=$$(find data/ -name "$*.pdb"); \
	emmap=$$(find sim/ -name "$*_$(res).mrc" -or -name "$*_$(res).ccp4"); \
	$(pymol) $${structure} $${emmap} -d "isomesh mesh, $*_$(res), 1"

# calculate the histogram for a given structure
hist/%: build/bin/hist
	@ structure=$$(find data/ -name "$*.pdb"); \
	$< $${structure} output/hist/$*/ ${options}
	make plot/output/hist/$*/

# flip the axes of an EM map
order := ""
rotate/%: build/bin/rotate_map
	$< data/$* ${order}

# main executable. primarily used for debugging purposes
main/%: build/bin/main
	$< $*

# Inspect the header of an EM map
inspect/%: build/bin/inspect_map
	@ emmaps=$$(find data/ -name "$*.map" -or -name "$*.ccp4" -or -name "$*.mrc" -or -name "$*.rec"); \
	for emmap in $${emmaps}; do\
		echo "Opening " $${emmap} " ...";\
		sleep 1;\
		$< $${emmap};\
	done

# Fit an EM map to a SAXS measurement file.  
# The wildcard should be the name of a measurement file. All EM maps in the same folder will then be fitted to the measurement.
em_fit/%: build/bin/em_fitter
	@ measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	folder=$$(dirname $${measurement}); \
	emmaps=$$(cat $${folder}/maps.txt); \
	for map in $${emmaps}; do\
		path=$$(find data/$${map}/ -name "*.map" -or -name "*.ccp4" -or -name "*.mrc"); \
		echo "Fitting " $${path} "..."; \
		sleep 1; \
		$< $${path} $${measurement} ${options}; \
		make plot_em/$*; \
		make plot/output/em_fitter/$*; \
	done

# Fit both an EM map and a PDB file to a SAXS measurement. 
# The wildcard should be the name of a measurement file. All EM maps in the same folder will then be fitted to the measurement. 
em/%: build/bin/em
	@ structure=$$(find data/ -name "$*.pdb"); \
	measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	emmap=$$(find data/ -name "*$*.map" -or -name "*$*.ccp4"); \
	$< $${emmap} $${structure} $${measurement}

optimize_radius/%: build/source/scripts/optimize_radius
	$< data/$*.pdb output/

# Perform a rigid-body optimization of the input structure. 
# The wildcard should be the name of both a measurement file and an associated PDB structure file. 
rigidbody/%: build/bin/rigidbody
	@ structure=$$(find data/ -name "$*.pdb"); \
	measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	echo "$< $${structure} $${measurement}";\
	$< $${structure} $${measurement} ${options}
	make plot/output/rigidbody/$*
	make plot_fits/$*

# perform a fit with crysol
crysol/%: 
	@ mkdir -p temp/crysol/
	@ measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	folder=$$(dirname $${measurement}); \
	structure=$$(find $${folder}/ -name "*.pdb"); \
	crysol $${measurement} $${structure} --prefix="temp/crysol/out" --constant ${options}
	@ mv temp/crysol/out.fit output/intensity_fitter/$*/crysol.fit

# perform a fit with pepsi-saxs
pepsi/%:
	@ mkdir -p temp/pepsi/
	@ measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	folder=$$(dirname $${measurement}); \
	structure=$$(find $${folder}/ -name "*.pdb"); \
	~/tools/Pepsi-SAXS/Pepsi-SAXS $${structure} $${measurement} -o "temp/pepsi/pepsi.fit"
	@ mv temp/pepsi/pepsi.fit output/intensity_fitter/$*/pepsi.fit
	

# Perform a fit of a structure file to a measurement. 
# All structure files in the same location as the measurement will be fitted. 
intensity_fit/%: build/bin/intensity_fitter
	@ measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat" -or -name "$*.xvg"); \
	folder=$$(dirname $${measurement}); \
	structure=$$(find $${folder}/ -name "*.pdb"); \
	for pdb in $${structure}; do\
		echo "Fitting " $${pdb} " ...";\
		sleep 1;\
		$< $${measurement} $${pdb} ${options};\
		make plot/output/intensity_fitter/$*;\
		make plot_fits/$*;\
	done

# Check the consistency of the program. 
# The wildcard should be the name of an EM map. A number of SAXS measurements will be simulated from the map, and then fitted to it. 
consistency/%: build/bin/consistency
	@ measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	folder=$$(dirname $${measurement}); \
	emmap=$$(find $${folder}/ -name "*.map" -or -name "*.ccp4" -or -name "*.mrc"); \
	for map in $${emmap}; do\
		echo "Testing " $${map} " ...";\
		sleep 1;\
		$< $${map};\
	done

# usage: make fit_consistency/2epe res=10
# Check the consistency of the program. 
# The wildcard should be the name of both a measurement file and an associated PDB structure file. 
# A simulated EM map must be available. The resolution can be specified with the "res" argument. 
fit_consistency/%: build/bin/fit_consistency
	@ structure=$$(find data/ -name "$*.pdb"); \
	measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	emmap=$$(find sim/ -name "$*_${res}.ccp4" -or -name "$*_${res}.mrc"); \
	$< $${emmap} $${structure} $${measurement}

# usage: make fit_consistency2/2epe res=10
# Check the consitency of the program.
# The wildcard should be the name of a PDB structure file. 
# A simulated EM map must be present available with the given resolution. 
# A measurement will be simulated from the PDB structure, and fitted to the EM map. 
fit_consistency2/%: build/bin/fit_consistency2
	@ structure=$$(find data/ -name "$*.pdb"); \
	emmap=$$(find sim/ -name "$*_${res}.ccp4" -or -name "$*_${res}.mrc"); \
	echo "$< $${emmap} $${structure}"; \
	$< $${emmap} $${structure}

map_consistency/%: build/bin/em_pdb_fitter
	@ map=$$(find data -name "$*.ccp4" -or -name "$*.map" -or -name "$*.mrc"); \
	folder=$$(dirname $${map}); \
	pdb=$$(find $${folder} -name "*.ent"); \
	$< $${map} $${pdb} $${options}
	make plot/output/em_pdb_fitter/$*

# Rebin a SAXS measurement file. This will dramatically reduce the number of data points. 
# The wildcard should be the name of a SAXS measurement file. 
rebin/%: build/bin/rebin
	@ measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	$< $${measurement}

# Calculate a unit cell and write it to the file as a CRYST1 record. 
# The wildcard should be the name of a PDB structure file. 
unit_cell/%: build/bin/unit_cell
	@ structure=$$(find data/ -name "$*.pdb"); \
	$< $${structure}


####################################################################################
###			     SIMULATIONS					 ###
####################################################################################

# eval "$(~/tools/eman2/bin/conda shell.bash hook)"
simulate/%: 
	@ structure=$$(find data/ -name "$*.pdb"); \
	python3 ~/tools/eman2/bin/e2pdb2mrc.py $${structure} sim/$*_$(res).mrc res=$(res) het

#simulate/%: 
#	@ structure=$(shell find data/ -name "$*.pdb"); \
#	$(simprog) $${structure} sim/$*_$(res).mrc res=$(res) het center

simfit/%: build/bin/fit_consistency
	@ structure=$$(find data/ -name "$*.pdb"); \
	measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	$(simprog) $${structure} sim/$*_$(res).mrc res=$(res) het center; \
	echo "./fit_consistency $${structure} $${measurement} sim/$*_$(res).mrc\n"; \
	$< sim/$*_$(res).mrc $${structure} $${measurement}

old_simulate/%: 
	@ structure=$$(find data/ -name "$*.pdb"); \
	phenix.fmodel $${structure} high_resolution=$(res) generate_fake_p1_symmetry=True;\
	phenix.mtz2map mtz_file=$(*F).pdb.mtz pdb_file=$${structure} labels=FMODEL,PHIFMODEL output.prefix=$(*F);\
	rm $(*F).pdb.mtz;\
	mv $(*F)_fmodel.ccp4 sim/$(*F)_$(res).ccp4;\

stuff/%: build/bin/stuff
	@ structure=$$(find data/ -name "$*.pdb"); \
	$< $${structure}

####################################################################################
###				TESTS						 ###
####################################################################################
tags := ""
exclude_tags := "~[broken] ~[manual] ~[slow] ~[disable]"
test_files = $(addprefix test/, $(shell find test/ -wholename "*.cpp" -printf "%P "))

memtest/%: $$(shell find test/ -wholename "%.cpp")
	@ make -C build "test_$*" -j${cmake_threads}
	valgrind --track-origins=yes --log-file="valgrind.txt" build/test/bin/test_$* ~[slow] ~[broken] ${tags}

debug_tests: $(source) $(test_files)
	@ make -C build tests -j${cmake_threads}
	@ for test in $$(find build/test/bin/test_*); do\
		$${test} $(exclude_tags);\
	done

debug:
	echo $(test_files)

tests: $(source_files) $(test_files)
	@ make -C build tests -j${cmake_threads}
	@ mkdir -p build/test/reports
	@ for test in $$(find build/test/bin/test_*); do\
		$${test} $(exclude_tags) --reporter junit --out build/test/reports/$$(basename $${test}).xml;\
	done

test/%: $$(shell find test/ -wholename "%.cpp")
	@ make -C build "test_$*" -j${cmake_threads}
	build/test/bin/test_$* ~[slow] ~[broken] ${tags}

# special build target for our tests since they obviously depend on themselves, which is not included in $(source_files)
build/test/%: $$(shell find source/ -print) $(shell find test -name *%.cpp) build/Makefile
	@ cmake --build build --target $* -j${cmake_threads}


####################################################################################
###				BUILD						 ###
####################################################################################
.PHONY: build winbuild
build: 
	@ mkdir -p build; 
	@ cd build; cmake ../ $(ARGS)

buildstatic: 
	@ mkdir -p build; 
	@ cd build; cmake -DBUILD_SHARED_LIBS=OFF ../ 

winbuild: 
	@ mkdir -p winbuild;
	@ cd winbuild; cmake -DCMAKE_TOOLCHAIN_FILE=cmake/TC-mingw.cmake -DBUILD_SHARED_LIBS=OFF ../ 

winbuild/bin/%: $(source) $(include) executable/%.cpp
	@ cmake --build winbuild/ --target $(*F) -j${cmake_threads}

winbuild/%: $(source) $(include)
	@ cmake --build winbuild/ --target $(*F) -j${cmake_threads} 

build/bin/%: $(source) $(include) executable/%.cpp
	@ cmake --build build/ --target $(*F) -j${cmake_threads} 

build/%: $(source) $(include)
	@ cmake --build build/ --target $(*F) -j${cmake_threads}
	
build/Makefile: CMakeLists.txt
	@ mkdir -p build
	@ cd build; cmake ../
	
clean/build: 
	@ rmdir -f build

clean/winbuild:
	@ rmdir -f winbuild
