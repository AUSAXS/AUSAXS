# generate lists of files for easy use as dependencies
pymol := pymol
simprog := pdb2mrc

cmake_threads := 6

source := $(addprefix source/, $(shell find source/ -printf "%P "))
include := $(addprefix include/, $(shell find include/ -printf "%P "))

# all targets passes the options string to the executable
options :=


.SECONDARY:
#################################################################################
###				PYTHON PLOTS				      ###
#################################################################################
# Plot a SAXS dataset along with any accompanying fits
plot_fits/%: scripts/compare_fit.py
	@ measurement=$(shell find figures/intensity_fitter/ -name "$*.dat"); \
	python3 $< $${measurement}

plot_em/%: scripts/plot_fit.py
	@ measurement=$$(find figures/em_fitter/ -name "$*.RSR" -or -name "$*.dat"); \
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

# calculate the scattering from a pdb structure
scatter/%: build/executable/scattering
	@ structure=$(shell find data/ -name "$*.pdb"); \
	$< $${structure} ${options}
	make plot/figures/scattering/$*

# hydrate a structure and show it in pymol
hydrate/%: build/executable/new_hydration
	@ structure=$(shell find data/ -name "$*.pdb"); \
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
hist/%: build/executable/hist
	@structure=$(shell find data/ -name "$*.pdb"); \
	$< $${structure} figures/hist/$*/ ${options}
	make plot/figures/hist/$*/

# flip the axes of an EM map
order := ""
rotate/%: build/executable/rotate_map
	$< data/$* ${order}

# main executable. primarily used for debugging purposes
main/%: build/executable/main
	$< $*

# Inspect the header of an EM map
inspect/%: build/executable/inspect_map
	@ emmaps=$$(find data/ -name "$*.map" -or -name "$*.ccp4" -or -name "$*.mrc"); \
	for emmap in $${emmaps}; do\
		echo "Opening " $${emmap} " ...";\
		sleep 1;\
		$< $${emmap};\
	done

# Fit an EM map to a SAXS measurement file.  
# The wildcard should be the name of a measurement file. All EM maps in the same folder will then be fitted to the measurement.
em_fit/%: build/executable/em_fitter
	@ measurement=$$(find data/ -name "$*.RSR" -or -name "$*.dat"); \
	folder=$$(dirname $${measurement}); \
	emmaps=$$(cat $${folder}/maps.txt); \
	for map in $${emmaps}; do\
		path=$$(find data/$${map}/ -name "*.map" -or -name "*.ccp4" -or -name "*.mrc"); \
		echo "Fitting " $${path} "..."; \
		sleep 1; \
		$< $${path} $${measurement} ${options}; \
		make plot_em/$*; \
		make plot/figures/em_fitter/$*; \
	done

# Fit both an EM map and a PDB file to a SAXS measurement. 
# The wildcard should be the name of a measurement file. All EM maps in the same folder will then be fitted to the measurement. 
em/%: build/executable/em
	@ structure=$(shell find data/ -name "$*.pdb"); \
	measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
	emmap=$(shell find data/ -name "*$*.map" -or -name "*$*.ccp4"); \
	$< $${emmap} $${structure} $${measurement}

optimize_radius/%: build/source/scripts/optimize_radius
	$< data/$*.pdb figures/

# Perform a rigid-body optimization of the input structure. 
# The wildcard should be the name of both a measurement file and an associated PDB structure file. 
rigidbody/%: build/executable/rigidbody
	@ structure=$(shell find data/ -name "$*.pdb"); \
	measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
	$< $${structure} $${measurement} figures/ ${options}
	make plot/output
	make plot_fits/$*

# perform a fit with crysol
crysol/%: 
	@ mkdir -p temp/crysol/
	@ measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
	folder=$$(dirname $${measurement}); \
	structure=$$(find $${folder}/ -name "*.pdb"); \
	crysol $${measurement} $${structure} --prefix="temp/crysol/out" --constant ${options}
	@ mv temp/crysol/out.fit figures/intensity_fitter/$*/crysol.fit

# perform a fit with pepsi-saxs
pepsi/%:
	@ mkdir -p temp/pepsi/
	@ measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
	folder=$$(dirname $${measurement}); \
	structure=$$(find $${folder}/ -name "*.pdb"); \
	~/tools/Pepsi-SAXS/Pepsi-SAXS $${structure} $${measurement} -o "temp/pepsi/pepsi.fit"
	@ mv temp/pepsi/pepsi.fit figures/intensity_fitter/$*/pepsi.fit
	

# Perform a fit of a structure file to a measurement. 
# All structure files in the same location as the measurement will be fitted. 
intensity_fit/%: build/executable/intensity_fitter
	@ measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat" -or -name "$*.xvg"); \
	folder=$$(dirname $${measurement}); \
	structure=$$(find $${folder}/ -name "*.pdb"); \
	for pdb in $${structure}; do\
		echo "Fitting " $${pdb} " ...";\
		sleep 1;\
		$< $${measurement} $${pdb} ${options};\
		make plot/figures/intensity_fitter/$*;\
		make plot_fits/$*;\
	done

# Check the consistency of the program. 
# The wildcard should be the name of an EM map. A number of SAXS measurements will be simulated from the map, and then fitted to it. 
consistency/%: build/executable/consistency
	@ measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
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
fit_consistency/%: build/executable/fit_consistency
	@ structure=$(shell find data/ -name "$*.pdb"); \
	measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
	emmap=$(shell find sim/ -name "$*_${res}.ccp4" -or -name "$*_${res}.mrc"); \
	$< $${emmap} $${structure} $${measurement}

# usage: make fit_consistency2/2epe res=10
# Check the consitency of the program.
# The wildcard should be the name of a PDB structure file. 
# A simulated EM map must be present available with the given resolution. 
# A measurement will be simulated from the PDB structure, and fitted to the EM map. 
fit_consistency2/%: build/executable/fit_consistency2
	@ structure=$(shell find data/ -name "$*.pdb"); \
	emmap=$$(find sim/ -name "$*_${res}.ccp4" -or -name "$*_${res}.mrc"); \
	echo "$< $${emmap} $${structure}"; \
	$< $${emmap} $${structure}

map_consistency/%: build/executable/em_pdb_fitter
	@ map=$$(find data -name "$*.ccp4" -or -name "$*.map" -or -name "$*.mrc"); \
	folder=$$(dirname $${map}); \
	pdb=$$(find $${folder} -name "*.ent"); \
	$< $${map} $${pdb} $${options}
	make plot/figures/em_pdb_fitter/$*

# Rebin a SAXS measurement file. This will dramatically reduce the number of data points. 
# The wildcard should be the name of a SAXS measurement file. 
rebin/%: build/executable/rebin
	@ measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
	$< $${measurement}

# Calculate a unit cell and write it to the file as a CRYST1 record. 
# The wildcard should be the name of a PDB structure file. 
unit_cell/%: build/executable/unit_cell
	@ structure=$(shell find data/ -name "$*.pdb"); \
	$< $${structure}


####################################################################################
###			     SIMULATIONS					 ###
####################################################################################

# eval "$(~/tools/eman2/bin/conda shell.bash hook)"
simulate/%: 
	@ structure=$(shell find data/ -name "$*.pdb"); \
	python3 ~/tools/eman2/bin/e2pdb2mrc.py $${structure} sim/$*_$(res).mrc res=$(res) het

#simulate/%: 
#	@ structure=$(shell find data/ -name "$*.pdb"); \
#	$(simprog) $${structure} sim/$*_$(res).mrc res=$(res) het center

simfit/%: build/executable/fit_consistency
	@ structure=$(shell find data/ -name "$*.pdb"); \
	measurement=$(shell find data/ -name "$*.RSR" -or -name "$*.dat"); \
	$(simprog) $${structure} sim/$*_$(res).mrc res=$(res) het center; \
	echo "./fit_consistency $${structure} $${measurement} sim/$*_$(res).mrc\n"; \
	$< sim/$*_$(res).mrc $${structure} $${measurement}

old_simulate/%: 
	@ structure=$(shell find data/ -name "$*.pdb"); \
	phenix.fmodel $${structure} high_resolution=$(res) generate_fake_p1_symmetry=True;\
	phenix.mtz2map mtz_file=$(*F).pdb.mtz pdb_file=$${structure} labels=FMODEL,PHIFMODEL output.prefix=$(*F);\
	rm $(*F).pdb.mtz;\
	mv $(*F)_fmodel.ccp4 sim/$(*F)_$(res).ccp4;\

stuff/%: build/executable/stuff data/%.pdb
#	@$< data/$*.pdb sim/native_20.ccp4 sim/native_21.ccp4 sim/native_22.ccp4 sim/native_23.ccp4
	@$< data/$*.pdb $(shell find sim/ -name "$**" -printf "%p\n" | sort | awk '{printf("%s ", $$0)}')


####################################################################################
###				TESTS						 ###
####################################################################################
tags := ""
exclude_tags := "~[broken] ~[manual] ~[slow] ~[disable]"
memtest/%: $(shell find source/ -print) test/%.cpp	
	valgrind --track-origins=yes --log-file="valgrind.txt" build/test [$(*F)] ${tags}

tests: $(shell find source/ -print) $(shell find test/ -print)
	@ make -C build tests -j${cmake_threads}
	@ for test in $$(find build/test/*); do\
		$${test} $(exclude_tags);\
	done

test/%: test/%.cpp
	@ make -C build "test_$*" -j${cmake_threads}
	build/test/test_$* ~[slow] ~[broken] ${tags}

# special build target for our tests since they obviously depend on themselves, which is not included in $(source_files)
build/test/%: $(shell find source/ -print) $(shell find test -name *%.cpp) build/Makefile
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

winbuild/executable/%: $(source) $(include) executable/%.cpp
	@ cmake --build winbuild/ --target $(*F) -j${cmake_threads}

winbuild/%: $(source) $(include)
	@ cmake --build winbuild/ --target $(*F) -j${cmake_threads} 

build/executable/%: $(source) $(include) executable/%.cpp
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
