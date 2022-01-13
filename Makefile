# generate lists of files for easy use as dependencies
pymol := pymol

cmake_threads := 6

test_files := $(basename $(shell find source/tests/ -maxdepth 1 -name "*.cpp" -printf "%P "))
source_files := $(addprefix source/, $(shell find source/ -type f -not -wholename "source/tests/*" -printf "%P "))

#################################################################################
###				EXECUTABLES					 ###
#################################################################################
.SECONDARY:

.PHONY:
gui: build/source/gui/gui
	build/gui

gwidth := 1
bwidth := 1
ra := 2.4
rh := 1.5
ps := Radial
hydrate/%: build/source/scripts/new_hydration
	$< data/$*.pdb output/$*.pdb --grid_width ${gwidth} --radius_a ${ra} --radius_h ${rh} --placement_strategy ${ps}
#	$(pymol) output/$*.pdb -d "show spheres; color orange, hetatm"
	$(pymol) output/$*.pdb -d "hide all; show spheres, hetatm; color orange, hetatm"

hist/%: build/source/scripts/hist
	$< data/$*.pdb figures/ --grid_width ${gwidth} --radius_a ${ra} --radius_h ${rh} --bin_width ${bwidth} --placement_strategy ${ps}

main/%: build/source/scripts/main
	$< data/$*.pdb output/filled_volume.pdb

optimize_radius/%: build/source/scripts/optimize_radius
	$< data/$*.pdb figures/
	
qlow := 0
qhigh := 1000
center := center
intensity_fit/%: build/source/scripts/intensity_fitter
	$< data/$*.pdb data/$*.RSR figures/ --qlow ${qlow} --qhigh ${qhigh} --${center} --radius_a ${ra} --radius_h ${rh} --grid_width ${gwidth} --bin_width ${bwidth} --placement_strategy ${ps}

#################################################################################
###				TESTS						 ###
#################################################################################
test/%: $(shell find source/ -print) test/%.cpp
	@ make -C build test
	build/test [$(*F)] ~[broken] ~[manual]

tests: $(shell find source/ -print) $(shell find test/ -print) build/Makefile
	@ make -C build test
	build/test ~[broken] ~[manual]

# special build target for our tests since they obviously depend on themselves, which is not included in $(source_files)
build/source/tests/%: $(shell find source/ -print) build/Makefile
	@ make -C build $* -j${cmake_threads}

#################################################################################
###				BUILD						 ###
#################################################################################
.PHONY: build
build: 
	@ mkdir -p build; 
	@ cd build; cmake ../ -j${cmake_threads}

build/%: $(source_files) build/Makefile
	@ cmake --build build/ --target $(*F) -j${cmake_threads}
	
build/Makefile: $(shell find -name "CMakeLists.txt" -printf "%P ")
	@ mkdir -p build
	@ cd build; cmake ../ -j${cmake_threads}
	
clean/build: 
	@ rmdir -f build

