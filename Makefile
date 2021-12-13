# generate lists of files for easy use as dependencies
pymol := pymol

test_files := $(basename $(shell find source/tests/ -maxdepth 1 -name "*.cpp" -printf "%P "))
source_files := $(addprefix source/, $(shell find source/ -type f -not -wholename "source/tests/*" -printf "%P "))

#################################################################################
###				EXECUTABLES					 ###
#################################################################################
.SECONDARY:

.PHONY:
gui: build/source/scripts/gui
	$<

.PHONY:
hydrate/%: build/source/scripts/new_hydration
	$< data/$*.pdb output/$*.pdb
	$(pymol) output/$*.pdb -d "show spheres; color orange, hetatm"

.PHONY:
hist/%: build/source/scripts/hist
	$< data/$*.pdb figures/

.PHONY:
main/%: build/source/scripts/main
	$< data/$*.pdb output/filled_volume.pdb

.PHONY:
optimize_radius/%: build/source/scripts/optimize_radius
	$< data/$*.pdb figures/
	
.PHONY:
intensity_fit/%: build/source/scripts/intensity_fitter
	$< data/$*.pdb data/$*.RSR figures/

#################################################################################
###				TESTS						 ###
#################################################################################
tests: $(addprefix build/source/tests/, $(test_files))
	for program in $^ ; do \
	    $$program ; \
	done
	
test/%: build/source/tests/% source/tests/%.cpp
	$<

# special build target for our tests since they obviously depend on themselves, which is not included in $(source_files)
build/source/tests/%: $(shell find source/ -print) build/Makefile
	@ make -C build $*	

#################################################################################
###				BUILD						 ###
#################################################################################
build/%: $(source_files) build/Makefile
	@ cmake --build build/ --target $(*F)
	
build/Makefile: $(shell find -name "CMakeLists.txt" -printf "%P ")
	@ mkdir -p build
	@ cd build; cmake ../
	
clean/build: 
	@ rmdir -f build
