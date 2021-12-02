# generate lists of files for easy use as dependencies
pymol := ~/tools/pymol/pymol

test_files := $(basename $(shell find source/tests/ -maxdepth 1 -name "*.cpp" -printf "%P "))
source_files := $(addprefix source/, $(shell find source/ -type f -not -wholename "source/tests/*" -printf "%P "))

#################################################################################
###				EXECUTABLES					 ###
#################################################################################
.phony:
hydrate/%: build/source/scripts/new_hydration
	$< data/$* output/$*
	$(pymol) output/$* -d "show spheres; color orange, hetatm"

.phony:
hist/%: build/source/scripts/hist
	$< data/$* figures/

.phony:
main/%: build/source/scripts/main
	$< data/$* output/filled_volume.pdb

.phony:
optimize_radius/%: build/source/scripts/optimize_radius
	$< data/$* figures/
	
.phony:
intensity_fit: build/source/scripts/intensity_fitter
	$< data/2epe.pdb data/Lyz_03m.RSR figures/

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
