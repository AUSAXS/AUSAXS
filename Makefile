# generate lists of files for easy use as dependencies
pymol := ~/tools/pymol/pymol

test_files := $(basename $(shell find source/tests/ -maxdepth 1 -name "*.cpp" -printf "%P "))
source_files := $(addprefix source/, $(shell find source/ -type f -not -wholename "source/tests/*" -printf "%P "))

.phony:
hydrate/%: build/source/scripts/new_hydration
	$< data/$* output/$*
	$(pymol) output/$* -d "show spheres; color orange, hetatm"

.phony:
hist: build/source/scripts/hist
	$< data/$* output/$*

tests: $(addprefix build/source/tests/, $(test_files))
	for program in $^ ; do \
	    $$program ; \
	done
	
test/%: build/source/tests/%
	$<
	
build/%: $(source_files) build/Makefile
	@ cmake --build build/ --target $(*F)
	
build: $(shell find -name "CMakeLists.txt" -printf "%P ") build/Makefile
	cd build; cmake ../
	
build/Makefile:
	@ mkdir -p build
	@ cd build; cmake ../
