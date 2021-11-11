# generate lists of files for easy use as dependencies
test_files := $(basename $(shell find source/tests/ -maxdepth 1 -name "*.cpp" -printf "%P "))
source_files := $(addprefix source/, $(shell find source/ -type f -not -wholename "source/tests/*" -printf "%P "))

.phony:
hydrate/%: build/source/scripts/new_hydration
	@ echo $(source_files)
	$< data/$* output/$*

run/%: build/main main.cpp source/* source/io/*
	build/main data/$*

build/main: main.cpp source/* build/Makefile
	make -C build main

tests: $(addprefix build/source/tests/, $(test_files))
	for program in $^ ; do \
	    $$program ; \
	done
	
test/%: build/source/tests/%
	$<
	
build/source/tests/%: $(shell find source/ -print) build/Makefile
	@ make -C build $*

build/%: $(source_files)
	@ make -C $(@D) $(*F)
	
build: $(shell find -name "CMakeLists.txt" -printf "%P ") build/Makefile
	cd build; cmake ../
	
build/Makefile:
	@ mkdir -p build
	@ cd build; cmake ../
