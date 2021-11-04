run/%: build/main main.cpp source/* source/io/*
	build/main data/$*

build/main: main.cpp source/* build/Makefile
	make -C build main

test_files := $(basename $(shell find source/tests/ -maxdepth 1 -name "*.cpp" -printf "%P "))
tests: $(addprefix build/source/tests/, $(test_files))
	for program in $^ ; do \
	    $$program ; \
	done
	
build/source/tests/%: source/tests/%.cpp build/Makefile
	@ make -C build $*

build/Makefile: CMakeLists.txt source/tests/CMakeLists.txt
	@ mkdir -p build
	@ cd build; cmake ../
