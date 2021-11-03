run: build/main
	@ build/main

build/main: main.cpp source/*
	mkdir -p build
	cd build; cmake ../; make
