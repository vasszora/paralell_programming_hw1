.PHONY: build run test benchmark

build:
	mkdir -p build
	cd build; cmake .. -DCMAKE_BUILD_TYPE=Release; make -j

run: build
	./build/acm
