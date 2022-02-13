.PHONY: build run test benchmark clean

build:
	cmake -E make_directory $(CURDIR)/release
	cmake -S $(CURDIR) -B $(CURDIR)/release -DCMAKE_BUILD_TYPE=Release -DUSE_WARNINGS=ON -DENABLE_TESTING=OFF
	cmake --build $(CURDIR)/release --parallel

run: build
	./release/acm 64

clean:
	$(MAKE) -C release clean

tbuild:
	cmake -E make_directory $(CURDIR)/debug
	cmake -S $(CURDIR) -B $(CURDIR)/debug -DCMAKE_BUILD_TYPE=Debug -DUSE_WARNINGS=ON -DENABLE_TESTING=ON
	cmake --build $(CURDIR)/debug --parallel

test: tbuild
	ctest -V --test-dir debug
