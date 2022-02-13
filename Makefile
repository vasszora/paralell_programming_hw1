.PHONY: build run test benchmark clean

build:
	cmake -E make_directory $(CURDIR)/release
	cmake -S $(CURDIR) -B $(CURDIR)/release -DCMAKE_BUILD_TYPE=Release -DUSE_WARNINGS=ON
	cmake --build $(CURDIR)/release --parallel

run: build
	./release/acm 64

clean:
	$(MAKE) -C release clean
