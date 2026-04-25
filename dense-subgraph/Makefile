CXX      := g++
CXXFLAGS := -std=c++17 -O2

INPUT  ?= input.txt
OUTPUT ?= output.txt

.PHONY: all clean core_exact exact exact_parallel \
        run-core_exact run-exact run-exact_parallel

all: core_exact exact exact_parallel

core_exact: core_exact.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

exact: exact.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

exact_parallel: exact_parallel.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run-core_exact: core_exact
	./core_exact $(INPUT)

run-exact: exact
	./exact $(INPUT) $(OUTPUT)

run-exact_parallel: exact_parallel
	./exact_parallel $(INPUT) $(OUTPUT)

clean:
	rm -f core_exact exact exact_parallel
