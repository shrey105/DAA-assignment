CXX      := g++
CXXFLAGS := -std=c++17 -O2

BINDIR_DS      := dense-subgraph
BINDIR_FL      := flowless

CORE_EXACT     := $(BINDIR_DS)/core_exact
EXACT          := $(BINDIR_DS)/exact
EXACT_PARALLEL := $(BINDIR_DS)/exact_parallel
GREEDY_PP      := $(BINDIR_FL)/greedy_plus_plus

INPUT   ?= input.txt
OUTPUT  ?= output.txt
T       ?= 10

.PHONY: all clean \
        core_exact exact exact_parallel greedy_plus_plus \
        run-core_exact run-exact run-exact_parallel run-greedy_plus_plus

all: core_exact exact exact_parallel greedy_plus_plus

core_exact: $(CORE_EXACT)
$(CORE_EXACT): $(BINDIR_DS)/core_exact.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

exact: $(EXACT)
$(EXACT): $(BINDIR_DS)/exact.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

exact_parallel: $(EXACT_PARALLEL)
$(EXACT_PARALLEL): $(BINDIR_DS)/exact_parallel.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

greedy_plus_plus: $(GREEDY_PP)
$(GREEDY_PP): $(BINDIR_FL)/greedy_plus_plus.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run-core_exact: core_exact
	$(CORE_EXACT) $(INPUT)

run-exact: exact
	$(EXACT) $(INPUT) $(OUTPUT)

run-exact_parallel: exact_parallel
	$(EXACT_PARALLEL) $(INPUT) $(OUTPUT)

run-greedy_plus_plus: greedy_plus_plus
	$(GREEDY_PP) $(INPUT) $(T) $(OUTPUT)

clean:
	rm -f $(CORE_EXACT) $(EXACT) $(EXACT_PARALLEL) $(GREEDY_PP)
