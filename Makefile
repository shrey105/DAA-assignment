CXX      := g++
CXXFLAGS := -std=c++17 -O2

BINDIR_DS      := dense-subgraph
BINDIR_FL      := flowless

CORE_EXACT     := $(BINDIR_DS)/core_exact
EXACT          := $(BINDIR_DS)/exact
GREEDY 		   := $(BINDIR_FL)/greedy
GREEDY_PP      := $(BINDIR_FL)/greedy_plus_plus

INPUT   ?= input.txt
OUTPUT  ?= output.txt
T       ?= 10

.PHONY: all clean \
        core_exact exact greedy greedy_plus_plus \
        run-core_exact run-exact run-greedy run-greedy_plus_plus

all: core_exact exact greedy greedy_plus_plus
run: core_exact exact greedy greedy_plus_plus
	@base=$$(basename $(INPUT) .txt); \
	echo "Running on $(INPUT)"; \
	$(CORE_EXACT) $(INPUT) $${base}_core_exact.txt; \
	$(EXACT) $(INPUT) $${base}_exact.txt; \
	$(GREEDY) $(INPUT) $${base}_greedy.txt; \
	$(GREEDY_PP) $(INPUT) $(T) $${base}_greedy_pp.txt; \
	echo "Outputs:"; \
	echo "  $${base}_core_exact.txt"; \
	echo "  $${base}_exact.txt"; \
	echo "  $${base}_greedy.txt"; \
	echo "  $${base}_greedy_pp.txt"

core_exact: $(CORE_EXACT)
$(CORE_EXACT): $(BINDIR_DS)/core_exact.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

exact: $(EXACT)
$(EXACT): $(BINDIR_DS)/exact.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

greedy: $(GREEDY)
$(GREEDY): $(BINDIR_FL)/greedy.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

greedy_plus_plus: $(GREEDY_PP)
$(GREEDY_PP): $(BINDIR_FL)/greedy_plus_plus.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run-core_exact: core_exact
	$(CORE_EXACT) $(INPUT)

run-exact: exact
	$(EXACT) $(INPUT) $(OUTPUT)

run-greedy: greedy
	$(GREEDY) $(INPUT) $(OUTPUT)

run-greedy_plus_plus: greedy_plus_plus
	$(GREEDY_PP) $(INPUT) $(T) $(OUTPUT)

clean:
	rm -f $(CORE_EXACT) $(EXACT) $(GREEDY) $(GREEDY_PP)
