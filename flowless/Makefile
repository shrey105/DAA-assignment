CXX      := g++
CXXFLAGS := -std=c++17 -O2

INPUT  ?= input.txt
OUTPUT ?= output.txt
T      ?= 10

.PHONY: all clean greedy_plus_plus run-greedy_plus_plus

all: greedy_plus_plus

greedy_plus_plus: greedy_plus_plus.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run-greedy_plus_plus: greedy_plus_plus
	./greedy_plus_plus $(INPUT) $(T) $(OUTPUT)

clean:
	rm -f greedy_plus_plus
