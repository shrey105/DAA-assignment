# DAA Assignment — Densest Subgraph Algorithms

Four algorithms for finding the densest subgraph:

| Algorithm | Source | Paper |
|---|---|---|
| `core_exact` | `dense-subgraph/core_exact.cpp` | Dense Subgraph |
| `exact` | `dense-subgraph/exact.cpp` | Dense Subgraph |
| `exact_parallel` | `dense-subgraph/exact_parallel.cpp` | Dense Subgraph |
| `greedy_plus_plus` | `flowless/greedy_plus_plus.cpp` | Flowless |

## Prerequisites

- `g++` with C++17 support

## Build

Build all four algorithms at once:

```bash
make
```

Or build individually:

```bash
make core_exact
make exact
make exact_parallel
make greedy_plus_plus
```

## Run

Each algorithm accepts an input graph file. The input format is an edge list: one `u v` pair per line.

### core_exact

```bash
# compile + run
make run-core_exact INPUT=path/to/graph.txt

# or run the binary directly
./dense-subgraph/core_exact <input_file>
```

Output is written to `<input_basename>_result.txt` in the same directory as the input file.

### exact

```bash
make run-exact INPUT=path/to/graph.txt OUTPUT=path/to/output.txt

# binary directly
./dense-subgraph/exact <input_file> [output_file]
```

If `output_file` is omitted, results are printed to stdout.

### exact_parallel

```bash
make run-exact_parallel INPUT=path/to/graph.txt OUTPUT=path/to/output.txt

# binary directly
./dense-subgraph/exact_parallel <input_file> [output_file]
```

If `output_file` is omitted, results are printed to stdout.

### greedy_plus_plus

```bash
# T controls the number of iterations (default: 10)
make run-greedy_plus_plus INPUT=path/to/graph.txt T=10 OUTPUT=path/to/output.txt

# binary directly
./flowless/greedy_plus_plus <input_file> [T] [output_file]
```

## Makefile variables

| Variable | Default | Description |
|---|---|---|
| `INPUT` | `input.txt` | Path to the input graph file |
| `OUTPUT` | `output.txt` | Path to the output file |
| `T` | `10` | Iteration count for `greedy_plus_plus` |

## Clean

```bash
make clean
```
