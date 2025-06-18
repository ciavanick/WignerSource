# Scripts and Setup

This document describes the key scripts and setup files used for running and analyzing Wigner coalescence simulations using ROOT. It includes instructions for both local and Docker-based workflows.

---

## Environment Setup — `wignerenv.sh`

**Purpose**:  
Initializes the environment by configuring paths needed to run the simulation and plotting macros.

**What it does**:
- Adds the `bin/` directory to the `PATH`
- Adds the `lib/` directory to the `LD_LIBRARY_PATH`

**Usage**:
```bash
source wignerenv.sh
```

**Note**: This script is automatically sourced inside `simulation.sh`, but must be run manually when working interactively.

---

## Local Simulation — `simulation.sh`

**Purpose**:  
Runs parallel simulations by dividing the momentum (`k*`) range and executing ROOT macros in multiple jobs.

**Usage**:
```bash
./simulation.sh <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [config_file]
```

**Example**:
```bash
./simulation.sh 0.001 2.0 8 0.005 simres res input.txt
```

**What it does**:
- Splits the total range of `k*` into `n_jobs` intervals
- Each job executes `macros/wignersim.cpp` over a subrange
- Merges the resulting `.root` files using `hadd`
- Runs the plotting macro `macros/makeplots.cpp` on the merged file

---

## Docker-Based Execution — `rundocker.sh` and `Dockerfile`

### `Dockerfile`

**Purpose**:  
Defines the build for a Docker image containing ROOT, GCC-12, and the Wigner simulation tools.

**How to build**:
```bash
docker build -t wignerutils .
```

### `rundocker.sh`

**Purpose**:  
Wrapper for launching simulations inside a Docker container, with optional CPU, memory, and configuration bindings.

**Usage**:
```bash
./rundocker.sh  <image_name> <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [cpus] [memory] [config_file]
```

**Example**:
```bash
./rundocker.sh wignerutils 0.001 2.0 8 0.005 simres res 8 12g input.txt
```

**Details**:
- Automatically mounts `wigner_output` on host to container output folder
- Optional parameters for resource limits (CPUs, memory)
- Supports mounting a config `.txt` file into the container

---

## Notes

- The default configuration file is `config/default.txt` if none is specified.
- Scripts assume you are in the root of the repository.
- All ROOT macros are run using `wigneroot` in batch (`-b`) and quiet (`-q`) mode.
