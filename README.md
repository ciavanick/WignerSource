# WignerSource

**WignerSource** is a C++ module for computing and handling Wigner functions using the ROOT framework. It provides utilities for simulating quantum mechanical systems and visualizing phase space distributions.

---

## Features

- Utilities for Wigner function manipulation
- Seamless integration with ROOT
- Easy compilation with provided script

---

## Requirements

- [ROOT Framework](https://root.cern/)
  - Ensure `root-config` is available in your environment
- C++17 or higher
- cmake

---

## Compilation and Installation

Clone the repository:

```bash
git clone https://github.com/ciavanick/WignerSource.git
cd WignerSource
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=..
cmake --build .
cmake --install .
cd ..
```

## How to Run

```bash
root -l wignertest2.cpp
```
## Docker

```bash
docker build -t wignerutils .
docker run --rm -it wignerutils 
```

If you have an arm achitecture:

```bash
docker build --platform linux/amd64 -t wignerutils .
docker run --platform linux/amd64 --rm -it wignerutils 
```
