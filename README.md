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
- Bash (for running `run.sh`)

---

## Installation

Clone the repository:

```bash
git clone https://github.com/ciavanick/WignerSource.git
cd WignerSource
```

## How to Compile and Run

To run the code, you need to compile `CWignerSource.cpp` and `CWignerUtils.cpp`, and link them properly with ROOT.  
This is handled automatically via the `run.sh` script and the `rootlogon.C` file.

On the first run, execute:

```bash
source run.sh
```