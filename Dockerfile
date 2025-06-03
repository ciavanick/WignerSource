FROM rootproject/root:latest

# Install system packages (if needed for ROOT or build tools)
RUN apt-get update && \
    apt-get install -y cmake g++ make && \
    rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /wigner

# Copy all project files (CMakeLists.txt, src/, include/, etc.)
COPY . .

# Create and enter build directory
RUN mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/wigner .. && \
    cmake --build . -- -j$(nproc) && \
    cmake --install .

# Confirm shared libraries were built
RUN ls /wigner/lib

# Run ROOT test with macro
CMD root -l -q wignertest2.cpp++
