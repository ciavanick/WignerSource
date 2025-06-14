# syntax=docker/dockerfile:1.4
FROM --platform=linux/amd64 rootproject/root:latest

# Install a stable compiler (gcc-12) to avoid ICEs
RUN apt-get update && \
    apt-get install -y g++-12 cmake make time bc && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100 && \
    rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /wigner
COPY . .

# Build the project
RUN mkdir -p build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/wigner .. && \
    cmake --build . && \
    cmake --install .

# Make sure the script is executable
RUN chmod +x simulation.sh

# Run simulation by default
CMD ./simulation.sh 0.001 0.5 8 0.005 simtest data
