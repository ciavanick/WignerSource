##
# @addtogroup DockerSetup
# @{
##
##
# @file Dockerfile
# @brief Docker image definition for the Wigner simulation environment.
#
# This Dockerfile builds an image configured for running Wigner simulations
# and analysis macros using ROOT. It includes the required compiler (gcc-12),
# dependencies, and sets up the build environment for the Wigner codebase.
#
# Default execution (`simulation.sh`) runs a predefined simulation.
#
# Build command:
# @code
# docker build -t wignerutils .
# @endcode
#
# Example run command (direct Docker):
# @code
# docker run --cpus="8" --memory="12g" --rm -it \
#   -v $(pwd)/wigner_output:/wigner/simtestn \
#   -v $(pwd)/input_data.txt:/wigner/input_data.txt \
#   wignerutils 0.001 0.5 8 0.005 simtestn result input_data.txt
# @endcode
#
##

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
ENTRYPOINT ["./simulation.sh"]
CMD ["0.001", "0.5", "8", "0.005", "simtest", "data", "config/default.txt"]
##
# @}
##