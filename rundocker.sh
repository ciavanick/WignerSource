#!/bin/bash

# ------------------------------------------------------------------------------
# rundocker.sh - Wrapper script for executing Wigner simulations in Docker
#
# Usage:
#   ./rundocker.sh <image_name> <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [cpus] [memory] [config_file]
#
# Example:
#   ./rundocker.sh wignerutils 0.001 2.0 8 0.005 simres res 8 12g input.txt
# ------------------------------------------------------------------------------

set -euo pipefail
export LC_NUMERIC=C

if [ "$#" -lt 7 ] || [ "$#" -gt 10 ]; then
    echo "Usage: $0 <image_name> <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [cpus] [memory] [config_file]"
    exit 1
fi

IMAGE_NAME=$1
START=$2
END=$3
NJOBS=$4
INCREMENT=$5
OUTDIR=$6
PREFIX=$7
CPUS=${8:-}
MEMORY=${9:-}
CONFIG_FILE=${10:-}

# Ensure the output folder exists on host
mkdir -p wigner_output

# Start building the Docker command
DOCKER_CMD=(
  docker run
  --rm -it
  -v "$(pwd)/wigner_output:/wigner/$OUTDIR"
)

# Add CPU/memory only if non-empty
if [[ -n "${CPUS}" ]]; then
  DOCKER_CMD+=(--cpus="$CPUS")
fi

if [[ -n "${MEMORY}" ]]; then
  DOCKER_CMD+=(--memory="$MEMORY")
fi

# Add config file mount if provided
if [[ -n "$CONFIG_FILE" ]]; then
  if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file '$CONFIG_FILE' not found."
    exit 1
  fi
  DOCKER_CMD+=("-v" "$(pwd)/$CONFIG_FILE:/wigner/$CONFIG_FILE")
fi

# Append image and simulation args
DOCKER_CMD+=("$IMAGE_NAME" "$START" "$END" "$NJOBS" "$INCREMENT" "$OUTDIR" "$PREFIX")

# Append config file as the final argument if provided
if [[ -n "$CONFIG_FILE" ]]; then
  DOCKER_CMD+=("$CONFIG_FILE")
fi

# Show and run
echo "Running simulation with:"
printf '  %q\n' "${DOCKER_CMD[@]}"
"${DOCKER_CMD[@]}"
