# ------------------------------------------------------------------------------
# rundocker.sh - Wrapper script for executing Wigner simulations in Docker
#
# This script provides a convenient way to run Wigner simulations within a Docker
# container. It mounts necessary volumes, sets optional resource limits (CPU and memory),
# and optionally includes a configuration file.
#
# Usage:
#   ./rundocker.sh <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [cpus] [memory] [config_file]
#
# Example:
#   ./rundocker.sh 0.001 2.0 8 0.005 simres res 8 12g input.txt
#
# Arguments:
# - <start>         : Starting k* value.
# - <end>           : Ending k* value.
# - <n_jobs>        : Number of parallel jobs.
# - <increment>     : Step increment within jobs.
# - <output_folder> : Directory in container for outputs (mapped to wigner_output on host).
# - <file_prefix>   : Prefix for output files.
# - [cpus]          : (Optional) CPU limit.
# - [memory]        : (Optional) Memory limit.
# - [config_file]   : (Optional) Path to the configuration file.
# ------------------------------------------------------------------------------


#!/bin/bash
set -euo pipefail


if [ "$#" -lt 6 ] || [ "$#" -gt 9 ]; then
    echo "Usage: $0 <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [config_file] [cpus] [memory]"
    exit 1
fi

START=$1
END=$2
NJOBS=$3
INCREMENT=$4
OUTDIR=$5
PREFIX=$6
CPUS=${7:-}
MEMORY=${8:-}
CONFIG_FILE=${9:-}


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
DOCKER_CMD+=("wignerutils" "$START" "$END" "$NJOBS" "$INCREMENT" "$OUTDIR" "$PREFIX")

# Append config file as the final argument if provided
if [[ -n "$CONFIG_FILE" ]]; then
  DOCKER_CMD+=("$CONFIG_FILE")
fi

# Show and run
echo "Running simulation with:"
printf '  %q\n' "${DOCKER_CMD[@]}"
"${DOCKER_CMD[@]}"
