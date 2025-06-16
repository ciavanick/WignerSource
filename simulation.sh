#!/bin/bash

# ------------------------------------------------------------------------------
# simulation.sh - Runs parallel simulations and post-processing
#
# Usage:
#   ./simulation.sh <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [config_file]
#
# Example:
#   ./simulation.sh 0.001 2.0 8 0.005 simres res input.txt
#
# Explanation:
# - Runs n parallel jobs covering the k* range from 0.001 to 2.0.
# - Each job incrementally steps by 0.005 in the range.
# - Output ROOT files are saved into the directory `simres/`.
# - Files are named using the prefix `res`.
# - Simulation parameters are read from `input.txt`.
# ------------------------------------------------------------------------------



set -euo pipefail
source wignerenv.sh
export LC_NUMERIC=C



if [ "$#" -lt 6 ] || [ "$#" -gt 7 ]; then
    echo "Usage: $0 <start> <end> <n_jobs> <increment> <output_folder> <file_prefix> [config_file]"
    exit 1
fi

# Parse arguments
START=$1
END=$2
NJOBS=$3
INCREMENT=$4
OUTDIR=$5
PREFIX=$6
CONFIG_FILE=${7:-config/default.txt}  # Default value

echo "START=$START END=$END NJOBS=$NJOBS INCREMENT=$INCREMENT OUTDIR=$OUTDIR PREFIX=$PREFIX CONFIG_FILE=$CONFIG_FILE"

# Compute range step using awk
RANGE_STEP=$(awk "BEGIN { printf \"%.6f\", ($END - $START)/$NJOBS }")
echo "RANGE_STEP=$RANGE_STEP"

mkdir -p "$OUTDIR"

# Loop over jobs
for ((i=0; i<NJOBS; i++)); do
    JOB_START=$(awk "BEGIN { printf \"%.6f\", $START + $i * $RANGE_STEP }")
    JOB_END=$(awk "BEGIN { printf \"%.6f\", $JOB_START + $RANGE_STEP }")
    OUTFILE="$OUTDIR/${PREFIX}_part_$i.root"

    echo "Launching job $i: [$JOB_START, $JOB_END] -> $OUTFILE"

    if [[ -z "$JOB_START" || -z "$JOB_END" ]]; then
        echo "ERROR: Empty range in job $i"
        exit 1
    fi

    wigneroot -l -b -q "macros/wignersim.cpp($JOB_START, $JOB_END, $INCREMENT, \"$OUTFILE\", \"$CONFIG_FILE\")" &
done

wait

MERGED="$OUTDIR/${PREFIX}_merged.root"
echo "Merging to $MERGED"
hadd -f "$MERGED" "$OUTDIR/${PREFIX}_part_"*.root

echo "All done. Merged output: $MERGED"
echo "Making the plots"

wigneroot -l -q -b "macros/makeplots.cpp(\"$OUTDIR\",\"${PREFIX}_merged.root\")"
