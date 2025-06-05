#!/bin/bash
set -euo pipefail

# Usage:
# ./simulation.sh <start> <end> <n_jobs> <increment> <output_folder> <file_prefix>

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <start> <end> <n_jobs> <increment> <output_folder> <file_prefix>"
    exit 1
fi

# Parse arguments
START=$1
END=$2
NJOBS=$3
INCREMENT=$4
OUTDIR=$5
PREFIX=$6

# Validate input
echo "START=$START END=$END NJOBS=$NJOBS INCREMENT=$INCREMENT OUTDIR=$OUTDIR PREFIX=$PREFIX"

# Compute range per job
RANGE_STEP=$(echo "scale=6; ($END - $START)/$NJOBS" | bc -l)
echo "RANGE_STEP=$RANGE_STEP"

# Create output directory if not exists
mkdir -p "$OUTDIR"

# Loop over jobs
for ((i=0; i<NJOBS; i++)); do
    JOB_START=$(echo "scale=6; $START + $i * $RANGE_STEP" | bc -l)
    JOB_END=$(echo "scale=6; $JOB_START + $RANGE_STEP" | bc -l)
    OUTFILE="$OUTDIR/${PREFIX}_part_$i.root"

    echo "Launching job $i: [$JOB_START, $JOB_END] -> $OUTFILE"

    # Safety check
    if [[ -z "$JOB_START" || -z "$JOB_END" ]]; then
        echo "ERROR: Empty range in job $i"
        exit 1
    fi

    # Run the ROOT macro in the background
    root -l -b -q "wignersim.cpp($JOB_START, $JOB_END, $INCREMENT, \"$OUTFILE\")" &
done

# Wait for all background jobs to complete
wait

# Merge all part files
MERGED="$OUTDIR/${PREFIX}_merged.root"
echo "Merging to $MERGED"
hadd -f "$MERGED" "$OUTDIR/${PREFIX}_part_"*.root

echo "✅ All done. Merged output: $MERGED"
