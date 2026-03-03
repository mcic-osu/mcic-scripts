#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --output=slurm-cntfiles-%j.out

# This script counts the number of files in each top-level dir below the focal dir

# Parse args
focaldir=$1

# Report
echo "# Starting script cntfiles.sh"
date
echo -e "\n# Focal dir: $focaldir \n"

# Move to focal dir
cd "$focaldir" || exit

# Count files
find . -type f | cut -d/ -f2 | sort | uniq -c

# Report
echo -e "\n# Done."
date
