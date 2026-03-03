#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=dl-data
#SBATCH --output=slurm-dl-data-%j.out

# FTP Download Script
# Downloads CGVRLL128235 folder from FTP server

# FTP connection details
FTP_HOST="198.2.254.43"
FTP_PORT="2122"
FTP_USER="cdgenomics-ftp3"
FTP_PASS="XXX"
REMOTE_FOLDER="CGVRLL128235"

# Local download directory (optional - will download to current directory if not set)
LOCAL_DIR="./downloads"

# Create local directory if it doesn't exist
mkdir -p "$LOCAL_DIR"

echo "Starting FTP download..."
echo "Server: $FTP_HOST:$FTP_PORT"
echo "Remote folder: $REMOTE_FOLDER"
echo "Local directory: $LOCAL_DIR"

# Download using lftp
lftp -c "
set ssl:verify-certificate no
set net:timeout 10
set net:max-retries 3
open -u $FTP_USER,$FTP_PASS $FTP_HOST:$FTP_PORT
lcd $LOCAL_DIR
mirror --verbose --parallel=4 --continue $REMOTE_FOLDER
quit
"

# Check if download was successful
if [ $? -eq 0 ]; then
    echo "Download completed successfully!"
    echo "Files downloaded to: $LOCAL_DIR/$REMOTE_FOLDER"
else
    echo "Download failed!"
    exit 1
fi
