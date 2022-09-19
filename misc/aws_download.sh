## Load the aws-cli conda environment
module load miniconda3
source activate /fs/project/PAS0471/jelmer/conda/aswcli

## Configure credentials (the last one, 'output' can probably be left empty)
aws configure

## Download the entire bucket into a new directory 'BGI_data' (change as wanted) 
aws s3 sync s3://f22ftsusat0591/F22FTSUSAT0591_HUMfexdR/ BGI_data

## Download one file
# aws s3api get-object --bucket f22ftsusat0591 --key F22FTSUSAT0591_HUMfexdR/BGI_F22FTSUSAT0591_HUMfexdR.report.en.pdf BGI_F22FTSUSAT0591_HUMfexdR.report.en.pdf
