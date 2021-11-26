#!/bin/bash

#SBATCH --account=def-wyeth

## Mail Options
#SBATCH --mail-user=melsiddieg@cmmt.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
## CPU Usage
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-12:00
#SBATCH --ntasks=1
module load nextflow
nextflow run my_structural6.nf --input samples.csv --tfile transposons.txt --fasta Inputs/hs37d5.fa -profile slurm  -with-trace -with-report -c my_nextflow.config --catalog Inputs/extData/variant_catalog.json  --annotation_file hg19.genes.bed --rfile cr.bed.gz -with-dag flowchart.png -resume
