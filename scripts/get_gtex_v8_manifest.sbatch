#!/bin/bash

#SBATCH --job-name=get_gtex_v8_manifest
#SBATCH --output=get_gtex_v8_manifest.%j.out
#SBATCH --error=get_gtex_v8_manifest.%j.err
#SBATCH --time=00:15:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=pi-mstephens

###Run script
find ../data/cis_eqtl_analysis_ready/ -maxdepth 1 -type f | grep 'ENSG' > ../data/gtex-v8-manifest.txt
