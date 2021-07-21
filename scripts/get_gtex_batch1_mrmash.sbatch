#!/bin/bash

#SBATCH --job-name=get_gtex_batch1_mrmash
#SBATCH --output=get_gtex_batch1_mrmash.%j.out
#SBATCH --error=get_gtex_batch1_mrmash.%j.err
#SBATCH --time=00:05:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=pi-mstephens

head -n 1000 ../data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix.txt > ../data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-batch1.txt
