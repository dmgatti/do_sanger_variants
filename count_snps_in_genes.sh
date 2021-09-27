#!/bin/bash
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 6 # number of cores
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-4:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

CONTAINER=~/containers/bioconductor.sif
SCRIPT=/projects/compsci/dgatti/projects/do_sanger_variants/count_snps_in_genes.R

module load singularity

singularity exec ${CONTAINER} Rscript ${SCRIPT}
