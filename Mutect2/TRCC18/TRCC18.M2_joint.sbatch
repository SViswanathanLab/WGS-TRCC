#!/bin/bash
#SBATCH --mem=20G
#SBATCH --job-name=TRCC18.M2_joint
#SBATCH --array=1-300
#SBATCH --output=/czlab/cbao/TRCC/Mutect2_joint/TRCC18/hg38_scatter_300/TRCC18.M2_joint.%A_%a.out.log
#SBATCH --error=/czlab/cbao/TRCC/Mutect2_joint/TRCC18/hg38_scatter_300/TRCC18.M2_joint.%A_%a.err.log
sed -n "$SLURM_ARRAY_TASK_ID"p /czlab/cbao/TRCC/Mutect2_joint/TRCC18/TRCC18.M2_joint.s300.sh | bash

