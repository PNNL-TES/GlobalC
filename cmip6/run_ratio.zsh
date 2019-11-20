#!/bin/zsh
#SBATCH -A FORTEPROJECT
#SBATCH -t 90
#SBATCH -N 1
 
source /etc/profile.d/modules.sh >& /dev/null
module load R/3.4.3 
module load gcc/5.2.0
 
date
Rscript /pic/projects/GCAM/Dorheim/Dorheim/GlobalC/cmip6/1.cmip6_RStoGPP.R
date
