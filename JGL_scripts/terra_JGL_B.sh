#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=getNetworkTAM111up      #Set the job name to "JobExample2"
#SBATCH --time=48:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=2         #72 for hisat-build           #Request 1 node
#SBATCH --ntasks-per-node=2     #18 for hisat-build     #Request 8 tasks/cores per node
#SBATCH --mem=20G                   #60G for hisat-build  #Request 8GB per node 
#SBATCH --output=terra_out/getNetworkTAM111up_A_out.%j      #Send stdout/err to "Example2Out.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --account=122825952664             #Set billing account to 123456
#SBATCH --partition=knl
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=r.kapoor@tamu.edu    #Send all emails to email_address

#module load R/3.6.0-iomkl-2018b-recommended-mt
#module load R_tamu/4.0.3-foss-2020a-recommended-mt
#module load R_tamu/3.6.0-iomkl-2018b-recommended-mt
module load R_tamu/3.6.3-foss-2020a
Rscript terra_cluster_JGL_B.R
#Rscript Test.R