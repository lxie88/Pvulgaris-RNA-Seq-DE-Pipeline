#!/bin/bash
#SBATCH --job-name=Finalproject_momo# Job name
#SBATCH --partition=batch # Partition (queue) name
#SBATCH --ntasks=4 # Single task job
#SBATCH --cpus-per-task=4 # Number of cores per task
#SBATCH --mem=16gb # Total memory for job
#SBATCH --time=20:00:00 # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/lx83659/FinalProject/log.%j # Standard output and error log
#SBATCH --mail-user=lx83659@uga.edu # Where to send mail
#SBATCH --mail-type=END,FAIL # Mail events (BEGIN, END, FAIL, ALL)

# Download necessary modules 

module load kallisto/0.43.1-foss-2016b
module load R/3.4.4-foss-2016b-X11-20160819-GACRC
module load SRA-Toolkit/2.9.1-centos_linux64

# Downloads rna-seq reads using fastq-dump (SRA-Toolkit/2.9.1-centos_linux64)
# You need to figure out which sequence you would like to download
# download RNAseq data from https://www. ncbi.nlm.nih.gov/sra/PRJNA498535

for n in {10..21};
do  fastq-dump --split-files --gzip SRR81132${n}
done

# download Phaselous vulgaris cds sequence from Phytozome v10 https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Pvulgaris

# login to the Phytozome website 
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=momosan@uga.edu' --data-urlencode 'password=077338Frank' -c cookies > /dev/null

# download cds sequence 
curl "https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Pvulgaris/download/_JAMO/57fecb5f7ded5e3135bc3576/Pvulgaris_442_v2.1.cds.fa.gz" -b cookies > Pvulgaris_cds.fa.gz
gunzip Pvulgaris_cds.fa.gz 

# Makes a kallisto index for cds sequence 
kallisto index -i Pvulgaris_cds.index Pvulgaris_cds.fa    

# use kallisto quantify the reads 
mkdir output
for i in {10..21};
do 
	dir="kallisto_SRR81132${i}" 
	file1="SRR81132${i}_1.fastq.gz"
        file2="SRR81132${i}_2.fastq.gz"
  	kallisto quant -i Pvulgaris_cds.index -o $dir $file1 $file2 -b 100 -t 4
done

mv kallisto_SRR81132* output/

# create a directory to store slueth results 
mkdir results 

# start running slueth R commonbean script  

Rscript sleuth_commonbean.R




