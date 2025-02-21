#!/bin/bash
#SBATCH --job-name=RNAseq_Project        # Job name
#SBATCH --partition=batch                # Partition (queue) name
#SBATCH --ntasks=4                       # Number of tasks
#SBATCH --cpus-per-task=4               # Number of cores per task
#SBATCH --mem=16gb                       # Total memory for job
#SBATCH --time=20:00:00                  # Time limit hrs:min:sec
#SBATCH --output=/path/to/log.%j         # Standard output and error log
# #SBATCH --mail-user=your_email@domain.com  # Uncomment and set if you want email notifications
# #SBATCH --mail-type=END,FAIL               # Mail events (BEGIN, END, FAIL, ALL)

# Load necessary modules
module load kallisto/0.43.1-foss-2016b
module load R/3.4.4-foss-2016b-X11-20160819-GACRC
module load SRA-Toolkit/2.9.1-centos_linux64

# Download RNA-seq data with fastq-dump (SRA Toolkit)
# Data source: NCBI BioProject PRJNA498535
for n in {10..21}; do
  fastq-dump --split-files --gzip SRR81132${n}
done

# Download Phaseolus vulgaris CDS from the JGI Phytozome portal
# Make sure you have valid credentials to log in
curl 'https://signon.jgi.doe.gov/signon/create' \
  --data-urlencode 'login=YourLoginName' \
  --data-urlencode 'password=YourPassword' \
  -c cookies > /dev/null

curl "https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Pvulgaris/download/_JAMO/57fecb5f7ded5e3135bc3576/Pvulgaris_442_v2.1.cds.fa.gz" \
  -b cookies > Pvulgaris_cds.fa.gz

gunzip Pvulgaris_cds.fa.gz

# Build Kallisto index
kallisto index -i Pvulgaris_cds.index Pvulgaris_cds.fa

# Quantify reads using Kallisto
mkdir output

for i in {10..21}; do
  dir="kallisto_SRR81132${i}"
  file1="SRR81132${i}_1.fastq.gz"
  file2="SRR81132${i}_2.fastq.gz"
  kallisto quant -i Pvulgaris_cds.index -o "$dir" "$file1" "$file2" -b 100 -t 4
done

mv kallisto_SRR81132* output/

# Create a directory for Sleuth results
mkdir results

# Run the Sleuth R script
Rscript sleuth.R
