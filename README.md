# Genome_annotation
Genome annotation done with Demi (my script and notes, we were annotation Barbus barbus genome)

Anotace genomu

Raw file directory: /storage/praha1/home/eliasok1/fishevo/data_source/barbus
Working directory: /storage/brno12-cerit/home/konvickz/barbus/genome_annotation

Konverze složky SRA na fastq:

Jeden transkriptom byl stažen z NCBI databáze jako SRA, musí převeden na fastq a zazipován.
Vytvořit skript:
  nano convert.pbs

#!/bin/bash
#PBS -l select=1:ncpus=1:mem=12gb:scratch_local=12gb
#PBS -l walltime=01:00:00
cd /storage/praha1/home/eliasok1/fishevo/data_source/barbus/transcriptome_genbank/ 
module add conda-modules-py37
conda activate sra-tools-3.0.3
fastq-dump --split-3 ERR10123689.sra > barbus_spleen.fastq

Následně zazipovat:

    gzip barbus_spleen.fastq

Skript nakonec vytvořil složku barbus_spleen.fastq (tuto nepotřebujeme, obsahuje asi nějaký report), a dvě složky ERR10123689_1.fastq, ERR10123689_2.fastq, které obsahují naše forward a reverse reads.

Filtering and cleaning the reads before mapping

https://github.com/OpenGene/fastp
https://wiki.metacentrum.cz/wiki/Fastp


#!/bin/bash
#PBS -l select=1:ncpus=6:mem=12gb:scratch_local=12gb
#PBS -l walltime=04:00:00
cd /storage/praha1/home/eliasok1/fishevo/data_source/barbus/RNA_raw_data/
module add conda-modules-py37
conda activate fastp
output_dir=/storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files
for forward_read in *_1.fq.gz; do
    base=$(basename "$forward_read" _1.fq.gz)
    reverse_read="${base}_2.fq.gz"
    output_forward="$output_dir/${base}_fastp_1.fq.gz"
    output_reverse="$output_dir/${base}_fastp_2.fq.gz"
        fastp -i "$forward_read" -I "$reverse_read" \
          -o "$output_forward" -O "$output_reverse" -l 50 --threads 6\
          --html "$output_dir/${base}_fastp.html" --json "$output_dir/${base}_fastp.json"
          
    echo "Processed $base"
done

cd /storage/praha1/home/eliasok1/fishevo/data_source/barbus/RNA_raw_data/

output_dir=/storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files
for forward_read in *_1.fq.gz; do
    base=$(basename "$forward_read" _1.fq.gz)
    reverse_read="${base}_2.fq.gz"
    output_forward="$output_dir/${base}_fastp_1.fq.gz"
    output_reverse="$output_dir/${base}_fastp_2.fq.gz"
    echo "Output Forward File: $output_forward"
    echo "Output Reverse File: $output_reverse"
done

Test for my loop with file names:
fastp -i /path/to/test_1.fq.gz -I /path/to/test_2.fq.gz \
      -o /path/to/output/test_fastp_1.fq.gz -O /path/to/output/test_fastp_2.fq.gz
Script from Demi to process one sample first:

/storage/brno1-cerit/home/dburghez/BAROMBI/WGS/bin/fastp --in1 $path/RAW/$indiv.R1.fastq.gz --in2 $path/RAW/$indiv.R2.fastq.gz --ou
t1 $path/FILTERED/$indiv.R1.trimmed.fastq.gz --out2 $path/FILTERED/$indiv.R2.trimmed.fastq.gz --unpaired1 $path/FILTERED/$indiv.R1.
unpaired.fastq.gz --unpaired2 $path/FILTERED/$indiv.R2.unpaired.fastq.gz -l 50 --thread 6 -h $path/HTML/$indiv.html &> $path/HTML/$
indiv.log\n\n"


#!/bin/bash
#PBS -l select=1:ncpus=6:mem=12gb:scratch_local=12gb
#PBS -l walltime=04:00:00
cd /storage/praha1/home/eliasok1/fishevo/data_source/barbus/RNA_raw_data/
module add conda-modules-py37
conda activate fastp
fastp \
  --in1 Barb_104B7_1.fq.gz \
  --in2 Barb_104B7_2.fq.gz \
  --out1 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_1.fq.gz \
  --out2 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_2.fq.gz \
  -l 50 \
  --thread 6 \
  -h /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/HTML/Barb_104B7.html



#!/bin/bash
#PBS -l select=1:ncpus=1:mem=12gb:scratch_local=12gb
#PBS -l walltime=01:00:00

cd /storage/praha1/home/eliasok1/fishevo/data_source/barbus/transcriptome_genbank/ module add conda-modules-py37
conda activate fastp
fastp -i in.*1.fq.gz -I in.*2.fq.gz \
-l 50 –-threads 6 
-o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/out.*1.fq.gz \
-O /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/out.*2.fq.gz


#!/bin/bash

# Input directory containing raw fastq files
input_dir="/path/to/input_files"
# Output directory to store processed files
output_dir="/path/to/output_files"
for forward_read in "$input_dir"/*_1.fq.gz; do
    base=$(basename "$forward_read" _1.fq.gz)
    reverse_read="$input_dir/${base}_2.fq.gz"
    output_forward="$output_dir/${base}_fastp_1.fq.gz"
    output_reverse="$output_dir/${base}_fastp_2.fq.gz"
        fastp -i "$forward_read" -I "$reverse_read" \
          -o "$output_forward" -O "$output_reverse" \
          --html "$output_dir/${base}_fastp.html" --json "$output_dir/${base}_fastp.json"
          
    echo "Processed $base"
done

for f in in.*1.fq.gz; do
  r=${f%1.fq.gz}2.fq.gz  # This gets the corresponding reverse read file
  base=$(basename "$f" 1.fq.gz)  # Extract the base name
  fastp -i "$f" -I "$r" \
        -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/out."$base"1.fq.gz \
        -O /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/out."$base"2.fq.gz
done

Downloading HTML files after fastp:
scp -r konvickz@skirit.ics.muni.cz:/storage/brno12-cerit/home/konvickz/ barbus/genome_annotation/fastp_files/HTML .


Adaptovat původní skript pro první vzorek pro všechny vzorky (vytvoří samostatné skripty)
sed 's/104B7/14d10/g' fastp_1sample > fastp_14d10.pbs

Genome indexing (with HISAT2)

module load hisat2/2


$ hisat2-build -p 16 genome.fa genome

#!/bin/bash
#PBS -l select=1:ncpus=6:mem=24gb:scratch_local=24gb
#PBS -l walltime=24:00:00

cd /storage/praha1/home/eliasok1/fishevo/data_source/barbus/genome_raw_data
module add hisat2-2.2.1
hisat2-build -p 6 Barbus_barbusA_masked.fa /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/genome_indexed/Barbus_barbusA_index

#!/bin/bash
#PBS -l select=1:ncpus=6:mem=24gb:scratch_local=24gb
#PBS -l walltime=24:00:00

cd /storage/praha1/home/eliasok1/fishevo/data_source/barbus/genome_raw_data
module add hisat2-2.2.1
hisat2-build -p 6 Barbus_barbus_chrom_masked.fa /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/genome_indexed/Barbus_barbus_chrom_index



Barbus_barbusB_masked.fa
Barbus_barbus_chrom_masked.fa


/storage/brno12-cerit/home/konvickz/barbus/genome_annotation/genome_indexed


Mapping and conversion to bam file (mapping always produces sam files)

#!/bin/bash
#PBS -l select=1:ncpus=8:mem=24gb:scratch_local=24gb
#PBS -l walltime=24:00:00

export OMP_NUM_THREADS=7

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/genome_indexed

module add hisat2-2.2.1

hisat2 -x Barbus_barbusA_index -1 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_1.fq.gz -2 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_2.fq.gz -S /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/sam_mapped/Barbus_barbusA_104B7_mapped.sam --threads 7 --no-unal

module add samtools-1.11

samtools sort -O bam -@ 7 -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/sam_mapped/Barbus_barbusA_104B7_mapped.bam\
 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/sam_mapped/Barbus_barbusA_104B7_mapped.sam

Adapting script with sed command again:
sed 's/104B7/roz1/g' map_A_104B7.pbs > map_A_roz1.pbs


For B chromosomes and all chromosomes:

#!/bin/bash
#PBS -l select=1:ncpus=8:mem=24gb:scratch_local=24gb
#PBS -l walltime=24:00:00

export OMP_NUM_THREADS=7

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/genome_indexed

module add hisat2-2.2.1

hisat2 -x Barbus_barbus_chrom_index -1 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_1.fq.gz -2 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_2.fq.gz -S /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/sam_mapped/Barbus_barbus_chrom_104B7_mapped.sam --threads 7 --no-unal --dta

module add samtools-1.11

samtools sort -O bam -@ 7 -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/sam_mapped/Barbus_barbus_chrom_104B7_mapped.bam\
 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/sam_mapped/Barbus_barbus_chrom_104B7_mapped.sam 


StringTie

#!/bin/bash
#PBS -l select=1:ncpus=4:mem=12gb:scratch_local=12gb
#PBS -l walltime=2:00:00

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/sam_mapped

module add stringtie-2.1.0

stringtie -p 4 -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/stringtie/Barbus_barbus_chrom_104B7.gtf Barbus_barbus_chrom_104B7_mapped.bam

Adapting script with sed command again:
sed 's/104B7/roz1/g' str_chrom_104B7.pbs > str_chrom_roz1.pbs

Downloading Integrative Genomics Viewer
https://igv.org/doc/desktop/#DownloadPage/

Downloading gtf files into my laptop
•	open new terminal
•	go to the file where I want gtf files (in my laptop)
•	log into metacentrum using sftp konvickz@skirit.ics.muni.cz (instead of ssh)
•	then get *.gtf
•	
•	
Assessing mapping results in IGV
•	download genomes (masked) into my laptop
•	download gtf files produced in Stringtie step
•	download gff file produced in masking step (done by Demi before we started this protocol)
•	upload everything into IGV

Merging gtf files with Stringtie

Barbus_barbusA_104B7.gtf
Barbus_barbusA_21d10.gtf
Barbus_barbusA_28d2.gtf
Barbus_barbusA_6wk4.gtf
Barbus_barbusA_roz2.gtf
Barbus_barbusA_14d10.gtf
Barbus_barbusA_21d1.gtf
Barbus_barbusA_28d3.gtf
Barbus_barbusA_6wk5.gtf
Barbus_barbusA_roz3.gtf
Barbus_barbusA_14d1.gtf
Barbus_barbusA_21d2.gtf
Barbus_barbusA_28d4.gtf
Barbus_barbusA_91B8.gtf
Barbus_barbusA_roz5.gtf
Barbus_barbusA_14d2.gtf
Barbus_barbusA_21d3.gtf
Barbus_barbusA_6wk1.gtf
Barbus_barbusA_box1B2.gtf	
Barbus_barbusA_spleen.gtf
Barbus_barbusA_14d3.gtf
Barbus_barbusA_21d4.gtf
Barbus_barbusA_6wk2.gtf
Barbus_barbusA_box5H9.gtf
Barbus_barbusA_14d4.gtf
Barbus_barbusA_28d1.gtf
Barbus_barbusA_6wk3.gtf
Barbus_barbusA_roz1.gtf

#!/bin/bash
#PBS -l select=1:ncpus=4:mem=12gb:scratch_local=12gb
#PBS -l walltime=2:00:00

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/stringtie

module add stringtie-2.1.0

stringtie --merge -p 4 -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/merged/Barbus_barbus_chrom_merged.gtf /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/scripts/merge_str/chrom_merge_list.txt

sed 's/barbusA/barbusB/g' A_merge_list.txt > B_merge_list.txt
Concatenate A and B merged gtf files:

cat fileA fileB >> fileC

cat Barbus_barbus_A_merged.gtf Barbus_barbus_B_merged.gtf > Barbus_barbus_merged_ABcombined.gtf


cp /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/merged/Barbus_barbus_chrom_merged.gtf /storage/praha1/home/eliasok1/fishevo/data_source/barbus

Downloading files from online databases to Metacentrum
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/012/432/095/GCF_012432095.1_ASM1243209v1/GCF_012432095.1_ASM1243209v1_protein.faa.gz


Giving permission to work with my file to all Metacentrum users
chmod 777 GCF_012432095.1_ASM1243209v1_protein.faa
chmod 777 -R directory name (for directory and everything in it, R as recursive)

Remaping with HISAT2 to get even cleaner result:

#PBS -l select=1:ncpus=8:mem=24gb:scratch_local=24gb
#PBS -l walltime=24:00:00

export OMP_NUM_THREADS=7

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/genome_indexed

module add hisat2-2.2.1

hisat2 -x Barbus_barbus_chrom_index -1 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_1.fq.gz\
 -2 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/fastp_files/Barb_104B7_fastp_2.fq.gz\
 -S /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files/Barbus_barbus_chrom_104B7_map.sam\
 --threads 7 --no-unal --dta --no-mixed --no-discordant 

module add samtools-1.11

samtools sort -O bam -@ 7 -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files/Barbus_barbus_chrom_104B7_map.bam\
 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files/Barbus_barbus_chrom_104B7_map.sam

samtools view -F3844 -q 60 /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files/Barbus_barbus_chrom_104B7_map.bam -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files/Barbus_barbus_chrom_104B7_filt.bam



sed 's/104B7/roz1/g' map_filt_chrom_104B7.pbs > map_filt_chrom_roz1.pbs


Merging the bam files before submitting to Stringtie with samtools cat
We will merge all A and B files together and chrom files separately (i.e. two outputs)



#!/bin/bash
#PBS -l select=1:ncpus=8:mem=50gb:scratch_local=50gb
#PBS -l walltime=24:00:00


module add samtools-1.11

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files

samtools cat -b list_AB_filt.txt -o Barbus_AB_samtools_merged.bam


#!/bin/bash
#PBS -l select=1:ncpus=8:mem=50gb:scratch_local=50gb
#PBS -l walltime=24:00:00


module add samtools-1.11

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files

samtools cat -b list_chrom_filt.txt -o Barbus_chrom_samtools_merged.bam

StringTie

#!/bin/bash
#PBS -l select=1:ncpus=8:mem=24gb:scratch_local=24gb
#PBS -l walltime=24:00:00

cd /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/map_filtered_files

module add stringtie-2.1.0

stringtie -p 8 -o /storage/brno12-cerit/home/konvickz/barbus/genome_annotation/stringtie/Barbus_chrom_samtools_merged.gtf Barbus_chrom_samtools_merged.bam


Making list:
Barbus_barbusA_104B7_filt.bam
Barbus_barbusA_14d4_filt.bam
Barbus_barbusA_21d4_filt.bam
Barbus_barbusA_6wk1_filt.bam
Barbus_barbusA_91B8_filt.bam	
Barbus_barbusA_roz3_filt.bam
Barbus_barbusA_14d10_filt.bam
Barbus_barbusA_21d10_filt.bam
Barbus_barbusA_28d1_filt.bam
Barbus_barbusA_6wk2_filt.bam
Barbus_barbusA_box1B2_filt.bam
Barbus_barbusA_roz5_filt.bam
Barbus_barbusA_14d1_filt.bam
Barbus_barbusA_21d1_filt.bam
Barbus_barbusA_28d2_filt.bam
Barbus_barbusA_6wk3_filt.bam
Barbus_barbusA_box5H9_filt.bam
Barbus_barbusA_spleen_filt.bam
Barbus_barbusA_14d2_filt.bam
Barbus_barbusA_21d2_filt.bam
Barbus_barbusA_28d3_filt.bam
Barbus_barbusA_6wk4_filt.bam
Barbus_barbusA_roz1_filt.bam
Barbus_barbusA_14d3_filt.bam
Barbus_barbusA_21d3_filt.bam
Barbus_barbusA_28d4_filt.bam
Barbus_barbusA_6wk5_filt.bam
Barbus_barbusA_roz2_filt.bam

sed 's/barbusA/barbus_chrom/g' list_AB_filt.txt > list_chrom_filt.txt
sed 's/barbusA/barbusB/g' list_AB_filt.txt >> list_AB_filt.txt

