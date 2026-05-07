#!/bin/bash

# ==========================================
# RNA-seq preprocessing workflow
# E. coli pressure response project
# ==========================================

# ------------------------------------------
# Required files
# ------------------------------------------
# raw fastq files
# reference_genome.fa
# annotation.gtf

# ------------------------------------------
# 1. Create directories
# ------------------------------------------

mkdir -p raw_fastq
mkdir -p trimmed_fastq
mkdir -p star_index
mkdir -p star_alignment
mkdir -p counts

# ------------------------------------------
# 2. Download RNA-seq data from SRA
# ------------------------------------------

prefetch SRR7217923
prefetch SRR7217925
prefetch SRR7217927
prefetch SRR7217929

# ------------------------------------------
# 3. Convert SRA to FASTQ
# ------------------------------------------

for sample in SRR7217923 SRR7217925 SRR7217927 SRR7217929
do

fasterq-dump ${sample} \
-O raw_fastq

done

# ------------------------------------------
# 4. Quality trimming with Trimmomatic
# ------------------------------------------

for sample in SRR7217923 SRR7217925 SRR7217927 SRR7217929
do

java -jar trimmomatic.jar PE \
-threads 4 \
-phred33 \
raw_fastq/${sample}_1.fastq \
raw_fastq/${sample}_2.fastq \
trimmed_fastq/${sample}_paired_1.fastq \
trimmed_fastq/${sample}_unpaired_1.fastq \
trimmed_fastq/${sample}_paired_2.fastq \
trimmed_fastq/${sample}_unpaired_2.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:20 \
TRAILING:20 \
SLIDINGWINDOW:4:20 \
MINLEN:50

done

# ------------------------------------------
# 5. Build STAR genome index
# ------------------------------------------

STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir star_index \
--genomeFastaFiles reference_genome.fa

# ------------------------------------------
# 6. Map reads using STAR
# ------------------------------------------

for sample in SRR7217923 SRR7217925 SRR7217927 SRR7217929
do

STAR \
--runThreadN 4 \
--genomeDir star_index \
--readFilesIn \
trimmed_fastq/${sample}_paired_1.fastq \
trimmed_fastq/${sample}_paired_2.fastq \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix star_alignment/${sample}_

done

# ------------------------------------------
# 7. Generate count matrix using mmquant
# ------------------------------------------

for sample in SRR7217923 SRR7217925 SRR7217927 SRR7217929
do

mmquant \
-s U \
-f BAM \
-t 1 \
-a annotation.gtf \
-r star_alignment/${sample}_Aligned.sortedByCoord.out.bam \
-o counts/${sample}_mmquant.out

done

# ------------------------------------------
# Final output
# ------------------------------------------

# tagcount.out

echo "RNA-seq preprocessing completed"
