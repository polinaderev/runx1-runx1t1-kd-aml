#Concatenate fastq files from 2 sequencing runs

# cat A_01_EKDL220016748-1A_HKKNCDSX5_L3_1.fq.gz A_01_EKDL220016748-1A_HKKNFDSX5_L4_1.fq.gz > siRE_rep1_R1.fastq.gz
# cat A_01_EKDL220016748-1A_HKKNCDSX5_L3_2.fq.gz A_01_EKDL220016748-1A_HKKNFDSX5_L4_2.fq.gz > siRE_rep1_R2.fastq.gz
# cat A_02_EKDL220016749-1A_HKJK3DSX5_L2_1.fq.gz A_02_EKDL220016749-1A_HKKNFDSX5_L2_1.fq.gz > siRE_rep2_R1.fastq.gz
# cat A_02_EKDL220016749-1A_HKJK3DSX5_L2_2.fq.gz A_02_EKDL220016749-1A_HKKNFDSX5_L2_2.fq.gz > siRE_rep2_R2.fastq.gz
# cat A_03_EKDL220016750-1A_HJTWJDSX5_L2_1.fq.gz A_03_EKDL220016750-1A_HKKNFDSX5_L3_1.fq.gz > siRE_rep3_R1.fastq.gz
# cat A_03_EKDL220016750-1A_HJTWJDSX5_L2_2.fq.gz A_03_EKDL220016750-1A_HKKNFDSX5_L3_2.fq.gz > siRE_rep3_R2.fastq.gz
# cat A_04_EKDL220016751-1A_HJTWJDSX5_L3_1.fq.gz A_04_EKDL220016751-1A_HKKTGDSX5_L1_1.fq.gz > siMM_rep1_R1.fastq.gz
# cat A_04_EKDL220016751-1A_HJTWJDSX5_L3_2.fq.gz A_04_EKDL220016751-1A_HKKTGDSX5_L1_2.fq.gz > siMM_rep1_R2.fastq.gz
# cat A_05_EKDL220016752-1A_HJTH2DSX5_L3_1.fq.gz > siMM_rep2_R1.fastq.gz
# cat A_05_EKDL220016752-1A_HJTH2DSX5_L3_2.fq.gz > siMM_rep2_R2.fastq.gz
# cat A_06_EKDL220016753-1A_HJTH2DSX5_L4_1.fq.gz A_06_EKDL220016753-1A_HJTWJDSX5_L4_1.fq.gz > siMM_rep3_R1.fastq.gz
# cat A_06_EKDL220016753-1A_HJTH2DSX5_L4_2.fq.gz A_06_EKDL220016753-1A_HJTWJDSX5_L4_2.fq.gz > siMM_rep3_R2.fastq.gz



#Trim adapters

module load Trimmomatic/0.39-Java-11

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 \
fastq/siRE_rep1_R1.fastq.gz fastq/siRE_rep1_R2.fastq.gz \
fastq/siRE_rep1_trimmed_paired_R1.fastq.gz fastq/siRE_rep1_trimmed_unpaired_R1.fastq.gz \
fastq/siRE_rep1_trimmed_paired_R2.fastq.gz fastq/siRE_rep1_trimmed_unpaired_R2.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 \
fastq/siRE_rep2_R1.fastq.gz fastq/siRE_rep2_R2.fastq.gz \
fastq/siRE_rep2_trimmed_paired_R1.fastq.gz fastq/siRE_rep2_trimmed_unpaired_R1.fastq.gz \
fastq/siRE_rep2_trimmed_paired_R2.fastq.gz fastq/siRE_rep2_trimmed_unpaired_R2.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 \
fastq/siRE_rep3_R1.fastq.gz fastq/siRE_rep3_R2.fastq.gz \
fastq/siRE_rep3_trimmed_paired_R1.fastq.gz fastq/siRE_rep3_trimmed_unpaired_R1.fastq.gz \
fastq/siRE_rep3_trimmed_paired_R2.fastq.gz fastq/siRE_rep3_trimmed_unpaired_R2.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 \
fastq/siMM_rep1_R1.fastq.gz fastq/siMM_rep1_R2.fastq.gz \
fastq/siMM_rep1_trimmed_paired_R1.fastq.gz fastq/siMM_rep1_trimmed_unpaired_R1.fastq.gz \
fastq/siMM_rep1_trimmed_paired_R2.fastq.gz fastq/siMM_rep1_trimmed_unpaired_R2.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 \
fastq/siMM_rep2_R1.fastq.gz fastq/siMM_rep2_R2.fastq.gz \
fastq/siMM_rep2_trimmed_paired_R1.fastq.gz fastq/siMM_rep2_trimmed_unpaired_R1.fastq.gz \
fastq/siMM_rep2_trimmed_paired_R2.fastq.gz fastq/siMM_rep2_trimmed_unpaired_R2.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 \
fastq/siMM_rep3_R1.fastq.gz fastq/siMM_rep3_R2.fastq.gz \
fastq/siMM_rep3_trimmed_paired_R1.fastq.gz fastq/siMM_rep3_trimmed_unpaired_R1.fastq.gz \
fastq/siMM_rep3_trimmed_paired_R2.fastq.gz fastq/siMM_rep3_trimmed_unpaired_R2.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50


#Align

mkdir alignments

module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.12-GCC-10.3.0

bowtie2 -p 8 --very-sensitive-local -x Genome_indexes/Human/hg38/Bowtie2/hg38 -1 fastq/siRE_rep1_trimmed_paired_R1.fastq.gz -2 fastq/siRE_rep1_trimmed_paired_R2.fastq.gz | \
samtools sort -@ 8 -T siRE_rep1 -o alignments/siRE_rep1.bam

bowtie2 -p 8 --very-sensitive-local -x Genome_indexes/Human/hg38/Bowtie2/hg38 -1 fastq/siRE_rep2_trimmed_paired_R1.fastq.gz -2 fastq/siRE_rep2_trimmed_paired_R2.fastq.gz | \
samtools sort -@ 8 -T siRE_rep2 -o alignments/siRE_rep2.bam

bowtie2 -p 8 --very-sensitive-local -x Genome_indexes/Human/hg38/Bowtie2/hg38 -1 fastq/siRE_rep3_trimmed_paired_R1.fastq.gz -2 fastq/siRE_rep3_trimmed_paired_R2.fastq.gz | \
samtools sort -@ 8 -T siRE_rep3 -o alignments/siRE_rep3.bam

bowtie2 -p 8 --very-sensitive-local -x Genome_indexes/Human/hg38/Bowtie2/hg38 -1 fastq/siMM_rep1_trimmed_paired_R1.fastq.gz -2 fastq/siMM_rep1_trimmed_paired_R2.fastq.gz | \
samtools sort -@ 8 -T siMM_rep1 -o alignments/siMM_rep1.bam

bowtie2 -p 8 --very-sensitive-local -x Genome_indexes/Human/hg38/Bowtie2/hg38 -1 fastq/siMM_rep2_trimmed_paired_R1.fastq.gz -2 fastq/siMM_rep2_trimmed_paired_R2.fastq.gz | \
samtools sort -@ 8 -T siMM_rep2 -o alignments/siMM_rep2.bam

bowtie2 -p 8 --very-sensitive-local -x Genome_indexes/Human/hg38/Bowtie2/hg38 -1 fastq/siMM_rep3_trimmed_paired_R1.fastq.gz -2 fastq/siMM_rep3_trimmed_paired_R2.fastq.gz | \
samtools sort -@ 8 -T siMM_rep3 -o alignments/siMM_rep3.bam

#Removeduplicates

module load picard/2.21.1-Java-11

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=alignments/siRE_rep1.bam O=alignments/siRE_rep1_rmdup.bam M=alignments/siRE_rep1_metrics.txt AS=true REMOVE_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=alignments/siRE_rep2.bam O=alignments/siRE_rep2_rmdup.bam M=alignments/siRE_rep2_metrics.txt AS=true REMOVE_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=alignments/siRE_rep3.bam O=alignments/siRE_rep3_rmdup.bam M=alignments/siRE_rep3_metrics.txt AS=true REMOVE_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=alignments/siMM_rep1.bam O=alignments/siMM_rep1_rmdup.bam M=alignments/siMM_rep1_metrics.txt AS=true REMOVE_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=alignments/siMM_rep2.bam O=alignments/siMM_rep2_rmdup.bam M=alignments/siMM_rep2_metrics.txt AS=true REMOVE_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=alignments/siMM_rep3.bam O=alignments/siMM_rep3_rmdup.bam M=alignments/siMM_rep3_metrics.txt AS=true REMOVE_DUPLICATES=true

samtools index alignments/siRE_rep1_rmdup.bam
samtools index alignments/siRE_rep2_rmdup.bam
samtools index alignments/siRE_rep3_rmdup.bam
samtools index alignments/siMM_rep1_rmdup.bam
samtools index alignments/siMM_rep2_rmdup.bam
samtools index alignments/siMM_rep3_rmdup.bam

#Genome browser files

module purge; module load bluebear
module load deepTools/3.5.0-foss-2020a-Python-3.8.2

bamCoverage -b alignments/siRE_rep1_rmdup.bam -o alignments/siRE_rep1.bw --normalizeUsing CPM
bamCoverage -b alignments/siRE_rep2_rmdup.bam -o alignments/siRE_rep2.bw --normalizeUsing CPM
bamCoverage -b alignments/siRE_rep3_rmdup.bam -o alignments/siRE_rep3.bw --normalizeUsing CPM
bamCoverage -b alignments/siMM_rep1_rmdup.bam -o alignments/siMM_rep1.bw --normalizeUsing CPM
bamCoverage -b alignments/siMM_rep2_rmdup.bam -o alignments/siMM_rep2.bw --normalizeUsing CPM
bamCoverage -b alignments/siMM_rep3_rmdup.bam -o alignments/siMM_rep3.bw --normalizeUsing CPM

#Peakcalling

module purge; module load bluebear
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

mkdir Peaks

macs2 callpeak -t alignments/siRE_rep1_rmdup.bam -n siRE_rep1 --outdir Peaks -q 0.0005 -B --trackline --nomodel --shift -100 --extsize 200

macs2 callpeak -t alignments/siRE_rep2_rmdup.bam -n siRE_rep2 --outdir Peaks -q 0.0005 -B --trackline --nomodel --shift -100 --extsize 200

macs2 callpeak -t alignments/siRE_rep3_rmdup.bam -n siRE_rep3 --outdir Peaks -q 0.0005 -B --trackline --nomodel --shift -100 --extsize 200

macs2 callpeak -t alignments/siMM_rep1_rmdup.bam -n siMM_rep1 --outdir Peaks -q 0.0005 -B --trackline --nomodel --shift -100 --extsize 200

macs2 callpeak -t alignments/siMM_rep2_rmdup.bam -n siMM_rep2 --outdir Peaks -q 0.0005 -B --trackline --nomodel --shift -100 --extsize 200

macs2 callpeak -t alignments/siMM_rep3_rmdup.bam -n siMM_rep3 --outdir Peaks -q 0.0005 -B --trackline --nomodel --shift -100 --extsize 200

#Filter for blacklist peaks

module purge; module load bluebear
module load BEDTools/2.29.2-GCC-9.3.0

bedtools intersect -a Peaks/siRE_rep1_peaks.narrowPeak -b Annotations/hg38-blacklist.v3.bed.gz -v > Peaks/siRE_rep1_filtered.bed
bedtools intersect -a Peaks/siRE_rep2_peaks.narrowPeak -b Annotations/hg38-blacklist.v3.bed.gz -v > Peaks/siRE_rep2_filtered.bed
bedtools intersect -a Peaks/siRE_rep3_peaks.narrowPeak -b Annotations/hg38-blacklist.v3.bed.gz -v > Peaks/siRE_rep3_filtered.bed

bedtools intersect -a Peaks/siMM_rep1_peaks.narrowPeak -b Annotations/hg38-blacklist.v3.bed.gz -v > Peaks/siMM_rep1_filtered.bed
bedtools intersect -a Peaks/siMM_rep2_peaks.narrowPeak -b Annotations/hg38-blacklist.v3.bed.gz -v > Peaks/siMM_rep2_filtered.bed
bedtools intersect -a Peaks/siMM_rep3_peaks.narrowPeak -b Annotations/hg38-blacklist.v3.bed.gz -v > Peaks/siMM_rep3_filtered.bed



