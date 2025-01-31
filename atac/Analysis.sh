#Make union and annotate to compare replicates

cat Peaks/siRE_rep1_filtered.bed Peaks/siRE_rep2_filtered.bed Peaks/siRE_rep3_filtered.bed Peaks/siMM_rep1_filtered.bed Peaks/siMM_rep2_filtered.bed Peaks/siMM_rep3_filtered.bed | bedtools merge -i stdin > Peaks/Union.bed

module load HOMER/4.11-foss-2019b

annotatePeaks.pl Peaks/Union.bed hg38 -bedGraph Peaks/siRE_rep1_treat_pileup.bdg Peaks/siRE_rep2_treat_pileup.bdg Peaks/siRE_rep3_treat_pileup.bdg Peaks/siMM_rep1_treat_pileup.bdg Peaks/siMM_rep2_treat_pileup.bdg Peaks/siMM_rep3_treat_pileup.bdg -size 200 > Peaks/Union_annotated.bed

#RScript comparereplicates.r

#Make peak set from only peaks shared in all 3 replicates

module load BEDTools/2.29.2-GCC-9.3.0

bedtools intersect -a Peaks/siRE_rep1_filtered.bed -b Peaks/siRE_rep2_filtered.bed -u | bedtools intersect -a stdin -b Peaks/siRE_rep3_filtered.bed -u > Peaks/siRE_shared_peakset.bed

bedtools intersect -a Peaks/siMM_rep1_filtered.bed -b Peaks/siMM_rep2_filtered.bed -u | bedtools intersect -a stdin -b Peaks/siMM_rep3_filtered.bed -u > Peaks/siMM_shared_peakset.bed

cat Peaks/siRE_shared_peakset.bed Peaks/siMM_shared_peakset.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > Peaks/Union_filteredforreplicates.bed

mkdir Counts

awk 'OFS="\t" {print $1"."$2+1"."$3, $1, $2+1, $3, "."}' Peaks/Union_filteredforreplicates.bed > Counts/Union.saf

#Get counts for all peaks

module load Subread/2.0.1-GCC-8.3.0

featureCounts -a Counts/Union.saf -F SAF -T 8 -p -o Counts/Counts.tsv alignments/siRE_rep1_rmdup.bam alignments/siRE_rep2_rmdup.bam alignments/siRE_rep3_rmdup.bam alignments/siMM_rep1_rmdup.bam alignments/siMM_rep2_rmdup.bam alignments/siMM_rep3_rmdup.bam

#RScript DApeaks.r

#Annotate to get distance from TSS

module load HOMER/4.11-foss-2019b

annotatePeaks.pl Analysis/Ordered_allpeaks.txt hg38 > Analysis/Ordered_allpeaks_annnotated.tsv
annotatePeaks.pl Analysis/Lost_allpeaks.txt hg38 > Analysis/Lost_allpeaks_annnotated.tsv
annotatePeaks.pl Analysis/Gained_allpeaks.txt hg38 > Analysis/Gained_allpeaks_annnotated.tsv

#Separate out distal peaks

#Motif enrichment analysis

findMotifsGenome.pl Analysis/Lost_distalpeaks.txt hg38 Analysis/LostMotifs/ -size 200 -noknown -preparsedDir ../motifs
findMotifsGenome.pl Analysis/Gained_distalpeaks.txt hg38 Analysis/GainedMotifs/ -size 200 -noknown -preparsedDir ../motifs

#Make density plots

mkdir CDTfiles

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -bedGraph Peaks/siMM_rep1_treat_pileup.bdg -size 2000 -hist 10 -ghist > CDTfiles/siMM_rep1_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -bedGraph Peaks/siMM_rep2_treat_pileup.bdg -size 2000 -hist 10 -ghist > CDTfiles/siMM_rep2_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -bedGraph Peaks/siMM_rep3_treat_pileup.bdg -size 2000 -hist 10 -ghist > CDTfiles/siMM_rep3_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -bedGraph Peaks/siRE_rep1_treat_pileup.bdg -size 2000 -hist 10 -ghist > CDTfiles/siRE_rep1_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -bedGraph Peaks/siRE_rep2_treat_pileup.bdg -size 2000 -hist 10 -ghist > CDTfiles/siRE_rep2_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -bedGraph Peaks/siRE_rep3_treat_pileup.bdg -size 2000 -hist 10 -ghist > CDTfiles/siRE_rep3_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/ap1.motif -size 2000 -hist 10 -ghist > CDTfiles/ap1motif_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/ets.motif -size 2000 -hist 10 -ghist > CDTfiles/etsmotif_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/cebp.motif -size 2000 -hist 10 -ghist > CDTfiles/cebpmotif_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/runx.motif -size 2000 -hist 10 -ghist > CDTfiles/runxmotif_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/pu1.motif -size 2000 -hist 10 -ghist > CDTfiles/pu1motif_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/myb.motif -size 2000 -hist 10 -ghist > CDTfiles/mybmotif_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/nfe2.motif -size 2000 -hist 10 -ghist > CDTfiles/nfe2motif_distalordered.cdt

annotatePeaks.pl Analysis/Ordered_distalpeaks.txt hg38 -m ../motifs/gata.motif -size 2000 -hist 10 -ghist > CDTfiles/gatamotif_distalordered.cdt


#Compare lost and gained distal peaks to healthy cells

module load BEDTools/2.29.2-GCC-9.3.0

cat Corces_et_al_ATACseq/peak_union/Distal_sites.bed Analysis/Lost_distalpeaks.txt Analysis/Gained_distalpeaks.txt | sort -k1,1 -k2,2n | bedtools merge -i stdin > Analysis/Distalsites_Corcesunion.tsv

module purge; module load bluebear
module load HOMER/4.11-foss-2019b

annotatePeaks.pl Analysis/Distalsites_Corcesunion.tsv hg38 -bedGraph Peaks/siMM_rep1_treat_pileup.bdg Peaks/siMM_rep2_treat_pileup.bdg Peaks/siMM_rep3_treat_pileup.bdg Peaks/siRE_rep1_treat_pileup.bdg Peaks/siRE_rep2_treat_pileup.bdg Peaks/siRE_rep3_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/CLP_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/CMP_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/GMP_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/HSC_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/LMPP_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/MEP_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/Monocyte_treat_pileup.bdg Corces_et_al_ATACseq/CellType_peaks/MPP_treat_pileup.bdg -size 200 > Analysis/Distalsites_Corcesunion_annotated.bed

#RScript corcesheatmap.r

#Intersect cell type specific peaks for density plots
module load BEDTools/2.29.2-GCC-9.3.0

bedtools intersect -a Analysis/Ordered_distalpeaks.txt -b Corces_et_al_ATACseq/CellType_specific_peaks/Monocyte_specific.bed -c > Analysis/Ordered_distal_monocytehits.tsv

bedtools intersect -a Analysis/Ordered_distalpeaks.txt -b Corces_et_al_ATACseq/CellType_specific_peaks/GMP_specific.bed -c > Analysis/Ordered_distal_GMPhits.tsv

bedtools intersect -a Analysis/Ordered_distalpeaks.txt -b Corces_et_al_ATACseq/CellType_specific_peaks/MPP_specific.bed -c > Analysis/Ordered_distal_MPPhits.tsv

bedtools intersect -a Analysis/Ordered_distalpeaks.txt -b Corces_et_al_ATACseq/CellType_specific_peaks/HSC_specific.bed -c > Analysis/Ordered_distal_HSChits.tsv


#Make above files into cdt for JavaTreeView - 4th column becomes 3rd = value, 1st column peakID numbers 1:inf, 2nd column GWEIGHT = all 1

####Priming analysis####

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR175/009/SRR1752199/SRR1752199.fastq.gz -o SRR1752199_GSM1581820_CD34_PBSC_DNAse-Seq_Homo_sapiens_DNase-Hypersensitivity.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR639/003/SRR6391633/SRR6391633.fastq.gz -o SRR6391633_GSM2893606_ITD-1_DHS_Homo_sapiens_DNase-Hypersensitivity.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR639/004/SRR6391634/SRR6391634.fastq.gz -o SRR6391634_GSM2893607_ITD-2_DHS_Homo_sapiens_DNase-Hypersensitivity.fastq.gz

#Align hg38, remove duplicates, per Align.sh

samtools merge siMM_rmdup_merged.bam siMM_rep1_rmdup.bam siMM_rep2_rmdup.bam siMM_rep3_rmdup.bam
samtools merge siRE_rmdup_merged.bam ../alignments/siRE_rep1_rmdup.bam siRE_rep2_rmdup.bam siRE_rep3_rmdup.bam
samtools merge ITD-1-2.bam ITD-1.bam ITD-2.bam
samtools index siMM_rmdup_merged.bam
samtools index siRE_rmdup_merged.bam
samtools index ITD-1-2.bam

bamCoverage -b PBSC.bam -o PBSC.bw --normalizeUsing CPM
bamCoverage -b ITD-1-2.bam -o ITD-1-2.bw --normalizeUsing CPM
bamCoverage -b siMM_rmdup_merged.bam -o siMM.bw --normalizeUsing CPM
bamCoverage -b siRE_rmdup_merged.bam -o siRE.bw --normalizeUsing CPM

#Filter all sites associated with genes from Supercluster lists from Analysis/Ordered_allpeaks_annnotated.tsv

computeMatrix reference-point -S siMM.bw siRE.bw ITD-1-2.bw PBSC.bw -R Supercluster1sites.bed -a 1500 -b 1500 --referencePoint center -o Supercluster1.mat.gz
computeMatrix reference-point -S siMM.bw siRE.bw ITD-1-2.bw PBSC.bw -R Supercluster2sites.bed -a 1500 -b 1500 --referencePoint center -o Supercluster2.mat.gz
plotHeatmap -m Supercluster1.mat.gz -o Supercluster1.pdf
plotHeatmap -m Supercluster2.mat.gz -o Supercluster2.pdf

###CEBPE average profile####
#Make bed file of CEBPE enhancer of interest from Analysis/Ordered_allpeaks_annnotated.tsv

computeMatrix reference-point -S siMM.bw siRE.bw -R CEBPEenhancer.bed -a 1500 -b 1500 --referencePoint center -o temp.mat.gz --outFileNameMatrix cebpeenhancerforggplot.tsv

R

library(ggplot2)

enh <- read.delim("cebpeenhancerforggplot.tsv", sep="\t", header=TRUE)

ggplot(enh, aes(x=bin, y=Value, color=Cond)) + 
  geom_line(linewidth = 1.6) + 
  theme_classic() +
  scale_color_manual(values = c('#000000', '#f50000')) +
  labs(x = 'Distance from ATAC-seq peak summit (bp)', y = 'Counts Per Million (CPM)') +
  scale_x_discrete(limits = c(0,60,120), labels = c('-1.5kb', '0', '1.5kb')) +
  scale_y_continuous(limits = c(0,max(enh$Value))) + 
  theme(axis.text.x = element_text(size = 32), axis.title.x = element_text(size = 32)) +
  theme(axis.text.y = element_text(size = 32), axis.title.y = element_text(size = 32)) +
  theme(axis.line = element_line(size = 1.2), axis.ticks = element_line(size = 1.2), axis.ticks.length = unit(0.3,'cm'))
