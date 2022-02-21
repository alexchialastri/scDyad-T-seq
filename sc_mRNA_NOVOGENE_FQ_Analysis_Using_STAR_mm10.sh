#!/bin/bash                                                                                                                                                                                         
#PBS -l nodes=1:ppn=6                                                                                                                                                                           
#PBS -l walltime=48:00:00


#User Input
START_DIR="/home/achialastri/scSerum2i_RNA-TE/scTE_testing"



#Standard Usage per person
CEL_BARCODES="/home/achialastri/perlscripts/mRNA_Mapping/cel-seq_barcode_whitelist"
ANNO_GENOME_DIR="/home/achialastri/Genomes/mm10_STAR_Aligner_Genes"
ANNO_GENOME_DIR_TE="/home/achialastri/Genomes/mm10_STAR_Aligner_TE"
PERL_DIR="/home/achialastri/perlscripts/mRNA_Mapping"
HOME_DIR="/home/achialastri"
MAPPER_PATH="/home/achialastri/STAR/bin/Linux_x86_64_static"
QualiMap_Path="/home/achialastri/qualimap_v2.2.1/"
GTF_FILE_Path="/home/achialastri/Genomes/mm10_Anno/Mus_musculus.GRCm38.99_Formated.gtf"
FASTQC_PATH="/home/achialastri/FastQC/"
scTE_idx="/home/achialastri/Genomes/mm10_scTE/mm10.exclusive.idx"

#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2
LOG=${FASTQ_R1%??????}
LOG2=${FASTQ_R2%??????}

TRIMR1="_val_1.fq.gz"
TRIMR2="_val_2.fq.gz"
TRIMED_R1=$OUT_NAME$TRIMR1
TRIMED_R2=$OUT_NAME$TRIMR2

BAM_ENDING="Aligned.sortedByCoord.out.bam"
BAM_NAME=$OUT_NAME$BAM_ENDING
GENECOUNT_ENDING="ReadsPerGene.out.tab"
GENECOUNT_NAME=$OUT_NAME$GENECOUNT_ENDING


#Fastqc for inital check (Uncomment optionally to perform FastQC before triming)
#mkdir $START_DIR/FastQC/
#mkdir $START_DIR/FastQC/PreTrim/
#$FASTQC_PATH/fastqc $START_DIR/*_1.fq.gz  $START_DIR/*_2.fq.gz -o $START_DIR/FastQC/PreTrim/


#Copy fq files into my notation and so we do not touch the raw files at any time.
cp $START_DIR/*_1.fq.gz $START_DIR/$FASTQ_R1.all.gz
cp $START_DIR/*_2.fq.gz $START_DIR/$FASTQ_R2.all.gz

#Use Trim Galore default settings
perl /home/achialastri/BisulfiteTools/TrimGalore-0.6.5/trim_galore --path_to_cutadapt /home/achialastri/miniconda3/bin/cutadapt --length 20 --paired --basename $OUT_NAME -o $START_DIR $START_DIR/$FASTQ_R1.all.gz $START_DIR/$FASTQ_R2.all.gz

#Fastqc on Trimed Data (Uncomment optionally to perform FastQC after triming)
#mkdir $START_DIR/FastQC/PostTrim/
#$FASTQC_PATH/fastqc $START_DIR/$TRIMED_R1 $START_DIR/$TRIMED_R2 -o $START_DIR/FastQC/PostTrim/


#Align to Genome using STAR, the genome has already been indexed using the corresponding refseq .gft file so it knows where exons and introns are in the genome. Warning STAR can use alot of memory ~32GB
#Step 1 aligns just for Genes, Step 2 aligns for TEs resulting in 2 separate sparce transcript matrixes 
$MAPPER_PATH/STAR --genomeDir $ANNO_GENOME_DIR --readFilesIn $START_DIR/$TRIMED_R2 $START_DIR/$TRIMED_R1 --readFilesCommand zcat --runThreadN 6 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $START_DIR/alignments_Genes/$OUT_NAME --soloType CB_UMI_Simple --soloCBwhitelist $CEL_BARCODES --soloUMIdedup Exact --outSAMattributes NH HI AS nM CR CY UR UY --soloBarcodeReadLength 0 --soloCBstart 1   --soloCBlen 8   --soloUMIstart 9   --soloUMIlen 4 --soloFeatures Gene Velocyto  --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outMultimapperOrder Random --runRNGseed 777 --outSAMmultNmax 1

$MAPPER_PATH/STAR --genomeDir $ANNO_GENOME_DIR_TE --readFilesIn $START_DIR/$TRIMED_R2 $START_DIR/$TRIMED_R1 --readFilesCommand zcat --runThreadN 6 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $START_DIR/alignments_TE/$OUT_NAME --soloType CB_UMI_Simple --soloCBwhitelist $CEL_BARCODES --soloUMIdedup Exact --outSAMattributes NH HI AS nM CR CY UR UY --soloBarcodeReadLength 0 --soloCBstart 1   --soloCBlen 8   --soloUMIstart 9   --soloUMIlen 4 --soloFeatures Gene --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outMultimapperOrder Random --runRNGseed 777 --outSAMmultNmax 1


#Run Qualimap. This will create a .html report on the aligment (Uncomment to run Qualimap)
#$QualiMap_Path/qualimap rnaseq -pe --java-mem-size=2000M -bam $START_DIR/alignments/$BAM_NAME -gtf $GTF_FILE_Path -outdir $START_DIR/alignments/QC -p non-strand-specific


# Remove unimportant files
rm $START_DIR/$FASTQ_R1.all.gz
rm $START_DIR/$FASTQ_R2.all.gz
rm $START_DIR/$TRIMED_R1
rm $START_DIR/$TRIMED_R2

