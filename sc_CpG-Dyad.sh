
#!/bin/bash
                             
#PBS -l nodes=1:ppn=6 
                           
#PBS -l walltime=48:00:00

#User Input    Make OUT_NAME the same as the fastq files before _L00#_R#_001.fastq
START_DIR="/home/achialastri/scSerum2i_CpG/Si7L1/Si7L1_CpGDiad"




#Standard Usage, no input required
#BARCODES="/home/achialastri/perlscripts/CpG_Diad/5mCpGDiadBulkAdapters.txt"
BARCODES="/home/achialastri/perlscripts/CpG_Diad/5mCpGDiadSingleCellBarcodes.txt"
#GENOME="/home/achialastri/Genomes/hg19_Zymo_LambdaPhage/hg19_Zymo_LambdaPhage.fa"
#GENOME_BISULFITE="/home/achialastri/Genomes/hg19_Zymo_LambdaPhage"
GENOME="/home/achialastri/Genomes/mm10/mm10.fa"
GENOME_BISULFITE="/home/achialastri/Genomes/mm10"
PERL_DIR="/home/achialastri/perlscripts/CpG_Diad"
HOME_DIR="/home/achialastri/"
BisMark_Dir="/home/achialastri/Bismark-master"
SAMTOOLS_Dir="/home/cwangsanuwat/src/samtools/"
BOWTIE_Dir="/home/achialastri/bowtie2-2.3.5-linux-x86_64"

#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2

intermediateR1=${FASTQ_R1%??????}
intermediateR2=${FASTQ_R2%??????}
OUT_NAME_R1=$intermediateR1-CpGDiad
OUT_NAME2_R1=$intermediateR1-CpGDiad_bismark_bt2
OUT_NAME3_R1=$intermediateR1-CpGDiad_bismark_bt2_SE_report
OUT_NAME_R2=$intermediateR2-MSPJI

UMILength=4


#Trim Files to desired length (86 is historical from typically having 76 mapping bases in scMspJI-seq)
perl $PERL_DIR/MakeFastqShorter.pl $START_DIR/*_1.fq $START_DIR/$FASTQ_R1 86

#Deduplicate the fastq file
/home/achialastri/BisulfiteTools/bbmap/clumpify.sh in=$START_DIR/$FASTQ_R1 out=$START_DIR/$intermediateR1-dedup.fastq dedupe subs=0

rm $START_DIR/$FASTQ_R1

#Extract only read with a correct barcode and add it to the end of the read name for later useage
perl $PERL_DIR/scCpGDiad_Trim_and_ExtractBarcodes_se-HammingCorrection.pl $BARCODES $START_DIR/$intermediateR1-dedup.fastq $UMILength

#Use Trim Galore
perl /home/achialastri/BisulfiteTools/TrimGalore-0.6.5/trim_galore --path_to_cutadapt /home/achialastri/miniconda3/bin/cutadapt --length 20 -o $START_DIR $START_DIR/$intermediateR1-dedup-CpGDiad.fastq

mv $START_DIR/$intermediateR1-dedup-CpGDiad_trimmed.fq $START_DIR/$OUT_NAME_R1.fastq


#Map to bisulfite converted genome using bismark
$BisMark_Dir/bismark --genome $GENOME_BISULFITE --path_to_bowtie2 $BOWTIE_Dir --multicore 2 --samtools_path $SAMTOOLS_Dir -o $START_DIR $START_DIR/$OUT_NAME_R1.fastq

#create sam file from bismark output bam file
$SAMTOOLS_Dir/samtools view -h  $START_DIR/$OUT_NAME2_R1.bam > $START_DIR/$OUT_NAME2_R1.sam

#deduplicate sam file using bismark
$BisMark_Dir/deduplicate_bismark --barcode --samtools_path $SAMTOOLS_Dir $START_DIR/$OUT_NAME2_R1.sam
OUT_NAME3_R1=$OUT_NAME2_R1.deduplicated
mv  $HOME_DIR/$OUT_NAME3_R1.bam $START_DIR
mv  $HOME_DIR/$OUT_NAME2_R1.deduplication_report.txt $START_DIR
$SAMTOOLS_Dir/samtools view -h  $START_DIR/$OUT_NAME3_R1.bam > $START_DIR/$OUT_NAME3_R1.sam

#evaluate sam file for correct C position (MSPJI) and evaluate corresponding dyad status
perl $PERL_DIR/process_scCpGDiadwithQC.pl $GENOME $START_DIR/$OUT_NAME3_R1.sam $BARCODES $UMILength

#Extract Overall Maintenance per cell
perl $PERL_DIR/Pullout-scMaintainedCGandCHG.pl $START_DIR/$OUT_NAME3_R1.faba

#Take Faba Format into a Maintained and Unmaintained CpG txt file with corrected locations
perl $PERL_DIR/FabaToCGtxt.pl $START_DIR/$OUT_NAME3_R1.faba

#Get Total Faba Read Counts per Cell as an easy QC to look at
perl $PERL_DIR/Reads_Per_Cell_Count.pl $START_DIR/$OUT_NAME3_R1.faba

#Remove Unneeded Items
rm $START_DIR/$OUT_NAME_R1.fastq
rm $START_DIR/$intermediateR1-dedup-CpGDiad.fastq
rm $START_DIR/$intermediateR1-dedup.fastq
rm $START_DIR/$OUT_NAME2_R1.sam
rm $START_DIR/$OUT_NAME2_R1.bam
rm $START_DIR/$OUT_NAME3_R1.sam
rm $START_DIR/*.raba  $START_DIR/*.CGraba  $START_DIR/*.CHGraba #Already heavily Deduplicated so raba is uless except maybe for spike ins 

#Move reports to different area
#Move all reports to a reports folder
mkdir $START_DIR/Diad_Reports
mv $START_DIR/*.M-bias.txt $START_DIR/Diad_Reports
mv $START_DIR/*_report.txt $START_DIR/Diad_Reports
mv $START_DIR/*.QC $START_DIR/Diad_Reports
mv $START_DIR/*.QC2 $START_DIR/Diad_Reports