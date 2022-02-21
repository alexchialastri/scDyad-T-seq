
#!/bin/bash
                             
#PBS -l nodes=1:ppn=6 
                           
#PBS -l walltime=48:00:00

#User Input    Make OUT_NAME the same as the fastq files before _L00#_R#_001.fastq
START_DIR="/home/achialastri/Expanding_CpG_Diad_BulkSL_2i_M_H_Rep1/S_M_H_rep1/BS_Test"




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
OUT_NAME_R1=$intermediateR1-BS_of_CpGDiad
OUT_NAME2=$OUT_NAME-dedup-BS_of_CpGDiad_bismark_bt2
OUT_NAME2_Report=$OUT_NAME-dedup-BS_of_CpGDiad_bismark_bt2_report
OUT_NAME_R2=$intermediateR2-MSPJI

#Other Parameters
UMILength=4
EnzymeClipLength=20
RandomPrimerLength=9
BCLength=10
HardClipThreePrimeR2Val=$(($UMILength + $BCLength + $EnzymeClipLength))


#Cat Fastq Files if needed first before running.

#perl /home/achialastri/perlscripts/MakeFastqShorter.pl $START_DIR/*_1.fq $START_DIR/$FASTQ_R1 86
cp $START_DIR/*_1.fq $START_DIR/$FASTQ_R1
cp $START_DIR/*_2.fq $START_DIR/$FASTQ_R2

#/home/achialastri/BisulfiteTools/bbmap/clumpify.sh in=$START_DIR/$FASTQ_R1 out=$START_DIR/$intermediateR1-dedup-R1only.fastq dedupe subs=0
#Deduplicate the fastq file
/home/achialastri/BisulfiteTools/bbmap/clumpify.sh in=$START_DIR/$FASTQ_R1 in2=$START_DIR/$FASTQ_R2 out=$START_DIR/$intermediateR1-dedup.fastq out2=$START_DIR/$intermediateR2-dedup.fastq dedupe subs=0

rm $START_DIR/$FASTQ_R1
rm $START_DIR/$FASTQ_R2

#Extract only read with a correct barcode and add it to the end of the read name for later useage
perl $PERL_DIR/scCpGDiad_Trim_and_ExtractBarcodes_se-HammingCorrection_TransferToR2.pl $BARCODES $START_DIR/$intermediateR1-dedup.fastq $UMILength $START_DIR/$intermediateR2-dedup.fastq

rm $START_DIR/$intermediateR1-dedup.fastq 
rm $START_DIR/$intermediateR2-dedup.fastq

#Use Trim Galore
perl /home/achialastri/BisulfiteTools/TrimGalore-0.6.5/trim_galore --path_to_cutadapt /home/achialastri/miniconda3/bin/cutadapt --length 20 --paired --clip_R1 $EnzymeClipLength --clip_R2 $RandomPrimerLength --three_prime_clip_R2 $HardClipThreePrimeR2Val --three_prime_clip_R1 $RandomPrimerLength --adapter2 CCACATCACCCAAACC -o $START_DIR $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad.fastq $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad.fastq
   #output contains same names as in put but remove the.fastq and end with *_val_1.fq for R1 and *_val_2.fq


rm $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad.fastq 
rm $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad.fastq

#Map to bisulfite converted genome using bismark
$BisMark_Dir/bismark --genome $GENOME_BISULFITE --path_to_bowtie2 $BOWTIE_Dir --multicore 2 --samtools_path $SAMTOOLS_Dir -o $START_DIR $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1.fq
$BisMark_Dir/bismark --genome $GENOME_BISULFITE --path_to_bowtie2 $BOWTIE_Dir --multicore 2 --samtools_path $SAMTOOLS_Dir --pbat -o $START_DIR $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2.fq

rm $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1.fq
rm $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2.fq


#deduplicate sam file using bismark
 #R1
$BisMark_Dir/deduplicate_bismark --barcode --samtools_path $SAMTOOLS_Dir $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1_bismark_bt2.bam
mv $HOME_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1_bismark_bt2.deduplicated.bam $START_DIR
mv $HOME_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1_bismark_bt2.deduplication_report.txt $START_DIR
OUT_NAME3_R1=$OUT_NAME2.deduplicated_R1
OUT_NAME3_Report_R1=$OUT_NAME2.deduplication_report_R1
mv $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1_bismark_bt2.deduplicated.bam $START_DIR/$OUT_NAME3_R1.bam
mv $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1_bismark_bt2.deduplication_report.txt $START_DIR/$OUT_NAME3_Report_R1.txt
$SAMTOOLS_Dir/samtools view -h  $START_DIR/$OUT_NAME3_R1.bam > $START_DIR/$OUT_NAME3_R1.sam

 #R2
$BisMark_Dir/deduplicate_bismark --barcode --samtools_path $SAMTOOLS_Dir $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2_bismark_bt2.bam
mv $HOME_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2_bismark_bt2.deduplicated.bam $START_DIR
mv $HOME_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2_bismark_bt2.deduplication_report.txt $START_DIR
OUT_NAME3_R2=$OUT_NAME2.deduplicated_R2
OUT_NAME3_Report_R2=$OUT_NAME2.deduplication_report_R2
mv $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2_bismark_bt2.deduplicated.bam $START_DIR/$OUT_NAME3_R2.bam
mv $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2_bismark_bt2.deduplication_report.txt $START_DIR/$OUT_NAME3_Report_R2.txt
$SAMTOOLS_Dir/samtools view -h  $START_DIR/$OUT_NAME3_R2.bam > $START_DIR/$OUT_NAME3_R2.sam


#extract traditional bisulfite data from the reads
$BisMark_Dir/bismark_methylation_extractor --samtools_path $SAMTOOLS_Dir --merge_non_CpG --report --multicore 2 -s -o $START_DIR $START_DIR/$OUT_NAME3_R1.sam $START_DIR/$OUT_NAME3_R2.sam

#extract Methylation concordance score (are large strings of 5mC fully methylated or not?) ARVGs = barcodes, UMI Length, SpikeINFLAG (2=just genome), Num of CGs needed in the read, Threshold to call fully meth, threshold to call fully unmeth
perl $PERL_DIR/MethylationConcordanceScoring.pl $START_DIR/$OUT_NAME3_R1.sam $BARCODES $UMILength 2 5 90 10
perl $PERL_DIR/MethylationConcordanceScoring.pl $START_DIR/$OUT_NAME3_R2.sam $BARCODES $UMILength 2 5 90 10

#Turn Extracted bisulfite reads into our faba format
   #CpG
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/CpG_OB_*txt $BARCODES $UMILength
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/CpG_CTOB_*txt $BARCODES $UMILength
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/CpG_OT_*txt $BARCODES $UMILength
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/CpG_CTOT_*txt $BARCODES $UMILength
   #Non_CpG
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/Non_CpG_OB_*txt $BARCODES $UMILength
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/Non_CpG_CTOB_*txt $BARCODES $UMILength
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/Non_CpG_OT_*txt $BARCODES $UMILength
perl $PERL_DIR/Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl $START_DIR/Non_CpG_CTOT_*txt $BARCODES $UMILength

#Merge 5mCpG from R1 and R2, 
cat $START_DIR/CpG_OB_*txtFormat $START_DIR/CpG_CTOB_*txtFormat > $START_DIR/CpG_All_OB.txtFormat
sort -u $START_DIR/CpG_All_OB.txtFormat >  $START_DIR/CpG_All_OB.txtFormatu
awk '{ print $1, $2, $3, $4, $5, $6 }' $START_DIR/CpG_All_OB.txtFormatu > $START_DIR/$OUT_NAME2-CpG_All_OB.txt
perl $PERL_DIR/Reads_Per_Cell_Count_BismarkZorz.pl $START_DIR/$OUT_NAME2-CpG_All_OB.txt 
rm $START_DIR/CpG_OB_*txtFormat $START_DIR/CpG_CTOB_*txtFormat $START_DIR/CpG_All_OB.txtFormat $START_DIR/CpG_All_OB.txtFormatu
rm $START_DIR/CpG_OB_*txt $START_DIR/CpG_CTOB_*txt

cat $START_DIR/CpG_OT_*txtFormat $START_DIR/CpG_CTOT_*txtFormat > $START_DIR/CpG_All_OT.txtFormat
sort -u $START_DIR/CpG_All_OT.txtFormat >  $START_DIR/CpG_All_OT.txtFormatu
awk '{ print $1, $2, $3, $4, $5, $6 }' $START_DIR/CpG_All_OT.txtFormatu > $START_DIR/$OUT_NAME2-CpG_All_OT.txt
perl $PERL_DIR/Reads_Per_Cell_Count_BismarkZorz.pl $START_DIR/$OUT_NAME2-CpG_All_OT.txt 
rm $START_DIR/CpG_OT_*txtFormat $START_DIR/CpG_CTOT_*txtFormat $START_DIR/CpG_All_OT.txtFormat $START_DIR/CpG_All_OT.txtFormatu
rm $START_DIR/CpG_OT_*txt $START_DIR/CpG_CTOT_*txt

#Do samething but for Non_CpG
cat $START_DIR/Non_CpG_OB_*txtFormat $START_DIR/Non_CpG_CTOB_*txtFormat > $START_DIR/Non_CpG_All_OB.txtFormat
sort -u $START_DIR/Non_CpG_All_OB.txtFormat >  $START_DIR/Non_CpG_All_OB.txtFormatu
awk '{ print $1, $2, $3, $4, $5, $6 }' $START_DIR/Non_CpG_All_OB.txtFormatu > $START_DIR/$OUT_NAME2-Non_CpG_All_OB.txt
perl $PERL_DIR/Reads_Per_Cell_Count_Bismark_NonCGMethCount.pl $START_DIR/$OUT_NAME2-Non_CpG_All_OB.txt 
rm $START_DIR/Non_CpG_OB_*txtFormat $START_DIR/Non_CpG_CTOB_*txtFormat $START_DIR/Non_CpG_All_OB.txtFormat $START_DIR/Non_CpG_All_OB.txtFormatu
rm $START_DIR/Non_CpG_OB_*txt $START_DIR/Non_CpG_CTOB_*txt

cat $START_DIR/Non_CpG_OT_*txtFormat $START_DIR/Non_CpG_CTOT_*txtFormat > $START_DIR/Non_CpG_All_OT.txtFormat
sort -u $START_DIR/Non_CpG_All_OT.txtFormat >  $START_DIR/Non_CpG_All_OT.txtFormatu
awk '{ print $1, $2, $3, $4, $5, $6 }' $START_DIR/Non_CpG_All_OT.txtFormatu > $START_DIR/$OUT_NAME2-Non_CpG_All_OT.txt
perl $PERL_DIR/Reads_Per_Cell_Count_Bismark_NonCGMethCount.pl $START_DIR/$OUT_NAME2-Non_CpG_All_OT.txt 
rm $START_DIR/Non_CpG_OT_*txtFormat $START_DIR/Non_CpG_CTOT_*txtFormat $START_DIR/Non_CpG_All_OT.txtFormat $START_DIR/Non_CpG_All_OT.txtFormatu
rm $START_DIR/Non_CpG_OT_*txt $START_DIR/Non_CpG_CTOT_*txt

#Create PerCell 5mC Percent Files from BS data (order is OB_Meth, OB_Unmeth, OT_Meth, OT_Unmeth, OutputName)
perl $PERL_DIR/Meth_Percent_PerCell_BS.pl $START_DIR/$OUT_NAME2-CpG_All_OB_MethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_OB_UnmethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_OT_MethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_OT_UnmethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_MethPercentPerCell.txt

perl $PERL_DIR/Meth_Percent_PerCell_BS.pl $START_DIR/$OUT_NAME2-Non_CpG_All_OB_MethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-Non_CpG_All_OB_UnmethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-Non_CpG_All_OT_MethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-Non_CpG_All_OT_UnmethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-Non_CpG_All_MethPercentPerCell.txt

#Combine OB OT Files and make them unique
sort -u $START_DIR/$OUT_NAME2-CpG_All_OB.txt > $START_DIR/$OUT_NAME2-CpG_All_OB_Unique.txt
sort -u $START_DIR/$OUT_NAME2-CpG_All_OT.txt > $START_DIR/$OUT_NAME2-CpG_All_OT_Unique.txt
perl $PERL_DIR/Binning_MouseGenomeForStrandBias_BS.pl $START_DIR/$OUT_NAME2-CpG_All_OB_Unique.txt 10000000
perl $PERL_DIR/Binning_MouseGenomeForStrandBias_BS.pl $START_DIR/$OUT_NAME2-CpG_All_OB_Unique.txt 9999999999
perl $PERL_DIR/Binning_MouseGenomeForStrandBias_BS.pl $START_DIR/$OUT_NAME2-CpG_All_OT_Unique.txt 10000000
perl $PERL_DIR/Binning_MouseGenomeForStrandBias_BS.pl $START_DIR/$OUT_NAME2-CpG_All_OT_Unique.txt 9999999999

perl $PERL_DIR/Reads_Per_Cell_Count_BismarkZorz.pl $START_DIR/$OUT_NAME2-CpG_All_OB_Unique.txt
perl $PERL_DIR/Reads_Per_Cell_Count_BismarkZorz.pl $START_DIR/$OUT_NAME2-CpG_All_OT_Unique.txt

perl $PERL_DIR/Meth_Percent_PerCell_BS.pl $START_DIR/$OUT_NAME2-CpG_All_OB_Unique_MethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_OB_Unique_UnmethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_OT_Unique_MethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_OT_Unique_UnmethylatedCountsPerCell.txt $START_DIR/$OUT_NAME2-CpG_All_Unique_MethPercentPerCell.txt

#Combine the OB and OT unique files
cat $START_DIR/$OUT_NAME2-CpG_All_OB_Unique.txt $START_DIR/$OUT_NAME2-CpG_All_OT_Unique.txt  > $START_DIR/$OUT_NAME2-CpG_All_BothStrands_Unique.txt
perl $PERL_DIR/Binning_MouseGenomeForStrandBias_BS.pl $START_DIR/$OUT_NAME2-CpG_All_BothStrands_Unique.txt 10000000
perl $PERL_DIR/Binning_MouseGenomeForStrandBias_BS.pl $START_DIR/$OUT_NAME2-CpG_All_BothStrands_Unique.txt 9999999999

#Remove more unneeded files
rm $START_DIR/$OUT_NAME3_R1.sam
rm $START_DIR/$OUT_NAME3_R2.sam
rm $START_DIR/$intermediateR1-dedup-BS_of_CpGDiad_val_1_bismark_bt2.bam
rm $START_DIR/$intermediateR2-dedup-BS_of_CpGDiad_val_2_bismark_bt2.bam

#Move all reports to a reports folder
mkdir $START_DIR/Reports
mv $START_DIR/*.M-bias.txt $START_DIR/Reports
mv $START_DIR/*_report.txt $START_DIR/Reports
mv $START_DIR/*_report_R1.txt $START_DIR/Reports
mv $START_DIR/*_report_R2.txt $START_DIR/Reports

#Move Per Cell Counts to a PerCell Folder
mkdir $START_DIR/PerCell
mv $START_DIR/*Unique_*ethylatedCountsPerCell.txt $START_DIR/PerCell
mv $START_DIR/*CpG_All_Unique_MethPercentPerCell.txt $START_DIR/PerCell

#Move non Unique to a NonUnique Folder
mkdir $START_DIR/NonUnique
mv $START_DIR/$OUT_NAME2-CpG_All_OB.txt $START_DIR/NonUnique
mv $START_DIR/$OUT_NAME2-CpG_All_OT.txt $START_DIR/NonUnique
mv $START_DIR/$OUT_NAME2-Non_CpG_All_OB.txt $START_DIR/NonUnique
mv $START_DIR/$OUT_NAME2-Non_CpG_All_OT.txt $START_DIR/NonUnique
mv $START_DIR/*ethylatedCountsPerCell.txt $START_DIR/NonUnique
mv $START_DIR/*CpG_All_MethPercentPerCell.txt $START_DIR/NonUnique

#Created BinedBS folder
mkdir $START_DIR/BinnedBS
mv $START_DIR/$OUT_NAME2-CpG_All_*_Unique*_10000000bp.txt $START_DIR/BinnedBS
mv $START_DIR/$OUT_NAME2-CpG_All_*_Unique*_Full_Chr.txt $START_DIR/BinnedBS