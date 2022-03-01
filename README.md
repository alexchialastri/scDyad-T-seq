# scDyad-T-seq
Scripts for investigation of Dyad-seq methodologies inculding scDyad&amp;T-seq

### See the accompanying manuscript for more details: 




## To perform the basic data analysis of Dyad-seq (any version) you must run two separate shell scripts.
* sc_BS_CpG_Diad_SE_BothEnds.sh and sc_CpG-Dyad.sh are for Dyad-seq verions using MspJI (M-M, M-H and scDyad-seq).
* sc_BS_5hmCpG_Diad_SE_BothEnds.sh and sc_5hmCpG-Dyad.sh are for Dyad-seq version using AbaSI (H-M and H-H).
* sc_mRNA_NOVOGENE_FQ_Analysis_Using_STAR_mm10.sh is used to quanitify single-cell RNA-seq detection in scDyad&T-seq.

* These scripts also require a Linux operating system, Perl, Clumpify from BBTools, Bismark mapper, Bowtie 2, Trim Galore, Samtools, and STAR (STARsolo)

#### 
#### Detecting absolute 5mC or 5hmC levels
sc_BS_CpG_Diad_SE_BothEnds.sh & sc_BS_5hmCpG_Diad_SE_BothEnds.sh perform the same steps using a different barcode set.
1. Fastq files are demultiplexed using Clumpify
2. The cell barcode and UMI are appended to the Read 1 and Read 2 Fastq read names
3. Trim Galore is used to remove any bias from enzymatic cutting (by AbaSI or MspJI) and random priming
4. Read 1 and Read 2 are mapped independently to the genome using Bismark mapper. 
5. Bismark mapper is used to deduplicate the mapped bam file based on mapping location and CellBarcode/UMI. The bam file is converted to a sam file using samtools
6. Bismark mapper is used to extract called 5mC and unmethylated cytosine locations (similar result if 5hmC conversion if performed but methylated calls refer to 5hmC)
7. Methylation Concordance is scored for reference
8. Extracted 5mC and C calls are reformated, then Linux command line functions are used to concatinate calls from Read 1 and read 2. Then these calls are further reformated and deduplicated based on the UMI and cell barcode (for instance when read 1 and read 2 are overlapping). T
9. The 5mC and C counts per cell for each strand and for both strands combined are recorded for reference.
10. Step 8 and 9 are performed for CpG and non-CpG methylation
11. The top and bottom strand files are concatinated into all 5mC and C sites detected. The results for each cell are then binned by a user defined bin size.
12. Files are cleaned up and removed. Other files are moved into different folders based on similarities.
13. The resulting "*-CpG_All_BothStrands_Unique.txt" files list detected methylated or unmethylated CpG sites. Column 1 is the cell number (in bulk this is the sample replicate indicator), column 2 is the chromosome (with chr23 indicating the X chromosome, chr24 indicating the Y chromosome, and chr25 indicating a mitocondrial based detection), column 3 is the genomic loci, column 4 is the the strand of the detected mark, column 5 is the UMI detected and column 6 indicates if the cytosine was methylated (Z) or unmethylated (z).

#### 
#### Detecting 5mCpG Maintenance (sc_CpG-Dyad.sh) or 5hmCpG Maintenance (sc_5hmCpG-Dyad.sh)
sc_CpG-Dyad.sh & sc_5hmCpG-Dyad.sh perform similar steps using different barcode sets and different scripts to account for differences between the cutsites of MspJI and AbaSI. This detection only requires read 1
1. The read 1 Fastq file is trimmed to 86 bases.
2. The trimmed Fastq file is demultiplexed using Clumpify
3. The cell barcode and UMI are appended to the Read 1 Fastq read names
4. Trim Galore is used to remove Illumina adapters and improve mapping
5. Read 1 is mapped to the genome using Bismark mapper.
6. Bismark mapper is used to deduplicate the mapped bam file based on mapping location and CellBarcode/UMI. The bam file is converted to a sam file using samtools.
7. The inferred detection of 5mC (MspJI based methods) or 5hmC (AbaSI based methods) is performed, using the known cutsite of the enzyme used. Probing the direct 5mC calls from the sequencing allows the calling of fully or hemi-methylated sites.
8. A summary of detection in each cell is produced. The number of fully or hemi-methylated sites found in a CpG for (AbaSI and MspJI methods) and for CpHpG (for MspJI methods only). The detection of indirect calls of 5mC (MspJI based methods) or 5hmC (AbaSI based methods) are recorded. The final fully or hemi-methylated detection files are created.
9. The resulting 5mCpG Dyad maintenance detection files [noted in the naming of the file as *_CGMaintained.txt (for loci detected as maintained {e.g. fully methylated dyads for M-M-Dyad-seq}) and *_CGNonMaintained.txt (for loci detected as not maintained  {e.g. hemi-methylated dyads for M-M-Dyad-seq})] list detected methylated or unmethylated CpG maintenance sites. Column 1 is the cell number (in bulk this is the sample replicate indicator), column 2 is the chromosome (with chr23 indicating the X chromosome, chr24 indicating the Y chromosome, and chr25 indicating a mitocondrial based detection), column 3 is the genomic loci, column 4 is the UMI detected and column 4 is the strand of the detected mark.

#### 
#### Detecting RNA transcripts in scDyad&T-seq
RNA transcripts are mapped to the genome twice. First to look for genes and then to look for expression of transposable element
1. Trim Galore is used to remove Illumina adapters and improve mapping
2. The resulting reads are mapped using STARsolo to the genome that has been annotated using the gene annotation from Ensembl.
3. The resulting reads are mapped using STARsolo to the genome that has been annotated using the transposable elements annotation file described in TEtranscripts (Jin, Bioinformatics 2015)
4. Resulting sparce matrixes are the same as produced in the 10x Cell Ranger, and can be analyzed using Seurat.



#### 
## How to use each script
scCpGDiad_Trim_and_ExtractBarcodes_se-HammingCorrection_TransferToR2.pl - Used to identify reads with the proper PCR sequence "GGTGTAGTGGGTTTGG" and proper barcodes. Also concatenates the UMI and cell barcode to the read name for read1 and 2
* Argument 0 = Barcode list
* Argument 1 = Read 1 Input Fastq file
* Argument 2 = UMI length (Typically 4 for most experiments performed in manuscript)
* Argument 3 = Read 2 Input Fastq file

MethylationConcordanceScoring.pl -
*Argument 0 = Mapped deduplicated sam file from Bismark
*Argument 1 = Barcode list
*Argument 2 = UMI length (typically = 4)
*Argument 3 = SpikeInFlag (typically = 2, genome only)
*Argument 4 = Number of CpGs needed on a single read (typically = 5)
*Argument 5 = Percent threshold to call read as fully methylated (typically = 90)
*Argument 6 = Percent threshold to call read as fully unmethylated (typically = 10)

Pullout-scMethylationFromBismarkMethylationExtractor_wReadNames.pl - Reformats the Bismark Methylation Extractor format to demultiplex by cell barcode and UMI.
Argument 0 = Bismark Methylation Extractor text file. There are 8 options. [CpG_OB_*txt, CpG_CTOB_*txt, CpG_OT_*txt, CpG_CTOT_*txt, Non_CpG_OB_*txt, Non_CpG_CTOB_*txt, Non_CpG_OT_*txt, Non_CpG_CTOT_*txt]. Run this script on each of the Bismark Methylation Extractor output files.
Argument 1 = Barcode list
Argument 2 = UMI length (typically = 4)

Reads_Per_Cell_Count_BismarkZorz.pl - summarizes the count detection per cell as two output files in order by cell number. Number of methylated 5mCpGs found or number of unmethylated CpGs found.
Argument 0 = Text file containing 6 columns space delimited (cell chr loci strand UMI MethCall{Z,z})

Reads_Per_Cell_Count_Bismark_NonCGMethCount.pl - summarizes the count detection per cell as two output files in order by cell number. Number of methylated non CpG 5mCs found or number of unmethylated non CpG cytosine found.
Argument 0 = Text file containing 6 columns space delimited (cell chr loci strand UMI MethCall{XH,xh})

Meth_Percent_PerCell_BS.pl - Summarizes all strands together. Inputs are *MethylatedCountsPerCell.txt files. Can be used on CpG or nonCpG files. Results in a file with two columns describing the total detection of sites and percent methylated of sites per cell, in order of cell number.
Argument 0 = *_OB_MethylatedCountsPerCell.txt
Argument 1 = *OB_UnmethylatedCountsPerCell.txt
Argument 2 = *OT_MethylatedCountsPerCell.txt
Argument 3 = *OT_UnmethylatedCountsPerCell.txt
Argument 4 = Desired output file name.

Binning_MouseGenomeForStrandBias_BS.pl - Takes 5mC calls for all cells and demultiplexes the cells as well as counts detection in each bin of the genome of a user defined size. Results in 3 files, the number of methylated calls in the bin by cell, the number of unmethylated calls in the bin by cell, and the resulting division of these two tables aka the percent methylation of each bin by cell.
Argument 0 = Text file containing 6 columns space delimited (cell chr loci strand UMI MethCall{Z,z}) typically all sites detected combined on all strands after deduplication.
Argument 1 = Bin size

MakeFastqShorter.pl - Used to hard clip Fastq files
Argument 0 = Input Fastq file name
Argument 1 = Desired Output Fastq file name
Argument 2 = Desired Fastq file trimmed length

scCpGDiad_Trim_and_ExtractBarcodes_se-HammingCorrection.pl - Used to identify reads with the proper PCR sequence "GGTGTAGTGGGTTTGG" and proper barcodes.
Argument 0 = Barcode list
Argument 1 = Input Fastq file
Argument 2 = UMI length (Typically 4 for most experiments performed in manuscript)

process_scCpGDiadwithQC.pl - Searches the genome location of mapped reads to identify fully or hemi-methylated CpG or CpH sites. Also detects all 5mC calls inferred by the cutting action of MspJI.
Argument 0 = Genome used for mapping (non bisulfite converted)
Argument 1 = Mapped deduplicated sam file from Bismark
Argument 2 = Cell barcodes used for M-M-Dyad-seq, M-H-Dyad-seq or scDyad-seq.
Argument 3 = UMI length used (typically = 4)

Pullout-scMaintainedCGandCHG.pl - Summarizes the detection per cell. Resulting .fabaMaintained text file gives the cell number, the number of fully methylated detections, the number of hemimethylated detections, and the percent fully methylated for CpG sites followed by CpHpG sites in the next 3 columns.
Argument 0 = .faba file generated by the process_scCpGDiadwithQC.pl script.

FabaToCGtxt.pl - Takes .faba file and creates a CpG maintained (fully methylated) and non-maintained files (hemi-methylated) text files. The exact loci of the CpG site is modified slightly from the faba file to correct for which strand it was detected on. These resulting files are used in downstream analysis and are submitted on GEO.
Argument 0 = .faba file generated by the process_scCpGDiadwithQC.pl script.

Reads_Per_Cell_Count.pl - summarizes the number of 5mC detections inferred MspJI. This value is always greater than the number of CpG sites detected with maintenance information, and accounts for detections in non CpG context as well.
Argument 0 = .faba file generated by the process_scCpGDiadwithQC.pl script.

process_sc5hmC_CpGDiad.pl - Searches the genome location of mapped reads to identify fully or hemi-hydroxymethylated CpG sites. Also detects all 5hmC calls inferred by the cutting action of AbaSI.
Argument 0 = Genome used for mapping (non bisulfite converted)
Argument 1 = Mapped deduplicated sam file from Bismark
Argument 2 = Cell barcodes used for H-M-Dyad-seq, H-H-Dyad-seq.
Argument 3 = UMI length used (typically = 4)

Pullout-scMaintainedCG_ForAbaSI.pl - Summarizes the detection per cell. Resulting .fabaMaintained text file gives the cell number, the number of fully hydroxymethylated detections, the number of hemi-hydroxymethylated detections, and the percent fully hydroxymethylated for CpG sites.
Argument 0 = .faba file generated by the process_sc5hmC_CpGDiad.pl script.

AbaSIFabaToCGtxt.pl - Takes .faba file and creates a CpG hydroxymaintained (fully hydroxymethylated) and non-maintained files (hemi-hydroxymethylated) text files. These resulting files are used in downstream analysis and are submitted on GEO.
Argument 0 = .faba file generated by the process_sc5hmC_CpGDiad.pl script.
