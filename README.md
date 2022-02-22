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

