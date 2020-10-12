#!/usr/bin/env bash

##Mike Mariani UVM 2020 Frietze Lab
##Adapted from https://github.com/broadinstitute/Drop-seq/blob/v2.1.0/doc/Drop-seq_Alignment_Cookbook.pdf

##10/02/2020
##Let's try to piece together the DropSeq pipeline
##and try it out on the Drayman dataset first
##Note there is a mix of my notes and verbatim
##text from the dropseq cook book linked to above.

##Below command is for submission of this script
##to our SGE galaxy server, the uncommented code
##can be pasted in the terminal.
##qsub -V -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/logs /slipstream/home/mmariani/scripts/hsv1_scrna_scripts/drop_seq/drop_seq_pipeline_mm.bash

##Overview of Alignment 
##The raw reads from the sequencer must be converted into a Picard-queryname-sorted 
##BAM file for each library in the sequencer run.    Since there are many sequencers 
##and pipelines available to do this, we leave this step to the user.  For example, 
##we use either Picard ​IlluminaBasecallsToSam​ ​​(preceded by Picard ​ExtractIlluminaBarcodes​ 
##for a library with sample barcodes);​ ​​or Illumina’s ​bcl2fastq​ ​​followed by Picard ​FastqToSam​.  
##​​Once you have an unmapped, queryname-sorted BAM, you can follow this set of steps to 
##align your raw reads and create a BAM file that is suitable to produce digital gene expression 
##(DGE) results. 

cd /slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/data/wt

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526684.1_1.fastq \
F2=SRR8526684.1_2.fastq \
O=SRR8526684.1.bam \
SM=SRR8526684.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"
 
/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526682_1.fastq \
F2=SRR8526682_2.fastq \
O=SRR8526682.bam \
SM=SRR8526682 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526681.1_1.fastq \
F2=SRR8526681.1_2.fastq \
O=SRR8526681.1.bam \
SM=SRR8526681.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526683.1_1.fastq \
F2=SRR8526683.1_2.fastq \
O=SRR8526683.1.bam \
SM=SRR8526683.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

cd /slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/data/mock

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526684.1_1.fastq \
F2=SRR8526684.1_2.fastq \
O=SRR8526684.1.bam \
SM=SRR8526684.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526683.1_1.fastq \
F2=SRR8526683.1_2.fastq \
O=SRR8526683.1.bam \
SM=SRR8526683.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526682_1.fastq \
F2=SRR8526682_2.fastq \
O=SRR8526682.bam \
SM=SRR8526682 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526681.1_1.fastq \
F2=SRR8526681.1_2.fastq \
O=SRR8526681.1.bam \
SM=SRR8526681.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526679.1_1.fastq \
F2=SRR8526679.1_2.fastq \
O=SRR8526679.1.bam \
SM=SRR8526679.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526678.1_1.fastq \
F2=SRR8526678.1_2.fastq \
O=SRR8526678.1.bam \
SM=SRR8526678.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

/slipstream/home/mmariani/programs/miniconda3/bin/picard FastqToSam \
F1=SRR8526677.1_1.fastq \
F2=SRR8526677.1_2.fastq \
O=SRR8526677.1.bam \
SM=SRR8526677.1 \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

######################################## METADATA CREATION ###############################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

##MetaData Creation Programs 
##A few files and required to generate meta data: a GTF file, and a fastq file.  
##From these two files we can derive various other files needed by the the Drop-seq software.  

##cd /slipstream/home/mmariani/references/drop_seq
cd /slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1

##CreateSequenceDictionary 
##The first file needed is the sequence dictionary.  This is a list of the contigs 
##in the fastq file and their lengths. 

##Downloaded files from BROAD:
##/slipstream/home/mmariani/references/drop_seq/GSM1629193_hg19_ERCC.dict.txt
##/slipstream/home/mmariani/references/drop_seq/GSM1629193_hg19_ERCC.fasta
##/slipstream/home/mmariani/references/drop_seq/GSM1629193_hg19_ERCC.gtf
##/slipstream/home/mmariani/references/drop_seq/GSM1629193_hg19_ERCC.refFlat.txt

##Or make my own using my combined GRCH38+HSV1 reference:

##Note might need to edit gtf file to include "transcipt_name" field:
##sed 's/\r//g' genes.gtf \
##|  awk -F $'\t' '{if(index($9,"transcript_id")!=0 && $1=="NC_001806.2")\
##{split($9,a,";");changeVar=a[1];};\
##gsub("transcript_id","transcript_name",changeVar);\
##print $0"; "changeVar;}' > genes.mm.edit.for.dropseq.gtf

##Define variables related to reference creation 
##and alignment with the STAR aligner

ref_fasta="/slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1/fasta/genome.fa"
ref_gtf="/slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1/genes/genes.mm.edit.for.dropseq.gtf"
##ref_dict="/slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1/dropseq/hg38_and_hsv1.dict"
ref_dict="/slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1/fasta/hg38_and_hsv1.dict"
ref_flat="/slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1/dropseq/hg38_and_hsv1.refFlat"
ref_reduced_gtf="/slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1/dropseq/genes.reduced.gtf"
ref_star="/slipstream/home/mmariani/references/cellranger/hsv1_10x_references/mm_cellranger_combined_reference_make_with_cell_ranger/hg38_and_hsv1/star_dropseq"

##Create the sequence dictionary

/slipstream/home/mmariani/programs/miniconda3/bin/picard CreateSequenceDictionary \
REFERENCE=$ref_fasta \
OUTPUT=$ref_dict \
SPECIES=hg38_and_hsv1

##ConvertToRefFlat 
##The next file is the refFlat file, which is generated using the sequence 
##dictionary generated above. 

/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/ConvertToRefFlat \
ANNOTATIONS_FILE=$ref_gtf \
SEQUENCE_DICTIONARY=$ref_dict \
OUTPUT=$ref_flat

##ReduceGTF 
##The may be useful if you need an easy to parse version of your annotations 
##in a language like R, and is also used to generate the other metadata. 

/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/ReduceGtf \
SEQUENCE_DICTIONARY=$ref_dict \
GTF=$ref_gtf \
OUTPUT=$ref_reduced_gtf \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

####CreateIntervalsFiles 
####As a last step, we create interval files needed for various programs 
####in the Drop-seq pipeline.  This program generates a number of interval 
####files for genes, exons, consensus introns, rRNA, and mt.  
####The example below uses the human MT contig name, but if you use 
####a different organism you should set that argument appropriately. 

/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/CreateIntervalsFiles \
SEQUENCE_DICTIONARY=$ref_dict \
REDUCED_GTF=$ref_reduced_gtf \
PREFIX=hg38_and_hsv1 \
OUTPUT=/slipstream/home/mmariani/references/drop_seq/intervals_files ##\
##MT_SEQUENCE=MT 

##The above command was not working with "MT_SEQUENCE=MT" so 
##I removed it, perhaps the mitochondiral features are not defined?

##MetaData Generation Pipeline 
##We’ve provided a shell script to generate new meta data sets for single organism 
##data in the distribution.  This script is called create_Drop-seq_reference_metadata.sh, 
##and the options for the program can be accessed by running with the -h option: 

##/path/to/dropseq_tools/create_Drop-seq_reference_metadata.sh -h 

######################################################## ALIGNMENT #################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

##Alignment Pipeline Programs 
##On the Drop-seq website you will find a zipfile containing the programs described below.  
##The zipfile also contains a script Drop-seq_alignment.sh that executes the process described below.  
##Because of differences in computing environments, this script is not guaranteed to work for all users.  
##However, we hope it will serve as an example of how the various programs should be invoked. 

##TagBamWithReadSequenceExtended 
##This Drop-seq program extracts bases from the cell/molecular barcode encoding read 
##(BARCODED_READ), and creates a new BAM tag with those bases on the ​genome read​.  
##By default, we use the BAM tag XM for molecular barcodes, and XC for cell barcodes, 
##using the TAG_NAME parameter. 
##This program is run once per barcode extraction to add a tag.  On the first iteration, 
##the cell barcode is extracted from bases 1-12.  This is controlled by the BASE_RANGE option.  
##On the second iteration, the molecular barcode is extracted from bases 13-20 of the barcode read.  
##This program has an option to drop a read (DISCARD_READ), which we use after both barcodes 
##have been extracted, which makes the output BAM have unpaired reads with additional tags.  
##Additionally, this program has a BASE_QUALITY option, which is the minimum ​base quality​ of 
##all bases of the barcode being extracted.  If more than NUM_BASES_BELOW_QUALITY bases falls 
##below this quality, the read pair is discarded.  

cd /slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/output

##Select the data dir you want to draw from and the output dir to 
##send alignment files to. I comment out the ones that I am not 
##using.  You can switch the commented variables and resubmit 
##on a separate node to have both datasets crunching in parallel

##data_dir="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/data/mock"
data_dir="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/data/wt"

##output_dir="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/output/mock"
output_dir="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/output/wt"

cd $data_dir

arrayNow=($(ls | cat | grep ".bam"))

for i in "${arrayNow[@]}"
do

##mkdir $(basename $i ".bam")
##cd $(basename $i ".bam")

prefix=$output_dir"/"$(basename $i ".bam")

##Example Cell Barcode: 
/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/TagBamWithReadSequenceExtended \
INPUT=$i \
OUTPUT=$prefix".unaligned_tagged_Cell.bam" \
SUMMARY=$prefix".unaligned_tagged_Cellular.bam_summary.txt" \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1 

##Example Molecular Barcode: 
/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/TagBamWithReadSequenceExtended \
INPUT=$prefix".unaligned_tagged_Cell.bam" \
OUTPUT=$prefix".unaligned_tagged_CellMolecular.bam" \
SUMMARY=$prefix".unaligned_tagged_Molecular.bam_summary.txt" \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1  

##FilterBam: 
##This Drop-seq program is used to remove reads where the cell or molecular barcode 
##has low quality bases.  During the run of TagBamWithReadSequenceExtended, 
##an XQ tag is added to each read to represent the number of bases that have 
##quality scores below the BASE_QUALITY threshold. These reads are then 
##removed from the BAM. 
 
/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/FilterBam \
TAG_REJECT=XQ \
INPUT=$prefix".unaligned_tagged_CellMolecular.bam" \
OUTPUT=$prefix".unaligned_tagged_filtered.bam" 

##TrimStartingSequence 
##This Drop-seq program is one of two sequence cleanup programs designed to trim away 
##any extra sequence that might have snuck it’s way into the reads.  In this case, 
##we trim the SMART Adapter that can occur 5’ of the read.  In our standard run, 
##we look for at least 5 contiguous bases (NUM_BASES) of the SMART adapter (SEQUENCE) 
##at the 5’ end of the read with no errors (MISMATCHES) , and hard clip those bases off the read. 

/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/TrimStartingSequence \
INPUT=$prefix".unaligned_tagged_filtered.bam" \
OUTPUT=$prefix".unaligned_tagged_trimmed_smart.bam" \
OUTPUT_SUMMARY=$prefix".adapter_trimming_report.txt" \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5  

##PolyATrimmer 
##This Drop-seq program is the second sequence cleanup program designed to trim away 
##trailing polyA tails from reads.  It searches for at least 6 (NUM_BASES) contiguous 
##A’s in the read with 0 mismatches (MISMATCHES), and hard clips the read to remove 
##these bases and all bases 3’ of the polyA run. 
 
/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/PolyATrimmer \
INPUT=$prefix".unaligned_tagged_trimmed_smart.bam" \
OUTPUT=$prefix".unaligned_mc_tagged_polyA_filtered.bam" \
OUTPUT_SUMMARY=$prefix".polyA_trimming_report.txt" \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true 

##SamToFastq 
##Now that your data has had the cell and molecular barcodes extracted, 
##the reads have been cleaned of SMARTSeq primer and polyA tails, 
##and the data is now unpaired reads, it’s time to align.  
##To do this, we extract the FASTQ files using Picard’s ​SamToFastq​ program.  
 
/slipstream/home/mmariani/programs/miniconda3/bin/picard SamToFastq \
INPUT=$prefix".unaligned_mc_tagged_polyA_filtered.bam" \
FASTQ=$prefix".unaligned_mc_tagged_polyA_filtered.fastq" \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

##Alignment - STAR 
##We use ​STAR​ as our RNA aligner.  The manual for STAR can be found ​here​. 
##There are many potential aligners one could use at this stage, and it’s 
##possible to substitute in your lab’s favorite.  We haven’t tested other 
##aligners in methodical detail, but all should produce valid BAM files 
##that can be plugged into the rest of the process detailed here. 
##If you’re unsure how to create an indexed reference for STAR, 
##please read the STAR manual. 
##Below is a minimal invocation of STAR.  Since STAR contains a huge number 
##of options to tailor alignment to a library and trade off sensitivity 
##vs specificity, you can alter the default settings of the algorithm 
##to your liking, but we find the defaults work reasonably well for Drop-seq.  
##Be aware that STAR requires roughly 30 gigabytes of memory to align 
##a single human sized genome, and 60 gigabytes for our human/mouse reference. 
 
/slipstream/home/mmariani/programs/miniconda3/bin/STAR \
--genomeDir $ref_star \
--readFilesIn $prefix".unaligned_mc_tagged_polyA_filtered.fastq" \
--outFileNamePrefix $prefix 

##SortSam 
##This ​picard program​ is invoked after alignment, to guarantee that the output 
##from alignment is sorted in queryname order.  As a side bonus, the output file 
##is a BAM (compressed) instead of SAM (uncompressed.) 
 
/slipstream/home/mmariani/programs/miniconda3/bin/picard SortSam \
I=$prefix"Aligned.out.sam" \
O=$prefix"Aligned.sorted.bam"  \
SO=queryname \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"

##MergeBamAlignment 
##This Picard program merges the sorted alignment output from STAR (ALIGNED_BAM) 
##with the unaligned BAM that had been previously tagged with molecular/cell 
##barcodes (UNMAPPED_BAM).  This recovers the BAM tags that were “lost” during alignment.  
##The REFERENCE_SEQUENCE argument refers to the fasta metadata file.  
##We ignore secondary alignments, as we want only the best alignment from STAR 
##(or another aligner), instead of assigning a single sequencing read to 
##multiple locations on the genome. 

########### Note: .dict file needs to be in same folder and with same prefix as .fa 
########### for MergeBamAlignment to work, e.g. genome.fa and genome.dict in same folder.
/slipstream/home/mmariani/programs/miniconda3/bin/picard MergeBamAlignment \
REFERENCE_SEQUENCE=$ref_fasta \
UNMAPPED_BAM=$prefix".unaligned_mc_tagged_polyA_filtered.bam" \
ALIGNED_BAM=$prefix"Aligned.sorted.bam" \
OUTPUT=$prefix".merged.bam" \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false \
TMP_DIR="/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/tmp"
 
##/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/TagReadWithGeneExon \
##I=$prefix".merged.bam" \
##O=$prefix".star_gene_exon_tagged.bam" \
##ANNOTATIONS_FILE=$ref_flat \
##TAG=GE 

####Invocation (​​The call to TagReadWithGeneFunction is the same as TagReadWithGeneExon) 
/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/TagReadWithGeneFunction \
I=$prefix".merged.bam" \
O=$prefix".star_gene_exon_tagged.bam" \
ANNOTATIONS_FILE=$ref_flat  

##The below command does not finish completely 
##and produces the output and error:

##INFO	2020-10-08 14:21:02	UMIIterator	Sorting finished.
##INFO	2020-10-08 14:21:02	DetectBeadSubstitutionErrors	Gathering UMI counts per cell and filtering out UMI biased barcodes as appropriate
##INFO	2020-10-08 14:21:02	DetectBeadSubstitutionErrors	Finished gathering a list of cell barcodes to collapse
##WARNING	2020-10-08 14:21:02	DetectBeadSubstitutionErrors	Was not able to determine the polyTPosition from the input data.  This may be due to not having any properly prepared reads to evaluate.
##WARNING	2020-10-08 14:21:02	DetectBeadSubstitutionErrors	No barcodes found for collapse.  This means you have no cell barcodes with at least [20] transcripts and aren't UMI biased at the last base. You might have a problem in your input!
##INFO	2020-10-08 14:21:02	DetectBeadSubstitutionErrors	Starting Barcode Collapse of [0] barcodes
##INFO	2020-10-08 14:21:02	DetectBeadSubstitutionErrors	Barcode Collapse Complete - [0] barcodes collapsed

/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/DetectBeadSubstitutionErrors \
I=$data_dir"/"$(basename $prefix)".bam" \
O=$prefix".clean_subtitution.bam" \
OUTPUT_REPORT=$prefix".clean.substitution_report.txt" 
I=$prefix".bam" \

##Thus we cannot run the below command
##because a $prefix".clean_subtitution.bam"
##file is not produced

##/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/DetectBeadSynthesisErrors \
##I=$prefix".clean_subtitution.bam" \
##O=$prefix".clean.bam" \
##REPORT=$prefix".clean.indel_report.txt" \
##OUTPUT_STATS=$prefix".synthesis_stats.txt" \
##SUMMARY=$prefix".synthesis_stats.summary.txt" \
##PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC 

##DGE
##In this example, we extract the DGE for the top 100 most 
##commonly occurring cell barcodes in the
##aligned BAM, using CODING+UTR regions on the SENSE strand.

/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/DigitalExpression \
I=$prefix".star_gene_exon_tagged.bam" \
O=$prefix".star_gene_exon_tagged.dge.txt.gz" \
SUMMARY=$prefix".star_gene_exon_tagged.dge.summary.txt" \
NUM_CORE_BARCODES=100

##Cell Selection 
##A key question to answer for your data set is how many cells you want to extract from your BAM.  
##One way to estimate this is to extract the number of reads per cell, then plot the cumulative 
##distribution of reads and select the “knee” of the distribution.  
##We provide a tool to extract the reads per cell barcode in the Drop-seq software 
##called BAMTagHistogram.  This extracts the number of reads for any BAM tag in a BAM file, 
##and is a general purpose tool you can use for a number of purposes.  For this purpose, we extract the cell tag “XC”: 

/slipstream/home/mmariani/programs/Drop-seq_tools-2.4.0/BamTagHistogram \
I=$prefix".star_gene_exon_tagged.bam" \
O=$prefix".star_cell_readcounts.txt.gz" \
TAG=XC 

##I=$prefix".out_gene_exon_tagged.bam" \
##O=$prefix".out_cell_readcounts.txt.gz" \

##Once we run this program, a little bit of R code can create a cumulative distribution plot.  
##Here’s an example using the 100 cells data from the Drop-seq initial publication (Figures 3C and 3D): 

##R code from the drop-seq manual
##a=read.table("100cells_numReads_perCell_XC_mq_10.txt.gz", header=F, stringsAsFactors=F) 
##x=cumsum(a$V1) 
##x=x/max(x) 
##plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,500)) 

##I generated separate R files for producing the kneeplots and 
##connecting the counts files to Seurat.

done
