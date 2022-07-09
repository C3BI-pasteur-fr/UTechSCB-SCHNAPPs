library(scPipe)
library(SingleCellExperiment)
library(Rsubread)
data_dir = tempdir()

# To process the data, we need the genome fasta file, gff3 exon annotation and a cell barcode annotation. The barcode annotation should be a `.csv` file with at least two columns, where the first column has the cell id and the second column contains the barcode sequence. We use ERCC spike-in genes for this demo. All files can be found in the `extdata` folder of the `scPipe` package:
  
# file path:
ERCCfa_fn = "/Volumes/CBUtechsZeus/bernd/levraud/index/fasta/marsseq.grcz11.fa"
ERCCanno_fn = "/Volumes/CBUtechsZeus/bernd/levraud/index/GRCz11.intervals.sorted.bed"
barcode_annotation_fn = system.file("extdata", "barcode_anno.csv", package = "scPipe")

# The read structure for this example is paired-ended, with the first longer read mapping to transcripts and second shorter read consisting of 6bp UMIs followed by 8bp cell barcodes. **NOTE**: by convention, `scPipe` always assumes `read1` refers to the read with the transcript sequence, which is usually the longer read. These data are also available in the in `extdata` folder:
  
fq_R1 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB26/Undetermined_S0_L001_R1_001ac.fastq.gz"
fq_R2 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB26/Undetermined_S0_L001_R2_001ac.fastq.gz"


# The pipeline starts with fastq file reformatting. We move the barcode and UMI sequences to the read name and leave the transcript sequence as is. This outputs a read name that looks like `@[barcode_sequence]*[UMI_sequence]#[readname]` ... The read structure of read 2 in our example dataset has the 8 bp long cell barcode starting at position 6 and the 6 bp long UMI sequence starting at the first position. So the read structure will be : `list(bs1=-1, bl1=0, bs2=6, bl2=8, us=0, ul=6)`. `bs1=-1, bl1=0` means we don't have an index in read 1 so we set a negative value as its start position and give it zero length. `bs2=6, bl2=8` means we have an index in read two which starts at position 6 in the read and is 8 bases in length. `us=0, ul=6` means we have a UMI at the start of read 2 which is 6 bases long. **NOTE**: we use a zero based index system, so the indexing of the sequence starts at zero.

# Seq_batch_ID	Run_name	Date	R1_design	I5_design	R2_design
# B6167	2107CSFCNsinfD39	120321	5I.4P.46M	NA	7W.8R


fq_R1 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB26/Undetermined_S0_L001_R1_001ac.fastq.gz"
fq_R2 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB26/Undetermined_S0_L001_R2_001ac.fastq.gz"
sc_trim_barcode(file.path(data_dir, "001ac.fastq"),
                fq_R1,
                fq_R2,
                read_structure = list(bs1=0, bl1=5, bs2=0, bl2=7, us=7, ul=8))

fq_R1 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB26/Undetermined_S0_L001_R1_001ae.fastq.gz"
fq_R2 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB26/Undetermined_S0_L001_R2_001ac.fastq.gz"

sc_trim_barcode(file.path(data_dir, "001ae.fastq"),
                fq_R1,
                fq_R2,
                read_structure = list(bs1=0, bl1=5, bs2=0, bl2=7, us=7, ul=8))


fq_R1 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB28/Undetermined_S0_L002_R1_001aa.fastq"
fq_R2 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB28/Undetermined_S0_L002_R2_001aa.fastq"

sc_trim_barcode(file.path(data_dir, "001aa.fastq"),
                fq_R1,
                fq_R2,
                read_structure = list(bs1=0, bl1=5, bs2=0, bl2=7, us=7, ul=8))

fq_R1 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB28/Undetermined_S0_L002_R1_001ae.fastq"
fq_R2 = "/Volumes/CBUtechsZeus/bernd/levraud/MARS-seq2.0_pipeline/scdb/raw_reads/SB28/Undetermined_S0_L002_R2_001ae.fastq"

sc_trim_barcode(file.path(data_dir, "001ab.fastq"),
                fq_R1,
                fq_R2,
                read_structure = list(bs1=0, bl1=5, bs2=0, bl2=7, us=7, ul=8))

system2("cat", args = paste(file.path(data_dir, "001aa.fastq"), file.path(data_dir, "001ab.fastq"),
                            file.path(data_dir, "001ac.fastq"), file.path(data_dir, "001ae.fastq"), sep=" "), 
        stdout = file.path(data_dir, "001aAll.fastq"))

# Next we align reads to the genome. This example uses `Rsubread` but any aligner that support RNA-seq alignment and gives standard BAM output can be used here.


Rsubread::buildindex(basename=file.path(data_dir, "ERCC_index"), reference=ERCCfa_fn)

Rsubread::align(index=file.path(data_dir, "ERCC_index"),
                readfile1=file.path(data_dir, "001aa.fastq"),
                nTrim5 = 4, minFragLength=25,type = "rna",
                nthreads = 8, sortReadsByCoordinates=T,
                output_file=file.path(data_dir, "out.aln.bam"), phredOffset=33)
