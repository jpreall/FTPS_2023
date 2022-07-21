# FTPS_2022
## CSHL Plant Science Course, Single Cell Lectures, 2022

### [Link to pre-baked data](https://www.dropbox.com/s/izwly580v4r0nvc/FTPS22_data.tar.gz?dl=1) (637MB .tar.gz file)

### [Lecture Slides](https://www.dropbox.com/s/vezs9ppkwpudg97/Preall_FTPS22.pptx?dl=0) (158MB .pptx file)
-------

This tutorial will be a guide through the first few steps of primary data analysis:
1. [FASTQ generation with `cellranger mkfastq`](#section1)
2. [Making a custom genome reference with `cellranger mkref`](#section2)
3. [Mapping and count matrix generation with `cellranger count`](#section3)
4. [Combining two samples into a shared, normalized matrix with `cellranger aggr`](#section4)

-------
## Single Cell Lab: Demonstration Experiment

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/be/Root-tip-tag.png/440px-Root-tip-tag.png" width="300">

**Root Meristem**

In the lab, the groups prepared some 10X Genomics 3' Gene Expression libraries. We attempted to prepare 4 samples: 2 replicates each of treatment vs control.  Only one replicate of each passed QC during sample preparation, so we have a total of 2 sequenced datasets to analyze.

**cDNA Quality Control**
After the cDNA amplification step, we always check libraries on a BioAnalyzer or Tapestation
<img src="https://github.com/jpreall/FTPS_2022/blob/main/images/cDNA_traces.png" width="800">

Ideally, these traces should contain a broad peak ranging from 1,000-10,000 base pairs, representing the full assortment of expressed full-length genes in the sample. We can tell already that the Treated sample generated less overall cDNA than the Control, despite identical PCR conditions. This could suggest that we either captured fewer cells than intended, or that something about that sample didn't go as planned (for example, the silane cleanup beads were over-dried for that sample, leading to less efficient elution of the cDNA during cleanup. We will get a closer look when the data comes back.

**Loading the Library**
These were pooled and sequenced at the CSHL NGS Core on a NextSeq2000, P3 flow cell, 100 cycle kit. Using the read format:

| Read1 | Index1 | Index2 | Read2 |
|---|---|---|---|
|28bp|10bp|10bp|90bp|

Here's what we roughly expect from the yield of this flow cell:

**1.1B read pairs / (2 samples * 6000 cells each) = 91,667 reads per cell**

Let's dive in:

### Illumina sequencing output
*(This is taken care of for you this year.  But here is some useful information about this step anyway, in case it's ever your responsibility to do the FASTQ generation step.):*

**Congratulations!  You didn't screw up an experiment.  Now you might have data.**

You submitted your 10X libraries to your sequencing core and the run completed successfully.  If your core is nice enough to provide an Illumina quality score plot as part of your data delivery, it might look something like this:

![QC images](https://github.com/jpreall/FTPS_2022/blob/main/images/Maize_Seq_QC.png "QC data from NextSeq2000")

Don't let those ugly spikes in the "% Base" (right panel) at the end of R1 and going on through the beginning of R4 worry you.  This is very typical.  R2 and R3 are the two index reads, which contain the sample barcodes, which in the case of our experiment is just a pool of 2 sequences.  Thus, it's totally expected that they have non-uniform base percentages. The region at the beginning of "R4" (which is actually Read 2 of the genomic insert) that has a stretch of non-uniform base utilization is also quite normal to see in 10X Genomics libraries, and it comes from some common but tolerable artifacts of the library prep and sequencing.  If you look closely (or use a tool like [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), you'll even be able to determine the sequence of an abundant "contaminating" signature at the start of Read2:

```bash
AAGCAGTGGTATCAACGCAGAGTACATGGG ## Template Switch Oligo Sequence
```
As it turns out, this is the sequence of the 10X Template Switch Oligo (TSO).  Lots of reads contain with the TSO, due to artifacts of the library chemistry.  [Don't Panic.](https://en.wikipedia.org/wiki/Don%27t_Panic_(The_Hitchhiker%27s_Guide_to_the_Galaxy))

## Installing Cellranger

**Cellranger** is 10X Genomic's free-to-use software that carries out mapping, and primary data analysis for all 10X Genomics sequencing-based pipelines. *This is not meant to be run on your personal laptop.* It is an intensive, memory-hungry Linux application that is built to run on high-performance workstations or clusters. For the purposes of this course, the Cellranger steps will be carried out ahead of time before our interactive session, partly because it takes several hours to complete (and that is if there is no queue on our cluster).  But if you would like to know more about setting up Cellranger at your home institution, 10X has some helpful instructions:  [Cellranger Installation Help](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in)

## <a name="section1">Creating 10X-compatible FASTQ files with `cellranger mkfastq`</a>
[Cellranger mkfastq instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq)

Cellranger expects its input FASTQ files to be in a very specific format, which is ever-so-slightly different than the default format produced by Illumina's FASTQ generation pipeline.  For this reason, they bundled their own version of Illumina's `bcl2fastq` into a program called `mkfastq` that spits out files in a Cellranger-compatible format. You may never have to do this part yourself, since it is likely that NGS core will be the one generating your 10X- FASTQ files that will plug nicely into the subsequent `count` pipeline.  

```bash
cellranger mkfastq \
	--localcores=12 \
	--run=/path/to/basecalls/ \
	--samplesheet=/path/to/SampleSheet.csv \
```
**What the heck is a SampleSheet?**
A sample sheet tells the FASTQ generation pipeline how to break the reads out into separate folders based on their 8bp sample indices, and then how to name these files and organize them into folders.  Your Illumina sequencer will generate a SampleSheet.csv as part of the data generation process, but you may need to modify it in order to build a FASTQ file for 10X Genomics workflows:

**Example SampleSheet.csv**

```bash
[Header],,,,,,,
IEMFileVersion,4,,,,,,
Date,7/15/22,,,,,,
Workflow,GenerateFASTQ,,,,,,
Application,NextSeq FASTQ Only,,,,,,
Assay,TruSeq HT,,,,,,
Description,,,,,,,
Chemistry,Amplicon,,,,,,
,,,,,,,
[Reads],,,,,,,
28,,,,,,,
90,,,,,,,
[Settings],,,,,,,
,,,,,,,
[Data],,,:,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
300180,FTPS22_Control,,,SI-TT-A9,SI-TT-A9,FTPS22,
300183,FTPS22_Control,,,SI-TT-A10,SI-TT-A10,FTPS22,
```

### Dual-indexed libraries

If you use `cellranger mkfastq` to generate your FASTQs, it will translate `SI-TT-A9` into the pair of i5/i7 barcodes based on this [table provided by 10X Genomics.](https://cdn.10xgenomics.com/raw/upload/v1655151897/support/in-line%20documents/Dual_Index_Kit_TT_Set_A.csv) For example:

```bash
index_name,index(i7),index2_workflow_a(i5),index2_workflow_b(i5)
SI-TT-A9,AAGTGGAGAG,TTCCTGTTAC,GTAACAGGAA
SI-TT-A10,CGTGACATGC,ATGGTCTAAA,TTTAGACCAT
```
Here, 'workflow a' and 'workflow b' refer to the specific chemistry of the Illumina sequencer being used.  The NextSeq500 and NextSeq2000 read their i5 indices in opposite directions, so this barcode could be read out directly or as a the reverse complement, depending on the sequencer.  Cellranger will detect this automatically, so you needn't worry any further about it. Just remember to correcly specify your index in the SampleSheet using `SI-TT-xx` according to whichever well of the barcoding plate you used, and all will be magically taken care of.


### Single-indexed libraries (soon to be discontinued)
Originally, 10X used only a single index read, added at the final PCR step, to distinguish samples. This works fine for older sequencing instruments like the NextSeq500, but on the newer patterned flow-cell instruments like the NextSeq2000 and NovaSeq began to suffer from the problem of [index hopping](https://www.10xgenomics.com/blog/sequence-with-confidence-understand-index-hopping-and-how-to-resolve-it). While newer libraries use dual indexing, you may still encounter single-indexed strategies from time to time, so I will include here a description of how they work:


**Single-indexing using Pooled Barcodes:** 10X Genomics uses a clever trick to make barcoding simple.  Each sample barcode is actually a carefully chosen pool of four unique 8bp indices.  This is why you won't see anything that looks like `AAGTCTGA` in these sample sheets, but rather something that looks like `SI-GA-C7`

If you use `cellranger mkfastq` to generate your FASTQs, it will translate `SI-GA-A2` into a set of four 8bp barcodes based on this [table provided by 10X Genomics.](https://s3-us-west-2.amazonaws.com/10x.files/supp/cell-exp/chromium-shared-sample-indexes-plate.csv).   For example:

```bash
SI-GA-A1,GGTTTACT,CTAAACGG,TCGGCGTC,AACCGTAA
SI-GA-A2,TTTCATGA,ACGTCCCT,CGCATGTG,GAAGGAAC
SI-GA-A3,CAGTACTG,AGTAGTCT,GCAGTAGA,TTCCCGAC
...etc
```



-------

These are the files produced by Cellranger mkfastq from a NextSeq2000 sequencing run:

```bash
$ ls /path/to/fastqs/

FTPS22_Ctrl_S1_L001_I1_001.fastq.gz
FTPS22_Ctrl_S1_L001_I2_001.fastq.gz
FTPS22_Ctrl_S1_L001_R1_001.fastq.gz
FTPS22_Ctrl_S1_L001_R2_001.fastq.gz
FTPS22_Ctrl_S1_L002_I1_001.fastq.gz
FTPS22_Ctrl_S1_L002_I2_001.fastq.gz
FTPS22_Ctrl_S1_L002_R1_001.fastq.gz
FTPS22_Ctrl_S1_L002_R2_001.fastq.gz

```

-------
### What are all those files?

#### ...I1... = 10bp Sample index1 (i7) read:
```bash
>zcat FTPS22_Ctrl_S1_L001_I1_001.fastq.gz | head -n 4
@VH00553:6:AAALMHYHV:1:1101:18383:1000 1:N:0:AAGTGGAGAG+GTAACAGGAA
AAGTGGAGAG
+
CCCCCCCCCC

```

#### ...I2... = 10bp Sample index2 (i5) read:
```bash
>zcat FTPS22_Ctrl_S1_L001_I2_001.fastq.gz | head -n 4
@VH00553:6:AAALMHYHV:1:1101:18383:1000 2:N:0:AAGTGGAGAG+GTAACAGGAA
GTAACAGGAA
+
CCCCCCCCCC

```
These files contain the Sample Index read(s) information.  Note that the Cellranger mkfastq pipeline also writes the error-corrected sample index into the name line of the corresponding R1 and R2 reads.  See below.  

-------

#### ...R1... = Illumina Read 1 (aka the 'forward read').  This contains the cell barcode and UMI:
```bash
>zcat FTPS22_Ctrl_S1_L001_R1_001.fastq.gz | head -n 4
@VH00553:6:AAALMHYHV:1:1101:18383:1000 1:N:0:AAGTGGAGAG+GTAACAGGAA
NGCGTATAGGCTGGATGAAGTTAGTCGG
+
#CCC;CCCCCCCCCCCCCCCCCCCC-CC

```
In the current chemistry, this will be a 28bp read:
  * Bases 1-16: Cell Barcode
  * Bases 17-28: UMI
  
**Note:** Older (version 2) libraries used a 26bp read (16bp Cell barcode, 10bp UMI)

```bash
CATCAAGGTCAGATAAGGTCGATCCGTT
< Cell Barcode >
                <    UMI   >
```

The cell barcodes are only accepted if they are a close match to the curated list of possible barcodes generated by the synthesis chemistry. This list is buried within the Cellranger directory:

```bash
$ >head -n 4 /fake/path/cellranger-7.0.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCG
```

-------



#### ...R2... = Illumina Read 2.  This contains the gene body read (or Cell Hash / CITE-seq tags):

Our current sequencing method of choice is the NextSeq500 HighOutput SE75 flow cell.  These kits allow us to do 56 base pair reads for the gene body.  This is because the 75 cycle kit actually ships enough extra reagent to sequence the Illumina indices, plus a little more to spare.  In reality, we can squeeze out 92bp from a SE75 kit. 

**56 (gene body) + 28 (cell barcode / UMI) + 8 (Sample Index) = 92 bases**
 
```bash
>zcat FTPS22_Ctrl_S1_L001_R2_001.fastq.gz | head -n 4
@VH00553:6:AAALMHYHV:1:1101:18383:1000 2:N:0:AAGTGGAGAG+GTAACAGGAA
AAGCAGTGGTATCAACGCAGAGTACATGGCCAAGTACTACCTGGACGACACGGTGGACGTGGTCAAGATGCTGGACGGCCTGGCCAGCGC
+
CCCCCCCCCCCC-CCCCC;CCCCCCC-C;CCC-CCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCC-CCCCCCCCCCC;CCCCCCCCC-
```

-------
## <a name="section2">Building a custom genome reference with `cellranger mkref`</a>

Cellranger will map the reads to the reference genome that you specify and count digital gene expression according to the transcriptome model that was used during building of the reference.  For most users, we recommend downloading the pre-built references for the [human and mouse genomes provided by 10X](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).  But you weirdos study plants, so there is no convenient pre-built reference genome for you.  You'll have to make your own.

### [Cellranger mkref instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references)

## <a name="section3">Primary data analysis with `cellranger count`</a>

For each sample in your experiment, you'll need to run `cellranger count`.  [(detailed map of the pipeline)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count)



`cellranger count` is widely used to generate the primary count table, and has become a de facto standard despite there being a few alternatives.  It requires a linux machine with 8GB memory per CPU (to hold the human reference genome in memory), and is quite astoundingly slow and prone to crashing.  Mapping and counting a single sample typically takes 4-8 hours on a single multi-core node with 12-16 CPU cores and 128GB memory. 

If you have a HPC cluster for distributed computing, good for you.  You'll probably get to know your sysadmin really well as you repeatedly crash jobs, hog memory, and wreak havoc on your HPC cluster.  

Unfortunately, Cellranger is still the best option for generating a robust, industry-standard count matrix plus easily shareable results files and a handsome QC summary HTML page.  If you are a power user and simply want an accurate count matrix as quickly as possible, consider switching to [STAR](https://github.com/alexdobin/STAR).  Cellranger actually relies on an older verion of STAR for the basic transcript alignment, but the newest versions of STAR supports generating count matrices from single cell libraries such as 10X Genomics, or any custom barcoding format of your own design.  You won't get the nice shareable Loupe file, but you will get your data after lunch instead of tomorrow morning.  

Here is how to run Cellranger in local mode, if you have only a single workstation computer or if you want to restrict the run to a single node on your cluster.  I find this is slow but reliable:

```bash
#!/bin/sh

SAMPLE=FTPS22_Ctrl
TRANSCRIPTOME=/seq/CellRanger/references/Zea_Mays_v3_Mar2019


cellranger count \
  --id=$SAMPLE \
  --jobmode=local \
  --localcores=12 \
  --transcriptome=$TRANSCRIPTOME \
	--fastqs=/path/to/folder/containing/your/fastqs/ \
	--sample=$SAMPLE \
  
  ```
In the FASTQ example above, the sample ID specified was `FTPS22_Ctrl`.  `cellranger count` will search through your FASTQ folders to find files whose names match this sample ID, and automatically recognize Read1 vs Read2 based on the file name.  If you ran multiple flowcells and wish to combine the data for additional depth, provide two paths to each fastq folder,separated by a comma with no spaces:
  
```bash
  	--fastqs=/path/to/fastq/folder1/,/path/to/fastq/folder2/
```
Make sure you used the same sample ID when preparing the fastq file from both flowcells.

### Congratulations, you successfully ran cellranger count!
How many tries did it take you?

Your results will now be stored in a folder called `Sample_ID/outs`.  There will also be a bunch of other files containing diagnostic information about the run that you can dig through if you are a masochist.  

For instance, this is what the BAM (possorted_genome_bam.bam) file looks like:
```bash
# use samtools to view the first read on chromosome 1 of this indexed BAM file:
$samtools view possorted_genome_bam.bam 1 | head -n 1
VH00553:6:AAALMHYHV:1:1208:32471:20576	0	1	1	255	90M	*	0	0	GAATTCCAAAGCCAAAGATTGCATCAGTTCTGCTGCTATTTCCTCCTATCATTCTTTCTGATGTTGAAAATGATATTAAGCCTAGGATTC	CCCCCCCCCCCCCCCCC;CCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCC	NH:i:1	HI:i:1	AS:i:88	nM:i:0	RG:Z:FTPS22_Ctrl:0:1:AAALMHYHV:1	RE:A:I	xf:i:0	CR:Z:TCATGAGTCGTGAGAG	CY:Z:CCCCCCCCCCCCCCCC	CB:Z:TCATGAGTCGTGAGAG-1	UR:Z:TGTTCTAGTGCC	UY:Z:CCCCCCCCCCCC	UB:Z:TGTTCTAGTGCC

```
*Note: The pre-baked data that I linked above doesn't include the BAM file, since it's about 31GB in size.*

If you want to learn about what all these columns and tags mean, check out [this guide](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam)

In your `outs/` folder, you should see these files/folders:

```console
analysis/
cloupe.cloupe
filtered_feature_bc_matrix/
filtered_feature_bc_matrix.h5
metrics_summary.csv
molecule_info.h5
possorted_genome_bam.bam
possorted_genome_bam.bam.bai
raw_feature_bc_matrix/
raw_feature_bc_matrix.h5
web_summary.html
```
The first thing you should look at is the `web_summary.html`:
![Web Summary](https://github.com/jpreall/FTPS_2022/blob/main/images/maize_websumm.png "Web Summary Preview")

   
Who am I kidding, the first thing you did was download and view the pretty Loupe file
![Loupe snapshot](https://github.com/jpreall/FTPS_2022/blob/main/images/maize_loupe.png "Your awesome Loupe file")

That's ok, we all do it.  But seriously, let's look at the web summaries.  We're going to talk over what all those values mean in class. 

[Web Summary: Control](https://github.com/jpreall/FTPS_2022/blob/main/files/web_summary_Control.html)

[Web Summary: Treat](https://github.com/jpreall/FTPS_2022/blob/main/files/web_summary_Treat.html)

*Instrumental Break*

#### Sequencing Saturation
The single most useful piece of information stored in this summary is the estimate of sequencing saturation.  This will tell you how deeply you have sequenced these libraries, and whether it would be worth your time and money to add additional lanes of sequencing to identify new transcripts improve count numbers for differential expression.  

Total saturation is listed on the summary page, with a more thorough view in the `analysis` tab:

<img src="https://github.com/jpreall/FTPS_2022/blob/main/images/SeqTech_Saturation.png" width="500">

To a first approximation, we can see Control are pretty decent, but the Treated sample is not so good:
| Sample | Estimated Cells | Median UMIs/cell | Median genes/cell | Seq. saturation | Genome Mapping Rate |
| Control | 4,399 | 13,795 | 3,643 | 55.0% | 64.5% |
| Treat | 1,005 | 1,982 | 1,131 | 72.3% | 28.2% |

For both samples, we attempted to load ~10,000 cells.  Based on 10X's reported recovery rate of ~60%, we would have expected a yield of about 6,000 cells after barcoding. The yield on the Control sample is tolerable, but the Treated sample is quite poor. Not only that, but we see an awful low mapping rate.  Likely, the majority of these non-mappers are coming from adapter-dimer artifacts that carried through the library prep because there simply wasn't enough cDNA coming from successfully captured cells to make a good library.


#### The Data Matrix
Two other files that will be of extreme value to you are the actual data matrices.  Cellranger packages what would otherwise be an enormous data file into a clean, compressed hierarchical data (HDF5) file format.  It creates two versions: one that has been filtered of "empty" cells based on its [filtering algorithm](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview), and the raw matrix containing all 3M+ barcodes and their associated gene counts, regardless if they are likely to contain valid cells or not.

 * filtered_feature_bc_matrix.h5
 * raw_feature_bc_matrix.h5

The filtered and raw matrices are both also stored under a separate matrix market exchange (.mtx) file format along with separate .csv files listing the cell barcodes and gene names, which can be stiched together into a unified data matrix.  Look for these folders under `filtered_feature_bc_matrix` and `raw_feature_bc_matrix`, respectively.  You can open any of these files with [Seurat](https://satijalab.org/seurat/), or [Scanpy](https://scanpy.readthedocs.io/en/latest/), or potentially other 3rd party analysis packages.   


## <a name="section4"> Combining samples with `cellranger aggr`</a>

Let's combine both the control and the lard diet samples into a unified data matrix.  
Be careful not to accidentally bait a computation scientist into discussing the relative merits of the many different strategies for aggregating multiple data sets.  You will have to gnaw your foot off before they finally get to the punchline: 
*there is no single best way to jointly analyze multiple datasets*

10X Genomics has decided to sidestep the issue by providing a simple data aggregation pipeline that takes a conservative approach as a first step, and leaving the more sophisticated steps in your capable hands.  `cellranger aggr` combines two or more datasets by randomly downsampling (discarding) reads from the richer datasets, until all samples have approximately the same median number of reads per cell.  This helps mitigate one of the simplest and easy to fix batch-effects caused by sequencing depth, but will not correct for the zillions of other variables that injected unintended variation to your samples. 


#### Create an aggr.csv file:

First, tell cellranger which samples to aggregate by creating an aggr.csv file formatted thusly:

```bash
sample_id,molecule_h5
FTPS22_Ctrl,/fake/path/FTPS22/count/FTPS22_Ctrl/outs/molecule_info.h5
FTPS22_Treat,/fake/path/FTPS22/count/FTPS22_Treat/outs/molecule_info.h5
```
`cellranger aggregate` uses the `molecule_info.h5` file as the primary data source to do its downsampling.  This file contains rich data about each unique cDNA detected, including the number of duplicated or redundant reads mapping to a common UMI.  It is cleaner to downsample sequencing data based on this data rather than a simplified count matrix, which has discarded any information about the library complexity, PCR duplications, etc.  Cellranger uses this richer data source, but other tools seem to work with the final count matrix just fine.  Again, don't ask a bioinformation about it if you have children to feed some time today.

#### Run cellranger aggr:

```bash
cellranger aggr --id=FTPS22 \
	--jobmode=local \
	--csv=$BASEDIR/aggr.csv \
	--normalize=mapped \
	--localcores=16 \
	--localmem=64 
```
`cellranger aggr` is significantly less memory and cpu intensive than `cellranger count`.  If you are aggregating only a few samples, this should take less than an hour.  

Once it's done, you can view the web summary to see what was done to normalize the libraries:

<img src="https://github.com/jpreall/FTPS_2022/blob/main/images/maize_aggr_websumm.png" width="800">

Let's take a look at that aggr Loupe file.  Each sample is now stored as a separate category under "LibraryID":

<img src="https://github.com/jpreall/FTPS_2022/blob/main/images/maize_aggr.png" width="500">
