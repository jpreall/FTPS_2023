# FTPS_2023
## CSHL Plant Science Course, Single Cell Lectures, 2023

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

In the lab, the groups prepared eight PIPseq libraries from Maize root meristem, four treated with a peptide (FCP) and four control.

**cDNA Quality Control**
cDNA was checked by Jack and the students in the lab.

Ideally, these traces should contain a broad peak ranging from 1,000-10,000 base pairs, representing the full assortment of expressed full-length genes in the sample.

**Loading the Library**
These were pooled and sequenced at the CSHL NGS Core on a NextSeq2000, P2 flow cell, 100 cycle kit. Using the read format:

| Read1 | Index1 | Index2 | Read2 |
|---|---|---|---|
|54bp|8bp|8bp|67bp|

Here's what we would naively expect from the yield of this flow cell: based on how the samples were loaded

**500M read pairs / (8 samples * 6000 cells each) = 91,667 reads per cell**

Let's dive in:

### Illumina sequencing output
*(This is taken care of for you this year.  But here is some useful information about this step anyway, in case it's ever your responsibility to do the FASTQ generation step.):*

**Congratulations!  You didn't screw up an experiment.  Now you might have data.**

You submitted your libraries to your sequencing core and the run completed successfully.  If your core is nice enough to provide an Illumina quality score plot as part of your data delivery, it might look something like this:

![QC images](https://github.com/jpreall/FTPS_2023/blob/main/images/PIP_FASTQ_qual.png "Quality scores from NextSeq2000")
![QC images](https://github.com/jpreall/FTPS_2023/blob/main/images/PIP_FASTQ_base.png "Base composition from NextSeq2000")

Don't let those ugly spikes in the "% Base" (right panel) at the end of R1 and going on through the beginning of R4 worry you.  This is very typical.  R2 and R3 are the two index reads, which contain the sample barcodes, which in the case of our experiment is just a pool of 8 sequences.  Thus, it's totally expected that they have non-uniform base percentages. The region at the beginning of "R4" (which is the genomic insert) that has a stretch of non-uniform base utilization is also quite normal to see in libraries using TSOs in the RT chemistry (eg. 10X and PIPseq), and it comes from some common but tolerable artifacts of the library prep and sequencing. 

**Pipseeker** is Fluent Bio's free-to-use software that carries out mapping, and simple data analysis. *This is not meant to be run on your personal laptop.* For the purposes of this course, the Pipseeker steps will be carried out ahead of time before our interactive session, partly because it takes some time to complete (and that is if there is no queue on our cluster).  


-------

These are the files produced by Pipseeker:

```bash
$ ls -F1 pipseeker/Control-1-results

barcodes/
cell_calling/
clustering/
filtered_matrix/
logs/
metrics/
raw_matrix/
report.html
run_config.json

```

-------
### What are all those files?

The first thing you should look at is the `report.html`:

### START EDITING HERE
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
