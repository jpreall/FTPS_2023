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


[Report: Control-1](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-1_report.html) 
[Report: Control-2](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-2_report.html) 
[Report: Control-3](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-3_report.html) 
[Report: Control-4](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-4_report.html) 

[Report: FCP-1](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-1_report.html) 
[Report: FCP-2](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-2_report.html) 
[Report: FCP-3](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-3_report.html) 
[Report: FCP-4](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-4_report.html) 

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

