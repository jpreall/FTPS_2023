# FTPS_2023
## CSHL Plant Science Course, Single Cell Lectures, 2023

### [Link to pre-baked data](https://www.dropbox.com/scl/fi/fybbd2xbwbycju90l96e2/FTPS23_data.tar.gz?rlkey=gu68wxv32k3pcripjo8nk23b1&dl=0) (~400MB .tar.gz file)

### [Lecture Slides](https://www.dropbox.com/scl/fi/p3f86aalhgnwsjtjtgjk9/Preall_FTPS23.pptx?rlkey=5skruntfaco7w0fco8rgbfigw&dl=0) (158MB .pptx file)

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

**500M read pairs / (8 samples * 2000 cells targeted/sample) = 31,250 reads per cell**


### Illumina sequencing output
**Congratulations!  Some of your experiments sort of worked, maybe.**
*(This is taken care of for you this year.  But here is some useful information about this step anyway, in case it's ever your responsibility to do the FASTQ generation step.):*

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
### QC summary reports for all 8 samples:

[Report: Control-1](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-1_report.html)  
[Report: Control-2](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-2_report.html)  
[Report: Control-3](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-3_report.html)  
[Report: Control-4](https://github.com/jpreall/FTPS_2023/blob/main/files/files/Control-4_report.html)  
  
[Report: FCP-1](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-1_report.html)  
[Report: FCP-2](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-2_report.html)  
[Report: FCP-3](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-3_report.html)  
[Report: FCP-4](https://github.com/jpreall/FTPS_2023/blob/main/files/files/FCP-4_report.html)  

