# miniQuant: optimal and fast gene isoform abundance estimation
**miniQuant** features: 

1. **Optimal** use of long and/or short RNA-seq reads: transcript abundance estimation that can be applied to different data scenarios: long-read-alone and hybrid (long reads + short reads) integrating the strengths of both technologies.
2. **Fast** RNA-seq quantification: less than 15 minutes to analyze 40 million paired-end short reads + 5 million long reads on a standard laptop computer.
3. Calculate novel **K-value** metric: a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study.

To reproduce the results of our paper, download old version from [miniQuant v1.0](https://github.com/Augroup/miniQuant/tree/miniQuant_NBT_results_reproduce)

## Table of contents
  * [Dependency](#dependency)
  * [Installation](#installation)
  * [Usage](#usage)
    + [Data Preparation](#data-preparation-if-start-from-fastafastq-file)
    + [gene isoform quantification by miniQuant](#isoform-quantification-by-miniquant)
      - [1. If quantify using long reads data alone](#1-if-quantify-using-long-reads-data-alone)
      - [2. If quantify using short and long reads data in hybrid mode](#2-if-quantify-using-short-and-long-reads-data-in-hybrid-mode)
    + [Calculate K-value by miniQuant](#calculate-k-value-by-miniquant)

## Dependency
```
Linux operating system
```
## Installation
Simply download binary executable from [miniQuant]()

## Usage
*miniQuant* provides two options for gene isoform quantification:
1. quantify by long reads data alone.
2. quantify using short and long reads data in hybrid mode. <br>
A toy dataset example is provided in `example/`. Please following example command below for instruction.

#### 1. If quantify using long reads data alone

<details>
 <summary>Click me</summary>

 
**miniQuant** requires reference transcripts sequences in `FASTA` format (`-r`) and long-read RNA-seq sequences in plain or gzipped `FASTA/FASTQ` format (`-l`) as the input.
##### Example: quantify using long reads data (`example/LR.fasta.gz`) with reference transcripts sequences (e.g. `example/reference.fa.gz`), results in `miniQuant_LR_alone_res` folder
```
miniQuant quant -r example/reference.fa.gz -l example/LR.fasta.gz -t NUM_THREADS -o miniQuant_LR_alone_res
```
##### Available parameters
```
Required arguments:
  -r, --reference arg           Reference sequence file in plain or gzipped
                                FASTA format
  -l, --long_reads arg          Input long reads file in plain or gzipped
                                FASTA/FASTQ format. Leave blank if using only
                                short reads.
  -o arg, --output arg          The path of output folder.

Optional arguments:
  --long_reads_library_prep arg The library preparation for long reads.
                                Choices:[cDNA-ONT,dRNA-ONT,cDNA-PacBio]
                                (default: cDNA-ONT)
  -t arg, --threads arg         Number of threads. Default is 1.
```
##### Results explanation 
`miniQuant_LR_alone_res/abundance.tsv`
```
Transcript_id	TPM	Expected_num_long_reads
ENST00000379080.5	0	0
ENST00000379081.5	24154.8	8.57495
ENST00000379084.5	0	0
ENST00000379087.5	2.07187	0.000735514
ENST00000379089.5	0.00517853	1.83838e-06
ENST00000651358.1	1.23203	0.000437371
ENST00000445726.5	1.24601	0.000442334
ENST00000297620.8	37370.7	13.2666
ENST00000422409.5	0.000763295	2.7097e-07
ENST00000379078.1	9091.41	3.22745
ENST00000294244.9	812632	288.484
ENST00000540893.1	57424.8	20.3858
ENST00000535820.1	59322	21.0593
```
* `Transcript ID`: transcript ID given in the reference sequences
* `TPM`:  transcript abundance in TPM (Transcripts Per Kilobase Million) <br>
* `Expected_num_long_reads`: expected counts calculated by TPM and total number of long reads<br>
The result is a TSV file showing the abundance of each transcript, one transcript per line.
</details>


#### 2. If quantify using short and long reads data in hybrid mode

<details>
 <summary>Click me</summary>


* Integrates short and long reads RNA-seq reads from the same organism for better quantification performance. <br>
* In hybrid mode, **miniQuant** requires reference transcripts sequences in `FASTA` format (`-r`), long-read RNA-seq sequences in plain or gzipped `FASTA/FASTQ` format (`-l`), and short-read paired-end RNA-seq sequences in plain or gzipped `FASTA/FASTQ` format (`-1` and `-2`) as the input. <br>
##### Example: quantify using short reads (e.g. `example/SR_R1.fasta.gz` and `example/SR_R2.fasta.gz`) and long reads (e.g. `example/LR.fasta.gz`) with reference transcripts sequences (e.g. `example/reference.fa.gz`), results in `miniQuant_hybrid_res` folder
```
miniQuant quant -r example/reference.fa.gz -l example/LR.fasta.gz -1 example/SR_R1.fasta.gz -2 example/SR_R2.fasta.gz -t NUM_THREADS -o miniQuant_hybrid_res
```
##### Available parameters
```
Required arguments:
  -r, --reference arg           Reference sequence file in plain or gzipped
                                FASTA format
  -l, --long_reads arg          Input long reads file in plain or gzipped
                                FASTA/FASTQ format. Leave blank if using only
                                short reads.
  -1, --short_reads_pair_1 arg  Input short reads pair 1 in plain or
                                gzipped FASTA/FASTQ format. Leave blank if using
                                only long reads. (default: "")
  -2, --short_reads_pair_2 arg  Input short reads pair 2 in plain or
                                gzipped FASTA/FASTQ format. Leave blank if using
                                only long reads. (default: "")
  -o arg, --output arg          The path of output folder.

Optional arguments:
  --long_reads_library_prep arg The library preparation for long reads. Choices:[cDNA-ONT,dRNA-ONT,cDNA-PacBio] (default: cDNA-ONT)
  --short_reads_strandness arg  The strandness of short reads.          Choices:[unstranded,fr-stranded,rf-stranded] (default: unstranded)

                                *fr-strandred: Strand specific reads, first
                                read forward
                                *rf-stranded: Strand specific reads, first
                                read reverse
                                 
  -t arg, --threads arg         Number of threads. Default is 1.
```
##### Results explanation 
`miniQuant_res_hybrid/abundance.tsv`
```
Transcript_id	TPM	Expected_num_long_reads	Expected_num_short_read_pairs
ENST00000379080.5	0	0	0
ENST00000379081.5	11322.6	4.01951	46.1734
ENST00000379084.5	0	0	0
ENST00000379087.5	11140.4	3.95485	45.4306
ENST00000379089.5	0.00844047	2.99637e-06	3.44202e-05
ENST00000651358.1	0	0	0
ENST00000445726.5	0	0	0
ENST00000297620.8	26191.2	9.29789	106.808
ENST00000422409.5	4659.63	1.65417	19.002
ENST00000379078.1	12534	4.44955	51.1135
ENST00000294244.9	833640	295.942	3399.58
ENST00000540893.1	50781.5	18.0274	207.087
ENST00000535820.1	49731.1	17.6545	202.803
```

* `Transcript ID`: transcript ID given in the reference sequences
* `TPM`:  transcript abundance in TPM (Transcripts Per Kilobase Million) <br>
* `Expected_num_long_reads`: expected counts calculated by TPM and total number of long reads<br>
* `Expected_num_short_read_pairs`: expected counts calculated by TPM and total number of short read pairs<br>
The result is a TSV file showing the abundance of each transcript, one transcript per line.

</details>

### Calculate K-value by miniQuant

<details>
 <summary>Click me</summary>

 
**K-value** is a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study. K-value can be calculated given a gene isoforms annotation in GTF format
<div style="display: inline-flex; align-items: center;">
  <!-- Video Thumbnail -->
  <a href="https://www.youtube.com/watch?v=h9xTFpaJFgs" target="_blank" style="display: inline-block;">
    <img src="https://github.com/Augroup/miniQuant/blob/1d39426bc1820f0ca8bc42c3cee0d0e192a3ce4b/kvalue_intro_Qi1.png" style="width: 80%; display: block;">
  </a>

  <!-- Play Button -->
  <a href="https://www.youtube.com/watch?v=h9xTFpaJFgs" target="_blank" style="display: inline-block;">
    <img src="https://upload.wikimedia.org/wikipedia/commons/b/b8/YouTube_play_button_icon_%282013%E2%80%932017%29.svg" 
         style="width: 50px; height: auto; margin-left: 5px;">
  </a>
</div>

##### Example: calculate K-value given annotation in GTF format (e.g. `example/annotation.gtf`)
```
```
##### Results explanation 
K-value for each gene<br>
`miniQuant_kvalue/kvalues.out`
```
Gene	Chr	Num_isoforms	Kvalue
ENSG00000000003.15	chrX	5	14.263027941780145
```
* `Gene`: gene ID
* `Chr`: chromsosome ID
* `Num_isoforms`: number of isoforms in the gene
* `Kvalue`: K value <br>

*For gene that consists only short isoforms (i.e. all isoforms with length <150 bp), K-value will not be calculated and a `NA` value will be given.


</details>
