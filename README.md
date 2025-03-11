# miniQuant: fast and optimal gene isoform abundance quantification
**miniQuant** features: 

1. **Optimal** use of long and/or short RNA-seq reads: transcript abundance estimation that can be applied to different data scenarios: long-read-alone and hybrid (long reads + short reads) integrating the strengths of both technologies.
2. **Fast** RNA-seq quantification: less than 15 minutes to analyze **unaligned** 40 million paired-end short reads + 5 million long reads on a standard laptop computer.
3. Calculate novel **K-value** metric: a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study.

Our [newest version](https://github.com/Augroup/miniQuant/releases) is recommended. However, to reproduce the results on our [Nature Biotechnology paper](), download old version from [miniQuant v1.0](https://github.com/Augroup/miniQuant/tree/miniQuant_NBT_results_reproduce).

## Table of contents
  * [Dependency](#dependency)
  * [Installation](#installation)
  * [Usage](#usage)
    + [Gene isoform quantification by miniQuant](#1-if-quantify-using-long-reads-data-alone)
      - [1. If quantify using long reads data alone](#1-if-quantify-using-long-reads-data-alone)
      - [2. If quantify using short and long reads data in hybrid mode](#2-if-quantify-using-short-and-long-reads-data-in-hybrid-mode)
    + [Calculate K-value by miniQuant](#calculate-k-value-by-miniquant)

## Dependency
```
Linux operating system (tested on Red Hat 8.8)
```
## Installation
1. Download [latest binary executable](https://github.com/Augroup/miniQuant/releases) (`miniQuant_linux-v1.1.tar.gz`) and decompress by `tar -zxvf miniQuant_linux-v1.1.tar.gz`.
2. `cd miniQuant_linux && chmod +x miniQuant`
3. (Optional if you want to directly call `miniQuant`) `cp ./miniQuant /usr/local/bin; cp ./miniQuant ~/.local/bin`
3. Run `./miniQuant`

## Usage
*miniQuant* provides two options for gene isoform quantification:
1. quantify by long reads data alone.
2. quantify using short and long reads data in hybrid mode. <br>
A toy dataset example is provided in `example/`. Please following example command below for instruction.

### 1. If quantify using long reads data alone

<details>
 <summary>Click me</summary>

 
**miniQuant** requires reference transcripts sequences in `FASTA` format (`-r`) and long-read RNA-seq sequences in plain or gzipped `FASTA/FASTQ` format (`-l`) as the input.
##### Example: quantify using long reads data (`example/LR.fasta.gz`) with reference transcripts sequences (e.g. `example/reference.fa`), results in `miniQuant_LR_alone_res` folder
```
miniQuant quant -r example/reference.fa -l example/LR.fasta.gz -t 1 -o miniQuant_LR_alone_res
```
#### Available parameters
```
Required arguments:
  -r, --reference arg           Reference sequence file in plain or gzipped
                                FASTA format
  -l, --long_reads arg          Input long reads file in plain or gzipped
                                FASTA/FASTQ format.(default: "")

Optional arguments:
  -o arg, --output arg          The path of output folder. (default: ./miniQuant_res/)
  --long_reads_library_prep arg The library preparation for long reads.
                                Choices:[cDNA-ONT,dRNA-ONT,cDNA-PacBio]
                                (default: cDNA-ONT)
  -t arg, --threads arg         Number of threads. Default is 1.
```
#### Results explanation 
The result will be in TSV format (`miniQuant_LR_alone_res/abundance.tsv`) showing the abundance of each transcript, one transcript per line, with following columns:
* `Transcript ID`: transcript ID provided in the reference sequences(`--reference`)
* `TPM`: transcript relative abundance in TPM ([Transcripts Per Kilobase Million](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/)). <br>
* `Expected_num_long_reads`: expected counts of long reads, corresponding to the total number of long reads of input. <br>
<details>
<summary>Click me for example</summary>

| Transcript_id | TPM | Expected_num_long_reads |
| --- | --- | --- |
| ENST00000379080.5 | 0 | 0 |
| ENST00000379081.5 | 24154.8 | 8.57495 |
| ENST00000379084.5 | 0 | 0 |
| ENST00000379087.5 | 2.07187 | 0.000735514 |
| ENST00000379089.5 | 0.00517853 | 0.00000183838 |
| ENST00000651358.1 | 1.23203 | 0.000437371 |
| ENST00000445726.5 | 1.24601 | 0.000442334 |
| ENST00000297620.8 | 37370.7 | 13.2666 |
| ENST00000422409.5 | 0.000763295 | 0.00000027097 |
| ENST00000379078.1 | 9091.41 | 3.22745 |
| ENST00000294244.9 | 812632 | 288.484 |
| ENST00000540893.1 | 57424.8 | 20.3858 |
| ENST00000535820.1 | 59322 | 21.0593 |

</details>
</details>


### 2. If quantify using short and long reads data in hybrid mode

<details>
 <summary>Click me</summary>


* Integrates short and long reads RNA-seq reads from the same organism for better quantification performance. <br>
* In hybrid mode, **miniQuant** requires reference transcripts sequences in `FASTA` format (`-r`), long-read RNA-seq sequences in plain or gzipped `FASTA/FASTQ` format (`-l`), and short-read paired-end RNA-seq sequences in plain or gzipped `FASTA/FASTQ` format (`-1` and `-2`) as the input. <br>
#### Example: quantify using short reads (e.g. `example/SR_R1.fasta.gz` and `example/SR_R2.fasta.gz`) and long reads (e.g. `example/LR.fasta.gz`) with reference transcripts sequences (e.g. `example/reference.fa`), results in `miniQuant_hybrid_res` folder
```
miniQuant quant -r example/reference.fa -l example/LR.fasta.gz -1 example/SR_R1.fasta.gz -2 example/SR_R2.fasta.gz -t 1 -o miniQuant_hybrid_res
```
#### Available parameters
```
Required arguments:
  -r, --reference arg           Reference sequence file in plain or gzipped
                                FASTA format
  -l, --long_reads arg          Input long reads file in plain or gzipped
                                FASTA/FASTQ format.(default: "")
  -1, --short_reads_pair_1 arg  Input short reads pair 1 in plain or
                                gzipped FASTA/FASTQ format. Leave blank if using
                                only long reads. (default: "")
  -2, --short_reads_pair_2 arg  Input short reads pair 2 in plain or
                                gzipped FASTA/FASTQ format. Leave blank if using
                                only long reads. (default: "")

Optional arguments:
  -o arg, --output arg          The path of output folder. (default: ./miniQuant_res/)
  --long_reads_library_prep arg The library preparation for long reads. Choices:[cDNA-ONT,dRNA-ONT,cDNA-PacBio] (default: cDNA-ONT)
  --short_reads_strandness arg  The strandness of short reads.          Choices:[unstranded,fr-stranded,rf-stranded] (default: unstranded)

                                *fr-stranded: Strand specific reads, first
                                read forward
                                *rf-stranded: Strand specific reads, first
                                read reverse
                                 
  -t arg, --threads arg         Number of threads. Default is 1.
```
#### Results explanation 
The result will be in TSV format (`miniQuant_res_hybrid/abundance.tsv`) showing the abundance of each transcript, one transcript per line, with following columns:
* `Transcript ID`: transcript ID provided in the reference sequences(`--reference`)
* `TPM`: transcript relative abundance in TPM ([Transcripts Per Kilobase Million](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/)). It is calculated by integrating both short and long reads. <br>
* `Expected_num_long_reads`: expected counts of long reads. It is calculated by integrating both short and long reads, corresponding to the total number of long reads of input. <br>
* `Expected_num_short_read_pairs`: expected counts of short read pairs, corresponding to the total number of short read pairs of input. <br>
* `Effective_length`: effective length of each transcript.
<details>
<summary>Click me for example</summary>

| Transcript_id | TPM | Expected_num_long_reads | Expected_num_short_read_pairs | Effective_length |
| --- | --- | --- | --- | --- |
| ENST00000379080.5 | 0.000615982 | 0.000000218674 | 0.00000251197 | 3357 |
| ENST00000379081.5 | 11328.6 | 4.02167 | 46.1982 | 3309 |
| ENST00000379084.5 | 0 | 0 | 0 | 659 |
| ENST00000379087.5 | 9749.14 | 3.46095 | 39.757 | 3339 |
| ENST00000379089.5 | 1604.44 | 0.569576 | 6.54291 | 3390 |
| ENST00000651358.1 | 31.5397 | 0.0111966 | 0.128619 | 3411 |
| ENST00000445726.5 | 34.3598 | 0.0121977 | 0.140119 | 3410 |
| ENST00000297620.8 | 27464.6 | 9.74994 | 112.001 | 3354 |
| ENST00000422409.5 | 2634.61 | 0.935286 | 10.7439 | 564 |
| ENST00000379078.1 | 11953.4 | 4.24347 | 48.7461 | 536 |
| ENST00000294244.9 | 842092 | 298.943 | 3434.05 | 1720 |
| ENST00000540893.1 | 45946.1 | 16.3109 | 187.368 | 320 |
| ENST00000535820.1 | 47160.9 | 16.7421 | 192.322 | 495 |

</details>
</details>

### Calculate K-value by miniQuant

<details>
<summary>Click me</summary>

 
**K-value** is a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study. K-value can be calculated given a gene isoforms annotation in GTF/GFF3/genePred format.
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

#### Example: calculate K-value given annotation in GTF/GFF3/genePred format (e.g. `example/annotation.gtf`)
```
miniQuant kvalue -a example/annotation.gtf -o miniQuant_kvalue -t 1
```
#### Available parameters
```
Required arguments:
  -a, --annotation arg  Gene isoform annotation file in GTF, GFF or 
                        genePred format

 Optional arguments:
  -o, --output arg      The path of output folder (default: ./miniQuant_kvalue/)
  -t, --threads arg             Num of threads (default: 1)
      --short_reads_mean_fragment_length arg
                                Mean value of short reads fragment lengths 
                                (default: 235.0)
```
#### Results explanation 
The result will be in TSV format (`miniQuant_kvalue/kvalues.tsv`) showing the K-value of each gene, one gene per line, with following columns:
* `Gene`: gene ID
* `K-value`: K-value. Larger K-value indicates higher quantification error. <br>

<details>
<summary>Click me for example</summary>

| Gene_id | K-value |
| --- | --- |
| ENSG00000164970.15 | 331.422233 |
| ENSG00000168005.9 | 1.320074 |

</details>



*For gene that consists only short isoforms (i.e. all isoforms with length less than `--short_reads_mean_fragment_length`), K-value will not be calculated and a `NA` value will be given.


</details>
