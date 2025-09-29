# Improving gene isoform quantification with miniQuant
**miniQuant** features: 

1. **Optimal** use of long and/or short RNA-seq reads: transcript abundance estimation that can be applied to different data scenarios: long-read-alone and hybrid (long reads + short reads) integrating the strengths of both technologies.
2. **Fast** RNA-seq quantification: less than 15 minutes to analyze **unaligned** 40 million paired-end short reads + 5 million long reads on a standard laptop computer.
3. Calculate novel **K-value** metric: a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study.

Our [newest version](https://github.com/Augroup/miniQuant/releases) is recommended with faster speed and better performance. However, to reproduce the results on our [Nature Biotechnology paper](https://doi.org/10.1038/s41587-025-02633-9), download old version from [miniQuant v1.0](https://github.com/Augroup/miniQuant/tree/miniQuant_NBT_results_reproduce).
Feel free to run miniQuant [online](https://doi.org/10.24433/CO.9449101.v1) without installation!

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
1. Download [latest binary executable](https://github.com/Augroup/miniQuant/releases) (`wget https://github.com/Augroup/miniQuant/releases/download/latest/miniQuant_linux_latest.tar.gz`) and decompress by `tar -zxvf miniQuant_linux_latest.tar.gz`.
2. `cd miniQuant_linux && chmod +x miniQuant`
3. (<b>Optional</b>. Only if you want to directly call `miniQuant` in command line) `cp ./miniQuant /usr/local/bin; cp ./miniQuant ~/.local/bin`
3. Run `./miniQuant` 
### If your operating system doesn't have GLIBC 2.28 or later
1. Download [latest binary executable](https://github.com/Augroup/miniQuant/releases) (`wget https://github.com/Augroup/miniQuant/releases/download/latest/miniQuant_linux_latest.tar.gz`) and decompress by `tar -zxvf miniQuant_linux_latest.tar.gz`.
2. `cd miniQuant_linux && chmod +x miniQuant`
3. Run miniQuant under Docker or Singularity container with following commands: <br>
`docker run -i -t tidesun/miniquant:latest ./miniQuant` 
<br><br> OR <br><br> 
`singularity run docker://tidesun/miniquant:latest ./miniQuant`
<br><br> OR <br><br> 
`singularity run https://miniquant.s3.us-east-2.amazonaws.com/miniQuant_latest.sif ./miniQuant`

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
  --mem arg                     Max RAM usage in GB allowed when aligning 
                                the reads (default: 20.0)
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
| ENST00000379081.5 | 30326.9 | 9.85623 |
| ENST00000379084.5 | 0 | 0 |
| ENST00000379087.5 | 0.000665636 | 0.000000216332 |
| ENST00000379089.5 | 0 | 0 |
| ENST00000651358.1 | 2.68181 | 0.00087159 |
| ENST00000445726.5 | 2.76325 | 0.000898056 |
| ENST00000297620.8 | 31447 | 10.2203 |
| ENST00000422409.5 | 0.0521039 | 0.0000169338 |
| ENST00000379078.1 | 12294.7 | 3.99577 |
| ENST00000294244.9 | 807593 | 262.468 |
| ENST00000540893.1 | 56604.9 | 18.3966 |
| ENST00000535820.1 | 61728.4 | 20.0617 |

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
  --mem arg                     Max RAM usage in GB allowed when aligning 
                                the reads (default: 20.0)
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
| ENST00000379080.5 | 0.000143216 | 0.0000000465451 | 0.00000118102 | 3357 |
| ENST00000379081.5 | 10983.9 | 3.56975 | 89.2826 | 3309 |
| ENST00000379084.5 | 0 | 0 | 0 | 659 |
| ENST00000379087.5 | 11371 | 3.69557 | 93.2673 | 3339 |
| ENST00000379089.5 | 283.145 | 0.0920222 | 2.35789 | 3390 |
| ENST00000651358.1 | 9.77862 | 0.00317805 | 0.081936 | 3411 |
| ENST00000445726.5 | 10.9382 | 0.00355493 | 0.0916257 | 3410 |
| ENST00000297620.8 | 27378.1 | 8.89788 | 225.57 | 3354 |
| ENST00000422409.5 | 4603.86 | 1.49626 | 6.37848 | 564 |
| ENST00000379078.1 | 11825 | 3.84313 | 15.5698 | 536 |
| ENST00000294244.9 | 841246 | 273.405 | 3554.4 | 1720 |
| ENST00000540893.1 | 44896.1 | 14.5912 | 35.2918 | 320 |
| ENST00000535820.1 | 47392.2 | 15.4025 | 57.6272 | 495 |

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
      --not_normalize_entry     Whether NOT normalize region-isoform matrix
                                (A matrix) before calculating K-value
      --kvalue_entry_type arg   What kind of entry to use for
                                region-isoform matrix (A matrix).
                                Choices:[effective_length,binary] (default:
                                effective_length)
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
