# miniQuant
**M͟i͟**xed Bayesian **n̲**etwork for **i̲**soform quantification (**miniQuant**) provides a highly-accurate bioinformatics tool for transcript abundance estimation.
**miniQuant** features: 
1. Novel **K-value** metric: a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study
2. **Mixed Bayesian network**: a novel mixed Bayesian network model for transcript abundance estimation that can be applied to different data scenarios: long-read-alone and hybrid (i.e., long reads plus short reads) integrating the strengths of both long reads and short reads.
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

## Table of contents
  * [Dependency](#dependency)
  * [Installation](#installation)
    + [[Recommended] Use Docker](#recommended-docker)
    + [[Recommended] Use Singularity](#recommended-singularity)
    + [[Not Recommended] Install from source](#not-recommended-install-from-source)
  * [Usage](#usage)
    + [Data Preparation](#data-preparation-if-start-from-fastafastq-file)
    + [Isoform quantification by miniQuant](#isoform-quantification-by-miniquant)
      - [1. If quantify using long reads data alone](#1-if-quantify-using-long-reads-data-alone)
      - [2. If quantify using short and long reads data in hybrid mode](#2-if-quantify-using-short-and-long-reads-data-in-hybrid-mode)
    + [Calculate K-value by miniQuant](#calculate-k-value-by-miniquant)

## Dependency
```
Linux operating system
```
The software has been tested with following software version
```
Python==3.9.7
minimap2==2.24
bowtie2==2.4.1
```
## Installation
It is recommended to use a Docker or Singularity to run the software.
### [Recommended] Docker

<details>
 <summary>Click me</summary>

#### Use [Docker](https://docs.docker.com/engine/install/)
```
# download and load docker image
wget https://miniquant.s3.us-east-2.amazonaws.com/miniQuant.tar.gz && docker load --input miniQuant.tar.gz && rm miniQuant.tar.gz
# run inside container
docker run -it --rm tidesun/miniquant:1.0 bash
cd / && source miniQuant/base/bin/activate
```

</details>


### [Recommended] Singularity

<details>
 <summary>Click me</summary>
 
#### Use [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#)
```
wget https://miniquant.s3.us-east-2.amazonaws.com/miniQuant.sif && singularity build --sandbox miniQuant_singularity miniQuant.sif && rm miniQuant.sif
singularity run -C --writable miniQuant_singularity bash
cd / && source miniQuant/base/bin/activate
```

</details>


### [Not Recommended] Install from source

<details>
 <summary>Click me</summary>
 
#### Dependency
```
Linux operating system with conda installed
```
```
conda create -n miniQuant python=3.8 openblas
conda activate miniQuant
wget -qO- https://miniquant.s3.us-east-2.amazonaws.com/miniQuant-1.0.0.tar.gz | tar xvz --one-top-level=miniQuant --strip-components 1
cd miniQuant
python -m venv base
source base/bin/activate
wget -qO- https://miniquant.s3.us-east-2.amazonaws.com/pretrained_models.tar.gz | tar xvz
pip install --upgrade pip
pip install setuptools==57.4.0
pip install -r requirements.txt
```
### Optional: install pretrained models for SIRV set-4 real data
```
cd miniQuant
wget -qO- https://miniquant.s3.us-east-2.amazonaws.com/SIRV_pretrained_models.tar.gz | tar xvz
```

</details>


## Usage
MiniQuant starts from Sequence Alignment/Map (SAM) file. For fasta/fastq file input, refer to [Data Preparation](#data-preparation-if-start-from-fastafastq-file) section, otherwise, refer to [Isoform quantification by miniQuant](#isoform-quantification-by-miniquant) section.
### Data Preparation (if start from fasta/fastq file)

<details>
 <summary>Click me</summary>
<br>
 
<b>Preparation:</b>
* install minimap2(v2.24) and bowtie2(v2.4.1)
<br>
<b>Required:</b>

* long reads alignment data mapped to reference genome in SAM format, example data can be found in `miniQuant/example/LR.sam`
* gene isoform annotation in GTF format, example data can be found in `miniQuant/example/annotation.gtf`
<br>

<b>Optional:</b>

* short reads alignment data mapped to reference transcriptome in SAM format, example data can be found in `miniQuant/example/SR.sam`
<br>
<b>Sequence alignment recommendation:</b>

##### use `minimap2` to map long reads data (e.g. [ENCFF714YOZ.fastq.gz](https://www.encodeproject.org/files/ENCFF714YOZ/@@download/ENCFF714YOZ.fastq.gz)) to reference genome (e.g. [GRCh38.primary_assembly.genome.fa](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz))
###### For dRNA-ONT data
```
minimap2 -a --MD -t 10 -N 0 -u f -x splice -o LR.sam 
GRCh38.primary_assembly.genome.fa ENCFF714YOZ.fastq.gz
```
###### For cDNA-ONT or cDNA-PacBio data
```
minimap2 -a --MD -t 10 -N 0 -x splice -o LR.sam 
GRCh38.primary_assembly.genome.fa ENCFF714YOZ.fastq.gz
```
##### use `Bowtie2` to map short reads data (e.g. paired end reads: [ENCFF892WVN.fastq.gz](https://www.encodeproject.org/files/ENCFF892WVN/@@download/ENCFF892WVN.fastq.gz) and [ENCFF481BLH.fastq.gz](https://www.encodeproject.org/files/ENCFF481BLH/@@download/ENCFF481BLH.fastq.gz) to reference transcriptome (e.g. [gencode.v39.transcripts.fa](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.transcripts.fa.gz))
```
bowtie2-build -f 
gencode.v39.transcripts.fa bowtie2_index

bowtie2 -f --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -p 10 -k 200 \
-x bowtie2_index -1 ENCFF892WVN.fastq.gz -2 ENCFF481BLH.fastq.gz > SR.sam

```

</details>

### Isoform quantification by miniQuant
*miniQuant* provides two options for isoform quantification:

1. quantify by long reads data alone.
2. quantify using short and long reads data in hybrid mode. <br>

An toy dataset example is provided in `miniQuant/example/`. Please following example command below for instruction.

#### 1. If quantify using long reads data alone

<details>
 <summary>Click me</summary>

 
miniQuant requires gene isoform annotation in `GTF` format (`-gtf`) and long reads sequence alignment mapped to the reference genome in `SAM` format (`-lrsam`) as the input.
##### Example: quantify using long reads data (`miniQuant/example/LR.sam`) with annotation (e.g. `miniQuant/example/annotation.gtf`), results in `miniQuant_res` folder
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
-gtf miniQuant/example/annotation.gtf \
-lrsam miniQuant/example/LR.sam \
-t 1 \
-o miniQuant_res

arguments:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of isoform annotation file in GTF format
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file mapping to reference genome.
  -t THREADS, --threads THREADS
                        Number of threads. Default is 1.
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory
```
##### Results explanation 
Isoform quantification abundance <br>
`miniQuant_res/Isoform_abundance.out`
```
Isoform	Gene	TPM
ENST00000373020.9	ENSG00000000003.15	710234.9711212328
ENST00000494424.1	ENSG00000000003.15	0.06848555891537092
ENST00000496771.5	ENSG00000000003.15	103773.58490566035
ENST00000612152.4	ENSG00000000003.15	3.2608726820945185e-20
ENST00000614008.4	ENSG00000000003.15	181274.39435547238
```
* `Isoform`: isoform ID
* `Gene`: gene ID
* `TPM`: isoform TPM <br>
The result is a TSV file showing the abundance of each gene isoform, one isoform per line.
</details>


#### 2. If quantify using short and long reads data in hybrid mode

<details>
 <summary>Click me</summary>


* Integrates short and long reads sequencing data from the same organism for better quantification performance. <br>
* In hybrid mode, miniQuant requires gene isoform annotation in `GTF` format (`-gtf`), long reads sequence alignment mapped to the reference genome in `SAM` format (`-lrsam`), and short reads sequence alignment mapped to reference transcriptome in `SAM` format (`-srsam`) as the input. <br>
* A pretrained machine learning model will be used for optimal intergration by simply set `-pretrained_model_path` to the long reads sequencing platform (i.e. `cDNA-ONT` for Oxford Nanopore cDNA sequencing, `cDNA-PacBio` for PacBio cDNA sequencing, and `dRNA-ONT` for Oxford Nanopore direct-RNA sequencing.
##### Example: quantify using short reads (e.g. `miniQuant/example/SR.sam`) and long reads data (e.g. `miniQuant/example/SR.sam`) by annotation (e.g. `miniQuant/example/annotation.gtf`), results in `miniQuant_res_hybrid` folder
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
-gtf miniQuant/example/annotation.gtf \
-lrsam miniQuant/example/LR.sam \
-srsam miniQuant/example/SR.sam \
--pretrained_model_path dRNA-ONT \
--EM_choice hybrid \
-t 1 \
-o miniQuant_res_hybrid

arguments:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of isoform annotation file in GTF format
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file mapping to reference genome.
  -srsam SHORT_READ_SAM_PATH, --short_read_sam_path SHORT_READ_SAM_PATH
                        The path of short read sam file mapping to reference transcriptome.
  --pretrained_model_path PRETRAINED_MODEL_PATH
                        The pretrained model path to identify the alpha. default: cDNA-ONT. \n
                        Can be one of the options [cDNA-ONT,dRNA-ONT,cDNA-PacBio] or file path of pretrained model.
  -t THREADS, --threads THREADS
                        Number of threads. Default is 1.
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory
```
##### Advanced parameters for hybrid quantification
```
optional arguments
  --eff_len_option EFF_LEN_OPTION
                        How to calculate the effective length [kallisto,RSEM]. Choose kallisto 
                        or RSEM to calculate the effective length in the same way as the 
                        corresponding method. Default is kallisto.
  --EM_SR_num_iters EM_SR_NUM_ITERS
                        Number of maximum iterations for EM algorithm. Default is 200.
```
##### Results explanation 
Isoform quantification abundance <br>
`miniQuant_res_hybrid/Isoform_abundance.out`
```
Isoform	Gene	Effective length	TPM
ENST00000373020.9	ENSG00000000003.15	3535.9141630901286	728571.217176296
ENST00000494424.1	ENSG00000000003.15	587.9141630901288	0.08438441577205032
ENST00000496771.5	ENSG00000000003.15	792.9141630901288	97566.37817660523
ENST00000612152.4	ENSG00000000003.15	3563.9141630901286	0.008307999564622307
ENST00000614008.4	ENSG00000000003.15	667.9141630901288	173862.31195468327
```
* `Isoform`: isoform ID
* `Gene`: gene ID
* `Effective length`: isoform effective length
* `TPM`: isoform TPM <br>
The result is a TSV file showing the abundance of each gene isoform, one isoform per line.

</details>

### Calculate K-value by miniQuant

<details>
 <summary>Click me</summary>

 
**K-value** is a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study. K-value can be calculated given a gene isoforms annotation in GTF format
##### Example: calculate K-value given annotation in GTF format (e.g. `miniQuant/example/annotation.gtf`)
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py cal_K_value \
-gtf miniQuant/example/annotation.gtf \
-t 1 \
-o miniQuant_kvalue

optional arguments:
  -t THREADS, --threads THREADS
                        Number of threads
  --sr_region_selection SR_REGION_SELECTION
                        SR region selection methods
                        [default:read_length][read_length,num_exons]
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
