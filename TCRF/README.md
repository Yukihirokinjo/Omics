#   TranscriptomeRarefaction ver. 1.1


## 0. Introduction



## 1. Prerequisites

TranscriptomeRarefaction depend on:


        seqtk                    (https://github.com/lh3/seqtk)
        Salmon                   ver. >1.3       (https://combine-lab.github.io/salmon/)
        Bowtie2 (Optional)       ver. >2.1       (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
        Samtools (Optional)      ver. >=1.9      (https://github.com/samtools/samtools)
        Trinity (Optional)       ver. >=2.8      (https://github.com/trinityrnaseq/trinityrnaseq)
        RSEM (Optional)          ver. >=1.3.3      (http://deweylab.github.io/RSEM)


        Add paths for the executables of these tools to your PATH.

The scripts work on the Linux environment. We have tested following environments:

        CentOS 7.2
        CentOS 8.0



## 2. Installation
No need to compile since TCSF and IMRA are bash wrapper scripts. Put them where you want to install them.
In the directory, change their permissions to be executable.

```
git clone https://github.com/Yukihirokinjo/Omics.git
cd Omics/TCRF

$ chmod u+x *.bash
```

Thereafter, add the path to your PATH.
For example..
```
$ echo 'export PATH=/path/to/Omics_dir/TCRF:$PATH' >> ~/.bashrc
$ source ~/.bashrc
```


## 3. Running TranscriptomeRarefaction.bash

### 3.1 Input data
```
Transcriptome assembly (fasta)
Illumina paired-end reads (fastq)
gene_trans_map file (Trinity.fasta.gene_trans_map)
```


### 3.2 Quick start

##### 
```
$ TranscriptomeRarefaction.bash   -i Trinity.fasta_file -m gene_trans_map_file -1 read_R1.fq -2 read_R2.fq -o output_directory 
```



### 3.3 Command line options
--------------------------------------------------------------------------------
##### 

        Mandatory:
        -i              <FILE>  Trinity assembly file (Trinity.fasta).

        -m              <FILE>  "gene_trans_map" file (Trinity.fasta.gene_trans_map). 

        -1              <FILE>  Input read file (forward).

        -2              <FILE>  Input read file (reverse).

        Optional:

        -Salmon         <BOOL>   Use Salmon for abundance estimation (default).

        -RSEM         <BOOL>   Use RSEM for abundance estimation.

        -o              <STR>   Output directory (default: "TCRF_Out_<current time>).

        -c              <INT>   Number of threads to be used for computation (default: 1).

        -s              <INT>   Random seed number for subsampling bam file  (default: 101).

        -t              <INT>   Threshold for the "expected read count" (default: 2).


--------------------------------------------------------------------------------

### 3.4 Output directories/files
--------------------------------------------------------------------------------
        ./OUTDIR/                        :Directory contains output files.

        ./OUTDIR/GeneCount.txt           :The result gene count file.

--------------------------------------------------------------------------------


