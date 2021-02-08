# OESannot ver. 0.8



## Requirements 

OESannot-1

  - prodigal
  - infernal
  - barnap
  - tRNAscan-SE

OESannot-2

- blast+
- diamond https://github.com/bbuchfink/diamond
- R
- seqinr(R)
- stringr(R)
- IRanges(R)

## Installation
```
tar -zxvf OESannot.tar.gz
cd OESannot
```

set an environment variable for OESannot
```
  echo "# OESannot" >> ~/.bashrc
  echo "export OESADIR=$(pwd)" >> ~/.bashrc

  ln -s $OESADIR/*.sh /path/to/$PATH_avilable_directory             
  ln -s $OESADIR/Scripts/*.{sh,R} /path/to/$PATH_avilable_directory 
```
or
```
  echo "export PATH=$PATH:$OESADIR" >> ~/.bashrc
  echo "export PATH=$PATH:$OESADIR/Scripts" >> ~/.bashrc
```

## USAGE 
```
OESannot-1.sh  -f input_genome.fna  [Optin: -c num_CPU ]

OESannot-2.sh  -f input_genome.fna -gff input_genome_merge.gff [Optin: -c num_CPU ]
```
NOTE: "input_genome_merge.gff" is the output of OESannot-1.sh.


##  Workflow 

OESannot-1 -> OESannot-2


### OESannot-1

  1-1 Initial ORF prediction (Prodigal)

  1-2 rRNA prediction (barrnap)

  1-3 tRNA prediction (tRNAscan-SE)

  1-4 Other non-coding RNA prediction (Infernal)

  1-5 Merge predictions


### OESannot-2

2-1 IGS ORF searech:
  Search additional ORFs from the region between predicted genes by using Blastx against COG database.
  The detected gene will be specified as "IGS_BLASThit" at 2nd colum on the corresponding row in output gff file.

2-2 ORF length evalation:
  Evaluate predicted ORFs length while referring to a length distribution of proteins in COG database which have same COGid.
  If a protein is evaluated as "short", "Short" flag will be put at the 6th colum on the corresponding row in xxx_OESa.gff file.

2-3 Pseudo-pseudo gene detection:
  Identify pseudo-pseudogenes.
  Information for identified pseudogenes will be stored with "PPgene" flag (this pipeline original) in the output gff file.
  The ORFs detected as pseudo-pseudogene will be modified by introducing artificial polymelase srippage to be a "complete ORF".
  The modified "complete ORF" an its translate will be written in .ffn and .faa files.

2-4 Domains completeness evaluation:
  Evaluate conserved domain completeness by using CDD search. 
  The information of completeness is stored at 6th colum.

2-5 Pseudogene detection:
  Identify pseudogenes based on the results from 2-2 and 2-4.
  Information for identified pseudogenes will be stored with "pseudo" flag in the output gff file.

## Output files
```
Intermediate
xxx_merge.gff : 
 An intermediate output from step 1-5. Gff files from ORF, rRNA, tRNA, and other ncRNAs predictions are merged into this file.

<Final>
xxx_OESA(or PGAP).gff :
 This gff file include all information about normal ORFs, pseudogens, pseudo-pseudogenes, and ncRNAs.
 NOTE: This gff file has the pipeline-specific format, especially at 3rd ("PPgene" flag), 6th (length|domain completeness information), and 9th colums.

xxx.ffn : 
 Multi fasta file contains nucleotide sequences of predicted CDSs and pseudo-pseudo genes. 
 Sequence of pseudo-pseudo gene is modifided to be complete ORF.

xxx.faa :
 Multi fasta file contains amino acids sequences of predicted CDSs and pseudo-pseudo genes.

xxx.frn :
 Multi fasta file contains nucleotide sequences of predicted RNAs.
```

## Reference

Under prep...



 2018.02.09
 Yukihiro Kinjo
