#!/bin/bash
#
# TranscriptomeRarefaction.bash 
#
# depends on:
#   Trinity
#   RSEM
#   bowtie2
#   Samtools
#   Salmon
#   Seqtk

InFasta=""
InMap=""
Read1=""
Read2=""
outDir=""
Seed=""
EstMethod=""
numCPU=""
Info_SamBam=""
version=2.0

usage_exit() {
cat << EOS
  # 
  # TranscriptomeRarefaction  version-${version} 
  #
  # Usage: TranscriptomeRarefaction.bash  -i Trinity.fasta_file -m gene_trans_map_file -1 read_R1.fq -2 read_R2.fq -o output_directory  
  #
  # Arguments (mandatory): 
  # -i Transcript sequence file (Trinity.fasta) 
  # -m gene_trans_map file (Trinity.fasta.gene_trans_map)
  # -1 read_R1   
  # -2 read_R2   
  # Arguments (Optional): 
  # --RSEM		use RSEM for abundance estimation (default)
  # --Salmon	use Salmon for abundance estimation
  # -c num_CPUs 
  # -s random seed number (default: 101)
  # -t read count threshold (default: 2)
  #
  # Requirements:
  # - Trinity 
  # - bowtie2
  # - RSEM
  # - Salmon
  # - Seqtk
  # - Samtools 
EOS
  exit 1
}

Error_Check() {
    if [ "$?" -ne 0 ]; then
      echo "[Error] $1 failed. Please check the Messages avobe" 1>&2
      exit 1
    fi
}


while [ "$#" -gt 0 ]
do
  case "$1" in
    '-v' | '--version' )
      echo "Ver._$version" 
      exit 1
      ;;
    '-h' | '--help' )
      usage_exit
      ;;
    '-i')
      if  [  -e "$2"  ]; then
        InFasta="$2" 
        shift 2
      else
        echo "[Error] The input assembly file is not found" 1>&2
        exit 1   
      fi
      ;;
    '-m')
      if  [  -e "$2"  ]; then
        InMap="$2" 
        shift 2
      else
        echo "[Error] The input gene-trans mapping file is not found" 1>&2
        exit 1   
      fi
      ;;
    '-1')
      if  [  -e "$2"  ]; then
        readR1="$2" 
        shift 2
      else
        echo "[Error] The input read file is not found" 1>&2
        exit 1   
      fi
      ;;
    '-2')
      if  [  -e "$2"  ]; then
        readR2="$2" 
        shift 2
      else
        echo "[Error] The input read file is not found" 1>&2
        exit 1   
      fi
      ;;
    '-o')
      if  [ ! -e "$2"  ]; then
        outDir="$2" 
        shift 2
      else
        echo "[Error] The output directory is already exist" 1>&2  
        exit 1   
      fi
      ;;
    '-t')
      if [ -z "$2" ]; then
        echo "PROGRAM: option requires an argument $1" 1>&2
        exit 1
      else
        if  [ `expr "$2" : "[0-9]*$"` -gt 0  ]; then
          gtCount="$2" 
          shift 2
        else
          echo " Argument with option $1 should be an integer " 1>&2
          exit 1
        fi
      fi
      ;;
    '-RSEM')
		EstMethod="RSEM"
		shift 1
	  ;;
     '-Salmon')
		EstMethod="Salmon"
		shift 1
	  ;;
    '-c')
      if [ -z "$2" ]; then
        echo "PROGRAM: option requires an argument $1" 1>&2
        exit 1
      else
        if  [ `expr "$2" : "[0-9]*$"` -gt 0  ]; then
          numCPU="$2" 
          shift 2
        else
          echo " Argument with option $1 should be an integer " 1>&2
          exit 1
        fi
      fi
      ;;
    '-s')
      if [ -z "$2" ]; then
        echo "PROGRAM: option requires an argument $1" 1>&2
        exit 1
      else
        Seed="$2"
        shift 2
      fi
      ;;
    '-Info')
      Info_SamBam="on"
      shift 1
      ;;
    *)
    echo "Invalid option "$1" " 1>&2 
    usage_exit
    ;;
  esac
done

if [ -z "$readR1" ]; then
  echo "Read file is not specified" 1>&2
  usage_exit
fi
if [ -z "$readR2" ]; then
  echo "Read file is not specified" 1>&2
  usage_exit
fi

# default output directory
[ -z $outDir ] && outDir="TCRF_OUT_`date +%Y%m%d`_`date +%H%M`"

if [ -z "$EstMethod" ]; then
  EstMethod="Salmon"
fi


## prep output directory and file
mkdir -p ${outDir}
Error_Check Create_OutDir
echo "Prop. NumGenes" > ${outDir}/GeneCount.txt

############################################################ if EstMethod == RSEM 

if [ ${EstMethod} == "RSEM" ] ; then
	##  RSEM  1st (for All data)
	$TRINITY_HOME/util/align_and_estimate_abundance.pl \
	              --transcripts ${InFasta}  \
	              --seqType fq \
	              --left ${readR1}  --right ${readR2}  \
	              --est_method RSEM   \
	              --aln_method bowtie2 \
	              --prep_reference \
	              --gene_trans_map  ${InMap} \
	              --thread_count ${numCPU:=1} \
	              --output_dir  ${outDir}/RSEM_All

	Error_Check  align_and_estimate_abundance

	mkdir ${outDir}/RSEM_All/tmp
	mv Trinity.fasta.{RSEM,bowtie2}.* ${outDir}/RSEM_All/tmp



	##  RSEM  2nd (for subsampling)
	for i in 0.01 0.05 0.1 0.2 0.3 0.5 0.75 ; do
	  # set iteration variables
	  p=`echo "scale=0; ${i} * 100" | bc`   # 0.01 -> 1.00
	  Int=${p%.*}                           # 1.00 -> 1 (%)
	  Frac=$(printf %02d $Int);             # 1 -> 01 (%)
	 
	  # subsample bam file ( -s "SeedNum" + "." + "Proportion" )
	  samtools view -s ${Seed:=101}.${Frac} -b ${outDir}/RSEM_All/bowtie2.bam > ${outDir}/p_${Frac}.sam  
	  Error_Check  SAMtools

	  # RSEM for each subset
	  rsem-calculate-expression --num-threads ${numCPU:=1} --sam ${outDir}/p_${Frac}.sam --paired-end ${outDir}/RSEM_All/tmp/Trinity.fasta.RSEM  ${outDir}/p_${Frac}
	  Error_Check  RSEM

	  # count number of expressed genes (above the TPM threshold)
	  NumGene=$( awk -v gtC=${gtCount:=2} '{if($5 > gtC) print $0}' ${outDir}/p_${Frac}.RSEM.genes.results | wc -l ) 
	  Error_Check GeneCount

	  # record results in output file
	  echo "${Frac} ${NumGene}" >>  ${outDir}/GeneCount.txt

	  if [ -z "$Info_SamBam" ]; then
	    rm  ${outDir}/p_${Frac}.{sam,bam}    
	  fi

	done

	# record original gene count (all)
	NumGene=$( awk -v gtC=${gtCount:=2} '{if($6 > gtC) print $0}' ${outDir}/RSEM_All/RSEM.genes.results | wc -l ) 
	Error_Check GeneCount
	echo "100 ${NumGene}" >>  ${outDir}/GeneCount.txt

############################################################ if EstMethod == Salmon
else
	##  Slmon  1st (for All data)
	salmon index -p ${numCPU:=1} -t ${InFasta} -i ref_idx
	Error_Check  SalmonInDex

	salmon quant -i ref_idx -l A -1 ${readR1} -2 ${readR2}  -p ${numCPU:=1} -o ${outDir}/Salmon_All --validateMappings
	Error_Check  SalmonQuant

	mkdir ${outDir}/Salmon_All/tmp

	##  Salmon  2nd (for subsampling)
	for i in 0.01 0.05 0.1 0.2 0.3 0.5 0.75 ; do
	  # set iteration variables
	  p=`echo "scale=0; ${i} * 100" | bc`   # 0.01 -> 1.00
	  Int=${p%.*}                           # 1.00 -> 1 (%)
	  Frac=$(printf %02d $Int);             # 1 -> 01 (%)
	
	  # subsample bam file ( -s "SeedNum" + "." + "Proportion" )
	  seqtk sample -s ${Seed:=101} ${readR1}  ${i}  > ${outDir}/p_${i}_R1.fastq 
	  seqtk sample -s ${Seed:=101} ${readR2}  ${i}  > ${outDir}/p_${i}_R2.fastq 
	  Error_Check  Seqtk

	  # salmon for each subset
	  salmon quant -i ref_idx -l A -1 ${outDir}/p_${i}_R1.fastq -2 ${outDir}/p_${i}_R2.fastq -p ${numCPU:=1} -o ${outDir}/p_${Frac} --validateMappings
	  Error_Check  SalmonQuant

	  # count number of expressed genes (above the TPM threshold)
	  NumGene=$( awk -v gtC=${gtCount:=2} '{if($5 > gtC) print $0}' ${outDir}/p_${Frac}/quant.sf | wc -l ) 
	  Error_Check GeneCount

	  # record results in output file
	  echo "${Frac} ${NumGene}" >>  ${outDir}/GeneCount.txt

	  rm ${outDir}/p_${i}_R{1,2}.fastq

	done

	NumGene=$( awk -v gtC=${gtCount:=2} '{if($5 > gtC) print $0}' ${outDir}/Salmon_All/quant.sf | wc -l ) 
	Error_Check GeneCount
	echo "100 ${NumGene}" >>  ${outDir}/GeneCount.txt


fi



#EOF
