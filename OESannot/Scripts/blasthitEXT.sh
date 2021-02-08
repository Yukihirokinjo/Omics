#!/bin/bash

#
# blasthitEXT.sh : Extract sequences which have confident blast hit aginst reference proteins
#
# Note: the IGS sequences must be extracted by IGSextract.bash. (sequence name structure (>ID:Start:End) is important to this script)


Query=""
Ref=""
name=""
QType=""
numCPU=""
Eval=""

usage_exit() {
	echo "Usage: blasthitEXT.sh  [ -q query_proteome/IGSs ] [ -g template_genome ] [ -r reference proteome ] [-t query_type (igs/prot)] [-o output_Prefix] [-n num._CPU]" 1>&2
	exit 1
}


while [ "$#" -gt 0 ]
do
	case "$1" in
		'-h' | '--help' )
			usage_exit
			;;
		'-q')
			if  [  -e "$2"  ]; then
				Query="$2"
				shift 2
			else
				echo "[Error] The input Query file is not found" 1>&2
				exit 1
			fi
			;;
		'-g')
			if  [  -e "$2"  ]; then
				Genome="$2"
				shift 2
			else
				echo "[Error] The input Genome file is not found" 1>&2
				exit 1
			fi
			;;
		'-r')
			if  [  -e "$2"  ]; then
				Ref="$2"
				shift 2
			else
				echo "[Error] The input Reference file is not found" 1>&2
				exit 1
			fi
			;;
		'-t')
			if  [ "$2" == "igs" -o  "$2" == "prot" ]; then
				Qtype="$2"
				shift 2
			else
				echo "[Error] The query sequence type (igs/prot) is not specifided" 1>&2
				exit 1
			fi
			;;
		'-o')
			if  [ -n  "$2" ]; then
				name="$2"
				shift 2
			else
				echo "[Error] The output prefix is not specifided" 1>&2
				exit 1
			fi
			;;
		'-e')
			if [ -z "$2" ]; then
				echo "PROGRAM: option requires an argument -- $1" 1>&2
				exit 1
			else
				if  	[ `expr "$2" : "[0-9]*$"` -gt 0  ]		||
					[ `expr "$2" : "[1-9][eE][-+][0-9]*$"` -gt 0 ]; then
					Eval=$2
					shift 2
				elif	[ `echo "$2 >= 0" | bc -l` -eq 1 ]; then
					Eval=`echo $2 | bc -l`
					shift 2
				else
					echo " Argument with option $1 should be a positive real number " 1>&2
					exit 1
				fi
			fi
			;;
                '-gc' | '--code' )
                        if [ -z "$2" ]; then
                                echo "PROGRAM: option requires an argument $1" 1>&2
                                exit 1
                        else
                                if  [ `expr "$2" : "[0-9]*$"` -gt 0  ]; then
                                        gCode="$2"
                                        shift 2
                                else
                                        echo " Argument with option $1 should be an integer " 1>&2
                                        exit 1
                                fi
                        fi
                        ;;
		'-n')
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
		*)
		echo "Invalid option "$1" " 1>&2
		usage_exit
		;;
	esac
done

#----------mandatory parameter check
if [ -z "$Query" ]; then
	echo "Input Query file is not specified" 1>&2
	usage_exit
fi
if [ -z "$Genome" ]; then
	echo "Input Genome file is not specified" 1>&2
	usage_exit
fi
if [ -z "$Ref" ]; then
	echo "Input Reference file is not specified" 1>&2
	usage_exit
fi
if [ -z "$Qtype" ]; then
	echo "Input Query type is not specified" 1>&2
	usage_exit
fi
if [ -z "$name" ]; then
	name="blahitEXT"
fi


#-----------------------------------

mkdir blaEXT
cp ${Genome} ${Query} blaEXT
cd blaEXT

pRefName=$(basename ${Ref})
RefName=${pRefName%.*}

# dual-blast (plastx/p and t/blastn)
# to get both position information (blasthit and genomic)

makeblastdb  -in ${Ref}    -dbtype prot -out tmp_ref  1> blastdb-ref.log
makeblastdb  -in ${Genome} -dbtype nucl -out tmp_genome  1> blastdb-genome.log

if [ "$Qtype"  == "igs" ]; then
  blastx -db tmp_ref -query ${Query} -query_gencode ${gCode:=11}  -evalue ${Eval:=1e-10} -outfmt 6  -out ${name}_${RefName}.blast  -num_threads ${numCPU:=1} -max_target_seqs 1 -max_hsps 1   1> blastx.log
  blastn -db tmp_genome -query ${Query} -evalue 1e-10 -outfmt 6 -out ${name}_GenomicPos.blast -num_threads ${numCPU:=1} -max_target_seqs 1 -max_hsps 1   1> blastn.log

  cat ${Query} > tmp_IGSs.fna
fi

if [ "$Qtype"  == "prot" ]; then
  blastp -db tmp_ref -query ${Query} -evalue ${Eval:=1e-10}  -outfmt 6  -out ${name}_${RefName}.blast  -num_threads ${numCPU:=1} -max_target_seqs 1 -max_hsps 1   1> blastp.log
  tblastn -db tmp_genome -query ${Query} -outfmt 6 -out ${name}_GenomicPos.blast -num_threads ${numCPU:=1} -max_target_seqs 1 -max_hsps 1   1> tblastn.log

  cat ${Query} > tmp_proteome.faa                                     
fi


if [ -s ${name}_${RefName}.blast ]; then

  # Extract ID and hit-positions
  awk '{print $1, $7, $8}' ${name}_${RefName}.blast > tmp_BLASThit.txt           # input for blasthitEXT.R
  # Extract genomic positions
  cat ${name}_GenomicPos.blast | awk '{print $1, $9, $10}' > tmp_TGBhit.txt  # input for blasthitEXT.R
  cat ${Genome}  >  tmp_genome.fna                                           # input for blasthitEXT.R

  # File initializing
  : > tmp_BLASThit.fna
  : > tmp_BLASThit.gff

   R --vanilla --slave --args  ${name}  ${Qtype} ${gCode:=11}  < $(which blasthitEXT.R)  

  # Sort output gff (Because the output gff/fasta entries are originally stored by BLAST hit score for each sequence)
  cat tmp_BLASThit.gff | sort -k 1,1 -k 4n,5  >  ${name}_BLASThit.gff
  cat tmp_BLASThit.txt | sort > ${name}_BLASThit.txt
  cat tmp_BLASThit.fna > ${name}_BLASThit.fna

cp ${name}_BLASThit.gff ..
else

  echo "##\n## No BLASThit was found.\n##\n"

fi
rm  tmp_*


cd ..
# END OF FILE
