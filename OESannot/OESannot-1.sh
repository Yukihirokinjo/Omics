#!/bin/bash
#################################
# Genome Annot pipeline
# For Bacteria 1st
#################################

# genome(fasta, nucleotide; .fna/fasta)

usage_exit() {
	echo "Usage: OESannot-1.sh [-f Input_fasta ]  < Options.. >" 1>&2
	echo ""
	echo "Input data:"
	echo "  -f  | --fasta     : input genome sequence data in a fasta format file"
	echo ""
	echo "Options:"
	echo "  -v  | --version  : prints version of this script"
	echo "  -h  | --help     : prints this usage"
	echo "  -gc | --code     : specify genetic code (default: 11)"
	echo "  -c  | --cpu      : specify number of cpus (default: 1)"
	echo "  -o  | --out      : specify output directory prefix. The output directories will be 'XXX_predict' and 'XXX_curation'. (default: basename of input fasta file)"
	exit 1
}


Error_Check() {
    if [ "$?" -ne 0 ]; then
      echo "[Error] $1 failed. Please check the Messages avobe" 1>&2
      exit 1
    fi
}

version="0.3"

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
		'-f' | '--fasta')
			if  [  -e "$2"  ]; then
				genome="$2" 
				shift 2
			else
				echo "[Error] The input fasta file is not found" 1>&2
				exit 1	 
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
		'-c' | '--cpu')
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
		'-o' | '--out')
			if [ -z "$2" ]; then
				echo "PROGRAM: option requires an argument $1" 1>&2
				exit 1
			else
				OutDir="${2}_predict" 
				if  [ ! -e "${OutDir}"  ]; then
					echo "Output directory = ${OutDir}" 
					shift 2
				elif [ -e "${OutDir}" ] || [ -d "${OutDir}" ]; then
					echo "The specified output directory is already exist"
						echo -n 'Overwrite ? [y/n]: '
						read ans
							case $ans in
								y)
									rm -r ./"${OutDir}"
									;;
								n)
									printf "\nExit out of the IMRA\n\n" 
									exit 1
									;;
								*)
									echo "Press y(yes) or n(no)"
									;;
							esac
							echo
				fi
			fi
			;;
		*)
		echo "Invalid option "$1" " 1>&2 
		usage_exit
		;;
	esac
done

[ -z "$genome" ] && usage_exit

[ -z "$OESADIR" ] && echo "The OESADIR environment variable is not set. Please check README file."


Genome=$(basename ${genome})
GenName="${Genome%.*}"

[ -z "${OutDir}" ] && OutDir=${GenName}_predict

if [ -e "${OutDir}" ] || [ -d "${OutDir}" ]; then
	echo "The specified output directory is already exist"
	echo -n 'Overwrite ? [y/n]: '
	read ans
	case $ans in
		y)
			rm -r ./"${OutDir}"
			;;
		n)
			printf "\nExit out of the IMRA\n\n" 
			exit 1
			;;
		*)
			echo "Press y(yes) or n(no)"
			;;
	esac
	echo
fi

Rfam=${OESADIR}/db/rfam/Rfam.cm

mkdir ./${OutDir}
cp ${genome}  ./${OutDir}/
cd  ./${OutDir}

#------------------------------------------------1. Gene Prediction
##[CDS]

echo -e "\n## ------------------------ Prodigal ORFs Search ------------------ ##\n\n"
echo "Genetic Code =  ${gCode:=11}"

 prodigal -i ${Genome} -o ${GenName}_prodi.gff  -f gff -g ${gCode:=11}  -a ${GenName}_prodi.faa -d ${GenName}_prodi_mrna.fasta -s ${GenName}_start

# Prodigi conf < 60 (output: ${GenName}_prodiHC.gff )
# ProdiLCORFrm.bash -gff ${GenName}_prodi.gff -faa ${GenName}_prodi.faa -g  ${Genome} -ref  ${RefGenome} -o  ${GenName} -e 1e-10

##[rRNA]
echo -e "\n## ------------------------ BARNAP rRNAs Search ------------------ ##\n\n"

 barrnap  ${Genome} > ${GenName}_rRNAs.gff
# NOTE: modify barnap binary file Line119 regarding gff output ("barnap:0.8" -> "barnap-0.8")

##[tRNA]
echo -e "\n## ------------------------ tRNAscan-SE tRNAs Search ------------------ ##\n\n"

 tRNAscan-SE -B ${Genome} -o ${GenName}_tRNASE.out -f tRNASE.ss -m tRNASE.sum -F tRNASE.fpos 

##[ncRNA(infernal)]
echo -e "\n## ------------------------ Infernal Other ncRNAs Search ------------------ ##\n\n"

 cmscan --incE 1e-10 --tblout ${GenName}_infernal.out --fmt 2 --oskip ${Rfam}  ${Genome}


##--------- Annotations merge

echo -e "\n## ------------------------ Generating Result Files ------------------ ##\n\n"

 RNAout2gff.sh  ${GenName}_infernal.out  ${GenName}_tRNASE.out  ${GenName}_rRNAs.gff  ${GenName}

cat ${GenName}_prodi.gff  ${GenName}_RNAoutMerge.gff |  sort -k 1,1 -k 4n,4 | sed 's/ /-/g' > ${GenName}_merge.gff

cp ${GenName}_merge.gff ..
cd ..

echo -e "\n## ------------------------ Initial predictions has been finished. ------------------ ##\n"
echo -e "\n## ------------------------ Please proceed to the second pipeline. ------------------ ##\n"
echo -e "\n## ------------------------ Thank You ! \n\n"
## END OF FILE
