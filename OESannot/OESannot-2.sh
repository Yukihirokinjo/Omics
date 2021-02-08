#!/bin/bash
#################################
# Genome Annot pipeline
# For Bacteria 2nd
#################################
#
# USAGE: OESannot-2.sh  [ -f input_genome.fna ] [ -gff BORF-1_output.gff ]  [ -d BORF directory ] [ -c numCPU ]
#
# genome(fasta, nucleotide; .fna/fasta)
# refgenome(fasta, amino acid; .faa)
# NOTE: File for length information of reference sequences must be prepared. -> <DBid> <Avg.Length> <sd> <MaxLength> <MinLength>

usage_exit() {
	echo "Usage: OESannot-2.sh [ -f input_genome.fna ] [ -gff BORF-1_output.gff ]   < Options.. >" 1>&2
	echo ""
	echo "Input data:"
	echo "  -f  | --fasta       : input genome sequence data in a fasta format file"
	echo "  -gff                : input gff file generated by BORF-1"
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
		'-f' )
			if  [  -e "$2"  ]; then
				genome="$2"
				shift 2
			else
				echo "[Error]The input fasta file is not found" 1>&2
				exit 1
			fi
			;;
		'-gff' | '--gff' )
			if  [  -e "$2"  ]; then
				gff="$2"
				shift 2
			else
				echo "[Error] The gff file is not found" 1>&2
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
		'-o' | '--out')
			if [ -z "$2" ]; then
				echo "PROGRAM: option requires an argument $1" 1>&2
				exit 1
			else
				OutDir="${2}_curation" 
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

[ -z "$OESADIR" ] && echo "'OESADIR' environment variable is not set. Please check README file."

[ -z "$gff" ] && usage_exit
Gff=`basename ${gff}`

CDD=${OESADIR}/db/cdd/CDDncbi_Pfam
CDDid=${OESADIR}/db/cdd/cddid.tbl
COGinfo=${OESADIR}/db/cog/cog2003-2014.csv
COGLen=${OESADIR}/db/cog/COGconsensus_LengthInfo2.txt
RefGenome=${OESADIR}/db/cog/prot2003-2014.fa 

Genome=$(basename ${genome})
GenName=${Genome%.*}

[ -z "${OutDir}" ] && OutDir=${GenName}_curation

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

mkdir ${OutDir}
cp ${genome}  ./${OutDir}/
cp ${gff}  ./${OutDir}/
cd  ./${OutDir}

#------------------------------------------------1. IGS ORF Searech

# IGS ORF(>30nt) (output: IGS_Seq_${GenName}.fna, IGS_${GenName}.gff)

printf "\n## ------------------------ IGS ORF Search ------------------ ##\n\n"
IGSextract.bash  ${Gff}  ${Genome}  ${GenName}

# Extract BLAST hit from IGS  (output: IGS_BLASThit.fna, IGS_BLASThit.gff)

blasthitEXT.sh  -q  IGS_Seq_${GenName}.fna  -g ${Genome} -r ${RefGenome}  -t igs  -o IGS -e 1e-5  -n ${numCPU:=1} --code ${gCode:=11}

# ORFwalk for IGS hit  (output:  ${GenName}_ORFw.gff)

ORFwalk.sh       IGS_BLASThit.gff  ${Genome}  60  IGS ${gCode:=11}

# Integrate IGShit_ORFw seq and remove large overlap with exist sequence  (output: ${GenName}_2nd_merge.gff)

GFFmerge.bash   ${Gff}  IGS_ORFw.gff  ${GenName}_IGShit

#------------------------------------------------2. ORF length evalation

# Short ORF detection (Length evalation): 1st
printf "\n## ------------------------ ORF Length Evaluation ------------------ ##\n\n"

ORFLengthEval.sh   ${GenName}_IGShit_merge.gff   ${Genome}  ${RefGenome}  ${COGinfo} ${COGLen}   1st  ${numCPU:=1}  0.7  ${gCode:=11}

# ORF-walk for detected short ORFs

ORFwalk.sh       ShortEval_1st.gff  ${Genome}  60  Short1st  ${gCode:=11}

# Short ORF detection (Length evalation): 2nd

ORFLengthEval.sh  Short1st_ORFw.gff  ${Genome}  ${RefGenome}   ${COGinfo} ${COGLen}  2nd  ${numCPU:=1} 0.8 ${gCode:=11}

GFFmerge.bash    ShortEval_2nd.gff  ${GenName}_IGShit_merge.gff  ${GenName}_IGShit_Short2nd

#------------------------------------------------3. Pseudo-pseudo/Pseudo gene Searech

# PPgene  (output: ${GenName}_PPgene.gff, ${GenName}_PPgene.txt, ${GenName}_PPrecode.fna, ${GenName}_PPrecode.faa)

printf "\n## ------------------------ Pseudo-pseudo Gene Detection ------------------ ##\n\n"

PPgeneDetect.sh  ${GenName}_IGShit_Short2nd_merge.gff  ${Genome}  ${RefGenome}  ${GenName}  ${numCPU:=1} ${COGinfo} ${gCode:=11}

ORFLengthEval.sh   ${GenName}_PPgeneDetect.gff   ${Genome}  ${RefGenome}  ${COGinfo} ${COGLen}   PPgene  ${numCPU:=1}  0.8 ${gCode:=11}

# Merge gffs again

GFFmerge.bash    ${GenName}_IGShit_Short2nd_merge.gff  ShortEval_PPgene.gff  ${GenName}_IGShit_Short2nd_PPgene

# Domain completeness evaluation

printf  "\n## ------------------------ Domains Completeness Evaluation ------------------ ##\n\n"

DomCompEval.sh  ${GenName}_IGShit_Short2nd_PPgene_merge.gff ${Genome}  ${CDD}  ${CDDid}  ${GenName}  ${numCPU:=1}  0.8

# colum6 in gff file will be "<SC>|<Short/Normal>|<x.x(%)>|<XxDomXXX>" for short ORFs

# Merge gffs again

GFFmerge.bash   ${GenName}_IGShit_Short2nd_PPgene_merge.gff   ${GenName}_DomEval.gff  ${GenName}_IGShit_PPgene_ShortEval

# PseudoGene detect

printf  "\n## ------------------------ PseudoGene Detection ------------------ ##\n\n"

PseudoGeneDetect.sh  ${GenName}_IGShit_PPgene_ShortEval_merge.gff ${GenName}

#------------------------------------------------4. Final
## Write out to gff, faa, fna(RNA)

printf  "\n## ------------------------ Generating Result Files ------------------ ##\n\n"

Rscript $(which gff2faaFin.R)  ${GenName}_PseudoEval_merge.gff ${Genome}  ${GenName}  ${gCode:=11} 3
mv  addLocusTag.gff  ${GenName}_OESA.gff
Rscript $(which gffSummary.R)  ${GenName}_OESA.gff ${Genome} ${GenName}

# cp ${GenName}_OESA.gff ${GenName}.faa ${GenName}.ffn ${GenName}.frn  ${GenName}_Summary.txt  ..

cd ..

printf  "\n## ------------------------ Pipeline Finished! Thank You!! ------------------ ##\n\n"

## END OF FILE