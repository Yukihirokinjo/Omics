#!/bin/bash
#
# ProdiLCORFrm.bash  [ prodi.gff ] [ prodi.faa ] [temp_genome.fasta ] [ Ref_prot.faa ] [ output_prefix ]
#
# 1.  with no conservation -> annot artifact
# 2.  with high conservation -> specific hypo


ProdiGFF=""
ProdiFAA=""
Genome=""
RefProteome=""
Name=""
Eval=""


usage_exit() {
	echo "Usage: ProdiLCORFrm.bash  [ -gff prodi_gff] [-faa prodi_faa] [ -g template_genome ] [ -ref reference proteome ]  [-o output_Prefix]  [-e E-value ] [-n num._CPU]" 1>&2
	exit 1
}


while [ "$#" -gt 0 ]
do
	case "$1" in
		'-h' | '--help' )
			usage_exit
			;;
		'-gff')
			if  [  -e "$2"  ]; then
				ProdiGFF="$2"
				shift 2
			else
				echo "[Error] The input GFF file is not found" 1>&2
				exit 1
			fi
			;;
		'-faa')
			if  [  -e "$2"  ]; then
				ProdiFAA="$2"
				shift 2
			else
				echo "[Error] The input Proteome(.faa) file is not found" 1>&2
				exit 1
			fi
			;;
		'-ref')
			if  [  -e "$2"  ]; then
				RefProteome="$2"
				shift 2
			else
				echo "[Error] The input Reference(.faa) file is not found" 1>&2
				exit 1
			fi
			;;
		'-g')
			if  [  -e "$2"  ]; then
				Genome="$2"
				shift 2
			else
				echo "[Error] The input Genome(.fasta/fna) file is not found" 1>&2
				exit 1
			fi
			;;
		'-o')
			if  [ -n  "$2" ]; then
				Name="$2"
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
if [ -z "$ProdiGFF" ]; then
	echo "Input GFF file is not specified" 1>&2
	usage_exit
fi
if [ -z "$Genome" ]; then
	echo "Input Genome file is not specified" 1>&2
	usage_exit
fi
if [ -z "$RefProteome" ]; then
	echo "Input Reference file is not specified" 1>&2
	usage_exit
fi
if [ -z "$ProdiFAA" ]; then
	echo "Input Proteome(.faa) is not specified" 1>&2
	usage_exit
fi
if [ -z "$Name" ]; then
	Name="LC"
fi

mkdir  LCcheck
cp ${ProdiGFF} ${ProdiFAA} ${Genome} ${RefProteome}  LCcheck
cd LCcheck

grep -v ^"#" ${ProdiGFF} | \
awk  'BEGIN{FS=";"}{ Conf=(substr($7, (index($7,"=") + 1) ))}{ if(int(Conf)< 70) print $0 }' > ProdigiLC.gff
diff ${ProdiGFF}  ProdigiLC.gff  | grep ^"<" | sed "s/< //g" > tmp_ProdigiLCrm.gff

  Rscript $(which gff2faa.R)  ProdigiLC.gff  ${Genome}  LowConf

cat gff2.faa > ${Name}_prodi_LowConf.faa
rm  gff2.faa

  blasthitEXT.sh  -q ${Name}_prodi_LowConf.faa  -g ${Genome} -r ${RefProteome}  -t prot -o LowConf -e ${Eval:=1e-10}

cat tmp_ProdigiLCrm.gff  ./blaEXT/LowConf_BLASThit.gff  | sort -k 1,1 -k 4n,4  > ${Name}_ProdiHC.gff
cp ${Name}_ProdiHC.gff ..

cd ..

## END OF FILE
