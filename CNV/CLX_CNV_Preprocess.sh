#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o log/
#$ -j y

## dir lib
dir=soft/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles
apt=soft/apt-1.16.1-x86_64-intel-linux/bin/apt-copynumber-workflow
files=soft/files


## ------------------------------------------------------------------------------------

## Genotipos
soft/apt-1.16.1-x86_64-intel-linux/bin/apt-probeset-genotype \
	-o Genotype \
	-c $dir/GenomeWideSNP_6.cdf \
	--set-gender-method cn-probe-chrXY-ratio \
	--chrX-probes $dir/GenomeWideSNP_6.chrXprobes \
	--chrY-probes $dir/GenomeWideSNP_6.chrYprobes \
	--special-snps $dir/GenomeWideSNP_6.specialSNPs \
	--read-models-birdseed $dir/GenomeWideSNP_6.birdseed-v2.models \
	-a birdseed-v2 \
	--cel-files $files/listfile

$apt \
--adapter-type-normalization true \
--reference-output $files/GenomeWideSNP_6.hapmap270.na34.r1.a5.ref \
--set-analysis-name clx \
--cdf-file $dir/GenomeWideSNP_6.cdf \
--chrX-probes $dir/GenomeWideSNP_6.chrXprobes \
--chrY-probes $dir/GenomeWideSNP_6.chrYprobes \
--special-snps $dir/GenomeWideSNP_6.specialSNPs \
--annotation-file $files/GenomeWideSNP_6.na34.annot.db \
--delete-files true \
--o Results \
--text-output true \
--delete-files false \
--cel-files $files/listfile

## ------------------------------------------------------------------------------------

ls Results | grep ".clx.CN5.CNCHP.txt" > list_file

while read file
do
   name="$(echo $file | cut -d'_' -f 3)"
   status="$(echo $file | cut -d'_' -f 4 | cut -d'.' -f 1)"
   sed '/#%/d' Results/$file > res_ind/$status/$name
done < list_file

rm -f list_file

## ------------------------------------------------------------------------------------



