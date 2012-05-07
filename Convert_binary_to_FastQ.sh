#!/bin/bash
#$ -N convert-MiSeq
#$ -P nne-790-ab
#$ -l h_rt=48:00:00
#$ -pe default 8
#$ -cwd

NSLOTS=8

outputDirectory=Convert_binary_to_FastQ

source /rap/nne-790-ab/software/CASAVA_v1.8.2/module-load.sh

configureBclToFastq.pl \
--input-dir Data/Intensities/BaseCalls \
--output-dir  $outputDirectory \
--use-bases-mask Y*,Y*,Y*,Y*

cd $outputDirectory

make -j $NSLOTS

