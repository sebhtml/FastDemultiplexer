#!/bin/bash
#$ -N convert-120321_M00212_0002_AMS2002881-00300
#$ -P nne-790-ab
#$ -l h_rt=48:00:00
#$ -pe default 8
#$ -cwd

sequenceWorld=/rap/nne-790-ab/Instruments/Illumina_MiSeq_Tombstone/
run=120321_M00212_0002_AMS2002881-00300
NSLOTS=8

outputDirectory=$sequenceWorld/$run/Sequences/Convert_binary_to_FastQ

source /rap/nne-790-ab/software/CASAVA_v1.8.2/module-load.sh

configureBclToFastq.pl \
--input-dir $sequenceWorld/$run/Data/Intensities/BaseCalls \
--output-dir  $sequenceWorld/$run/FastQ-Sequences/no-demul-Unaligned \
--use-bases-mask Y*,Y*,Y*,Y*

mkdir -p $outputDirectory
cd $outputDirectory

make -j $NSLOTS

