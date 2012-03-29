#!/bin/bash
#$ -N Demultiplex-with-CASAVA-2012-03-29.1
#$ -P nne-790-ab
#$ -l h_rt=48:00:00
#$ -pe default 8
#$ -cwd

sequenceWorld=/rap/nne-790-ab/Instruments/Illumina_MiSeq_Tombstone/
run=120321_M00212_0002_AMS2002881-00300
baseMask=Y*,I*,I*,Y*
NSLOTS=8

##############################

# charger gcc, les modules perl et CASAVA
source /rap/nne-790-ab/software/CASAVA_v1.8.2/module-load.sh

#  --mismatches 1,1 may increase the yield...

outputDirectory=$sequenceWorld/$run/Demultiplex-with-CASAVA-2012-03-29.1

configureBclToFastq.pl --input-dir $sequenceWorld/$run/Data/Intensities/BaseCalls \
 --output-dir  $outputDirectory \
 --sample-sheet  $sequenceWorld/$run/SampleSheet-selected.csv \
 --use-bases-mask $baseMask

cd $outputDirectory

make -j $NSLOTS

