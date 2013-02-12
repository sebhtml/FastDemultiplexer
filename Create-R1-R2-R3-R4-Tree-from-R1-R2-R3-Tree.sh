#!/bin/bash
# Sébastien Boisvert
# 2013-02-12

if test $# != 1
then
	echo "This tool takes a tree with R1, R2, and R3 files (sequences in R1+R3, 1 barcodes in R2)"
	echo "and create a tree with R1, R2, R3 and R4 files (sequences in R1+R4, 2 barcodes in R2+R3"
	echo "Known bug: the project must be a relative path (patches are welcomed!)"
	echo "Usage"
	echo "$0 Project_C1MMMACXX"
	exit
fi

project=$1

projectName=$(basename $project)

mkdir $projectName

for sample in $(ls $project|grep Sample_)
do
	mkdir $projectName/$sample

	for sequenceFile in $(ls $project/$sample/*|grep _R1_)
	do
		ln -s ../../$sequenceFile $projectName/$sample/$(basename $sequenceFile)
	done

	for sequenceFile in $(ls $project/$sample/*|grep _R2_)
	do
		ln -s ../../$sequenceFile $projectName/$sample/$(basename $sequenceFile)
		ln -s ../../$sequenceFile $projectName/$sample/$(echo $(basename $sequenceFile)|sed 's/_R2_/_R3_/g')
	done

	for sequenceFile in $(ls $project/$sample/*|grep _R3_)

	do
		ln -s ../../$sequenceFile $projectName/$sample/$(echo $(basename $sequenceFile)|sed 's/_R3_/_R4_/g')
	done

done
