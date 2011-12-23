#!/bin/bash

path=$1

projectName=$(basename $path)

mkdir $projectName

for sample in $(ls $path)
do
	mkdir $projectName/$sample
	cd $projectName/$sample
	for i in $(find ../../$path/$sample|grep fastq)
	do
		ln -s $i
	done
	cd ../../
done

