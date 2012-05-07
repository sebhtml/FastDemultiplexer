#!/bin/bash
#$ -N Demultiplex-123
#$ -P nne-790-ab
#$ -l h_rt=48:00:00
#$ -pe default 8
#$ -cwd

PATH=.:$PATH

for i in $(seq 1 8)
do
	FastDemultiplexer.py SampleSheet.csv Project/Sample_lane$i Demultiplexed > l$i.log &
done

wait
