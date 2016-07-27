#!/bin/bash
t=20
m=600
# Generate Data
rm -Rf dataset/
mkdir dataset
mkdir dataset/FASTA
# Run the algorithms
d=1
for l in 9 11 13 15 17 
do
	((d++))
	# Generate datasets for the current (l, d) values
	for((i=1; i<=t; i++)) do 
		filename="data-$l-$d-$i"
		java DataSetGenerator $l $d $i
	done
	# clear and create folder for the current (l, d)
	rm -Rf result-$l-$d/
	mkdir result-$l-$d
	mkdir result-$l-$d/EMS_GT
	mkdir result-$l-$d/qPMS9
	mkdir result-$l-$d/PMS8
	for((i=1; i<=t; i++)) do 
		filename="$l,$d,$i"
		./EMS_GT dataset/$filename >> result-$l-$d/EMS_GT/$filename
		./qpms9 -l $l -d $d dataset/FASTA/$filename &>> result-$l-$d/qPMS9/$filename
		./PMS8 -l $l -d $d dataset/FASTA/$filename &>> result-$l-$d/PMS8/$filename
	done
done





