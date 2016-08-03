#!/bin/bash
 
# experiments.sh
# Sets up EMS-GT and competitors PMS8 and qPMS9.
#  Generates (l,d) planted motif problem datasets.
#  Performs r experimental runs and formats results.

# @author Aia Sia
# @version 1.0 9/09/2015

r=50 # number of experimental runs
# javac -d bin src/*.java # compile java source files

# cd src/PMS8; make; mv Debug/PMS8  ../../bin/PMS8;  # uncomment these lines to compile PMS8 and qPMS9 from scratch.
# cd ../qPMS9; make; mv NoMpi/qpms9 ../../bin/qPMS9; # see Readme.md files in src/PMS8 and src/qPMS9 for dependency info.

d=1
rm -Rf datasets/
rm -Rf results/
mkdir results
mkdir datasets
mkdir datasets/FASTA

cd datasets
for l in 9 11 13 15 17 
do
	((d++))
	for((i=1; i <= r; i++)) do
		java -cp ../ DatasetGenerator  $l $d $i
		java -cp ../ DatasetConverter  $l,$d,$i

		../EMS_GT_NEW $l,$d,$i 	>> ../results/EMS_GT_NEW-$l,$d
	done
done
cd ..
