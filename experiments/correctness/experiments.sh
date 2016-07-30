#!/bin/bash
 
# experiments.sh
# Sets up EMS-GT and competitors PMS8 and qPMS9.
#  Generates (l,d) planted motif problem datasets.
#  Performs r experimental runs and formats results.

# @author Aia Sia
# @version 1.0 9/09/2015

r=20 # number of experimental runs
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
	# print headers in EMS_GT result files
	# start r runs
	for((i=1; i <= r; i++)) do
		# generate a unique (l,d) dataset for this run
		java -cp ../ DatasetGenerator  $l $d $i
		# convert this dataset to FASTA
		java -cp ../ DatasetConverter  $l,$d,$i
		# test all programs on this dataset
		# java -cp ../bin EMS_GT $l,$d,$i 	>> ../results/EMS_GT-$l,$d
		#java -cp ../bin EMS_GT_32 $l,$d,$i 	>> ../results/EMS_GT_32-$l,$d
		../EMS_GT_CLEAN $l,$d,$i 10	>> ../results/EMS_GT_CLEAN-$l,$d


		# candidate elimination only
		../EMS_GT_NEW $l,$d,$i 	>> ../results/EMS_GT_NEW-$l,$d
		
	done
done
cd ..
