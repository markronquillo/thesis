#!/bin/bash
 
# experiments.sh
# Sets up EMS-GT and competitors PMS8 and qPMS9.
#  Generates (l,d) planted motif problem datasets.
#  Performs r experimental runs and formats results.

# @author Aia Sia
# @version 1.0 9/09/2015

r=20 # number of experimental runs
# javac -d bin src/*.java # compile java source files

d=1
rm -Rf datasets/
rm -Rf results/
mkdir results
mkdir datasets

cd datasets
# (l,d): (9,2) (11,3) (13,4) (15,5) 
for l in 9 11 13 15 17 
do
	((d++))
	# print headers in EMS_GT result files
	# start r runs
	for((i=1; i <= r; i++)) do
		# generate a unique (l,d) dataset for this run
		java -cp ../bin DatasetGenerator  $l $d $i

		../bin/EMS_GT_CLEAN $l,$d,$i 	&>> ../results/EMS_GT-$l,$d
	done
done
cd ..