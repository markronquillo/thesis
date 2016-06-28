#!/bin/bash
 
# experiments.sh
# Sets up EMS-GT and competitors PMS8 and qPMS9.
#  Generates (l,d) planted motif problem datasets.
#  Performs r experimental runs and formats results.

# @author Aia Sia
# @version 1.0 9/09/2015

r=5 # number of experimental runs
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
		java -cp ../bin DatasetGenerator  $l $d $i
		# convert this dataset to FASTA
		java -cp ../bin DatasetConverter  $l,$d,$i
		# test all programs on this dataset
		# java -cp ../bin EMS_GT $l,$d,$i 	>> ../results/EMS_GT-$l,$d
		#java -cp ../bin EMS_GT_32 $l,$d,$i 	>> ../results/EMS_GT_32-$l,$d
		for((n=5; n <= 10; n++)) do
			../bin/EMS_GT_NEW $l,$d,$i 	&>> ../results/EMS_GT_NEW-$l,$d
		done
		# java -cp ../bin EMS_GT_64 $l,$d,$i 	>> ../results/EMS_GT_64-$l,$d
		# ../bin/PMS8 FASTA/$l,$d,$i $l $d 	&>> ../results/rawPMS8-$l,$d
		# ../bin/qPMS9 -l $l -d $d FASTA/$l,$d,$i 	&>> ../results/rawqPMS9-$l,$d
	done
	# print headers in PMS result files
	# echo "l,d,run,time(s)" > ../results/PMS8-$l,$d
	# echo "l,d,run,time(s)" > ../results/qPMS9-$l,$d
	# interpret rawPMS8 & rawqPMS9 files; note that outputReaders
	#  are specific to this directory structure & file naming convention
	# java -cp ../bin PMS8outputReader  ../results/rawPMS8-$l,$d  >> ../results/PMS8-$l,$d
	# java -cp ../bin qPMS9outputReader ../results/rawqPMS9-$l,$d >> ../results/qPMS9-$l,$d
done
cd ..