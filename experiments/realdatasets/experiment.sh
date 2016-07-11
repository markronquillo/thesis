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

rm -Rf results/
mkdir results

cd datasets
# print headers in EMS_GT result files
# start r runs
for((i=1; i <= r; i++)) do
	../bin/EMS_GT real_c-fos.txt 3 13 5 &>> ../results/real_c-fos-13,5
	../bin/EMS_GT real_c-fos.txt 3 12 4 &>> ../results/real_c-fos-12,4
	../bin/EMS_GT real_c-fos.txt 3 10 4 &>> ../results/real_c-fos-10,4
	../bin/EMS_GT real_c-fos.txt 3 10 3 &>> ../results/real_c-fos-10,3

	../bin/EMS_GT real_c-myc.txt 4 11 4 &>> ../results/real_c-myc-11,4
	../bin/EMS_GT real_c-myc.txt 4 11 3 &>> ../results/real_c-myc-11,3
	../bin/EMS_GT real_c-myc.txt 4 10 3 &>> ../results/real_c-myc-10,3
	../bin/EMS_GT real_c-myc.txt 4 8 3 &>> ../results/real_c-myc-8,3
	../bin/EMS_GT real_c-myc.txt 4 8 2 &>> ../results/real_c-myc-8,2
	../bin/EMS_GT real_c-myc.txt 4 8 1 &>> ../results/real_c-myc-8,1

	../bin/EMS_GT real_growth_hormone.txt 8 9 4 &>> ../results/real_growth_hormone-9,4
	../bin/EMS_GT real_growth_hormone.txt 8 9 3 &>> ../results/real_growth_hormone-9,3
	../bin/EMS_GT real_growth_hormone.txt 8 8 3 &>> ../results/real_growth_hormone-8,3

	../bin/EMS_GT real_histoneH1.txt 2 10 3 &>> ../results/real_histoneH1-10,3

	../bin/EMS_GT real_insulin.txt 4 10 4 &>> ../results/real_insulin-10,4
	../bin/EMS_GT real_insulin.txt 4 10 3 &>> ../results/real_insulin-10,3
	../bin/EMS_GT real_insulin.txt 4 8 2 &>> ../results/real_insulin-8,2

	../bin/EMS_GT real_interleukin-3.txt 3 10 4 &>> ../results/real_interleukin-3-10,4
	../bin/EMS_GT real_interleukin-3.txt 3 9 3 &>> ../results/real_interleukin-3-9,3
	../bin/EMS_GT real_interleukin-3.txt 3 8 2 &>> ../results/real_interleukin-3-8,2

	../bin/EMS_GT real_metallothionein.txt 13 11 5 &>> ../results/real_metallothionein-11,5
	../bin/EMS_GT real_metallothionein.txt 13 10 3 &>> ../results/real_metallothionein-10,3
	../bin/EMS_GT real_metallothionein.txt 13 9 3 &>> ../results/real_metallothionein-9,3
	../bin/EMS_GT real_metallothionein.txt 13 8 3 &>> ../results/real_metallothionein-8,3

	../bin/EMS_GT yeast_hse_hstf.txt 3 8 3 &>> ../results/yeast_hse_hstf-8,3
	../bin/EMS_GT yeast_hse_hstf.txt 3 8 2 &>> ../results/yeast_hse_hstf-8,2

	../bin/EMS_GT yeast_mcb.txt 3 6 2 &>> ../results/yeast_mcb-6,2

	../bin/EMS_GT yeast_pdr3.txt 4 8 2 &>> ../results/yeast_pdr3-8,2

	../bin/EMS_GT yeast_pho4.txt 2 6 2 &>> ../results/yeast_pho4-6,2
done
cd ..