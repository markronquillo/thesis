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
	../bin/EMS_GT_NEW real_c-fos.txt &>> ../results/real_c-fos
	../bin/EMS_GT_NEW real_c-myc.txt &>> ../results/real_c-myc
	../bin/EMS_GT_NEW real_growth_hormone.txt &>> ../results/real_growth_hormone
	../bin/EMS_GT_NEW real_histoneH1.txt &>> ../results/real_histoneH1
	../bin/EMS_GT_NEW real_insulin.txt &>> ../results/real_insulin
	../bin/EMS_GT_NEW real_interleukin-3.txt &>> ../results/real_interleukin-3
	../bin/EMS_GT_NEW real_metallothionein.txt &>> ../results/real_metallothionein

	../bin/EMS_GT_NEW yeast_hse_hstf.txt &>> ../results/yeast_hse_hstf
	../bin/EMS_GT_NEW yeast_mcb.txt &>> ../results/yeast_mcb
	../bin/EMS_GT_NEW yeast_pdr3.txt &>> ../results/yeast_pdr3
	../bin/EMS_GT_NEW yeast_pho4.txt &>> ../results/yeast_pho4
done
cd ..