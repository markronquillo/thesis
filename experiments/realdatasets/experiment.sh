#!/bin/bash
 
# experiments.sh
# Sets up EMS-GT and competitors PMS8 and qPMS9.
#  Generates (l,d) planted motif problem datasets.
#  Performs r experimental runs and formats results.

# @author Aia Sia
# @version 1.0 9/09/2015

r=1 # number of experimental runs
# javac -d bin src/*.java # compile java source files

# cd src/PMS8; make; mv Debug/PMS8  ../../bin/PMS8;  # uncomment these lines to compile PMS8 and qPMS9 from scratch.
# cd ../qPMS9; make; mv NoMpi/qpms9 ../../bin/qPMS9; # see Readme.md files in src/PMS8 and src/qPMS9 for dependency info.

rm -Rf results/
mkdir results

cd datasets
# print headers in EMS_GT result files
# start r runs
for((i=1; i <= r; i++)) do
	../bin/EMS_GT real_c-fos.txt 3 13 5 GTTCCCGTCAATC &>> ../results/real_c-fos-13,5-GTTCCCGTCAATC
	../bin/EMS_GT real_c-fos.txt 3 12 4 TACTCCAACCGC &>> ../results/real_c-fos-12,4-TACTCCAACCGC
	../bin/EMS_GT real_c-fos.txt 3 10 4 AGGACATCTG &>> ../results/real_c-fos-10,4-AGGACATCTG
	../bin/EMS_GT real_c-fos.txt 3 10 4 CACAGGATGT &>> ../results/real_c-fos-10,4-CACAGGATGT
	../bin/EMS_GT real_c-fos.txt 3 10 3 GAGTTGGCTG &>> ../results/real_c-fos-10,3-GAGTTGGCTG

	../bin/EMS_GT real_c-myc.txt 4 11 4 AGCAGAGGGCG &>> ../results/real_c-myc-11,4-AGCAGAGGGCG
	../bin/EMS_GT real_c-myc.txt 4 11 3 ATCTCCGCCCA &>> ../results/real_c-myc-11,3-ATCTCCGCCCA
	../bin/EMS_GT real_c-myc.txt 4 10 3 GGCGCGCAGT &>> ../results/real_c-myc-10,3-GGCGCGCAGT
	../bin/EMS_GT real_c-myc.txt 4 10 3 CAGCTGTTCC &>> ../results/real_c-myc-10,3-CAGCTGTTCC
	../bin/EMS_GT real_c-myc.txt 4 8 3  GGCGTGGG &>> ../results/real_c-myc-8,3-GGCGTGGG
	../bin/EMS_GT real_c-myc.txt 4 8 2 TTGCTGGG &>> ../results/real_c-myc-8,2-TTGCTGGG
	../bin/EMS_GT real_c-myc.txt 4 8 1 GTTTATTC &>> ../results/real_c-myc-8,1-GTTTATTC

	../bin/EMS_GT real_growth_hormone.txt 8 9 4 ATTATCCAT &>> ../results/real_growth_hormone-9,4-ATTATCCAT
	../bin/EMS_GT real_growth_hormone.txt 8 9 4 ATAAATGTA &>> ../results/real_growth_hormone-9,4-ATAAATGTA
	../bin/EMS_GT real_growth_hormone.txt 8 9 4 TCATGTTTT &>> ../results/real_growth_hormone-9,4-TCATGTTTT
	../bin/EMS_GT real_growth_hormone.txt 8 9 3 TATAAAAAG &>> ../results/real_growth_hormone-9,3-TATAAAAAG
	../bin/EMS_GT real_growth_hormone.txt 8 9 3 TTAGCACAA &>> ../results/real_growth_hormone-9,3-TTAGCACAA
	../bin/EMS_GT real_growth_hormone.txt 8 8 3 GGGAGGAG &>> ../results/real_growth_hormone-8,3-GGGAGGAG
	../bin/EMS_GT real_growth_hormone.txt 8 8 3 GTCAGTGG &>> ../results/real_growth_hormone-8,3-GTCAGTGG

	../bin/EMS_GT real_histoneH1.txt 2 10 3 CAATCACCAC &>> ../results/real_histoneH1-10,3-CAATCACCAC
	../bin/EMS_GT real_histoneH1.txt 2 10 3 AAACAAAAGT &>> ../results/real_histoneH1-10,3-AAACAAAAGT

	../bin/EMS_GT real_insulin.txt 4 10 4 AAGACTCTAA &>> ../results/real_insulin-10,4-AAGACTCTAA
	../bin/EMS_GT real_insulin.txt 4 10 3 GCCATCTGCC &>> ../results/real_insulin-10,3-GCCATCTGCC
	../bin/EMS_GT real_insulin.txt 4 8 2 CTATAAAG &>> ../results/real_insulin-8,2-CTATAAAG
	../bin/EMS_GT real_insulin.txt 4 8 2 GGGAAATG &>> ../results/real_insulin-8,2-GGGAAATG

	../bin/EMS_GT real_interleukin-3.txt 3 10 4 &>> ../results/real_interleukin-3-10,4
	../bin/EMS_GT real_interleukin-3.txt 3 9 3 &>> ../results/real_interleukin-3-9,3
	../bin/EMS_GT real_interleukin-3.txt 3 8 2 &>> ../results/real_interleukin-3-8,2

	../bin/EMS_GT real_metallothionein.txt 13 11 5 TAACTGATAAA &>> ../results/real_metallothionein-11,5-TAACTGATAAA
	../bin/EMS_GT real_metallothionein.txt 13 10 3 TTTGCACACG &>> ../results/real_metallothionein-10,3-TTTGCACACG
	../bin/EMS_GT real_metallothionein.txt 13 9 3 CATGCGCAG &>> ../results/real_metallothionein-9,3-CATGCGCAG
	../bin/EMS_GT real_metallothionein.txt 13 9 3 TACACTCAG &>> ../results/real_metallothionein-9,3-TACACTCAG
	../bin/EMS_GT real_metallothionein.txt 13 9 3 CAGGCACCT &>> ../results/real_metallothionein-9,3-CAGGCACCT
	../bin/EMS_GT real_metallothionein.txt 13 9 3 GTACATTGT &>> ../results/real_metallothionein-9,3-GTACATTGT
	../bin/EMS_GT real_metallothionein.txt 13 8 3 GCTATAAA &>> ../results/real_metallothionein-8,3-GCTATAAA

	../bin/EMS_GT yeast_hse_hstf.txt 3 8 3 &>> ../results/yeast_hse_hstf-8,3
	../bin/EMS_GT yeast_hse_hstf.txt 3 8 2 &>> ../results/yeast_hse_hstf-8,2

	../bin/EMS_GT yeast_mcb.txt 3 6 2 &>> ../results/yeast_mcb-6,2

	../bin/EMS_GT yeast_pdr3.txt 4 8 2 &>> ../results/yeast_pdr3-8,2

	../bin/EMS_GT yeast_pho4.txt 2 6 2 &>> ../results/yeast_pho4-6,2
done
cd ..