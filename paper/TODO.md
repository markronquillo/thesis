## TODO:

. Appendix
	. conference paper
	. Codes

-------------------------------------------------
## Targets: Task - Deadline


--------------------------------------------------
## Notes

## Manuscript Sections
. Acknowledgments
. Abstract
. ToC
. List of tables
. Chapter 1
	. Introduction
		1.1 Context of the Study
		1.2 Objectives of the Study
		1.3 Research Questions
		1.4 Significance of the Study
		1.5 Scope and Limitations
. Chapter 2: RRL
	. Introduction
	. Approximate algorithms
	. Exact algorithms
. Chapter 3: Methodology
	. Introduction
	. Datasets
	. ?
	. Parameter fine tuning
	. Evaluation
. Chapter 4: Results and Analysis
	. Introduction
	. Pruning 
	. Block flags
	. Fast candidate elimination
	. Pre-computation of pair bits for Hamming distance computation
	. [INCLUDES]
		. pseudocode
		. tables
		. figures
. Chapter 5: Conclusion

## Ideas
* Count how many repetitions per blockstarts are there
	- create an array of size (# of blockstarts possible)
	- initialize value with zero
	- increment the value if block start is encounterettd
* Other way to generate a candidate motif?
	- new data structure? or adjust the current?
	- new way to populate the search space.
	** use blockmask twice 7 5 5  