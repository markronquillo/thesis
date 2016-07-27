## Introduction
Motif finding and pattern search in biological sequences has many applications thus numerous research has been conducted in this area. Finding motif in a set of DNA sequences has many applications such as finding transcription site, gene regularization ....

This study mainly focuses on Planted Motif Search (PMS) problem, also known as (l, d) motif problem. The problem simply states that "Given a set of DNA sequences can we find an unknown motif of length l that appears at different positions in each of the sequences" (Pevsner). PMS problem takes an input n sequences of length m each and integers l and d that holds the value for the length of the motif planted and the number of allowed mutations.

----
## RRL 

A lot of algorithms have been proposed in this study. There are exact and approximate algorithms. An algorithm that always find all the motifs in the given input is called exact algorithm. Since PMS problem is known as NP-complete, all known exact algorithms have exponential worst case runtime.
(elaborate on algorithms on both) The algorithm that is improved in this study falls on the exact algorithm category. (read about pattern-based and the one that exhausts the 4**l lmers)

This study will propose speedup-techniques that will improve the EMS-GT with blockmasking algorithm.

Previous studies about the algorithm have proven that ... (qpms9 comparison)

-----
## Methodology

<!-- Introduction, about 2 sentences -->

### Datasets
(mention that you have updated the datasetgenerator to generate also in FASTA format)

### Implementation
(mention the port to c++ and its importance, majority of the PMS algorithm was implemented using c++)


### Evaluation
(comparison between EMS-GT with blockmasking) and the algorithm will be compared with qPMS9.

### Speedup Techniques

#### Hamming distance Pre-Computation 
The Test phase of the algorithm checks if a candidate motif, from the filtered search space, is a motif by checking if it exists in the remaining string sequences. The checking part of the algorithm is done by getting the Hamming Distance of the motif versus every lmer in remaining sequences. If a candidate motif is not present in a sequence it is already eliminated as a candidate motif. The EMS-GT algorithm uses exclusive OR (XOR) bitwise operation between the mappings of two lmers that will output nonzero pair of bits for every mismatch.
(explain how the speedup technique works) (thus improving the Test Phase)


#### Heuristic value for t-prime
The parameter t-prime holds the number of string sequences that will be used in the Generation of common neighborhood phase of the algorithm. In both EMS-GT and EMS-GT with blockmasking used the value 10 for the t-prime variable. Since in this study we improved the running time of the Test phase, it is ideal to set the value of the t-prime depending on the number of candidates motif left in the Generation phase.

---- 


## Results

<!-- (17,6),  -->