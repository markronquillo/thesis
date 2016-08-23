[1]
Our study focuses on solving the Motif Finding Problem. 
EMS-GT2 is an improvement of the previous EMS-GT, an exact enumerative algorithm that solves the motif finding problem.

[3]
In simple definition, the motif finding problem is a problem that is about looking for common substrings, which is called a motif, in a set of DNA sequences considering mutations since biological DNA sequences are subject to mutations. This property of the problem made it NP-Complete.

[4]
One variant of the motif finding problem is the (l, d)-planted motif search problem. And by definition, it is the a problem that accepts `n` number of DNA sequences with `m` length each. And given an l-mer to be planted in each sequence only once and an integer value `d` that is the number of mutations allowed.

[5] For a demonstration of the problem.
Given the [configuration]. 
In this demonstration the planted l-mer is `tactg` and substrings in each sequence in yellow background are the positions where the motif is planted. We can see that all these l-mers are at most d-distance with the motif `tactg`.

[14] Neighborhood of a Sequences
Basically it is the union of all the d-neighborhood of all l-mers in that sequence.

[17] EMS-GT
There are two common ways of solving the (l, d)-planted motif search problem that are used widely by exact algorithms ... First, it checks every possible combination of positions across different sequences and tests if it is the correct positions where the motif is planted. Second is to exhaustively search all 4^l possible l-mers if it is the motif. 

We introduced a novel algorithm, an exact enumerative algorithm that solves the (l, d)-planted motif search problem. EMS-GT exhaustively search all 4^l possible l-mers. Exact algorithms return the correct answer all the time at the expense of having an exponential runtime since the problem is NP-Complete. 

[18]
EMS-GT has two phases, the Generate phase and Test phase. The Generate phase quickly eliminates the l-mers in the 4^l by generating the d-neighborhood of a sequence and 

The Test phase, given the candidate motifs, each will be tested if it has at least one d-neighbor in each of the rest of the sequences.

[EMS-GT DEMO]

[33] Implementation
There are some notable techniques in our implementation of the algorithm.

First,
[34] Integer Mapping of l-mers
We convert the l-mers in its corresponding integer values by mapping each character with its corresponding two-bit representation
a=00, c=01, g=10, t=11
For example:

Next,
[35] Bit-based set representation and l-mer enumeration
The algorithm needs a way to track or maintain the set of d-neighborhood of an l-mer or a sequence. To do this, we maintain an array of bits that flags if an l-mer is a member of the set or not. We used the integer value of an l-mer as the index in array.

Furthermore
[36] Bit-array compression
We improved the space requirement of the array by compressing it using a 32-bit integer to represent flags. So instead of just one flag per row (or index), we can have 32 flags each.


[37] XOR-based hamming distance computations
For the hamming distance computations, since we are using integer representation of the l-mers we can use the boolean operator XOR to efficiently locate the mismatch positions.

For example:
The aacgt is mapped to this binary and the tacgc is mapped to this binary
Since the XOR operator results 0 if the bits are the same and 1 if otherwise, we can count each pair of bits in the resulting binary to get the hamming distance,

[38] Recursive generation of neighborhood
As demonstrated earlier we need to generate a set of d-neighbor of an l-mer, to do this, we recursively alter each character in each position. each node in this recursive tree represents an l-mer in the neighborhood.


[39] Block-based optimization for neighborhood generation.
The recursive generation of neighborhood sets the array bits one bit a time. We improved this by setting the bits by a pre-computed blocks of bits instead of just one bit at a time.

[44-45] Demonstration
For example, lets say we will generate the d-neighborhood of this 15-mer, since initially we know the size of our pre-computed block patterns, we know the value k, since it is 4^k. The k value represents the suffix length in the l-mer. Which means we only to recusively generate upto the l-k characters in the l-mer, which is 10 in this example. And for each prefix we apply the corresponding block mask based on the remaining number of allowed mutations.

[47] Additional Speed ups
Introduction:
This study introduces two new speedup techniques that is mainly used in the Test phase.

We introduced a faster way to eliminate a candidate motif through block processing. 
First, we list some observations and properties which helped in constructing the technique.
Read 1. 2. and 3.

To do this, we first partitioned the array into 4^k size blocks. We can derive two properties that is true in each block. 
1. Each l-mer in a block share the same prefix of length l-k.
2. Second the distance between any l-mer in a block is at most `k`

Illustrate:
For an illustration of the theorem we used.
lets say there's an l-mer x and y within a block, what we know is the hamming distance between them is at most `k`.
Lets say there is another l-mer `z` located in a distant block, and the distance between x and z
Then we know that the distance between y and z cannot be less than or equal to d.

We used this theorem in Test phase by filtering l-mers in a sequences.

The improved test phase is now processed by blocks,
The first step is to test candidate motifs in a block until it is eliminated.
Lets say gagct is eliminated at Sn'' sequence. then we can eliminate l-mers in Sn'' where the distance between gagct and that l-mer is greater than d + k), using the theorem explained earlier.
By doing this we can quickly eliminate a candidate motif if it is indeed not a candidate motif.


[65] Pre-computation
The second speedup technique that we implemented in our implementation is the pre-computation of mismatch values for the hamming distance computation.

To recall, the computation of hamming distance uses the binary representation of l-mers and uses XOR operator to get the mismatch position. XOR results another integer that contains the mismatch positions.
For every hamming distance computation in the original EMS-GT, it counts every pairs of bits and keep track of the non-zero ones since it represents mismatch position.

The precomputed array looks like this.
But pre-computation of all mismatch values for all integer values will cause overhead computation problems. What we did is, we only pre-compute upto a certain number of bits only and use it repeatedly

To illustrate:
Lets say we have 16-bit of pre-computed mismatch values, and given this binary XOR result. We first get the mismatch values of the first 16-bit. the total hamming distance in this case is 13

So instead of checking the 16 pairs of bits, we only need to lookup twice.

[75] Method
For the dataset, we used 20 different DNA sequence where each sequence 600 characters long and we evaluated the algorithms using (l, d) challenge instances

To be exact, we did the evaluation of EMS-GT, qpms9 and the improved EMS-GT2 using the .... challenge instances over 20 iterations and we get the average runtime of each.

[80] To conclude

EMS-GT2 is efficient in solving challenge instances where l < 18
EMS-GT2 proved its competitiveness over qPMS9, by beating the algorithm on challenge instances mentioned earlier.
Lastly, in practice, studies have proved that motifs are normally 10bps long so eventhough EMS-GT2 has a memory constraint for instances where l > 17, it is still significant.





