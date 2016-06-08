#include <string>
#include <vector>

using namespace std;

class Config {
public:
	Config(const DataSetParams *ds) {
		if (ds->numberOfSequences == 0) {
			cout << "DataSetParameters is not set." << endl;
			return;
		}

    	numberOfBlockRows = 1 << ( (2*(ds->lengthOfMotif)) - blockDegree - TABLE_WIDTH_LOG_2 );

    	// compute the preliminary values needed by the program
		mask =  (((long) 1) << ((2*ds->lengthOfMotif)-2)) - 1;  

		prefixMask = (((long) 1) << 2*(ds->lengthOfMotif-blockDegree-1)) - 1;

		suffixMask = (((long) 1) << 2*(blockDegree-1)) - 1;

		prefixShift = (2*blockDegree - TABLE_WIDTH_LOG_2);

    	numberOfPossibleLmersInSequence = ds->lengthOfSequence - ds->lengthOfMotif + 1;

		numberOfCandidateMotifsRows = (((long)1 << (2*ds->lengthOfMotif)) >> 5);
	}

	/**
	 *  This holds the value for the following
	 *  - length of suffix
	 *  - 4 ^ blockDegree = number of lmers per block
	 */
	const int blockDegree = 5;

	/**
	 *  The collection of candidate motifs will only 
	 *  be done up to this value
	 */
	const int tPrime = 10;

	/**
	 *	Usage of boolean flags in generation of neighborhood only starts
	 *	from this value
	 */
	const int tPrime_2 = 15;


	/**
	 *  Holds the MAX integer bit value
	 */
	const int TABLE_WIDTH_LOG_2 = 5; // 2^5 = 32

	/**
	 *  Holds the computed table width based on the TABLE_WIDTH_LOG_2 value
	 */
	const int TABLE_WIDTH = 1 << TABLE_WIDTH_LOG_2; 

	/**
	 *  Holds the total number of lmers in a block
	 *  Computed as: 4 ^ blockDegree
	 */	
	const int numberOfLmersInBlock = (1 << (blockDegree*2));

	/**
	 *  Holds the total number of rows in a block
	 *  Computed as: 4 ^ blockDegree / 32
	 */
	const int numberOfRowsInBlock = numberOfLmersInBlock / 32;
	

	// --------------------------------------------------------------------


	/**
	 *  Holds the total number of blocks
	 */
	int numberOfBlockRows;

	/**
	 *  Holds the total number of rows in a block
	 *  Computed as: ((2 * lengthOfMotif) - blockDegree)
	 */
	int prefixShift;

	/**
	 *	Removes exceeding bits
	 */
	long mask, suffixMask, prefixMask;

	/**
	 *  This holds the number of possible l-mers 
	 *  given a string sequence and the length of the motif.
	 *
	 *  numberOfSequences - length of motif + 1
	 */
	int numberOfPossibleLmersInSequence;

	/**
	 *  This holds the height value of the candidateMotifs array
	 *
	 *
	 *  numberOfCandidateMotifsRows = (int) (( (long)1 << (2*l) ) >>> 5); 
	 */
	long numberOfCandidateMotifsRows;


};