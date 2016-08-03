#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <bitset>
#include <climits>


using namespace std;
using namespace std::chrono;

// #if defined(INT_MAX)
// #else 
//     #define INT_MAX 32767
// #endif

class EMS_GT {
public:
	EMS_GT(const DataSetParams *ds, const Config *c): 
		ds(ds), config(c)
	{
    	currentNeighborhood = new int[config->numberOfCandidateMotifsRows];
	}

	void start() {
    	high_resolution_clock::time_point t1 = high_resolution_clock::now();
    	// steady_clock::time_point begin = std::chrono::steady_clock::now();

		generateMismatchesCount();

		generateBlockMasks();

		collectCandidateMotifs();


    	high_resolution_clock::time_point ts = high_resolution_clock::now();
		transformLmerSequences();
    	high_resolution_clock::time_point te = high_resolution_clock::now();
    	long duration = duration_cast<microseconds>( te - ts ).count();

    	ts = high_resolution_clock::now();
    	searchMotif();
		te = high_resolution_clock::now();
    	duration = duration_cast<microseconds>( te - ts ).count();

		printResults();
    	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    	duration = duration_cast<microseconds>( t2 - t1 ).count();
		float sec = (float) duration / 1000000;
    	cout << "Motif Found: " << motifFound << endl;

		// steady_clock::time_point end= std::chrono::steady_clock::now();
		// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<std::endl;
		// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() <<std::endl;
	}

private:
	const DataSetParams *ds;
	const Config *config;
	int* candidateMotifs;
	int* currentNeighborhood;

	// TODO: change this to a long
	int*** blockMasks;
	int** currentBlockMasks;
	int* mismatches;

	long prefix, suffix;
	long currentMapping;

	long currentBlockCol; 
	long currentBlockRow;

	string foundMotifs;

	// contains all lmers per sequence
	// from sequence_tPrime to sequence_m
	vector< vector<long> > lmerMappings;

	// holds the string sequence in the dataset
	// that contains all the filteredLmerMappings
	int filteredLmerSequence;

	// holds the filtered lmers that will 
	// be used in further testing
	// the other bits in a block
	vector<long> filteredLmerMappings;

	// flags that tells if a block has at least 
	// one bit set in it.
	int* blockFlags;

	// count of number of lmers (bit represented)
	// in a block
	int* numberOfBitsPerBlock;

	// holds the minimum hamming distance between
	// the current lmer (generation of neighborhood)
	// versus all the processed lmers in the current sequence
	int currentMinDistance;

	// these are the processed lmers in sequence (generation phase)
	// that is within d-1 distance versus the current lmer processing
	vector< vector<long> > pruneLmers;

	int motifFound = 0;


	// -----------------------------------------------------------------------------
	//	Private Methods
	// -----------------------------------------------------------------------------

	/**
	 *	Precomputes the nonzero-bits(per 2 bits)
	 *	for the mismatches
	 */
	void generateMismatchesCount() {
	    long size =  (1 << 18);
	    mismatches = new int[size];
	    for (int i=0; i < size; i++) {
	        int result = i;
	        int distance = 0;
	        for (int k=0; k < ds->lengthOfMotif; k++) {
	            if ((result & 3) != 0) {
	                distance++;
	            }
	            result = result >> 2;
	        }
	        mismatches[i] = distance;
	    }
	}

	/**
	 *	Precomputes the patterns for each 4^5 possible lmers
	 */
	void generateBlockMasks() {
		blockMasks =  new int**[config->numberOfLmersInBlock];

		for (int i=0; i < config->numberOfLmersInBlock; i++) {
		    blockMasks[i] = new int*[config->blockDegree-1];

		    // instantiate pattern block array
		    for (int x=0; x < config->blockDegree-1; x++)
		        blockMasks[i][x] = new int[config->numberOfRowsInBlock];

		    for (int row=0; row < config->numberOfRowsInBlock; row++) {
		        for (int col=31; col > -1; col--) {
		            int distance = computeHammingDistance(i, row*config->TABLE_WIDTH+col);
		
		            for (int k=0; k < config->blockDegree - 1; k++) {
		                if (distance <= k+1) {
		                    blockMasks[i][k][row]++;
		                } 
		                if (col > 0) {
		                    blockMasks[i][k][row] = blockMasks[i][k][row] << 1;
		                }
		            }
		        }
		    }
		}
	}
	
	void collectCandidateMotifs() {

		// Tracks the 1-bits in each block
		numberOfBitsPerBlock = new int[config->numberOfCandidateMotifsRows/32];
		for (int i=0; i < config->numberOfCandidateMotifsRows/32; i++)
		{
			numberOfBitsPerBlock[i] = 0;
		}

    	high_resolution_clock::time_point ts = high_resolution_clock::now();
		// generate neighborhood of first string sequence in the dataset
		generateNeighborhood(0);

		// copy the whole neighborhood to the candidate motifs array
		candidateMotifs = currentNeighborhood;
		// clear the current neighbordhood for the succeeding generations
    	currentNeighborhood = new int[config->numberOfCandidateMotifsRows];

    	high_resolution_clock::time_point te = high_resolution_clock::now();
    	long  duration = duration_cast<microseconds>( te - ts ).count();
        // cout << "     >    (1) Duration (s): " << duration / 1000000 << endl;
 	
		blockFlags = new int[config->numberOfBlockRows];

		int middle_T = config->tPrime_2;
		if (config->tPrime < config->tPrime_2) middle_T = config->tPrime;

    	for(int i=1; i < middle_T; i++) {
    		high_resolution_clock::time_point ts = high_resolution_clock::now();
        	generateNeighborhood(i);

        	for (int j=0; j < config->numberOfCandidateMotifsRows; j++) {
        	    candidateMotifs[j] &= currentNeighborhood[j];
        	    currentNeighborhood[j] = 0;

        	    if (i == middle_T-1) {
		    	    if (candidateMotifs[j] != 0) {
		    	    	// count bits, 
		    	    	int val = candidateMotifs[j];
		    	    	for (int k=0; k < 32; k++) {
		    	    	    if ((val & 1) != 0) {
		    	    	    	numberOfBitsPerBlock[(int)(j/32)]++;
		    	    	    }
		    	    	    val = val >> 1;
		    	    	}
		    	    	setBlockFlag(j);
		    	    }
        	    }
        	}

        	high_resolution_clock::time_point te = high_resolution_clock::now();
        	long  duration = duration_cast<microseconds>( te - ts ).count();
        	// cout << "     >    (" << i+1 << ") Duration (s): " << duration / 1000000 << endl;
    	}

		for(int i=middle_T; i < config->tPrime; i++) {
			// cout << middle_T << endl;
			high_resolution_clock::time_point ts = high_resolution_clock::now();
	    	generateNeighborhoodB(i);
			blockFlags = new int[config->numberOfBlockRows];
	    	for (int j=0; j < config->numberOfCandidateMotifsRows; j++) {
	    	    candidateMotifs[j] &= currentNeighborhood[j];
	    	    currentNeighborhood[j] = 0;

	    	    if (candidateMotifs[j] != 0)
	    	    	setBlockFlag(j);
	    	}
    		// cout << "BLOCK FLAGS: " << countBlockFlags() << endl;

	    	high_resolution_clock::time_point te = high_resolution_clock::now();
	    	long  duration = duration_cast<microseconds>( te - ts ).count();
	    	// cout << "     >    (" << i+1 << ") Duration (s): " << duration / 1000000 << endl;
		}
	}

	void generateNeighborhood(int s) {
	    string currentSequence = ds->stringSequences[s];

	    // get the mapping for the first l characters
	    prefix = 0;
	    suffix = 0;
	    currentMapping = 0;
	    for (int i=0; i < ds->lengthOfMotif; i++) {
	        char c = currentSequence[i];
	        int base = 0;
	        switch(c) {
	            case 'C': base = 1; break;
	            case 'G': base = 2; break;
	            case 'T': base = 3; break;
	        }
	        if ( i < ds->lengthOfMotif - config->blockDegree)
	            prefix = (prefix << 2) + base;
	        else
	            suffix = (suffix << 2) + base;
	        currentMapping = (currentMapping << 2) + base;
	    }

	    currentBlockRow =  (int) (suffix / config->TABLE_WIDTH);
	    currentBlockCol =  (int) (suffix % config->TABLE_WIDTH);
	    currentBlockMasks = blockMasks[(int)suffix];
	    long blockStart = (int) (prefix * config->numberOfRowsInBlock);

	    for(int offset=0; offset < config->numberOfRowsInBlock; offset++) {
	        if(ds->numberOfAllowedMutations >= config->blockDegree)
	            currentNeighborhood[blockStart+offset] = INT_MAX;
	        else
	            currentNeighborhood[blockStart+offset] 
	                |= currentBlockMasks[ds->numberOfAllowedMutations - 1][offset];
	    }
	    addNeighbors(prefix, 0, ds->numberOfAllowedMutations);

	    // adjust one character at a time and get the mapping
	    for (int i=ds->lengthOfMotif; i < ds->lengthOfSequence; i++) {
	        prefix = (prefix & config->prefixMask) << 2;
	        suffix = (suffix & config->suffixMask) << 2;
	        currentMapping = (currentMapping & config->mask) << 2;

	        char c = currentSequence[i-config->blockDegree];
	        switch(c) {
	            case 'C': prefix += 1; break;
	            case 'G': prefix += 2; break;
	            case 'T': prefix += 3; break;
	        }

	        c = currentSequence[i];
	        switch(c) {
	            case 'C': suffix += 1; currentMapping+=1; break;
	            case 'G': suffix += 2; currentMapping+=2; break;
	            case 'T': suffix += 3; currentMapping+=3; break;
	        }

	        // collectProcessedLmersToPrune(currentSequence, currentMapping, i);

	        // if (currentMinDistance > 0) {
		        currentBlockRow = (int) (suffix / config->TABLE_WIDTH);
		        currentBlockCol = (int) (suffix % config->TABLE_WIDTH);
		        currentBlockMasks = blockMasks[(int)suffix];
		        long blockStart = (int) (prefix * config->numberOfRowsInBlock);

		        for(int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		            if(ds->numberOfAllowedMutations >= config->blockDegree)
		                currentNeighborhood[blockStart+offset] = INT_MAX;
		            else
		                currentNeighborhood[blockStart+offset] 
		                    |= currentBlockMasks[ds->numberOfAllowedMutations - 1][offset];
		        }

		        // it is possible to create a new addNeighbor function that checks if 
		        // rather than comparing all the time
		        addNeighbors(prefix, 0, ds->numberOfAllowedMutations);
	        // }
	    }
	}

	void addNeighbors(long prefix, int start, int allowedMutations) {

		int shift = (ds->lengthOfMotif - config->blockDegree- start) * 2;

		for (int k=start; k < ds->lengthOfMotif-config->blockDegree; ++k) 
		{
		    shift -= 2;

		    long alt1 = prefix ^ (((long) 1) << shift);
		    long alt2 = prefix ^ (((long) 2) << shift);
		    long alt3 = prefix ^ (((long) 3) << shift);

		    long mapping1 = (alt1 << (config->blockDegree*2)) + suffix;
		    long mapping2 = (alt2 << (config->blockDegree*2)) + suffix;
		    long mapping3 = (alt3 << (config->blockDegree*2)) + suffix;

		    int blockStart1 =  (int) alt1 << config->prefixShift;
		    int blockStart2 =  (int) alt2 << config->prefixShift;
		    int blockStart3 =  (int) alt3 << config->prefixShift;


		    // now check if a blockStart exists in one

		    // int level = (ds->numberOfAllowedMutations - allowedMutations + 1);
		    // if (pruneLmers.size() > 0 && level < ds->numberOfAllowedMutations && !pruneLmers[level].empty()) {
		    // 	// cout << level << endl;
		    // 	vector<long> lmers = pruneLmers[level];

		    // 	if (find(lmers.begin(), lmers.end(), mapping1) != lmers.end()) {
		    // 		continue;
		    // 	}
		    // 	if (find(lmers.begin(), lmers.end(), mapping2) != lmers.end()) {
		    // 		continue;
		    // 	}
		    // 	if (find(lmers.begin(), lmers.end(), mapping3) != lmers.end()) {
		    // 		continue;
		    // 	}
		    // }

		    int allow_d = allowedMutations - 1;
		    if (allow_d >= config->blockDegree)  {
		    	// if (currentNeighborhood[blockStart1] == INT_MAX) {}
		    	// else {
			        for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
			            currentNeighborhood[blockStart1 + offset] = INT_MAX;
			            currentNeighborhood[blockStart2 + offset] = INT_MAX;
			            currentNeighborhood[blockStart3 + offset] = INT_MAX;
			        }
		    	// }
		    } 
		    else if (allow_d > 0) {
		        for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		            currentNeighborhood[blockStart1 + offset] |= currentBlockMasks[allow_d - 1][offset];
		            currentNeighborhood[blockStart2 + offset] |= currentBlockMasks[allow_d - 1][offset];
		            currentNeighborhood[blockStart3 + offset] |= currentBlockMasks[allow_d - 1][offset];
		        }
		    } else {
		    	// if (((currentNeighborhood[blockStart1+currentBlockRow] >> currentBlockCol) & 1) == 1) {}
		    	// else {
			        currentNeighborhood[blockStart1 + currentBlockRow]  |= 1 << (currentBlockCol);
			        currentNeighborhood[blockStart2 + currentBlockRow]  |= 1 << (currentBlockCol);
			        currentNeighborhood[blockStart3 + currentBlockRow]  |= 1 << (currentBlockCol);
		    	// }
		    }

		    if (allowedMutations > 1) {
	    		addNeighbors(alt1, k+1, allowedMutations-1);
	    		addNeighbors(alt2, k+1, allowedMutations-1);
	    		addNeighbors(alt3, k+1, allowedMutations-1);
		    }
		}
	}

	void generateNeighborhoodB(int s) {
	    string currentSequence = ds->stringSequences[s];

	    // get the mapping for the first l characters
	    prefix = 0;
	    suffix = 0;
	    currentMapping = 0;
	    for (int i=0; i < ds->lengthOfMotif; i++) {
	        char c = currentSequence[i];
	        int base = 0;
	        switch(c) {
	            case 'C': base = 1; break;
	            case 'G': base = 2; break;
	            case 'T': base = 3; break;
	        }
	        if ( i < ds->lengthOfMotif - config->blockDegree)
	            prefix = (prefix << 2) + base;
	        else
	            suffix = (suffix << 2) + base;

	        currentMapping = (currentMapping << 2) + base;
	    }

	    currentBlockRow = (int) (suffix / config->TABLE_WIDTH);
	    currentBlockCol = (int) (suffix % config->TABLE_WIDTH);
	    currentBlockMasks = blockMasks[(int)suffix];
	    int blockStart = (int) (prefix * config->numberOfRowsInBlock);

	    for(int offset=0; offset < config->numberOfRowsInBlock; offset++) {
	        if(ds->numberOfAllowedMutations >= config->blockDegree)
	            currentNeighborhood[blockStart+offset] = INT_MAX;
	        else
	            currentNeighborhood[blockStart+offset] 
	                |= currentBlockMasks[ds->numberOfAllowedMutations - 1][offset];
	    }
	    addNeighborsB(prefix, 0, ds->numberOfAllowedMutations);

	    // adjust one character at a time and get the mapping
	    for (int i=ds->lengthOfMotif; i < ds->lengthOfSequence; i++) {
	        prefix = (prefix & config->prefixMask) << 2;
	        suffix = (suffix & config->suffixMask) << 2;
	        currentMapping = (currentMapping & config->mask) << 2;

	        char c = currentSequence[i-config->blockDegree];
	        switch(c) {
	            case 'C': prefix += 1; break;
	            case 'G': prefix += 2; break;
	            case 'T': prefix += 3; break;
	        }

	        c = currentSequence[i];
	        switch(c) {
	            case 'C': suffix += 1; currentMapping += 1; break;
	            case 'G': suffix += 2; currentMapping += 2; break;
	            case 'T': suffix += 3; currentMapping += 3; break;
	        }

	        currentBlockRow = (int) (suffix / config->TABLE_WIDTH);
	        currentBlockCol = (int) (suffix % config->TABLE_WIDTH);
	        currentBlockMasks = blockMasks[(int)suffix];
	        int blockStart = (int) (prefix * config->numberOfRowsInBlock);

	        // collectProcessedLmersToPrune(currentSequence, currentMapping, i);
	        // if (currentMinDistance > 0) {
		        for(int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		            if(ds->numberOfAllowedMutations >= config->blockDegree)
		                currentNeighborhood[blockStart+offset] = INT_MAX;
		            else
		                currentNeighborhood[blockStart+offset] 
		                    |= currentBlockMasks[ds->numberOfAllowedMutations - 1][offset];
		        }
		        addNeighborsB(prefix, 0, ds->numberOfAllowedMutations);
	        // }
	    }
	}

	void addNeighborsB(long prefix, int start, int allowedMutations) {

	    int shift = (ds->lengthOfMotif - config->blockDegree- start) * 2;

	    for (int k=start; k < ds->lengthOfMotif-config->blockDegree; ++k) {
	        shift -= 2;

	        long alt1 = prefix ^ (((long) 1) << shift);
	        long alt2 = prefix ^ (((long) 2) << shift);
	        long alt3 = prefix ^ (((long) 3) << shift);

	        long mapping1 = (alt1 << (config->blockDegree*2)) + suffix;
		    long mapping2 = (alt2 << (config->blockDegree*2)) + suffix;
		    long mapping3 = (alt3 << (config->blockDegree*2)) + suffix;

	        int blockStart1 = (int) alt1 << config->prefixShift;
	        int blockStart2 = (int) alt2 << config->prefixShift;
	        int blockStart3 = (int) alt3 << config->prefixShift;

	        // int level = (ds->numberOfAllowedMutations - allowedMutations + 1);
	        // if (pruneLmers.size() > 0 && level < ds->numberOfAllowedMutations && !pruneLmers[level].empty()) {
	        // 	vector<long> lmers = pruneLmers[level];

	        // 	if (find(lmers.begin(), lmers.end(), mapping1) != lmers.end()) {
	        // 		continue;
	        // 	}
	        // 	if (find(lmers.begin(), lmers.end(), mapping2) != lmers.end()) {
	        // 		continue;
	        // 	}
	        // 	if (find(lmers.begin(), lmers.end(), mapping3) != lmers.end()) {
	        // 		continue;
	        // 	}
	        // }

	        int allow_d = allowedMutations - 1;
	        if (allow_d >= config->blockDegree)  {
	        	if (isSetBlockFlag(blockStart1) )
		            for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                currentNeighborhood[blockStart1 + offset] = INT_MAX;
		            }
	        	if (isSetBlockFlag(blockStart2) )
		            for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                currentNeighborhood[blockStart2 + offset] = INT_MAX;
		            }
	        	if (isSetBlockFlag(blockStart3) )
		            for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                currentNeighborhood[blockStart3 + offset] = INT_MAX;
		            }
	        } 
	        else if (allow_d > 0) {
	        	if (isSetBlockFlag(blockStart1) )
		            for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                currentNeighborhood[blockStart1 + offset] |= currentBlockMasks[allow_d - 1][offset];
		            }
	        	if (isSetBlockFlag(blockStart2) )
		            for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                currentNeighborhood[blockStart2 + offset] |= currentBlockMasks[allow_d - 1][offset];
		            }
	        	if (isSetBlockFlag(blockStart3) )
		            for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                currentNeighborhood[blockStart3 + offset] |= currentBlockMasks[allow_d - 1][offset];
		            }
	        } else {
	            currentNeighborhood[blockStart1 + currentBlockRow]  |= 1 << (currentBlockCol);
	            currentNeighborhood[blockStart2 + currentBlockRow]  |= 1 << (currentBlockCol);
	            currentNeighborhood[blockStart3 + currentBlockRow]  |= 1 << (currentBlockCol);
	        }


	        if (allowedMutations > 2 && false) {
	            int shift2 = (ds->lengthOfMotif - config->blockDegree - k-1) * 2;
	            for (int i=k+1; i < ds->lengthOfMotif-config->blockDegree; i++) {
	                shift2 -= 2;

	                long alt1Prime1 = alt1 ^ (((long) 1) << shift2);
	                long alt1Prime2 = alt1 ^ (((long) 2) << shift2);
	                long alt1Prime3 = alt1 ^ (((long) 3) << shift2);

	                long alt2Prime1 = alt2 ^ (((long) 1) << shift2);
	                long alt2Prime2 = alt2 ^ (((long) 2) << shift2);
	                long alt2Prime3 = alt2 ^ (((long) 3) << shift2);

	                long alt3Prime1 = alt3 ^ (((long) 1) << shift2);
	                long alt3Prime2 = alt3 ^ (((long) 2) << shift2);
	                long alt3Prime3 = alt3 ^ (((long) 3) << shift2);


	                int blockStart11 = (int) alt1Prime1 << config->prefixShift;
	                int blockStart12 = (int) alt1Prime2 << config->prefixShift;
	                int blockStart13 = (int) alt1Prime3 << config->prefixShift;

	                int blockStart21 = (int) alt2Prime1 << config->prefixShift;
	                int blockStart22 = (int) alt2Prime2 << config->prefixShift;
	                int blockStart23 = (int) alt2Prime3 << config->prefixShift;

	                int blockStart31 = (int) alt3Prime1 << config->prefixShift;
	                int blockStart32 = (int) alt3Prime2 << config->prefixShift;
	                int blockStart33 = (int) alt3Prime3 << config->prefixShift;

	                int allow_d2 = allowedMutations - 2;
	                if (allow_d2 >= config->blockDegree)  {
	                	if (isSetBlockFlag(blockStart11))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart11 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart12))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart12 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart13))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart13 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart21))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart21 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart22))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart22 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart23))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart23 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart31))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart31 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart32))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart32 + offset] = INT_MAX;
		                    }
	                	if (isSetBlockFlag(blockStart33))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart33 + offset] = INT_MAX;
		                    }
	                } 
	                else if (allow_d2 > 0) {
	                    int allow_d2_prime = allow_d2 - 1;
	                	if (isSetBlockFlag(blockStart11))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart11 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }
	                	if (isSetBlockFlag(blockStart12))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart12 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }
	                	if (isSetBlockFlag(blockStart13))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart13 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }
		                if (isSetBlockFlag(blockStart21))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart21 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }
		                if (isSetBlockFlag(blockStart22))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart22 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }
		                if (isSetBlockFlag(blockStart23))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart23 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }
		                if (isSetBlockFlag(blockStart31))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart31 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }

		                if (isSetBlockFlag(blockStart32))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart32 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }

		                if (isSetBlockFlag(blockStart33))
		                    for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		                        currentNeighborhood[blockStart33 + offset] |= currentBlockMasks[allow_d2_prime][offset];
		                    }
	                } else {
	                    currentNeighborhood[blockStart11 + currentBlockRow]  |= 1 << (currentBlockCol);
	                    currentNeighborhood[blockStart12 + currentBlockRow]  |= 1 << (currentBlockCol);
	                    currentNeighborhood[blockStart13 + currentBlockRow]  |= 1 << (currentBlockCol);

	                    currentNeighborhood[blockStart21 + currentBlockRow]  |= 1 << (currentBlockCol);
	                    currentNeighborhood[blockStart22 + currentBlockRow]  |= 1 << (currentBlockCol);
	                    currentNeighborhood[blockStart23 + currentBlockRow]  |= 1 << (currentBlockCol);

	                    currentNeighborhood[blockStart31 + currentBlockRow]  |= 1 << (currentBlockCol);
	                    currentNeighborhood[blockStart32 + currentBlockRow]  |= 1 << (currentBlockCol);
	                    currentNeighborhood[blockStart33 + currentBlockRow]  |= 1 << (currentBlockCol);
	                }

	                int iPrime = i+1;
	                addNeighborsB(alt1Prime1, iPrime, allow_d2);
	                addNeighborsB(alt1Prime2, iPrime, allow_d2);
	                addNeighborsB(alt1Prime3, iPrime, allow_d2);

	                addNeighborsB(alt2Prime1, iPrime, allow_d2);
	                addNeighborsB(alt2Prime2, iPrime, allow_d2);
	                addNeighborsB(alt2Prime3, iPrime, allow_d2);

	                addNeighborsB(alt3Prime1, iPrime, allow_d2);
	                addNeighborsB(alt3Prime2, iPrime, allow_d2);
	                addNeighborsB(alt3Prime3, iPrime, allow_d2);

	            }
	        }
	        else if (allowedMutations > 1) {
	            addNeighborsB(alt1, k+1, allowedMutations-1);
	            addNeighborsB(alt2, k+1, allowedMutations-1);
	            addNeighborsB(alt3, k+1, allowedMutations-1);
	        }
	    }
	}

	void collectProcessedLmersToPrune(string currentSequence, long mapping, int pos)
	{
		long p;
		int min = INT_MAX;

		pruneLmers.clear();
		pruneLmers.resize(ds->numberOfAllowedMutations);

		// process the first lmer
		for (int i=0; i < ds->lengthOfMotif; i++) {
	        char c = currentSequence[i];
	        int base = 0;
	        switch(c) {
	            case 'C': base = 1; break;
	            case 'G': base = 2; break;
	            case 'T': base = 3; break;
	        }
	        p = (p << 2) + base;
	    }

	    int hd = computeHammingDistance(mapping, p);
	    if (hd < ds->numberOfAllowedMutations) { 
	    	pruneLmers[hd].push_back(p);
	    }
	    if (min > hd) {
	    	min = hd;
	    }

	     // adjust one character at a time and get the mapping
	    for (int i=ds->lengthOfMotif; i < pos; i++) {
	        p = (p & config->mask) << 2;

	        char c = currentSequence[i-config->blockDegree];
	        switch(c) {
	            case 'C': p+= 1; break;
	            case 'G': p+= 2; break;
	            case 'T': p+= 3; break;
	        }

		    int hd = computeHammingDistance(mapping, p);
		    if (hd < ds->numberOfAllowedMutations) { 
		    	pruneLmers[hd].push_back(p);
		    }
		    if (min > hd)
		    {
		    	min = hd;
		    }
	    }

	    currentMinDistance = min;
	}

	void transformLmerSequences() {
	    // collect local mappings
	    vector<long> mappings;

	    for (int i=config->tPrime; i < ds->numberOfSequences; i++) {
	        mappings.clear();

	        string currentSequence = ds->stringSequences[i];
	        long mapping; 

	        for (int j=0; j < ds->lengthOfMotif; j++) {
	            char c = currentSequence[j];
	            int base = 0;
	            switch(c) {
	                case 'C': base=1; break;
	                case 'G': base=2; break;
	                case 'T': base=3; break;
	            }
	            mapping = (mapping << 2) + base;
	        }
	        mappings.push_back(mapping);

	        int k=0;
	        for(int j=ds->lengthOfMotif; j < ds->lengthOfSequence;j++) {
	            char c = currentSequence[j];
	            int base = 0;
	            switch(c) {
	                case 'C': base=1; break;
	                case 'G': base=2; break;
	                case 'T': base=3; break;
	            }
	            mapping = ((mapping & config->mask) << 2) + base;
	            mappings.push_back(mapping);
	        }
	        lmerMappings.push_back(mappings);
	    }
	}

	void searchMotif() {
	    int value, numMotifs = 0;
	    foundMotifs = "";

	    // cout << trimmedRows.size() << endl;

	    // for (int i=0; i < trimmedRows.size(); i++) {
	    for (int i=0; i < config->numberOfCandidateMotifsRows; i++) {
	    	// if this is the start of the block
	    	// initialize/clear the filtered lmer mappings
	    	// TODO: change this to i % 32 -- since this is dependent 
	    	// on the number of bits in the INTEGER
	    	if (i % (1 << config->blockDegree) == 0) {
	    		filteredLmerMappings.clear();
	    		filteredLmerSequence = -1;
	    	}

	    	// if value is zero then proceed to the next row	
	        if ( (value = candidateMotifs[i]) == 0 ) {
	            continue;
	        }

	        long base = ((long) i) << 5;

	        // int row = trimmedRows[i];
	        // value = candidateMotifs[row];
	        // long base = ((long) row) << 5;

	        for (int j=0; j < 32; j++) {
	            if ((value & 1) != 0) {
	                long candidate = base + j;

	                // if there are less than `2` lmers to be tested in a block
	                // theres no need for the filtering speedup since, it won't be used anyway
	                if (numberOfBitsPerBlock[(int)(i/32)] <= 1) {
	                	if (isMotif(candidate)) {
	                    	string motif = decode(candidate, ds->lengthOfMotif);
	                        foundMotifs += " " + motif;
	                        numMotifs++;

	                        if (motif == ds->plantedMotif) motifFound = 1;
		                }
	                }

	                // if this is the first lmer in the block to be tested
	                else if (filteredLmerSequence == -1) {
		                if (isMotifInitializeFilter(candidate)) {
	                    	string motif = decode(candidate, ds->lengthOfMotif);
	                        foundMotifs += " " + motif;
	                        numMotifs++;

	                        if (motif == ds->plantedMotif) motifFound = 1;
		                }
	                }

	                // else if there exists a filtered lmers and a sequences, use that
	                else {
	                	if (isMotifUseFilter(candidate)) {
	                    	string motif = decode(candidate, ds->lengthOfMotif);
	                        foundMotifs += " " + motif;
	                        numMotifs++;

	                        if (motif == ds->plantedMotif) motifFound = 1;
		                }
	                }
	            }
	            value = value >> 1;
	        }
	    }
	}

	/**
	 *  Given a mapping, check the rest of the sequences if that mapping
	 *  exists in all of it.
	 *
	 *	This isMotif function assumes that the filteredLmerMappings is not yet initialized
	 */
	bool isMotifInitializeFilter(long mapping) {
	    for (int i=0; i < lmerMappings.size(); i++) {
	        bool found = false;

	        filteredLmerSequence = i;
	    	filteredLmerMappings.clear(); 
	        for (int j=0; j < config->numberOfPossibleLmersInSequence; j++) {
	        	// current lmer for comparison
	            long lmer = lmerMappings[i][j];
	            // compute hammingdistance of current lmer vs the candidate motif
	            int hammingDistance = computeHammingDistance(mapping, lmer);

	            // collect all lmers in the sequence that is within
	            // d + blockDegree distance between the mapping
	            if (hammingDistance <= (config->blockDegree + ds->numberOfAllowedMutations)) 
	            {
	            	// if our current min is beaten
	            	// then clear the list, and add the current champ
	            	filteredLmerMappings.push_back(lmer);
	            } 

	            if (hammingDistance <= ds->numberOfAllowedMutations) {
	                found = true;
	                break;
	            }
	        }

	        // if there is a sequence where the mapping is not present
	        // return false
	        if ( !found ) {
	            return false;
	        }
	    }

	    // if the mapping is a motif we can disregard
	    // the collected filter lmer mappings
	    filteredLmerMappings.clear(); 
	    filteredLmerSequence = -1;
	    
	    return true;
	}

	/**
	 *  Given a mapping, check the rest of the sequences if that mapping
	 *  exists in all of it.
	 *
	 *	This isMotif function assumes that the filteredLmerMappings is already initialized 
	 *	and uses it.
	 */
	bool isMotifUseFilter(long mapping) {
		// check if there is a d-distance between the mapping versus all the lmers 
		// in the filtered lmer mappings set
		bool found = false;
		
		// cout << filteredLmerMappings.size() << endl;
		for (const auto &lmer : filteredLmerMappings)
		{
	        int hammingDistance = computeHammingDistance(mapping, lmer);
	        if (hammingDistance <= ds->numberOfAllowedMutations) {
	        	found = true;
	        	break;
	        }
		}

		// if there are no d-distance lmer in the filtered set
		// we can say that the mapping is not a motif.
		if (found == false) 
		{
			return false;
		}

		// if there is a d-distance lmer in the filtered set, 
		// check if there is for the rest of the string sequences
	    for (int i=0; i < lmerMappings.size(); i++) {
	        bool found = false;

	        // move to the next sequence
	        if (i == filteredLmerSequence) continue;

	        for (int j=0; j < config->numberOfPossibleLmersInSequence; j++) {
	        	// current lmer for comparison
	            long lmer = lmerMappings[i][j];
	            // compute hammingdistance of current lmer vs the candidate motif
	            int hammingDistance = computeHammingDistance(mapping, lmer);

	            if (hammingDistance <= ds->numberOfAllowedMutations) {
	                found = true;
	                break;
	            }
	        }

	        // if there is a sequence where the mapping is not present
	        // return false
	        if ( !found ) {
	            return false;
	        }
	    }
	    return true;
	}

	/**
	 *  Given a mapping, check the rest of the sequences if that mapping
	 *  exists in all of it.
	 *
	 *	This isMotif function assumes that the filteredLmerMappings is already initialized 
	 *	and uses it.
	 */
	bool isMotif(long mapping) {

		// if there is a d-distance lmer in the filtered set, 
		// check if there is for the rest of the string sequences
	    for (int i=0; i < lmerMappings.size(); i++) {
	        bool found = false;

	        for (int j=0; j < config->numberOfPossibleLmersInSequence; j++) {
	        	// current lmer for comparison
	            long lmer = lmerMappings[i][j];
	            // compute hammingdistance of current lmer vs the candidate motif
	            int hammingDistance = computeHammingDistance(mapping, lmer);

	            if (hammingDistance <= ds->numberOfAllowedMutations) {
	                found = true;
	                break;
	            }
	        }

	        // if there is a sequence where the mapping is not present
	        // return false
	        if ( !found ) {
	            return false;
	        }
	    }
	    return true;
	}

	// -----------------------------------------------------------------------------
	//	Utility Methods
	// -----------------------------------------------------------------------------
	int computeHD(long lmer1, long lmer2) {
	    int distance = 0;
	    long result = lmer1 ^ lmer2;

	    for (int i=0; i < ds->lengthOfMotif; i++) {
	        if ((result & 3) != 0) {
	            distance++;
	        }
	        result = result >> 2;
	    }
	    return distance;
	}

	int computeHammingDistance(long lmer1, long lmer2) {
	    int distance = 0;
	    long result = lmer1 ^ lmer2;

	    int c = 2;
	    while (c--) {
	        int i = (result & ((1 << 18)-1));
	        distance += mismatches[i];
	        result = result >> 18;
	    }
	    return distance;
	}


	string decode(long mapping, int strlen) {
	    string decoding = "";

	    for (int i=0; i < strlen; i++) {
	        int base = (int) mapping & 3;

	        switch(base) {
	            case 0: decoding = "A" + decoding; break;
	            case 1: decoding = "C" + decoding; break;
	            case 2: decoding = "G" + decoding; break;
	            case 3: decoding = "T" + decoding; break;
	        }
	        mapping = mapping >> 2;
	    }
	    return decoding;
	}

	void setBlockFlag(int blockStart) {
	    blockFlags[(int)blockStart/32] = 1;
	}

	int isSetBlockFlag(int blockRow) {
		int bs = (int)(blockRow/32);
		if (bs < config->numberOfBlockRows-1 && blockRow % 32 != 0)
			return blockFlags[bs] || blockFlags[bs+1];
	    return blockFlags[bs];
	}

	int countBlockFlags() {
	    int count = 0;
	    for (int i=0; i < (config->numberOfBlockRows); i++) {
	        if (blockFlags[i] > 0) count++;
	    }
	    return count;
	}

	void printResults() {
	    cout << "--------------------------------------------------------" << endl;
	    cout << " RESULTS: " << endl;
	    // cout << "Dataset: " << inputFileHeader << inputFileName << endl;
	    cout << "   > Planted Motif: " << ds->plantedMotif << endl;
	    cout << "   > Found Motifs: " << foundMotifs << endl;
	    cout << "-------------------------------------------------------" << endl;
	}

	void countSurvivingCandidateMotifs() {
		long count = 0;
		for (int i=0; i < config->numberOfCandidateMotifsRows; i++) {
			int row = candidateMotifs[i];

			if (row == 0) continue;
			for (int j=0; j < 32; j++) {
				int val = row & 1;
				if (val == 1) count++;

				row = row >> 1;
			}
		}
		cout << count << endl;
	}

};