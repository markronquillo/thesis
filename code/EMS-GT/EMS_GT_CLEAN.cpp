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


		// cout << "======================================================" << endl;
		// cout << " EMS_GT v2.0" << endl;
		// cout << "======================================================" << endl;
		// cout << " CONFIGURATION" << endl;
		// cout << " tPrime:\t\t " << config->tPrime << endl;
		// cout << " tPrime_2:\t\t " << config->tPrime_2 << endl;
		// cout << " l:\t\t " << ds->lengthOfMotif << endl;
		// cout << " d:\t\t " << ds->numberOfAllowedMutations << endl;
		// cout << "======================================================" << endl;

		// cout << "" << endl;
		// cout << " GENERATE CANDIDATE MOTIFS" << endl;
		// cout << "------------------------------------------------------" << endl;

		// generateMismatchesCount();
		// cout << "   > Hamming Distance count generated" << endl;

		generateBlockMasks();
		// cout << "   > BlockMasks generated" << endl;

		collectCandidateMotifs();
		// cout << "   > Candidate Motifs generated" << endl;

		// cout << "" << endl;
		// cout << " TEST CANDIDATE MOTIFS" << endl;
		// cout << "------------------------------------------------------" << endl;

    	high_resolution_clock::time_point ts = high_resolution_clock::now();
		transformLmerSequences();
    	high_resolution_clock::time_point te = high_resolution_clock::now();
    	long duration = duration_cast<microseconds>( te - ts ).count();
		// cout << "   > Lmer Sequences generated: " << duration / 1000000 << " (s)" << endl;

    	ts = high_resolution_clock::now();
    	searchMotif();
		te = high_resolution_clock::now();
    	duration = duration_cast<microseconds>( te - ts ).count();
		// cout << "   > Search Motif is done: " << duration / 1000000 << " (s)" << endl;

		printResults();
    	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    	duration = duration_cast<microseconds>( t2 - t1 ).count();
		float sec = (float) duration / 1000000;
    	cout << "Duration: " << duration << endl;
    	cout << "Duration (s): " <<  sec << endl;

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


	// -----------------------------------------------------------------------------
	//	Private Methods
	// -----------------------------------------------------------------------------

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
		            int distance = computeHD(i, row*config->TABLE_WIDTH+col);
		
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

    	for(int i=1; i < config->tPrime; i++) {
    		high_resolution_clock::time_point ts = high_resolution_clock::now();
        	generateNeighborhood(i);

        	for (int j=0; j < config->numberOfCandidateMotifsRows; j++) {
        	    candidateMotifs[j] &= currentNeighborhood[j];
        	    currentNeighborhood[j] = 0;

        	}

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

		    int allow_d = allowedMutations - 1;
		    if (allow_d >= config->blockDegree)  {
		        for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		            currentNeighborhood[blockStart1 + offset] = INT_MAX;
		            currentNeighborhood[blockStart2 + offset] = INT_MAX;
		            currentNeighborhood[blockStart3 + offset] = INT_MAX;
		        }
		    } 
		    else if (allow_d > 0) {
		        for (int offset=0; offset < config->numberOfRowsInBlock; offset++) {
		            currentNeighborhood[blockStart1 + offset] |= currentBlockMasks[allow_d - 1][offset];
		            currentNeighborhood[blockStart2 + offset] |= currentBlockMasks[allow_d - 1][offset];
		            currentNeighborhood[blockStart3 + offset] |= currentBlockMasks[allow_d - 1][offset];
		        }
		    } else {
		        currentNeighborhood[blockStart1 + currentBlockRow]  |= 1 << (currentBlockCol);
		        currentNeighborhood[blockStart2 + currentBlockRow]  |= 1 << (currentBlockCol);
		        currentNeighborhood[blockStart3 + currentBlockRow]  |= 1 << (currentBlockCol);
		    }

		    if (allowedMutations > 1) {
	    		addNeighbors(alt1, k+1, allowedMutations-1);
	    		addNeighbors(alt2, k+1, allowedMutations-1);
	    		addNeighbors(alt3, k+1, allowedMutations-1);
		    }
		}
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

	    for (int i=0; i < config->numberOfCandidateMotifsRows; i++) {
	        if ( (value = candidateMotifs[i]) == 0 ) {
	            continue;
	        }
	        long base = ((long) i) << 5;
	        for (int j=0; j < 32; j++) {
	            if ((value & 1) != 0) {
	                long candidate = base + j;
	                if (isMotif(candidate)) {
	                    foundMotifs += " " + decode(candidate, ds->lengthOfMotif);
	                    numMotifs++;
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
	            int hammingDistance = computeHD(mapping, lmer);

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

	void printResults() {
	    cout << "--------------------------------------------------------" << endl;
	    cout << " RESULTS: " << endl;
	    // cout << "Dataset: " << inputFileHeader << inputFileName << endl;
	    cout << "   > Planted Motif: " << ds->plantedMotif << endl;
	    cout << "   > Found Motifs: " << foundMotifs << endl;
	    cout << "-------------------------------------------------------" << endl;
	}
};