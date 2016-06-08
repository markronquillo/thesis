#include <vector>

using namespace std;

struct DataSetParams {
	int numberOfSequences;
	int lengthOfMotif;
	int lengthOfSequence;
	int numberOfAllowedMutations;

	vector<int> plantedAlignments;
	string plantedMotif;
	vector<string> stringSequences;
};