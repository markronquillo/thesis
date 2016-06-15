#include <iostream>
#include <vector>
#include "DataSetParams.cpp"
#include "Config.cpp"
#include "EMS_GT_CLEAN.cpp"
#include "utils.cpp"

using namespace std;

// TODO:
// make a generic EMS_GT, using long and int (32 vs 64)
// 
int main(int argc, char* argv[]) {
	// string inputFileHeader = "../datasets/";
	string inputFileName;
	int tPrime = -1;

    DataSetParams dsParams;

    if (argc > 1)   
        inputFileName = argv[1];
    else
        inputFileName = "17,6,0"; // hardcoded muna

    if (argc > 2)
    	tPrime = atoi(argv[2]);

    // put the dataset config and put it in the dsConfig struct
	readInput(inputFileName, &dsParams);

	// other config details are based on the dataset parameters
	// so we need to pass it to
    Config config = Config(&dsParams);

    // change if argument is passed
    if (tPrime != -1)
    	config.tPrime = tPrime;

	EMS_GT ems = EMS_GT(&dsParams, &config);
	ems.start();

	return 0;
}