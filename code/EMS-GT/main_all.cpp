#include <iostream>
#include <vector>
#include "DataSetParams.cpp"
#include "Config.cpp"
#include "EMS_GT.cpp"
#include "utils.cpp"

using namespace std;

// TODO:
// make a generic EMS_GT, using long and int (32 vs 64)
//
int main(int argc, char* argv[]) {
    // string inputFileHeader = "../datasets/";
    string inputFileName;
    int tPrime = -1;
    string plantedMotif = "";

    DataSetParams dsParams;

    if (argc > 1)
        inputFileName = argv[1];
    else
        inputFileName = "17,6,0"; // hardcoded muna

    // put the dataset config and put it in the dsConfig struct
    readInput(inputFileName, &dsParams);

    if (argc > 3)
        dsParams.lengthOfMotif = atoi(argv[3]);
    if (argc > 4)
        dsParams.numberOfAllowedMutations = atoi(argv[4]);

    // other config details are based on the dataset parameters
    // so we need to pass it to
    Config config = Config(&dsParams);

    // adjust tPrime (nPrime) value depending on the dataset
    // hardcoded
    switch(dsParams.lengthOfMotif)  {
        case 13: config.tPrime = 9; break;
        case 15: config.tPrime = 8; break;
        case 17: config.tPrime = 7; break;
    }

    // override change if argument is passed
    if (argc > 2)
        config.tPrime = atoi(argv[2]);

    if (argc > 5)
        plantedMotif = argv[5];
    
    config.tPrime_2 = 6;

    // EMS_GT ems = EMS_GT(&dsParams, &config, plantedMotif);
    EMS_GT ems = EMS_GT(&dsParams, &config);
    ems.start();

    return 0;
}
