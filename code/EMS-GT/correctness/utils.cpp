#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>


using namespace std;


/**
 *	Stores the data in the DataSetConfig
 */
int readInput(string filename, DataSetParams *ds) {
 	ifstream inputFile(filename);

    if (!inputFile) {
        cout << "Can't find input file." << endl;
        // exit(1);
        return 0;
    }

    string line;

    // read number of sequences
    getline(inputFile, line);
    ds->numberOfSequences = atoi(line.c_str());

    // read length of each sequences
    getline(inputFile, line);
    ds->lengthOfSequence = atoi(line.c_str());

    // read length of planted motif
    getline(inputFile, line);
    ds->lengthOfMotif = atoi(line.c_str());

    // read number of allowed mutation
    getline(inputFile, line);
    ds->numberOfAllowedMutations = atoi(line.c_str());

    // planted motif alignments
    getline(inputFile, line);
    istringstream iss(line);
    int n;
    while (iss >> n) {
        ds->plantedAlignments.push_back(n);
    }

    // planted motif
    getline(inputFile, ds->plantedMotif);

    // read the sequences
    for (int i=0; i < ds->numberOfSequences; i++) {
        getline(inputFile, line);
        ds->stringSequences.push_back(line);
    }

    return 1;
}


void printConfig(Config *c) {
    cout << "-------------------------------------------------------------------" << endl;
    cout << "PRINTING CONFIG:" << endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << "Number of Lmers in Block: " << c->numberOfLmersInBlock << endl;
    cout << "Number of Rows in block: " << c->numberOfRowsInBlock << endl;
    cout << "-------------------------------------------------------------------" << endl;
}


