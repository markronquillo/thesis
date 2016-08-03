// 540682650946920004

#include <iostream>

using namespace std;

int* mismatches;


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
        for (int k=0; k < 15; k++) {
            if ((result & 3) != 0) {
                distance++;
            }
            result = result >> 2;
        }
        mismatches[i] = distance;
    }
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

int computeHD(long lmer1, long lmer2) {
    int distance = 0;
    long result = lmer1 ^ lmer2;

    for (int i=0; i < 15; i++) {
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
        cout << mismatches[i] << endl;
        result = result >> 18;
    }
    return distance;
}

long encode(string lmer, int strlen) {
	long mapping = 0;
	for (int x=0; x< strlen; x++) {
		char c = lmer[x];
		cout << c << " ";
		int base = 0;
		switch(c) {
			case 'C': base = 1; break;
			case 'G': base = 2; break;
			case 'T': base = 3; break;
		}
		mapping = (mapping << 2) + base;
	}
	return mapping;
}


int main() {
	generateMismatchesCount();

	string l1 = "GACTGGACCTGCACA";
	string l2 = "GATTGGATGCGCACG";
	cout << encode(l1, 15) << endl;
	cout << encode(l2, 15) << endl;

	cout << decode(568876612, 15) << endl;
	cout << decode(602465862, 15) << endl;

	cout << computeHammingDistance(568876612, 602465862) << endl;
	// cout << decode(540682650946920004, 15) << endl;
}
