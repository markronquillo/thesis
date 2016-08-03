#include<iostream>

using namespace std;

int main() {

    for (int x=0; x < 6; x++) {
        for (int y =0; y < 32; y++)
        {
            int r = rand() % 2;
            cout << r << ",";
        }
        cout << endl;
    }
}
