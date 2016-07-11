#include <iostream>

using namespace std;
int main() {
    char c = 'c';

    switch (c) {
        case 'C':
            cout << "uppercase"; break;
        case 'c':
            cout << "lowercase"; break;
    }
}
