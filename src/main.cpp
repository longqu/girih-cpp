#include <exception>
#include <iostream>

#include "utilities/affinity/affinity.hpp"

using namespace std;

int main() {
    string str;
    girih::Affinity affinity;

    affinity.get_string(str,0);
    cout << str.c_str() << endl;

//    affinity.set_openmp(3);

//    affinity.get_string(str,0);
 //   cout << str.c_str() << endl;

//    affinity.get_string(str,0);
 //   cout << str.c_str() << endl;

//    affinity.set_openmp_nested(2,2);

//    affinity.get_string(str,0);
//    cout << str.c_str() << endl;

//    affinity.get_string(str,0);
//    cout << str.c_str() << endl;
}