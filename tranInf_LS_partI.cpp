//
//  main.cpp
//  transInf_LS_partI
//
//  Created by YZhan on 5/16/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#define Q 101    // max test total number of qubits = Q - 1
// #define NN 500 // estimation of number of possibilities

// int b[NN][Q] = {0}; // branching parameters b[total number][i-th possibility][branching para]
// int d[NN] = {0};   // depth of tree d[total number][i-th possibility]
int branPara[Q] = {0};  // record current set of branching parameters for output
ofstream fout[Q];   // for output of branching parameters
int testQBegin = 7;    // starting point of testQ
long long Index = 1;  // number of possible tree structures for a certain Q
void Factorize(int totNum, int avaiNum, int parentLevPara, int level);

int main() {
    // all possible branching parameters for given Q range
    long long possibility[Q] = {0};   // record the number of possibilities for all Q
    
//#pragma omp parallel for
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        Index = 1;  // reset # of possible structures
        // factorize to find all possibilities
        int rootPara = 1;   // only one root
        branPara[0] = 1;    // level 0 has one photon
        // open a file for output
//        string filename, fn1, fn2, fn3;
//        fn1 = "/Users/zhanyuan/Desktop/treeBranPara/Q";
//        fn2 = to_string(testQ);
//        fn3 = ".txt";
//        filename = fn1 + fn2 + fn3;
//        fout[testQ - testQBegin].open(filename);
        
        Factorize(testQ, testQ, rootPara, 1);
        Index--;
//        fout[testQ - testQBegin].close();
        possibility[testQ - testQBegin] = Index;
        cout<< "Number of possible branching parameter sets is "<< Index<< " for a total number of qubits Q="<< testQ<< endl;
    }
//    // output number of possibilities to a file
//    string ffilename = "/Users/zhanyuan/Desktop/treeBranPara/numPossibility.txt";
//    ofstream ffout(ffilename);
//    for(int i = 0; i < Q - testQBegin; ++i) {
//        ffout<< possibility[i]<< " ";
//    }
//    ffout.close();
    return 0;
}

void Factorize(int totNum, int avaiNum, int parentLevPara, int level) {
    int Inte;   // the integer that needs factorization
    Inte = avaiNum/parentLevPara - 1;
    // cout<< Inte<< endl;
    for(int i = 1; i <= Inte; ++i) {
        if(Inte%i == 0) {
            branPara[level] = i;
            if(i == Inte) {  // we got a new set
                //b[Index][level] = i;
                //d[Index] = level;   // record the depth
                //b[Index][level + 1] = 0;  // mark the end
                //cout<< "The "<< Index<< "-th structure has branching parameters {";
//                for(int j = 0; j < level; ++j) {
//                    //cout<< branPara[j]<< " ";
//                    fout[totNum - testQBegin]<< branPara[j]<< " ";  // output to a file
//                }
//                //cout<< branPara[level]<< "}"<< endl;
//                fout[totNum - testQBegin]<< branPara[level]<< endl;
                Index++;
            }
            else {
//                for(int j = Index; j < NN; ++j) {   // cover all parent nodes
//                    b[j][level] = i;
//                }
                Factorize(totNum, Inte, i, level + 1);
            }
        }
    }
}

