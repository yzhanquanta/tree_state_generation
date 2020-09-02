//
//  main.cpp
//  tranInf_LS_partIII
//
//  Created by YZhan on 5/18/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;
#define Q 190    // max test total number of qubits = Q - 1
int testQBegin = 7;

int main() {
    // all possible branching parameters for Q=7-101
    double ep = 0.1;  // single-qubit total loss rate
    double cohT = pow(10, 5.5);   // (coherence time)/(bottom-layer time interval)
    
    double bestTotInfid[Q] = {0}; // total best infidelity for each Q
    for(int i = 0; i < Q; ++i)
        bestTotInfid[i] = 1;
    
    // read the number of possibilities for each Q
    int possibility[Q] = {0};    // number of possibilities for all Q
    string ffilename = "/Users/zhanyuan/Desktop/treeBranPara/numPossibility.txt";
    ifstream ffin(ffilename);
    for(int i = 0; i < Q - testQBegin; ++i) {
        ffin>> possibility[i];
    }
    ffin.close();
    
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        int Index = possibility[testQ - testQBegin];  // # of possible structures
        cout<< "Number of possible branching parameter sets is "<< Index<< " for a total number of qubits Q="<< testQ<< endl;
        // open corresponding branching parameter file
        string filename, fn1, fn2, fn3;
        fn1 = "/Users/zhanyuan/Desktop/treeBranPara/Q";
        fn2 = to_string(testQ);
        fn3 = ".txt";
        filename = fn1 + fn2 + fn3;
        ifstream fin(filename);
        
        for(int i = 0; i < Index; ++i) {
            // read in the branching parameters of the i-th possibility for testQ
            string paraLine;
            getline(fin, paraLine); // each line represents a set of branching parameters
            istringstream iss(paraLine);
            int b[Q] = {0};
            int inputBP = 0, inputBPInd = 0;
            while(iss>> inputBP) {
                b[inputBPInd] = inputBP;
                //cout<< b[inputBPInd]<< " ";
                inputBPInd++;
            }
            //cout<< endl;
            inputBPInd--;   // depth of this tree
            
            // calculate effective loss rate for the i-th tree
            double effLossRate = 0;    // effective loss rates
            double R[Q + 1] = {0};  // success prob. of indirect Z-measurement in i-th level
            for(int j = inputBPInd - 1; j > 0; --j) {
                R[j] = 1 - pow(1 - (1 - ep)*pow(1 - ep + ep*R[j + 2], b[j + 2]), b[j + 1]);
            }
            effLossRate = 1 - (pow(1 - ep + ep*R[1], b[1]) - pow(ep*R[1], b[1]))*pow(1 - ep + ep*R[2], b[2]);
            
            // calculate infidelity due to limited coherence time
            double totTime = 0; // total time to generate the best tree
            double ETime = 0;   // time for all E gates
            double CZTime = 0;  // time for all CZ gates
            double probSuccCoh = 0; // success probability due to limited coherence time
            double probSuccTree = 0;    // success probability from tree encoding
            probSuccTree = 1 - effLossRate;
            ETime = testQ*1;    // take the time for a E gate (bottom-layer time interval) as the unit
            int delT[Q] = {0};  // time interval in level-j
            double timeCZ[Q] = {0}; // time for a CZ gate in level-j, j = 1~depth
            delT[inputBPInd] = 1;  // bottom layer
            timeCZ[inputBPInd] = (b[inputBPInd] - 1)*delT[inputBPInd] + 1;
            delT[inputBPInd - 1] = timeCZ[inputBPInd] + 1;  // second bottom most layer, it is different from upper layers because photons in bottom layer are not exactly equally distributed
            timeCZ[inputBPInd - 1] = (b[inputBPInd - 1] - 1)*delT[inputBPInd - 1] + 1;
            for(int j = inputBPInd - 2; j > 0; --j) {  // recursive calculation
                delT[j] = timeCZ[j + 1] + delT[j + 1] - 1;
                timeCZ[j] = (b[j] - 1)*delT[j] + 1;
            }
            int nodeNum[Q] = {0};   // number of nodes in level-j
            nodeNum[0] = 1;
            for(int j = 1; j <= inputBPInd; ++j) {
                nodeNum[j] = nodeNum[j - 1]*b[j];
            }
            for(int j = 1; j <= inputBPInd; ++j) {
                CZTime += timeCZ[j]*nodeNum[j - 1];
            }
            
            // calculate total infidelity
            totTime = ETime + CZTime;
            probSuccCoh = exp(-totTime/cohT);
            double totInfidSort = 0;
            totInfidSort = 1 - probSuccCoh*probSuccTree;
            if(totInfidSort < bestTotInfid[testQ]) {    // update the minimum infidelity
                bestTotInfid[testQ] = totInfidSort;
            }
        }
        
        
        fin.close();
    }
    
    
    // output to a txt file
    string fileName, fn1, fn2, fn3;
    fn1 = "/Users/zhanyuan/Desktop/data/totInfid2DSG01_1055Test";
    fn2 = to_string(Q);
    fn3 = ".txt";
    fileName = fn1 + fn2 + fn3;
    ofstream fout(fileName);
    if(!fout) {
        cout<< "Error when opening file "<< fileName<< " !"<< endl;
        return 1;
    }
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        cout<< "Q="<< testQ<< ", infidelity="<< bestTotInfid[testQ]<< " ";
        fout<< bestTotInfid[testQ]<< " ";
    }
    fout.close();
    
    
    return 0;
}


