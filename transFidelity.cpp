//
//  main.cpp
//  transFidelity
//
//  Created by YZhan on 5/10/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#define Q 102    // max test total number of qubits = Q - 1
#define NN 90000 // estimation of number of possibilities

int b[NN][Q] = {0}; // branching parameters
int d[NN] = {0};   // depth of tree
int Index = 1;  // number of possible tree structures for a certain Q
void Factorize(int avaiNum, int parentLevPara, int level);

int main() {
    // all possible branching parameters for Q=5-103
    // double currentMin = 0.5;
    // int currentIndex = 0;
    double ep = 0.1;   // single-qubit total loss rate
    double cohTime = pow(10, 5);  // (coherence time)/(bottom-layer time interval)
    double bestBestTotInF = 1;    // the optimal result for all Q
    int bestBestQ = 0;  // corresponding Q for the optimal result
    int bestBestDep = 0;    // corresponding depth
    double bestBestTree[Q] = {0};   // corresponding branching parameters
    double bestTotInfid[Q] = {0}; // total best infidelity for each Q
    int bestDep[Q] = {0}; // corresponding depth
    int bestBranPara[Q][Q] = {0}; // corresponding branching parameters
    int testQBegin = 7;
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        // clear all arrays
        Index = 1;  // reset # of possible structures
        for(int i = 1; i < NN; ++i) {
            if(d[i] == 0) {  // all clear
                break;
            }
            d[i] = 0;
            for(int j = 1; j < Q; ++j) {
                if(b[i][j] == 0) {   // have reached the tail
                    break;
                }
                b[i][j] = 0;
            }
        }
        
        // factorize to find all possibilities
        int rootPara = 1;   // only one root
        Factorize(testQ, rootPara, 1);
        Index--;
        cout<< "Number of possible branching parameter sets is "<< Index<< " for a total number of qubits Q="<< testQ<< endl;
    
        // calculate total infidelity for the whole set with testQ qubits
        double totInfidSort[NN] = {0};  // record all total infidelities for later sorting
        for(int i = 1; i <= Index; ++i) {
            // calculate effective loss rate for the i-th tree
            double effLossRate = 0;    // effective loss rates
            double R[Q + 1] = {0};  // success prob. of indirect Z-measurement in i-th level
            for(int j = d[i] - 1; j > 0; --j) {
                R[j] = 1 - pow(1 - (1 - ep)*pow(1 - ep + ep*R[j + 2], b[i][j + 2]), b[i][j + 1]);
            }
            effLossRate = 1 - (pow(1 - ep + ep*R[1], b[i][1]) - pow(ep*R[1], b[i][1]))*pow(1 - ep + ep*R[2], b[i][2]);
            
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
            delT[d[i]] = 1;  // bottom layer
            timeCZ[d[i]] = (b[i][d[i]] - 1)*delT[d[i]] + 1;
            delT[d[i] - 1] = timeCZ[d[i]] + 1;  // second bottom most layer, it is different from upper layers because photons in bottom layer are not exactly equally distributed
            timeCZ[d[i] - 1] = (b[i][d[i] - 1] - 1)*delT[d[i] - 1] + 1;
            for(int j = d[i] - 2; j > 0; --j) {  // recursive calculation
                delT[j] = timeCZ[j + 1] + delT[j + 1] - 1;
                timeCZ[j] = (b[i][j] - 1)*delT[j] + 1;
            }
            int nodeNum[Q] = {0};   // number of nodes in level-j
            nodeNum[0] = 1;
            for(int j = 1; j <= d[i]; ++j) {
                nodeNum[j] = nodeNum[j - 1]*b[i][j];
            }
            for(int j = 1; j <= d[i]; ++j) {
                CZTime += timeCZ[j]*nodeNum[j - 1];
            }
//            cout<< "CZTime="<< CZTime<< " for branching parameters={";
//            for(int j = 1; j <= d[i]; ++j) {
//                cout<< b[i][j]<< " ";
//            }
//            cout<< "}"<< endl;
            
            // calculate total infidelity
            totTime = ETime + CZTime;
            probSuccCoh = exp(-totTime/cohTime);
            totInfidSort[i] = 1 - probSuccCoh*probSuccTree;
        }
        
        // sort total infidelities
        int oriIndex[NN] = {0};   // keep track of original index
        for(int i = 1; i <= Index; ++i) {
            oriIndex[i] = i;
        }
        for(int i = 1; i < Index; ++i) { // bubble sort
            for(int j = 1; j < Index - i; ++j) {
                if(totInfidSort[j] > totInfidSort[j + 1]) {
                    double temp = totInfidSort[j];
                    totInfidSort[j] = totInfidSort[j + 1];
                    totInfidSort[j + 1] = temp;
                    int tempIndex = oriIndex[j];
                    oriIndex[j] = oriIndex[j + 1];
                    oriIndex[j + 1] = tempIndex;
                }
            }
        }
        
//        // test output for check
//        for(int i = 1; i < Index; ++i) {
//            cout<< "This is the "<< i<< "-th best case, whose total infidelity="<< totInfidSort[i]<< " and the branching parameters are {";
//            for(int j = 1; j <= d[oriIndex[i]]; ++j) {
//                cout<< b[oriIndex[i]][j]<< " ";
//            }
//            cout<< "}"<< endl;
//        }
        
        bestTotInfid[testQ] = totInfidSort[1];
        bestDep[testQ] = d[oriIndex[1]];
        for(int i = 1; i <= bestDep[testQ]; ++i) {
            bestBranPara[testQ][i] = b[oriIndex[1]][i];
        }
    }
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        cout<< "For Q="<< testQ<< ", the best total infidelity="<< bestTotInfid[testQ]<< ", with a tree whose branching parameters are {";
        for(int i = 1; i < bestDep[testQ]; ++i) {
            cout<< bestBranPara[testQ][i]<< " ";
        }
        cout<< bestBranPara[testQ][bestDep[testQ]]<< "};"<< endl;
        // find the best best result
        if(bestTotInfid[testQ] < bestBestTotInF) {
            bestBestTotInF = bestTotInfid[testQ];
            bestBestQ = testQ;
            bestBestDep = bestDep[testQ];
            for(int i = 1; i <= bestBestDep; ++i) {
                bestBestTree[i] = bestBranPara[testQ][i];
            }
        }
    }
    
    // output for plotting
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        cout<< bestTotInfid[testQ]<< " ";
    }
    cout<< endl<< endl;
    cout<< "For total loss rate="<< ep<< ", coherence time="<< cohTime<< ", the best best total infidelity="<< bestBestTotInF<< ", with Q="<< bestBestQ<< ", and branching parameters={";
    for(int i = 1; i < bestBestDep; ++i) {
        cout<< bestBestTree[i]<< " ";
    }
    cout<< bestBestTree[bestBestDep]<< "}."<< endl<< endl;
    
    // output to a txt file
//    string fileName, fn1, fn2, fn3;
//    fn1 = "/Users/zhanyuan/Desktop/data/totInfid4";
//    fn2 = to_string(int(cohTime));
//    fn3 = ".txt";
//    fileName = fn1 + fn2 + fn3;
    string fileName = "/Users/zhanyuan/Desktop/data/totInfidIndSG3.txt";
    ofstream fout(fileName);
    if(!fout) {
        cout<< "Error when opening file "<< fileName<< " !"<< endl;
        return 1;
    }
    for(int testQ = testQBegin; testQ < Q - 1; ++testQ) {
        fout<< bestTotInfid[testQ]<< " ";
    }
    fout<< bestTotInfid[Q - 1];
    fout.close();
    return 0;
}

void Factorize(int avaiNum, int parentLevPara, int level) {
    int Inte;   // the integer that needs factorization
    Inte = avaiNum/parentLevPara - 1;
    // cout<< Inte<< endl;
    for(int i = 1; i <= Inte; ++i) {
        if(Inte%i == 0) {
            if(i == Inte) {  // we got a new set
                b[Index][level] = i;
                d[Index] = level;   // record the depth
                b[Index++][level + 1] = 0;  // mark the end
            }
            else {
                for(int j = Index; j < NN; ++j) {   // cover all parent nodes
                    b[j][level] = i;
                }
                Factorize(Inte, i, level + 1);
            }
        }
    }
}
