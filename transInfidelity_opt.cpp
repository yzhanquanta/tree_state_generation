//
//  main.cpp
//  transInfidelity_opt
//
//  Created by YZhan on 5/13/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#define Q 52    // max test total number of qubits = Q - 1
#define NN 90000 // estimation of number of possibilities
#define TestN 15    // 2D precision

int b[NN][Q] = {0}; // branching parameters b[total number][i-th possibility][branching para]
int d[NN] = {0};   // depth of tree d[total number][i-th possibility]
int Index = 1;  // number of possible tree structures for a certain Q
void Factorize(int avaiNum, int parentLevPara, int level);

int main() {
    // all possible branching parameters for Q=7-101
    double ep2D[TestN] = {0};  // single-qubit total loss rate
    double cohT2D[TestN] = {0};   // (coherence time)/(bottom-layer time interval)
    double epStep = (-0.5 - (-4))/(TestN - 1);
    double cohStep = (10.0 - 1.0)/(TestN - 1);
//        cout<< epStep<< " "<< cohStep;
    for(int i = 0; i < TestN; ++i) {
        ep2D[i] = pow(10, -4 + epStep*i);
        cohT2D[i] = pow(10, 1 + cohStep*i);
    }
    int TestNep = 1;
    int TestNcoh = 1;
    ep2D[0] = 0.15;
    cohT2D[0] = 5000;
//    ep2D[1] = 0.15;
//    ep2D[2] = 0.2;
//    cohT2D[0] = 1000;
//    cohT2D[1] = 5000;
//    cohT2D[2] = 10000;
    double optimalInF[TestN][TestN] = {0};    // optimal result for different ep and cohTime
    int optimalQ[TestN][TestN] = {0};
    double bestTotInfid[Q][TestN][TestN] = {0}; // total best infidelity for each Q
    
    int testQBegin = 7;
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        Index = 1;  // reset # of possible structures
        // factorize to find all possibilities
        int rootPara = 1;   // only one root
        Factorize(testQ, rootPara, 1);
        Index--;
        cout<< "Number of possible branching parameter sets is "<< Index<< " for a total number of qubits Q="<< testQ<< endl;
        
        // for each set of parameters, find the best choice of Q and corresponding infidelity
        for(int epInd = 0; epInd < TestNep; ++epInd) {
            for(int cohInd = 0; cohInd < TestNcoh; ++cohInd) {
                // calculate total infidelity for the whole set with testQ qubits
                double totInfidSort[NN] = {0};  // record all total infidelities for later sorting
                int oriIndex[NN] = {0};   // keep track of original index
                for(int i = 1; i <= Index; ++i) {
                    // calculate effective loss rate for the i-th tree
                    oriIndex[i] = i;
                    double effLossRate = 0;    // effective loss rates
                    double R[Q + 1] = {0};  // success prob. of indirect Z-measurement in i-th level
                    for(int j = d[i] - 1; j > 0; --j) {
                        R[j] = 1 - pow(1 - (1 - ep2D[epInd])*pow(1 - ep2D[epInd] + ep2D[epInd]*R[j + 2], b[i][j + 2]), b[i][j + 1]);
                    }
                    effLossRate = 1 - (pow(1 - ep2D[epInd] + ep2D[epInd]*R[1], b[i][1]) - pow(ep2D[epInd]*R[1], b[i][1]))*pow(1 - ep2D[epInd] + ep2D[epInd]*R[2], b[i][2]);
                    
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
                    
                    // calculate total infidelity
                    totTime = ETime + CZTime;
                    probSuccCoh = exp(-totTime/cohT2D[cohInd]);
                    totInfidSort[i] = 1 - probSuccCoh*probSuccTree;
                }
                
                // sort total infidelities
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
                bestTotInfid[testQ][epInd][cohInd] = totInfidSort[1];
            }
        }
    }
    
    // for each set of parameters, find the best choice of Q and corresponding infidelity
    for(int epInd = 0; epInd < TestNep; ++epInd) {
        for(int cohInd = 0; cohInd < TestNcoh; ++cohInd) {
//            // output to a txt file
//            string fileName, fn1, fn21, fn22, fn3;
//            fn1 = "/Users/zhanyuan/Desktop/data/totInfidInd";
//            fn21 = to_string(epInd);
//            fn22 = to_string(cohInd);
//            fn3 = "_opt.txt";
//            fileName = fn1 + fn21 + fn22 + fn3;
//            ofstream fout(fileName);
//            if(!fout) {
//                cout<< "Error when opening file "<< fileName<< " !"<< endl;
//                return 1;
//            }
//            for(int testQ = testQBegin; testQ < Q - 1; ++testQ) {
//                fout<< bestTotInfid[testQ][epInd][cohInd]<< " ";
//            }
//            fout<< bestTotInfid[Q - 1][epInd][cohInd];
//            fout.close();
            
            // find the best best result
            cout<< "For total loss rate="<< ep2D[epInd]<< ", coherence time="<< cohT2D[cohInd]<< ", the total infidelity for each Q= ";
            double bestBestTotInF = 1;    // the optimal result for all Q
            int bestBestQ = 0;  // corresponding Q for the optimal result
            for(int testQ = testQBegin; testQ < Q; ++testQ) {
                //        cout<< "For Q="<< testQ<< ", the best total infidelity="<< bestTotInfid[testQ]<< ", with a tree whose branching parameters are {";
                //        for(int i = 1; i < bestDep[testQ]; ++i) {
                //            cout<< bestBranPara[testQ][i]<< " ";
                //        }
                //        cout<< bestBranPara[testQ][bestDep[testQ]]<< "};"<< endl;
                if(bestTotInfid[testQ][epInd][cohInd] < bestBestTotInF) {
                    bestBestTotInF = bestTotInfid[testQ][epInd][cohInd];
                    bestBestQ = testQ;
                }
                cout<< bestTotInfid[testQ][epInd][cohInd]<< " ";
            }
            cout<< endl;
            
            // output for plotting
            optimalInF[epInd][cohInd] = bestBestTotInF; // record
            optimalQ[epInd][cohInd] = bestBestQ;
            cout<< "For total loss rate="<< ep2D[epInd]<< ", coherence time="<< cohT2D[cohInd]<< ", the best best total infidelity="<< bestBestTotInF<< ", with Q="<< bestBestQ<< "."<< endl;
        }
    }
    
    // output to a txt file
    string fileName, fn1, fn2, fn3;
    fn1 = "/Users/zhanyuan/Desktop/data/totInfidSCTest";
    fn2 = to_string(Q);
    fn3 = ".txt";
    fileName = fn1 + fn2 + fn3;
    ofstream fout(fileName);
    if(!fout) {
        cout<< "Error when opening file "<< fileName<< " !"<< endl;
        return 1;
    }
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        fout<< bestTotInfid[testQ][0][0]<< " ";
    }
    fout.close();
//
//    // output to a txt file
//    string fileNameQ, fn1Q, fn2Q, fn3Q;
//    fn1Q = "/Users/zhanyuan/Desktop/data/totInfid2DSLQ";
//    fn2Q = to_string(Q);
//    fn3Q = ".txt";
//    fileNameQ = fn1Q + fn2Q + fn3Q;
//    ofstream foutQ(fileNameQ);
//    if(!fout) {
//        cout<< "Error when opening file "<< fileNameQ<< " !"<< endl;
//        return 1;
//    }
//    for(int epInd = 0; epInd < TestNep; ++epInd) {
//        for(int cohInd = 0; cohInd < TestN - 1; ++cohInd) {
//            foutQ<< optimalQ[epInd][cohInd]<< " ";
//        }
//        foutQ<< optimalQ[epInd][TestN - 1]<< endl;
//    }
//    foutQ.close();
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

