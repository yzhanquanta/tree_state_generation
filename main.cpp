//
//  main.cpp
//  SC_maxTest2
//
//  Created by YZhan on 5/24/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;
#define TestNum 400
int testQBegin = 7; // starting test point
int testQEnd = 325;   // ending test point
#define Q 1000    // max test total number of qubits = Q - 1
int branPara[Q] = {0};  // record current set of branching parameters for calculation
long long Index = 1;  // number of possible tree structures for a certain Q
void Factorize(int totNum, int avaiNum, int parentLevPara, int level);
#define TestN 15    // 2D precision
int TestNep = 1;
int TestNcoh = 1;
double ep2D[TestN] = {0};  // single-qubit total loss rate
double cohT2D[TestN] = {0};   // (coherence time)/(bottom-layer time interval)
double bestTotInfid[TestNum][TestN][TestN] = {0}; // total best infidelity for each Q
int bestQ[TestN][TestN] = {0};  // best choice of Q for each set of paras
int bestBranPara[Q] = {0};
int bestDepth = 0;

int main(int argc, char *argv[]) {
    // parameter initialization
//    testQBegin = atoi(argv[1]);
//    testQEnd = atoi(argv[2]);
    double epStep = (-0.5 - (-4))/(TestN - 1);
    double cohStep = (10.0 - 1.0)/(TestN - 1);
    //        cout<< epStep<< " "<< cohStep;
    for(int i = 0; i < TestN; ++i) {
        ep2D[i] = pow(10, -4 + epStep*i);
        cohT2D[i] = pow(10, 1 + cohStep*i);
    }
    ep2D[0] = 0.1;
    cohT2D[0] = pow(10, 9);
    //    ep2D[1] = 0.15;
    //    ep2D[2] = 0.2;
    //    cohT2D[0] = 1000;
    //    cohT2D[1] = 5000;
    //    cohT2D[2] = 10000;
    double optimalInF[TestN][TestN] = {0};    // optimal result for different ep and cohTime
    int optimalQ[TestN][TestN] = {0};   // corresponding Q
    
    for(int i = 0; i < TestN; ++i)
        for(int j = 0; j < TestN; ++j)
            optimalInF[i][j] = 1;

    for(int i = 0; i < TestNum; ++i)
        for(int j = 0; j < TestN; ++j)
            for(int k = 0; k < TestN; ++k)
                bestTotInfid[i][j][k] = 1;
    
    // test each Q and keep updating the best choice of tree structure for given parameters
    for(int testQ = testQBegin; testQ < testQEnd; ++testQ) {
        Index = 1;  // reset # of possible structures
        bestDepth = 1;
        // factorize to find all possibilities
        int rootPara = 1;   // only one root
        branPara[0] = 1;    // level 0 has one photon
        Factorize(testQ, testQ, rootPara, 1);
        Index--;
        //cout<< "Number of possible branching parameter sets is "<< Index<< " for a total number of qubits Q="<< testQ<< ", "<< endl;
        // update the optimal choice of Q and the optimal infidelity
        for(int epInd = 0; epInd < TestNep; ++epInd) {
            for(int cohInd = 0; cohInd < TestNcoh; ++cohInd) {
                //cout<< "with infidelity: "<< bestTotInfid[testQ - testQBegin][epInd][cohInd]<< endl;
                if(bestTotInfid[testQ - testQBegin][epInd][cohInd] < optimalInF[epInd][cohInd]) {
                    optimalInF[epInd][cohInd] = bestTotInfid[testQ - testQBegin][epInd][cohInd];
                    optimalQ[epInd][cohInd] = testQ;
                    cout<< "Update: Q="<< testQ<< ", and the best tree sturcture out of "<< Index<< " choices is:{";
                    for(int k = 0; k < bestDepth; ++k) {
                        cout<< bestBranPara[k]<< " ";
                    }
                    cout<< bestBranPara[bestDepth]<< "}, with depth= "<< bestDepth<< " and total infidelity= "<< bestTotInfid[testQ - testQBegin][epInd][cohInd]<< endl;
                }
            }
        }
    }
    
    for(int epInd = 0; epInd < TestNep; ++epInd) {
        for(int cohInd = 0; cohInd < TestNcoh; ++cohInd) {
            cout<< "For single-photon loss rate="<< ep2D[epInd]<< " and coherence time="<< cohT2D[cohInd]<< ", the best transmission infidelity for different Q's are: ";
            for(int testQ = testQBegin; testQ < testQEnd; ++testQ) {
                cout<< bestTotInfid[testQ - testQBegin][epInd][cohInd]<< " ";
            }
            cout<< endl<< "The optimal infidelity is: "<< optimalInF[epInd][cohInd]<< ", with Q="<< optimalQ[epInd][cohInd]<< "."<< endl;
        }
    }
    
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
                branPara[level + 1] = 0;    // mind the tail for recursive calculation
                branPara[level + 2] = 0;
                for(int epInd = 0; epInd < TestNep; ++epInd) {
                    for(int cohInd = 0; cohInd < TestNcoh; ++cohInd) {  // for each set of given parameters, calculate infidelity for this new structure
                        // calculate effective loss rate for this new set
                        double effLossRate = 0;    // effective loss rates
                        double R[Q + 1] = {0};  // success prob. of indirect Z-measurement in i-th level
                        for(int j = level - 1; j > 0; --j) {
                            R[j] = 1 - pow(1 - (1 - ep2D[epInd])*pow(1 - ep2D[epInd] + ep2D[epInd]*R[j + 2], branPara[j + 2]), branPara[j + 1]);
                        }
                        effLossRate = 1 - (pow(1 - ep2D[epInd] + ep2D[epInd]*R[1], branPara[1]) - pow(ep2D[epInd]*R[1], branPara[1]))*pow(1 - ep2D[epInd] + ep2D[epInd]*R[2], branPara[2]);
                        
                        // calculate infidelity due to limited coherence time
                        double totTime = 0; // total time to generate the best tree
                        double ETime = 0;   // time for all E gates
                        double CZTime = 0;  // time for all CZ gates
                        double probSuccCoh = 0; // success probability due to limited coherence time
                        double probSuccTree = 0;    // success probability from tree encoding
                        probSuccTree = 1 - effLossRate;
                        ETime = totNum*1;    // take the time for a E gate (bottom-layer time interval) as the unit
                        int delT[Q] = {0};  // time interval in level-j
                        double timeCZ[Q] = {0}; // time for a CZ gate in level-j, j = 1~depth
                        delT[level] = 1;  // bottom layer
                        timeCZ[level] = (branPara[level] - 1)*delT[level] + 1;
                        delT[level - 1] = timeCZ[level] + 1;  // second bottom most layer, it is different from upper layers because photons in bottom layer are not exactly equally distributed
                        timeCZ[level - 1] = (branPara[level - 1] - 1)*delT[level - 1] + 1;
                        for(int j = level - 2; j > 0; --j) {  // recursive calculation
                            delT[j] = timeCZ[j + 1] + delT[j + 1] - 1;
                            timeCZ[j] = (branPara[j] - 1)*delT[j] + 1;
                        }
                        int nodeNum[Q] = {0};   // number of nodes in level-j
                        nodeNum[0] = 1;
                        for(int j = 1; j <= level; ++j) {
                            nodeNum[j] = nodeNum[j - 1]*branPara[j];
                        }
                        for(int j = 1; j <= level; ++j) {
                            CZTime += timeCZ[j]*nodeNum[j - 1];
                        }
                        
                        // calculate total infidelity
                        totTime = ETime + CZTime;
                        probSuccCoh = exp(-totTime/cohT2D[cohInd]);
                        double totInfidSort = 0;
                        totInfidSort = 1 - probSuccCoh*probSuccTree;
                        if(totInfidSort < bestTotInfid[totNum - testQBegin][epInd][cohInd]) {    // update the minimum infidelity
                            bestTotInfid[totNum - testQBegin][epInd][cohInd] = totInfidSort;
                            for(int k = 0; k <= level; ++k) {
                                bestBranPara[k] = branPara[k];
                                bestDepth = level;
                            }
                        }
                    }
                }
                Index++;
            }
            else {
                Factorize(totNum, Inte, i, level + 1);
            }
        }
    }
}

