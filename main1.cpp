//
//  main.cpp
//  tranInf_depSearch
//
//  Created by YZhan on 5/27/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <math.h>
using namespace std;
#define TestN 1

int main(int argc, char *argv[]) {
    // parameter initialization
    int TestNep = TestN;    // # of tested single-photon loss
    int TestNcoh = TestN;   // # of tested coherence time
    double ep2D[TestN] = {0};  // single-qubit total loss rate
    double cohT2D[TestN] = {0};   // (coherence time)/(bottom-layer time interval)
    double epStep = (-0.5 - (-4))/(TestN - 1);
    double cohStep = (10.0 - 1.0)/(TestN - 1);
    for(int i = 0; i < TestN; ++i) {
        ep2D[i] = pow(10, -4 + epStep*i);
        cohT2D[i] = pow(10, 1 + cohStep*i);
    }
    ep2D[0] = 0.4;
    cohT2D[0] = 5000;
//    cin>> ep2D[0]>> cohT2D[0];
    double optimalInF[TestN][TestN] = {0};    // optimal result for different ep and cohTime
    int optimalQ[TestN][TestN] = {0};   // corresponding Q
    int optimalD[TestN][TestN] = {0};   // corresponding depth
    int optimalBran[TestN][TestN][10] = {0};  // corresponding branching parameter
    for(int i = 0; i < TestN; ++i)
        for(int j = 0; j < TestN; ++j)
            optimalInF[i][j] = 1;

    // explore possible tree structures for depth=3, 4, 5
    cout<< "The optimal results of total transmission infidelity and corresponding trees, when testing only structures with depth of 3, 4, 5:"<< endl;
    for(int epInd = 0; epInd < TestNep; ++epInd) {
        for(int cohInd = 0; cohInd < TestNcoh; ++cohInd) {
            int testWidth = 50; // test range for each level
            int branPara[10] = {0}; // branching parameters
            double R[10] = {0}; // success prob. of indirect Z-measurement in i-th level (i<depth)
            for(int i = 0; i < testWidth; ++i) {    // level-5, 0 means depth<=4
                branPara[5] = i;    // b_4
                R[5] = 0;
                for(int j = 0; j < testWidth; ++j) {    // level-4, 0 means depth<=3
                    branPara[4] = j;
                    if((i > 0) && (j == 0)) { // b_4!=0, so b_3 cannot be 0
                        continue;;
                    }
                    R[4] = 1 - pow(1 - (1 - ep2D[epInd])*pow(1 - ep2D[epInd] + ep2D[epInd]*R[6], branPara[6]), branPara[5]);
                    for(int k = 1; k < testWidth; ++k) {    // level-3
                        branPara[3] = k;
                        R[3] = 1 - pow(1 - (1 - ep2D[epInd])*pow(1 - ep2D[epInd] + ep2D[epInd]*R[5], branPara[5]), branPara[4]);
                        for(int l = 1; l < testWidth; ++l) {    // level-2
                            branPara[2] = l;
                            R[2] = 1 - pow(1 - (1 - ep2D[epInd])*pow(1 - ep2D[epInd] + ep2D[epInd]*R[4], branPara[4]), branPara[3]);
                            for (int m = 1; m < testWidth; ++m) {   // level-1
                                branPara[1] = m;
                                R[1] = 1 - pow(1 - (1 - ep2D[epInd])*pow(1 - ep2D[epInd] + ep2D[epInd]*R[3], branPara[3]), branPara[2]);
                                
                                // effective loss rate
                                double effLossRate = 0;
                                effLossRate = 1 - (pow(1 - ep2D[epInd] + ep2D[epInd]*R[1], branPara[1]) - pow(ep2D[epInd]*R[1], branPara[1]))*pow(1 - ep2D[epInd] + ep2D[epInd]*R[2], branPara[2]);
                                double probSuccTree = 0;    // success probability from tree encoding
                                probSuccTree = 1 - effLossRate;
                                
                                // tree information
                                int depth = 0;  // depth of the tree
                                if(j == 0) {
                                    depth = 3;
                                }
                                else if(i == 0) {
                                    depth = 4;
                                }
                                else {
                                    depth = 5;
                                }
                                int totNumQ = 1;    // total # of qubits in the tree
                                int nodeNum[10] = {0};   // number of nodes in each level
                                nodeNum[0] = 1;
                                for(int dd = 1; dd <= depth; ++dd) {
                                    nodeNum[dd] = nodeNum[dd - 1]*branPara[dd];
                                    totNumQ += nodeNum[dd];
                                }
                                
                                // infidelity due to limited coherence time
                                int delT[10] = {0};  // time interval in level-j
                                double timeCZ[10] = {0}; // time for a CZ gate in level-j, j = 1~depth
                                double CZTime = 0;  // time for all CZ gates
                                delT[depth] = 1;  // bottom layer
                                timeCZ[depth] = (branPara[depth] - 1)*delT[depth] + 1;
                                delT[depth - 1] = timeCZ[depth] + 1;  // second bottom most layer, it is different from upper layers because photons in bottom layer are not exactly equally distributed
                                timeCZ[depth - 1] = (branPara[depth - 1] - 1)*delT[depth - 1] + 1;
                                for(int dd = depth - 2; dd > 0; --dd) {  // recursive calculation
                                    delT[dd] = timeCZ[dd + 1] + delT[dd + 1] - 1;
                                    timeCZ[dd] = (branPara[dd] - 1)*delT[dd] + 1;
                                }
                                for(int dd = 1; dd <= depth; ++dd) {
                                    CZTime += timeCZ[dd]*nodeNum[dd - 1];
                                }
                                
                                double ETime = 0;   // time for all E gates
                                ETime = totNumQ*1;    // take the time for a E gate (bottom-layer time interval) as the unit
                                double totTime = 0; // total time to generate the best tree
                                double probSuccCoh = 0; // success probability due to limited coherence time
                                totTime = ETime + CZTime;
                                probSuccCoh = exp(-totTime/cohT2D[cohInd]);
                                double totInfidSort = 0;
                                totInfidSort = 1 - probSuccCoh*probSuccTree;
                                // totInfidSort = effLossRate;
                                if(totInfidSort < optimalInF[epInd][cohInd]) {    // update the minimum infidelity
                                    optimalInF[epInd][cohInd] = totInfidSort;
                                    optimalQ[epInd][cohInd] = totNumQ;
                                    optimalD[epInd][cohInd] = depth;
                                    for(int dd = 1; dd <= depth; ++dd) {
                                        optimalBran[epInd][cohInd][dd] = branPara[dd];
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            // output the result
            cout<< "For single-photon loss probability = "<< ep2D[epInd]<< ", and coherence time = "<< cohT2D[cohInd]<< ", the minimum total infidelity = "<< optimalInF[epInd][cohInd]<< ", with a tree whose size Q = "<< optimalQ[epInd][cohInd]<< ", depth = "<< optimalD[epInd][cohInd]<< ", and branching parameters = {";
            for(int dd = 1; dd < optimalD[epInd][cohInd]; ++dd) {
                cout<< optimalBran[epInd][cohInd][dd]<< " ";
            }
            cout<< optimalBran[epInd][cohInd][optimalD[epInd][cohInd]]<< "};"<< endl;
        }
    }
    
    // output for plot
    for(int epInd = 0; epInd < TestNep; ++epInd) {
        for(int cohInd = 0; cohInd < TestNcoh - 1; ++cohInd) {
            cout<< optimalInF[epInd][cohInd]<< " ";
        }
        cout<< optimalInF[epInd][TestNcoh - 1]<< endl;
    }
    
    return 0;
}
