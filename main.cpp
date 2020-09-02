//
//  main.cpp
//  tranInf_depSearch
//
//  Created by YZhan on 5/27/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//
//  This version is modified to calculate the total error probability in presence of CZ gate error, for fixed single-photon loss and coherence time. The CZ gate error information is input from the results calculated using Mathematica.

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
#define TestN 31   // 400 different photon bandwidth to test, 400=2/0.005

int main(int argc, char *argv[]) {
    // parameter initialization
    double ep2D = 0.1;
    double cohT2D = pow(10, 5); // t_coh*gamma_R
    double optimalInF[TestN] = {0};    // optimal result for different photon bw
    int optimalQ[TestN] = {0};   // corresponding Q
    int optimalD[TestN] = {0};   // corresponding depth
    int optimalBran[TestN][10] = {0};  // corresponding branching parameter
    for(int i = 0; i < TestN; ++i)
        optimalInF[i] = 1;
    
    double gammaLOverR[TestN] = {0};    // gamma_L/gamma_R
    for(int i = 0; i < TestN; ++i) {
        gammaLOverR[i] = 0.005*i;
    }
    
    double CZError[TestN] = {0};
    string ffilename = "/Users/zhanyuan/Desktop/Fig4Data/CZError.txt";  // input CZ error data from file
    ifstream ffin(ffilename);
    if(!ffin) {
        cout<< "Error when opening file "<< ffilename<< " !"<< endl;
    }
    for(int i = 0; i < TestN; ++i) {
        ffin>> CZError[i];
    }
    ffin.close();
    
    // explore possible tree structures for depth=3, 4, 5
    cout<< "The optimal results of total transmission infidelity and corresponding trees, when testing only structures with depth of 3, 4, 5, and for single-photon loss (CZ error NOT included) = "<< ep2D<< ", and t_coh*gammaR = "<<  cohT2D<< ":"<< endl;
    for(int bwInd = 0; bwInd < TestN; ++bwInd) {
        ep2D = 0.1;
        cohT2D = pow(10, 5);    // reset
        ep2D = 1 - (1 - ep2D)*(1 - CZError[bwInd]); // the new single-photon loss including the CZ error
        cohT2D = cohT2D*gammaLOverR[bwInd]/2;   // R=t_coh/tau=t_coh*gamma_L/2=t_coh*gamma_R*gammaLOverR/2, for tau=2/gamma_L, use 1/gammaR as the unit for all the times, including t_coh and tau
        int testWidth = 15; // test range for each level
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
                R[4] = 1 - pow(1 - (1 - ep2D)*pow(1 - ep2D + ep2D*R[6], branPara[6]), branPara[5]);
                for(int k = 1; k < testWidth; ++k) {    // level-3
                    branPara[3] = k;
                    R[3] = 1 - pow(1 - (1 - ep2D)*pow(1 - ep2D + ep2D*R[5], branPara[5]), branPara[4]);
                    for(int l = 1; l < testWidth; ++l) {    // level-2
                        branPara[2] = l;
                        R[2] = 1 - pow(1 - (1 - ep2D)*pow(1 - ep2D + ep2D*R[4], branPara[4]), branPara[3]);
                        for (int m = 1; m < testWidth; ++m) {   // level-1
                            branPara[1] = m;
                            R[1] = 1 - pow(1 - (1 - ep2D)*pow(1 - ep2D + ep2D*R[3], branPara[3]), branPara[2]);
                            
                            // effective loss rate
                            double effLossRate = 0;
                            effLossRate = 1 - (pow(1 - ep2D + ep2D*R[1], branPara[1]) - pow(ep2D*R[1], branPara[1]))*pow(1 - ep2D + ep2D*R[2], branPara[2]);
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
                            
                            // infidelity due to CZ gate
                            double probSuccCZ = 0;
                            probSuccCZ = pow(1 - CZError[bwInd], double(totNumQ - 1));
                            
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
                            probSuccCoh = exp(-totTime/cohT2D);
                            double totInfidSort = 0;
                            probSuccCZ = 1;
                            totInfidSort = 1 - probSuccCoh*probSuccTree*probSuccCZ;
                            // totInfidSort = effLossRate;
                            if(totInfidSort < optimalInF[bwInd]) {    // update the minimum infidelity
                                optimalInF[bwInd] = totInfidSort;
                                optimalQ[bwInd] = totNumQ;
                                optimalD[bwInd] = depth;
                                for(int dd = 1; dd <= depth; ++dd) {
                                    optimalBran[bwInd][dd] = branPara[dd];
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // output the result
        cout<< "For gamma_L/gamma_R = "<< gammaLOverR[bwInd]<< ", CZ gate error = "<< CZError[bwInd]<< ", the minimum total infidelity = "<< optimalInF[bwInd]<< ", with a tree whose size Q = "<< optimalQ[bwInd]<< ", depth = "<< optimalD[bwInd]<< ", and branching parameters = {";
        for(int dd = 1; dd < optimalD[bwInd]; ++dd) {
            cout<< optimalBran[bwInd][dd]<< " ";
        }
        cout<< optimalBran[bwInd][optimalD[bwInd]]<< "};"<< endl;
    }
    
    
    // output for plot
    cout<< endl<< "Minimum total error probability for plot:"<< endl;
    for(int bwInd = 0; bwInd < TestN; ++bwInd) {
        cout<< optimalInF[bwInd]<< " ";
    }
    cout<< endl;
    
    return 0;
}
