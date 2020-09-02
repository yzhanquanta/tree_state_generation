//
//  main.cpp
//  tranInf_LS_partII_BackUp
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
#define Q 102    // max test total number of qubits = Q - 1
#define TestN 15    // 2D precision
int testQBegin = 7;

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
    //ep2D[0] = 0.1;
    cohT2D[0] = pow(10,5);
    int TestNep = TestN;
    int TestNcoh = 1;
    //    ep2D[1] = 0.15;
    //    ep2D[2] = 0.2;
    //    cohT2D[0] = 1000;
    //    cohT2D[1] = 5000;
    //    cohT2D[2] = 10000;
    double optimalInF[TestN][TestN] = {0};    // optimal result for different ep and cohTime
    int optimalQ[TestN][TestN] = {0};
    double bestTotInfid[Q][TestN][TestN] = {0}; // total best infidelity for each Q
    for(int i = 0; i < Q; ++i)
        for(int j = 0; j < TestN; ++j)
            for(int k = 0; k < TestN; ++k)
                bestTotInfid[i][j][k] = 1;
    int origInd[Q][TestN] = {0};
    
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
            
            // for each set of parameters, find the corresponding infidelity for each set and find the best for testQ
            for(int epInd = 0; epInd < TestNep; ++epInd) {
                for(int cohInd = 0; cohInd < TestNcoh; ++cohInd) {
                    // calculate effective loss rate for the i-th tree
                    double effLossRate = 0;    // effective loss rates
                    double R[Q + 1] = {0};  // success prob. of indirect Z-measurement in i-th level
                    for(int j = inputBPInd - 1; j > 0; --j) {
                        R[j] = 1 - pow(1 - (1 - ep2D[epInd])*pow(1 - ep2D[epInd] + ep2D[epInd]*R[j + 2], b[j + 2]), b[j + 1]);
                    }
                    effLossRate = 1 - (pow(1 - ep2D[epInd] + ep2D[epInd]*R[1], b[1]) - pow(ep2D[epInd]*R[1], b[1]))*pow(1 - ep2D[epInd] + ep2D[epInd]*R[2], b[2]);
                    
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
                    probSuccCoh = exp(-totTime/cohT2D[cohInd]);
                    double totInfidSort = 0;
                    totInfidSort = 1 - probSuccCoh*probSuccTree;
                    if(totInfidSort < bestTotInfid[testQ][epInd][cohInd]) {    // update the minimum infidelity
                        bestTotInfid[testQ][epInd][cohInd] = totInfidSort;
                        origInd[testQ][epInd] = i;
                    }
                }
                
            }
        }
        fin.close();
    }
    
    // for different coherence times, the best tree structure for each Q
    for(int testQ = testQBegin; testQ < Q; ++testQ) {
        cout<< "For Q="<< testQ;
        for(int epInd = 0; epInd < TestNep; ++epInd) {
            cout<< ", when single-photon loss rate="<< ep2D[epInd]<< ", the best tree structure is the"<< origInd[testQ][epInd]<< "-th one";
        }
        cout<< endl;
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
            }
            
            // output for plotting
            optimalInF[epInd][cohInd] = bestBestTotInF; // record
            optimalQ[epInd][cohInd] = bestBestQ;
            cout<< "For total loss rate="<< ep2D[epInd]<< ", coherence time="<< cohT2D[cohInd]<< ", the best best total infidelity="<< bestBestTotInF<< ", with Q="<< bestBestQ<< "."<< endl;
        }
    }
    
//        // output to a txt file
//        string fileName, fn1, fn2, fn3;
//        fn1 = "/Users/zhanyuan/Desktop/data/totInfid2DSLTest3BackUp";
//        fn2 = to_string(Q);
//        fn3 = ".txt";
//        fileName = fn1 + fn2 + fn3;
//        ofstream fout(fileName);
//        if(!fout) {
//            cout<< "Error when opening file "<< fileName<< " !"<< endl;
//            return 1;
//        }
//        for(int epInd = 0; epInd < TestNep; ++epInd) {
//            for(int cohInd = 0; cohInd < TestN - 1; ++cohInd) {
//                fout<< optimalInF[epInd][cohInd]<< " ";
//            }
//            fout<< optimalInF[epInd][TestN - 1]<< endl;
//        }
//        fout.close();
    
//        // output to a txt file
//        string fileNameQ, fn1Q, fn2Q, fn3Q;
//        fn1Q = "/Users/zhanyuan/Desktop/data/totInfid2DSL_BU";
//        fn2Q = to_string(Q);
//        fn3Q = ".txt";
//        fileNameQ = fn1Q + fn2Q + fn3Q;
//        ofstream foutQ(fileNameQ);
//        if(!fout) {
//            cout<< "Error when opening file "<< fileNameQ<< " !"<< endl;
//            return 1;
//        }
//        for(int epInd = 0; epInd < TestNep; ++epInd) {
//            for(int cohInd = 0; cohInd < TestN - 1; ++cohInd) {
//                foutQ<< optimalQ[epInd][cohInd]<< " ";
//            }
//            foutQ<< optimalQ[epInd][TestN - 1]<< endl;
//        }
//        foutQ.close();
    return 0;
}


