//
//  main.cpp
//  effLoss_depSrch
//
//  Created by YZhan on 5/27/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <math.h>
using namespace std;
#define TestNum 10000

int main() {
    // parameter initialization
    double ep = 0.49;
    double effLoss[TestNum] = {0};
//    for(int i = 0; i < TestNum; ++i) {
//        effLoss[i] = 1;
//    }
    int minQ[TestNum] = {0};
    // int bestBran[TestNum][10] = {0};
    int Index = 0;
    int mark = 0;
    int flag = 0;
    double soFarBest = 1;
    
    // explore possible tree structures for depth=3, 4, 5
    int testWidth = 100; // test range for each level
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
            R[4] = 1 - pow(1 - (1 - ep)*pow(1 - ep + ep*R[6], branPara[6]), branPara[5]);
            for(int k = 1; k < testWidth; ++k) {    // level-3
                branPara[3] = k;
                R[3] = 1 - pow(1 - (1 - ep)*pow(1 - ep + ep*R[5], branPara[5]), branPara[4]);
                for(int l = 1; l < testWidth; ++l) {    // level-2
                    branPara[2] = l;
                    R[2] = 1 - pow(1 - (1 - ep)*pow(1 - ep + ep*R[4], branPara[4]), branPara[3]);
                    for (int m = 1; m < testWidth; ++m) {   // level-1
                        branPara[1] = m;
                        R[1] = 1 - pow(1 - (1 - ep)*pow(1 - ep + ep*R[3], branPara[3]), branPara[2]);
                        
                        // effective loss rate
                        double effLossRate = 0;
                        effLossRate = 1 - (pow(1 - ep + ep*R[1], branPara[1]) - pow(ep*R[1], branPara[1]))*pow(1 - ep + ep*R[2], branPara[2]);
                        
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
                        
                        // record good choices;
                        if(effLossRate < soFarBest) {
                            effLoss[Index] = effLossRate;
                            soFarBest = effLossRate;
                            minQ[Index] = totNumQ;
                            if((effLoss[Index] < pow(10, -8.5)) && (flag == 0)) {
                                mark = Index;
                                flag = 1;
                            }
                            Index++;
                        }
                        
                        
                        
                    }
                }
            }
        }
    }
    mark = Index;
    // output the result
    cout<< "For single-photon loss probability = "<< ep<< ", there are "<< Index<< " good choices of trees as follow:"<< endl;
    cout<< "mark: "<< mark<< endl;
    for(int i = 0; i < mark; ++i) {
        if(effLoss[i] < 1) {
            cout<< minQ[i]<< " ";
        }
    }
    cout<< endl<< ", with effective loss rates:"<< endl;
    for(int i = 0; i < mark; ++i) {
        if(effLoss[i] < 1) {
            cout<< effLoss[i]<< " ";
        }
    }
    cout<< endl;
    cout<< "output to excel"<< endl;
    for(int i = 0; i < mark; ++i) {
        if(effLoss[i] < 1) {
            cout<< minQ[i]<< " "<< effLoss[i]<< endl;
        }
    }
    return 0;
}

