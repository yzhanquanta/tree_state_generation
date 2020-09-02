//
//  main.cpp
//  incCZ
//
//  Created by YZhan on 6/7/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
#define TestN 101   // 400 different photon bandwidth to test, 400=2/0.005

int main(int argc, char *argv[]) {
    // parameter initialization
    double ep2Dori[15] = {0};  // single-qubit total loss rate
    double cohT2Dori[18] = {0};   // (coherence time)/(bottom-layer time interval)
    double epStep = (-0.5 - (-4))/(15 - 1);
    double cohStep = (10.0 - 1.0)/(15 - 1);
    for(int i = 0; i < 15; ++i) {
        ep2Dori[i] = pow(10, -4 + epStep*i);
        cohT2Dori[i + 3] = pow(10, 1 + cohStep*i);
    }
    cohT2Dori[0] = pow(10, 1 - 9.0/14.0*3.0);
    cohT2Dori[1] = pow(10, 1 - 9.0/14.0*2.0);
    cohT2Dori[2] = pow(10, 1 - 9.0/14.0*1.0);
    
    double optimalInF[15][15] = {0};    // optimal result for different photon bw
//    int optimalQ[15][15] = {0};   // corresponding Q
//    int optimalD[15][15] = {0};   // corresponding depth
//    int optimalBran[15][15][10] = {0};  // corresponding branching parameter
    for(int i = 0; i < 15; ++i)
        for(int j = 0; j < 15; ++j)
            optimalInF[i][j] = 1;
    
    // input CZ error data from file
    double gammaLOverR[TestN] = {0};    // gamma_L/gamma_R
    for(int i = 0; i < TestN; ++i) {
        gammaLOverR[i] = 0.005*i;
    }
    double CZError[TestN] = {0};
    string ffilename = "/Users/zhanyuan/Desktop/Fig4Data/CZError.txt";
    ifstream ffin(ffilename);
    if(!ffin) {
        cout<< "Error when opening file "<< ffilename<< " !"<< endl;
    }
    for(int i = 0; i < TestN; ++i) {
        ffin>> CZError[i];
    }
    ffin.close();
    
    // input total error probabilties for give mu and R
    double totErrProb[15][18] = {0};
    string fffilename = "/Users/zhanyuan/Desktop/depSrc/depSrc_data.txt";
    ifstream fffin(fffilename);
    if(!fffin) {
        cout<< "Error when opening file "<< fffilename<< " !"<< endl;
    }
    cout<< "The original minimum total error probability:"<< endl;
    for(int i = 0; i < 15; ++i) {
        for(int j = 0; j < 18; ++j) {
            fffin>> totErrProb[i][j];
            cout<< totErrProb[i][j]<< " ";
        }
        cout<< endl;
    }
    fffin.close();
    
    // explore possible tree structures for depth=3, 4, 5
    for(int epInd = 0; epInd < 15; ++epInd) {
        for(int cohInd = 0; cohInd < 15; ++cohInd) {
            // cout<< "The optimal results of total transmission infidelity and corresponding trees, when testing only structures with depth of 3, 4, 5, and for single-photon loss (CZ error NOT included) = "<< ep2Dori[epInd]<< ", and t_coh*gammaR = "<<  cohT2Dori[cohInd]<< ":"<< endl;
            for(int bwInd = 0; bwInd < TestN; ++bwInd) {
                double ep2D = ep2Dori[epInd];
                double cohT2D = cohT2Dori[cohInd + 3];    // reset
                ep2D = 1 - (1 - ep2D)*(1 - CZError[bwInd]); // the new single-photon loss including the CZ error
                cohT2D = cohT2D*gammaLOverR[bwInd]/2;   // R=t_coh/tau=t_coh*gamma_L/2=t_coh*gamma_R*gammaLOverR/2, for tau=2/gamma_L, use 1/gammaR as the unit for all the times, including t_coh and tau
                
                // search depSrc_data for the closest result
                int closestEpInd = 14;
                int closestCohInd = 17;
                for(int closestEpSrc = 0; closestEpSrc < 15; ++closestEpSrc) {  // search ep
                    if(ep2Dori[closestEpSrc] > ep2D) {
                        if(fabs(ep2Dori[closestEpSrc] - ep2D) < fabs(ep2Dori[closestEpSrc - 1] - ep2D)) {
                            closestEpInd = closestEpSrc;
                        }
                        else {
                            closestEpInd = closestEpSrc - 1;
                        }
                        break;
                    }
                }
                for(int closestCohSrc = 0; closestCohSrc < 18; ++closestCohSrc) {  // search cohT
                    if(cohT2Dori[0] > cohT2D) {
                        closestCohInd = 0;
                        break;
                    }
                    if(cohT2Dori[closestCohSrc] > cohT2D) {
                        if(fabs(cohT2Dori[closestCohSrc] - cohT2D) < fabs(cohT2Dori[closestCohSrc - 1] - cohT2D)) {
                            closestCohInd = closestCohSrc;
                        }
                        else {
                            closestCohInd = closestCohSrc - 1;
                        }
                        break;
                    }
                }
                // cout<< closestEpInd<< " "<< closestCohInd<< endl;
                // update the optimalInf
                if(totErrProb[closestEpInd][closestCohInd] < optimalInF[epInd][cohInd]) {
                    optimalInF[epInd][cohInd] = totErrProb[closestEpInd][closestCohInd];
                }
            }
        }
    }
    
    // output for plot
    cout<< endl<< "Minimum total error probability for plot:"<< endl;
    for(int epInd = 0; epInd < 15; ++epInd) {
        for(int cohInd = 0; cohInd < 14; ++cohInd) {
            cout<< optimalInF[epInd][cohInd]<< " ";
        }
        cout<< optimalInF[epInd][14]<< endl;
    }
    
    return 0;
}
