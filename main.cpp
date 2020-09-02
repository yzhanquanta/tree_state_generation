//
//  main.cpp
//  tranInf_depSrc_Sort
//
//  Created by YZhan on 5/27/20.
//  Copyright © 2020 YuanZhan. All rights reserved.
//

#include <iostream>
#include <fstream>
using namespace std;
#define num 270

int main() {
    string filename = "/Users/zhanyuan/Desktop/049.txt";
    ifstream fin(filename);
    
    int outQ[num] = {0};
    double outEffLoss[num] = {0};
    double bestSoFar = 1;
    int Index = 0;
    
    for(int i = 0; i < num; ++i) {
        int Q;
        double effLoss;
        fin>> Q>> effLoss;
        if(effLoss < bestSoFar) {
            outQ[Index] = Q;
            outEffLoss[Index] = effLoss;
            bestSoFar = effLoss;
            Index++;
        }
    }
    fin.close();
    
    cout<< "There are "<< Index<< " trees worth outputing:"<< endl;
    for(int i = 0; i < Index; ++i) {
        cout<< outQ[i]<< " ";
    }
    cout<< endl<< endl;
    for(int i = 0; i < Index; ++i) {
        cout<< outEffLoss[i]<< " ";
    }
    cout<< endl;
    return 0;
}
