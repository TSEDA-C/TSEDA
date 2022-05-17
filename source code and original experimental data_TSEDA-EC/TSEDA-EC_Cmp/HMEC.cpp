//
// Created by qqq on 2021/12/13.
//

#include "common.h"
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"

double runHMEC(string XmlFile, string RscAlcFile, double& SchTime) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    CalculateLevelList();
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b,cc,ww);
    chromosome Chrom_HMEC_b = GnrChr_HMEC(Rank_b);
    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    return Chrom_HMEC_b.EnergyConsumption;
}