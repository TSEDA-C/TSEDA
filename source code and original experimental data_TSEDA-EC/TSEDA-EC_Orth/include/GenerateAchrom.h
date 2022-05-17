//
// Created by qqq on 2021/9/15.
//

#ifndef FRAME_GENERATEACHROM_H
#include "common.h"
#define FRAME_GENERATEACHROM_H

void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime);
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
chromosome GnrChr_HEFT(vector<double> Rank_b);
chromosome GnrChr_HMEC(vector<double> Rank_b);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
double GnrML_Evl_EFT(chromosome& ch);
double GnrML_Evl_MEC(chromosome& ch);
double DcdEvl(chromosome& ch, bool IsFrw);
void SelectRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
double CalculateEnergy(chromosome& Chrom);
double CalculateECByDelta2(chromosome& Chrom, int& NumBer);
double IFBD(chromosome& ch);
void LBCA(chromosome& chrom);
void InitProModelOfResAlc(vector<vector<double>>& PMR);
void InitProModelOfTskSch(vector<vector<double>>& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
chromosome GnrTskLstOfChr(vector<vector<double> >& PMS);
void UpdatePMR(vector<vector<double>>& PMR, vector<chromosome>& Pop);
void UpdatePMS(vector<vector<double>>& PMS, vector<chromosome>& Pop);

#endif //FRAME_GENERATEACHROM_H
