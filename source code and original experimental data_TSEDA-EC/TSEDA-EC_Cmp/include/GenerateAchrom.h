//
// Created by qqq on 2021/9/15.
//

#ifndef FRAME_GENERATEACHROM_H
#include "common.h"
#define FRAME_GENERATEACHROM_H

//void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime);
void UpdateITL(set<double>& ITLofRscId, double& StartTime, double& EndTime);
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w);
void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
void ModifyRscAlcLstByCode_RK(chromosome& chrom);
void GnrCode_RK(chromosome& chrom);
void GnrRscAlcTskSchLstFromCode_RK(chromosome& chrom);
vector<int> GnrSS_TS();
vector<int> GnrSS_Lvl();
chromosome GnrChr_HEFT_t(vector<double> rank_t);
chromosome GnrChr_HEFT_b_t(vector<double> rank_b_t);
chromosome GnrChr_HEFT(vector<double> Rank_b);
chromosome GnrChr_HEFT_b_ADBRKGA(vector<double> Rank_b);
chromosome GnrChr_HEFT_t_ADBRKGA(vector<double> Rank_t);
chromosome GnrChr_DIHEFT(vector<double>& Rank_b);
chromosome GnrChr_IHEFT3_b(vector<double> Rank_b);
chromosome GnrChr_IHEFT3_t(vector<double> Rank_t);
double IHEFT3(chromosome& ch);
double ClcAvrReadyTime(int TskId, chromosome& chrom);
vector<double> GnrDecimalsByAscend();
chromosome GnrChr_Lvl_Ran();
chromosome GnrChr_HMEC(vector<double> Rank_b);
chromosome GnrChr_HEFT_Baseline();
chromosome GnrPrtByRank_Rnd(vector<double>& Rank);
chromosome GnrPrtByRank_EFT(vector<double>& Rank);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
double GnrML_Evl_EFT(chromosome& ch);
double HrsDcd_EFT(chromosome& ch);
double NrmDcd(chromosome& ch, bool IsFrw);
double HrsDcd_CTP(chromosome& ch);
double GnrML_Evl_MEC(chromosome& ch);
double DcdEvl(chromosome& ch, bool IsFrw);
void AdpDcd (chromosome&chrom ,double& CurTime,double& TotalTime);
void SelectRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
double CalculateEnergy(chromosome& Chrom);
double CalculateECByDelta2(chromosome& Chrom, int& NumBer);
void RepairPriorityAndGnrSchOrd(chromosome& chrom);
void RepairMapAndGnrRscAlcLst(chromosome& ch);
int FindNearestRscId(int TaskId, double value);
void UpdateParticle(chromosome& ch,chromosome& Pbest, chromosome& Gbest, double& runtime, double& SchTime);
double IFBD(chromosome& ch);
double IFBS(chromosome& ch);
void LBCA(chromosome& chrom);
void LBCA_IFBS(chromosome& ch);
void InitProModelOfResAlc(vector<vector<double>>& PMR);
void InitProModelOfTskSch(vector<vector<double>>& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
chromosome GnrTskLstOfChr(vector<vector<double> >& PMS);
void UpdatePMR(vector<vector<double>>& PMR, vector<chromosome>& Pop);
void UpdatePMS(vector<vector<double>>& PMS, vector<chromosome>& Pop);

#endif //FRAME_GENERATEACHROM_H
