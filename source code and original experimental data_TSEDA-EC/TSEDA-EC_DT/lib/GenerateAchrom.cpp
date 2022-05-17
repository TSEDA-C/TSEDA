//
// Created by qqq on 2021/9/15.
//

#include "GenerateAchrom.h"
#include <cstdlib>
#include "tools.h"
#include "common.h"
#include "GenOperator.h"
//{calculate the average execution time of tasks}
void W_Cal_Average(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        int RscSize = Tasks[i].ElgRsc.size();
        for (int j = 0; j < RscSize; ++j)
            sum += 1.0 / Rscs[Tasks[i].ElgRsc[j]].pc;
        w[i] = Tasks[i].length * sum / RscSize;
    }
}
//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

//calculate the rank of tasks based on independent IO using transfer time C[i][j]
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = w[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId=TskLstInLvl[i][j];
            double ChildMaxRankc = 0;
            for (int k = 0; k < Tasks[TaskId].children.size(); ++k) {
                int tem = Tasks[TaskId].children[k];
                double CompareObject = RankList[tem] + c[TaskId][tem];
                if(ChildMaxRankc  < CompareObject ){
                    ChildMaxRankc = CompareObject;
                }
            }
            RankList[TaskId] =w[TaskId] + ChildMaxRankc;
        }
    }
}

void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w) {
    for(int i =1 ;i < TskLstInLvl.size(); ++i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId = TskLstInLvl[i][j];
            for (int k = 0; k < Tasks[TaskId].parents.size(); ++k) {
                int tem = Tasks[TaskId].parents[k];
                double re = w[tem] + c[tem][TaskId] + RankList[tem];
                if (RankList[TaskId] < re) {
                    RankList[TaskId] = re;
                }
            }
        }
    }
}

chromosome GnrChr_HEFT_t(vector<double> Rank_t) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[i];
    }
    GnrML_Evl_EFT(TemChrom);
    CalculateEnergy(TemChrom);
    return TemChrom;
}

//{in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of rank_b_t}
chromosome GnrChr_HEFT_b_t(vector<double> Rank_b_t) {
    chromosome TemChrom;
    IntChr(TemChrom);
    int cur = 0;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        vector<double> TemList(TskLstInLvl[i].size());
        for (int j = 0; j < TskLstInLvl[i].size(); ++j)
            TemList[j] = Rank_b_t[TskLstInLvl[i][j]];
        vector<int> ind(TskLstInLvl[i].size());
        IndexSortByValueOnAscend(ind, TemList);
        for (int j = TskLstInLvl[i].size() - 1; j > -1; j--)
            TemChrom.TskSchLst[cur++] = TskLstInLvl[i][ind[j]];
    }
    GnrML_Evl_EFT(TemChrom);
    CalculateEnergy(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrML_Evl_EFT(TemChrom);
    CalculateEnergy(TemChrom);
    return TemChrom;
}

chromosome GnrChr_Lvl_Ran() {
    chromosome chrom;
    IntChr(chrom);
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()] + (LevelIdOfTask[i] + (rand() % 10000) / 10000.0) / TskLstInLvl.size();
    }
    return chrom;
}

chromosome GnrChr_HMEC(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrML_Evl_MEC(TemChrom);
    return TemChrom;
}

double GnrML_Evl_MEC(chromosome& ch){
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    ch.MakeSpan = 0;
    ch.EnergyConsumption = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1,TaskIndex = ch.TskSchLst[i];
        double FinalEC = 9999999999, FinalEndTime = 9999999999, FinalStartTime = 0;
        for (int v: Tasks[TaskIndex].ElgRsc) {
            double ReadyTime = 0;
            for (int ParentIndex : Tasks[TaskIndex].parents) {
                int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                double ftt = ch.EndTime[ParentIndex];
                if(v != ParentRscIndex){
                    double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                    ftt += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                }
                if (ReadyTime + PrecisionValue < ftt){
                    ReadyTime = ftt;
                }
            }
            double ExeTime = Tasks[TaskIndex].length  / Rscs[v].pc;
            ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[v], ExeTime, ReadyTime);
            ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExeTime;
            ch.RscAlcLst[TaskIndex] = v;
            double CurEC = CalculateECByDelta2(ch,i);
            //{find/record the earliest finish time}
            if (CurEC + PrecisionValue < FinalEC) {
                FinalStartTime = ch.StartTime[TaskIndex];
                FinalEndTime = ch.EndTime[TaskIndex];
                FinalEC = CurEC;
                RscIndex = v;
            }
        }
        ch.EnergyConsumption = FinalEC;
        ch.StartTime[TaskIndex] = FinalStartTime;
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        ch.MakeSpan = XY_MAX(ch.MakeSpan, FinalEndTime);
        UpdateITL(ITL[RscIndex],FinalStartTime,FinalEndTime);
    }
    return ch.EnergyConsumption;
}

// I/O independent
double DcdEvl(chromosome& ch, bool IsFrw) {
    ch.MakeSpan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);  a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i], RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the task and resource (Rsc) allocated to this task
        double ReadyTime = 0;
        if(IsFrw) {                              //forward-loading
            if (Tasks[TaskIndex].parents.size() != 0) {
                for (int j = 0; j < Tasks[TaskIndex].parents.size(); ++j) {
                    int ParentTask = Tasks[TaskIndex].parents[j], ParentRsc = ch.RscAlcLst[ParentTask];
                    double TransferTime = 0;
                    if(RscIndex != ParentRsc) {
                        TransferTime = ParChildTranFileSizeSum[ParentTask][TaskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw)); // -xy
                    }
                    double sum = ch.EndTime[ParentTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        } else {                                //backward-loading
            if (Tasks[TaskIndex].children.size() != 0) {
                for (int j = 0; j < Tasks[TaskIndex].children.size(); ++j) {
                    int ChildTask = Tasks[TaskIndex].children[j], ChildRsc = ch.RscAlcLst[ChildTask];
                    double TransferTime = 0;
                    if(RscIndex != ChildRsc) {
                        TransferTime = ParChildTranFileSizeSum[TaskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                    }
                    double sum = ch.EndTime[ChildTask] + TransferTime;
                    if (ReadyTime + PrecisionValue < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        }
        double ExecutionTime = Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[RscIndex],ExecutionTime,ReadyTime);
        ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExecutionTime;
        if (ch.MakeSpan + PrecisionValue < ch.EndTime[TaskIndex]) {
            ch.MakeSpan = ch.EndTime[TaskIndex];
        }
        //{update ITL}
        UpdateITL(ITL[RscIndex],ch.StartTime[TaskIndex],ch.EndTime[TaskIndex]);
    }
    return ch.MakeSpan;
}

void AdpDcd (chromosome& chrom, double& CurTime, double& TotalTime) {
    vector<double> SP(3);
    SP[0] = pow(CurTime/TotalTime,Parameter_ADBRKGA.alpha);
    SP[1] = Parameter_ADBRKGA.beta * (1-pow(CurTime/TotalTime,Parameter_ADBRKGA.alpha));
    SP[2] = (1-Parameter_ADBRKGA.beta) * (1-pow(CurTime/TotalTime,Parameter_ADBRKGA.alpha));
    double RandNum = double (rand()%1000) / 1000;
    if (RandNum < SP[0]){
        NrmDcd(chrom, true);
    }
    if (RandNum >= SP[0] && RandNum < (SP[0]+SP[1])){
        HrsDcd_EFT(chrom);
//        GnrRscAlcTskSchLstFromCode_RK(chrom);
//        GnrML_Evl_MEC(chrom);
//        ModifyRscAlcLstByCode_RK(chrom);
    }
    if (RandNum >= (SP[0]+SP[1])){
        HrsDcd_CTP(chrom);
    }
}


//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk, -1);
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk, -1);
    chrom.EndTime.resize(comConst.NumOfTsk);
    chrom.StartTime.resize(comConst.NumOfTsk);
    chrom.Code_TD.resize(comConst.NumOfRsc);
    chrom.TskSchPart.resize(comConst.NumOfTsk);
    chrom.RscAlcPart.resize(comConst.NumOfTsk);
    chrom.VTskSchPart.resize(comConst.NumOfTsk,0);
    chrom.VRscAlcPart.resize(comConst.NumOfTsk,0);
}

chromosome GnrChr_HEFT_Baseline() {
    chromosome TemChrom;
    IntChr(TemChrom);
    int ScheduleOrder = 0;
    for (int j = 0; j < TskLstInLvl.size(); ++j) {
        if (TskLstInLvl[j].size() < 2) {
            TemChrom.TskSchLst[ScheduleOrder++]=TskLstInLvl[j][0];
            continue;
        }
        vector<int> SonTaskNum;
        for (int i = 0; i < TskLstInLvl[j].size(); ++i)
            SonTaskNum.push_back(Tasks[TskLstInLvl[j][i]].children.size());

        vector<int> ind(TskLstInLvl[j].size());
        IndexSortByValueOnAscend(ind, SonTaskNum);
        for (int i = TskLstInLvl[j].size() - 1; i >= 0; i--) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][ind[i]];
        }
    }
    GnrML_Evl_EFT(TemChrom);
    CalculateEnergy(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT_b_ADBRKGA(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<double> Decimals = GnrDecimalsByAscend();
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.Code_RK[ind[comConst.NumOfTsk - i - 1]] = Decimals[i];
    }
    HrsDcd_EFT(TemChrom);
    CalculateEnergy(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT_t_ADBRKGA(vector<double> Rank_t) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<double> Decimals = GnrDecimalsByAscend();
    IndexSortByValueOnAscend(ind, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.Code_RK[ind[i]] = Decimals[i];
    }
    HrsDcd_EFT(TemChrom);
    CalculateEnergy(TemChrom);
    return TemChrom;
}

chromosome GnrChr_DIHEFT(vector<double>& Rank_b) {
    chromosome chrom;
    IntChr(chrom);
    vector<int> upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }
    vector<set<double>> ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    int IndexCount = 0;
    vector<double > AvrRtSet(comConst.NumOfTsk,0.0);
    while (!RTI.empty()){
        double max_RANK = 0;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double priority = Rank_b[*lit] + AvrRtSet[*lit];
            if (priority - PrecisionValue > max_RANK) {
                pit = lit; max_RANK = priority;
            }
        }
        int CrnTskId = *pit;
        RTI.erase(pit);
        chrom.TskSchLst[IndexCount] = CrnTskId;
        IndexCount++;
        int FinalRscForCrnTask = -1, NeedProcessChildTask = -1;
        double FinalStartTimeOfCrnTask = 0, FinalEndTimeOfCrnTask = InfiniteValue, TemMaxRnk = -1;
        vector<int> ReadyChildTaskSet;
        for (int m = 0; m < Tasks[CrnTskId].children.size(); ++m) {
            int childId = Tasks[CrnTskId].children[m];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                if (Rank_b[childId] - PrecisionValue > TemMaxRnk) {
                    NeedProcessChildTask = childId;
                    TemMaxRnk = Rank_b[NeedProcessChildTask];
                }
                ReadyChildTaskSet.push_back(childId);
            }
        }
        if(NeedProcessChildTask == -1){ //只需要处理当前任务
            SelectRsc_EFT(chrom,ITL,CrnTskId,FinalRscForCrnTask,FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
            chrom.StartTime[CrnTskId] = FinalStartTimeOfCrnTask;
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.MakeSpan = XY_MAX(chrom.MakeSpan,chrom.EndTime[CrnTskId]);
            UpdateITL(ITL[FinalRscForCrnTask],FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
        } else{  //需要和可以处理的子任务一起考虑分配资源
            int FinalRscForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = InfiniteValue;
            vector<set<double>> ITL_CrnTaskScheduled;
            for (int j = 0; j < Tasks[CrnTskId].ElgRsc.size(); ++j) {
                double ReadyTimeOfCrnTask = 0;
                int RscForCrnTsk = Tasks[CrnTskId].ElgRsc[j];
                chrom.RscAlcLst[CrnTskId] = RscForCrnTsk;
                for (int i2 = 0; i2 < Tasks[CrnTskId].parents.size(); ++i2) {
                    int ParentTask = Tasks[CrnTskId].parents[i2];
                    int RscOfParentTask = chrom.RscAlcLst[ParentTask];
                    double ftt = chrom.EndTime[ParentTask];
                    if(RscForCrnTsk != RscOfParentTask){
                        ftt += ParChildTranFileSizeSum[ParentTask][CrnTskId] / VALUE * 8 / (XY_MIN(Rscs[RscForCrnTsk].bw,Rscs[RscOfParentTask].bw));
                    }
                    if (ReadyTimeOfCrnTask + PrecisionValue < ftt){
                        ReadyTimeOfCrnTask = ftt;
                    }
                }
                double ExeTimeOfCrnTask = Tasks[CrnTskId].length / Rscs[RscForCrnTsk].pc;
                double CrnTaskStartTime = FindIdleTimeSlot(ITL[RscForCrnTsk],ExeTimeOfCrnTask,ReadyTimeOfCrnTask);
                double CrnTaskEndTime = CrnTaskStartTime + ExeTimeOfCrnTask;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL[RscForCrnTsk],CrnTaskStartTime,CrnTaskEndTime);
//                chrom.StartTime[CrnTskId] = CrnTaskStartTime;
                chrom.EndTime[CrnTskId] = CrnTaskEndTime;
                for (int j1 = 0; j1 < Tasks[NeedProcessChildTask].ElgRsc.size(); ++j1) {
                    double ChildReadyTime = 0;
                    int RscOfChildTsk = Tasks[NeedProcessChildTask].ElgRsc[j1];
                    for (int i4 = 0; i4 < Tasks[NeedProcessChildTask].parents.size(); ++i4) {
                        int TemParentTask = Tasks[NeedProcessChildTask].parents[i4];
                        int RscOfTemParTsk = chrom.RscAlcLst[TemParentTask];
                        double ftt = chrom.EndTime[TemParentTask];
                        if (RscOfChildTsk != RscOfTemParTsk) {
                            ftt += ParChildTranFileSizeSum[TemParentTask][NeedProcessChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscOfChildTsk].bw, Rscs[RscOfTemParTsk].bw));
                        }
                        if (ChildReadyTime + PrecisionValue < ftt) {
                            ChildReadyTime = ftt;
                        }
                    }
                    double ChildExeTime = Tasks[NeedProcessChildTask].length / Rscs[RscOfChildTsk].pc;
                    double ChildStartTime = FindIdleTimeSlot(TemITL[RscOfChildTsk],ChildExeTime,ChildReadyTime);
                    double ChildEndTime = ChildStartTime + ChildExeTime;
                    if (FinalEndTimeOfChildTask > ChildEndTime + PrecisionValue ) {
                        FinalRscForChildTask = RscOfChildTsk;
                        FinalStartTimeOfChildTask = ChildStartTime;
                        FinalEndTimeOfChildTask = ChildEndTime;
                        FinalStartTimeOfCrnTask = CrnTaskStartTime;
                        FinalEndTimeOfCrnTask = CrnTaskEndTime;
                        FinalRscForCrnTask = RscForCrnTsk;
                        ITL_CrnTaskScheduled = TemITL;
                    }
                }
            }
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.StartTime[CrnTskId] = FinalStartTimeOfCrnTask;
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.TskSchLst[IndexCount] = NeedProcessChildTask;
            IndexCount++;
            chrom.RscAlcLst[NeedProcessChildTask] = FinalRscForChildTask;
            chrom.StartTime[NeedProcessChildTask] = FinalStartTimeOfChildTask;
            chrom.EndTime[NeedProcessChildTask] = FinalEndTimeOfChildTask;
            chrom.MakeSpan = XY_MAX(chrom.MakeSpan,chrom.EndTime[NeedProcessChildTask]);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL[FinalRscForChildTask],FinalStartTimeOfChildTask,chrom.EndTime[NeedProcessChildTask]);
            for (int q = 0; q < ReadyChildTaskSet.size(); ++q) { //把所有还没有处理的就绪子任务添加到RTI中；
                int TskId = ReadyChildTaskSet[q];
                if (TskId != NeedProcessChildTask) {
                    RTI.push_back(TskId);
                    AvrRtSet[TskId] = ClcAvrReadyTime(TskId, chrom);
                }
            }
            for (int m = 0; m < Tasks[NeedProcessChildTask].children.size(); ++m) {//把就绪的子任务的子任务添加到RTI中；
                int childId = Tasks[NeedProcessChildTask].children[m];
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    RTI.push_back(childId);
                    AvrRtSet[childId] = ClcAvrReadyTime(childId, chrom);
                }
            }
        }
    }
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[chrom.TskSchLst[i]] = chrom.RscAlcLst[chrom.TskSchLst[i]] + Decimals[i];
    }
    CalculateEnergy(chrom);
    return chrom;
}

double ClcAvrReadyTime(int TskId, chromosome& chrom) {
    int j = Tasks[TskId].ElgRsc.size();
    double AvrRt = 0;
    for (int i = 0; i < Tasks[TskId].parents.size(); ++i) {
        int parent = Tasks[TskId].parents[i];
        int parentRsc = chrom.RscAlcLst[parent];
        double transfertime = 0;
        for (int k = 0; k < j; ++k) {
            int v = Tasks[TskId].ElgRsc[k];
            if (v != parentRsc) {
                transfertime += ParChildTranFileSizeSum[parent][TskId] / VALUE / (XY_MIN(Rscs[v].bw, Rscs[parentRsc].bw) / 8);
            }
        }
        double fft = transfertime / j + chrom.EndTime[parent];
        if (AvrRt < fft) {
            AvrRt = fft;
        }
    }
    return AvrRt;
}

vector<double> GnrDecimalsByAscend() {
    vector<double> decimals(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        decimals[i] = (i + RandomDouble(0,1)) / comConst.NumOfTsk;
    }
    return decimals;
}

chromosome GnrChr_IHEFT3_b(vector<double> Rank_b) {
    chromosome Chrom_Rank_b;
    IntChr(Chrom_Rank_b);
    IndexSortByValueOnDescend(Chrom_Rank_b.TskSchLst, Rank_b);
    IHEFT3(Chrom_Rank_b);
    CalculateEnergy(Chrom_Rank_b);
    return Chrom_Rank_b;
}

chromosome GnrChr_IHEFT3_t(vector<double> Rank_t) {
    chromosome Chrom_Rank_t;
    IntChr(Chrom_Rank_t);
    IndexSortByValueOnAscend(Chrom_Rank_t.TskSchLst, Rank_t);
    IHEFT3(Chrom_Rank_t);
    CalculateEnergy(Chrom_Rank_t);
    return Chrom_Rank_t;
}

double IHEFT3(chromosome& ch) {
    list <int> TemTskSchLst;
    TemTskSchLst.assign(ch.TskSchLst.begin(),ch.TskSchLst.end());
    ch.RscAlcLst.resize(comConst.NumOfTsk,-1);
    ch.TskSchLst.resize(comConst.NumOfTsk,-1);
    vector<set<double>> ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    vector<int> upr(comConst.NumOfTsk,0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
    }
    int IndexCount = 0;
    while (!TemTskSchLst.empty()){
        int CrnTask = TemTskSchLst.front();
        ch.TskSchLst[IndexCount] = CrnTask;
        IndexCount++;
        TemTskSchLst.erase(TemTskSchLst.begin());
        int FinalRscForCrnTask = -1;
        double FinalStartTimeOfCrnTask = 0;
        double FinalEndTimeOfCrnTask = InfiniteValue;
        vector<int> NeedProcessChildTaskSet;
        for (int m = 0; m < Tasks[CrnTask].children.size(); ++m) {
            int childId = Tasks[CrnTask].children[m];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                NeedProcessChildTaskSet.push_back(childId);
            }
        }
        if(NeedProcessChildTaskSet.empty()){
            SelectRsc_EFT(ch,ITL,CrnTask,FinalRscForCrnTask,FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
            ch.StartTime[CrnTask] = FinalStartTimeOfCrnTask;
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscForCrnTask;
            ch.MakeSpan = XY_MAX(ch.MakeSpan,ch.EndTime[CrnTask]);
            UpdateITL(ITL[FinalRscForCrnTask],FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
        } else{
            int FinalChildTask = -1 , FinalRscForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = InfiniteValue;
            vector<set<double> > ITL_CrnTaskScheduled;
            for (int j = 0; j < Tasks[CrnTask].ElgRsc.size(); ++j) {
                double ReadyTimeOfCrnTask = 0;
                int CrnTaskRsc = Tasks[CrnTask].ElgRsc[j];
                ch.RscAlcLst[CrnTask] = CrnTaskRsc;
                for (int i2 = 0; i2 < Tasks[CrnTask].parents.size(); ++i2) {
                    int ParentTask = Tasks[CrnTask].parents[i2];
                    int RscOfParentTask = ch.RscAlcLst[ParentTask];
                    double max = ch.EndTime[ParentTask];
                    if(CrnTaskRsc != RscOfParentTask){
                        max += ParChildTranFileSizeSum[ParentTask][CrnTask] / VALUE * 8 / (XY_MIN(Rscs[CrnTaskRsc].bw,Rscs[RscOfParentTask].bw));
                    }
                    if (ReadyTimeOfCrnTask + PrecisionValue < max){
                        ReadyTimeOfCrnTask = max;
                    }
                }
                double ExeTimeOfCrnTask = Tasks[CrnTask].length / Rscs[CrnTaskRsc].pc;
                double CrnTaskStartTime = FindIdleTimeSlot(ITL[CrnTaskRsc],ExeTimeOfCrnTask,ReadyTimeOfCrnTask);
                double CrnTaskEndTime = CrnTaskStartTime + ExeTimeOfCrnTask;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL[CrnTaskRsc],CrnTaskStartTime,CrnTaskEndTime);
//                ch.StartTime[CrnTask] = CrnTaskStartTime;
                ch.EndTime[CrnTask] = CrnTaskEndTime;
                for (int i3 = 0; i3 < NeedProcessChildTaskSet.size(); ++i3) {
                    int TemChildTask = NeedProcessChildTaskSet[i3];
                    for (int j1 = 0; j1 < Tasks[TemChildTask].ElgRsc.size(); ++j1) {
                        double TemChildReadyTime = 0;
                        int TemChildRsc = Tasks[TemChildTask].ElgRsc[j1];
                        for (int i4 = 0; i4 < Tasks[TemChildTask].parents.size(); ++i4) {
                            int TemParentTask = Tasks[TemChildTask].parents[i4];
                            int TemParRsc = ch.RscAlcLst[TemParentTask];
                            double max = ch.EndTime[TemParentTask];
                            if (TemChildRsc != TemParRsc) {
                                max += ParChildTranFileSizeSum[TemParentTask][TemChildTask] / VALUE * 8 / (XY_MIN(Rscs[TemChildRsc].bw, Rscs[TemParRsc].bw));
                            }
                            if (TemChildReadyTime + PrecisionValue < max) {
                                TemChildReadyTime = max;
                            }
                        }
                        double TemChildExeTime = Tasks[TemChildTask].length / Rscs[TemChildRsc].pc;
                        double TemChildStartTime = FindIdleTimeSlot(TemITL[TemChildRsc],TemChildExeTime,TemChildReadyTime);
                        double TemChildEndTime = TemChildStartTime + TemChildExeTime;
                        if (FinalEndTimeOfChildTask > TemChildEndTime + PrecisionValue ) {
                            FinalChildTask = TemChildTask;
                            FinalRscForChildTask = TemChildRsc;
                            FinalStartTimeOfChildTask = TemChildStartTime;
                            FinalEndTimeOfChildTask = TemChildEndTime;
                            FinalStartTimeOfCrnTask = CrnTaskStartTime;
                            FinalEndTimeOfCrnTask = CrnTaskEndTime;
                            FinalRscForCrnTask = CrnTaskRsc;
                            ITL_CrnTaskScheduled = TemITL;
                        }
                    }
                }
            }
            ch.StartTime[CrnTask] = FinalStartTimeOfCrnTask;
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscForCrnTask;
            ch.TskSchLst[IndexCount] = FinalChildTask;
            IndexCount++;
            ch.RscAlcLst[FinalChildTask] = FinalRscForChildTask;
            ch.StartTime[FinalChildTask] = FinalStartTimeOfChildTask;
            ch.EndTime[FinalChildTask] = FinalEndTimeOfChildTask;
            ch.MakeSpan = XY_MAX(ch.MakeSpan,ch.EndTime[FinalChildTask]);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL[FinalRscForChildTask],FinalStartTimeOfChildTask,ch.EndTime[FinalChildTask]);
            TemTskSchLst.erase(find(TemTskSchLst.begin(),TemTskSchLst.end(),FinalChildTask));
            for (int m = 0; m < Tasks[FinalChildTask].children.size(); ++m) {
                int childId = Tasks[FinalChildTask].children[m];
                upr[childId] = upr[childId] - 1;
            }
        }
    }
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.Code_RK[ch.TskSchLst[i]] = ch.RscAlcLst[ch.TskSchLst[i]] + Decimals[i];
    }
    return ch.MakeSpan;
}

//void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime){
//    if(ITL[RscIndex].find(StartTime) != ITL[RscIndex].end()) {
//        ITL[RscIndex].erase(StartTime);
//    } else {
//        ITL[RscIndex].insert(StartTime);
//    }
//    if(ITL[RscIndex].find(EndTime) != ITL[RscIndex].end()) {
//        ITL[RscIndex].erase(EndTime);
//    } else {
//        ITL[RscIndex].insert(EndTime);
//    }
//}

void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime){
    if(ITLofRscId.find(StartTime) != ITLofRscId.end()) {
        ITLofRscId.erase(StartTime);
    } else {
        ITLofRscId.insert(StartTime);
    }
    if(ITLofRscId.find(EndTime) != ITLofRscId.end()) {
        ITLofRscId.erase(EndTime);
    } else {
        ITLofRscId.insert(EndTime);
    }
}


double GnrML_Evl_EFT(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1, TaskId = ch.TskSchLst[i];
        double FinalEndTime = 9999999999, FinalStartTime = 0;
        SelectRsc_EFT(ch,ITL,TaskId,RscId,FinalStartTime,FinalEndTime);  //Find the resource that can finish the task earliest
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        UpdateITL(ITL[RscId],FinalStartTime,FinalEndTime);              //update ITL
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.MakeSpan = makespan;
    return makespan;
}

double HrsDcd_EFT(chromosome& ch) {
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk,0.0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0)  RTI.push_back(i);
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1;
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; //小数部分最小的那个优先调度
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);
        double FinalEndTime = InfiniteValue;
        double FinalStartTime = 0;
        SelectRsc_EFT(ch,ITL,ch.TskSchLst[i],RscId,FinalStartTime,FinalEndTime);
        ch.StartTime[ch.TskSchLst[i]] = FinalStartTime;
        ch.EndTime[ch.TskSchLst[i]] = FinalEndTime;
        ch.Code_RK[ch.TskSchLst[i]] = RscId + tmp;
        ch.RscAlcLst[ch.TskSchLst[i]] = RscId;
        UpdateITL(ITL[RscId],FinalStartTime,FinalEndTime); //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
        for (int l = 0; l < Tasks[ch.TskSchLst[i]].children.size(); ++l) {
            int ChildId = Tasks[ch.TskSchLst[i]].children[l];
            upr[ChildId] = upr[ChildId] - 1;
            if (upr[ChildId]==0)   RTI.push_back(ChildId);
        }
    }
    ch.MakeSpan = makespan;
    return makespan;
}  //ADBRKGA

double NrmDcd(chromosome& ch, bool IsFrw) {
    vector<int > upr(comConst.NumOfTsk,-1);
    list<int> RTI;
    if(IsFrw)
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].parents.size();
            if (upr[i]==0)  RTI.push_back(i);
        }
    else
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].children.size();
            if (upr[i]==0)  RTI.push_back(i);
        }
    //generate resource allocation list and task scheduling order list
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.RscAlcLst[i] = floor(ch.Code_RK[i]);
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            int decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; //小数部分最小的那个优先调度
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);
        //更新RTI;
        if (IsFrw)
            for (int l = 0; l < Tasks[ch.TskSchLst[i]].children.size(); ++l) {
                int childId = Tasks[ch.TskSchLst[i]].children[l];
                upr[childId] = upr[childId] - 1;
                if (upr[childId]==0)   RTI.push_back(childId);
            }
        else
            for (int l = 0; l < Tasks[ch.TskSchLst[i]].parents.size(); ++l) {
                int parentId = Tasks[ch.TskSchLst[i]].parents[l];
                upr[parentId] = upr[parentId] - 1;
                if (upr[parentId]==0)  RTI.push_back(parentId);
            }
    }
    DcdEvl(ch, IsFrw);
    return ch.MakeSpan;
}

double HrsDcd_CTP(chromosome& ch) {
    vector<double> w(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        if(Tasks[i].parents.size() !=  0){
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {
                int parent = Tasks[i].parents[j];
                int ParRsc = floor(ch.Code_RK[parent]);
                if(ParRsc != RscIndex){
                    TransferTime[parent][i] = ParChildTranFileSizeSum[parent][i] / VALUE * 8 / XY_MIN(Rscs[RscIndex].bw,Rscs[ParRsc].bw) ;
                }
            }
        }
        w[i] = Tasks[i].length / Rscs[RscIndex].pc;
    }
    Calculate_Rank_b(Rank_b,TransferTime, w);
    IndexSortByValueOnAscend(ind, Rank_b);
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.TskSchLst[i] = TaskIndex;
        ch.Code_RK[TaskIndex] = floor(ch.Code_RK[TaskIndex]) + Decimals[i];
    }
    NrmDcd(ch, true);
    return ch.MakeSpan;
}

void SelectRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime) {
    for (int RscIdOfCrnTsk : Tasks[TaskIndex].ElgRsc) {
        double ReadyTime = 0;
        for (int ParentIndex: Tasks[TaskIndex].parents) { //calculate the ready time of the task
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            double FTT = ch.EndTime[ParentIndex];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                FTT += TransferData / VALUE * 8 / (XY_MIN(Rscs[RscIdOfCrnTsk].bw,Rscs[RscIdOfPrnTsk].bw));
            }
            if (ReadyTime + PrecisionValue < FTT){
                ReadyTime = FTT;
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc;
        double StartTime = FindIdleTimeSlot(ITL[RscIdOfCrnTsk],ExeTime,ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = RscIdOfCrnTsk;
        }
    }
}

double FindIdleTimeSlot(set<double>& ITLofRscId,double& ExeTime,double& ReadyTime){
    set<double>::iterator pre  = ITLofRscId.begin();
    set<double>::iterator post = ITLofRscId.begin();
    ++post;
    while(post != ITLofRscId.end()) {
        if((*post - *pre) > ExeTime - PrecisionValue && ReadyTime - PrecisionValue < (*post)-ExeTime) {
            return  XY_MAX(*pre, ReadyTime);
        } else {
            ++pre; ++pre; ++post; ++post;
        }
    }
}

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }

    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}

vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
}

void InitProModelOfResAlc(vector<vector<double> >& PMR) {
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int j : Tasks[i].ElgRsc) {
            PMR[i][j] =  1.0 / Tasks[i].ElgRsc.size();
        }
    }
}

void InitProModelOfTskSch(vector<vector<double>>& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants) {
    vector<int> STS(comConst.NumOfTsk, 0);
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int k = NumOfAncestors[i]; k < NumOfNonDescendants[i]; ++k) {
            PMS[i][k] = 1;
            ++STS[k];
        }
    }//PMS[i][k] represents the probability that the k-th scheduled task is task i

    for(int k = 0; k < comConst.NumOfTsk; ++k) {
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            PMS[i][k] = PMS[i][k] / STS[k];
        }
    }
}

chromosome GnrTskLstOfChr(vector<vector<double> >& PMS) {
    chromosome chrom;
    IntChr(chrom);
    vector<int > upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0)  RTI.push_back(i);
    }
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double sum = 0;
        for(int k : RTI){
            sum += PMS[k][i];
        }
        vector<double> SltProb(comConst.NumOfTsk);
        for (int k : RTI) {
            SltProb[k] = PMS[k][i] / sum;
        }
        double rnd = double(rand()%100) / 100;
        double ProbSum = 0;
        int taskIndex;
        for (int k : RTI) {
            ProbSum += SltProb[k];
            if (rnd + PrecisionValue < ProbSum) {
                taskIndex = k;
                break;
            }
        }
        chrom.TskSchLst[i] = taskIndex;
        RTI.erase(find(RTI.begin(), RTI.end(), taskIndex));
        for (int k = 0; k < Tasks[taskIndex].children.size(); ++k) {
            upr[Tasks[taskIndex].children[k]]--;
            if (upr[Tasks[taskIndex].children[k]] == 0){
                RTI.push_back(Tasks[taskIndex].children[k]);
            }
        }
    }
    return chrom;
}

void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR) {
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double rnd = double(rand()%100) / 100;
        double sum = 0;
        for(int j: Tasks[i].ElgRsc){
            sum += PMR[i][j];
            if(rnd < sum) {
                chrom.RscAlcLst[i] = j;
                break;
            }
        }
    }
}

chromosome GnrPrtByRank_Rnd(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        chrom.RscAlcPart[i] =RandomDouble2(0,comConst.NumOfRsc-1); //RandomDouble2(0,comConst.NumOfRsc);//rand() % comConst.NumOfRsc + rand() % 1000 / 1000.0 - 0.5;
    }
    RepairMapAndGnrRscAlcLst(chrom); //GnrRscAlcLst(chrom); //
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    DcdEvl(chrom, true);
    CalculateEnergy(chrom);
    return chrom;
}

chromosome GnrPrtByRank_EFT(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    GnrML_Evl_EFT(chrom);
    CalculateEnergy(chrom);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.RscAlcPart[i] = chrom.RscAlcLst[i] - 0.5 + (rand() % 10000) / 10000.0;
    }
    return chrom;
}

void RepairPriorityAndGnrSchOrd(chromosome& chrom) {
    vector<int> V, Q;
    vector<int> N(comConst.NumOfTsk, -1);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        N[i] = round(chrom.TskSchPart[i]);
    }
    vector<int> upr(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0) {
            Q.push_back(i);
        }
    }
    int MaxV = -1;
    while (V.size() != comConst.NumOfTsk) {
        for (int i = 0; i < Q.size(); ++i) {
            int TaskId = Q[i];
            int MaxP = -1;
            for (int i1 = 0; i1 < Tasks[TaskId].parents.size(); ++i1) {
                if (MaxP < N[Tasks[TaskId].parents[i1]]) {
                    MaxP = N[Tasks[TaskId].parents[i1]];
                }
            }
            if (N[TaskId] <= MaxP) {
                N[TaskId] = MaxP + 1;
            }
            for (int i1 = 0; i1 < V.size(); ++i1) {
                if (N[TaskId] == N[V[i1]])  {
                    N[TaskId] = MaxV + 1;
                    MaxV += 1;
                    break;
                }
            }
            MaxV = XY_MAX(N[TaskId], MaxV);
            V.push_back(TaskId);
        }
        vector<int> TemQ;
        for (int i = 0; i < Q.size(); ++i) {
            int taskId = Q[i];
            for (int i2 = 0; i2 < Tasks[taskId].children.size(); ++i2) {
                int childId = Tasks[taskId].children[i2];
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    TemQ.push_back(childId);
                }
            }
        }
        Q = TemQ;
    }
//    chrom.TskSchPart = N;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.TskSchPart[i] = N[i];
    }
    IndexSortByValueOnAscend(chrom.TskSchLst, N);
}

void RepairMapAndGnrRscAlcLst(chromosome& ch) {
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        int RscId = round(ch.RscAlcPart[i]);
        if(RscId < Tasks[i].ElgRsc[0]) { //超出下限的处理
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            continue;
        }
        if(RscId > Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1]) { //超出上限的处理
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1];
            continue;
        }
        if(find(Tasks[i].ElgRsc.begin(), Tasks[i].ElgRsc.end(), RscId) == Tasks[i].ElgRsc.end()){ //不存在的处理
            if(Tasks[i].ElgRsc.size() == 1) {
                ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
                ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            } else {
                int TemRscId = FindNearestRscId(i, ch.RscAlcPart[i]);
                ch.RscAlcPart[i] = TemRscId;
                ch.RscAlcLst[i] = TemRscId;
            }
            continue;
        }
        ch.RscAlcLst[i] = RscId;
    }
}

int FindNearestRscId(int TaskId, double value) {
    for (int j = 0; j < Tasks[TaskId].ElgRsc.size()-1; ++j ){
        if (Tasks[TaskId].ElgRsc[j] < value && value < Tasks[TaskId].ElgRsc[j+1] ) {
            if ( Tasks[TaskId].ElgRsc[j+1] - value < value - Tasks[TaskId].ElgRsc[j] ) {
                return Tasks[TaskId].ElgRsc[j+1];
            } else {
                return Tasks[TaskId].ElgRsc[j];
            }
        }
    }
}

void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime){
    Parameter_HPSO.InertiaWeight = 0.1 * (1-(runtime / SchTime)) + 0.9;
    Parameter_HPSO.c1 = 2 * (1-(runtime / SchTime));
    Parameter_HPSO.c2 = 2 * (runtime / SchTime);
    double r1 = RandomDouble(0,1);
    double r2 = RandomDouble(0,1);
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        ch.VTskSchPart[i] = Parameter_HPSO.InertiaWeight * ch.VTskSchPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.TskSchPart[i] - ch.TskSchPart[i])
                            + Parameter_HPSO.c2 * r2 * (Gbest.TskSchPart[i] - ch.TskSchPart[i]);
        ch.TskSchPart[i] += ch.VTskSchPart[i];

        ch.VRscAlcPart[i] = Parameter_HPSO.InertiaWeight * ch.VRscAlcPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.RscAlcPart[i] - ch.RscAlcPart[i])
                            + Parameter_HPSO.c2 * r2 * (Gbest.RscAlcPart[i] - ch.RscAlcPart[i]);
        ch.RscAlcPart[i] += ch.VRscAlcPart[i];
    }
    RepairMapAndGnrRscAlcLst(ch);   //GnrRscAlcLst(ch); //
    RepairPriorityAndGnrSchOrd(ch);
}

double IFBD(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSortByValueOnAscend(ind, OldChrom.EndTime);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.TskSchLst[comConst.NumOfTsk - 1 - i] = ind[i];
        }
        DcdEvl(NewChrom, IsFrw);
        CalculateEnergy(NewChrom);
        IsFrw = !IsFrw;
    } while (NewChrom.EnergyConsumption + PrecisionValue < OldChrom.EnergyConsumption);
    if (IsFrw) {
        ch = OldChrom;
    } else {
        ch = NewChrom;
    }
    return ch.EnergyConsumption;
}

double IFBS(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    vector<double> Decimals = GnrDecimalsByAscend();
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSortByValueOnAscend(ind, OldChrom.EndTime);
        //vector<double> Decimals = GnrDecimalsByAscend();
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]] = floor(NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]]) + Decimals[i];
        }
        NrmDcd(NewChrom, IsFrw);
        CalculateEnergy(NewChrom);
        IsFrw = !IsFrw;
    } while (NewChrom.EnergyConsumption + PrecisionValue < OldChrom.EnergyConsumption);
    if (IsFrw) { //the last is backward
        ch = OldChrom;
    } else {
        ch = NewChrom;
    }
    return ch.EnergyConsumption;
} //ADBRKGA

//{ Load Balancing with Communication Reduction Improvement (LBCRI)}
void LBCA(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Id(comConst.NumOfRsc,0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        Id[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc;
        TSK[RscIndex].push_back(i);
    }
    vector<int> ind(comConst.NumOfRsc);
    IndexSortByValueOnAscend(ind, Id);         //sorting according to loads
    int RscWithMinLd = ind[0];          //find out the resource (Rsc) with the lowest load;
    set<int> ST;
    if (abs(Id[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    } else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            int TaskIndex = TSK[RscWithMinLd][i];
            ST.insert(Tasks[TaskIndex].children.begin(),Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(),Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==
                Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if(ST.empty()){
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Id[ch.RscAlcLst[s]]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    ch.RscAlcLst[t[0].first] = RscWithMinLd;
    DcdEvl(ch, true);
    CalculateEnergy(ch);
    IFBD(ch);
    if (OldCh.EnergyConsumption+ PrecisionValue < ch.EnergyConsumption) {
        ch = OldCh;
    }
}

void LBCA_IFBS(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Ld(comConst.NumOfRsc,0.0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        Ld[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc;
        TSK[RscIndex].push_back(i);
    }
    int RscWithMinLd = 0;
    double TemLd = Ld[0];
    for (int j = 1; j < comConst.NumOfRsc; ++j) { //find out the resource (Rsc) with the lowest load -new-xy;
        if (TemLd > Ld[j]) {
            TemLd = Ld[j]; RscWithMinLd = j;
        }
    }
    set<int> ST;
    if (abs(Ld[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    } else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            int TaskIndex = TSK[RscWithMinLd][i];
            ST.insert(Tasks[TaskIndex].children.begin(),Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(),Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==  Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if(ST.empty()){//-w
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Ld[floor(ch.Code_RK[s])]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    double decimal = ch.Code_RK[t[0].first] - floor(ch.Code_RK[t[0].first]);
    ch.Code_RK[t[0].first] = RscWithMinLd + decimal;
    NrmDcd(ch, true);
    CalculateEnergy(ch);
    IFBS(ch);
    if (OldCh.EnergyConsumption + PrecisionValue < ch.EnergyConsumption) {
        ch = OldCh;
    }
}



double CalculateEnergy(chromosome &Chrom){
    Chrom.EnergyConsumption = 0;
    set<double>AllTime;
    AllTime.insert( Chrom.EndTime.begin(),Chrom.EndTime.end());
    AllTime.insert(Chrom.StartTime.begin(),Chrom.StartTime.end());
    vector<double>T;
    T.assign(AllTime.begin(),AllTime.end());
    for(int x = 0; x < T.size()-1; ++x) {
        if (fabs(T[x+1]-T[x]) < 1e-5) {
            continue;
        }
        vector<double> HstLd(HstSet.size(),0);
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            int taskId = Chrom.TskSchLst[i];
            int RscId = Chrom.RscAlcLst[taskId];
            if(Chrom.StartTime[taskId] - PrecisionValue < T[x] && T[x+1]  - PrecisionValue < Chrom.EndTime[taskId]) {
                int taskHT = Rscs[RscId].Hostid;
                double TemLd = Rscs[RscId].pc / HstSet[taskHT].pc;
                HstLd[taskHT] += TemLd;
            }
        }
        for(int k = 0; k < HstSet.size(); ++k) {
            Chrom.EnergyConsumption += (T[x+1] - T[x]) * CalculatePowerByLoad(HstLd[k],k);
        }
    }
    return Chrom.EnergyConsumption;
}

double CalculateECByDelta2(chromosome &ch, int &Number){
    double DeltaEC = 0;
    vector<int>STn;
    int CurTask = ch.TskSchLst[Number];
    int CurHT = Rscs[ch.RscAlcLst[CurTask]].Hostid;
    for (int i = 0; i < Number; ++i) {
        int taskId = ch.TskSchLst[i];
        int RscId = ch.RscAlcLst[taskId];
        int taskHT = Rscs[RscId].Hostid;
        if((ch.EndTime[taskId] > ch.StartTime[CurTask] + PrecisionValue && ch.EndTime[taskId] - PrecisionValue < ch.EndTime[CurTask] && taskHT == CurHT) ||
           (ch.StartTime[taskId] > ch.StartTime[CurTask] - PrecisionValue && ch.StartTime[taskId] + PrecisionValue < ch.EndTime[CurTask] && taskHT == CurHT) ||
           (ch.StartTime[taskId] + PrecisionValue < ch.StartTime[CurTask] && ch.EndTime[taskId] > ch.EndTime[CurTask] + PrecisionValue && taskHT == CurHT)){ //加精度控制-xy-已改-qmq
            STn.push_back(taskId);
        }
    }
    set<double>AllTime;
    vector<double>T;
    AllTime.insert(ch.StartTime[CurTask]);
    AllTime.insert(ch.EndTime[CurTask]);
    for (int i = 0; i < STn.size(); ++i) {
        if (ch.StartTime[STn[i]] > ch.StartTime[CurTask] + PrecisionValue ){
            AllTime.insert(ch.StartTime[STn[i]]);
        }
        if (ch.EndTime[STn[i]] + PrecisionValue < ch.EndTime[CurTask]){
            AllTime.insert(ch.EndTime[STn[i]]);
        }
    }
    T.assign(AllTime.begin(),AllTime.end());
    for(int x = 0; x < T.size()-1; ++x) {
        if (fabs(T[x+1]-T[x]) < 1e-5) {
            continue;
        }
        double HstLd = 0;
        for(int i = 0; i < STn.size(); ++i) {
            int taskId = STn[i];
            int RscId = ch.RscAlcLst[taskId];
            if(ch.StartTime[taskId] - PrecisionValue < T[x] && T[x+1]  - PrecisionValue < ch.EndTime[taskId]) {
                double TemLd = Rscs[RscId].pc / HstSet[CurHT].pc;
                HstLd += TemLd;
            }
        }
        double CurLd = HstLd + Rscs[ch.RscAlcLst[CurTask]].pc / HstSet[CurHT].pc;
        double TemEc = (T[x+1] - T[x]) * (CalculatePowerByLoad(CurLd,CurHT) - CalculatePowerByLoad(HstLd,CurHT));
        DeltaEC += TemEc;
    }
    if (ch.EndTime[CurTask] > ch.MakeSpan + PrecisionValue ){
        for (int i = 0; i < HstSet.size(); ++i) {
            DeltaEC += (ch.EndTime[CurTask] - ch.MakeSpan) * CalculatePowerByLoad(0,i);
        }
    }
    return ch.EnergyConsumption + DeltaEC;
}

void UpdatePMR(vector<vector<double>>& PMR, vector<chromosome>& Pop){
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            int count = 0;
            for (int n = 0; n < Parameter_TSEDA.NumOfEliteOfPop; ++n) {
                if (Pop[n].RscAlcLst[i] == j) {
                    ++count;
                }
            }
            PMR[i][j] = (1 - Parameter_TSEDA.theta1) * PMR[i][j] + Parameter_TSEDA.theta1 * (1.0 * count / Parameter_TSEDA.NumOfEliteOfPop);
        }
    }
}

void UpdatePMS(vector<vector<double>>& PMS, vector<chromosome>& Pop){
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int j = 0; j < comConst.NumOfTsk; ++j) {
            int count = 0;
            for(int n = 0; n < Parameter_TSEDA.NumOfEliteOfPop; ++n) {
                if(Pop[n].TskSchLst[i] == j) {
                    ++count;
                }
            }
            PMS[j][i] = (1-Parameter_TSEDA.theta2)*PMS[j][i] + Parameter_TSEDA.theta2*(1.0*count/Parameter_TSEDA.NumOfEliteOfPop);
        }
    }
}

