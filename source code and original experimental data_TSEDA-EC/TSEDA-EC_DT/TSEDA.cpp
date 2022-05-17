//
// Created by qqq on 2021/10/15.
//

#include "TSEDA.h"
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"

double runTSEDA(string XmlFile, string RscAlcFile,double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_TSEDA();
    CalculateLevelList();
    CalculateDescendants();
    CalculateAncestors();
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b, cc, ww);
    vector<int> NumOfAncestors(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfAncestors[i] = Ancestors[i].size();
    }
    vector<int> NumOfDescendants(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfDescendants[i] = Descendants[i].size();
    }
    vector<int> NumOfNonDescendants(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfNonDescendants[i] = comConst.NumOfTsk - NumOfDescendants[i];
    }

    vector<vector<double>> PMR(comConst.NumOfTsk, vector<double>(comConst.NumOfRsc, 0));
    InitProModelOfResAlc(PMR); //PMR[i][j] represents the probability that task i is assigned to resource j;
    vector<vector<double>> PMS(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    InitProModelOfTskSch(PMS, NumOfAncestors, NumOfNonDescendants); //PMS[i][k] represents the probability that the k-th scheduled task is task i

    vector<chromosome> population, NewPopulation(Parameter_TSEDA.NumOfChromPerPop);

    NewPopulation[0] = GnrChr_HMEC(Rank_b);         //a chromosome according to HMEC is seeded into the population;
    NewPopulation[1] = GnrChr_HEFT(Rank_b);         //a chromosome according to HEFT is seeded into the population;
    int num = 2;

    double bestFitness = NewPopulation[0].EnergyConsumption;
    int terminationNum = ceil(50 * sqrt(ModelScale)/Parameter_TSEDA.NumOfChromPerPop);
    int NumOfNoImpGen = 0;

    while (NumOfNoImpGen < terminationNum * Parameter_TSEDA.RunTimeRatioOfStg1) {
        #pragma omp parallel for
        for (int n = num; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            NewPopulation[n] = GnrTskLstOfChr(PMS);
            GnrML_Evl_MEC(NewPopulation[n]);
        }

        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());
        sort(population.begin(), population.end(), SortByEnergyConsumption);
        vector<chromosome>::iterator Iter = population.begin();
        advance(Iter, Parameter_TSEDA.NumOfChromPerPop);
        population.assign(population.begin(), Iter);

        UpdatePMR(PMR, population);
        UpdatePMS(PMS, population);
        ++iteration; num = 0;
        ++NumOfNoImpGen;
        if( population[0].EnergyConsumption + PrecisionValue < bestFitness ){
            bestFitness = population[0].EnergyConsumption;
            NumOfNoImpGen = 0;
        }
    }

    while (NumOfNoImpGen < terminationNum ){
        #pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            NewPopulation[n] = GnrTskLstOfChr(PMS);
            GnrRscLstOfChr(NewPopulation[n],PMR);
            DcdEvl(NewPopulation[n],true);
            CalculateEnergy(NewPopulation[n]);
        }

        sort(NewPopulation.begin(), NewPopulation.end(), SortByEnergyConsumption);
        #pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfImproveOfPop; ++n) {
            IFBD(NewPopulation[n]);
            LBCA(NewPopulation[n]);
        }

        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());
        sort(population.begin(), population.end(), SortByEnergyConsumption);
        vector<chromosome>::iterator Iter = population.begin();
        advance(Iter, Parameter_TSEDA.NumOfChromPerPop);
        population.assign(population.begin(), Iter);

        UpdatePMR(PMR, population);
        UpdatePMS(PMS, population);
        ++iteration;
        ++NumOfNoImpGen;
        if( population[0].EnergyConsumption + PrecisionValue < bestFitness ){
            bestFitness = population[0].EnergyConsumption;
            NumOfNoImpGen = 0;
        }
    }
    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    return population[0].EnergyConsumption;
}