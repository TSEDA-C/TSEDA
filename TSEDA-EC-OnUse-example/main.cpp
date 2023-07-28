#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "NGA.h"
#include "CGA.h"
#include "HGA.h"
#include "TSEDA.h"
#include "LWSGA.h"
#include "HEFT.h"
#include "HMEC.h"
#include "HPSO.h"
#include "ADBRKGA.h"
#include "Exhaustive.h"
#include "tools.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"

using namespace std;

int main() {
    srand((int) time(0));
    map<string, double> SchTime;
    SchTime["CyberShake8_1.0"] = 2.0;

    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();
    ofstream PrintChromFile("../PrintChrom.txt", ios::out);
    PrintChromFile.close();

//    ReadFile("CyberShake_8.xml","8_1.0_0.txt");
//    CalculateLevelList();                           //calculate the levels of tasks
//    chromosome tstChrom;
//    IntChr(tstChrom);

//    //测试HGA是否能找到最优解
//    tstChrom.RscAlcLst = {1,1,1,1,1,1,2,1};
//    tstChrom.RscAlcLst = {2,2,2,2,2,2,1,2};
//    GnrTskSchLst_HGA_S(tstChrom);
//    DcdEvl_S(tstChrom);
//    CalculateEnergy(tstChrom);
//    PrintChrom(tstChrom);
//    //测试NGA是否能找到最优解
//    tstChrom.TskSchLst = {0,1,2,4,3,5,6,7};
//    tstChrom.TskSchLst = {0,1,2,4,3,6,5,7};
//    tstChrom.TskSchLst = {0,1,2,4,3,5,7,6};
//    tstChrom.TskSchLst = {0,1,2,4,6,3,5,7};
//    tstChrom.TskSchLst = {0,2,1,4,3,5,6,7};
//    tstChrom.TskSchLst = {0,2,1,4,3,6,5,7};
//    tstChrom.TskSchLst = {0,2,1,4,3,5,7,6};
//    tstChrom.TskSchLst = {0,2,1,4,6,3,5,7};
//    GnrML_Evl_EFT_S(tstChrom);
//    CalculateEnergy(tstChrom);
//    PrintChrom(tstChrom);

//    tstChrom.TskSchLst = {0,1,2,4,3,6,5,7}; tstChrom.RscAlcLst = {1,1,1,1,1,0,2,0};
//    DcdEvl_S(tstChrom);
//    CalculateEnergy(tstChrom);
////    LBCA_S(tstChrom);
//    PrintIFBD_S(tstChrom);
//    exit(3);

    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << endl << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + ".xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";

        double Exh_SchTime  = 0 ;
        chromosome bestChrom = runExhaustive(XmlFile, RscAlcFile,Exh_SchTime);
        ClearALL();
        double HEFT_SchTime  = 0 ;
        chromosome HEFT_Result = runHEFT(XmlFile, RscAlcFile,HEFT_SchTime);
        ClearALL();
        double HGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HGA_Iteration = 0;
        chromosome HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);
        ClearALL();
        double NGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int NGA_Iteration = 0;
        chromosome NGA_Result = runNGA(XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration);
        ClearALL();
        double LWSGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int LWSGA_Iteration = 0;
        chromosome LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);
        ClearALL();
        double ADBRKGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int ADBRKGA_Iteration = 0;
        chromosome ADBRKGA_Result = runADBRKGA(XmlFile, RscAlcFile, ADBRKGA_SchTime, ADBRKGA_Iteration);
        ClearALL();
        double HPSO_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HPSO_Iteration = 0;
        chromosome HPSO_Result = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);
        ClearALL();
        double TSEDA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int TSEDA_Iteration = 0;
        chromosome TSEDA_Result = runTSEDA(XmlFile, RscAlcFile, TSEDA_SchTime, TSEDA_Iteration);
        ClearALL();
//        double HMEC_SchTime  = 0 ;
//        chromosome HMEC_Result = runHMEC(XmlFile, RscAlcFile,HMEC_SchTime);
//        ClearALL();
//        double CGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int CGA_Iteration = 0;
//        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration);
//        ClearALL();

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << bestChrom.MakeSpan << " " << bestChrom.EnergyConsumption << " " << Exh_SchTime << " "
                << HEFT_Result.MakeSpan << " " << HEFT_Result.EnergyConsumption << " " << HEFT_SchTime << " "
                << HGA_Result.MakeSpan << " " << HGA_Result.EnergyConsumption << " " << HGA_SchTime << " "
                << NGA_Result.MakeSpan << " " << NGA_Result.EnergyConsumption << " " << NGA_SchTime << " "
                << LWSGA_Result.MakeSpan << " " << LWSGA_Result.EnergyConsumption << " " << LWSGA_SchTime << " "
                << ADBRKGA_Result.MakeSpan << " " << ADBRKGA_Result.EnergyConsumption << " " << ADBRKGA_SchTime << " "
                << HPSO_Result.MakeSpan << " " << HPSO_Result.EnergyConsumption << " " << HPSO_SchTime << " "
                << TSEDA_Result.MakeSpan << " " << TSEDA_Result.EnergyConsumption << " " << TSEDA_SchTime << " "
//                << HMEC_Result.MakeSpan << " " << HMEC_Result.EnergyConsumption << " " << HMEC_SchTime << " "
//                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "
                << endl;
        PrintChrom(bestChrom);
        PrintChrom(HEFT_Result);
        PrintChrom(HGA_Result);
        PrintChrom(NGA_Result);
        PrintChrom(LWSGA_Result);
        PrintChrom(ADBRKGA_Result);
        PrintChrom(HPSO_Result);
        PrintChrom(TSEDA_Result);
//        PrintChrom(HMEC_Result);

        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}
