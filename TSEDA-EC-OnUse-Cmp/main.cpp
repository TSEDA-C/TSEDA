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

using namespace std;

int main() {
    srand((int) time(0));
    map<string, double> SchTime;
    //CDT=20
    SchTime["CyberShake30_1.0"] = 1.389;
    SchTime["CyberShake50_1.0"] = 2.334;
    SchTime["CyberShake100_1.0"] = 7.944;
    SchTime["Epigenomics24_1.0"] = 0.818;
    SchTime["Epigenomics47_1.0"] = 2.599;
    SchTime["Epigenomics100_1.0"] = 8.112;
    SchTime["Ligo30_1.0"] = 1.069;
    SchTime["Ligo50_1.0"] = 1.936;
    SchTime["Ligo100_1.0"] = 6.165;
    SchTime["Montage25_1.0"] = 0.794;
    SchTime["Montage50_1.0"] = 1.805;
    SchTime["Montage100_1.0"] = 5.258;
    SchTime["Sipht29_1.0"] = 1.438;
    SchTime["Sipht58_1.0"] = 4.940;
    SchTime["Sipht97_1.0"] = 15.177;

//    SchTime["CyberShake30_1.0"] = 3.039;
//    SchTime["CyberShake50_1.0"] = 16.348;
//    SchTime["CyberShake100_1.0"] = 293.321;
//    SchTime["Epigenomics24_1.0"] = 1.177;
//    SchTime["Epigenomics47_1.0"] = 6.459;
//    SchTime["Epigenomics100_1.0"] = 126.112;
//    SchTime["Ligo30_1.0"] = 1.789;
//    SchTime["Ligo50_1.0"] = 8.966;
//    SchTime["Ligo100_1.0"] = 95.097;
//    SchTime["Montage25_1.0"] = 1.113;
//    SchTime["Montage50_1.0"] = 21.875;
//    SchTime["Montage100_1.0"] = 344.133;
//    SchTime["Sipht29_1.0"] = 5.193;
//    SchTime["Sipht58_1.0"] = 41.285;
//    SchTime["Sipht97_1.0"] = 345.426;

//    SchTime["CyberShake30_1.0"] = 3.039;
//    SchTime["CyberShake50_1.0"] = 10.050;
//    SchTime["CyberShake100_1.0"] = 72.763;
//    SchTime["Epigenomics24_1.0"] = 1.177;
//    SchTime["Epigenomics47_1.0"] = 5.020;
//    SchTime["Epigenomics100_1.0"] = 24.991;
//    SchTime["Ligo30_1.0"] = 1.621;
//    SchTime["Ligo50_1.0"] = 3.736;
//    SchTime["Ligo100_1.0"] = 16.868;
//    SchTime["Montage25_1.0"] = 1.086;
//    SchTime["Montage50_1.0"] = 9.103;
//    SchTime["Montage100_1.0"] = 75.345;
//    SchTime["Sipht29_1.0"] = 5.193;
//    SchTime["Sipht58_1.0"] = 26.810;
//    SchTime["Sipht97_1.0"] = 124.423;

    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();
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
        double HMEC_SchTime  = 0;
        chromosome HMEC_Result = runHMEC(XmlFile, RscAlcFile, HMEC_SchTime);
        ClearALL();
        double HEFT_SchTime  = 0;
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
//        << NGA_Result.EnergyConsumption << " "
        << HEFT_Result.EnergyConsumption << " "
        << HGA_Result.EnergyConsumption << " " << NGA_Result.EnergyConsumption << " "
        << LWSGA_Result.EnergyConsumption << " " << ADBRKGA_Result.EnergyConsumption << " " << HPSO_Result.EnergyConsumption << " "
        << TSEDA_Result.EnergyConsumption << " "
        << HMEC_Result.EnergyConsumption << " "
        << HEFT_Result.MakeSpan << " "
        << HGA_Result.MakeSpan << " " << NGA_Result.MakeSpan << " "
        << LWSGA_Result.MakeSpan << " " << ADBRKGA_Result.MakeSpan << " " << HPSO_Result.MakeSpan << " "
        << TSEDA_Result.MakeSpan << " "
        << HMEC_Result.MakeSpan << " "
        << HEFT_SchTime << " "
        << HGA_SchTime << " " << NGA_SchTime << " "
        << LWSGA_SchTime << " " << ADBRKGA_SchTime << " " << HPSO_SchTime << " "
        << TSEDA_SchTime << " "
        << HMEC_SchTime << " "
        << HGA_Iteration << " " << NGA_Iteration << " "
        << LWSGA_Iteration << " " << ADBRKGA_Iteration << " " << HPSO_Iteration << " "
        << TSEDA_Iteration << " "
        << endl;

        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}
