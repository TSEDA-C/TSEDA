#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "NGA.h"
#include "CGA.h"
#include "HGA.h"
#include "TSEDA.h"
#include "LWSGA.h"
#include "HPSO.h"
#include "ADBRKGA.h"

using namespace std;

int main() {
    srand((int) time(0));
    CDT = 20;

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
        double HGA_SchTime  = 0;
        int HGA_Iteration = 0;
        chromosome HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);
        ClearALL();
        double NGA_SchTime = 0;
        int NGA_Iteration = 0;
        chromosome NGA_Result = runNGA(XmlFile, RscAlcFile, NGA_SchTime, NGA_Iteration);
        ClearALL();
        double LWSGA_SchTime = 0;
        int LWSGA_Iteration = 0;
        chromosome LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);
        ClearALL();
        double ADBRKGA_SchTime = 0;
        int ADBRKGA_Iteration = 0;
        chromosome ADBRKGA_Result = runADBRKGA(XmlFile, RscAlcFile, ADBRKGA_SchTime, ADBRKGA_Iteration);
        ClearALL();
        double HPSO_SchTime = 0;
        int HPSO_Iteration = 0;
        chromosome HPSO_Result = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);
        ClearALL();
        double TSEDA_SchTime = 0;
        int TSEDA_Iteration = 0;
        chromosome TSEDA_Result = runTSEDA(XmlFile, RscAlcFile, TSEDA_SchTime, TSEDA_Iteration);
        ClearALL();
//        double CGA_SchTime = 0;
//        int CGA_Iteration = 0;
//        chromosome CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration);
//        ClearALL();

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(3);
//        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
//                << HGA_Result.EnergyConsumption << " " << HGA_SchTime << " " << HGA_Iteration << " "
//                << NGA_Result.EnergyConsumption << " " << NGA_SchTime << " " << NGA_Iteration << " "
//                << LWSGA_Result.EnergyConsumption << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
//                << ADBRKGA_Result.EnergyConsumption << " " << ADBRKGA_SchTime << " " << ADBRKGA_Iteration << " "
//                << HPSO_Result.EnergyConsumption << " " << HPSO_SchTime << " " << HPSO_Iteration << " "
////                << TSEDA_Result.EnergyConsumption << " " << TSEDA_SchTime << " " << TSEDA_Iteration << " "
////                << CGA_Result.EnergyConsumption << " " << CGA_SchTime << " " << CGA_Iteration << " "
//                << endl;
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << HGA_SchTime << " " << NGA_SchTime << " " << LWSGA_SchTime << " " << ADBRKGA_SchTime << " " << HPSO_SchTime << " "
                << TSEDA_SchTime << " "
                << endl;

        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}
