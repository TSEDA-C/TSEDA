#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "TSEDA.h"

using namespace std;

int main() {
    srand((int) time(0));

    //{set the runtime (termination) time of the algorithm} -xy4
    map<string, double> SchTime;

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

    string strLine;
    ifstream iFile("../ExpParaSet.txt");
    if (!iFile) {
        cout << "filelist open failed!\n";
        exit(1);
    }
    while(getline(iFile,strLine)){
        istringstream is(strLine);
        Orthogonal TemOrthogonal;
        is >> TemOrthogonal.PopSizeFactor >> TemOrthogonal.theta1 >> TemOrthogonal.theta2
           >> TemOrthogonal.fdhi >> TemOrthogonal.ImprovementRate >> TemOrthogonal.RunTimeRatioOfStg1;
        orthogonal.push_back(TemOrthogonal);
    }
    iFile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string strLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, strLine);
        if (strLine.size() < 1) {
            cout << endl << "Empty input file(fileList)" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(strLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + ".xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        int index =0;
        for(Orthogonal TemOrthogonal: orthogonal){
            ++index;
            ofstream outfile("../OrthResultOutput/result"+to_string(index)+".txt", ios::app);
            if (!outfile) {
                cout << "Open the result file failure...\n";
                exit(0);
            }
            outfile.setf(ios::fixed, ios::floatfield);
            outfile.precision(3);
            cout <<endl<< "Parameter" + to_string(index) << " " << Model << " " << NumOfTask << " " << RscAvlRatio << " ";
            for (int times = 0; times < 10; ++times) {
                double TSEDA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
                int TSEDA_Iteration = 0;
                double TSEDA_Result = runTSEDA(XmlFile, RscAlcFile, TemOrthogonal, TSEDA_SchTime, TSEDA_Iteration);
                ClearALL();
                outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                        << TSEDA_Result << " " << TSEDA_SchTime << " " << TSEDA_Iteration
                        << endl;
            }
            outfile.close();
        }
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}


