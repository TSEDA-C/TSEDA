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

    SchTime["Epigenomics24_1.0"] = 1.235;
    SchTime["Epigenomics47_1.0"] = 4.139;
    SchTime["Epigenomics100_1.0"] = 22.057;

    SchTime["Ligo30_1.0"] = 1.346;
    SchTime["Ligo50_1.0"] = 3.032;
    SchTime["Ligo100_1.0"] = 11.675;

    SchTime["Montage25_1.0"] = 0.901;
    SchTime["Montage50_1.0"] = 3.367;
    SchTime["Montage100_1.0"] = 11.309;

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
           >> TemOrthogonal.EliteRate >> TemOrthogonal.ImprovementRate >> TemOrthogonal.RunTimeRatioOfStg1;
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
            cout << "Empty input file(fileList)" << endl;
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


