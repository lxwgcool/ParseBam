#ifndef CLSCOMPARISON_H
#define CLSCOMPARISON_H
#include <string>
#include <vector>
#include "../../../../ShareLibrary/clsvcf1000genome.h"
using namespace std;

struct St_BreakPoint
{
    string strChromID;
    int iStart;
    int iEnd;
};

class ClsComparison
{
public:
    ClsComparison();

public:
    void ParsePindel(string strFilePath, vector<St_BreakPoint>& vBP);
    void CompareStdSVWithPindel(string strFilePath, vector<St_SV>& vSvDEL);
    void CompareStdSVWithLumpy(string strFilePath, vector<St_SV>& vSvDEL);
};

#endif // CLSCOMPARISON_H
