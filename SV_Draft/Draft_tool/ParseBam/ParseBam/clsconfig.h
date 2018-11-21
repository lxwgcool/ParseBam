#ifndef CLSCONFIG_H
#define CLSCONFIG_H

#include <string>
using namespace std;

struct St_Config
{
    //General

    //Debug
    string strVcf; //VCF
    string strBamFile; //BAMFILE
    string strRef; // Reference File
    string strMultiAlignSeq;
    string strPindelPath;
    string strLumpyPath;

    St_Config():strVcf(""), strBamFile(""), strRef(""), strMultiAlignSeq(""),
                strPindelPath(""), strLumpyPath("")
    {}
};

class ClsConfig
{
public:
    ClsConfig();

public:
    void ReadConfig(St_Config& stConfig, char* cpIniPath);
    bool CheckConfig(St_Config& stConfig);
};

#endif // CLSCONFIG_H
