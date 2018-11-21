#include "clsconfig.h"
#include "../../../../ShareLibrary/clsreadconfigini.h"

ClsConfig::ClsConfig()
{

}

void ClsConfig::ReadConfig(St_Config& stConfig, char* cpIniPath)
{
    ClsReadConfigIni* pIni = new ClsReadConfigIni();
    pIni->ReadIni(cpIniPath); // The first valid parameters

    for(vector<St_Section>::iterator itr = pIni->GetConfigInfo().vSect.begin();
        itr != pIni->GetConfigInfo().vSect.end(); itr++)
    {
        if(itr->strName == "General")
        {
        }
        else if(itr->strName == "Debug")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "VCF")
                    stConfig.strVcf = itrmp->second;
                else if(itrmp->first == "BamFile")
                    stConfig.strBamFile = itrmp->second.c_str();
                else if(itrmp->first == "MultiAlignSeq")
                    stConfig.strMultiAlignSeq = itrmp->second.c_str();
                else if(itrmp->first == "Ref")
                    stConfig.strRef = itrmp->second.c_str();
                else if(itrmp->first == "PindelPath")
                    stConfig.strPindelPath = itrmp->second.c_str();
                else if(itrmp->first == "LumpyPath")
                    stConfig.strLumpyPath = itrmp->second.c_str();
                else{}
            }
        }
        else{}
    }

    delete pIni;
    pIni = NULL;
}

bool ClsConfig::CheckConfig(St_Config& stConfig)
{
    return true;
}
