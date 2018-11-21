#include "clscomparison.h"
#include <fstream>
#include "../../../../ShareLibrary/clsbasealgorithm.h"

ClsComparison::ClsComparison()
{

}

void ClsComparison::ParsePindel(string strFilePath, vector<St_BreakPoint>& vBP)
{
    fstream ifsPindel;
    ifsPindel.open(strFilePath.c_str());
    string strLine = "";
    vBP.clear();
    St_BreakPoint stBP;
    while(!ifsPindel.eof())
    {
        getline(ifsPindel, strLine);

        if(strLine.find("##########") == string::npos) // do not find
            continue;

        //If find it
        getline(ifsPindel, strLine);
        int iStart = 0;
        int iLen = 0;
        int iEnd = 0;

        //For Chrom ID
        iStart = strLine.find("ChrID");
        iStart += 5;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.strChromID = strLine.substr(iStart, iLen);
        trim(stBP.strChromID);

        //For BP Start position
        iStart = strLine.find("BP", iEnd);
        iStart = iStart + 2;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.iStart = atoi(strLine.substr(iStart, iLen).c_str());

        //For ending position
        iStart = iEnd + 1;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.iEnd = atoi(strLine.substr(iStart, iLen).c_str());

        if(stBP.iStart > stBP.iEnd)
        {
            int iTmp = stBP.iStart;
            stBP.iStart = stBP.iEnd;
            stBP.iEnd = iTmp;
        }

        vBP.push_back(stBP);
    }
}

void ClsComparison::CompareStdSVWithPindel(string strFilePath, vector<St_SV>& vSvDEL)
{
     vector<St_BreakPoint> vBP;
     ParsePindel(strFilePath, vBP);
     const int MAXDIFF = 15;

     //Compare to check how many covered
     int iMiss = 0;
     for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
     {
        bool bFind =  false;
        for(vector<St_BreakPoint>::iterator subItr = vBP.begin(); subItr != vBP.end(); subItr++)
        {
            if(abs(itr->iPos - subItr->iStart) <= MAXDIFF ||
               abs(itr->iEnd - subItr->iEnd) <= MAXDIFF)
            {
                bFind = true;
            }
        }
        if(!bFind)
        {
            cout << "Miss SV: " << " --- " << "Len: " << to_string(itr->iEnd - itr->iPos) << " --- "
                 << "(" << to_string(itr->iPos) << ", "
                 << to_string(itr->iEnd) << ")" << endl;
            iMiss++;
        }
     }
     cout << endl << "-----------" << endl
          << "Miss: " << to_string(iMiss) << endl << "-----------" << endl;
     cout << "STD SV(DEL) Total   : " << to_string(vSvDEL.size()) << endl;
     cout << "Pindel SV(DEL) Tatal: " << to_string(vBP.size()) << endl;
}

void ClsComparison::CompareStdSVWithLumpy(string strFilePath, vector<St_SV>& vSvDEL)
{
    ClsVcf1000Genome* pClsVcf = new ClsVcf1000Genome();
    pClsVcf->ParseVcf(strFilePath);
    vector<St_SV> vSvLumpyDEL; // only pick out deletion
    pClsVcf->GetDeletion(vSvLumpyDEL);
    cout << "vSvDEL size: " << to_string(vSvLumpyDEL.size()) << endl;
    delete pClsVcf;
    pClsVcf = NULL;

    //Compare
    //Compare to check how many covered
    const int MAXDIFF = 15;
    int iMiss = 0;
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
       bool bFind =  false;
       for(vector<St_SV>::iterator subItr = vSvLumpyDEL.begin(); subItr !=vSvLumpyDEL.end(); subItr++)
       {
           if(abs(itr->iPos - subItr->iPos) <= MAXDIFF ||
              abs(itr->iEnd - subItr->iEnd) <= MAXDIFF)
           {
               bFind = true;
           }
       }
       if(!bFind)
       {
           cout << "Miss SV: " << " --- " << "Len: " << to_string(itr->iEnd - itr->iPos) << " --- "
                << "(" << to_string(itr->iPos) << ", "
                << to_string(itr->iEnd) << ")" << endl;
           iMiss++;
       }
    }
    cout << endl << "-----------" << endl
         << "Miss: " << to_string(iMiss) << endl << "-----------" << endl;
    cout << "STD SV(DEL) Total   : " << to_string(vSvDEL.size()) << endl;
    cout << "Lumpy SV(DEL) Tatal: " << to_string(vSvLumpyDEL.size()) << endl;
}


