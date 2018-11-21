#ifndef CLSDEBUG_H
#define CLSDEBUG_H

#include "../../../../ShareLibrary/clsmuscle.h"
#include "../../../../ShareLibrary/clsbasealgorithm.h"
#include "clsparsebam.h"

struct St_SCRSet
{
    vector<St_BorderSCR> vSet;
    unsigned int iAvgClipPos;
    unsigned int iAvgMapPos;
    bool bLeftYes;
    bool bRightYes;

    St_SCRSet()
    {
        Clear();
    }

    void Clear()
    {
        vSet.clear();
        iAvgClipPos = 0;
        iAvgMapPos = 0;
        bLeftYes = false;
        bRightYes = false;
    }
};

struct St_ClusterPair
{
    St_SCRSet stLeftSCRSet;
    St_SCRSet stRightSCRSet;
    int iCount;

    St_ClusterPair(): iCount(0)
    {
        Clear();
    }

    void Clear()

    {
        stLeftSCRSet.Clear();
        stRightSCRSet.Clear();
        iCount = 0;
    }
};

struct St_OrgReads
{
    string strName;
    string strSeq;
    int iMatePos;
    int iPos;
    int iInsertSize;

    //-->for debug
    int iDiscCountNum;
    //<--

    St_OrgReads()
    {
        Clear();
    }

    void Clear()
    {
       strName = "";
       strSeq = "";
       iMatePos = -1;
       iInsertSize = -1;
       iPos = -1;
       iDiscCountNum = 0;
    }
};

struct St_SvDiff
{
    int iSvStart;
    int iSvEnd;

    //Start Part
    St_OrgReads stSvStartReads;
    int iSvStartMinDiff;

    //End Part
    St_OrgReads stSvEndReads;
    int iSvEndMinDiff;

    St_SvDiff()
    {
        Clear();
    }

    void Clear()
    {
        iSvStartMinDiff = 1000000;
        iSvEndMinDiff = 1000000;
        iSvStart = -1;
        iSvEnd = -1;
        stSvStartReads.Clear();
        stSvStartReads.Clear();
    }

    string GetStrand()
    {
        if(iSvStart == -1 || iSvEnd == -1)
            return "*";
        else if(iSvStart == iSvEnd)
            return "=";
        else if(iSvStart < iSvEnd)
            return "+";
        else if(iSvStart > iSvEnd)
            return "-";
        else
            return "?";
    }

    int GetDiff()
    {
        if(iSvStartMinDiff > iSvEndMinDiff)
            return iSvEndMinDiff;
        else
            return iSvStartMinDiff;
    }

    int GetDiffSum() //For discordant reads
    {
        return iSvStartMinDiff + iSvEndMinDiff;
    }

    int GetInsertSize()
    {
        if(iSvStartMinDiff > iSvEndMinDiff)
            return stSvEndReads.iInsertSize;
        else
            return stSvStartReads.iInsertSize;
    }

    string GetSvHitSide()
    {
        if(iSvStartMinDiff < iSvEndMinDiff)
            return "Start";
        else
            return "End";
    }

    int GetDiscCountNum()
    {
        if(iSvStartMinDiff > iSvEndMinDiff)
            return stSvEndReads.iDiscCountNum;
        else
            return stSvStartReads.iDiscCountNum;
    }
};

struct St_SvReadsSet
{
    St_SV stSv;
    vector<St_BorderSCR> vLeftBoard; //small position
    vector<St_BorderSCR> vRightBoard; //large position

    St_SvReadsSet()
    {
        Clear();
    }

    void Clear()
    {
        stSv.Clear();
        vLeftBoard.clear();
        vRightBoard.clear();
    }
};

struct St_DiscClipGroup
{
    St_BorderSCR stClipReads;
    //vector<St_DiscordantReads*> vDiscdReads;
    int iCount;

    void Clear()
    {
        stClipReads.Clear();
        iCount = 0;
        //vDiscdReads.clear();
    }
};

class ClsDebug
{
public:
    ClsDebug();

    //void GetDELRelatedReads(vector<St_SV>& vSvDEL, vector<St_BorderSCR>& vBorderSCR);
    void CollectRawSvCandidate(string strBamFilePath, vector<St_SV>& vSvDEL);
    void DebugCheckSVByDiscordantReads(vector<St_DiscordantReads>& vDiscdRreads, vector<St_SV>& vSvDEL);
    void DebugCheckSVByClipReads(vector<St_ClipReads>& vClipReads, vector<St_SV>& vSvDEL);

    void CheckHitStatusByClipReads(vector<St_SV>& vSvDEL, vector<St_BorderSCR>& vBorderSCR);
    void GetChipReadsBySV(vector<St_SV>& vSvDEL, vector<St_BorderSCR>& vBorderSCR, bool bBlastConfirm=false);
    void CheckSVClipReadsByBlast(vector<St_SV>& vSvDEL, vector<St_BorderSCR>& vBorderSCR);

    void CheckHitStatusByDiscordantReads(vector<St_SV>& vSvDEL, vector<St_DiscordantReads>& vDiscdRreads);

    void CheckStdSvCount(vector<St_DiscordantReads>& vDiscdRreads, vector<St_SV>& vSvDEL);

    //Check
    void GetDiscClipGroup(vector<St_BorderSCR>& vClipReads, vector<St_DiscordantReads>& vDiscdRreads);

    void FilterSCReadsByDiscReads(vector<St_BorderSCR>& vBorderSCR, vector<St_DiscordantReads>& vDiscdRreads);

    //Save the corresponding reads
    string SaveClipReads(vector<St_BorderSCR>& vBorderSCR, string strFilePath="");

    //Do local assembly
    void LocalAssemblyByVelvet(string strReadsPath);

    //Use K-mer to get the candidate group -->
    void ClusterClipReads(vector<St_BorderSCR>& vBorderSCR, vector<St_SV>& vSvDEL);
    //<--

private:
    void UpdateSCRSet(St_SCRSet& stSRCSetOrg, St_SCRSet& stSRCSetTmp);
};

#endif // CLSDEBUG_H
