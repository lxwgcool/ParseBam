#ifndef CLSPARSEBAM_H
#define CLSPARSEBAM_H

#include <string>
#include <vector>
#include "../../../../ShareLibrary/clsbasealgorithm.h"
#include "../../../../ShareLibrary/clsvcf1000genome.h"
#include "../../../../ShareLibrary/clsfastareader.h"
#include "../../../../ShareLibrary/clsblast.h"
#include "clsdrawimage.h"
#include <map>

using namespace std;

//enum En_ClipPart{cpLeft, cpRight, cpMax};
enum En_ClipType{ctSoft, ctHard, ctMax};
enum En_MapType{mtMat1End, mtMat2End, mtMat2Start, mtMat1Start, mtMax};

struct St_ClipReads
{
    int iReadsMapPos;
    string strName;
    string strQuerySeq;
    string strAlignSeq;
    vector<int> vPos;
    vector<string> vClipSeq;

    int iMateMapPos;
    string strMateName;
    string strMateQuerySeq;
    string strMateAlignSeq;
    vector<int> vMatePos;
    vector<string> vMateClipSeq;

    int iInsertSize;
    bool bFirstMate;
    bool bSecondMate;
    bool bMateMapped;

    St_ClipReads()
    {
        Clear();
    }

    void Clear()
    {
        iReadsMapPos = -1;
        strName = "";
        strQuerySeq = "";
        strAlignSeq = "";
        vPos.clear();
        vClipSeq.clear();

        iMateMapPos = -1;
        strMateName = "";
        strMateQuerySeq = "";
        strMateAlignSeq = "";
        vMatePos.clear();
        vMateClipSeq.clear();

        iInsertSize = -1;
        bFirstMate = false;
        bSecondMate = false;
        bMateMapped = false;
    }
};

struct St_BorderSCR //Border Soft-Cliped-Reads
{
    int iReadsMapPos;
    int iMatPos;
    int iClipPos;
    En_ClipPart enClipPart;
    string strClipSeq;
    string strQuerySeq;
    int iInsertSize;
    string strName;

    //For debug
    bool bSvContribute;
    int iDiscCountNum;

    En_MapType enMapType;
    bool bFirstMate;
    bool bSecondMate;
    bool bMateMapped;

    St_BorderSCR()
    {
        Clear();
    }

    void Clear()
    {
        iReadsMapPos = -1;
        iMatPos = -1;
        iClipPos = -1;
        enClipPart = cpMax;
        strClipSeq = " ";
        strQuerySeq = " ";
        iInsertSize = -1;
        strName = " ";

        //-->For debug
        bSvContribute = false;
        iDiscCountNum = 0;
        //<--

        enMapType = mtMax;
        bFirstMate = false;
        bSecondMate = false;
        bMateMapped = false;
    }

    bool IsLeftPair()
    {
        if(iReadsMapPos < iMatPos)
            return true;
        else
            return false;
    }

    bool IsRightPair()
    {
        if(iReadsMapPos > iMatPos)
            return true;
        else
            return false;
    }

    En_ClipPart GetClipPart()
    {
        if(this->iReadsMapPos == this->iClipPos)
            return cpLeft;
        else
            return cpRight;
    }
};

struct St_SvDelReads //support the deletion recorded in vcf
{
    St_SV stSV; // Which sv those reads related to
    vector<St_BorderSCR> vSVLeftReads; // all the clip are related to reads itself
    vector<St_BorderSCR> vSVRightReads;
    //int iSvLefAlignPos; //LEFT aligned pos should equal to SV RIGHT boundary position
    //int iSvRightAlignPos; //RIGHT aligned pos should equal to SV LEFT boundary position
    bool bSvLeftSupportRight;
    bool bSvRightSupportLeft;
    string strSvLeftAlignSeq;
    string strSvRightAlignSeq;

    void Clear()
    {
        stSV.Clear();
        vSVLeftReads.clear();
        vSVRightReads.clear();
        //iSvLefAlignPos = -1;
        //iSvRightAlignPos = -1;
        bSvLeftSupportRight = false;
        bSvRightSupportLeft = false;
        strSvLeftAlignSeq = "";
        strSvRightAlignSeq = "";
    }
};

//-->For discondent reads
struct  St_DiscordantReads
{
    //For left reads
    int iReadsPos;  //left side reads
    string strReadsAlignSeq;
    string strReadsClipSeq;
    bool bClip;

    //For right reads
    int iMatePos; //right side reads
    string strMateAlignSeq;
    string strMateClipSeq;
    bool bMateClip;

    //Discordant info
    int iInsertSize;

    St_DiscordantReads()
    {
        Clear();
    }

    void Clear()
    {
        //Left Reads
        iReadsPos = 0;
        strReadsAlignSeq = "";
        strReadsClipSeq = "";
        bClip = false;

        //Right reads
        iMatePos = 0;
        strMateAlignSeq = "";
        strMateClipSeq = "";
        bMateClip = false;

        //Discordant Info
        iInsertSize = 0;
    }
};


struct St_MultiMapReads
{
    string strName;
    vector<int> vPos;
    vector<En_ClipType> vClipType;
    bool bSoftClip;
    bool bHardClip;

    St_MultiMapReads()
    {
        Clear();
    }

    void Clear()
    {
        strName = "";
        vPos.clear();
        vClipType.clear();
        bSoftClip = false;
        bHardClip = false;
    }
};
//<--

// use the similar way of delly to make clusrering
struct St_DRGroup // Discordant reads group
{
    //-->
    vector<St_DiscordantReads> vDR;  //Left-most first, Right-most last
    vector<St_BorderSCR> vDRRangedSCR; // SCR dropped into DR range
    //<--

    // -->
    vector<St_BorderSCR> vStartSCR;  // Start Border of Clipped Reads
    vector<St_BorderSCR> vEndSCR; //End Border of Clipped Reads
    // <--
    int iStart; // valid range for SV
    int iEnd; // valid range for SV

    int iLeftBoundary; //For real left boundary
    int iRightBoundary; //For real right boundary

    St_DRGroup()
    {
        Clear();
    }

    void Clear()
    {
        vDR.clear();
        vDRRangedSCR.clear();
        vStartSCR.clear();
        vEndSCR.clear();
        iStart = 0;
        iEnd = 0;
        iLeftBoundary = 0;
        iRightBoundary = 0;
    }

    int GetLen()
    {
        return iRightBoundary - iLeftBoundary;
    }
};

struct St_GroupBound
{
    vector<St_BorderSCR> vSCR;
    bool bLeftSVBoundGood;
    bool bRightSVBoundGood;
    int iClipPos;

    St_GroupBound()
    {
        Clear();
    }

    void Clear()
    {
        vSCR.clear();
        bLeftSVBoundGood = false;
        bRightSVBoundGood = false;
        iClipPos = 0;
    }
};

class ClsParseBam
{
public:
    ClsParseBam();
    ~ClsParseBam();

public:
    //This is the raw bam file --> form the original raw reads -->
    void ReadBamFile(string strBamFilePath, vector<St_BorderSCR>& vBorderSCR);
    void ReadBamFileMergedAlignReads(string strBamFilePath,
                                     vector<St_SvDelReads>& vSvDelReads, vector<St_SV>& vSvDEL);

    void DebugSpeciSv(string strBamFilePath, St_SV& stSv);
    //<--
    void GetDELRelatedReads(vector<St_SV>& vSvDEL,
                            vector<St_BorderSCR>& vBorderSCR,
                            vector<St_SvDelReads>& vSvDelReads);
    void OutputSvDelReads(vector<St_SvDelReads>& vSvDelReads);

    //For extracting all pissible SV candidate
    void GetDiscordantReads(string strBamFilePath, vector<St_DiscordantReads>& vDiscdRreads);
    void GetClipReads(string strBamFilePath, vector<St_ClipReads>& vClipReads);
    void GetMultiPlacedMappedReads(string strBamFilePath);

    //clustering Discordant Reads ---> Cluster the reads which contains the same part of overlap -->
    void ClusterDS(vector<St_DiscordantReads>& vDiscdRreads, vector<St_SV>& vSvDEL); // clustering discordant reads    

    void ClusterDiscrodantReadsDelly(vector<St_DRGroup>& vDRGroup,
                                     vector<St_DiscordantReads>& vDiscdRreads,
                                     vector<St_BorderSCR>& vBorderSCR,
                                     vector<St_SV>& vSvDEL);

    void GetClusterBySCR(vector<St_BorderSCR>& vBorderSCR, vector<St_DRGroup>& vSCRGroup);
    void GetClusterBySCRTargetReads(vector<St_BorderSCR>& vBorderSCR, vector<St_DRGroup>& vSCRGroup);
    void UpgradeSCR(vector<St_BorderSCR>& vBorderSCR, vector<St_BorderSCR>& vNewSCR);

    void GetSCRGroup(vector<St_GroupBound>& vGroupBound, vector<St_DRGroup>& vSCRGroup);
    bool GetGroupFromBoundary(vector<St_BorderSCR>& vSCR, St_DRGroup& stSCRGroup);

    //<--

    void DrawPics(vector<St_DRGroup>& vDRGroup, string strRefPath, vector<St_SV>& vSvDEL);
    void DrawPicsForSingleGroup(St_DRGroup& stDRGroup, St_Fasta& vFasta, ClsDrawImage* pDrawImg, int iGroupIndex);
    void DrawPicsForSingleQueryRefSeq(string& strQuerySeq, string& strRefFa, ClsBlast* pBlast,
                                      cv::Mat& img, ClsDrawImage* pDrawImg, int iOffset, int iBaseLine);
};

#endif // CLSPARSEBAM_H
