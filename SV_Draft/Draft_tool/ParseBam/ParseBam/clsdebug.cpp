#include "clsdebug.h"
#include <algorithm>
#include "../../../../ShareLibrary/clsvelvet.h"
#include "../../../../ShareLibrary/clsblast.h"
#include <unordered_map>
#include "../../../../ShareLibrary/clsbasealgorithm.h"
#include "../../../../ShareLibrary/clskmeralgorithm.h"

const int OFFSET = 10;
const char* ArrayMapType[mtMax] = {"Mat1End", "Mat2End", "Mat2Start", "Mat1Start"};
const char* ArrayClipPart[cpMax] = {"Start", "End", "Any"};

//Cluster discordant reads --> For "for_each"
struct AddVal
{
    int iVal;
    AddVal(int v) : iVal(v)
    {}

    void operator()(int &elem) const
    {
        elem += iVal;
    }
};

ClsDebug::ClsDebug()
{
}

bool sort_svdiff(St_SvDiff stV1, St_SvDiff stV2)
{
    if(stV1.GetDiff() > stV2.GetDiff())
        return true;
    else
        return false;
}

void ClsDebug::CheckHitStatusByClipReads(vector<St_SV>& vSvDEL, vector<St_BorderSCR>& vBorderSCR)
{
//    //--> For debug: only keep the target SV -->
//    for(vector<St_SV>::iterator itr = vSvDEL.end() - 1; itr >= vSvDEL.begin(); itr--)
//    {
//        if( (itr->iPos  == 190950 && itr->iEnd == 191169) ||
//            (itr->iPos  == 27657262 && itr->iEnd == 27657596) )
//        {}
//        else
//        {
//            vSvDEL.erase(itr);
//            continue;
//        }
//    }
//    //<--

    //const int DIFFMAX = 50;
    //Check if vBorderSCR already contain all of standard SV -->
    int iMissNum = 0;
    cout << endl << "-----" << endl;
    //Sort first
    vector<St_SvDiff> vSvDiff;
    St_SvDiff stSvDiff;
    const int CLOSLEN = 60;
    int iStartNum = 0;
    int iEndNum = 0;

    //--> Check the real reads support each different SV
    vector<St_SvReadsSet> vSvReadsSet;
    St_SvReadsSet stSvReadsSet;
    //<--

    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
//        if(itr->iPos == 4972926 && itr->iEnd == 4973837)
//        {
//            cout << "SV Length: " << to_string(itr->iEnd - itr->iPos) << endl;
//            pClsParseBam->DebugSpeciSv(stConfig.strBamFile, *itr);
//        }

        bool bFind = false;
        stSvDiff.Clear();
        stSvReadsSet.Clear();
        stSvReadsSet.stSv = *itr;
        for(vector<St_BorderSCR>::iterator subItr = vBorderSCR.begin(); subItr != vBorderSCR.end(); subItr++)
        {
            //if(subItr->enClipPart == cpLeft) //这里是判断左边被截去的情况
            {
                if(abs(itr->iEnd - subItr->iReadsMapPos) < stSvDiff.iSvEndMinDiff) // this is right border
                {
                    stSvDiff.iSvEndMinDiff = abs(itr->iEnd - subItr->iReadsMapPos);
                    stSvDiff.stSvEndReads.iInsertSize = subItr->iInsertSize;
                    stSvDiff.stSvEndReads.iMatePos = subItr->iMatPos;
                    stSvDiff.stSvEndReads.iPos = subItr->iReadsMapPos;
                    stSvDiff.stSvEndReads.strSeq = subItr->strQuerySeq;
                    stSvDiff.stSvEndReads.strName = subItr->strName;
                    stSvDiff.stSvEndReads.iDiscCountNum = subItr->iDiscCountNum;
                }

                if(abs(itr->iEnd - subItr->iReadsMapPos) < CLOSLEN)
                {
                    bFind = true;
                    subItr->bSvContribute = true;
                    stSvReadsSet.vRightBoard.push_back(*subItr);
                }

//                if(abs(itr->iEnd - subItr->iReadsMapPos) < CLOSLEN)
//                {
//                    subItr->bSvContribute = true;
//                }

            }
            //else if(subItr->enClipPart == cpRight)
            {
                if(abs(itr->iPos - subItr->iClipPos) < stSvDiff.iSvStartMinDiff) // this is left border
                {
                    stSvDiff.iSvStartMinDiff = abs(itr->iPos - subItr->iClipPos);
                    stSvDiff.stSvStartReads.iInsertSize = subItr->iInsertSize;
                    stSvDiff.stSvStartReads.iMatePos = subItr->iMatPos;
                    stSvDiff.stSvStartReads.iPos = subItr->iReadsMapPos;
                    stSvDiff.stSvStartReads.strSeq = subItr->strQuerySeq;
                    stSvDiff.stSvStartReads.strName = subItr->strName;
                    stSvDiff.stSvStartReads.iDiscCountNum = subItr->iDiscCountNum;
                }

                if(abs(itr->iPos - subItr->iClipPos) < CLOSLEN)
                {
                    bFind = true;
                    subItr->bSvContribute = true;
                    stSvReadsSet.vLeftBoard.push_back(*subItr);
                }

//                if(abs(itr->iPos - subItr->iClipPos) < DIFFMAX)
//                {
//                    subItr->bSvContribute = true;
//                }
            }
            //if(bFind)
            //    break;
        }

        //--> Save reads set for each SV -->
        vSvReadsSet.push_back(stSvReadsSet);
        //<--

        if(!bFind)
            iMissNum++;
        else
        {
            if(stSvDiff.GetSvHitSide() == "End")
                iEndNum++;
            else if(stSvDiff.GetSvHitSide() == "Start")
                iStartNum++;
        }

        stSvDiff.iSvStart = itr->iPos;
        stSvDiff.iSvEnd = itr->iEnd;
        vSvDiff.push_back(stSvDiff);
    }

    //--> Output the supported reads for each different SV
    ofstream ofsReadSets;
    ofsReadSets.open("./SvReadsSet.txt");
    int iBothSupport = 0;
    for(vector<St_SvReadsSet>::iterator itr = vSvReadsSet.begin(); itr != vSvReadsSet.end(); itr++)
    {
        if(itr->vLeftBoard.empty() || itr->vRightBoard.empty())
        {
            continue;
        }
        iBothSupport++;

        ofsReadSets << "Sv: <" << IntToStr(itr->stSv.iPos) << ", " << IntToStr(itr->stSv.iEnd) << ">" << endl;

        //Left Support
        ofsReadSets << "----- Left Border -----" << endl;
        int iMinLeftDiff = 100000;
        for(vector<St_BorderSCR>::iterator subItr = itr->vLeftBoard.begin(); subItr != itr->vLeftBoard.end(); subItr++)
        {
            int iDiff = abs(itr->stSv.iPos - subItr->iClipPos);
            ofsReadSets << subItr->strQuerySeq
                        << " --- " << to_string(iDiff)
                        << " --- " << "Mapping Pos: " << IntToStr(subItr->iReadsMapPos)
                        << " --- " << "Clip Pos: " << IntToStr(subItr->iClipPos)
                        << " --- " << "Insert_Size: " << to_string(subItr->iInsertSize)
                        << " --- " << "Mate_Mapped: " << string(subItr->bMateMapped ? "Yes" : "No")
                        << " --- " << "Mate_Pos: " << to_string(subItr->iMatPos)
                        << " --- " << "ArrayClipPart: " << string(ArrayClipPart[subItr->enClipPart])
                        << " --- " << "Is First Mate: " << string(subItr->bFirstMate ? "Yes" : "No")
                        << " --- " << "Is Second Mate: " << string(subItr->bSecondMate ? "Yes" : "No") << endl;

            if(iMinLeftDiff > iDiff)
                iMinLeftDiff = iDiff;
        }
        ofsReadSets << "Min Left Diff: " << to_string(iMinLeftDiff) << endl;

        //Right Support
        int iMinRightDiff = 100000;
        ofsReadSets << "----- Right Border -----" << endl;
        for(vector<St_BorderSCR>::iterator subItr = itr->vRightBoard.begin(); subItr != itr->vRightBoard.end(); subItr++)
        {
            int iDiff = abs(itr->stSv.iEnd - subItr->iReadsMapPos);
            ofsReadSets << subItr->strQuerySeq
                        << "  --  " << to_string(iDiff)
                        << "  --  " << "Mapping Pos: " << IntToStr(subItr->iReadsMapPos)
                        << "  --  " << "Clip Pos: " << IntToStr(subItr->iClipPos)
                        << " -- " << "Insert_Size: " << to_string(subItr->iInsertSize)
                        << " -- " << "Mate_Mapped: " << string(subItr->bMateMapped ? "Yes" : "No")
                        << " -- " << "Mate_Pos: " << to_string(subItr->iMatPos)
                        << " -- " << "ArrayClipPart: " << string(ArrayClipPart[subItr->enClipPart])
                        << " -- " << "Is First Mate: " << string(subItr->bFirstMate ? "Yes" : "No")
                        << " -- " << "Is Second Mate: " << string(subItr->bSecondMate ? "Yes" : "No") << endl;

            if(iMinRightDiff > iDiff)
                iMinRightDiff = iDiff;
        }
        ofsReadSets << "Min Right Diff: " << to_string(iMinRightDiff) << endl;

        ofsReadSets << endl;
    }
    ofsReadSets << to_string(iBothSupport) << endl;
    ofsReadSets.close();
    //<--

    //Set who contributes to each SV -->
    for(vector<St_SvDiff>::iterator itr = vSvDiff.begin(); itr != vSvDiff.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_BorderSCR>::iterator subItr = vBorderSCR.begin(); subItr != vBorderSCR.end(); subItr++)
        {
            if(itr->GetSvHitSide() == "End")
            {
                //the end reads
                if(itr->stSvEndReads.iPos == subItr->iReadsMapPos &&
                   itr->stSvEndReads.iMatePos == subItr->iMatPos)
                {
                    subItr->bSvContribute = true;
                    bFind = true;
                }
            }
            else if(itr->GetSvHitSide() == "Start")
            {
                //for the start side reads
                //the end reads
                if(itr->stSvStartReads.iPos == subItr->iReadsMapPos &&
                   itr->stSvStartReads.iMatePos == subItr->iMatPos)
                {
                    subItr->bSvContribute = true;
                    bFind = true;
                }
            }
            if(bFind)
                break;
        }
    }
    //<--

    sort(vSvDiff.begin(), vSvDiff.end(), sort_svdiff);
    int iIndex = 1;
    int iStartGood = 0;
    int iEndGood = 0;
    int iBothGood = 0;
    for(vector<St_SvDiff>::iterator itr = vSvDiff.begin(); itr != vSvDiff.end(); itr++, iIndex++)
    {
        cout << to_string(iIndex) << "\t"
             << to_string(itr->GetDiff())
             << " ---- "
             << "(" << to_string(itr->iSvStart) << ", " << to_string(itr->iSvEnd) << ")"
             << " ---- " << to_string(itr->GetDiscCountNum())
             << " ---- " << itr->GetSvHitSide()
             << " ---- " << to_string(itr->GetInsertSize())
             << " .... (" << to_string(itr->iSvStartMinDiff) << ", " << to_string(itr->iSvEndMinDiff) << ")"
             << " .... " << "Start(" << to_string(itr->stSvStartReads.iInsertSize)
                                     << ", " << itr->stSvStartReads.strName << ")"
             << " .... " << "End(" << to_string(itr->stSvEndReads.iInsertSize)
                                   << ", " << itr->stSvEndReads.strName << ")" << endl;

        if(itr->iSvStartMinDiff < CLOSLEN)
            iStartGood++;
        if(itr->iSvEndMinDiff < CLOSLEN)
            iEndGood++;
        if(itr->iSvStartMinDiff < CLOSLEN && itr->iSvEndMinDiff < CLOSLEN)
            iBothGood++;
    }

    cout <<"-----" << endl;

    cout << "vSvDEL for Sample: " << to_string(vSvDEL.size()) << endl;
    cout << "vBorderSCR Size  : " << to_string(vBorderSCR.size()) << endl;
    cout << "Miss Num         : " << to_string(iMissNum) << endl;
    cout << "Hit Num          : " << to_string(vSvDEL.size() - iMissNum) << endl;
    cout << "  Support Start  : " << to_string(iStartNum) << endl;
    cout << "  Support End    : " << to_string(iEndNum) << endl;

    cout << endl << " ---- " << endl;
    cout << "Start Good : " << to_string(iStartGood) << endl;
    cout << "End Good   : " << to_string(iEndGood) << endl;
    cout << "Both Good  : " << to_string(iBothGood) << endl;

    //Cout the number of contribute
    int iContribute = 0;
    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr != vBorderSCR.end(); itr++)
    {
        if(itr->bSvContribute)
            iContribute++;
    }
    cout << endl << "------" << endl;
    cout << IntToStr(iContribute) << endl;
    cout << "------" << endl;
}

struct St_SVClipReads
{
    St_SV stSv;
    vector<St_BorderSCR> vSvStartClipReads;
    vector<St_BorderSCR> vSvEndClipReads;

    bool bBlastResult;

    void Clear()
    {
        stSv.Clear();
        vSvStartClipReads.clear();
        vSvEndClipReads.clear();
        bBlastResult = false;
    }
};

void ClsDebug::GetChipReadsBySV(vector<St_SV>& vSvDEL, vector<St_BorderSCR>& vBorderSCR,
                                bool bBlastConfirm)
{
    const int CLOSLEN = 80;
    vector<St_SVClipReads> vSvClipReads;
    St_SVClipReads stSvClipReads;

    for(vector<St_SV>::iterator itrSv = vSvDEL.begin(); itrSv != vSvDEL.end(); itrSv++)
    {
        stSvClipReads.Clear();
        stSvClipReads.stSv = *itrSv;
        for(vector<St_BorderSCR>::iterator itrReads = vBorderSCR.begin(); itrReads != vBorderSCR.end(); itrReads++)
        {
            //Check Sv Start
            if(abs(itrSv->iPos - itrReads->iClipPos) < CLOSLEN)
                stSvClipReads.vSvStartClipReads.push_back(*itrReads);

            //Check Sv End
            if(abs(itrSv->iEnd - itrReads->iReadsMapPos) < CLOSLEN)
                stSvClipReads.vSvEndClipReads.push_back(*itrReads);
        }
        vSvClipReads.push_back(stSvClipReads);
    }

    if(!bBlastConfirm)
    {}
    else
    {
        //do blast confirm
        //Step 1: build comparison group
        ClsBlast* pBlast = new ClsBlast();
        const int BLASTLEN = 60;
        int iBlastGoodNum = 0;
        for(vector<St_SVClipReads>::iterator itr = vSvClipReads.begin(); itr != vSvClipReads.end(); itr++)
        {
            if(itr->bBlastResult)
                continue;

            if(itr->vSvStartClipReads.empty() || itr->vSvEndClipReads.empty())
            {
                itr->bBlastResult = false;
                continue;
            }

            for(vector<St_BorderSCR>::iterator itrStartReads = itr->vSvStartClipReads.begin();
                itrStartReads != itr->vSvStartClipReads.end(); itrStartReads++)
            {
                //For start: Last 60
                string strSvStart = itrStartReads->strQuerySeq; //itrStartReads->strQuerySeq.substr(itrStartReads->strQuerySeq.length() - BLASTLEN, BLASTLEN);
                string strRefFa = pBlast->CreatFaFile("RefSeq", strSvStart, ftRef);

                for(vector<St_BorderSCR>::iterator itrEndReads = itr->vSvEndClipReads.begin();
                    itrEndReads != itr->vSvEndClipReads.end(); itrEndReads++)
                {
                    //For End: first 60
                    string strSvEnd = itrEndReads->strClipSeq; //itrEndReads->strQuerySeq.substr(0, BLASTLEN);
                    string strQueryFa = pBlast->CreatFaFile("QuerySeq", strSvEnd, ftQuery);

                    int iAlignNum = pBlast->TwoSeqBlast(strQueryFa, strRefFa, false);
                    if(iAlignNum > 15)
                    {
                        itr->bBlastResult = true;
                        break;
                    }
                }
                if(itr->bBlastResult)
                {
                    iBlastGoodNum++;
                    break;
                }
            }

            if(itr->bBlastResult)
                continue;

            //Do it again: this time the ending side of cliped reads i reference and the start side of reads is query sequence  -------->
            for(vector<St_BorderSCR>::iterator itrEndReads = itr->vSvEndClipReads.begin();
                itrEndReads != itr->vSvEndClipReads.end(); itrEndReads++)
            {
                //For start: Last 60
                string strSvEnd = itrEndReads->strQuerySeq; //itrStartReads->strQuerySeq.substr(itrStartReads->strQuerySeq.length() - BLASTLEN, BLASTLEN);
                string strRefFa = pBlast->CreatFaFile("RefSeq", strSvEnd, ftRef);

                for(vector<St_BorderSCR>::iterator itrStartReads = itr->vSvStartClipReads.begin();
                    itrStartReads != itr->vSvEndClipReads.end(); itrStartReads++)
                {
                    //For End: first 60
                    string strSvStart = itrEndReads->strClipSeq; //itrEndReads->strQuerySeq.substr(0, BLASTLEN);
                    string strQueryFa = pBlast->CreatFaFile("QuerySeq", strSvStart, ftQuery);

                    int iAlignNum = pBlast->TwoSeqBlast(strQueryFa, strRefFa, false);
                    if(iAlignNum > 15)
                    {
                        itr->bBlastResult = true;
                        break;
                    }
                }
                if(itr->bBlastResult)
                {
                    iBlastGoodNum++;
                    break;
                }
            }
            //<---------------
        }
        //output statistic result -->
        iBlastGoodNum = 0;
        for(vector<St_SVClipReads>::iterator itr = vSvClipReads.begin(); itr != vSvClipReads.end(); itr++)
        {
            if(itr->bBlastResult)
                iBlastGoodNum++;
        }
        cout << endl << "Blast_Good_Num: " << to_string(iBlastGoodNum) << endl;
        delete pBlast;
        pBlast = NULL;
    }
    //<--

    //Output --> Go
    ofstream ofsSvClipReads;
    ofsSvClipReads.open("./SvClipReads.txt");
    int iBoth = 0;
    int iNone = 0;
    int iStartOnly = 0;
    int iEndOnly = 0;
    for(vector<St_SVClipReads>::iterator itr = vSvClipReads.begin(); itr != vSvClipReads.end(); itr++)
    {
        //1: output current SV
        ofsSvClipReads << "Std SV: (" << to_string(itr->stSv.iPos) << ", "
                       << to_string(itr->stSv.iEnd) << ")" << " --- ";
        if(itr->bBlastResult)
            ofsSvClipReads << "Pass_Blast --- ";
        else
            ofsSvClipReads << "Fail_Blast --- ";

        if(!itr->vSvStartClipReads.empty() && !itr->vSvEndClipReads.empty())
        {
            ofsSvClipReads << "Both_Clip" << endl;
            iBoth++;
        }
        else if(!itr->vSvStartClipReads.empty() && itr->vSvEndClipReads.empty())
        {
            ofsSvClipReads << "Start_Only_Clip" << endl;
            iStartOnly++;
        }
        else if(itr->vSvStartClipReads.empty() && !itr->vSvEndClipReads.empty())
        {
            ofsSvClipReads << "End_Only_Clip" << endl;
            iEndOnly++;
        }
        else
        {
            ofsSvClipReads << "None_Clip" << endl;
            iNone++;
        }

        //2: output Support Reads -- Sv Start
        ofsSvClipReads << "--- Start --- " << endl;
        for(vector<St_BorderSCR>::iterator itrReads = itr->vSvStartClipReads.begin(); itrReads != itr->vSvStartClipReads.end(); itrReads++)
        {
            ofsSvClipReads << itrReads->strQuerySeq << endl;
        }
        //3: output Support Reads -- Sv End
        ofsSvClipReads << "--- End --- " << endl;
        for(vector<St_BorderSCR>::iterator itrReads = itr->vSvEndClipReads.begin(); itrReads != itr->vSvEndClipReads.end(); itrReads++)
        {
            ofsSvClipReads << itrReads->strQuerySeq << endl;
        }
        ofsSvClipReads << endl;
    }

    ofsSvClipReads << endl << "**********" << endl
         << "Both      : " << to_string(iBoth) << endl
         << "None      : " << to_string(iNone) << endl
         << "Start_Only: " << to_string(iStartOnly) << endl
         << "End_Only  : " << to_string(iEndOnly) << endl;

    ofsSvClipReads.close();
}

void ClsDebug::CheckSVClipReadsByBlast(vector<St_SV>& vSvDEL, vector<St_BorderSCR>& vBorderSCR)
{
    const int CLOSLEN = 80;
    vector<St_SVClipReads> vSvClipReads;
    St_SVClipReads stSvClipReads;

    for(vector<St_SV>::iterator itrSv = vSvDEL.begin(); itrSv != vSvDEL.end(); itrSv++)
    {
        stSvClipReads.Clear();
        stSvClipReads.stSv = *itrSv;
        for(vector<St_BorderSCR>::iterator itrReads = vBorderSCR.begin(); itrReads != vBorderSCR.end(); itrReads++)
        {
            //Check Sv Start
            if(abs(itrSv->iPos - itrReads->iClipPos) < CLOSLEN)
                stSvClipReads.vSvStartClipReads.push_back(*itrReads);

            //Check Sv End
            if(abs(itrSv->iEnd - itrReads->iReadsMapPos) < CLOSLEN)
                stSvClipReads.vSvEndClipReads.push_back(*itrReads);
        }
        vSvClipReads.push_back(stSvClipReads);
    }

    //
}

void ClsDebug::CheckHitStatusByDiscordantReads(vector<St_SV>& vSvDEL, vector<St_DiscordantReads>& vDiscdRreads)
{
    //Check if vBorderSCR already contain all of standard SV -->
    int iMissNum = 0;
    cout << endl << "-----" << endl;

    cout << "Total Discordant reads number: " << to_string(vDiscdRreads.size()) << endl;
    int iContailClipNum = 0;
    int iInsertSize0Num = 0;
    bool bDoOneTime = true;
    vector<St_SvDiff> vSvDiff;
    St_SvDiff stSvDiff;

    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        bool bFind = false;
        bool bContainClip = false;
        for(vector<St_DiscordantReads>::iterator itrDiscdReads = vDiscdRreads.begin();
            itrDiscdReads != vDiscdRreads.end(); itrDiscdReads++)
        {
            if(itrDiscdReads->iInsertSize == 0)
            {
                if(bDoOneTime)
                    iInsertSize0Num++;
                //continue;
            }
            //Check if it is in the range
            if( (itr->iPos > (itrDiscdReads->iReadsPos - OFFSET)) &&
                (itr->iEnd < (itrDiscdReads->iMatePos + OFFSET)))
            {
                bFind = true;
                if(!bContainClip)
                {
                    if(itrDiscdReads->bClip || itrDiscdReads->bMateClip)
                    {
                        iContailClipNum++;
                        bContainClip = true;
                    }
                }
                int iDiffSum = abs(itr->iPos - itrDiscdReads->iReadsPos) + abs(itr->iEnd - itrDiscdReads->iMatePos);
                //1: Try to Check if this SV has been recorded
                bool bSvRecord = false;
                for(vector<St_SvDiff>::iterator itrSvDiff = vSvDiff.begin();
                    itrSvDiff != vSvDiff.end(); itrSvDiff++)
                {
                    if(itr->iPos == itrSvDiff->iSvStart &&
                       itr->iEnd == itrSvDiff->iSvEnd)
                    {
                        //找到了  --> update
                        bSvRecord = true;
                        if(iDiffSum < itrSvDiff->GetDiffSum())
                        {
                            //update start reads
                            itrSvDiff->stSvStartReads.iPos = itrDiscdReads->iReadsPos;
                            itrSvDiff->stSvStartReads.iInsertSize = itrDiscdReads->iInsertSize;
                            itrSvDiff->iSvStartMinDiff = abs(itr->iPos - itrDiscdReads->iReadsPos);

                            //update end reads
                            itrSvDiff->stSvEndReads.iPos = itrDiscdReads->iMatePos;
                            itrSvDiff->iSvEndMinDiff = abs(itr->iPos - itrDiscdReads->iMatePos);
                        }
                        break;
                    }
                }

                if(!bSvRecord) //是一个新的
                {
                    stSvDiff.Clear();
                    stSvDiff.stSvStartReads.iPos = itrDiscdReads->iReadsPos;
                    stSvDiff.stSvStartReads.iInsertSize = itrDiscdReads->iInsertSize;
                    stSvDiff.iSvStartMinDiff = abs(itr->iPos - itrDiscdReads->iReadsPos);

                    stSvDiff.stSvEndReads.iPos = itrDiscdReads->iMatePos;
                    stSvDiff.iSvEndMinDiff = abs(itr->iPos - itrDiscdReads->iMatePos);

                    stSvDiff.iSvStart = itr->iPos;
                    stSvDiff.iSvEnd = itr->iEnd;

                    vSvDiff.push_back(stSvDiff);
                }
            }
        }

        if(!bFind)
        {
            cout <<"(" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ")" << endl;
            iMissNum++;
        }
        bDoOneTime = false;
    }

    cout << endl << "-----" << endl;
    cout << "Total SV        : " << to_string(vSvDEL.size()) << endl;
    cout << "Missing Number  : " << to_string(iMissNum) << endl;
    cout << "Hit Number      : " << to_string(vSvDEL.size() - iMissNum) << endl;
    cout << "Contail Clip Num: " << to_string(iMissNum) << endl;

    cout << endl << "----- Hit Sv Detial -----" << endl;

    int iAbnormalDiffNum = 0;
    for(vector<St_SvDiff>::iterator itr = vSvDiff.begin(); itr != vSvDiff.end(); itr++)
    {
        cout <<"SV: (" << to_string(itr->iSvStart) << ", " << to_string(itr->iSvEnd) << ")"
             << " --- SV Len: " << to_string(itr->iSvEnd - itr->iSvStart)
             << " --- Insert Size: " << to_string(itr->stSvStartReads.iInsertSize)
             << " --- Diff(" << to_string(itr->iSvStartMinDiff) << ", " << to_string(itr->iSvEndMinDiff) << ")"
             << " --- Reads Pos(" << to_string(itr->stSvStartReads.iPos)
                                  << ", " << to_string(itr->stSvEndReads.iPos) << ")" << endl;
        if(itr->GetDiffSum() > 2000)
            iAbnormalDiffNum++;
    }
    cout << endl << "Abnormal Diff Sum (>2000): " << to_string(iAbnormalDiffNum) << endl;
    cout << endl << "Insert Size 0 Num        : " << to_string(iInsertSize0Num) << endl;
}

void ClsDebug::CollectRawSvCandidate(string strBamFilePath, vector<St_SV>& vSvDEL)
{
    ClsParseBam* pClsParseBam = new ClsParseBam();

    //1: Collect based on discordant reads
    vector<St_DiscordantReads> vDiscdRreads; //DiscdRreads: discordant reads
    pClsParseBam->GetDiscordantReads(strBamFilePath, vDiscdRreads);
    pClsParseBam->ClusterDS(vDiscdRreads, vSvDEL);

    //2: Collect based on softclip reads
    vector<St_ClipReads> vClipReads; //We record their mates as well
    //pClsParseBam->GetClipReads(strBamFilePath, vClipReads);

    cout << "DiscordantReads Num: " << to_string(vDiscdRreads.size()) << endl;
    cout << "Clip Reads Num     : " << to_string(vClipReads.size()) << endl;

    //Debug
    //DebugCheckSVByDiscordantReads(vDiscdRreads, vSvDEL);
    //DebugCheckSVByClipReads(vClipReads, vSvDEL);

    //Clustering those reads


    //Get the final raw candidates (we do not require accurate breakpoint at this time)

    //Release Pointer
    delete pClsParseBam;
    pClsParseBam = NULL;
}

void ClsDebug::DebugCheckSVByDiscordantReads(vector<St_DiscordantReads>& vDiscdRreads, vector<St_SV>& vSvDEL)
{
    //Check if vBorderSCR already contain all of standard SV -->
    int iMissNum = 0;
    cout << endl << "-----" << endl;
    //Sort first
    vector<St_SvDiff> vSvDiff;
    St_SvDiff stSvDiff;

    cout << "Total Cliped reads number: " << to_string(vDiscdRreads.size()) << endl;

    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        bool bFind = false;
        stSvDiff.Clear();
        for(vector<St_DiscordantReads>::iterator subItr = vDiscdRreads.begin();
            subItr != vDiscdRreads.end(); subItr++)
        {
            //Check if it is in the range
            if(itr->iPos >= subItr->iReadsPos && itr->iEnd <= subItr->iMatePos)
            {
                if(abs(itr->iPos - (int)(subItr->iReadsPos + subItr->strReadsAlignSeq.length())) < stSvDiff.iSvStartMinDiff)
                {
                    stSvDiff.iSvStartMinDiff = abs(itr->iPos - (int)(subItr->iReadsPos + subItr->strReadsAlignSeq.length()));
                    stSvDiff.stSvStartReads.iInsertSize = subItr->iInsertSize;
                }

                if(abs(itr->iEnd - subItr->iMatePos) < stSvDiff.iSvEndMinDiff)
                {
                    stSvDiff.iSvEndMinDiff = abs(itr->iEnd - subItr->iMatePos);
                    stSvDiff.stSvEndReads.iInsertSize = subItr->iInsertSize;
                }

                bFind = true;
            }
        }

        if(!bFind)
            iMissNum++;

        stSvDiff.iSvStart = itr->iPos;
        stSvDiff.iSvEnd = itr->iEnd;
        vSvDiff.push_back(stSvDiff);
    }

    sort(vSvDiff.begin(), vSvDiff.end(), sort_svdiff);
    for(vector<St_SvDiff>::iterator itr = vSvDiff.begin(); itr != vSvDiff.end(); itr++)
    {
        cout << to_string(itr->GetDiff()) << "   ----   "
             <<"(" << to_string(itr->iSvStart) << ", " << to_string(itr->iSvEnd) << ")"
             << "   ----   SV_Len: " << to_string(itr->iSvEnd - itr->iSvStart + 1)
             << "   ----   " << to_string(itr->GetInsertSize())
             << "   ----   " << itr->GetStrand() << endl;
    }

    cout <<"-----" << endl;
    cout << "Total SV      : " << to_string(vSvDEL.size()) << endl;
    cout << "Missing Number: " << to_string(iMissNum) << endl;

}

void ClsDebug::DebugCheckSVByClipReads(vector<St_ClipReads>& vClipReads, vector<St_SV>& vSvDEL)
{
    //Check if vBorderSCR already contain all of standard SV -->
    int iMissNum = 0;
    cout << endl << "-----" << endl;
    //Sort first
    vector<St_SvDiff> vSvDiff;
    St_SvDiff stSvDiff;

    cout << "Total Cliped reads number: " << to_string(vClipReads.size()) << endl;

    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        bool bFind = false;
        stSvDiff.Clear();
        for(vector<St_ClipReads>::iterator subItr = vClipReads.begin(); subItr != vClipReads.end(); subItr++)
        {
            if(abs(itr->iPos - (int)(subItr->iReadsMapPos + subItr->strAlignSeq.length())) < stSvDiff.iSvStartMinDiff)
            {
                stSvDiff.iSvStartMinDiff = abs(itr->iPos - (int)(subItr->iReadsMapPos + subItr->strAlignSeq.length()));
                stSvDiff.stSvStartReads.iInsertSize = subItr->iInsertSize;
            }

            if(abs(itr->iEnd - subItr->iMateMapPos) < stSvDiff.iSvEndMinDiff)
            {
                stSvDiff.iSvEndMinDiff = abs(itr->iEnd - subItr->iMateMapPos);
                stSvDiff.stSvEndReads.iInsertSize = subItr->iInsertSize;
            }

            if(stSvDiff.GetDiff() < 100)
                bFind = true;
        }

        if(!bFind)
            iMissNum++;

        stSvDiff.iSvStart = itr->iPos;
        stSvDiff.iSvEnd = itr->iEnd;
        vSvDiff.push_back(stSvDiff);
    }

    sort(vSvDiff.begin(), vSvDiff.end(), sort_svdiff);
    for(vector<St_SvDiff>::iterator itr = vSvDiff.begin(); itr != vSvDiff.end(); itr++)
    {
        cout << to_string(itr->GetDiff()) << "   ----   "
             <<"(" << to_string(itr->iSvStart) << ", " << to_string(itr->iSvEnd) << ")"
             << "   ----   " << to_string(itr->GetInsertSize())
             << "   ----   " << itr->GetStrand() << endl;
    }

    cout <<"-----" << endl;
    cout << "Total SV      : " << to_string(vSvDEL.size()) << endl;
    cout << "Missing Number: " << to_string(iMissNum) << endl;
}

//结论，我们确实发现了这些真正的sv的区域确实都对应了很多的SV, 但是问题在于我们没有办法定一个非常准确的阈值去进行分组 -->
void ClsDebug::CheckStdSvCount(vector<St_DiscordantReads>& vDiscdRreads, vector<St_SV>& vSvDEL)
{
    //Find the min and max offset
    int iMaxPos = 0;
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin(); itr != vDiscdRreads.end(); itr++)
    {
        if(iMaxPos < itr->iMatePos)
            iMaxPos = itr->iMatePos;
    }

    if(iMaxPos == 0)
    {
        cout << "False Length 0" << endl;
        return;
    }

    cout << endl << "Max Pos: " << to_string(iMaxPos) << endl;

    //-->Check hit status
    vector<int> vWholeSeq;
    size_t iSize = iMaxPos+2;
    vWholeSeq.resize(iSize, 0);

    //count times
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin(); itr != vDiscdRreads.end(); itr++)
    {
        if(itr->iReadsPos >= 0 && itr->iMatePos >= 0)
        {
            //cout << "(" << to_string(itr->iReadsPos) << ", " << to_string(itr->iMatePos) << ")" << endl;
            for_each(vWholeSeq.begin()+itr->iReadsPos, vWholeSeq.begin()+itr->iMatePos, AddVal(1));
        }
    }

    cout << "Counting finished" << endl;

//    int i=0;
//    for(vector<int>::iterator itr = vWholeSeq.begin(); itr != vWholeSeq.end(); itr++, i++)
//    {
//        if(i < 10000)
//        cout << to_string(*itr) << endl;
//    }

//    return;

    //Try to check the hitting status of the std SV
    vector<int> vCount;
    vector<int> vLen;
    int iIndex = 1;
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        int iMaxCount = -1;
        int iLen = 0;
        if(itr->iEnd >= vWholeSeq.size())
        {
            vCount.push_back(iMaxCount);
            vLen.push_back(iLen);
            cout << "End too big!!!" << endl;
            continue;
        }

        //If sv dropped in discordant range
        for(int i = itr->iPos; i < itr->iEnd; i++)
        {
            if(iMaxCount < vWholeSeq[i])
            {
                iMaxCount = vWholeSeq[i];
                iLen = 1;
                for(i = i+1; i<itr->iEnd; i++)
                {
                    if(iMaxCount == vWholeSeq[i])
                        iLen++;
                    else
                    {
                        i = i-1;
                        break;
                    }
                }
            }
        }

        vCount.push_back(iMaxCount);
        vLen.push_back(iLen);
        iIndex++;
    }

    int i = 0;
    //output
    cout << "vSvDEL Size: " << to_string(vSvDEL.size()) << endl;
    cout << "vCount Size: " << to_string(vCount.size()) << endl;
    cout << "vLen Size  : " << to_string(vLen.size()) << endl;

    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++, i++)
    {
        cout << "(" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ")"
             << "  ---  " << "Max Count: " << to_string(vCount[i])
             << "  ---  " << "Len: " << to_string(vLen[i]) << endl;
    }

    //Check average hit number
    int64_t iValidLen = 0;
    int64_t iCountSum = 0;
    for(vector<int>::iterator itr = vWholeSeq.begin(); itr != vWholeSeq.end(); itr++)
    {
        iCountSum += *itr;
        if(*itr != 0)
           iValidLen++;
    }
    float fAvg = (float)iCountSum/iValidLen;
    cout << endl         
         << "Average Hit Number: " << FloatToStr(fAvg) << endl;
    cout << "OrgValue          : " << to_string(fAvg)
         << " --- iCount Num        : " << to_string(iCountSum)
         << " --- iValid Len        : " << to_string(iValidLen)
         << endl;
}

void ClsDebug::GetDiscClipGroup(vector<St_BorderSCR>& vClipReads, vector<St_DiscordantReads>& vDiscdRreads)
{
    vector<St_DiscClipGroup> vDiscClipGroup;
    St_DiscClipGroup stDCGroup;
    //我们是将discordant reads 环绕着 clip reads来
    for(vector<St_BorderSCR>::iterator itr = vClipReads.begin(); itr != vClipReads.end(); itr++)
    {
        stDCGroup.Clear();
        stDCGroup.stClipReads = *itr;

        for(vector<St_DiscordantReads>::iterator subItr = vDiscdRreads.begin();
            subItr != vDiscdRreads.end(); subItr++)
        {
            if(itr->iClipPos >= subItr->iReadsPos && itr->iClipPos <= subItr->iMatePos)
            {
                //stDCGroup.vDiscdReads.push_back(&(*subItr));
                stDCGroup.iCount++;
            }
        }
        vDiscClipGroup.push_back(stDCGroup);
    }

    //Let's Count
    int iMinCount = 100;
    int iPotential = 0;
    cout << endl << "************ Hit Discordant Clip Group *************" << endl;
    int iContributeCount = 0;
    int iContributeNum = 0;
    int64_t iSumAvg = 0;
    for(vector<St_DiscClipGroup>::iterator itr = vDiscClipGroup.end() - 1; itr >= vDiscClipGroup.begin(); itr--)
    {
        if(itr->stClipReads.bSvContribute)
        {
            cout << to_string(itr->iCount) << endl;
            iContributeCount += itr->iCount;
            iContributeNum++;
        }
        //if(itr->iCount <= iMinCount)
        //    vDiscClipGroup.erase(itr);
        //else
        //    iPotential++;

        iSumAvg += itr->iCount;
    }

    cout << "Original Clip Reads Num      : " << to_string(vClipReads.size()) << endl;
    cout << "Original Discordant Reads Num: " << to_string(vDiscdRreads.size()) << endl;

    cout << endl
         << "iContributeCount: " << to_string(iContributeCount) << endl;
    cout << "iContributeNum  : " << to_string(iContributeNum) << endl;
    cout << "Avg             : " << to_string((float)iContributeCount / iContributeNum) << endl;

    cout << endl
         << "iSumAvg       : " << to_string(iSumAvg) << endl;
    cout << "vDiscClipGroup: " << to_string(vDiscClipGroup.size()) << endl;
    cout << "Avg           : " << to_string((float)iSumAvg / vDiscClipGroup.size()) << endl;
}

bool sort_border_scr(St_BorderSCR stV1, St_BorderSCR stV2)
{
    if(stV1.iDiscCountNum < stV2.iDiscCountNum)
        return true;
    else
        return false;
}

void ClsDebug::FilterSCReadsByDiscReads(vector<St_BorderSCR>& vBorderSCR,
                                        vector<St_DiscordantReads>& vDiscdRreads)
{
    //--------------->Get Count Number
    //Find the min and max offset
    int iMaxPos = 0;
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin(); itr != vDiscdRreads.end(); itr++)
    {
        if(iMaxPos < itr->iMatePos)
            iMaxPos = itr->iMatePos;
    }

    if(iMaxPos == 0)
    {
        cout << "False Length 0" << endl;
        return;
    }

    cout << endl << "Max Pos: " << to_string(iMaxPos) << endl;

    //-->Check hit status
    vector<int> vWholeSeq;
    size_t iSize = iMaxPos+2;
    vWholeSeq.resize(iSize, 0);

    //count times
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin(); itr != vDiscdRreads.end(); itr++)
    {
        if(itr->iReadsPos >= 0 && itr->iMatePos >= 0)
        {
            //cout << "(" << to_string(itr->iReadsPos) << ", " << to_string(itr->iMatePos) << ")" << endl;
            for_each(vWholeSeq.begin()+itr->iReadsPos, vWholeSeq.begin()+itr->iMatePos, AddVal(1));
        }
    }

    cout << "Counting finished" << endl;
    //<---------------------

    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr != vBorderSCR.end(); itr++)
    {
        if(itr->enClipPart = cpRight)
        {
            if(itr->iClipPos + OFFSET < iMaxPos)
                itr->iDiscCountNum = vWholeSeq[itr->iClipPos + OFFSET];
        }
        else if(itr->enClipPart = cpLeft)
        {
            if(itr->iClipPos - OFFSET < iMaxPos &&
               itr->iClipPos - OFFSET > 0)
                itr->iDiscCountNum = vWholeSeq[itr->iClipPos - OFFSET];
        }
    }

    sort(vBorderSCR.begin(), vBorderSCR.end(), sort_border_scr);

    int i = 0;
    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr != vBorderSCR.end(); itr++, i++)
    {
        if(i > 1000)
            break;

        cout << to_string(itr->iDiscCountNum) << endl;
    }
    cout << "Median value: " << to_string(vBorderSCR[vBorderSCR.size()/2].iDiscCountNum) << endl;
}

string ClsDebug::SaveClipReads(vector<St_BorderSCR>& vBorderSCR, string strFilePath)
{
    if(vBorderSCR.empty())
    {
        cout << "No reads detected" << endl;
        return "Nil";
    }

    if(strFilePath == "")
    {
        strFilePath = "../TempFile/clippedReads.fq";
    }

    ofstream ofsFq;
    ofsFq.open(strFilePath.c_str());
    int iReadsLen = vBorderSCR.begin()->strQuerySeq.length();
    string strQuality(iReadsLen, '#');
    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr != vBorderSCR.end(); itr++)
    {
        ofsFq << "@" << itr->strName << endl;
        ofsFq << itr->strQuerySeq << endl;
        ofsFq << "+" << endl;
        ofsFq << strQuality << endl;
    }
    ofsFq.close();
    return strFilePath;
}

void ClsDebug::LocalAssemblyByVelvet(string strReadsPath)
{
    //step 1: Use velveth to create hash table --> velveth
    //step 2: Use velvetg for local assembly --> velvetg
    ClsVelvet::GetInstance().LocalAssembly(strReadsPath, 100, 11, false);

    //step 3: Save and output the result
}

void ClsDebug::ClusterClipReads(vector<St_BorderSCR>& vBorderSCR, vector<St_SV>& vSvDEL)
{
    const int KMERLEN = 9;
    const int CLUSTERRANGE = 30; // 50 bps
    const int SETDIFF = 30;
    //Create K-mer table by unordered map
    unordered_map<unsigned, vector<St_SCRSet> > mpKmerTable;

    int iIter = 0;
    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr != vBorderSCR.end(); itr++)
    {                
        iIter++;

        if(itr->bSvContribute)
            cout << IntToStr(iIter) << endl;
        else
            continue;

        //use original sequence of reads to create k-mer table
        for(int i=0; i<itr->strQuerySeq.length() - KMERLEN + 1; i++)
        {
            string strKmer = itr->strQuerySeq.substr(i, KMERLEN);
            unsigned int uiKmer = ConvertKmerToNum32(strKmer);
            bool bFindKey = false;
            for(unordered_map<unsigned int, vector<St_SCRSet> >::iterator itrMpKmer = mpKmerTable.begin();
                itrMpKmer != mpKmerTable.end(); itrMpKmer++)
            {
                if(itrMpKmer->first == uiKmer)
                {
                    //Check if could join the existed cluster
                    bool bClustExist = false;
                    for(vector<St_SCRSet>::iterator itrSCRSet = itrMpKmer->second.begin();
                        itrSCRSet != itrMpKmer->second.end(); itrSCRSet++)
                    {
                        for(vector<St_BorderSCR>::iterator itrSCR = itrSCRSet->vSet.begin(); itrSCR != itrSCRSet->vSet.end(); itrSCR++)
                        {
                            //if(itr->iClipPos == itr->iReadsMapPos)
                            //{
                                if(abs(itrSCR->iReadsMapPos - itr->iReadsMapPos) < CLUSTERRANGE) //Assume it is the ending part of SV
                                {
                                    itrSCRSet->vSet.push_back(*itr);
                                    bClustExist = true;
                                    break;
                                }
                            //}
                            else
                            {
                                if(abs(itrSCR->iClipPos - itr->iClipPos) < CLUSTERRANGE) // Assume it is the beginning part of SV
                                {
                                    itrSCRSet->vSet.push_back(*itr);
                                    bClustExist = true;
                                    break;
                                }
                            }
                        }
                        if(bClustExist)
                            break;
                    }
                    if(!bClustExist)
                    {
                        St_SCRSet stSCRSet;
                        stSCRSet.vSet.push_back(*itr);
                        itrMpKmer->second.push_back(stSCRSet);
                    }
                    //Find k-mer
                    bFindKey = true;
                }
                if(bFindKey)
                    break;
            }

            //Failed to find this k-mer
            if(!bFindKey)
            {
                vector<St_SCRSet> vSCRSet;
                St_SCRSet stSCRSet;
                stSCRSet.vSet.push_back(*itr);
                vSCRSet.push_back(stSCRSet);
                mpKmerTable.insert(std::pair<unsigned int, vector<St_SCRSet> >(uiKmer, vSCRSet));
                //mpKmerTable[uiKmer] = vSCRSet;
                //cout << IntToStr(mpKmerTable.size()) << endl;
            }
        }
    }

    //Already get all elements of clustering --> Now get the real cluster
    cout << endl << "Length of K-mer Table: " << IntToStr(mpKmerTable.size()) << endl;

    //Calculate average Clip Pos for each different SCR Set -->
    for(unordered_map<unsigned int, vector<St_SCRSet> >::iterator itrMpKmer = mpKmerTable.begin();
        itrMpKmer != mpKmerTable.end(); itrMpKmer++)
    {
        if(itrMpKmer->second.size() > 1)
        {
            cout << "K-mer: " << ::ConvertNum32ToKmer(itrMpKmer->first, KMERLEN) << endl;
        }

        for(vector<St_SCRSet>::iterator itrSCRSet = itrMpKmer->second.begin();
            itrSCRSet != itrMpKmer->second.end(); itrSCRSet++)
        {
            //-->Get Average position for each SCR Set
            bool bSamePosMapClip = false;
            int iSamePosMapClipCount = 0;
            for(vector<St_BorderSCR>::iterator itrSCR = itrSCRSet->vSet.begin();
                itrSCR != itrSCRSet->vSet.end(); itrSCR++)
            {
                itrSCRSet->iAvgClipPos += itrSCR->iClipPos;
                itrSCRSet->iAvgMapPos += itrSCR->iReadsMapPos;

                if(itrSCR->iReadsMapPos == itrSCR->iClipPos)
                {
                    bSamePosMapClip = true;
                    iSamePosMapClipCount++;
                }

                if(itrMpKmer->second.size() > 1)
                {
                    cout << "\t" << "<" << to_string(itrSCR->iReadsMapPos) << ", " << to_string(itrSCR->iClipPos) << "> ";
                }
            }            
            itrSCRSet->iAvgClipPos /= itrSCRSet->vSet.size();
            itrSCRSet->iAvgMapPos /= itrSCRSet->vSet.size();            

            if(bSamePosMapClip)
            {
                if(itrSCRSet->vSet.size() == iSamePosMapClipCount ||
                   (itrSCRSet->vSet.size() > iSamePosMapClipCount && iSamePosMapClipCount > 1))
                {
                    itrSCRSet->bLeftYes = false;
                    itrSCRSet->bRightYes = true;
                }
                else
                {
                    itrSCRSet->bLeftYes = true;
                    itrSCRSet->bRightYes = false;
                }
            }
            else
            {
                itrSCRSet->bLeftYes = true;
                itrSCRSet->bRightYes = false;
            }

            if(itrMpKmer->second.size() > 1)
            {
                //-->Cout Avg mapping position and clip position
                cout << " -- " << "Avg_Map: " << to_string(itrSCRSet->iAvgMapPos) << ", " << "Avg_Clip: " << to_string(itrSCRSet->iAvgClipPos);
                cout << endl;
            }
            //<--
        }
    }
    cout << "Done: Calculate average Clip Pos for each different SCR Set" << endl;
    //<--

    //Merge clusters affiliated with different k-mer in k-mer table
    vector<St_ClusterPair> vClusterPair;
    vector<St_ClusterPair> vTmpPairSet;
    for(unordered_map<unsigned int, vector<St_SCRSet> >::iterator itrMpKmer = mpKmerTable.begin();
        itrMpKmer != mpKmerTable.end(); itrMpKmer++)
    {
        if(itrMpKmer->second.size() <= 1)
            continue;

        //Get all of pairs --> For this round
        vTmpPairSet.clear();

        //Get the longest Set -->
        vector<St_SCRSet>::iterator itrAnchorSet = itrMpKmer->second.begin();
        for(vector<St_SCRSet>::iterator itrSCRSet = itrMpKmer->second.begin() + 1;
            itrSCRSet < itrMpKmer->second.end(); itrSCRSet++)
        {
            if(itrAnchorSet->vSet.size() < itrSCRSet->vSet.size())
                itrAnchorSet = itrSCRSet;
        }

        //<--
        for(vector<St_SCRSet>::iterator itrSCRSetPair = itrMpKmer->second.begin();
            itrSCRSetPair != itrMpKmer->second.end(); itrSCRSetPair++)
        {
            if(itrAnchorSet == itrSCRSetPair)
                continue;

            St_ClusterPair stCP;
            if(itrAnchorSet->iAvgClipPos < itrSCRSetPair->iAvgMapPos)
            {
               if(itrAnchorSet->bLeftYes && itrSCRSetPair->bRightYes)
               {
                   stCP.stLeftSCRSet = *itrAnchorSet;
                   stCP.stRightSCRSet = *itrSCRSetPair;
                   stCP.iCount = 1;
                   vTmpPairSet.push_back(stCP);
               }
            }
            else if(itrAnchorSet->iAvgMapPos > itrSCRSetPair->iAvgClipPos)
            {
               if(itrSCRSetPair->bLeftYes && itrAnchorSet->bRightYes)
               {
                   stCP.stLeftSCRSet = *itrSCRSetPair;
                   stCP.stRightSCRSet = *itrAnchorSet;
                   stCP.iCount = 1;
                   vTmpPairSet.push_back(stCP);
               }
            }
        }


//        for(vector<St_SCRSet>::iterator itrSCRSet = itrMpKmer->second.begin();
//            itrSCRSet != itrMpKmer->second.end(); itrSCRSet++)
//        {
////            if(itrSCRSet->vSet.size() <= 1)
////                continue;

//            for(vector<St_SCRSet>::iterator itrSCRSetPair = itrSCRSet + 1;
//                itrSCRSetPair != itrMpKmer->second.end(); itrSCRSetPair++)
//            {
////                if(itrSCRSetPair->vSet.size() <= 1)
////                    continue;

//                //Create Pair -->
//                St_ClusterPair stCP;
//                if(itrSCRSet->iAvgMapPos < itrSCRSetPair->iAvgMapPos)
//                {
//                    //if(itrSCRSet->bLeftYes && itrSCRSetPair->bRightYes)
//                    {
//                        stCP.stLeftSCRSet = *itrSCRSet;
//                        stCP.stRightSCRSet = *itrSCRSetPair;
//                        stCP.iCount = 1;
//                        vTmpPairSet.push_back(stCP);
//                    }
//                }
//                else
//                {
//                    //if(itrSCRSetPair->bLeftYes && itrSCRSet->bRightYes)
//                    {
//                        stCP.stLeftSCRSet = *itrSCRSetPair;
//                        stCP.stRightSCRSet = *itrSCRSet;
//                        stCP.iCount = 1;
//                        vTmpPairSet.push_back(stCP);
//                    }
//                }
//            }
//        }
//        cout << "vTmpPairSet Size: " << IntToStr(vTmpPairSet.size()) << endl;

//        for(vector<St_ClusterPair>::iterator itrTmp = vTmpPairSet.begin(); itrTmp != vTmpPairSet.end(); itrTmp++)
//        {
//            cout << "\t" << " --- " << "Left Boundary: " << to_string(itrTmp->stLeftSCRSet.iAvgClipPos)
//                         << " --- " << "Right Boundary: " << to_string(itrTmp->stRightSCRSet.iAvgMapPos) << endl;
//        }

        //Merge all of temporary pair set with the main container
        for(vector<St_ClusterPair>::iterator itrTmp = vTmpPairSet.begin(); itrTmp != vTmpPairSet.end(); itrTmp++)
        {
            bool bFind = false;
            for(vector<St_ClusterPair>::iterator itr = vClusterPair.begin(); itr != vClusterPair.end(); itr++)
            {
                if( (abs( int(itr->stLeftSCRSet.iAvgClipPos - itrTmp->stLeftSCRSet.iAvgClipPos) ) < SETDIFF ||
                     abs( int(itr->stLeftSCRSet.iAvgMapPos - itrTmp->stLeftSCRSet.iAvgMapPos) ) < SETDIFF)
                        &&
                    abs( int(itr->stRightSCRSet.iAvgMapPos - itrTmp->stRightSCRSet.iAvgMapPos)) < SETDIFF )
                {
                    UpdateSCRSet(itr->stLeftSCRSet, itrTmp->stLeftSCRSet);
                    UpdateSCRSet(itr->stRightSCRSet, itrTmp->stRightSCRSet);
                    itr->iCount++;
                    bFind = true;
                    break;
                }
            }

            if(!bFind)
            {
                vClusterPair.push_back(*itrTmp);
            }
        }
    }

    //进一步进行过滤 -->　左边吧clip和map

    cout << "Before iCount Filter: " << IntToStr(vClusterPair.size()) << endl;

    if(vClusterPair.empty())
    {}
    else
    {
        for(vector<St_ClusterPair>::iterator itr = vClusterPair.end() - 1; itr >= vClusterPair.begin(); itr--)
        {
            if(itr->iCount < 6) // 15bps should be covered
                vClusterPair.erase(itr);
//            else
//                cout << IntToStr(itr->iCount) << endl;
        }
    }

    //Output -> The cluster result
    cout << endl << "------- vClusterPair ---------" << endl;

    //Cout each cluster --> Go----------
    for(vector<St_ClusterPair>::iterator itr = vClusterPair.begin(); itr != vClusterPair.end(); itr++)
    {
        cout << "Cluster Count: " << to_string(itr->iCount)
             << " --- " << "Left_Avg_Map_Clip_Pos: " << to_string(itr->stLeftSCRSet.iAvgMapPos) << ", " << to_string(itr->stLeftSCRSet.iAvgClipPos)
             << " --- " << "Right_Avg_Map_Clip_Pos: " << to_string(itr->stRightSCRSet.iAvgMapPos) << ", " << to_string(itr->stRightSCRSet.iAvgClipPos)
             << endl;
        //-->Get Average position for each SCR Set
        for(vector<St_BorderSCR>::iterator itrSCR = itr->stLeftSCRSet.vSet.begin();
            itrSCR != itr->stLeftSCRSet.vSet.end(); itrSCR++)
        {
            cout << "\t" << "<" << to_string(itrSCR->iReadsMapPos) << ", " << to_string(itrSCR->iClipPos) << "> ";
        }
        cout << endl;

        for(vector<St_BorderSCR>::iterator itrSCR = itr->stRightSCRSet.vSet.begin();
            itrSCR != itr->stRightSCRSet.vSet.end(); itrSCR++)
        {
            cout << "\t" << "<" << to_string(itrSCR->iReadsMapPos) << ", " << to_string(itrSCR->iClipPos) << "> ";
        }
        cout << endl;
//        //-->Cout Avg mapping position and clip position
//        cout << " \t-- " << "Avg_Map: " << to_string(itrSCRSet->iAvgMapPos) << ", " << "Avg_Clip: " << to_string(itrSCRSet->iAvgClipPos);
//        cout << endl;
    }
    //<--

    cout << "After iCount Filter : " << IntToStr(vClusterPair.size()) << endl;

    //Check how many get real SV
    cout << endl << "------------- Check how many get real SV -------------" << endl;
    int iHitRealSVNum = 0;
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        for(vector<St_ClusterPair>::iterator itrCP = vClusterPair.begin(); itrCP != vClusterPair.end(); itrCP++)
        {
            //For left breakpoint
            bool bLeftHit = false;
            if(abs( int(itr->iPos - itrCP->stLeftSCRSet.iAvgClipPos) ) < SETDIFF)
            {
                bLeftHit = true;
            }

            //For right breakpoint
            bool bRightHit = false;
            if(abs( int(itr->iEnd - itrCP->stRightSCRSet.iAvgMapPos) ) < SETDIFF)
            {
                bRightHit = true;
            }

            if(bLeftHit && bRightHit)
            {
                iHitRealSVNum++;
                //output SV
                cout << "<" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ">" << endl;
                break;
            }
        }
    }
    cout << endl << "SV Sum: " << IntToStr(vSvDEL.size()) << endl;
    cout << "Hit SV: " << IntToStr(iHitRealSVNum) << endl;
}

void ClsDebug::UpdateSCRSet(St_SCRSet& stSRCSetOrg, St_SCRSet& stSRCSetTmp)
{
    //Merge SCR
    for(vector<St_BorderSCR>::iterator itr = stSRCSetTmp.vSet.begin(); itr != stSRCSetTmp.vSet.end(); itr++)
    {
        bool bExist = false;
        for(vector<St_BorderSCR>::iterator subItr = stSRCSetOrg.vSet.begin(); subItr != stSRCSetOrg.vSet.end(); subItr++)
        {
            if(itr->iReadsMapPos == subItr->iReadsMapPos &&
               itr->iClipPos == subItr->iClipPos)
            {
                bExist = true;
                break;
            }
        }

        if(!bExist)
        {
            stSRCSetOrg.vSet.push_back(*itr);
        }
    }

    //Update Average Pos
    stSRCSetOrg.iAvgMapPos = 0;
    stSRCSetOrg.iAvgClipPos = 0;
    for(vector<St_BorderSCR>::iterator itr = stSRCSetOrg.vSet.begin(); itr != stSRCSetOrg.vSet.end(); itr++)
    {
        stSRCSetOrg.iAvgMapPos += itr->iReadsMapPos;
        stSRCSetOrg.iAvgClipPos += itr->iClipPos;
    }
    stSRCSetOrg.iAvgMapPos /= stSRCSetOrg.vSet.size();
    stSRCSetOrg.iAvgClipPos /= stSRCSetOrg.vSet.size();
}

