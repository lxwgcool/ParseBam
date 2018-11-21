#include "clsparsebam.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "../../../../ShareLibrary/clsmuscle.h"
#include <algorithm>

#include "../../../../ShareLibrary/bamtools/include/api/BamReader.h"
using namespace BamTools;

string arryClipPart[cpMax] = {"Left", "Right", "Any"};
string arryBoundaySV[cpMax] = {"Right", "Left", "Any"};
const int MEDIAN_INSERT_SIZE = 386;
const int MEAN_INSERT_SIZE = 379;
const int STD_DEVIATION = 38;
const int MIN_DEL_SIZE = 40;

ClsParseBam::ClsParseBam()
{
    string strTmp = "Hello World";
    int i = atoi(strTmp.c_str());
    cout << IntToStr(i) << endl;
}

 ClsParseBam::~ClsParseBam()
{}

bool sort_scrClipPos_func(St_BorderSCR stScr1, St_BorderSCR stScr2)
{
    if(stScr1.iClipPos <= stScr2.iClipPos)
        return true;
    else
        return false;
}

void ClsParseBam::ReadBamFile(string strBamFilePath, vector<St_BorderSCR>& vBorderSCR)
{
    ofstream ofsClippedReads;
    ofsClippedReads.open("./ClippedReads.txt");

    //Go to check sorted bam file
    //1: parse bam file & output the expected reads
    BamReader* pBamReader = new BamReader();

    int iTotalNum = 0;
    int iMappedNum = 0;
    int iClippedSum = 0;
    int iSoftClipNum = 0;
    int iHardClipNum = 0;

    vector<St_ClipReads> vClipReads;
    St_ClipReads stClipReads;

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
//        //Cout << the specific reads
//        if(al.QueryBases == "TTAAAAAAAAGTCTGCACTAATGCCACAATATTGATAGGAATAAGCCATAATTACTCCTGAAAATCCTTCAAATGACTTTTGTACATACTTTACTTGTTGT" ||
//           al.QueryBases == "CACACTACTTCAAGTAGACAATATCTTTTAAAAAAAAGTCTGCACTAATGCCACAATATTGATAGGAATAAGCCATAATTACTCCTGAAAATCCTTCAAAT" ||
//           al.QueryBases == "TAATGCCACAATATTGATAGGAATAAGCCATAATTACTCCTGAAAATCCTTCAAATGACTTTTGTACATACTTTACTTGTTGTTCGCTAAACGTGCCTTGT" ||
//           al.QueryBases == "CCTGAAAATCCTTCAAATGACTTTTGTACATACTTTACTTGTTGTTCACTAAACGTGCCTTGTATATAGTCTACCTTTGCACTTTTTCTCTCTCTGATTCA" ||
//           al.QueryBases == "GAAGAAAGAGAGAGGGAGGGAGGGAAGGAAGGAGGGAAGGAAGGAAAGGGGAGAGAGAGAGAAAAGGAGATGGGAGGGAATTCATAAGGCATAAATGAAAA")
//        {
//            cout << "\t Is_Reverse_Strand: " << string(al.IsReverseStrand() ? "Yes" : "No") << endl;
//            cout << "\t Is_Mate_Reverse_Strand: " << string(al.IsMateReverseStrand() ? "Yes" : "No") << endl;
//        }
//        else
//            continue;

        iTotalNum++;
        //Case 1: do not map  ------->
        if(!al.IsMapped() || al.QueryBases == "")
            continue;

        if(al.QueryBases.find('N') != string::npos ||
           al.QueryBases.find('n') != string::npos)
            continue;

        int iRealInsertSize = al.InsertSize;
        if(iRealInsertSize == 0 && al.IsMateMapped())
        {
            iRealInsertSize = al.MatePosition - al.Position;
            if(iRealInsertSize > 0)
                iRealInsertSize += al.QueryBases.length();
            else
                iRealInsertSize -= al.QueryBases.length();
        }

        stClipReads.Clear();
        stClipReads.strName = al.Name;
        stClipReads.strQuerySeq = al.QueryBases;
        stClipReads.iReadsMapPos = al.Position;
        stClipReads.iInsertSize = iRealInsertSize;
        stClipReads.iMateMapPos = al.MatePosition;
        stClipReads.bFirstMate = al.IsFirstMate();
        stClipReads.bSecondMate = al.IsSecondMate();
        stClipReads.bMateMapped = al.IsMateMapped();

        //Case 2: If the reads maps fine
        ///1: Check if there is soft clip
        bool bFindClip = false;
        int iOffSet = 0;

//        //Print the status of each Cigar tag
//        for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
//            itr != al.CigarData.end(); itr++)
//        {
//            cout << itr->Type << " -- " << IntToStr(itr->Length) << endl;
//        }

        for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
            itr != al.CigarData.end(); itr++)
        {
            switch(itr->Type)
            {
                case 'M': // alignment match (can be a sequence match or mismatch)
                    iMappedNum++;
                    iOffSet += itr->Length;
                    break;
                case 'I': // insertion to the reference
                    break;
                case 'D': // deletion from the reference
                case 'N':  // skipped region from the reference
                    break;
                case 'S':  // soft clipping (clipped sequences present in SEQ)
                    iSoftClipNum++;
                    iClippedSum++;
//                    cout << IntToStr(al.QueryBases.length()) << " -- " << to_string(iOffSet) << " : " << to_string(itr->Length)
//                         << endl;
                    stClipReads.vClipSeq.push_back(al.QueryBases.substr(iOffSet, itr->Length));
                    stClipReads.vPos.push_back(iOffSet);                    
                    iOffSet += itr->Length;
                    bFindClip = true;
                    break;
                case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                    iHardClipNum++;
                    iClippedSum++;
                    break;
                case 'P': // padding (silent deletion from padded reference)
                case '=': // sequence match
                case 'X': // sequence mismatch
                    break;
            }
        }

        //cout << "-----Ciga Finished: Align VS Query: "
        //     << to_string(al.AlignedBases.length()) << " -- " << to_string(al.QueryBases.length()) << endl;

        if(bFindClip)
        {
//            //Record Two types of sequence  -->
//            string strTmp = "> Query:Align " + to_string(al.QueryBases.length()) + ":" +
//                            to_string(al.AlignedBases.length()) +
//                            " -- clipNum:" + to_string(stClipReads.vClipSeq.size());
//            ofsClippedReads << strTmp << endl;

//            //1: The original sequence
//            ofsClippedReads << "Query" << endl << al.QueryBases << endl;

//            //2: The cliped seqeunce
//            if(!stClipReads.vClipSeq.empty())
//            {
//                ofsClippedReads << "----Clip----" << endl;
//                for(vector<string>::iterator itr = stClipReads.vClipSeq.begin();
//                    itr != stClipReads.vClipSeq.end(); itr++)
//                {
//                    ofsClippedReads << *itr << endl;
//                }
//                ofsClippedReads << "------------" << endl;
//            }

//            //3: Output alignbased
//            ofsClippedReads << "Aligned" << endl << al.AlignedBases << endl << endl;

            //Save clipped info
            vClipReads.push_back(stClipReads);
        }
    }

    //Summarize
    ofsClippedReads << "Total Reads Num : " << to_string(iTotalNum) << endl;
    ofsClippedReads << "Mapped Reads Num: " << to_string(iMappedNum) << endl;
    ofsClippedReads << "     Clipped Sum: " << to_string(iClippedSum) << endl;
    ofsClippedReads << "       Soft Sum : " << to_string(iSoftClipNum) << endl;
    ofsClippedReads << "       Hard Sum : " << to_string(iHardClipNum) << endl;
    ofsClippedReads << "ClippedReads Sum: " << IntToStr(vClipReads.size()) << endl;

    ofsClippedReads.close();

    //---->Initial St_BorderSCR
    vBorderSCR.clear();
    St_BorderSCR stBorderSCR;
    //<----

    //Output the target clipe reads --->
    ofstream ofsTargetInfo;
    ofsTargetInfo.open("./TargetInfo.txt");
    ofstream ofsTargetReads;
    ofsTargetReads.open("./TargetReads.txt");

    En_ClipPart enClipPart = cpMax;
    int iSum = 0;
    for(vector<St_ClipReads>::iterator itr = vClipReads.begin(); itr != vClipReads.end(); itr++)
    {
        //Pick the clipped reads in the boundary
        bool bBorder = false;

        int iIndex = 0;
        string strClipSeq = "";
        int iClipPosInRef = -1;
        for(vector<int>::iterator itrSCPos = itr->vPos.begin(); itrSCPos != itr->vPos.end();
            itrSCPos++, iIndex++)
        {
            //Check if it is the last clip part from the right side
            int iRightOffSet = *itrSCPos - (itr->strQuerySeq.length() - itr->vClipSeq[iIndex].length());
            if(abs(*itrSCPos) <= 10)  //Left side
            {
                enClipPart = cpLeft;
                strClipSeq = itr->vClipSeq[iIndex];
                iClipPosInRef = itr->iReadsMapPos + *itrSCPos;
                bBorder = true;
                break;
            }
            else //if(abs(iRightOffSet) <= 10)
            {
                enClipPart = cpRight;
                strClipSeq = itr->vClipSeq[iIndex];
                iClipPosInRef = itr->iReadsMapPos + itr->strQuerySeq.length() -
                                itr->vClipSeq[iIndex].length();
                bBorder = true;
                break;
            }
//            else
//                continue;
        }

        /* Format of "TargetInfo"
         * >1
         * ReadsClipPart         : Left
         * Support Boundary of SV: Right
         * ReadsAlignPos in Ref  : 3234
         * ClipPos in Reference  : 3234
         * ClipLength            : 12
         * ClipSeq               : ATGCAAAATTTTCCCCC
         *
         *
         * Format of "TargetReads" (Name: index ReadsMap_pos Clip_Pos)
         * >1 3234 3234
         * AAAAAAAATTTTTTCCCCGGGGGGGGGGGGGGGGGAAAT
         */
        if(bBorder)
        {
            if( itr->strQuerySeq.find('N') != string::npos ||
                itr->strQuerySeq.find('n') != string::npos)
                continue;

            //For "TargetInfo"
            ofsTargetInfo << ">" << IntToStr(iSum) << endl;            
            ofsTargetInfo << "Name                  : " << itr->strName << endl;
            ofsTargetInfo << "ReadsClipPart         : " << arryClipPart[enClipPart] << endl;
            ofsTargetInfo << "Support Boundary of SV: " << arryBoundaySV[enClipPart] << endl;
            ofsTargetInfo << "ReadsAlignPos in Ref  : " << IntToStr(itr->iReadsMapPos) << endl;
            ofsTargetInfo << "ClipPos in Reference  : " << IntToStr(iClipPosInRef) << endl;
            ofsTargetInfo << "ClipLength            : " << IntToStr(itr->vClipSeq[iIndex].length()) << endl;
            ofsTargetInfo << "ClipSeq               : " << itr->vClipSeq[iIndex] << endl;

            //For "Target Sequence"
            ofsTargetReads << "> " << IntToStr(iSum) << " "
                           << IntToStr(itr->iReadsMapPos) << " "
                           << IntToStr(iClipPosInRef) << endl;
            ofsTargetReads << itr->strQuerySeq << endl;
            //Count
            iSum++;

            //Save Border reads
            stBorderSCR.enClipPart = enClipPart;
            stBorderSCR.iClipPos = iClipPosInRef;
            stBorderSCR.iReadsMapPos = itr->iReadsMapPos;
            stBorderSCR.strClipSeq = itr->vClipSeq[iIndex];
            stBorderSCR.strQuerySeq = itr->strQuerySeq;
            stBorderSCR.iInsertSize = itr->iInsertSize;
            stBorderSCR.strName = itr->strName;
            stBorderSCR.iMatPos = itr->iMateMapPos;
            stBorderSCR.bFirstMate = itr->bFirstMate;
            stBorderSCR.bSecondMate = itr->bSecondMate;
            stBorderSCR.bMateMapped = itr->bMateMapped;
            vBorderSCR.push_back(stBorderSCR);
            stBorderSCR.Clear();
        }
    }
    ofsTargetInfo.close();
    ofsTargetReads.close();
    //<---

    //sort soft clip reads by clip position  from small to large
    sort(vBorderSCR.begin(), vBorderSCR.end(), sort_scrClipPos_func);

    delete pBamReader;
    pBamReader = NULL;    
}

const int OFFSETSV =10;  //90% of org reads
void ClsParseBam::GetDELRelatedReads(vector<St_SV>& vSvDEL,
                                     vector<St_BorderSCR>& vBorderSCR,
                                     vector<St_SvDelReads>& vSvDelReads)
{
    //iterator SV deletion
    vSvDelReads.clear();
    St_SvDelReads stSvDelReads;

    for(vector<St_SV>::iterator itrSvDel = vSvDEL.begin(); itrSvDel != vSvDEL.end(); itrSvDel++)
    {
        //iterator border soft-clipped reads
        for(vector<St_BorderSCR>::iterator itrBorderScR = vBorderSCR.begin();
            itrBorderScR != vBorderSCR.end(); itrBorderScR++)
        {
            //For SV Left reads
            if(abs(itrBorderScR->iClipPos - itrSvDel->iPos) <= OFFSETSV)
            {
                //this reads is target, and we need to insert it into the container (vSvDelReads)
                bool bFind = false;
                for(vector<St_SvDelReads>::iterator itrSvDelReads = vSvDelReads.begin();
                    itrSvDelReads != vSvDelReads.end(); itrSvDelReads++)
                {
                    if(bFind)
                        break;
                    if(itrSvDelReads->stSV == *itrSvDel)
                    {
                        bFind = true;
                        //Inset Left clip vector or right clip vector
//                        switch(itrBorderScR->enClipPart)
//                        {
//                            case cpLeft:
//                                itrSvDelReads->vLeftClipReads.push_back(*itrBorderScR);
//                                break;
//                            case cpRight:
//                                itrSvDelReads->vRightClipReads.push_back(*itrBorderScR);
//                                break;
//                            default:
//                                break;
//                        }
                        itrSvDelReads->vSVLeftReads.push_back(*itrBorderScR);
                    }
                }
                if(!bFind) //New sv supported reads
                {
                    stSvDelReads.Clear();
                    stSvDelReads.stSV = *itrSvDel;
//                    switch(itrBorderScR->enClipPart)
//                    {
//                        case cpLeft:
//                            stSvDelReads.vLeftClipReads.push_back(*itrBorderScR);
//                            break;
//                        case cpRight:
//                            stSvDelReads.vRightClipReads.push_back(*itrBorderScR);
//                            break;
//                        default:
//                            break;
//                    }
                    stSvDelReads.vSVLeftReads.push_back(*itrBorderScR);
                    vSvDelReads.push_back(stSvDelReads);
                }
            }
            //For SV Right reads
            else if(abs(itrBorderScR->iClipPos - itrSvDel->iEnd) <= OFFSETSV)
            {
                //this reads is target, and we need to insert it into the container (vSvDelReads)
                bool bFind = false;
                for(vector<St_SvDelReads>::iterator itrSvDelReads = vSvDelReads.begin();
                    itrSvDelReads != vSvDelReads.end(); itrSvDelReads++)
                {
                    if(bFind)
                        break;
                    if(itrSvDelReads->stSV == *itrSvDel)
                    {
                        bFind = true;
                        itrSvDelReads->vSVRightReads.push_back(*itrBorderScR);
                    }
                }
                if(!bFind) //New sv supported reads
                {
                    stSvDelReads.Clear();
                    stSvDelReads.stSV = *itrSvDel;
                    stSvDelReads.vSVRightReads.push_back(*itrBorderScR);
                    vSvDelReads.push_back(stSvDelReads);
                }
            }
        }
    }
}

void ClsParseBam::OutputSvDelReads(vector<St_SvDelReads>& vSvDelReads)
{
    /* The output format should be like this:
     * Notice: "1" is the reads index.
     * ----------
     * >Info
     * "reads Index" "Chrom Index" "SV Type" "SoftClip Position in Vcf" "SoftClip Position in Bam" "Alignment Position in Ref"
     * >QuerySeq
     * query sequence (should almost equal to reads sequence)
     * >ClipSeq
     * clip sequence
     * ----------
     * For example
     * ----------
     * >Info
     * 12 1 DEL 23453 23459 23450
     * >QuerySeq
     * AATGCGCGCGCGCGCGCGCCCCTTTAAAAAA
     * >ClipSeq
     * GCGCCCCTTTAAAAAA
     * ----------
     */

    /* We also need to output the statistic result
     * ----------
     * >SV: SV_Index
     * Type                      : Type string
     * Chrom Index               : Chrom_Index
     * Position                  : Alignment_Pos
     * End Pos                   : Ending Position
     * Length                    : Length
     * Support Reads - Left Clip : Number
     * Support Reads - Right Clip: Number
     * Sequence                  : sequence
     * ----------
     * For example
     * ----------
     * >SV: 1
     * Type                      : DEL
     * Chrom Index               : 11
     * Position                  : 34523
     * End Pos                   : 34600
     * Length                    : 40
     * Support Reads - Left Clip : 12
     * Support Reads - Right Clip: 56
     * Sequence                  : ATGCCCCTTAAA
     *
     * ----------
     */

    //Use Muscle to align the cliped seq (belongs to one SV break point)
    ClsMuscle* pMuscle = new ClsMuscle();

    ofstream ofsSvDelReads;
    ofsSvDelReads.open("./SvDelReads.txt");

    ofstream ofsSvStatistic;
    ofsSvStatistic.open("./SvStatistic.txt");

    ///output the cliped reads -->
    ofstream ofsClipTemp;
    string strClipTempFile = "./ClipTemp.txt";

    ofstream ofsSvMergedClipSeq;
    ofsSvMergedClipSeq.open("./SvMergedClipSeq.txt"); //record the multiple alignment result

    int iSum = 0;
    int iSvIndex = 0;
    for(vector<St_SvDelReads>::iterator itr = vSvDelReads.begin(); itr != vSvDelReads.end(); itr++)
    {
        //Output the sv Statistic Info -->
        ofsSvStatistic << ">SV: " << to_string(iSvIndex) << endl;
        ofsSvStatistic << "Type                      : " << itr->stSV.strType << endl;
        ofsSvStatistic << "Chrom Index               : " << to_string(itr->stSV.iChrom) << endl;
        ofsSvStatistic << "Position                  : " << to_string(itr->stSV.iPos) << endl;
        ofsSvStatistic << "End Pos                   : " << to_string(itr->stSV.iEnd) << endl;
        ofsSvStatistic << "Length                    : " << to_string(itr->stSV.iLen) << endl;
        ofsSvStatistic << "Support Reads - SV Left   : " << to_string(itr->vSVLeftReads.size()) << endl;
        ofsSvStatistic << "Support Reads - SV Right  : " << to_string(itr->vSVRightReads.size()) << endl;
        ofsSvStatistic << "Sequence                  : " << itr->stSV.strRef << endl;
        iSvIndex++;
        //<--

        int iTempIndex = 0;
        bool bFindValidAlignedClip = false;
        string strLastClipSeq = "";
        ofsClipTemp.open(strClipTempFile.c_str()); //prepare for multiple alignment
        for(vector<St_BorderSCR>::iterator subItr = itr->vSVLeftReads.begin();
            subItr != itr->vSVLeftReads.end(); subItr++)
        {
            ofsSvDelReads << ">Info" << endl;
            ofsSvDelReads << to_string(iSum)
                          << " " << to_string(itr->stSV.iChrom)
                          << " " << itr->stSV.strType
                          << " " << to_string(itr->stSV.iPos)
                          << " " << to_string(subItr->iClipPos)
                          << " " << to_string(subItr->iReadsMapPos) << endl;
            ofsSvDelReads << ">QuerySeq" << endl;
            ofsSvDelReads << subItr->strQuerySeq << endl;
            ofsSvDelReads << ">ClipSeq" << endl;
            ofsSvDelReads << subItr->strClipSeq << endl;

            if(subItr->enClipPart == cpRight && subItr->strClipSeq.length() > 10)
            {
                iTempIndex++;
                ofsClipTemp << ">" << to_string(iTempIndex) << endl;
                ofsClipTemp << subItr->strClipSeq << endl;
                strLastClipSeq = subItr->strClipSeq;
                if(!bFindValidAlignedClip)
                    bFindValidAlignedClip = true;
            }            
            iSum++;
        }
        ofsClipTemp.close();
        //-->Output Cliped reads to prepare the calculation later
        /* Some requirement:
         * For SV Left Reads -->
         * (1) only keep the right part clipped reads
         * (2) only keep the clipped part >= 10 bps
         * (3) multiple align those clipped part
         * (4) Merge those multiple clipped part to the longer string
         */
        if(bFindValidAlignedClip)
        {            
            string strName = itr->stSV.strType + "_" +
                             to_string(itr->stSV.iPos) + "_" +
                             to_string(itr->stSV.iEnd) + "_" +
                             "SV_Left";
            if(iTempIndex > 1) // More than one --> do multiple alignment
            {
                string strAlignedSeqPath = pMuscle->Run(strClipTempFile);
                string strMergedSeq = pMuscle->MergeAlignedSeq(strAlignedSeqPath);
                //save the aligned merged reads
                // a) Save name
                ofsSvMergedClipSeq << ">" << strName << endl;
                // b) save sequence
                ofsSvMergedClipSeq << strMergedSeq << endl;
            }
            else //For the case only supported by less than 1 number of reads
            {
                if(strLastClipSeq != "")
                {
                    // a) Save name
                    ofsSvMergedClipSeq << ">" << strName << endl;
                    // b) save sequence
                    ofsSvMergedClipSeq << strLastClipSeq << endl;
                }
            }
        }
        //<--

        iTempIndex = 0;
        bFindValidAlignedClip = false;
        strLastClipSeq = "";
        ofsClipTemp.open(strClipTempFile.c_str()); //prepare for multiple alignment
        for(vector<St_BorderSCR>::iterator subItr = itr->vSVRightReads.begin();
            subItr != itr->vSVRightReads.end(); subItr++)
        {
            ofsSvDelReads << ">Info" << endl;
            ofsSvDelReads << to_string(iSum)
                          << " " << to_string(itr->stSV.iChrom)
                          << " " << itr->stSV.strType
                          << " " << to_string(itr->stSV.iPos)
                          << " " << to_string(subItr->iClipPos)
                          << " " << to_string(subItr->iReadsMapPos) << endl;
            ofsSvDelReads << ">QuerySeq" << endl;
            ofsSvDelReads << subItr->strQuerySeq << endl;
            ofsSvDelReads << ">ClipSeq" << endl;
            ofsSvDelReads << subItr->strClipSeq << endl;

            if(subItr->enClipPart == cpLeft && subItr->strClipSeq.length() > 10)
            {
                iTempIndex++;
                ofsClipTemp << ">" << to_string(iTempIndex) << endl;
                ofsClipTemp << subItr->strClipSeq << endl;
                strLastClipSeq = subItr->strClipSeq;
                if(!bFindValidAlignedClip)
                    bFindValidAlignedClip = true;
            }

            iSum++;
        }
        ofsClipTemp.close();
        //-->Output Cliped reads to prepare the calculation later
        /* Some requirement:
         * For SV Right Reads -->
         * (1) only keep the Left part clipped reads
         * (2) only keep the clipped part >= 10 bps
         * (3) multiple align those clipped part
         * (4) Merge those multiple clipped part to the longer string
         */
        if(bFindValidAlignedClip)
        {
            string strName = itr->stSV.strType + "_" +
                             to_string(itr->stSV.iPos) + "_" +
                             to_string(itr->stSV.iEnd) + "_" +
                             "SV_Right";
            if(iTempIndex > 1)
            {
                string strAlignedSeqPath = pMuscle->Run(strClipTempFile);
                string strMergedSeq = pMuscle->MergeAlignedSeq(strAlignedSeqPath);
                //save the aligned merged reads

                // a) Save name
                ofsSvMergedClipSeq << ">" << strName << endl;
                // b) save sequence
                ofsSvMergedClipSeq << strMergedSeq << endl;
            }
            else
            {
                if(strLastClipSeq != "")
                {
                    // a) Save name
                    ofsSvMergedClipSeq << ">" << strName << endl;
                    // b) save sequence
                    ofsSvMergedClipSeq << strLastClipSeq << endl;
                }
            }
        }
        //<--
    }
    ofsSvDelReads.close();
    ofsSvStatistic.close();
    ofsSvMergedClipSeq.close();
    delete pMuscle;
    pMuscle = NULL;
}

struct St_MergeSeqAlignStatus // The alignment status of Multiple alignment merged sequence
{
    int iSvLeftAlignPos;
    int iSvRightAlignPos;
    string strAlignSeqSvLeft;
    string strAlignSeqSvRight;

    St_MergeSeqAlignStatus()
    {
        iSvLeftAlignPos = 0;
        iSvRightAlignPos = 0;
        strAlignSeqSvLeft = "";
        strAlignSeqSvRight = "";
    }

    void Clean()
    {
        iSvLeftAlignPos = 0;
        iSvRightAlignPos = 0;
        strAlignSeqSvLeft = "";
        strAlignSeqSvRight = "";
    }
};

//Parse the bam file to check if the related reads hits with the other part of SV in reference
void ClsParseBam::ReadBamFileMergedAlignReads(string strBamFilePath,
                                 vector<St_SvDelReads>& vSvDelReads, vector<St_SV>& vSvDEL)
{
    //Parse Bam file -->
    /* Notice;
     * (1) each reads matches with a specific SV (Deletion)
     * (2) Check if those merged cliped reads really support the correspondant boundary of SV
     */
    map<string, St_MergeSeqAlignStatus> mpAlignStatus;
    St_MergeSeqAlignStatus stAlignStatus;

    //1: Read Bam file
    BamReader* pBamReader = new BamReader();

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    //Save to a map<name_string, Left&&Right_Alignment_Status>  -->Go
    // *** May need to create a new data structure
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(!al.IsMapped()) // Only consider the mapped reads
            continue;

        int iBP = al.Name.rfind('_');
        string strSvName = al.Name.substr(0, iBP);
        string strPosPart = al.Name.substr(iBP + 1, al.Name.length()-iBP-1);
        int iSvPos = -1;
        if(strPosPart == "Left")
        {
            iSvPos = al.Position;
        }
        else if(strPosPart == "Right")
        {
            iSvPos = al.Position + al.AlignedBases.length();
        }
        //Check if Sv has been recorde
        bool bFind = false;
        for(map<string, St_MergeSeqAlignStatus>::iterator itr = mpAlignStatus.begin();
            itr != mpAlignStatus.begin(); itr++)
        {
            if(itr->first == strSvName)
            {
                //update the old SV Align Record
                if(strPosPart == "Left")
                {
                    itr->second.iSvLeftAlignPos = iSvPos;
                    itr->second.strAlignSeqSvLeft = al.AlignedBases;
                }
                else if(strPosPart == "Right")
                {
                    itr->second.iSvRightAlignPos = iSvPos;
                    itr->second.strAlignSeqSvRight = al.AlignedBases;
                }
                bFind = true;
            }
            if(bFind)
                break;
        }
        if(!bFind)
        {
            //Add a new SV align record
            stAlignStatus.Clean();
            if(strPosPart == "Left")
            {
                stAlignStatus.iSvLeftAlignPos = iSvPos;
                stAlignStatus.strAlignSeqSvLeft = al.AlignedBases;
            }
            else if(strPosPart == "Right")
            {
                stAlignStatus.iSvRightAlignPos = iSvPos;
                stAlignStatus.strAlignSeqSvRight = al.AlignedBases;
            }
            mpAlignStatus[strSvName] = stAlignStatus;
        }
    }

    //Compare the map created in the last step with the original SvDelReads -->Go
    int iSupportNumOneside = 0;
    int iSupportNumBoth = 0;

    for(vector<St_SvDelReads>::iterator itr = vSvDelReads.begin(); itr != vSvDelReads.end(); itr++)
    {
        string strName = itr->stSV.strType + "_" +
                         to_string(itr->stSV.iPos) + "_" +
                         to_string(itr->stSV.iEnd) + "_SV";
        if(mpAlignStatus.find(strName) != mpAlignStatus.end())
        {
            //Find it
            if(abs(mpAlignStatus[strName].iSvLeftAlignPos - itr->stSV.iEnd) < 10)
            {
                itr->bSvLeftSupportRight = true;
                itr->strSvLeftAlignSeq = mpAlignStatus[strName].strAlignSeqSvLeft;
            }
            if(abs(mpAlignStatus[strName].iSvRightAlignPos - itr->stSV.iPos) < 10)
            {
                itr->bSvRightSupportLeft = true;
                itr->strSvRightAlignSeq = mpAlignStatus[strName].strAlignSeqSvRight;
            }

            if(itr->bSvLeftSupportRight && itr->bSvRightSupportLeft)
                iSupportNumBoth++;
            else if(itr->bSvLeftSupportRight || itr->bSvRightSupportLeft)
                iSupportNumOneside++;
        }
    }

    //Output the supporting result;
    ofstream ofsConfidentResult;
    ofsConfidentResult.open("./ConfidentSupport.txt");

    //*************************************
    //*****Output the statistic result*****
    //*************************************
    ofsConfidentResult << "*************************************" << endl;
    ofsConfidentResult << "Sum SV                : " << to_string(vSvDEL.size()) <<endl;
    ofsConfidentResult << "SV For Current Sample : " << to_string(vSvDelReads.size()) <<endl;
    ofsConfidentResult << "Sum SV Support Oneside: " << to_string(iSupportNumOneside) <<endl;
    ofsConfidentResult << "Sum SV Support Both   : " << to_string(iSupportNumBoth) <<endl;
    ofsConfidentResult << "*************************************" << endl;

    //*********************************************
    //*****Output Detail for each confident SV*****
    //*********************************************
    ofsConfidentResult << endl;
    int iIndex = 1;
    for(vector<St_SvDelReads>::iterator itr = vSvDelReads.begin(); itr != vSvDelReads.end(); itr++)
    {
        if(!itr->bSvLeftSupportRight && !itr->bSvRightSupportLeft) //all of them failed
            continue;

        //Output the support result -->
         ofsConfidentResult << ">SV: " << to_string(iIndex) << endl
                          << "Type                      : " << itr->stSV.strType << endl
                          << "Chrom Index               : " << to_string(itr->stSV.iChrom) << endl
                          << "Position                  : " << to_string(itr->stSV.iPos) << endl
                          << "End Pos                   : " << to_string(itr->stSV.iEnd) << endl
                          << "Length                    : " << to_string(itr->stSV.iLen) << endl
                          << "Support Reads - SV Left   : " << to_string(itr->vSVLeftReads.size()) << endl
                          << "Support Reads - SV Right  : " << to_string(itr->vSVRightReads.size()) << endl
                          << "SvLeftSupportRight        : " << (itr->bSvLeftSupportRight ? "Yes" : "No") << endl
                          << "SvLeftAlignSeq            : " << itr->strSvLeftAlignSeq << endl
                          << "SvRightSupportLeft        : " << (itr->bSvRightSupportLeft ? "Yes" : "No") << endl
                          << "SvRightAlignSeq            : " << itr->strSvRightAlignSeq << endl
                          << "---" << endl;
    }
    ofsConfidentResult.close();

    delete pBamReader;
    pBamReader = NULL;
}

struct St_DiscReads
{
    int iPos;
    int iMatPos;
    int iInsertSize;
    int iSvDiff;
    string strPart;

    St_DiscReads():iPos(-1),iMatPos(-1),iInsertSize(-1), iSvDiff(-1), strPart("Nil")
    {}

    void Clear()
    {
        iPos = -1;
        iMatPos = -1;
        iInsertSize = -1;
        iSvDiff = -1;
        strPart = "Nil";
    }
};

bool sort_discordantReads(St_DiscReads stV1, St_DiscReads stV2)
{
    if(abs(stV1.iInsertSize) > abs(stV2.iInsertSize))
        return true;
    else
        return false;
}

void ClsParseBam::DebugSpeciSv(string strBamFilePath, St_SV& stSv)
{
    //parse bam file
    BamReader* pBamReader = new BamReader();

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    const int iOffSet = 400;
    vector<St_DiscReads> vDiscordantReads;
    St_DiscReads stDiscordantReads;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(!al.IsMapped() || al.QueryBases == "")
            continue;

        //Left part
        if(abs(al.Position + (int)al.AlignedBases.length() - stSv.iPos) < iOffSet)
        {
            //Check if it is discordant
            if(abs(al.InsertSize) > MEDIAN_INSERT_SIZE * 1.2)
            {
                stDiscordantReads.Clear();
                stDiscordantReads.iPos = al.Position;
                stDiscordantReads.iMatPos = al.MatePosition;
                stDiscordantReads.iInsertSize = al.InsertSize;
                stDiscordantReads.iSvDiff = abs(al.Position + (int)al.AlignedBases.length() - stSv.iPos);
                stDiscordantReads.strPart = "Left";
                vDiscordantReads.push_back(stDiscordantReads);
            }
        }

        //right part
        if(abs(al.Position - stSv.iEnd) < iOffSet)
        {
            //Check if it is discordant
            if(abs(al.InsertSize) > MEDIAN_INSERT_SIZE * 1.2)
            {
                stDiscordantReads.Clear();
                stDiscordantReads.iPos = al.Position;
                stDiscordantReads.iMatPos = al.MatePosition;
                stDiscordantReads.iInsertSize = al.InsertSize;
                stDiscordantReads.iSvDiff = abs(al.Position - stSv.iEnd);
                stDiscordantReads.strPart = "Right";
                vDiscordantReads.push_back(stDiscordantReads);
            }
        }
    }
    //Check all the reads around
    cout << "*************************" << endl
         << "*******Support Reads*****" << endl
         << "***SV: " << "(" << to_string(stSv.iPos) << ", " << to_string(stSv.iEnd) << ")" << "***" << endl;

    cout << "Num: " << to_string(vDiscordantReads.size()) << endl;
    cout << "detials" << endl;

    sort(vDiscordantReads.begin(), vDiscordantReads.end(), sort_discordantReads);

    for(vector<St_DiscReads>::iterator itr = vDiscordantReads.begin();
        itr != vDiscordantReads.end(); itr++)
    {
        cout << to_string(itr->iSvDiff) << " ---- " << to_string(itr->iInsertSize) << "   ----   "
             <<"(" << to_string(itr->iPos) << ", " << to_string(itr->iMatPos) << ")"
             << "  ---  "<< stDiscordantReads.strPart << endl;
    }
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl << endl;
    delete pBamReader;
    pBamReader = NULL;
}

bool sort_discordantReads_large_small_func(St_DiscordantReads stDR1, St_DiscordantReads stDR2)
{
    if(stDR1.iReadsPos > stDR2.iReadsPos)
        return true;
    else
        return false;
}

void ClsParseBam::GetDiscordantReads(string strBamFilePath, vector<St_DiscordantReads>& vDiscdRreads)
{
    /* How to do it:
     * 1: Read Bam
     * 2: Check if mapped
     * 3: Check discordant --> absolution value (we allow inversion)
     * 4: Record itself and itsmate
     * 5: Check duplicate
     * 6: End
     */

    //Go!!!
    //parse bam file
    BamReader* pBamReader = new BamReader();

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    St_DiscordantReads stDR;
    vDiscdRreads.clear();

    int iMinThreshold = MEAN_INSERT_SIZE - 3*STD_DEVIATION;
    int iMaxThreshold = MEAN_INSERT_SIZE + 3*STD_DEVIATION;

    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(!al.IsMapped() || al.QueryBases == "")
            continue;
        //if mapped

        if(al.QueryBases.find('N') != string::npos ||
           al.QueryBases.find('n') != string::npos)
            continue;

        // skip concordant reads
        int iRealInsertSize = al.InsertSize;
//        if(iRealInsertSize == 0) //这个情况属于认为的把mate去掉了,但是在bam file里面mate position是有记录的
//        {
//            iRealInsertSize = al.MatePosition - al.Position;
//            if(iRealInsertSize > 0)
//                iRealInsertSize += al.QueryBases.length();
//            else
//                iRealInsertSize -= al.QueryBases.length();
//        }

        if(//abs(iRealInsertSize) > iMinThreshold &&
           abs(iRealInsertSize) < iMaxThreshold)
            continue;

        //skip the case with too large insertion size --> (this too large means abnormal) -->
        int iAbnormalLarge = 100000;
        if(abs(iRealInsertSize) > iAbnormalLarge)
            continue;
        //<--

//        if(abs(iRealInsertSize) > 500000) // it is unrealistic
//            continue;

        if(!al.IsMateMapped()) //both side should be both mapped
            continue;

        //Now they are both mapped discordant reads
        if(//iRealInsertSize >= 0 ||
           al.MatePosition > al.Position) // current reads is left, the mate is right
        {
            //First check if it has been recorded -->
            bool bFind = false;
            for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin();
                itr != vDiscdRreads.end(); itr++)
            {
                if(itr->iReadsPos == al.Position &&
                   itr->iMatePos == al.MatePosition)
                {
                    //update existed value
                    if(itr->strReadsAlignSeq == "")
                    {
                        itr->strReadsAlignSeq = al.AlignedBases;
                        //Check if it contains clip part -->
                        for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                            itrCigar != al.CigarData.end(); itrCigar++)
                        {
                            if(itrCigar->Type == 'S')
                            {
                                itr->bClip = true;
                                break;
                            }
                        }
                        //<--
                    }
                    bFind = true;
                    break;
                }
            }
            if(!bFind)
            {
                //If it is new
                stDR.Clear();
                stDR.iReadsPos = al.Position;
                stDR.iMatePos = al.MatePosition;
                stDR.strReadsAlignSeq = al.AlignedBases;
                stDR.iInsertSize = iRealInsertSize;
                //Check if it contains clip part -->
                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                    itrCigar != al.CigarData.end(); itrCigar++)
                {
                    if(itrCigar->Type == 'S')
                    {
                        stDR.bClip = true;
                        break;
                    }
                }
                //<--
                vDiscdRreads.push_back(stDR);
            }
        }
        else //negative insert size !!!
        {
            //First check if it has been recorded -->
            bool bFind = false;
            for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin();
                itr != vDiscdRreads.end(); itr++)
            {
                if(itr->iReadsPos == al.MatePosition &&
                   itr->iMatePos == al.Position)
                {
                    //update existed value
                    if(itr->strMateAlignSeq == "")
                    {
                        itr->strMateAlignSeq = al.AlignedBases;
                        //Check if it contains clip part -->
                        for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                            itrCigar != al.CigarData.end(); itrCigar++)
                        {
                            if(itrCigar->Type == 'S')
                            {
                                itr->bMateClip = true;
                                break;
                            }
                        }
                        //<--
                    }
                    bFind = true;
                    break;
                }
            }
            if(!bFind)
            {
                //If it is new
                stDR.Clear();
                stDR.iReadsPos = al.MatePosition;
                stDR.iMatePos = al.Position;
                stDR.strMateAlignSeq = al.AlignedBases;
                stDR.iInsertSize = abs(iRealInsertSize);
                //Check if it contains clip part -->
                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                    itrCigar != al.CigarData.end(); itrCigar++)
                {
                    if(itrCigar->Type == 'S')
                    {
                        stDR.bMateClip = true;
                        break;
                    }
                }
                //<--
                vDiscdRreads.push_back(stDR);
            }
        }
    }

    cout << "Discordant Reads Size: " << to_string(vDiscdRreads.size()) << endl;

    //Check how many discordant contain clip reads
    int iLeftClipNum = 0;
    int iRightClipNum = 0;
    int iBothClipNum = 0;
    int iTotalNum = 0;
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin(); itr != vDiscdRreads.end(); itr++)
    {
        if(itr->bClip || itr->bMateClip)
            iTotalNum++;
        if(itr->bClip)
            iLeftClipNum++;
        if(itr->bMateClip)
            iRightClipNum++;
        if(itr->bClip && itr->bMateClip)
            iBothClipNum++;
    }

    cout << "---- Left one Clipped Reads Num: " << IntToStr(iLeftClipNum) << endl;
    cout << "---- Mate Clipped Reads Num    : " << IntToStr(iRightClipNum) << endl;
    cout << "---- Both Clipped Reads Num    : " << IntToStr(iBothClipNum) << endl;
    cout << "---- Total Clipped Reads Num   : " << IntToStr(iTotalNum) << endl;

//    //try to filter  the clipped reads that didn't contain cliped reads --> Comments temporary
//    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.end() - 1;
//        itr >= vDiscdRreads.begin(); itr--)
//    {
//        if(itr->bClip || itr->bMateClip)
//            continue;
//        else
//            vDiscdRreads.erase(itr);
//    }
//    cout << "Filtered vDiscdRreads Size: " << to_string(vDiscdRreads.size()) << endl;

    delete pBamReader;
    pBamReader = NULL;
}

void ClsParseBam::GetClipReads(string strBamFilePath, vector<St_ClipReads>& vClipReads)
{
    /* Steps
     * 1: Parse bam file
     * 2: find softclip
     * 3: record softclip
     *    (1) itself
     *    (2) its mate
     *    (3) insertsize
     */

    //Go!!!
    //Go to check sorted bam file
    //1: parse bam file & output the expected reads
    BamReader* pBamReader = new BamReader();

    vClipReads.clear();
    St_ClipReads stClipReads;

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        //Case 1: do not map  ------->
        if(!al.IsMapped() || al.QueryBases == "")
            continue;

        if(!al.IsMateMapped()) // if the mate doesn't map
            continue;

        // skip concordant reads
        if(abs(al.InsertSize) < MEDIAN_INSERT_SIZE + MIN_DEL_SIZE)
            continue;

        stClipReads.Clear();

        //Case 2: If the reads maps fine
        ///1: Check if there is soft clip
        bool bFindClip = false;
        int iOffSet = 0;

        for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
            itr != al.CigarData.end(); itr++)
        {
            switch(itr->Type)
            {
                case 'M': // alignment match (can be a sequence match or mismatch)
                    iOffSet += itr->Length;
                    break;
                case 'I': // insertion to the reference
                    break;
                case 'D': // deletion from the reference
                case 'N':  // skipped region from the reference
                    break;
                case 'S':  // soft clipping (clipped sequences present in SEQ)
//                    cout << IntToStr(al.QueryBases.length()) << " -- " << to_string(iOffSet) << " : " << to_string(itr->Length)
//                         << endl;
                    if(al.InsertSize > 0) //Left reads
                    {
                        stClipReads.vClipSeq.push_back(al.QueryBases.substr(iOffSet, itr->Length));
                        stClipReads.vPos.push_back(iOffSet);
                    }
                    else //right reads (insert size is negative)
                    {
                        stClipReads.vMateClipSeq.push_back(al.QueryBases.substr(iOffSet, itr->Length));
                        stClipReads.vMatePos.push_back(iOffSet);
                    }
                    iOffSet += itr->Length;
                    bFindClip = true;
                    break;
                case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                    break;
                case 'P': // padding (silent deletion from padded reference)
                case '=': // sequence match
                case 'X': // sequence mismatch
                    break;
            }
        }

        //cout << "-----Ciga Finished: Align VS Query: "
        //     << to_string(al.AlignedBases.length()) << " -- " << to_string(al.QueryBases.length()) << endl;
        //We need to record the clip status of its mate -->
        if(bFindClip)
        {
            bool bFind  = false;
            //First check if it has been recorded
            for(vector<St_ClipReads>::iterator itr = vClipReads.begin(); itr != vClipReads.end(); itr++)
            {
                if(al.InsertSize > 0) //the left reads
                {
                    if(itr->iReadsMapPos == al.Position &&
                       itr->iMateMapPos == al.MatePosition &&
                       itr->iInsertSize == al.InsertSize)
                    {
                        //find it
                        if(itr->strName == "")
                        {
                            itr->strName = al.Name;
                            itr->strQuerySeq = al.QueryBases;
                            itr->strAlignSeq = al.AlignedBases;
                            itr->iReadsMapPos = al.Position;
                            itr->vClipSeq = stClipReads.vClipSeq;
                            itr->vPos = stClipReads.vPos;
                        }
                        bFind = true;
                        break;
                    }
                }
                else //the right reads
                {
                    if(itr->iReadsMapPos == al.MatePosition &&
                       itr->iMateMapPos == al.Position &&
                       itr->iInsertSize == abs(al.InsertSize))
                    {
                        //find it
                        if(itr->strMateName == "")
                        {
                            itr->strMateName = al.Name;
                            itr->strMateQuerySeq = al.QueryBases;
                            itr->strMateAlignSeq = al.AlignedBases;
                            itr->iMateMapPos = al.Position;
                            itr->vMateClipSeq = stClipReads.vMateClipSeq;
                            itr->vMatePos = stClipReads.vMatePos;
                        }
                        bFind = true;
                        break;
                    }
                }
            }

            if(!bFind) //This is a new clip
            {
                if(al.InsertSize > 0)
                {
                    stClipReads.strName = al.Name;
                    stClipReads.iReadsMapPos = al.Position;
                    stClipReads.strQuerySeq = al.QueryBases;
                    stClipReads.strAlignSeq = al.AlignedBases;
                    stClipReads.iInsertSize = al.InsertSize;
                    stClipReads.iMateMapPos = al.MatePosition;
                    vClipReads.push_back(stClipReads);
                }
                else // this reads is mate
                {
                    stClipReads.strMateName = al.Name;
                    stClipReads.iMateMapPos = al.Position;
                    stClipReads.strMateQuerySeq = al.QueryBases;
                    stClipReads.strMateAlignSeq = al.AlignedBases;
                    stClipReads.iInsertSize = abs(al.InsertSize);
                    stClipReads.iReadsMapPos = al.MatePosition;
                    vClipReads.push_back(stClipReads);
                }
            }
        }
    }
    cout << "CR: " << to_string(vClipReads.size()) << endl;
}

void ClsParseBam::GetMultiPlacedMappedReads(string strBamFilePath)
{
    //Try to use hardclip and softclip to identify the related deletion part -->
    vector<St_MultiMapReads> vMultiMapReads;

    BamReader* pBamReader = new BamReader();

    vMultiMapReads.clear();
    St_MultiMapReads stMiltMapReads;

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        //Case 1: do not map  ------->
        if(!al.IsMapped() || al.QueryBases == "")
            continue;

        //Check if it is clip --> Go!!!
        stMiltMapReads.Clear();
        for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
            itr != al.CigarData.end(); itr++)
        {
            switch(itr->Type)
            {
                case 'M': // alignment match (can be a sequence match or mismatch)
                    break;
                case 'I': // insertion to the reference
                    break;
                case 'D': // deletion from the reference
                case 'N':  // skipped region from the reference
                    break;
                case 'S':  // soft clipping (clipped sequences present in SEQ)
                    stMiltMapReads.bSoftClip = true;
                    break;
                case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                    stMiltMapReads.bHardClip = true;
                    break;
                case 'P': // padding (silent deletion from padded reference)
                case '=': // sequence match
                case 'X': // sequence mismatch
                    break;
            }

            if(stMiltMapReads.bSoftClip || stMiltMapReads.bHardClip)
                break;
        }

        if(!stMiltMapReads.bSoftClip && !stMiltMapReads.bHardClip)
            continue; // do another loop diractely

        //Check if it is existed
        bool bExist = false;
        for(vector<St_MultiMapReads>::iterator itr = vMultiMapReads.begin(); itr != vMultiMapReads.end(); itr++)
        {
            //Check the name
            if(itr->strName == stMiltMapReads.strName)
            {
                if(stMiltMapReads.bHardClip)
                {
                    itr->vPos.push_back(al.Position);
                    itr->vClipType.push_back(ctHard);
                    itr->bHardClip = true;
                }

                if(stMiltMapReads.bSoftClip)
                {
                    itr->vPos.push_back(al.Position);
                    itr->vClipType.push_back(ctSoft);
                    itr->bSoftClip = true;
                }
                bExist = true;
                break;
            }
        }
        //If we cannot find it
        if(!bExist)
        {
            if(stMiltMapReads.bHardClip)
            {
                stMiltMapReads.vPos.push_back(al.Position);
                stMiltMapReads.vClipType.push_back(ctHard);
                vMultiMapReads.push_back(stMiltMapReads);
            }

            if(stMiltMapReads.bSoftClip)
            {
                stMiltMapReads.vPos.push_back(al.Position);
                stMiltMapReads.vClipType.push_back(ctSoft);
                vMultiMapReads.push_back(stMiltMapReads);
            }
        }
    }

    delete pBamReader;
    pBamReader = NULL;

    //Check some indicators
    //(1) how many contain multiple mapping positions
    //(2) how many contain both soft and hard clip
    int iMultiMapReadsNum = 0;
    int iBothHardSoftClipReadsNum = 0;

     for(vector<St_MultiMapReads>::iterator itr = vMultiMapReads.begin(); itr != vMultiMapReads.end(); itr++)
     {
         if(itr->vPos.size() > 1)
            iMultiMapReadsNum++;
         if(itr->vPos.size() > 1 && itr->bHardClip && itr->bSoftClip)
            iBothHardSoftClipReadsNum++;
     }
     //Output the result
     cout << endl
          << "Multi_Map_Reads_Num           : " << to_string(iMultiMapReadsNum) << endl
          << "Both_Hard_Soft_Clip_Reads_Num : " << to_string(iBothHardSoftClipReadsNum) << endl;
}

//Cluster discordant reads
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

struct St_DiscdSegment //discordant segment
{
    int iTimes;
    int iStartPos;
    int iEndPos;
    vector<int> vTimes;

    St_DiscdSegment():iTimes(0), iStartPos(-1), iEndPos(-1)
    {}

    St_DiscdSegment(int iV1, int iV2, int iV3)
    {
        iTimes = iV1;
        iStartPos = iV2;
        iEndPos = iV3;
    }

    void Clear()
    {
        iTimes = 0;
        iStartPos = -1;
        iEndPos = -1;
        vTimes.clear();
    }
};

void ClsParseBam::ClusterDS(vector<St_DiscordantReads>& vDiscdRreads, vector<St_SV>& vSvDEL)
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

    //-->Check hit status
    vector<int> vWholeSeq;
    vWholeSeq.resize(iMaxPos+1, 0);

    //count times
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.begin(); itr != vDiscdRreads.end(); itr++)
    {
        if(itr->iReadsPos >= 0 && itr->iMatePos >= 0)
        {
            for_each(vWholeSeq.begin()+itr->iReadsPos, vWholeSeq.begin()+itr->iMatePos, AddVal(1));
        }
    }

    //Get Clustering
    int iCurIndex = 0;
    int iPreValue = -1;
    int iPreIndex = 0;
    int iNonZero = 0;
    vector<St_DiscdSegment> vDiscdSegment;

    for(vector<int>::iterator itr = vWholeSeq.begin(); itr != vWholeSeq.end(); itr++, iCurIndex++)
    {
        if(iPreValue < 0)
        {
            iPreValue = *itr;
            continue;
        }
        else if(iPreValue != *itr) //--> this is new
        {
            //if(iPreValue > 1)
                iNonZero++;

            //if(iPreValue > 1)
                vDiscdSegment.push_back(St_DiscdSegment(iPreValue, iPreIndex, iCurIndex - 1));
            iPreValue = *itr;
            iPreIndex = iCurIndex;
        }
        else if(itr + 1 == vWholeSeq.end()) // if it reaches to the end
        {
            //if(iPreValue > 1)
                vDiscdSegment.push_back(St_DiscdSegment(iPreValue, iPreIndex, iCurIndex));
            //if(iPreValue > 1)
                iNonZero++;
        }
    }

    //Merge the same cluster -->
    const int MINSEGNUM = 5;
    vector<St_DiscdSegment> vMergedDiscdSegment;
    bool bIncrease = false;
    bool bDecrease = false;
    int iPreTimes = 0;

    St_DiscdSegment stDS;

    for(vector<St_DiscdSegment>::iterator itr = vDiscdSegment.begin(); itr != vDiscdSegment.end(); itr++)
    {
        if(itr->iTimes == 0) // this is breakpoint
        {
            if(bIncrease && bDecrease && stDS.vTimes.size() >= MINSEGNUM)
            {
                vMergedDiscdSegment.push_back(stDS);
                stDS.Clear();
                bIncrease = false;
                bDecrease = false;
                iPreTimes = 0;
            }
        }
        else // iTimes != 0
        {
            if(stDS.vTimes.empty()) // this is the first of new time
            {
                stDS.iStartPos = itr->iStartPos;
                iPreTimes = itr->iTimes;
                stDS.vTimes.push_back(itr->iTimes);
                if(stDS.iTimes < itr->iTimes)
                    itr->iTimes = itr->iTimes;
            }
            else //if it isn't the first one
            {
                if(itr->iTimes > iPreTimes) // going up
                {
                    if(!bIncrease)
                    {
                        if(!bDecrease)
                        {
                            bIncrease = true;
                            stDS.vTimes.push_back(itr->iTimes);
                            stDS.iEndPos = itr->iEndPos;
                            iPreTimes = itr->iTimes;
                            if(stDS.iTimes < itr->iTimes)
                                itr->iTimes = itr->iTimes;
                        }
                        else //Decrease is true; this is wrong case
                        {
                            stDS.Clear();
                            bIncrease = false;
                            bDecrease = false;
                            iPreTimes = itr->iTimes;
                        }
                    }
                    else // is increase
                    {
                        if(!bDecrease)
                        {
                            stDS.vTimes.push_back(itr->iTimes);
                            stDS.iEndPos = itr->iEndPos;
                            iPreTimes = itr->iTimes;
                            if(stDS.iTimes < itr->iTimes)
                                itr->iTimes = itr->iTimes;
                        }
                        else // decrease is true;
                        {
                            //this breakpoint --> need to check
                            if(bIncrease && bDecrease && stDS.vTimes.size() >= MINSEGNUM) //Check if could be record
                            {
                                vMergedDiscdSegment.push_back(stDS);
                            }
                            stDS.Clear();
                            stDS.iStartPos = itr->iStartPos;
                            bIncrease = false;
                            bDecrease = false;
                            iPreTimes = itr->iTimes;
                        }
                    }
                }
                else //going down. (Current Times smaller than previous one)
                {
                    if(!bIncrease)
                    {
                        if(!bDecrease)
                        {
                            stDS.Clear();
                            stDS.iStartPos = itr->iStartPos;
                            bIncrease = false;
                            bDecrease = false;
                            iPreTimes = itr->iTimes;
                        }
                        else
                        {
                            stDS.Clear();
                            stDS.iStartPos = itr->iStartPos;
                            bIncrease = false;
                            bDecrease = false;
                            iPreTimes = itr->iTimes;
                        }
                    }
                    else
                    {
                        if(!bDecrease)
                        {
                            stDS.vTimes.push_back(itr->iTimes);
                            stDS.iEndPos = itr->iEndPos;
                            iPreTimes = itr->iTimes;
                            bDecrease = true;
                            if(stDS.iTimes < itr->iTimes)
                                itr->iTimes = itr->iTimes;
                        }
                        else
                        {
                            stDS.vTimes.push_back(itr->iTimes);
                            stDS.iEndPos = itr->iEndPos;
                            iPreTimes = itr->iTimes;
                            if(stDS.iTimes < itr->iTimes)
                                itr->iTimes = itr->iTimes;
                        }
                    }
                }
            }
        }
    }
    //<--

    //Output Merged discordnat segment --> Go !!
    cout << endl << "*****************************" << endl
                 << "vMergedDiscdSegment Size: " << to_string(vMergedDiscdSegment.size()) << endl
                 << "*****************************" << endl;

    for(vector<St_DiscdSegment>::iterator itr = vMergedDiscdSegment.begin(); itr != vMergedDiscdSegment.end(); itr++)
    {
        cout << "Times: " << to_string(itr->iTimes) << "  ---  "
             << "(" << to_string(itr->iStartPos) << ", " << to_string(itr->iEndPos) << ")" << endl;

        for(vector<int>::iterator subItr = itr->vTimes.begin(); subItr != itr->vTimes.end(); subItr++)
        {
            cout << to_string(*subItr) << ", ";
        }
        cout << endl;
    }

    //return;
    //****************************************************************************************

    //Out put result:
    //Statistic
    cout << endl << "****************Statisic*******************" << endl;
    cout << "Total Num: " << to_string(vMergedDiscdSegment.size()) << endl;
    cout << "Non zero : " << to_string(iNonZero) << endl;
    cout << "****************End*******************" << endl;

    //Check the overlap with standard SV
    cout << endl << "****************Compare with SV*******************" << endl;
    int iFindNum = 0;
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        bool bFind = false;
        int iLeftDiff = 1000000;
        int iRightDiff = 1000000;
        int iTimes = 0;

        for(vector<St_DiscdSegment>::iterator subItr = vMergedDiscdSegment.begin(); subItr != vMergedDiscdSegment.end(); subItr++)
        {
            //--> Check if we it can be contained
            if(abs(itr->iPos - subItr->iStartPos) < 400 &&
               abs(itr->iEnd - subItr->iEndPos) < 400)
            {
                bFind = true;

                if(iLeftDiff > abs(itr->iPos - subItr->iStartPos))
                {
                    iLeftDiff = abs(itr->iPos - subItr->iStartPos);
                    iTimes = subItr->iTimes;
                }

                if(iRightDiff > abs(itr->iEnd - subItr->iEndPos))
                {
                    iRightDiff = abs(itr->iEnd - subItr->iEndPos);
                    iTimes = subItr->iTimes;
                }
            }
        }

        if(bFind) //Find it
        {
            cout << "(" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ")"
                 << "  ---  " << "Times: " << to_string(iTimes)
                 << "  ---  " << "Left Diff : " << to_string(iLeftDiff)
                 << "  ---  " << "Right Diff: " << to_string(iRightDiff) << endl;
            iFindNum++;
        }
        else //Fail to find it
        {
            cout << "(" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ")" << "  ---  " << endl;
        }
    }

    cout << "*****************************" << endl
         << "Total SV: " << to_string(vSvDEL.size()) << endl
         << "Find SV : " << to_string(iFindNum) << endl;

    cout << endl << "*************Detail of St_DiscdSegment****************" << endl;

    return;

    //output details
    int iTooShort = 0;
    for(vector<St_DiscdSegment>::iterator itr = vDiscdSegment.begin(); itr != vDiscdSegment.end(); itr++)
    {
        if(itr->iTimes > 0)
        {
            cout  << to_string(itr->iTimes) <<  "  ---  "
                  << "(" << to_string(itr->iStartPos) << ", " << to_string(itr->iEndPos) << ")" <<  "  ---  "
                  << to_string(itr->iEndPos - itr->iStartPos) << endl;
            if((itr->iEndPos - itr->iStartPos) < 100)
                iTooShort++;
        }
    }
    cout << ">>>>>>>>>>>>>>>>>>>>>> Too Short: " << to_string(iTooShort) << endl;
}

bool sort_SV_len_func(St_SV stSv1, St_SV stSv2)
{
    if((stSv1.iEnd - stSv1.iPos) < (stSv2.iEnd - stSv2.iPos))
        return true;
    else
        return false;
}

void ClsParseBam::ClusterDiscrodantReadsDelly(vector<St_DRGroup>& vDRGroup, vector<St_DiscordantReads>& vDiscdRreads,
                                              vector<St_BorderSCR>& vBorderSCR,
                                              vector<St_SV>& vSvDEL)
{       
    //1: Sort current discordant reads data se
    sort(vDiscdRreads.begin(), vDiscdRreads.end(), sort_discordantReads_large_small_func);

//    //Filter the discordant reads which do not contain any split reads --> Go
//    vector<St_BorderSCR>::iterator itrTmpScr = vBorderSCR.begin();
//    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.end() - 1; itr >= vDiscdRreads.begin(); itr--)
//    {
//        for(vector<St_BorderSCR>::iterator subItr = itrTmpScr; subItr < vBorderSCR.end(); subItr++)
//        {
//            if(itr->iReadsPos <= subItr->iClipPos)
//            {
//                if(itr->iMatePos >= subItr->iClipPos)
//                {
//                    //Keep it
//                    itrTmpScr = subItr;
//                }
//                else
//                {
//                    vDiscdRreads.erase(itr);
//                    itrTmpScr = subItr;
//                }
//                break;
//            }
//        }
//    }
//    cout << "After clip reads filter: " << to_string(vDiscdRreads.size()) << endl;
//    //<--

    //2: Make clustering
    if(vDiscdRreads.size() <= 1)
        return;

    int iAllowedDiff = 101; // Now we set this diff as Reads Length    
    St_DRGroup stDRGroup;
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.end() - 1; itr >= vDiscdRreads.begin(); itr--)
    {
        vector<St_DiscordantReads>::iterator tmpItr = itr;
        stDRGroup.vDR.push_back(*itr);
        stDRGroup.iEnd = itr->iMatePos;
        stDRGroup.iStart = itr->iReadsPos;
        stDRGroup.iLeftBoundary = itr->iReadsPos;
        stDRGroup.iRightBoundary = itr->iMatePos;

        for(vector<St_DiscordantReads>::iterator subItr = itr - 1; subItr >= vDiscdRreads.begin(); subItr--)
        {
            if(abs(int(subItr->iReadsPos - itr->iReadsPos)) < iAllowedDiff)
            {
                stDRGroup.vDR.push_back(*subItr);
                tmpItr = subItr;

                //Start should be the largest position, while, End should be the smallest position
                if(stDRGroup.iStart < tmpItr->iReadsPos)
                    stDRGroup.iStart = tmpItr->iReadsPos;

                if(stDRGroup.iEnd > tmpItr->iMatePos)
                    stDRGroup.iEnd = tmpItr->iMatePos;

                //For real boundary -->
                if(stDRGroup.iLeftBoundary > tmpItr->iReadsPos)
                    stDRGroup.iLeftBoundary = tmpItr->iReadsPos;

                if(stDRGroup.iRightBoundary < tmpItr->iMatePos)
                    stDRGroup.iRightBoundary = tmpItr->iMatePos;
                //<--
            }
            else
                break;
        }
        itr = tmpItr;

        //Add current group into group set
        if(stDRGroup.vDR.size() > 1) // erase isolated group (only contributed by 1 PE reads)
        {
            //we only keep the group contributed by more than one PE reads
            vDRGroup.push_back(stDRGroup);
        }

        //Clean current temporary DRGroup
        stDRGroup.Clear();
    }

    //Now we have done !!
    cout << "The Number of DRGroup: " << to_string(vDRGroup.size()) << endl;

    //Combine the overlapped group -->
    cout << endl <<  "---------" << endl;
    vector<St_DRGroup> vNewDRGroup;
    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
    {
//        cout << "( " << to_string(itr->iStart) << ",\t" << to_string(itr->iEnd) << ")" << " --- " << to_string(itr->vDR.size()) << endl;

        stDRGroup.Clear();
        stDRGroup = *itr;

        for(vector<St_DRGroup>::iterator subItr = itr+1; subItr < vDRGroup.end(); subItr++)
        {
            //dropped in the range
            if(subItr->iStart > stDRGroup.iLeftBoundary && subItr->iStart < stDRGroup.iRightBoundary)
            {
                stDRGroup.vDR.insert(stDRGroup.vDR.end(), subItr->vDR.begin(), subItr->vDR.end());
                //Update start and end
                if(subItr->iStart > stDRGroup.iStart)
                    stDRGroup.iStart = subItr->iStart;
                if(subItr->iEnd< stDRGroup.iEnd)
                    stDRGroup.iEnd = subItr->iEnd;

                //update real boundary --> we need to keep the shortest boundary range
                if(subItr->iLeftBoundary > stDRGroup.iLeftBoundary)
                    stDRGroup.iLeftBoundary = subItr->iLeftBoundary;
                if(subItr->iRightBoundary < stDRGroup.iRightBoundary)
                    stDRGroup.iRightBoundary = subItr->iRightBoundary;

                itr = subItr;

                if(subItr == vDRGroup.end())
                {
                     vNewDRGroup.push_back(stDRGroup);
                     break;
                }

            }
            else
            {
                vNewDRGroup.push_back(stDRGroup);
                break;
            }
        }
    }

    cout << "After merge: " << to_string(vNewDRGroup.size()) << endl;

    //Filter the discordant reads which do not contain any split reads --> Go
    vector<St_BorderSCR>::iterator itrTmpScr = vBorderSCR.end() - 1;
    for(vector<St_DRGroup>::iterator itr = vNewDRGroup.end() - 1; itr >= vNewDRGroup.begin(); itr--)
    {
        for(vector<St_BorderSCR>::iterator subItr = itrTmpScr; subItr >= vBorderSCR.begin(); subItr--)
        {
            if(itr->iEnd >= subItr->iClipPos)
            {
                if(itr->iStart <= subItr->iClipPos)
                {                    
                    //Keep it
                    itrTmpScr = subItr;
                    itr->vDRRangedSCR.push_back(*subItr);
                    //Put all of SCR with in the range of current group into the data structure
                    for(vector<St_BorderSCR>::iterator subItr = itrTmpScr - 1; subItr >= vBorderSCR.begin(); subItr--)
                    {
                        if(itr->iStart <= subItr->iClipPos)
                        {
                            itrTmpScr = subItr;
                            itr->vDRRangedSCR.push_back(*subItr);
                        }
                        else
                            break;
                    }
                }
                else
                {
                    vNewDRGroup.erase(itr);
                    itrTmpScr = subItr;
                }
                break;
            }
        }
    }
    cout << "After clip reads filter: " << to_string(vNewDRGroup.size()) << endl;
    //<--

//    //Here --> we can add additional situation to make a better cluser --> Go!! *************************************
//    vector<St_DRGroup> vSCRGroup;
//    GetClusterBySCR(vBorderSCR, vSCRGroup);
//    cout << "spilte reads new group: " << to_string(vSCRGroup.size()) << endl;
//

    vector<St_DRGroup> vSCRGroup;
    GetClusterBySCRTargetReads(vBorderSCR, vSCRGroup);
    //vNewDRGroup.clear();
    vNewDRGroup.insert(vNewDRGroup.end(), vSCRGroup.begin(), vSCRGroup.end());
    //--> add new group
    //<-- ***************************

    for(vector<St_DRGroup>::iterator itr = vNewDRGroup.begin(); itr != vNewDRGroup.end(); itr++)
    {
        cout << "( " << to_string(itr->iStart) << ",\t" << to_string(itr->iEnd) << ")" << " --- "
             << to_string(itr->vDR.size())
             << " --- " << to_string(itr->iEnd - itr->iStart)
             << endl;

    }

    cout <<"---------"<< endl << endl;

    cout << endl << " ============== SV =================" << endl;
    //Check the deletion length of Standard SV
    sort(vSvDEL.begin(), vSvDEL.end(), sort_SV_len_func);
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        cout << "<" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ">"
             << "\t ... SV_Len: " << to_string(itr->iEnd - itr->iPos) << endl;
    }
    cout << "============================" << endl;

    //Check how many hit with SV
    int iHitNum = 0;
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        cout << "<" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ">"
             << " ... SV Len: " << to_string(itr->iEnd - itr->iPos) << endl;
        unsigned int iStartDiff = 500;
        int iRealStartDiff = 0;
        unsigned int iEndDiff = 100000;
        int iRealEndDiff = 0;
        bool bFind = false;
        unsigned int iStart = 0;
        unsigned int iEnd = 0;
        for(vector<St_DRGroup>::iterator subItr = vNewDRGroup.begin(); subItr != vNewDRGroup.end(); subItr++)
        {
            if(itr->iPos >= subItr->iStart-50 && itr->iPos <= subItr->iEnd+50)
            {
                if(abs(int(itr->iPos - subItr->iStart)) < iStartDiff)
                {
                    iStart = subItr->iStart;
                    iEnd = subItr->iEnd;
                    iStartDiff = abs(int(itr->iPos - subItr->iStart));
                    iRealStartDiff = itr->iPos - subItr->iStart;
                    if(iEndDiff > abs(int(itr->iEnd - subItr->iEnd)))
                    {
                        iEndDiff = abs(int(itr->iEnd - subItr->iEnd));
                        iRealEndDiff = subItr->iEnd - itr->iEnd;
                    }
                    bFind = true;
                }
            }
        }

        if(bFind)
        {
            iHitNum++;
            cout << "\tFind --> " << "Start Diff: " << to_string(iRealStartDiff) << " --- "
                 << "End Diff: " << to_string(iRealEndDiff) << " --- "
                 << "(" << to_string(iStart) << ", " << to_string(iEnd) << ")" << endl;
        }
        else
            cout << "\t --- None ---" << endl;
    }

    cout << "Total Number of SV: " << to_string(vSvDEL.size()) << endl;
    cout << "Number of Hit Sv  : " << to_string(iHitNum) << endl;
}

struct St_SCRBoundary
{
    int iClipPos;
    En_ClipPart enClipPart;
    int iMatePos;

    St_SCRBoundary()
    {
        iClipPos = 0;
        enClipPart = cpMax;
        iMatePos = 0;
    }

    void Clear()
    {
        iClipPos = 0;
        enClipPart = cpMax;
        iMatePos = 0;
    }
};

//This logic doesn't work
void ClsParseBam::GetClusterBySCR(vector<St_BorderSCR>& vBorderSCR, vector<St_DRGroup>& vSCRGroup)
{
    //First of all group the bounday of each potential cluster -->
    vector<St_SCRBoundary> vBoundary;
    St_SCRBoundary stBoundary;
    int iRightClipNum = 0;
    int iLeftClipNum = 0;
    int iNilNum = 0;
    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr < vBorderSCR.end(); itr++)
    {
        if(itr->GetClipPart() == cpRight && itr->IsLeftPair()) //
        {
            iRightClipNum++;
            stBoundary.Clear();
            stBoundary.iClipPos = itr->iClipPos;
            stBoundary.enClipPart = itr->GetClipPart();
            stBoundary.iMatePos = itr->iMatPos;
            //Combine the closed and same clip part -->
            for(vector<St_BorderSCR>::iterator subItr = itr+1; subItr < vBorderSCR.end(); subItr++)
            {
                if( subItr->GetClipPart() == cpRight &&
                    subItr->IsLeftPair() &&
                    abs(subItr->iClipPos - stBoundary.iClipPos) < 101)
                {
                    //update boundary info -->
                    //1: update clip pos
                    if(stBoundary.iClipPos < subItr->iClipPos)
                    {
                        stBoundary.iClipPos = subItr->iClipPos;
                    }
                    //2: update mate pos
                    if(stBoundary.iMatePos > subItr->iMatPos)
                    {
                        stBoundary.iMatePos = subItr->iMatPos;
                    }
                    itr = subItr;
                    //<--
                }
                else
                {
                    vBoundary.push_back(stBoundary);
                    break;
                }
            }
        }
        else if(itr->GetClipPart() == cpLeft && itr->IsRightPair())
        {
            iLeftClipNum++;
            stBoundary.Clear();
            stBoundary.iClipPos = itr->iClipPos;
            stBoundary.enClipPart = itr->GetClipPart();
            stBoundary.iMatePos = itr->iMatPos;
            //Combine the closed and same clip part -->
            for(vector<St_BorderSCR>::iterator subItr = itr+1; subItr < vBorderSCR.end(); subItr++)
            {
                if( subItr->GetClipPart() == cpLeft &&
                    subItr->IsRightPair() &&
                    abs(subItr->iClipPos - stBoundary.iClipPos) < 101)
                {
                    //update boundary info -->
                    //1: update clip pos
                    if(stBoundary.iClipPos > subItr->iClipPos)
                    {
                        stBoundary.iClipPos = subItr->iClipPos;
                    }
                    //2: update mate pos
                    if(stBoundary.iMatePos < subItr->iMatPos)
                    {
                        stBoundary.iMatePos = subItr->iMatPos;
                    }
                    itr = subItr;
                    //<--
                }
                else
                {
                    vBoundary.push_back(stBoundary);
                    break;
                }
            }
        }
        else
        {
            iNilNum++;
            continue;
        }
    }

    cout << "Raw Group Number: " << to_string(vBoundary.size()) << endl;
    cout << "\t iRightClipNum: " << to_string(iRightClipNum) << endl;
    cout << "\t  iLeftClipNum: " << to_string(iLeftClipNum) << endl;
    cout << "\t       iNilNum: " << to_string(iNilNum) << endl;

    vSCRGroup.clear();
    St_DRGroup stSCRGroup;
    bool bLeftEdge = false;
    int iLeftEdgeIndex = 0;
    for(vector<St_SCRBoundary>::iterator itr = vBoundary.begin(); itr != vBoundary.end(); itr++)
    {
        if(itr->enClipPart == cpRight)
        {
            if(iLeftEdgeIndex == 0)
            {
                bLeftEdge = true;
                iLeftEdgeIndex = 1;
                stSCRGroup.iStart = itr->iClipPos;
                stSCRGroup.iRightBoundary = itr->iMatePos;
            }
            else if(iLeftEdgeIndex == 1)
            {
                bLeftEdge = true;
                iLeftEdgeIndex = 1;
                stSCRGroup.iStart = itr->iClipPos;
                stSCRGroup.iRightBoundary = itr->iMatePos;
            }
        }
        else if(itr->enClipPart == cpLeft)
        {
            if(bLeftEdge)
            {

                stSCRGroup.iEnd = itr->iClipPos;
                stSCRGroup.iLeftBoundary = itr->iMatePos;
                if(stSCRGroup.GetLen() > MEDIAN_INSERT_SIZE)
                    vSCRGroup.push_back(stSCRGroup);
                bLeftEdge = false;
                stSCRGroup.Clear();
            }
            else
                continue;
        }
        else
            continue;
    }
}

void ClsParseBam::GetClusterBySCRTargetReads(vector<St_BorderSCR>& vBorderSCR, vector<St_DRGroup>& vSCRGroup)
{
    vector<St_BorderSCR> vNewSCR;
    UpgradeSCR(vBorderSCR, vNewSCR);

    //Use target clip reads for clustering complementary --> now we start to make clustering
    //Step 1: Sort first: from small to large based on clip position
    sort(vNewSCR.begin(), vNewSCR.end(), sort_scrClipPos_func);

    //step 2: Start grouping
    int iClipDiff = 150; // the maximum length of
    int iInsetSizeDiff = 200;

    ///1: we first group the boundary reads --> Go
    vector<St_GroupBound> vGroupBound;
    St_GroupBound stGroupBound;
    cout << endl << "Group Info **************" << endl;
    for(vector<St_BorderSCR>::iterator itr = vNewSCR.begin(); itr < vNewSCR.end(); itr++)
    {
        stGroupBound.Clear();
        stGroupBound.vSCR.push_back(*itr);
        En_ClipPart enClipPart = itr->GetClipPart();

        for(vector<St_BorderSCR>::iterator subItr = itr+1; subItr < vNewSCR.end(); subItr++)
        {
            if(abs(itr->iClipPos - subItr->iClipPos) <= iClipDiff &&
               abs(abs(itr->iInsertSize) - abs(subItr->iInsertSize)) < iInsetSizeDiff)
                    //&& subItr->GetClipPart() == enClipPart
            {
                stGroupBound.vSCR.push_back(*subItr);
                itr = subItr;
            }
            else
                break;
        }
        vGroupBound.push_back(stGroupBound);

//        //output --> each group
//        cout << to_string(stGroupBound.vSCR.size()) << "\t";
//        for(vector<St_BorderSCR>::iterator itrSCR = stGroupBound.vSCR.begin();
//            itrSCR != stGroupBound.vSCR.end(); itrSCR++)
//        {
//            cout << to_string(itrSCR->iClipPos) << ", ";
//        }
//        cout << endl;

    }
    cout << "**********************" << endl << endl;

    cout << endl << "vGroupBound Size: " << to_string(vGroupBound.size()) << endl;


    //Filter isolate clip
    for(vector<St_GroupBound>::iterator itr = vGroupBound.end() - 1; itr >= vGroupBound.begin(); itr--)
    {
        if(itr->vSCR.size() <= 1)
            vGroupBound.erase(itr);
    }
    cout << "After Filter: " << to_string(vGroupBound.size()) << endl << endl;

    GetSCRGroup(vGroupBound, vSCRGroup);

//    ///2: calculate the feature for each group --> Update the info of each boundary group
//    for(vector<St_GroupBound>::iterator itr = vGroupBound.begin(); itr != vGroupBound.end(); itr++)
//    {
//        if(itr->vSCR.empty())
//            continue;

//        bool bMate1Start = false;
//        bool bMate1End = false;
//        bool bMate2Start = false;
//        bool bMate2End = false;
//        int iSmallPosClip = itr->vSCR.begin()->iClipPos;
//        int iLargePosClip = itr->vSCR.begin()->iClipPos;
//        for(vector<St_BorderSCR>::iterator subItr = itr->vSCR.begin(); subItr != itr->vSCR.end(); subItr++)
//        {
//            switch(subItr->enMapType)
//            {
//                case mtMat1End:
//                    bMate1End = true;
//                    break;
//                case mtMat1Start:
//                    bMate1Start = true;
//                    break;
//                case mtMat2End:
//                    bMate2End = true;
//                    break;
//                case mtMat2Start:
//                    bMate2Start = true;
//                    break;
//                default:
//                    break;
//            }
//            if(subItr->iClipPos < iSmallPosClip)
//                iSmallPosClip = subItr->iClipPos;
//            if(subItr->iClipPos > iLargePosClip)
//                iLargePosClip = subItr->iClipPos;
//        }

//        if(bMate1End || bMate2End)
//        {
//            itr->bLeftSVBoundGood = true;
//            itr->iClipPos = iLargePosClip;
//        }

//        if(bMate2Start || bMate1Start)
//        {
//            itr->bRightSVBoundGood = true;
//            itr->iClipPos = iSmallPosClip;
//        }
//    }

//    ///3: Do final group
//    vSCRGroup.clear();
//    St_DRGroup stSCRGroup;
//    int iPotentialGroup = 0;
//    for(vector<St_GroupBound>::iterator itr = vGroupBound.begin(); itr != vGroupBound.end(); itr++)
//    {
//        if(itr + 1 == vGroupBound.end())
//            break;

//        if(!itr->bLeftSVBoundGood)
//            continue;

//        if(itr->bLeftSVBoundGood)
//        {
//            if((itr+1)->bRightSVBoundGood) // Left Right --> We get the potential SV
//            {
//                iPotentialGroup++;
//                stSCRGroup.Clear();
//                stSCRGroup.iStart = itr->iClipPos;
//                stSCRGroup.iLeftBoundary = itr->iClipPos;

//                stSCRGroup.iEnd = (itr+1)->iClipPos;
//                stSCRGroup.iRightBoundary = (itr+1)->iClipPos;

//                itr++;

//                if(stSCRGroup.iEnd - stSCRGroup.iStart > 0 &&
//                   stSCRGroup.iEnd - stSCRGroup.iStart <= iMaxSVSize)
//                {
//                    vSCRGroup.push_back(stSCRGroup);
//                }
//                continue;
//            }
//        }
//    }
//    cout << "iPotentialGroup     : " << to_string(iPotentialGroup) << endl;
    cout << "Real Potential Group: " << to_string(vSCRGroup.size()) << endl;
    //<--
}

void ClsParseBam::UpgradeSCR(vector<St_BorderSCR>& vBorderSCR, vector<St_BorderSCR>& vNewSCR)
{
    //Only keep 4 types of softclip and set the type value for each specific one
    vNewSCR.clear();
    int iMat1End = 0;
    int iMat2End = 0;
    int iMat1Start = 0;
    int iMat2Start = 0;
    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr != vBorderSCR.end(); itr++)
    {
        //Type 1
        //Type 3
        if(itr->bFirstMate)
        {
            //Type 1
            if(itr->GetClipPart() == cpRight &&
               itr->bMateMapped &&
               abs(itr->iInsertSize) > MEDIAN_INSERT_SIZE &&
               abs(itr->iInsertSize) < MEDIAN_INSERT_SIZE + 3*STD_DEVIATION) //small deletion
            {
                itr->enMapType = mtMat1End;
                vNewSCR.push_back(*itr);
                iMat1End++;
            }

            //Type 2
            if(itr->GetClipPart() == cpLeft &&
               itr->bMateMapped &&
               abs(itr->iInsertSize) > MEDIAN_INSERT_SIZE - 3*STD_DEVIATION &&
               abs(itr->iInsertSize) < MEDIAN_INSERT_SIZE + 3*STD_DEVIATION)
            {
                itr->enMapType = mtMat1Start;
                vNewSCR.push_back(*itr);
                iMat1Start++;
            }
        }
        //Type 2
        //Type 4
        else if(itr->bSecondMate)
        {
            //Type 2
            if(itr->GetClipPart() == cpLeft &&
               itr->bMateMapped &&
               abs(itr->iInsertSize) > MEDIAN_INSERT_SIZE &&
               abs(itr->iInsertSize) < MEDIAN_INSERT_SIZE + 3*STD_DEVIATION) //small deletion
            {
                itr->enMapType = mtMat2Start;
                vNewSCR.push_back(*itr);
                iMat2Start++;
            }

            //Type 4
            if(itr->GetClipPart() == cpRight &&
               itr->bMateMapped &&
               abs(itr->iInsertSize) > MEDIAN_INSERT_SIZE - 3*STD_DEVIATION &&
               abs(itr->iInsertSize) < MEDIAN_INSERT_SIZE + 3*STD_DEVIATION)
            {
                itr->enMapType = mtMat2End;
                vNewSCR.push_back(*itr);
                iMat2End++;
            }
        }
        else{}
    }
    cout << "vBorderSCR Size: " << to_string(vBorderSCR.size()) << endl;
    cout << "vNewSCR Size   : " << to_string(vNewSCR.size()) << endl;
    cout << "\t iMate1End: " << to_string(iMat1End) << endl;
    cout << "\t iMate1Start: " << to_string(iMat1Start) << endl;
    cout << "\t iMate2End: " << to_string(iMat2End) << endl;
    cout << "\t iMate2Start: " << to_string(iMat2Start) << endl;
}

//Get in current group --> Go
//This time --> we cluster the closed clipped reads together based on clip position
//          and then --> check if the clipped reads dropped in each group could represent any potntial deletion  --> Go!!
void ClsParseBam::GetSCRGroup(vector<St_GroupBound>& vGroupBound, vector<St_DRGroup>& vSCRGroup)
{
    vSCRGroup.clear();

    //--> Go!!
    //1: get all of left side reads and right side reads from the boundary data set
    St_DRGroup stSCRGroup;
    for(vector<St_GroupBound>::iterator itr = vGroupBound.begin(); itr != vGroupBound.end(); itr++)
    {
        stSCRGroup.Clear();
        if(GetGroupFromBoundary(itr->vSCR, stSCRGroup))
        {
            vSCRGroup.push_back(stSCRGroup);
        }
    }
    cout << "vSCRGroup: " << IntToStr(vSCRGroup.size()) << endl;
}

bool ClsParseBam::GetGroupFromBoundary(vector<St_BorderSCR>& vSCR, St_DRGroup& stSCRGroup) // record all of reads which made contribution -->
{
    if(vSCR.size() <= 2)
        return false;

    //Get real left boundary
    //First separate all the reads into 2 groups
    vector<St_BorderSCR> vStartBound;
    vector<St_BorderSCR> vEndBound;
    for(vector<St_BorderSCR>::iterator itr = vSCR.begin(); itr != vSCR.end(); itr++)
    {
        switch(itr->enMapType)
        {
            case mtMat1End:
                vStartBound.push_back(*itr);
                break;
            case mtMat1Start:
                vEndBound.push_back(*itr);
                break;
            case mtMat2End:
                vStartBound.push_back(*itr);
                break;
            case mtMat2Start:
                vEndBound.push_back(*itr);
                break;
            default:
                break;
        }
    }

    cout << "vStartBound_1: " << IntToStr(vStartBound.size()) << endl;
    cout << "vEndBound_1  : " << IntToStr(vEndBound.size()) << endl;

    if(vEndBound.empty() || vStartBound.empty())
        return false;

    //Filter the abnormal clip reads
//    sort(vStartBound.begin(), vStartBound.end(), sort_scrClipPos_func);
//    sort(vEndBound.begin(), vEndBound.end(), sort_scrClipPos_func);
    cout << "sort both StartBound and EndBound Finished!" << endl;
    int iIndex = vStartBound.size()/2;
    int iMidClipPos = vStartBound[iIndex].iClipPos;
    int iOffSet = 100;
    for(vector<St_BorderSCR>::iterator itr = vStartBound.end() - 1; itr >=vStartBound.begin(); itr--)
    {
        if(itr->iClipPos >= iMidClipPos + iOffSet ||
           itr->iClipPos <= iMidClipPos - iOffSet)
            vStartBound.erase(itr);
    }


    iIndex = vEndBound.size()/2;
    if(vEndBound.size() % 2 == 0)
        iIndex--;
    iMidClipPos = vEndBound[iIndex].iClipPos;
    for(vector<St_BorderSCR>::iterator itr = vEndBound.end() - 1; itr >=vEndBound.begin(); itr--)
    {
        if(itr->iClipPos >= iMidClipPos + iOffSet ||
           itr->iClipPos <= iMidClipPos - iOffSet)
            vEndBound.erase(itr);
    }

    cout << "vStartBound_2: " << IntToStr(vStartBound.size()) <<  " >>> ";
    for(vector<St_BorderSCR>::iterator itr = vStartBound.begin(); itr != vStartBound.end(); itr++)
    {
        cout << "<" << to_string(itr->iReadsMapPos) << ", " << to_string(itr->iClipPos) << "> ";
    }
    cout << endl;

    cout << "vEndBound_2  : " << IntToStr(vEndBound.size()) << " >>> ";
    for(vector<St_BorderSCR>::iterator itr = vEndBound.begin(); itr != vEndBound.end(); itr++)
    {
        cout << "<" << to_string(itr->iReadsMapPos) << ", " << to_string(itr->iClipPos) << "> ";
    }
    cout << endl;


    if(vStartBound.size() + vEndBound.size() <= 2)
        return false;

    //Get Clip Pos
    //1: Start
    iIndex = vStartBound.size()/2;
    stSCRGroup.iStart = vStartBound[iIndex].iClipPos;//(vStartBound.end() - 1)->iClipPos;
    stSCRGroup.iLeftBoundary = vStartBound.begin()->iClipPos;

    //2: End
    iIndex = vEndBound.size()/2;
    //if(vEndBound.size() % 2 == 0)
    //    iIndex--;
    stSCRGroup.iEnd = vEndBound[iIndex].iClipPos;//vEndBound.begin()->iClipPos;
    stSCRGroup.iRightBoundary = (vEndBound.end() - 1)->iClipPos;

    cout << "stSCRGroup.iStart: " << IntToStr(stSCRGroup.iStart) << endl;
    cout << "stSCRGroup.iEnd  : " << IntToStr(stSCRGroup.iEnd) << endl;

    if(stSCRGroup.iEnd > stSCRGroup.iStart &&
       abs(stSCRGroup.iEnd - stSCRGroup.iStart) > 5 &&
       abs(stSCRGroup.iEnd - stSCRGroup.iStart) < 150)
    {
        cout << "Yes, I am good" << endl;
        //Save boundary infor into Group Info -->
        stSCRGroup.vStartSCR.insert(stSCRGroup.vStartSCR.end(), vStartBound.begin(), vStartBound.end());
        stSCRGroup.vEndSCR.insert(stSCRGroup.vEndSCR.end(), vEndBound.begin(), vEndBound.end());
        //<--
        return true;
    }
    else
         return false;
}

void ClsParseBam::DrawPics(vector<St_DRGroup>& vDRGroup, string strRefPath, vector<St_SV>& vSvDEL)
{
    if(vDRGroup.empty())
        return;

    cout << endl << " ********************** Draw Pics *************************" << endl;
    //Get the Deletion related Group first
    vector<St_DRGroup> vTargetGroup;
    int iHitNum = 0;

    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        cout << "<" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ">"
             << " ... SV Len: " << to_string(itr->iEnd - itr->iPos) << endl;
        unsigned int iStartDiff = 500;
        int iRealStartDiff = 0;
        unsigned int iEndDiff = 100000;
        int iRealEndDiff = 0;
        bool bFind = false;
        unsigned int iStart = 0;
        unsigned int iEnd = 0;
        for(vector<St_DRGroup>::iterator subItr = vDRGroup.begin(); subItr != vDRGroup.end(); subItr++)
        {
            if(itr->iPos >= subItr->iStart-50 && itr->iPos <= subItr->iEnd+50)
            {
                if(abs(int(itr->iPos - subItr->iStart)) < iStartDiff)
                {
                    iStart = subItr->iStart;
                    iEnd = subItr->iEnd;
                    iStartDiff = abs(int(itr->iPos - subItr->iStart));
                    iRealStartDiff = itr->iPos - subItr->iStart;
                    if(iEndDiff > abs(int(itr->iEnd - subItr->iEnd)))
                    {
                        iEndDiff = abs(int(itr->iEnd - subItr->iEnd));
                        iRealEndDiff = subItr->iEnd - itr->iEnd;
                    }
                    vTargetGroup.push_back(*subItr);
                    bFind = true;
                    break;
                }
            }
        }
        if(bFind)
        {
            iHitNum++;
            cout << "\tFind --> " << "Start Diff: " << to_string(iRealStartDiff) << " --- "
                 << "End Diff: " << to_string(iRealEndDiff) << " --- "
                 << "(" << to_string(iStart) << ", " << to_string(iEnd) << ")" << endl;
        }
        else
            cout << "\t --- None ---" << endl;
    }


    //Read Reference File -->
    ClsFastaReader* pFaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFaReader->ReadFastaRegular(strRefPath, vFasta);
    cout << "Fasta Size: " << to_string(vFasta.size()) << endl;
    //<--

    //Let's Draw Pics
    string strCmd = "mkdir -p ./Image";
    system(strCmd.c_str());
    strCmd = "rm ./Image/*";
    system(strCmd.c_str());

    ClsDrawImage* pDrawImg = new ClsDrawImage();
    int iGroupIndex = 0;

    for(vector<St_DRGroup>::iterator itr = vTargetGroup.begin(); itr != vTargetGroup.end(); itr++)
    {
        DrawPicsForSingleGroup(*itr, *vFasta.begin(), pDrawImg, iGroupIndex);
        iGroupIndex++;
    }

    delete pDrawImg;
    pDrawImg = NULL;
}

void ClsParseBam::DrawPicsForSingleGroup(St_DRGroup& stDRGroup, St_Fasta& stFasta, ClsDrawImage* pDrawImg, int iGroupIndex)
{
    //Get the related part from reference seuqence --> Create Reference Seq
    string strSubRef = stFasta.strSeq.substr(stDRGroup.iLeftBoundary, stDRGroup.GetLen());
    if("" == strSubRef)  // Do nothing if sub-reference is empty
        return;

    ClsBlast* pBlast = new ClsBlast();
    string strRefFa = pBlast->CreatFaFile("RefSeq", strSubRef, ftRef);

    //Now we can draw it -->
    int iRefLen = stDRGroup.iRightBoundary - stDRGroup.iLeftBoundary;
    cv::Mat img = cv::Mat::zeros(iRefLen*2, iRefLen*2, CV_8UC3);

    int iN = 0;
    int iOffset = 0;

    for(vector<St_DiscordantReads>::iterator itr = stDRGroup.vDR.begin(); itr != stDRGroup.vDR.end(); itr++)
    {
        iOffset = pow(-1, iN) * iN;
        DrawPicsForSingleQueryRefSeq(itr->strReadsAlignSeq, strRefFa, pBlast, img, pDrawImg, iOffset, stDRGroup.iLeftBoundary);
        iN++;
    }

    //Left-most first, Right-most last. vector<St_BorderSCR> vDRRangedSCR; // SCR dropped into DR range
    for(vector<St_BorderSCR>::iterator itr = stDRGroup.vDRRangedSCR.begin(); itr != stDRGroup.vDRRangedSCR.end(); itr++)
    {
        iOffset = pow(-1, iN) * iN;
        DrawPicsForSingleQueryRefSeq(itr->strQuerySeq, strRefFa, pBlast, img, pDrawImg, iOffset, stDRGroup.iLeftBoundary);
        iN++;
    }
    //<--

    // -->
    //vector<St_BorderSCR> vStartSCR;  // Start Border of Clipped Reads
    for(vector<St_BorderSCR>::iterator itr = stDRGroup.vStartSCR.begin(); itr != stDRGroup.vStartSCR.end(); itr++)
    {
        iOffset = pow(-1, iN) * iN;
        DrawPicsForSingleQueryRefSeq(itr->strQuerySeq, strRefFa, pBlast, img, pDrawImg, iOffset, stDRGroup.iLeftBoundary);
        iN++;
    }

    //vector<St_BorderSCR> vEndSCR; //End Border of Clipped Reads
    for(vector<St_BorderSCR>::iterator itr = stDRGroup.vEndSCR.begin(); itr != stDRGroup.vEndSCR.end(); itr++)
    {
        iOffset = pow(-1, iN) * iN;
        DrawPicsForSingleQueryRefSeq(itr->strQuerySeq, strRefFa, pBlast, img, pDrawImg, iOffset, stDRGroup.iLeftBoundary);
        iN++;
    }
    //<--

    //Save Image
    string strImagePath = "./Image/Image_Group_" + to_string(iGroupIndex) + ".jpg";
    pDrawImg->SaveImage(img, strImagePath);

    delete pBlast;
    pBlast = NULL;
}

void ClsParseBam::DrawPicsForSingleQueryRefSeq(string& strQuerySeq, string& strRefFa, ClsBlast* pBlast,
                                               cv::Mat& img, ClsDrawImage* pDrawImg, int iOffset, int iBaseLine)
{
    //1: Create Fa for Query Seqeunce (from the reads sequence)
    string strQueryFa = pBlast->CreatFaFile("QuerySeq", strQuerySeq, ftQuery);

    //2: Get Alignment Result
    vector<St_BlastAlignUnit> vBlastAlign;
    pBlast->GetTopNResultsFromTwoSeqBlast(strQueryFa, strRefFa, vBlastAlign, 2);

    //Draw Pics
    St_ImgPos stImgPos;
    stImgPos.iAmplify = 2;
    stImgPos.bPrimery = true;
    stImgPos.iVerticalOffset = iOffset;
    stImgPos.iBaseLine = iBaseLine;

    for(vector<St_BlastAlignUnit>::iterator itr = vBlastAlign.begin(); itr != vBlastAlign.end(); itr++)
    {
        //For the first Line --> Primary Part
        stImgPos.iCurRefStart = itr->iStartPos;
        stImgPos.iCurRefEnd = itr->iEndPos;

        stImgPos.iCurQueryStart = itr->iStartPosQuery;
        stImgPos.iCurQueryEnd = itr->iEndPosQuery;

        pDrawImg->DrawDiagnal(img, stImgPos);
    }
}
