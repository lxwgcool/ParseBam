#include <iostream>
#include "clsparsebam.h"
#include "clsconfig.h"
#include <algorithm>
#include "../../../../ShareLibrary/clsmuscle.h"
#include "../../../../ShareLibrary/clsbwa.h"
#include "clscomparison.h"
#include "clsdebug.h"
#include <map>

using namespace std;
const char* SAMPLENAME = "NA12878";

int main(int argc, char *argv[])
{
    if(argc > 2)
    {
        cout << "To many arguments" << endl;
        return 1;
    }

    unsigned int uiTmp1 = 1000;
    unsigned int uiTmp2 = 1100;
    cout << to_string(int(uiTmp1 - uiTmp2)) << endl;

//    vector<int> vTemp;
//    vTemp.resize(100, 1);
//    for_each(vTemp.begin() + 10, vTemp.begin() + 20, AddVal(1));

    //Read Config file
    ClsConfig* pConfig = new ClsConfig();
    St_Config stConfig;
    pConfig->ReadConfig(stConfig, argv[1]);
    delete pConfig;
    pConfig = NULL;

    //Read VCF
    ClsVcf1000Genome* pClsVcf = new ClsVcf1000Genome();
    pClsVcf->ParseVcf(stConfig.strVcf, SAMPLENAME);
    vector<St_SV> vSvDEL; // only pick out deletion
    pClsVcf->GetDeletion(vSvDEL, SAMPLENAME);
    cout << "vSvDEL size: " << to_string(vSvDEL.size()) << endl;
    delete pClsVcf;
    pClsVcf = NULL;

//    ClsComparison* pComparison = new ClsComparison();
//    //pComparison->CompareStdSVWithPindel(stConfig.strPindelPath, vSvDEL);
//    pComparison->CompareStdSVWithLumpy(stConfig.strLumpyPath, vSvDEL);
//    delete pComparison;
//    pComparison = NULL;
//    return 0;


    //Get Raw SV Candidates
    //CollectRawSvCandidate(stConfig.strBamFile, vSvDEL);

    //return 0;

    //Read Bam file
    ClsParseBam* pClsParseBam = new ClsParseBam();
    ClsDebug* pDebug = new ClsDebug();

//    // (0) Get Multiple Placed Mapped Reads
//    pClsParseBam->GetMultiPlacedMappedReads(stConfig.strBamFile);

    //(1) Get Clipped reads
    vector<St_BorderSCR> vBorderSCR;
    pClsParseBam->ReadBamFile(stConfig.strBamFile, vBorderSCR);    
    cout << "vBorderSCR size: " << to_string(vBorderSCR.size()) << endl;

    //***********************************
    ///Let's do the right thing -->
    //(1) Get discordant reads
    vector<St_DiscordantReads> vDiscdRreads; //DiscdRreads: discordant reads
    pClsParseBam->GetDiscordantReads(stConfig.strBamFile, vDiscdRreads);

    //(2) Sort Discordant reads
    //(3) Make clustering
    //(4) Check how many contain the potential standard SV
    vector<St_DRGroup> vDRGroup;
    pClsParseBam->ClusterDiscrodantReadsDelly(vDRGroup, vDiscdRreads, vBorderSCR, vSvDEL);

    //Now let's draw pics -->
    pClsParseBam->DrawPics(vDRGroup, stConfig.strRef, vSvDEL);

    /// <-----
    //***********************************
    //vector<St_BorderSCR> vNewSCR;
    //pClsParseBam->UpgradeSCR(vBorderSCR, vNewSCR);
    //pDebug->CheckHitStatusByClipReads(vSvDEL, vBorderSCR);

//    pDebug->ClusterClipReads(vBorderSCR, vSvDEL);

    //pDebug->GetChipReadsBySV(vSvDEL, vBorderSCR,  true);

//    string strReadsPath = pDebug->SaveClipReads(vBorderSCR);
//    //--> Go Velvet
//    pDebug->LocalAssemblyByVelvet(strReadsPath);

//    //(2) Get Discordant reads
//    vector<St_DiscordantReads> vDiscdRreads; //DiscdRreads: discordant reads
//    pClsParseBam->GetDiscordantReads(stConfig.strBamFile, vDiscdRreads);
//    //pDebug->CheckHitStatusByDiscordantReads(vSvDEL, vDiscdRreads);

//    /* (3) Score each discordant range
//     * a) Score the range by each discordant reads
//     * b) Check the real SV range to see the max score
//     */
//    cout << endl << "----------------" << endl << "*************** CheckStdSvCount ***************" << endl;

//    //pDebug->CheckStdSvCount(vDiscdRreads, vSvDEL);
//    pDebug->FilterSCReadsByDiscReads(vBorderSCR, vDiscdRreads);
//    pDebug->CheckHitStatusByClipReads(vSvDEL, vBorderSCR);

    //pDebug->GetDiscClipGroup(vBorderSCR, vDiscdRreads);

//    //(3) Get the multiple alignment alignment and check the result
//    ClsMuscle* pMuscle = new ClsMuscle();
//    string strAlignedSeqPath = pMuscle->Run(stConfig.strMultiAlignSeq);
//    string strMergedSeq = pMuscle->MergeAlignedSeq(strAlignedSeqPath);
//    cout << strMergedSeq << endl;
//    delete pMuscle;
//    pMuscle = NULL;

//    //Check The Deletion related reads --> Go!!
//    vector<St_SvDelReads> vSvDelReads;
//    pClsParseBam->GetDELRelatedReads(vSvDEL, vBorderSCR, vSvDelReads);
//    cout << "vSvDelReads size: " << to_string(vSvDelReads.size()) << endl;

//    //output the SvDelReads
//    pClsParseBam->OutputSvDelReads(vSvDelReads);
//    //<--

//    //Align the merged clipped parts back to reference to identify the other side of SV
//    string strReadsPath = "./SvMergedClipSeq.txt";
//    string strBamPath = ClsBWA::GetInstance().CreateBamBySingleReads(stConfig.strRef, strReadsPath,
//                                                 "", "", false, true, true, 12, true); //use 12 threads do it

//    //Let's try to parse this bam file --> Go!!!
//    pClsParseBam->ReadBamFileMergedAlignReads(strBamPath, vSvDelReads, vSvDEL);

    //Release memory
//    delete pClsVcf;
//    pClsVcf = NULL;
    delete pClsParseBam;
    pClsParseBam = NULL;

    delete pDebug;
    pDebug = NULL;

    cout << "Hello World!" << endl;
    return 0;
}


