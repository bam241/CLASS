/* 

Simple code used to convert root output file produced by CLASS to simple data file

Remains : 

First Cycle Time

Authors:

BaL
Nico. T.
ZaK

*/

#include "CLASSHeaders.hxx"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
using namespace std;

// Get the output of a linux command...
string exec(const char* cmd)
{
    char buffer[128];
    string result = "";
    shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
}

int main(int argc, char** argv)
{
//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------TO ADAPT----------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    // Size limits for crashed and correct jobs
    string s_SizeLimit = "14M";

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------VARIABLES---------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    string s_tmp; 
    Long64_t TimeSecond = 0;

    IsotopicVector v_Pu;
    v_Pu.Add(94,238,0,1);
    v_Pu.Add(94,239,0,1);
    v_Pu.Add(94,240,0,1);
    v_Pu.Add(94,241,0,1);
    v_Pu.Add(94,242,0,1);
    v_Pu.Add(95,241,0,1);

    IsotopicVector v_U;
    v_U.Add(92,238,0,1);
    v_U.Add(92,235,0,1);

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------BRANCHES ---------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    TFile *FileScenario = new TFile("Scenario.root","RECREATE");
    TTree *TreeScenario = new TTree("TreeScenario","TreeScenario"); 

    double  BU_UOX  = 0;    TreeScenario->Branch("BU_UOX",&BU_UOX,"BU_UOX/D");
    double  BU_MOX  = 0;    TreeScenario->Branch("BU_MOX",&BU_MOX,"BU_MOX/D");
    double  BU_SFR  = 0;    TreeScenario->Branch("BU_SFR",&BU_SFR,"BU_SFR/D");

    double  Fr_UOX  = 0;    TreeScenario->Branch("Fr_UOX",&Fr_UOX,"Fr_UOX/D");
    double  Fr_MOX  = 0;    TreeScenario->Branch("Fr_MOX",&Fr_MOX,"Fr_MOX/D");
    double  Fr_SFR  = 0;    TreeScenario->Branch("Fr_SFR",&Fr_SFR,"Fr_SFR/D");

    double  CT_UOX  = 0;    TreeScenario->Branch("CT_UOX",&CT_UOX,"CT_UOX/D");
    double  CT_MOX  = 0;    TreeScenario->Branch("CT_MOX",&CT_MOX,"CT_MOX/D");
    double  CT_SFR  = 0;    TreeScenario->Branch("CT_SFR",&CT_SFR,"CT_SFR/D");

    int     SM_MOX  = 0;    TreeScenario->Branch("SM_MOX",&SM_MOX,"SM_MOX/I");

    double  Fr_SPu  = 0;    TreeScenario->Branch("Fr_SPu",&Fr_SPu,"Fr_SPu/D");

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------INPUT ------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    string PathToROOTFiles = string(argv[1]);

    // Number and Names of correct ROOT Files
    string CMD = string("find -L ") + PathToROOTFiles+ string(" -type f -size +") + s_SizeLimit + string(" -name \"OUT.root\"");

    string s_NumberOfElements = exec((CMD + string(" | wc -l")).c_str()); int NumberOfElements = atoi(s_NumberOfElements.c_str()); 
    system((CMD + string(" > ROOTFileList.txt")).c_str());

    // Number and Names of crashed ROOT Files
    CMD = string("find -L ") + PathToROOTFiles+ string(" -type f -size -") + s_SizeLimit + string(" -name \"OUT.root\"");
    string s_NumberOfElementsCrashed = exec((CMD + string(" | wc -l")).c_str()); int NumberOfElementsCrashed = atoi(s_NumberOfElementsCrashed.c_str());
    system((CMD + string(" > ROOTFileListCrashed.txt")).c_str());

    //Get first line for root file example
    string s_OneFileForBranches = exec("sed -n 1p ROOTFileList.txt | tr -d '\040\011\012\015'");
    
    // Load ROOT file number f and load TTree
    TFile *TFileName = new TFile(s_OneFileForBranches.c_str());
    TTree *fData = new TTree(); 
    fData = (TTree*) gDirectory->Get(TFileName->GetListOfKeys()->At(TFileName->GetNkeys() - 1)->GetName());
    //Time Steps
    Long64_t NTime = fData->GetEntries();
    TFileName->Close();

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------WRITING ----------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    cout<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"--------- EXPERIMENT INFORMATIONS ----------------------"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;
    cout<<" Number of Good ROOT files : "<<NumberOfElements<<endl;
    cout<<" Number of Crashed ROOT Files : "<<NumberOfElementsCrashed<<endl;
    cout<<" Number of Time Steps / Simulation : "<<NTime<<endl;
    cout<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------LOAD BRANCHES-----------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    cout<<endl<<"#########################"<<endl;
    cout<<"Number Of ROOT Files : "<<NumberOfElements<<endl;
    cout<<"Number Of Time Step / Files : "<<NTime<<endl<<endl;
    cout<<"Progression : "<<endl;

    string s_ROOTFileName;
 
    // Store Bad Scenarios   
    ifstream f_ROOTFileListCrashed("ROOTFileListCrashed.txt"); 
    for (int f=1; f<=NumberOfElementsCrashed; f++)
    {
       // Get input arguments
        getline(f_ROOTFileListCrashed, s_ROOTFileName);
        
        BU_UOX = atof(s_ROOTFileName.substr(6,9).c_str());
        BU_MOX = atof(s_ROOTFileName.substr(16,9).c_str());
        BU_SFR = atof(s_ROOTFileName.substr(26,9).c_str());
        Fr_UOX = atof(s_ROOTFileName.substr(36,9).c_str());
        Fr_MOX = atof(s_ROOTFileName.substr(46,9).c_str());
        Fr_SFR = atof(s_ROOTFileName.substr(56,9).c_str());
        CT_UOX = atof(s_ROOTFileName.substr(66,9).c_str());
        CT_MOX = atof(s_ROOTFileName.substr(76,9).c_str());
        CT_SFR = atof(s_ROOTFileName.substr(86,9).c_str());
        SM_MOX = atof(s_ROOTFileName.substr(96,9).c_str());
        Fr_SPu = atof(s_ROOTFileName.substr(106,9).c_str());
        SIOK = 0;
        TreeScenario->Fill();

    }f_ROOTFileListCrashed.close();

    // Store Bad Scenarios   
    ifstream f_ROOTFileList("ROOTFileList.txt");
    for (int f=1; f<=NumberOfElements; f++)
    {
        int UOX_NLOAD_Theoric = 0; int MOX_NLOAD_Theoric = 0; int SFR_NLOAD_Theoric = 0;
        UOX_NLOAD = 0; MOX_NLOAD = 0; MOX_MLOAD = 0; UNAT = 0;

        // Get input arguments
        getline(f_ROOTFileList, s_ROOTFileName);

        BU_UOX = atof(s_ROOTFileName.substr(6,9).c_str());
        BU_MOX = atof(s_ROOTFileName.substr(16,9).c_str());
        BU_SFR = atof(s_ROOTFileName.substr(26,9).c_str());
        Fr_UOX = atof(s_ROOTFileName.substr(36,9).c_str());
        Fr_MOX = atof(s_ROOTFileName.substr(46,9).c_str());
        Fr_SFR = atof(s_ROOTFileName.substr(56,9).c_str());
        CT_UOX = atof(s_ROOTFileName.substr(66,9).c_str());
        CT_MOX = atof(s_ROOTFileName.substr(76,9).c_str());
        CT_SFR = atof(s_ROOTFileName.substr(86,9).c_str());
        SM_MOX = atof(s_ROOTFileName.substr(96,9).c_str());
        Fr_SPu = atof(s_ROOTFileName.substr(106,9).c_str());
        SIOK = 1;

        // Load ROOT file number f and load TTree
        cout<<s_ROOTFileName.c_str()<<endl;
        TFile *TFileName = new TFile(s_ROOTFileName.c_str());
        if (TFileName->IsZombie()) {TFileName->Close(); continue;}
        if (!TFileName->IsOpen()) {continue;}
        TTree *fData = new TTree(); 
        fData = (TTree*) gDirectory->Get(TFileName->GetListOfKeys()->At(TFileName->GetNkeys() - 1)->GetName());
        //fData->Print();
        fData->SetBranchStatus("*", 0); // All branches are unbranched

        // CONNECT BRANCHES
        string Branchname, ActiveBranchName;

        fData->SetBranchStatus("AbsTime", 1);       fData->SetBranchAddress("AbsTime", &TimeSecond);
        fData->SetBranchStatus("ParcPower", 1);     fData->SetBranchAddress("ParcPower", &Power);

        IsotopicVector *IV_TOTAL=0; Branchname = "TOTAL";
        fData->SetBranchStatus((Branchname+"*").c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &IV_TOTAL);

        Reactor* PWR_UOX = new Reactor();
        Branchname = "R_PWR_UOx";
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &PWR_UOX);

        Reactor* PWR_MOX = new Reactor();
        Branchname = "R_PWR_MOX";
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &PWR_MOX);

        Reactor* SFR_MOX = new Reactor();
        Branchname = "R_SFR_MOX";
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &SFR_MOX);

        Storage* Stock_UOX = new Storage();
        Branchname = "S_StockUOx"; 
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &Stock_UOX);

        Storage* Stock_MOX = new Storage();
        Branchname = "S_StockMOx"; 
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &Stock_MOX);

        //---------------------------------------------------------------------------------------------------------------------------------------------
        //---------------------------------------------------------------LOOP ON EVENTS AND FILE WRITING-----------------------------------------------
        //---------------------------------------------------------------------------------------------------------------------------------------------

        cout<<"\r =====> "<<(int)((double)f/NumberOfElements*100. +1)<<" % "<<flush;

        for (Long64_t  t = 1; t < NTime; t++)   //loop over scenario time
        {
            fData->GetEntry(t);     //Update all branched object to the new CLASS time step j

            Time = double(TimeSecond)/double(cYear); //The time (in year) at this time step

 
            TreeScenario->Fill();
        }
        TFileName->Close();
    }
    f_ROOTFileList.close();
    cout<<endl<<"END OF ..."<<endl<<"#########################"<<endl;

    system("rm -f ROOTFileList.txt");
    system("rm -f ROOTFileListCrashed.txt");

    FileScenario->Write();
    FileScenario->Close();
}

/*
 g++ -std=c++11 -o CLASS_R2D CLASS_ROOT2DAT.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
*/
