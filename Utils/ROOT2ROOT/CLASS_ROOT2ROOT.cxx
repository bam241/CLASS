/* 

   Simple code used to convert multi root output file in one root file 

Needs : 

Select some element of interest given a Z list

If you are interested in a specific observable you will need to implement the method in Scenar_t.hxx
There is no need to recompile the BIG root file

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
#include "TROOT.h"
//
#include "SSENAR/Scenar_t.hxx"

using namespace std;

// Get the output of a linux command...
string exec(string s_cmd)
{
    char *cmd = new char[s_cmd.length()+1]; strcpy(cmd,s_cmd.c_str());
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

// Show a progress bar ...
    template <class T>
void exec_progress(T i, T Ni)
{
    cout<<"\r =====> "<< (T)((double)i/Ni*100. +1) <<" % "<<flush;
}

//
// Retrieve quantities of interest
//

// // Obtain inventories from facilities or from a isotopic vector pointer
// double GetInvOfTime (CLASSFacility *object, IsotopicVector isovec)
// {
//     return object->GetInsideIV().GetThisComposition(isovec).GetTotalMass();
// }
// //
// double GetInvOfTime (IsotopicVector *object, IsotopicVector isovec)
// {
//     return object->GetThisComposition(isovec).GetTotalMass();
// }
// //
double GetInvAtBOC (Reactor *object, IsotopicVector isovec)
{
    return double(object->GetIVBeginCycle().GetThisComposition(isovec).GetTotalMass());
}
// //
// double GetInvAtEOC (Reactor *object, IsotopicVector isovec)
// {
//     return object->GetIVOutCycle().GetThisComposition(isovec).GetTotalMass();
// }
//
// //
// double GetDiffFlux (CLASSFacility *object, IsotopicVector isovec)
// {
//     double cinflux, coutflux;
//
//     cinflux  = object->GetCumulativeIVIn().GetThisComposition(isovec).GetTotalMass();
//     coutflux = object->GetCumulativeIVOut().GetThisComposition(isovec).GetTotalMass();
//
//     return (cinflux - coutflux);
// }
// //
//

int main(int argc, char** argv)
{

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------INIT---------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    if (argc != 2){
        cout << " Usage : ./ex ROOTFILESDIR " << endl;
        return EXIT_FAILURE;
    }
    //
    string s_ROOTDIR    = argv[1];
    string s_dirtag     = s_ROOTDIR + "/OUT_";
    size_t dirtagsize   = s_dirtag.length();
    size_t numformat    = 9; // nber of bytes used to print the values in the full root files dir names

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------VARIABLES---------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    //  Time and Power will Get dumped from original OUT.root files : respect the TYPE
    Long64_t TimeSecond = 0;
    double ParcPower    = 0.0;
    //


    //  Z list of interest : Only these elements will be dumped in the root file
    vector<unsigned short> ZLIST;//{92,94,95,96};
    ZLIST.push_back(92);
    ZLIST.push_back(94);
    ZLIST.push_back(95);
    ZLIST.push_back(96);

    // Isotopes of Interest
    //
    vector<IsotopicVector> ILIST;//{U,Pu,Am,Cm};
    //
    IsotopicVector v_U,v_U5, v_U8;
    ZAI U8 = ZAI(92,238,0);
    ZAI U5 = ZAI(92,235,0);
    v_U    = 1*U8 + 1*U5;
    ILIST.push_back(v_U);
    //
    v_U5 = 1*U5;
    v_U8 = 1*U8;

    IsotopicVector v_Pu, v_Pu_Fis, v_Pu8, v_Pu9, v_Pu0, v_Pu1, v_Pu2;
    ZAI Pu8 = ZAI(94,238,0);
    ZAI Pu9 = ZAI(94,239,0);
    ZAI Pu0 = ZAI(94,240,0);
    ZAI Pu1 = ZAI(94,241,0);
    ZAI Pu2 = ZAI(94,242,0);
    v_Pu    = 1*Pu8 + 1*Pu9 + 1*Pu0 + 1*Pu1 + 1*Pu2;
    ILIST.push_back(v_Pu);
    //
    v_Pu_Fis = 1*Pu9 + 1*Pu1;
    v_Pu8    = 1*Pu8;
    v_Pu9    = 1*Pu9;
    v_Pu0    = 1*Pu0;
    v_Pu1    = 1*Pu1;
    v_Pu2    = 1*Pu2;

    IsotopicVector v_Am, v_Am1, v_Am3;
    ZAI Am1  = ZAI(95,241,0);
    ZAI Am2x = ZAI(95,242,1);
    ZAI Am3  = ZAI(95,243,0);
    v_Am     = 1*Am1 + 1*Am2x + 1*Am3;
    ILIST.push_back(v_Am);
    //
    v_Am1 = 1*Am1;
    v_Am3 = 1*Am3;

    IsotopicVector v_Cm, v_Cm2, v_Cm4;
    ZAI Cm2 = ZAI(96,242,0);
    ZAI Cm3 = ZAI(96,243,0);
    ZAI Cm4 = ZAI(96,244,0);
    v_Cm    = 1*Cm2 +1*Cm3 + 1*Cm4;
    ILIST.push_back(v_Cm);
    //
    v_Cm2 = 1*Cm2;
    v_Cm4 = 1*Cm4;

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------OUTPUT -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    TFile *FileScenario = new TFile("SSenario.root","RECREATE");
    TTree *TreeScenario = new TTree("TreeScenario","TreeScenario");	

    //


    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------INPUT ------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    // Names of ROOT Files
    string s_Elements   = "find " + s_ROOTDIR + " -type f -name \"*.root\""             ;
    string s_Candidates = "find " + s_ROOTDIR + " -type f ! -size -20M -name \"*.root\"";
    string s_Crashed    = "find " + s_ROOTDIR + " -type f   -size -20M -name \"*.root\"";

    // Number and Names of ROOT Files
    string s_NumberOfElements   = exec(s_Elements   + " | wc -l");
    string s_NumberOfCandidates = exec(s_Candidates + " | wc -l");
    string s_NumberOfCrashed    = exec(s_Crashed    + " | wc -l");


    // Writing in corresponding files
    s_Elements   += " > ROOTFileList.txt";
    s_Candidates += " > ROOTFileCandidates.txt";
    s_Crashed    += " > ROOTFileCrashed.txt";
    system(s_Elements.c_str()  );
    system(s_Candidates.c_str());
    system(s_Crashed.c_str()   );

    unsigned short NumberOfElements = atoi(s_NumberOfElements.c_str());
    unsigned short NumberOfCandidates = atoi(s_NumberOfCandidates.c_str());
    unsigned short NumberOfCrashed = atoi(s_NumberOfCrashed.c_str());

    // Only candidates root files will be analysed
    string s_OneFileForBranches = exec("sed -n 1p ROOTFileCandidates.txt | tr -d '\040\011\012\015'");

    // Load ROOT file number f and load TTree
    TFile *TFileName = new TFile(s_OneFileForBranches.c_str());
    TTree *fData = new TTree();	
    fData = (TTree*) gDirectory->Get(TFileName->GetListOfKeys()->At(TFileName->GetNkeys() - 1)->GetName());
    //Time Steps
    Long64_t NTime = fData->GetEntries();
    TFileName->Close();

    // Branching to Tree
    //Scenar_t* Scen = new Scenar_t(NTime -1);
    Scenar_t* Scen = new Scenar_t();
    TreeScenario->Branch("ssenar","Scenar_t", &Scen);

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------LOAD BRANCHES-----------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    cout<< endl << "#########################" << endl;
    cout<<"Number Of ROOT Files :            " << NumberOfElements   <<endl;
    cout<<"Number Of ROOT Candidates Files : " << NumberOfCandidates <<endl;
    cout<<"Number Of ROOT Crashed Files    : " << NumberOfCrashed    <<endl;
    cout<<"Number Of Time Step / Files     : " << NTime <<endl<<endl;
    cout<<"Progression                     : " << endl;

    string s_ROOTFileName;

    // Failed simulations are dealt by a null score in the observables
    // the Tree is filled by default value which are set to zero
    ifstream f_ROOTFileCrashed("ROOTFileCrashed.txt"); 
    for (int f=1; f<=NumberOfCrashed; f++)
    {
        // Get input arguments
        getline(f_ROOTFileCrashed, s_ROOTFileName);

        Scen->BU_UOx  = atof(s_ROOTFileName.substr(dirtagsize                     ,numformat).c_str());
        Scen->TC_UOx  = atof(s_ROOTFileName.substr(dirtagsize +  1*(numformat + 1),numformat).c_str());
        Scen->BU_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  2*(numformat + 1),numformat).c_str());
        Scen->TC_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  3*(numformat + 1),numformat).c_str());
        Scen->TF_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  4*(numformat + 1),numformat).c_str());
        Scen->IsMOxAm = atoi(s_ROOTFileName.substr(dirtagsize +  5*(numformat + 1),numformat).c_str());
        Scen->Fr_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  6*(numformat + 1),numformat).c_str());
        Scen->NB_Fuel = atof(s_ROOTFileName.substr(dirtagsize +  7*(numformat + 1),numformat).c_str());
        Scen->Ks_Fuel = atof(s_ROOTFileName.substr(dirtagsize +  8*(numformat + 1),numformat).c_str());
        Scen->LF_Fuel = atof(s_ROOTFileName.substr(dirtagsize +  9*(numformat + 1),numformat).c_str());
        Scen->StMMOx  = atoi(s_ROOTFileName.substr(dirtagsize + 10*(numformat + 1),numformat).c_str());

        Scen->TimeStep  = 1.0/12.0; // each month

        TreeScenario->Fill();
        Scen->Clear();

    }f_ROOTFileCrashed.close();

    //
    ifstream f_ROOTFileCandidates("ROOTFileCandidates.txt");
    ofstream f_ROOTFileSelected("ROOTFileSelected.txt");
    //
    Scen->NTimeStep = (unsigned short) (NTime - 1); // total nber of time step, just a caution: should be set by ctor
    Scen->TimeStep  = 1.0/12.0; // each month

    //
    NumberOfCandidates = 1001;
    for (unsigned short f=1; f<=NumberOfCandidates; f++)
    {
        unsigned short UOx_NLOAD_Theoric = 0, MOx_NLOAD_Theoric = 0;
        Scen->UOx_NLOAD = 0;
        Scen->MOx_NLOAD = 0;
        Scen->MOx_MLOAD = 0;

        // Get input arguments
        getline(f_ROOTFileCandidates, s_ROOTFileName);

        Scen->BU_UOx  = atof(s_ROOTFileName.substr(dirtagsize                     ,numformat).c_str());
        Scen->TC_UOx  = atof(s_ROOTFileName.substr(dirtagsize +  1*(numformat + 1),numformat).c_str());
        Scen->BU_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  2*(numformat + 1),numformat).c_str());
        Scen->TC_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  3*(numformat + 1),numformat).c_str());
        Scen->TF_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  4*(numformat + 1),numformat).c_str());
        Scen->IsMOxAm = atoi(s_ROOTFileName.substr(dirtagsize +  5*(numformat + 1),numformat).c_str());
        Scen->Fr_MOx  = atof(s_ROOTFileName.substr(dirtagsize +  6*(numformat + 1),numformat).c_str());
        Scen->NB_Fuel = atof(s_ROOTFileName.substr(dirtagsize +  7*(numformat + 1),numformat).c_str());
        Scen->Ks_Fuel = atof(s_ROOTFileName.substr(dirtagsize +  8*(numformat + 1),numformat).c_str());
        Scen->LF_Fuel = atof(s_ROOTFileName.substr(dirtagsize +  9*(numformat + 1),numformat).c_str());
        Scen->StMMOx  = atoi(s_ROOTFileName.substr(dirtagsize + 10*(numformat + 1),numformat).c_str());

        //
        Scen->S_IsOk = true;

        // Load ROOT file number f and load TTree
        cout<<s_ROOTFileName.c_str()<<endl;
        TFile *TFileName = new TFile(s_ROOTFileName.c_str());

        // avoid corrupted files
        if (TFileName->IsZombie()) {TFileName->Close(); Scen->S_IsOk = false; NumberOfCrashed++; continue;}
        if (!TFileName->IsOpen()) {Scen->S_IsOk = false; NumberOfCrashed++; continue;}

        // Record Selected root files
        f_ROOTFileSelected << s_ROOTFileName << endl;

        //
        TTree *fData = new TTree();	
        fData = (TTree*) gDirectory->Get(TFileName->GetListOfKeys()->At(TFileName->GetNkeys() - 1)->GetName());
        //fData->Print();
        fData->SetBranchStatus("*", 0); // All branches are unbranched

        // CONNECT BRANCHES
        string Branchname, ActiveBranchName;

        fData->SetBranchStatus("AbsTime", 1); 		fData->SetBranchAddress("AbsTime", &TimeSecond);
        fData->SetBranchStatus("ParcPower", 1);		fData->SetBranchAddress("ParcPower", &ParcPower);

        // TOTAL IV
        IsotopicVector *IV_TOTAL=0; Branchname = "TOTAL";
        fData->SetBranchStatus((Branchname+"*").c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &IV_TOTAL);

        // WASTE IV
        IsotopicVector *IV_WASTE=0; Branchname = "WASTE";
        fData->SetBranchStatus((Branchname+"*").c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &IV_WASTE);

        // INCYCLE IV
        IsotopicVector *IV_INCYCLE=0; Branchname = "INCYCLE";
        fData->SetBranchStatus((Branchname+"*").c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &IV_INCYCLE);

        // FACILITIES
        Reactor* PWR_UOx = new Reactor();
        Branchname = "R_PWR_UOx";
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &PWR_UOx);

        Reactor* PWR_MOx = new Reactor();
        Branchname = "R_PWR_MOx";
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &PWR_MOx);

        Storage* StockUOx = new Storage();
        Branchname = "S_StockUOx"; 
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &StockUOx);

        Storage* StockMOx = new Storage();
        Branchname = "S_StockMOx"; 
        ActiveBranchName = Branchname + "*"; 
        fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
        fData->SetBranchAddress((Branchname+".").c_str(), &StockMOx);



        //---------------------------------------------------------------------------------------------------------------------------------------------
        //---------------------------------------------------------------LOOP ON EVENTS AND FILE WRITING-----------------------------------------------
        //---------------------------------------------------------------------------------------------------------------------------------------------

        exec_progress(f,NumberOfCandidates);

        // PROCEEDING
        for (Long64_t  t = 1; t < NTime; t++)	//loop over scenario time
        {
            fData->GetEntry(t);		//Update all branched object to the new CLASS time step t

            double tTime = double(TimeSecond)/double(cYear);
            Scen->Time.push_back(tTime); //The time (in year) at this time step
            Scen->Power.push_back((double) ParcPower);

            //
            Scen->UOx_CT = double(PWR_UOx->GetCycleTime())/double(cYear);
            Scen->MOx_CT = double(PWR_MOx->GetCycleTime())/double(cYear);


            //
            // Cumulated Number of Load
            //

            // New Reactor load - UOx
            if (tTime >= double(PWR_UOx->GetCreationTime())/double(cYear) + Scen->UOx_CT*(UOx_NLOAD_Theoric))
            {
                UOx_NLOAD_Theoric++;
                if (GetInvAtBOC(PWR_UOx, v_U5) > 0)
                {
                    Scen->UOx_NLOAD += 1;
                }
            }

            // New Reactor load - MOx
            if (tTime >= double(PWR_MOx->GetCreationTime())/double(cYear) + Scen->MOx_CT*(MOx_NLOAD_Theoric))
            {
                MOx_NLOAD_Theoric++;
                if (GetInvAtBOC(PWR_MOx, v_Pu) > 0) Scen->MOx_NLOAD += 1;
                else Scen->MOx_MLOAD += 1;
            }
            //

            //
            //
            IsotopicVector* PWR_UOx_BOC  = new IsotopicVector();
            IsotopicVector* PWR_MOx_BOC  = new IsotopicVector();

            IsotopicVector* PWR_UOx_EOC  = new IsotopicVector();
            IsotopicVector* PWR_MOx_EOC  = new IsotopicVector();

            IsotopicVector* StockUOx_CIN = new IsotopicVector();
            IsotopicVector* StockUOx_COU = new IsotopicVector();
            IsotopicVector* StockUOx_INV = new IsotopicVector();

            IsotopicVector* StockMOx_CIN = new IsotopicVector();
            IsotopicVector* StockMOx_COU = new IsotopicVector();
            IsotopicVector* StockMOx_INV = new IsotopicVector();

            IsotopicVector* TOTAL_INV    = new IsotopicVector();
            IsotopicVector* WASTE_INV    = new IsotopicVector();
            IsotopicVector* INCYCLE_INV  = new IsotopicVector();

            //// By ZLIST
            //for (unsigned short i = 0; i < ZLIST.size(); i++){

            //    PWR_UOx_BOC->Add(PWR_UOx->GetIVBeginCycle().GetSpeciesComposition(ZLIST[i]));
            //    PWR_MOx_BOC->Add(PWR_MOx->GetIVBeginCycle().GetSpeciesComposition(ZLIST[i]));

            //    PWR_UOx_EOC->Add(PWR_UOx->GetIVOutCycle().GetSpeciesComposition(ZLIST[i]));
            //    PWR_MOx_EOC->Add(PWR_MOx->GetIVOutCycle().GetSpeciesComposition(ZLIST[i]));

            //    StockUOx_CIN->Add(StockUOx->GetCumulativeIVIn().GetSpeciesComposition(ZLIST[i]));
            //    StockUOx_COU->Add(StockUOx->GetCumulativeIVOut().GetSpeciesComposition(ZLIST[i]));
            //    StockUOx_INV->Add(StockUOx->GetInsideIV().GetSpeciesComposition(ZLIST[i]));

            //    StockMOx_CIN->Add(StockMOx->GetCumulativeIVIn().GetSpeciesComposition(ZLIST[i]));
            //    StockMOx_COU->Add(StockMOx->GetCumulativeIVOut().GetSpeciesComposition(ZLIST[i]));
            //    StockMOx_INV->Add(StockMOx->GetInsideIV().GetSpeciesComposition(ZLIST[i]));

            //      TOTAL_INV->Add(IV_TOTAL->GetSpeciesComposition(ZLIST[i]));
            //      WASTE_INV->Add(IV_WASTE->GetSpeciesComposition(ZLIST[i]));
            //    INCYCLE_INV->Add(IV_INCYCLE->GetSpeciesComposition(ZLIST[i]));

            //}

            // By ILIST
            for (unsigned short i = 0; i < ILIST.size(); i++){

                PWR_UOx_BOC->Add(PWR_UOx->GetIVBeginCycle().GetThisComposition(ILIST[i]));
                PWR_MOx_BOC->Add(PWR_MOx->GetIVBeginCycle().GetThisComposition(ILIST[i]));

                PWR_UOx_EOC->Add(PWR_UOx->GetIVOutCycle().GetThisComposition(ILIST[i]));
                PWR_MOx_EOC->Add(PWR_MOx->GetIVOutCycle().GetThisComposition(ILIST[i]));

                StockUOx_CIN->Add(StockUOx->GetCumulativeIVIn().GetThisComposition(ILIST[i]));
                StockUOx_COU->Add(StockUOx->GetCumulativeIVOut().GetThisComposition(ILIST[i]));
                StockUOx_INV->Add(StockUOx->GetInsideIV().GetThisComposition(ILIST[i]));

                StockMOx_CIN->Add(StockMOx->GetCumulativeIVIn().GetThisComposition(ILIST[i]));
                StockMOx_COU->Add(StockMOx->GetCumulativeIVOut().GetThisComposition(ILIST[i]));
                StockMOx_INV->Add(StockMOx->GetInsideIV().GetThisComposition(ILIST[i]));

                  TOTAL_INV->Add(IV_TOTAL->GetThisComposition(ILIST[i]));
                  WASTE_INV->Add(IV_WASTE->GetThisComposition(ILIST[i]));
                INCYCLE_INV->Add(IV_INCYCLE->GetThisComposition(ILIST[i]));

            }

            //
            Scen->v_PWR_UOx_BOC.push_back(PWR_UOx_BOC);
            Scen->v_PWR_MOx_BOC.push_back(PWR_MOx_BOC);

            Scen->v_PWR_UOx_EOC.push_back(PWR_UOx_EOC);
            Scen->v_PWR_MOx_EOC.push_back(PWR_MOx_EOC);

            Scen->v_StockUOx_CIN.push_back(StockUOx_CIN);
            Scen->v_StockUOx_COU.push_back(StockUOx_COU);
            Scen->v_StockUOx.push_back(StockUOx_INV);

            Scen->v_StockMOx_CIN.push_back(StockMOx_CIN);
            Scen->v_StockMOx_COU.push_back(StockMOx_COU);
            Scen->v_StockMOx.push_back(StockMOx_INV);

            Scen->v_IVTOTAL.push_back(TOTAL_INV);
            Scen->v_IVINCYCLE.push_back(WASTE_INV);
            Scen->v_IVWASTE.push_back(INCYCLE_INV);
            //

            //
        }
        TreeScenario->Fill();
        Scen->Clear();
        TFileName->Close();
    }
    f_ROOTFileCandidates.close();
    f_ROOTFileSelected.close();
    //
    cout<<"ROOT Crashed Files Updated    : " << NumberOfCrashed    <<endl;
    cout<<endl<<"END OF ..."<<endl<<"#########################"<<endl;

    // system("rm -f ROOTFileList.txt");
    // system("rm -f ROOTFileListCrashed.txt");

    FileScenario->Write();
    FileScenario->Close();
}

/*
   g++ -std=c++11 -o CLASS_ROOT2ROOT.ex CLASS_ROOT2ROOT.cxx -I $CLASS_include -L $CLASS_lib -L SCENAR -lCLASSpkg -lScenar `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
   */
