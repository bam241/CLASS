/* 

Code to use in Sensitivity Analysis
Read and store in a TTree N Scenario information (Input and Output)

STEPS :

1 - Get All ROOT Files names and Store Tree Structure

2 - Build branches In Scenario.root output file

3 - Loop on Input file and Fill output variable. Fill the output tree.



REMAIN

INPUT VAR NOT GOOD
NLOAD and MLOAD NOT GOOD
UNAT Not implemented...

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
#include <regex>


using namespace std;

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------- Convzersion to double ------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
string dtoa(double num)
{
        ostringstream os(ostringstream::out);
        os<<setprecision(6)<<num;
        return os.str();
}
string itoa(int num)
{
        ostringstream os(ostringstream::out);
        os<<num;
        return os.str();
}
//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------- Get the output of a linux command ------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
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

//---------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------- MAIN ------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------USE AND DOC ------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------VARIABLES---------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    string s_tmp; 
    Long64_t TimeSecond = 0; Long64_t Time = 0;
    double Power = 0;
    int NStocks=0; int NPools=0; int NReactors=0; int NFabPlants=0;
    vector <string> v_Branches; // vector that will contain all the branches stored in the TTree

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------DATA TO CHANGE----------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    // Size limits for crashed and correct jobs
    string s_SizeLimit = "14M";

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------IV ---------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    // List Of Isotopes or Elements to store in the final root files
    vector <string> v_Isotopes = {"U","Np","Pu","Am","Cm","MA","Pu8","Pu9","Pu0","Pu1","Pu2","U5","U8","Np7","Am1","Am3","Cm4","Cm5"};
    int NumberOfIsotopes = v_Isotopes.size();

    vector<IsotopicVector> v_IsotopesList;
    // List Of Isotopes or Elements to store in the final root files, same order than above
    IsotopicVector v_U;
    v_U.Add(92,234,0,1);
    v_U.Add(92,235,0,1);
    v_U.Add(92,236,0,1);
    v_U.Add(92,237,0,1);
    v_U.Add(92,238,0,1);
    v_U.Add(92,239,0,1);
    v_IsotopesList.push_back(v_U);

    IsotopicVector v_Np;
    v_Np.Add(93,236,0,1);
    v_Np.Add(93,237,0,1);
    v_Np.Add(93,238,0,1);
    v_Np.Add(93,239,0,1);
    v_IsotopesList.push_back(v_Np);

    IsotopicVector v_Pu;
    v_Pu.Add(94,238,0,1);
    v_Pu.Add(94,239,0,1);
    v_Pu.Add(94,240,0,1);
    v_Pu.Add(94,241,0,1);
    v_Pu.Add(94,242,0,1);
    v_IsotopesList.push_back(v_Pu);

    IsotopicVector v_Am;
    v_Am.Add(95,241,0,1);
    v_Am.Add(95,242,0,1);
    v_Am.Add(95,242,1,1);
    v_Am.Add(95,243,0,1);
    v_Am.Add(95,244,0,1);
    v_IsotopesList.push_back(v_Am);

    IsotopicVector v_Cm;
    v_Cm.Add(96,242,0,1);
    v_Cm.Add(96,243,0,1);
    v_Cm.Add(96,244,0,1);
    v_Cm.Add(96,245,0,1);
    v_Cm.Add(96,246,0,1);
    v_Cm.Add(96,247,0,1);
    v_Cm.Add(96,248,0,1);
    v_IsotopesList.push_back(v_Cm);

    IsotopicVector v_MA = v_Np + v_Am + v_Cm;
    v_IsotopesList.push_back(v_MA);

    IsotopicVector Pu8;
    Pu8.Add(94,238,0,1);
    v_IsotopesList.push_back(Pu8);
    
    IsotopicVector Pu9;
    Pu9.Add(94,239,0,1);
    v_IsotopesList.push_back(Pu9);
    
    IsotopicVector Pu0;
    Pu0.Add(94,240,0,1);
    v_IsotopesList.push_back(Pu0);
    
    IsotopicVector Pu1;
    Pu1.Add(94,241,0,1);
    v_IsotopesList.push_back(Pu1);
    
    IsotopicVector Pu2;
    Pu2.Add(94,242,0,1);
    v_IsotopesList.push_back(Pu2);
    
    IsotopicVector U5;
    U5.Add(92,235,0,1);
    v_IsotopesList.push_back(U5);
    
    IsotopicVector U8;
    U8.Add(92,238,0,1);
    v_IsotopesList.push_back(U8);
    
    IsotopicVector Np7;
    Np7.Add(93,237,0,1);
    v_IsotopesList.push_back(Np7);
    
    IsotopicVector Am1;
    Am1.Add(95,241,0,1);
    v_IsotopesList.push_back(Am1);
    
    IsotopicVector Am3;
    Am3.Add(95,243,0,1);
    v_IsotopesList.push_back(Am3);
    
    IsotopicVector Cm4;
    Cm4.Add(96,244,0,1);
    v_IsotopesList.push_back(Cm4);
    
    IsotopicVector Cm5;
    Cm5.Add(96,245,0,1);
    v_IsotopesList.push_back(Cm5);

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------- CLASS BRANCHES INFO ---------------------------------------------------------
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
    fData->SetBranchStatus("*", 0); // All branches are unbranched
    int NBranches = fData->GetListOfBranches()->GetEntries();

    cout<<"--------------------------------------------------------"<<endl;
    cout<<"--------- TTree " << fData->GetName()<<" Loaded..."<<endl;
    cout<<"--------------------------------------------------------"<<endl;
        cout<<"List of existing Branches : "<<endl<<endl;
    for(int i=0; i<NBranches; i++)
    {
        s_tmp = fData->GetListOfBranches()->At(i)->GetName();
        if(s_tmp[s_tmp.size()-1]=='.') s_tmp = s_tmp.substr(0, s_tmp.size()-1);
        v_Branches.push_back(s_tmp);
        cout<<v_Branches[i]<<endl;
    } 
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    
    // Number of Reactors
    for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="R_") NReactors++;
    Reactor* B_Reactor[NReactors]; int IndiceReactor=0;
    // Number of Pools
    for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="P_") NPools++;
    Pool* B_Pool[NPools]; int IndicePool=0;
    // Number of Stocks
    for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="S_") NStocks++;
    Storage* B_Stock[NStocks]; int IndiceStock=0;
    // Number of Fabrication Plants
    for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="F_") NFabPlants++;
    FabricationPlant* B_FabPlant[NFabPlants]; int IndiceFabPlant=0;
    // Others...
    int NGlobalIV = 7;
    IsotopicVector* B_GlobalIV[NGlobalIV]; int IndiceGlobalIV=0;

    //Time Steps
    Long64_t NTime = fData->GetEntries();
    TFileName->Close();

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------ROOT Files info --------------------------------------------------------------
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
    //---------------------------------------------------------------INPUT BRANCHES TO ADAPT ------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    string s_FileScenarioName = "Scenario.root";

    TFile *FileScenario = new TFile(s_FileScenarioName.c_str(),"RECREATE");
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
    //--------------------------------------------------------------- CREATE BRANCHES IN Scenario FILE and FILL Output Info File ------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    int SIOK = 0;           TreeScenario->Branch("SIOK",&SIOK,"SIOK/I");

    int NumObs = 0;
    vector< vector <double>> Obs;
    vector< double> Obs_t;
    vector<string> NameObs;

    // File.root becomes File.info and is open
    size_t SPos = s_FileScenarioName.find(".root");    s_FileScenarioName.replace(SPos, std::string(".root").length(), ".info");
    ofstream DataFileName(s_FileScenarioName.c_str()); DataFileName<<scientific<<setprecision(5);

    DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
    DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
    DataFileName<<endl;
    DataFileName<<"Structure and Names of Branches Stored in The final TTree"<<endl;
    DataFileName<<"Each Branch is a vector of the time"<<endl<<endl;
    DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
    DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
    DataFileName<<endl<<endl;
    DataFileName<<setw(30)<<"AbsTime        T"<<endl;
    DataFileName<<setw(30)<<"ParcPower      P"<<endl;
    DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
    DataFileName<<setw(30)<<" "<<setw(5);
    for (int i=0; i<NumberOfIsotopes; i++) DataFileName<<v_Isotopes[i]<<setw(5); DataFileName<<endl;
    DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
    
    // Initialisation to zero
    for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 

    // Time and Power
    Obs.push_back(Obs_t);
    Obs.push_back(Obs_t);

    // Fill the Obs Vector with 0, And Fill the info File
    // int NTotalBranches = (NReactors*4 + NPools*3 + NStocks*3 + NFabPlants*3 + 10)*NumberOfIsotopes;
    // for(int b<=2; b<999; b++) Obs.push_back(Obs_t);
    for(int b=2; b<NBranches; b++) 
    {
        //It's a reactor
        if (v_Branches[b].substr(0,2)=="R_")
        {
            DataFileName<<endl;
            DataFileName<<setw(30)<<v_Branches[b]+string(" Inventory");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;

            DataFileName<<setw(30)<<v_Branches[b]+string(" Inv. BOC");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
            
            DataFileName<<setw(30)<<v_Branches[b]+string(" N LOAD"); Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs); DataFileName<<endl;

            DataFileName<<setw(30)<<v_Branches[b]+string(" M LOAD"); Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs); DataFileName<<endl;
        }
        //It's a Pool
        else if (v_Branches[b].substr(0,2)=="P_")
        {
            DataFileName<<endl;
            DataFileName<<setw(30)<<v_Branches[b]+string(" Inventory");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
/*
            DataFileName<<setw(30)<<v_Branches[b]+string(" Cumul In");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;

            DataFileName<<setw(30)<<v_Branches[b]+string(" Cumul Out");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
*/
        }
        //It's a Storage
        else if (v_Branches[b].substr(0,2)=="S_")
        {
            DataFileName<<endl;
            DataFileName<<setw(30)<<v_Branches[b]+string(" Inventory");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
/*
            DataFileName<<setw(30)<<v_Branches[b]+string(" Cumul In");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;

            DataFileName<<setw(30)<<v_Branches[b]+string(" Cumul Out");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
*/
        }
        //It's a FP
        else if (v_Branches[b].substr(0,2)=="F_")
        {
            DataFileName<<endl;
            DataFileName<<setw(30)<<v_Branches[b]+string(" Inventory");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
/*
            DataFileName<<setw(30)<<v_Branches[b]+string(" Cumul In");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;

            DataFileName<<setw(30)<<v_Branches[b]+string(" Cumul Out");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
*/
        }
        // It's global Observable (TOTAL, WASTE, etc...)
        else 
        {    // TOTAL Inventory
            DataFileName<<setw(30)<<v_Branches[b]+string(" Inventory");
            for(int i=0; i<NumberOfIsotopes; i++) {Obs.push_back(Obs_t); NumObs++; DataFileName<<setw(5)<<string("B")+itoa(NumObs);} DataFileName<<endl;
        }
        // Remind the line definition in the matrix
        if(NumObs%5==0)
        {
            DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
            DataFileName<<setw(30)<<" "<<setw(5);
            for (int i=0; i<NumberOfIsotopes; i++) DataFileName<<v_Isotopes[i]<<setw(5); DataFileName<<endl;
            DataFileName<<"C ----------------------------------------------------------------------------------------------------------------------"<<endl;
        }
    }

    // Create Branches in the tree
    TreeScenario->Branch("T",&Obs[0]);
    TreeScenario->Branch("P",&Obs[1]);
    for(int i=1; i<=NumObs; i++) TreeScenario->Branch((string("B")+itoa(i)).c_str(),&Obs[i+1]);
    //for(int b<=2; b<999; b++) TreeScenario->Branch((string("B")+itoa(b-1)).c_str(),&Obs[b]);
//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------START PROCESSING ROOT FILES --------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    cout<<endl<<"#########################"<<endl;
    cout<<"Number Of ROOT Files : "<<NumberOfElements<<endl;
    cout<<"Number Of Time Step / Files : "<<NTime<<endl<<endl;
    cout<<"Progression : "<<endl;

    string s_ROOTFileName;

    // Get the position of the first number...
    int PosFirstData = 0;
    for (int i=0; i<s_OneFileForBranches.length(); i++) if (isdigit(s_OneFileForBranches[i])) {PosFirstData=i; break;}
    int LengthData = 10;

    // Store Bad Scenarios   
    ifstream f_ROOTFileListCrashed("ROOTFileListCrashed.txt"); 
    for (int f=1; f<=NumberOfElementsCrashed; f++)
    {
       // Get input arguments
        getline(f_ROOTFileListCrashed, s_ROOTFileName);
        
        BU_UOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*0 ,LengthData-1).c_str());
        BU_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*1 ,LengthData-1).c_str());
        BU_SFR = atof(s_ROOTFileName.substr(PosFirstData+LengthData*2 ,LengthData-1).c_str());
        Fr_UOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*3 ,LengthData-1).c_str());
        Fr_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*4 ,LengthData-1).c_str());
        Fr_SFR = atof(s_ROOTFileName.substr(PosFirstData+LengthData*5 ,LengthData-1).c_str());
        CT_UOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*6 ,LengthData-1).c_str());
        CT_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*7 ,LengthData-1).c_str());
        CT_SFR = atof(s_ROOTFileName.substr(PosFirstData+LengthData*8 ,LengthData-1).c_str());
        SM_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*9 ,LengthData-1).c_str());
        Fr_SPu = atof(s_ROOTFileName.substr(PosFirstData+LengthData*10,LengthData-1).c_str());
        SIOK = 0;
        TreeScenario->Fill();

    }f_ROOTFileListCrashed.close();

    // Store Good Scenarios   
    ifstream f_ROOTFileList("ROOTFileList.txt");
    for (int f=1; f<=NumberOfElements; f++)
    {
        // int UOX_NLOAD_Theoric = 0; int MOX_NLOAD_Theoric = 0; int SFR_NLOAD_Theoric = 0;
        // UOX_NLOAD = 0; MOX_NLOAD = 0; MOX_MLOAD = 0; UNAT = 0;

        // Get input arguments
        getline(f_ROOTFileList, s_ROOTFileName);

        BU_UOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*0 ,LengthData-1).c_str());
        BU_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*1 ,LengthData-1).c_str());
        BU_SFR = atof(s_ROOTFileName.substr(PosFirstData+LengthData*2 ,LengthData-1).c_str());
        Fr_UOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*3 ,LengthData-1).c_str());
        Fr_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*4 ,LengthData-1).c_str());
        Fr_SFR = atof(s_ROOTFileName.substr(PosFirstData+LengthData*5 ,LengthData-1).c_str());
        CT_UOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*6 ,LengthData-1).c_str());
        CT_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*7 ,LengthData-1).c_str());
        CT_SFR = atof(s_ROOTFileName.substr(PosFirstData+LengthData*8 ,LengthData-1).c_str());
        SM_MOX = atof(s_ROOTFileName.substr(PosFirstData+LengthData*9 ,LengthData-1).c_str());
        Fr_SPu = atof(s_ROOTFileName.substr(PosFirstData+LengthData*10,LengthData-1).c_str());
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

        //---------------------------------------------------------------------------------------------------------------------------------------------
        //---------------------------------------------------------------CONNECT BRANCHES--------------------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------------------------------------
        string Branchname, ActiveBranchName;

        fData->SetBranchStatus("AbsTime", 1);       fData->SetBranchAddress("AbsTime", &TimeSecond);
        fData->SetBranchStatus("ParcPower", 1);     fData->SetBranchAddress("ParcPower", &Power);
    
        for(int i=2; i<NBranches; i++)
        {
            if (v_Branches[i].substr(0,2)=="S_")
            {
                B_Stock[IndiceStock] = new Storage();
                fData->SetBranchStatus((v_Branches[i] + "*").c_str(), 1);
                fData->SetBranchAddress((v_Branches[i] + ".").c_str(), &B_Stock[IndiceStock]);
                IndiceStock++;
            }
            else if (v_Branches[i].substr(0,2)=="P_")
            {
                B_Pool[IndicePool] = new Pool();
                fData->SetBranchStatus((v_Branches[i] + "*").c_str(), 1);
                fData->SetBranchAddress((v_Branches[i] + ".").c_str(), &B_Pool[IndicePool]);
                IndicePool++;
            }
            else if (v_Branches[i].substr(0,2)=="R_")
            {
            
                B_Reactor[IndiceReactor] = new Reactor();
                fData->SetBranchStatus((v_Branches[i] + "*").c_str(), 1);
                fData->SetBranchAddress((v_Branches[i] + ".").c_str(), &B_Reactor[IndiceReactor]);
                IndiceReactor++;
            }
            else if (v_Branches[i].substr(0,2)=="F_")
            {
                B_FabPlant[IndiceFabPlant] = new FabricationPlant();
                fData->SetBranchStatus((v_Branches[i] + "*").c_str(), 1);
                fData->SetBranchAddress((v_Branches[i] + ".").c_str(), &B_FabPlant[IndiceFabPlant]);
                IndiceFabPlant++;
            }
            else
            {
                B_GlobalIV[IndiceGlobalIV] = new IsotopicVector();
                fData->SetBranchStatus((v_Branches[i] + "*").c_str(), 1);                                 //Branch activation
                fData->SetBranchAddress((v_Branches[i] + ".").c_str(), &B_GlobalIV[IndiceGlobalIV]);      //Connection between variable and Branches
                IndiceGlobalIV++;
            }
        }

        //---------------------------------------------------------------------------------------------------------------------------------------------
        //---------------------------------------------------------------LOOP ON EVENTS AND FILE WRITING-----------------------------------------------
        //---------------------------------------------------------------------------------------------------------------------------------------------

        cout<<"\r =====> "<<(int)((double)f/NumberOfElements*100. +1)<<" % "<<flush;

        // Time and Energy
        for (Long64_t  t = 1; t < NTime; t++)   //loop over scenario time
        {
            fData->GetEntry(t);     //Update all branched object to the new CLASS time step j
            Obs[0][t] = double(TimeSecond)/double(cYear); //The time (in year) at this time step
            Obs[1][t] = Power;

            NumObs=2; IndiceReactor=0; IndicePool=0; IndiceStock=0; IndiceFabPlant=0; IndiceGlobalIV=0;

            for(int b=2; b<NBranches; b++)
            {
                //It's a reactor
                bool IsNewCycle=false;
                double R_CT = B_Reactor[IndiceReactor]->GetCycleTime()/double(cYear);
                double R_ST = B_Reactor[IndiceReactor]->GetCreationTime()/double(cYear);
                int NLOAD_Theoric_0 = int((Obs[0][t-1]-R_ST)/R_CT +1);
                int NLOAD_Theoric_1 = int((Obs[0][t]-R_ST)/R_CT +1);
                if (v_Branches[b].substr(0,2)=="R_")
                {
                    // Inventory
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Reactor[IndiceReactor]->GetInsideIV().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
                    // Inv. BOC
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Reactor[IndiceReactor]->GetIVBeginCycle().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
                    // N LOAD

                    Obs[NumObs][t] = Obs[NumObs][t-1];
                    Obs[NumObs+1][t] = Obs[NumObs+1][t-1];

                    if(Obs[0][t] >= R_ST && (NLOAD_Theoric_0<NLOAD_Theoric_1 || t==1)) IsNewCycle = true;

                    if (IsNewCycle)
                    {
                        if (B_Reactor[IndiceReactor]->GetInsideIV().GetThisComposition(v_IsotopesList[0]).GetTotalMass() > 0) Obs[NumObs][t] = Obs[NumObs][t-1] +1;
                        else Obs[NumObs+1][t] = Obs[NumObs+1][t-1] +1;
                    }
                    NumObs+=2;
                    IndiceReactor++;
                }
                //It's a Pool
                else if (v_Branches[b].substr(0,2)=="P_")
                {
                    // Inventory
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Pool[IndicePool]->GetInsideIV().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
/*
                    // Cumul IN
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Pool[IndicePool]->GetCumulativeIVIn().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
                    // Cumul Out
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Pool[IndicePool]->GetCumulativeIVOut().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
*/
                    IndicePool++;
                }
                //It's a Storage
                else if (v_Branches[b].substr(0,2)=="S_")
                {
                    // Inventory
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Stock[IndiceStock]->GetInsideIV().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
 /*
                    // Cumul IN
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Stock[IndiceStock]->GetCumulativeIVIn().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
                    // Cumul Out
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_Stock[IndiceStock]->GetCumulativeIVOut().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
*/
                    IndiceStock++;
                }
                //It's a FP
                else if (v_Branches[b].substr(0,2)=="F_")
                {
                    // Inventory
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_FabPlant[IndiceFabPlant]->GetInsideIV().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
/*
                    // Cumul IN
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_FabPlant[IndiceFabPlant]->GetCumulativeIVIn().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
                    // Cumul Out
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_FabPlant[IndiceFabPlant]->GetCumulativeIVOut().GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
*/    
                    IndiceFabPlant++;
                }
                // It's global Observable (TOTAL, WASTE, etc...)
                else 
                {    
                    for(int i=0; i<NumberOfIsotopes; i++) {Obs[NumObs][t] = B_GlobalIV[IndiceGlobalIV]->GetThisComposition(v_IsotopesList[i]).GetTotalMass(); NumObs++;}
                    IndiceGlobalIV++;
                }
            }
        }

        TreeScenario->Fill();
        TFileName->Close();

        // Free momory...
        // Delete pointers
        for(int i=0; i<NReactors; i++)  delete B_Reactor[i];
        for(int i=0; i<NPools; i++)     delete B_Pool[i];
        for(int i=0; i<NStocks; i++)    delete B_Stock[i];
        for(int i=0; i<NFabPlants; i++) delete B_FabPlant[i];
        for(int i=0; i<NGlobalIV; i++)  delete B_GlobalIV[i];
        // Indices to zero...
        IndiceReactor=0; IndicePool=0; IndiceStock=0; IndiceFabPlant=0; IndiceGlobalIV=0;
    }
    f_ROOTFileList.close();
    cout<<endl<<"END OF ..."<<endl<<"#########################"<<endl;

    system("rm -f ROOTFileList.txt");
    system("rm -f ROOTFileListCrashed.txt");

    FileScenario->Write();
    FileScenario->Close();

}

/*
g++ -std=c++11 -o CLASS_R2R CLASS_ROOT2ROOT.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
*/
