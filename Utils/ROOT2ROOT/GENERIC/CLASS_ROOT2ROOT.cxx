/* 

Code to use in Sensitivity Analysis
Read and store in a TTree N Scenario information (Input and Output)

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
    vector<IsotopicVector *> IV_Branch;
    int NStocks=0; int NPools=0; int NReactors=0; int NFabPlants=0;
    vector <string> v_Branches; // vector that will contain all the branches stored in the TTree

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------DATA TO CHANGE----------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    // Size limits for crashed and correct jobs
    string s_SizeLimit = "14M";

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

    //Time Steps
    Long64_t NTime = fData->GetEntries();
    TFileName->Close();

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------INPUT BRANCHES TO ADAPT ------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    string s_FileScenarioName = "Scenario.root";
    string s_TreeScenarioName = "TT";

    TFile *FileScenario = new TFile(s_FileScenarioName.c_str(),"RECREATE");
    TTree *TreeScenario = new TTree(s_TreeScenarioName.c_str(),s_TreeScenarioName.c_str()); 

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
    //---------------------------------------------------------------IV ---------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    // List Of Isotopes or Elements to store in the final root files
    vector <string> v_Isotopes = {"U","Np","Pu","Am","Cm","MA","PF","Pu8","Pu9","Pu0","Pu1","Pu2","U5","U8"};
    int NumberOfIsotopes = v_Isotopes.size();

    //---------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------- CREATE BRANCHES IN Scenario FILE --------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    int SIOK = 0;           TreeScenario->Branch("SIOK",&SIOK,"SIOK/I");

    vector<double> v_Time;  TreeScenario->Branch("T",&v_Time);
    vector<double> v_Power; TreeScenario->Branch("P",&v_Power);


    int NumObs = 1;

    vector< vector <double>> Obs;
    vector< double> Obs_t;
    vector<string> NameObs;

    // Branches Vector Initialization
    for(int b=2; b<NBranches; b++)
    {
        //It's a reactor
        if (v_Branches[b].substr(0,2)=="R_")
        {
            // Inventory in rector
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
            // Inventory @ BOC in rector
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
            // Fissile Enrichment in rector
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
            // Number of Load in rector
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
            // Number of Missed Load in rector
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
        }
        //It's a Storage, A Pool or a FP
        else if (v_Branches[b].substr(0,2)=="P_" || v_Branches[b].substr(0,2)=="S_" || v_Branches[b].substr(0,2)=="F_")
        {
            // Inventory in Facility
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
            // Cumul In In facility 
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
            // Cumul Out In facility 
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
        }
        // It's global Observable (TOTAL, WASTE, etc...)
        else 
        {    // Inventory
            for(int i=0; i<NumberOfIsotopes; i++)
            {
                for(int t=0; t<NTime; t++) Obs_t.push_back(0.); 
                Obs.push_back(Obs_t);
                NameObs.push_back(string("O") + itoa(NumObs)); NumObs++;
            }
        }
    }


cout<<Obs.size()<<endl;
cout<<NameObs.size()<<endl;
exit(1);









/*

    for(int i=2; i<NBranches; i++)
    {
        for(int e=0; e<NumberOfIsotopes; e++)
        {
            // ADD also BOC, NLOAD, MLOAD and Enr if Reactor
            // ADD also CumulIn and CumulOut if Others
            // 
            NumObs++;
            string s_Obs = string("O") + itoa(NumObs);

            if (v_Branches[i].substr(0,2)=="R_")
            {

                TreeScenario->Branch(s_Obs.c_str(),&Obs[NumObs]);



               
                vector<double> R_Inv;       TreeScenario->Branch("O1",&R_Inv);
                vector<double> R_InvBOC;    TreeScenario->Branch("O2",&R_InvBOC);
                vector<int> R_NLOAD;        TreeScenario->Branch("O3",&R_NLOAD);
                vector<int> R_MLOAD;        TreeScenario->Branch("O4",&R_MLOAD);
                vector<int> R_ENR;          TreeScenario->Branch("O5",&R_ENR);

                R_Inv[i][e] = vector<double> doit etre chargé

           
            }



        }
    }
*/
    //---------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------- Fill OUTPUT MATRIX In Info File----------------------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    size_t SPos = s_FileScenarioName.find(".root");    s_FileScenarioName.replace(SPos, std::string(".root").length(), ".info");
    ofstream DataFileName(s_FileScenarioName.c_str()); DataFileName<<scientific<<setprecision(5);
    
    DataFileName<<"C AbsTime        T"<<endl;
    DataFileName<<"C ParcPower      P"<<endl;
    DataFileName<<"C ------------------------------------------------------------------------------------------"<<endl;
    DataFileName<<"C "<<setw(20)<<" "<<setw(5)<<"U"<<setw(5)<<"Np"<<setw(5)<<"Pu"<<setw(5)<<"Am"<<setw(5)<<"Cm"<<setw(5)<<"MA"<<setw(5)<<"PF"<<setw(5)<<"Pu8"<<setw(5)<<"Pu9"<<setw(5)<<"Pu0"<<setw(5)<<"Pu1"<<setw(5)<<"Pu2"<<setw(5)<<"U5"<<setw(5)<<"U8"<<endl;
    DataFileName<<"C ------------------------------------------------------------------------------------------"<<endl;
    for(int i=2; i<NBranches; i++)
    {
        DataFileName<<"C "<<setw(20)<<v_Branches[i]; 
        //DataFileName.seekp((i-1) * 100 + 20);
        for(int e=0; e<NumberOfIsotopes; e++)
        {
            DataFileName<<setw(5)<<2 + NumberOfIsotopes*(i-2) + e;
        }DataFileName<<endl;
        if(i%6==0)
        {
            DataFileName<<"C ------------------------------------------------------------------------------------------"<<endl;
            DataFileName<<"C "<<setw(20)<<" "<<setw(5)<<"U"<<setw(5)<<"Np"<<setw(5)<<"Pu"<<setw(5)<<"Am"<<setw(5)<<"Cm"<<setw(5)<<"MA"<<setw(5)<<"PF"<<setw(5)<<"Pu8"<<setw(5)<<"Pu9"<<setw(5)<<"Pu0"<<setw(5)<<"Pu1"<<setw(5)<<"Pu2"<<setw(5)<<"U5"<<setw(5)<<"U8"<<endl;
            DataFileName<<"C ------------------------------------------------------------------------------------------"<<endl;
        }
    }
    DataFileName<<"C ------------------------------------------------------------------------------------------"<<endl;
    DataFileName<<"C"<<endl;
    DataFileName<<"C"<<endl;
    DataFileName<<"C"<<endl;

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
    //--------------------------------------------------------------- IV needed for getting info on rooti files -----------------------------------
    //---------------------------------------------------------------------------------------------------------------------------------------------

    IsotopicVector v_U;
    v_U.Add(92,234,0,1);
    v_U.Add(92,235,0,1);
    v_U.Add(92,236,0,1);
    v_U.Add(92,237,0,1);
    v_U.Add(92,238,0,1);
    v_U.Add(92,239,0,1);

    IsotopicVector v_Np;
    v_Np.Add(93,236,0,1);
    v_Np.Add(93,237,0,1);
    v_Np.Add(93,238,0,1);
    v_Np.Add(93,239,0,1);

    IsotopicVector v_Pu;
    v_Pu.Add(94,238,0,1);
    v_Pu.Add(94,239,0,1);
    v_Pu.Add(94,240,0,1);
    v_Pu.Add(94,241,0,1);
    v_Pu.Add(94,242,0,1);

    IsotopicVector v_Am;
    v_Am.Add(95,241,0,1);
    v_Am.Add(95,242,0,1);
    v_Am.Add(95,242,1,1);
    v_Am.Add(95,243,0,1);
    v_Am.Add(95,244,0,1);

    IsotopicVector v_Cm;
    v_Cm.Add(96,242,0,1);
    v_Cm.Add(96,243,0,1);
    v_Cm.Add(96,244,0,1);
    v_Cm.Add(96,245,0,1);
    v_Cm.Add(96,246,0,1);
    v_Cm.Add(96,247,0,1);
    v_Cm.Add(96,248,0,1);

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------START PROCESSING ROOT FILES --------------------------------------------------
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
        // int UOX_NLOAD_Theoric = 0; int MOX_NLOAD_Theoric = 0; int SFR_NLOAD_Theoric = 0;
        // UOX_NLOAD = 0; MOX_NLOAD = 0; MOX_MLOAD = 0; UNAT = 0;

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

        //---------------------------------------------------------------------------------------------------------------------------------------------
        //---------------------------------------------------------------CONNECT BRANCHES--------------------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------------------------------------
        string Branchname, ActiveBranchName;

        fData->SetBranchStatus("AbsTime", 1);       fData->SetBranchAddress("AbsTime", &TimeSecond);
        fData->SetBranchStatus("ParcPower", 1);     fData->SetBranchAddress("ParcPower", &Power);

        for(int i=0; i<NBranches; i++) IV_Branch.push_back(0);
    
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
                fData->SetBranchStatus((v_Branches[i] + "*").c_str(), 1);                   //Branch activation
                fData->SetBranchAddress((v_Branches[i] + ".").c_str(), &IV_Branch[i]);      //Connection between variable and Branches
            }
        }

        //---------------------------------------------------------------------------------------------------------------------------------------------
        //---------------------------------------------------------------LOOP ON EVENTS AND FILE WRITING-----------------------------------------------
        //---------------------------------------------------------------------------------------------------------------------------------------------

        cout<<"\r =====> "<<(int)((double)f/NumberOfElements*100. +1)<<" % "<<flush;

        for (Long64_t  t = 1; t < NTime; t++)   //loop over scenario time
        {
            fData->GetEntry(t);     //Update all branched object to the new CLASS time step j

            v_Time.push_back(double(TimeSecond)/double(cYear)); //The time (in year) at this time step
            v_Power.push_back(Power);







        }
        TreeScenario->Fill();
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
g++ -std=c++11 -o CLASS_R2R CLASS_ROOT2ROOT.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
*/
