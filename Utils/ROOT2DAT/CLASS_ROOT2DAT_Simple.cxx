/* 

Simple code used to convert root output file produced by CLASS to simple data file

Remains : 
	- Take into account all CLASSES (SepPlant, etc.)
	- Calculate Natural Uranium consumption


Authors:

BaL
Nico. T.

*/
#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv)
{

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------VARIABLES---------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    vector <string> v_Branches; // vector that will contain all the branches stored in the TTree

    string s_tmp; 
	Long64_t Time = 0;
	double Power = 0;
	vector<IsotopicVector *> IV_Branch;
	int NStocks=0; int NPools=0; int NReactors=0; int NFabPlants=0;
	vector<ZAI > v_ZAI;


//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------FILES AND USAGE---------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

	// If bad execution command
	if(!argc%3)
	{
		cout<<endl<<"##################"<<endl<<endl;
		cout<<"USAGE:"<<endl<<endl;
		cout<<"CLASS_R2D_S FileName.root Z1 A1 I1 Z2 A2 I3 etc..."<<endl;
		cout<<endl<<"##################"<<endl<<endl;
		exit(1);
	}

	// If bad ROOT file
	string FileName = argv[1];
	TFile *TFileName = new TFile(FileName.c_str());
	if (TFileName == 0)
	{
		cout << "can't read/open file " << FileName << endl;
		exit(0);
	}

	// Get wanted ZAI
	int NI = (int)((argc-2)/3); int Z=0; int A=0; int I=0;
	for(int i=0; i<NI; i++)
	{
		Z = atoi(argv[2+3*i]); A = atoi(argv[3+3*i]); I = atoi(argv[4+3*i]);
		v_ZAI.push_back(ZAI(Z,A,I));
	}
	//for(int i=0; i<NI; i++) cout<<v_ZAI[i].Z()<<" "<<v_ZAI[i].A()<<" "<<v_ZAI[i].I()<<" "<<endl;

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------LOAD BRANCHES-----------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

	// ROOT file information and List of Branches
	TTree *fData = new TTree();	
	fData = (TTree*) gDirectory->Get(TFileName->GetListOfKeys()->At(TFileName->GetNkeys() - 1)->GetName());
	fData->SetBranchStatus("*", 0); // All branches are unbranched
    int NBranches = fData->GetListOfBranches()->GetEntries();
	
	cout<<endl<<endl<<"################################################"<<endl;
	cout<<"TTree " << fData->GetName()<<" Loaded "<<endl;
	cout<<"################################################"<<endl;
	cout<<"List of existing Branches : "<<endl<<endl;
    for(int i=0; i<NBranches; i++)
    {
    	s_tmp = fData->GetListOfBranches()->At(i)->GetName();
    	if(s_tmp[s_tmp.size()-1]=='.') s_tmp = s_tmp.substr(0, s_tmp.size()-1);
    	v_Branches.push_back(s_tmp);
    	cout<<v_Branches[i]<<endl;
    } 
	cout<<"################################################"<<endl;
	cout<<"################################################"<<endl<<endl;
	
	// Number of Stocks
	for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="S_") NStocks++;
	Storage* B_Stock[NStocks]; int IndiceStock=0;
	// Number of Pools
	for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="P_") NPools++;
	Pool* B_Pool[NPools]; int IndicePool=0;
	// Number of Reactors
	for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="R_") NReactors++;
	Reactor* B_Reactor[NReactors]; int IndiceReactor=0;
	// Number of Fabrication Plants
	for(int i=0; i<NBranches; i++) if (v_Branches[i].substr(0,2)=="F_") NFabPlants++;
	FabricationPlant* B_FabPlant[NFabPlants]; int IndiceFabPlant=0;

	//Time Steps
	Long64_t NTime = fData->GetEntries();

	// Output data ascii file name
    size_t SPos = FileName.find(".root");    FileName.replace(SPos, std::string(".root").length(), ".dat");
    ofstream DataFileName(FileName.c_str()); DataFileName<<scientific<<setprecision(5);

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------FILE MATRIX-------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

int NumberOfElements = NI; // Number of group to be printed (i.e. MA, U, FP, Unat, etc...)

DataFileName<<"C AbsTime		0"<<endl;
DataFileName<<"C -------------------------------------------------------------------------------------------"<<endl;

DataFileName<<"C "<<setw(20)<<" "<<setw(5);
for (int i=0; i<NI; i++) DataFileName<<setw(4)<<v_ZAI[i].Z()<<"_"<<v_ZAI[i].A()<<"_"<<v_ZAI[i].I();
DataFileName<<endl;

DataFileName<<"C -------------------------------------------------------------------------------------------"<<endl;
for(int i=2; i<NBranches; i++)
{
	DataFileName<<"C "<<setw(20)<<v_Branches[i]; 
	//DataFileName.seekp((i-1) * 100 + 20);
	for(int e=0; e<NumberOfElements; e++)
	{
		DataFileName<<setw(10)<<2 + NumberOfElements*(i-2) + e;
	}DataFileName<<endl;
}
DataFileName<<"C -------------------------------------------------------------------------------------------"<<endl;
DataFileName<<"C"<<endl;
DataFileName<<"C"<<endl;
DataFileName<<"C"<<endl;


//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------CONNECT BRANCHES--------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    // Time
	fData->SetBranchStatus("AbsTime", 1);		//Branch activation
	fData->SetBranchAddress("AbsTime", &Time);	//Connection between variable and Branches
	// Thermal Power
	fData->SetBranchStatus("ParcPower", 1);
	fData->SetBranchAddress("ParcPower", &Power);

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
			fData->SetBranchStatus((v_Branches[i] + "*").c_str(), 1);					//Branch activation
			fData->SetBranchAddress((v_Branches[i] + ".").c_str(), &IV_Branch[i]);		//Connection between variable and Branches
		}
    }

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------LOOP ON EVENTS AND FILE WRITING-----------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------

    cout<<endl<<"#########################"<<endl<<"Progression : "<<endl;
	for (Long64_t  t = 0; t < NTime; t++)	//loop over scenario time
	{
		cout<<"\r =====> "<<(int)((double)t/NTime*100. +1)<<" %"<<flush;

		IndiceStock = 0; IndicePool = 0; IndiceReactor = 0; IndiceFabPlant = 0;
		fData->GetEntry(t);		//Update all branched object to the new CLASS time step j

		// Time 
		double Time_year = double(Time)/double(cYear); //The time (in year) at this time step

		// Calculate UNat

		// File Writing
		DataFileName<<Time_year<<"\t"<<Power<<"\t";

	    for(int i=2; i<NBranches; i++)
	    {
			if (v_Branches[i].substr(0,2)=="S_")
			{
				for(int j=0; j<NI; j++) DataFileName<<B_Stock[IndiceStock]->GetInsideIV().GetIsotopicQuantity()[v_ZAI[j]]<<"\t";
				IndiceStock++;
			}
			else if (v_Branches[i].substr(0,2)=="P_")
			{
				for(int j=0; j<NI; j++) DataFileName<<B_Pool[IndicePool]->GetInsideIV().GetIsotopicQuantity()[v_ZAI[j]]<<"\t";
				IndicePool++;
			}
			else if (v_Branches[i].substr(0,2)=="R_")
			{
				for(int j=0; j<NI; j++) DataFileName<<B_Reactor[IndiceReactor]->GetInsideIV().GetIsotopicQuantity()[v_ZAI[j]]<<"\t";
				IndiceReactor++;
			}
			else if (v_Branches[i].substr(0,2)=="F_")
			{
				for(int j=0; j<NI; j++) DataFileName<<B_FabPlant[IndiceFabPlant]->GetInsideIV().GetIsotopicQuantity()[v_ZAI[j]]<<"\t";
				IndiceFabPlant++;
			}
			else
			{
				for(int j=0; j<NI; j++) DataFileName<<IV_Branch[i]->GetIsotopicQuantity()[v_ZAI[j]]<<"\t";
			}
	    }
		DataFileName<<endl;
	}

    cout<<endl<<"END OF ..."<<endl<<"#########################"<<endl;
	DataFileName.close();
	delete fData;
}

/*
g++ -o CLASS_R2D CLASS_ROOT2DAT.cxx -I$CLASS_PATH/source/include -I$CLASS_PATH/source/external -I$CLASS_PATH/source/Model/Equivalence -I$CLASS_PATH/source/Model/Irradiation -I$CLASS_PATH/source/Model/XS -L$CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs`
*/
