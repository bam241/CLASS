/*************************************************
             DESCRIPTION

 This example explained how to get human readable
outputs from CLASS .root output file
Please read carefully the comments.

This example is able to read ouput of 
ExampleParc_MLP_MOX.cxx
author : BLG
/*************************************************/
#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <vector>
#include "StringLine.hxx"


using namespace std;

TTree *fData;
//_____________________________________________________________________		
void OpenRootFile(TString filename)
{	
	fData = new TTree();
	cout<<"Opening "<<filename<<" ...."<<endl;
	TFile *FileIn;
	FileIn = TFile::Open(filename);
	if (FileIn == 0)
	{
    	// if we cannot open the file, print an error message and return immediatly
    	cout <<"can'tread/open file "<<filename<<endl;;
    	exit(0);
   	}
	fData = (TTree*)gDirectory->Get(FileIn->GetListOfKeys()->At(FileIn->GetNkeys()-1)->GetName() ) ;
	//cout<<"TTree Loaded "<<endl;
}
//________________________________________________________________________

int main(int argc, char** argv)
{

	if(argc != 2)
	{
		cout<<" Usage :"<<endl;
		cout<<"ReadRoot ExampleParc_MLP.root"<<endl;
		cout<<" open this .cxx file and read comments for help and learning"<<endl;
		cout<<"This example is only able to read ouput of ExampleParc_MLP_MOX.cxx"<<endl;
		exit(1);
	}
	/**********************************/
	/**	Opening Root file *************/
	/**********************************/

	string FileName = argv[1]; //the name of the CLASS output file

	OpenRootFile(FileName);

	ofstream ouput("ReadRootOuput.dat"); // we save extracted values of interest in this file
	/**********************************/
	/**	BRANCHING INPUT ***************/
	/**********************************/
	/*
		In this section we connect the CLASS objects present in the
		CLASS output file (.root)  to new objects of same type. 
		Doing that we can interact with the CLASS output using the
		methods of the CLASS libraries. To connect object you need to 
		know their name you defined in your .cxx .

		********************************
		NAMING BRANCHES
		*********************************
					NAMING FACILITIES
				*********************************
		Branchname is related with the name of a CLASS object you defined in your CLASS input
		Let assume you define a storage in the CLASS input (.cxx) like :

		Storage *Stock = new Storage(); //Definition of the stock
		Stock->SetName("Stock_UOX"); //Its name
		gCLASS->Add(Stock); //Adding the stock to the Scenario 

		In this case the Branchname variable will be
		S_Stock_UOX (note the "S_" added before the name I gave to this Storage)

		The name format to branch object is :

		For Reactor 			: R_NameYouGave
		For Pool    			: P_NameYouGave
		For Storage 			: S_NameYouGave
		For FabricationPlant	: F_NameYouGave
		For SeparationPlant		: C_NameYouGave

		If you didn't give a name to your object the default names are the following :

		For Reactor 			: R_ReactorID
		For Pool    			: P_PoolID
		For Storage 			: S_StorageID
		For FabricationPlant	: F_FabricationPlantID
		For SeparationPlant		: C_SepPlantID

		ID is an integer. 0 is the first OBJECT of the kind you add to the scenario with  
		gCLASS->Add(OBJECT);

					NAMING OTHER BRANCHES
				*********************************
		NAME        = AbsTime
		EXPLANATION = Scenario time (is in seconds) : 
		
		NAME        = OUTINCOME
		EXPLANATION = IsotopicVector containing all the materials comming from 
		outside of the scenario  

		NAME        = INCYCLE
		EXPLANATION = IsotopicVector containing all the materials in the cycle
		(i.e) ALL except Reactors, Wastes, and outincome 

		NAME        = WASTE	
		EXPLANATION = IsotopicVector containing all the materials in the waste
		
		NAME        =  STOCK;
		EXPLANATION = IsotopicVector containing all the materials contained in all
				  stocks
		
		NAME        = COOLING
		EXPLANATION = IsotopicVector containing all the materials contained in all
				  Pool

		NAME        = FUELFABRICATION
		EXPLANATION = IsotopicVector containing all the materials contained in all
				  FabricationPlants
		
		NAME        = REACTOR
		EXPLANATION = IsotopicVector containing all the materials contained in all
				  Reactors

		
	*/
		// All the branch are "unbranched" ;) (to be sure)
		// fData is the loaded TTRee contained in the CLASS output (.root) 
		fData->SetBranchStatus("*", 0); 
	//Branching the scenario time
		Long64_t Time = 0; //We create an integer to be connected to the scenario time
		fData->SetBranchStatus("AbsTime", 1);//We activate this branch
		fData->SetBranchAddress("AbsTime", &Time);//we "connect" our local variable Time to the one in the root file

	//Storages 

		Storage* Stock_UOX = new Storage();//We create a Storage to be connected to the one in the .root file 
		string Branchname = "S_StockUOX"; 
		string ActiveBranchName = Branchname + "*"; 
		fData->SetBranchStatus(ActiveBranchName.c_str(), 1);//activating the branch
		fData->SetBranchAddress((Branchname+".").c_str(), &Stock_UOX);//we "connect" our local object Stock_UOX to the one in the root file

	
	/* You can declare and connect as many Storage as you defined in your CLASS .cxx
		
		Storage* Stock_MOX = new Storage();
		 Branchname = "S_StockMOX";
		 ActiveBranchName = Branchname + "*"; 
		fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
		fData->SetBranchAddress((Branchname+".").c_str(), &Stock_MOX);
	/*

	/* Adding a Pool 

		Pool* Pool_UOX = new Pool();
		 Branchname = "P_Pool_UOX";
		 ActiveBranchName = Branchname + "*"; 
		fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
		fData->SetBranchAddress((Branchname+".").c_str(), &Pool_UOX);
		
	*/

	/* Adding a Reactor*/

		Reactor* reactor = new Reactor();
		Branchname = "R_Dampierre_MOX";
	    ActiveBranchName = Branchname + "*"; 
		fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
		fData->SetBranchAddress((Branchname+".").c_str(), &reactor);
		


	//If you want to connect many Reactors with default names
	/*	int NumOfReactor = 10;

		Reactor* reactors[41]; //declaring all the Reactors

		for(int r = 0 ; r<NumOfReactor ; ++)
		{	reactors[i] = 0; //Initialize reactors
			
			stringstream ssss;
			ssss<<"R_FBR_"<<r;
		    Branchname = ssss.str();
			ActiveBranchName = Branchname + "*"; 
			fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
			fData->SetBranchAddress((Branchname+".").c_str(), &reactors[r]);
		}
	*/

	/*  Adding a Fabrication Plant

		FabricationPlant *FP_MOX = new FabricationPlant();
		Branchname = "F_Fab_MOX";
		ActiveBranchName = Branchname + "*"; 
		fData->SetBranchStatus(ActiveBranchName.c_str(), 1);
		fData->SetBranchAddress((Branchname+".").c_str(), &FP_MOX);
	*/

	//BRANCHING ALL
		IsotopicVector* ALL_Total=0;
		Branchname = "TOTAL";
		fData->SetBranchStatus((Branchname+"*").c_str(), 1);
		fData->SetBranchAddress((Branchname+".").c_str(), &ALL_Total);


	//Isotopic vector to extract wanted nuclei
		IsotopicVector PlutoniumVector;
		PlutoniumVector.Add(94,238,0,1);
		PlutoniumVector.Add(94,239,0,1);
		PlutoniumVector.Add(94,240,0,1);
		PlutoniumVector.Add(94,241,0,1);
		PlutoniumVector.Add(94,242,0,1);
		PlutoniumVector.Add(95,241,0,1);

		IsotopicVector UraniumVector;
		UraniumVector.Add(92,238,0,1);
		UraniumVector.Add(92,235,0,1);

		/*	WANTED OUTPUTS
			We choose here to extract this quantities over time :
			The plutonium content loaded in the MOX reactor	
			The total amount of Uranium neptunium plutonium americium and curium in the scenario

			Fill free to look at the $CLASS_PATH/documentation/doxygen to find methods to access the value you want
		*/

		ouput<< "TIME(y)  PuContent (%w)  All_U(t)    All_Np(t)    All_Pu(t)    All_Am(t)    All_Cm(t)     "<<endl;

		Long64_t nentries = fData->GetEntries(); //Number of CLASS time step

		for (Long64_t  j = 0; j < nentries; j++)	//loop over scenario time
		{

			fData->GetEntry(j);//Update all branched object to the new CLASS time step j
			//Temps 
			double Time_year = double(Time)/double(cYear); //The time (in year) at this time step

			//Plutonium entering in MOX

			IsotopicVector Pu_At_Loading = reactor->GetIVBeginCycle().GetThisComposition(PlutoniumVector);//GetPlutonium vector
			IsotopicVector U_At_Loading  = reactor->GetIVBeginCycle().GetThisComposition(UraniumVector);//GetUranium vector
//
			double U_Totmass = U_At_Loading.GetTotalMass();
			double Pu_Totmass = Pu_At_Loading.GetTotalMass();
			double	WeighPuContent  = Pu_Totmass/(Pu_Totmass+U_Totmass)*100 ; //The plutonium content loaded in reactor

			//	//Total inventories

			double 	ALL_Pu		= ALL_Total->GetSpeciesComposition(94).GetTotalMass();
			double 	ALL_U		= ALL_Total->GetSpeciesComposition(92).GetTotalMass();
			double 	ALL_Am		= ALL_Total->GetSpeciesComposition(95).GetTotalMass();
			double 	ALL_Cm		= ALL_Total->GetSpeciesComposition(96).GetTotalMass();
			double 	ALL_Np		= ALL_Total->GetSpeciesComposition(93).GetTotalMass();
			
			ouput<<setprecision(5)<<Time_year <<" "<<setprecision(5)<<WeighPuContent<<" "<<setprecision(5)<<ALL_U<<" "<<setprecision(5)<<ALL_Np<<" "<<setprecision(5)<<ALL_Pu<<" "<<setprecision(5)<<ALL_Am<<" "<<setprecision(5)<<ALL_Cm<<" "<<endl;

		}

		delete fData;
	
cout << "Results have been written in ReadRootOuput.dat"<<endl;

}


//==========================================================================================
// Compilation
//==========================================================================================
/*
 
 g++ -o ReadOutput ReadOutput.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
 
 
 */
