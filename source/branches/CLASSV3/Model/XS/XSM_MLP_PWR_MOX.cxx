
#include "XSModel.hxx"
#include "XSM_MLP_PWR_MOX.hxx"
#include "LogFile.hxx"
#include "StringLine.hxx"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include <TGraph.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include <dirent.h>
#include <errno.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

//________________________________________________________________________
//
//		XSM_MLP_PWR_MOX
//
//
//
//
//________________________________________________________________________
XSM_MLP_PWR_MOX::XSM_MLP_PWR_MOX(LogFile* Log,string TMVA_Weight_Directory,string InformationFile, bool IsTimeStep)
{

	SetLog(Log);

	fIsStepTime=IsTimeStep;
	fTMVAWeightFolder = TMVA_Weight_Directory;
	if(InformationFile=="")
		fMLPInformationFile=TMVA_Weight_Directory+"/Data_Base_Info.nfo";
	else
		fMLPInformationFile=InformationFile;

	GetMLPWeightFiles();
	SetDataBaseInformation();

	if(IsLog())
	{
		// Warning
		cout	<< "!!INFO!! !!!XSM_MLP_PWR_MOX!!! A EvolutionData has been define :" << endl;
		cout	<< "\t His TMVA folder is : \"" << fTMVAWeightFolder << "\"" << endl;

		
		GetLog()->fLog 	<< "!!INFO!! !!!XSM_MLP_PWR_MOX!!! A EvolutionData has been define :" << endl;
		GetLog()->fLog	<<"\t His TMVA folder is : \"" << fTMVAWeightFolder << "\"" << endl;
	}
}

//________________________________________________________________________
void XSM_MLP_PWR_MOX::SetDataBaseInformation()
{
	ifstream FILE(fMLPInformationFile.c_str());

	if(FILE.good())
	{
		double HM_Mass_tonne=0;
		double Power_watt=0;

		FILE>>HM_Mass_tonne;
		FILE>>Power_watt;

		fDataBaseHMMass=HM_Mass_tonne;
		fDataBasePower=Power_watt;

		while(FILE.eof())
		{
			double TIME=-1;
			FILE>>TIME;

			if(TIME!=-1)
				fMLP_Time.push_back(TIME);
		}
	}
	else
	{
		cout<<"Can't find/open file "<<fMLPInformationFile<<endl;
		exit(0);
	}

}
//________________________________________________________________________
void XSM_MLP_PWR_MOX::GetMLPWeightFiles()
{

	/**********Get All TMVA weight files*******************/

	//check if the folder containing weights exists
	DIR* rep = NULL;
	struct dirent* fichierLu = NULL;
	rep = opendir(fTMVAWeightFolder.c_str());
	if (rep == NULL)
	{
		perror("");
		exit(1);
	}

	/**Save file names of TMVA  weights*/
	fWeightFiles.resize(0);
	while ((fichierLu = readdir(rep)) != NULL)
	{
		string FileName= fichierLu->d_name ;
		if(FileName != "." && FileName != ".." )
		{
			if(FileName[FileName.size()-3]=='x'  &&  FileName[FileName.size()-2]=='m' && FileName[FileName.size()-1]=='l' )
				fWeightFiles.push_back(FileName);

		}
	}

}
//________________________________________________________________________
//________________________________________________________________________
//
//	Time  (MLP take time as parameter)
//________________________________________________________________________
//________________________________________________________________________
void XSM_MLP_PWR_MOX::ReadWeightFile(string Filename, int &Z, int &A, int &I, int &Reaction)
{
	Z=-1;
	A=-1;
	I=-1;
	Reaction=-1;

	size_t found = Filename.find("XS");

	string NameJOB;
	NameJOB=Filename.substr(found);
	int pos=0;

	StringLine::NextWord(NameJOB, pos, '_');

	Z = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );

	A = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );

	I = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );

   	string  sReaction = (StringLine::NextWord(NameJOB,pos,'_') ).c_str() ;
	if(sReaction=="fis")
		Reaction=0;
	if(sReaction=="cap")
		Reaction=1;
	if(sReaction=="n2n")
		Reaction=2;

}
//________________________________________________________________________
void XSM_MLP_PWR_MOX::CreateTMVAInputTree(IsotopicVector isotopicvector,int TimeStep)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/

	TFile*   InputFile = new TFile("./InputTMP.root","RECREATE");
	TTree*   InputTree = new TTree("InTMP", "InTMP");
	float teneur = 0;
	float Pu8    = 0;
	float Pu9    = 0;
	float Pu10   = 0;
	float Pu11   = 0;
	float Pu12   = 0;
	float Am1    = 0;
	float U5_enrichment 	 = 0;
	float Time   = 0;

	InputTree->Branch(	"teneur",&teneur,"teneur/F"	);
	InputTree->Branch(	"Pu8"	,&Pu8	,"Pu8/F"	);
	InputTree->Branch(	"Pu9"	,&Pu9	,"Pu9/F"	);
	InputTree->Branch(	"Pu10"	,&Pu10	,"Pu10/F"	);
	InputTree->Branch(	"Pu11"	,&Pu11	,"Pu11/F"	);
	InputTree->Branch(	"Pu12"	,&Pu12	,"Pu12/F"	);
	InputTree->Branch(	"Am1"	,&Am1	,"Am1/F"	);
	InputTree->Branch(	"U5_enrichment"	,&U5_enrichment	,"U5_enrichment/F"	);
	InputTree->Branch(	"Time"	,&Time	,"Time/F"	);



	float U8     = isotopicvector.GetZAIIsotopicQuantity(92,238,0);
	float U5     = isotopicvector.GetZAIIsotopicQuantity(92,235,0);
	float U4     = isotopicvector.GetZAIIsotopicQuantity(92,234,0);

	float UTOT = U8 + U5 + U4;

	Pu8    	   = isotopicvector.GetZAIIsotopicQuantity(94,238,0);
	Pu9    	   = isotopicvector.GetZAIIsotopicQuantity(94,239,0);
	Pu10   	   = isotopicvector.GetZAIIsotopicQuantity(94,240,0);
	Pu11   	   = isotopicvector.GetZAIIsotopicQuantity(94,241,0);
	Pu12   	   = isotopicvector.GetZAIIsotopicQuantity(94,242,0);
	Am1        = isotopicvector.GetZAIIsotopicQuantity(95,241,0);

	teneur = (Pu8+Pu9+Pu10+Pu11+Pu12+Am1)/(Pu8+Pu9+Pu10+Pu11+Pu12+Am1+U8+U5+U4); //prop mol. Pu

	double TOTPU=(Pu8+Pu9+Pu10+Pu11+Pu12+Am1);

	Pu8 = Pu8  / TOTPU;
	Pu9 = Pu9  / TOTPU;
	Pu10= Pu10 / TOTPU;
	Pu11= Pu11 / TOTPU;
	Pu12= Pu12 / TOTPU;
	Am1 = Am1  / TOTPU;

	U5_enrichment = U5 / UTOT;

	Time=fMLP_Time[TimeStep];

	if(Pu8 + Pu9 + Pu10 + Pu11 + Pu12 + Am1 > 1.00001 )//?????1.00001??? I don't know it! goes in condition if =1 !! may be float/double issue ...
	{
		cout<<"!!!!!!!!!!!ERRORR!!!!!!!!!!!!"<<endl;
		cout<<Pu8<<" "<<Pu9<<" "<<Pu10<<" "<<Pu11<<" "<<Pu12<<" "<<Am1<<endl;
		exit(0);
	}
	// All value are molar (!weight)

	InputTree->Fill();

	InputFile->Write();
	delete InputTree;
	InputFile-> Close();
	delete InputFile;
}
///________________________________________________________________________
double XSM_MLP_PWR_MOX::ExecuteTMVA(string WeightFile)
{
	// --- Create the Reader object

	TMVA::Reader *reader = new TMVA::Reader( "Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	Float_t Pu8,Pu9,Pu10,Pu11,Pu12,Am1,Time,teneur;
	reader->AddVariable("teneur",&teneur);
	reader->AddVariable( "Pu8"  ,&Pu8 );
	reader->AddVariable( "Pu9"  ,&Pu9 );
	reader->AddVariable( "Pu10" ,&Pu10);
	reader->AddVariable( "Pu11" ,&Pu11);
	reader->AddVariable( "Pu12" ,&Pu12);
	reader->AddVariable( "Am1"  ,&Am1 );
	reader->AddVariable( "Time" ,&Time);

	// --- Book the MVA methods

	string dir    = fTMVAWeightFolder;
	if(dir[dir.size()-1]!='/')
   		dir+="/";

	// Book method MLP
	TString methodName = "MLP method";
	TString weightpath = dir + WeightFile ;
	reader->BookMVA( methodName, weightpath );

	// Prepare input tree
	TFile *input(0);
	if (!gSystem->AccessPathName( "./InputTMP.root" )) {
		input = TFile::Open( "./InputTMP.root" ); // check if file in local directory exists
	}
	if (!input) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}

	TTree* theTree = (TTree*)input->Get("InTMP");

	theTree->SetBranchAddress("teneur",&teneur);
	theTree->SetBranchAddress( "Pu8"  ,&Pu8   );
	theTree->SetBranchAddress( "Pu9"  ,&Pu9   );
	theTree->SetBranchAddress( "Pu10" ,&Pu10  );
	theTree->SetBranchAddress( "Pu11" ,&Pu11  );
	theTree->SetBranchAddress( "Pu12" ,&Pu12  );
	theTree->SetBranchAddress( "Am1"  ,&Am1   );
	theTree->SetBranchAddress( "Time" ,&Time  );

	theTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;
	delete theTree;
	input->Close();
	delete input;
	//cout<<"....done"<<endl;

	//cout<<val<<endl;
	return (double)val;
}
//________________________________________________________________________
EvolutionData XSM_MLP_PWR_MOX::GetCrossSectionsTime(IsotopicVector IV)
{

	//cout<<"=====Building Evolution Data From TMVA MLP====="<<endl;

	EvolutionData EvolutionDataFromMLP = EvolutionData();

	map<ZAI,TGraph*> ExtrapolatedXS[3];
	/*************DATA BASE INFO****************/
	EvolutionDataFromMLP.SetReactorType("PWR");
	EvolutionDataFromMLP.SetFuelType("MOX");
	EvolutionDataFromMLP.SetPower(fDataBasePower);
	EvolutionDataFromMLP.SetHeavyMetalMass(fDataBaseHMMass);
	/************* The Cross sections***********/
	for(int i=0;i<int(fWeightFiles.size());i++)
	{
		int Z=-2;
		int A=-2;
		int I=-2;
		int Reaction=-2;
		ReadWeightFile( fWeightFiles[i], Z, A, I, Reaction);

		ZAI zaitmp = ZAI(Z,A,I);

		for(int TimeStep=0;TimeStep<int(fMLP_Time.size());TimeStep++)
		{
			std::ifstream ifs ("InputTMP.root");
	 		if (ifs.is_open())
	 		{  ifs.close();
	 			system( "rm InputTMP.root" );
	 		}
			CreateTMVAInputTree(IV,TimeStep);

			pair< map<ZAI, TGraph*>::iterator, bool> IResult;

			IResult = ExtrapolatedXS[Reaction].insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph() ) );


			if(IResult.second )
			{
				(IResult.first)->second->SetPoint(0, (double)GetMLPTime()[TimeStep], ExecuteTMVA(fWeightFiles[i]) );
			}
			else
			{
				(IResult.first)->second->SetPoint( (IResult.first)->second->GetN(), (double)GetMLPTime()[TimeStep], ExecuteTMVA(fWeightFiles[i]) );
			}

  			system( "rm InputTMP.root" );
  		}

	}

	/**********Sorting TGraph*********/
	for(int x=0;x<3;x++)
	{	map<ZAI,TGraph*>::iterator it;
		for(it = ExtrapolatedXS[x].begin(); it != ExtrapolatedXS[x].end(); it++)
			it->second->Sort();
	}
	/**********Filling Matrices*/
	EvolutionDataFromMLP.SetFissionXS(ExtrapolatedXS[0]);
	EvolutionDataFromMLP.SetCaptureXS(ExtrapolatedXS[1]);
	EvolutionDataFromMLP.Setn2nXS(ExtrapolatedXS[2]);

	//cout<<"=====Evolution Data Built====="<<endl;
	return EvolutionDataFromMLP;
}
//________________________________________________________________________
//________________________________________________________________________
//
//	Time step (1 MLP per time step)
//________________________________________________________________________
//________________________________________________________________________
void XSM_MLP_PWR_MOX::ReadWeightFileStep(string Filename, int &Z, int &A, int &I, int &Reaction, int &TimeStep)
{
	Z=-1;
	A=-1;
	I=-1;
	Reaction=-1;

	size_t found = Filename.find("XS");

	string NameJOB;
	NameJOB=Filename.substr(found);
	int pos=0;

	StringLine::NextWord(NameJOB, pos, '_');

	Z = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );

	A = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );

	I = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );

   	string  sReaction = (StringLine::NextWord(NameJOB,pos,'_') ).c_str() ;
	if(sReaction=="fis")
		Reaction=0;
	if(sReaction=="cap")
		Reaction=1;
	if(sReaction=="n2n")
		Reaction=2;

    TimeStep = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );

}
//________________________________________________________________________
void XSM_MLP_PWR_MOX::CreateTMVAInputTreeStep(IsotopicVector isotopicvector)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/

	TFile*   InputFile = new TFile("./InputTMP.root","RECREATE");
	TTree*   InputTree = new TTree("InTMP", "InTMP");
	float teneur = 0;
	float Pu8    = 0;
	float Pu9    = 0;
	float Pu10   = 0;
	float Pu11   = 0;
	float Pu12   = 0;
	float Am1    = 0;
	float U5_enrichment 	 = 0;

	InputTree->Branch(	"teneur",&teneur,"teneur/F"	);
	InputTree->Branch(	"Pu8"	,&Pu8	,"Pu8/F"	);
	InputTree->Branch(	"Pu9"	,&Pu9	,"Pu9/F"	);
	InputTree->Branch(	"Pu10"	,&Pu10	,"Pu10/F"	);
	InputTree->Branch(	"Pu11"	,&Pu11	,"Pu11/F"	);
	InputTree->Branch(	"Pu12"	,&Pu12	,"Pu12/F"	);
	InputTree->Branch(	"Am1"	,&Am1	,"Am1/F"	);
	InputTree->Branch(	"U5_enrichment"	,&U5_enrichment	,"U5_enrichment/F"	);


	float U8     = isotopicvector.GetZAIIsotopicQuantity(92,238,0);
	float U5     = isotopicvector.GetZAIIsotopicQuantity(92,235,0);
	float U4     = isotopicvector.GetZAIIsotopicQuantity(92,234,0);

	float UTOT = U8 + U5 + U4;

	Pu8    	   = isotopicvector.GetZAIIsotopicQuantity(94,238,0);
	Pu9    	   = isotopicvector.GetZAIIsotopicQuantity(94,239,0);
	Pu10   	   = isotopicvector.GetZAIIsotopicQuantity(94,240,0);
	Pu11   	   = isotopicvector.GetZAIIsotopicQuantity(94,241,0);
	Pu12   	   = isotopicvector.GetZAIIsotopicQuantity(94,242,0);
	Am1        = isotopicvector.GetZAIIsotopicQuantity(95,241,0);

	teneur = (Pu8+Pu9+Pu10+Pu11+Pu12+Am1)/(Pu8+Pu9+Pu10+Pu11+Pu12+Am1+U8+U5+U4); //prop mol. Pu

	double TOTPU=(Pu8+Pu9+Pu10+Pu11+Pu12+Am1);

	Pu8 = Pu8  / TOTPU;
	Pu9 = Pu9  / TOTPU;
	Pu10= Pu10 / TOTPU;
	Pu11= Pu11 / TOTPU;
	Pu12= Pu12 / TOTPU;
	Am1 = Am1  / TOTPU;

	U5_enrichment = U5 / UTOT;

	if(Pu8 + Pu9 + Pu10 + Pu11 + Pu12 + Am1 > 1.00001 )//?????1.00001??? I don't know it! goes in condition if =1 !! may be float/double issue ...
	{
		cout<<"!!!!!!!!!!!ERRORR!!!!!!!!!!!!"<<endl;
		cout<<Pu8<<" "<<Pu9<<" "<<Pu10<<" "<<Pu11<<" "<<Pu12<<" "<<Am1<<endl;
		exit(0);
	}
	// All value are molar (!weight)

	InputTree->Fill();

	InputFile->Write();
	delete InputTree;
	InputFile-> Close();
	delete InputFile;
}
///________________________________________________________________________
double XSM_MLP_PWR_MOX::ExecuteTMVAStep(string WeightFile)
{
	// --- Create the Reader object

	TMVA::Reader *reader = new TMVA::Reader( "Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	Float_t Pu8,Pu9,Pu10,Pu11,Pu12,Am1,teneur;
	reader->AddVariable("teneur",&teneur);
	reader->AddVariable( "Pu8"  ,&Pu8 );
	reader->AddVariable( "Pu9"  ,&Pu9 );
	reader->AddVariable( "Pu10" ,&Pu10);
	reader->AddVariable( "Pu11" ,&Pu11);
	reader->AddVariable( "Pu12" ,&Pu12);
	reader->AddVariable( "Am1"  ,&Am1 );
	// --- Book the MVA methods

	string dir    = fTMVAWeightFolder;
	if(dir[dir.size()-1]!='/')
   		dir+="/";

	// Book method MLP
	TString methodName = "MLP method";
	TString weightpath = dir + WeightFile ;
	reader->BookMVA( methodName, weightpath );

	// Prepare input tree
	TFile *input(0);
	if (!gSystem->AccessPathName( "./InputTMP.root" )) {
		input = TFile::Open( "./InputTMP.root" ); // check if file in local directory exists
	}
	if (!input) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}

	TTree* theTree = (TTree*)input->Get("InTMP");

	theTree->SetBranchAddress("teneur",&teneur);
	theTree->SetBranchAddress( "Pu8"  ,&Pu8   );
	theTree->SetBranchAddress( "Pu9"  ,&Pu9   );
	theTree->SetBranchAddress( "Pu10" ,&Pu10  );
	theTree->SetBranchAddress( "Pu11" ,&Pu11  );
	theTree->SetBranchAddress( "Pu12" ,&Pu12  );
	theTree->SetBranchAddress( "Am1"  ,&Am1   );

	theTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;
	delete theTree;
	input->Close();
	delete input;
	//cout<<"....done"<<endl;

	//cout<<val<<endl;
	return (double)val;
}
//________________________________________________________________________
EvolutionData XSM_MLP_PWR_MOX::GetCrossSectionsStep(IsotopicVector IV)
{

	std::ifstream ifs ("InputTMP.root");
	 if (ifs.is_open())
	 {  ifs.close();
	 	system( "rm InputTMP.root" );
	 }
	CreateTMVAInputTreeStep(IV);
	//cout<<"=====Building Evolution Data From TMVA MLP====="<<endl;

	EvolutionData EvolutionDataFromMLP = EvolutionData();

	map<ZAI,TGraph*> ExtrapolatedXS[3];
	/*************DATA BASE INFO****************/
	EvolutionDataFromMLP.SetReactorType("PWR");
	EvolutionDataFromMLP.SetFuelType("MOX");
	EvolutionDataFromMLP.SetPower(fDataBasePower);
	EvolutionDataFromMLP.SetHeavyMetalMass(fDataBaseHMMass);

	/************* The Cross sections***********/

	for(int i=0;i<int(fWeightFiles.size());i++)
	{
		int Z=-2;
		int A=-2;
		int I=-2;
		int Reaction=-2;
		int TimeStep=-2;
		ReadWeightFileStep( fWeightFiles[i], Z, A, I, Reaction, TimeStep);

		ZAI zaitmp = ZAI(Z,A,I);

		pair< map<ZAI, TGraph*>::iterator, bool> IResult;

		IResult = ExtrapolatedXS[Reaction].insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph() ) );

		if( IResult.second )
		{
			(IResult.first)->second->SetPoint(0, (double)GetMLPTime()[TimeStep], ExecuteTMVAStep(fWeightFiles[i]) );
		}
		else
		{
			(IResult.first)->second->SetPoint( (IResult.first)->second->GetN(), (double)GetMLPTime()[TimeStep], ExecuteTMVAStep(fWeightFiles[i]) );
		}

	}

	system( "rm InputTMP.root" );
	/**********Sorting TGraph*********/
	for(int x=0;x<3;x++)
	{	map<ZAI,TGraph*>::iterator it;
		for(it = ExtrapolatedXS[x].begin(); it != ExtrapolatedXS[x].end(); it++)
			it->second->Sort();
	}
	/**********Filling Matrices*/
	EvolutionDataFromMLP.SetFissionXS(ExtrapolatedXS[0]);
	EvolutionDataFromMLP.SetCaptureXS(ExtrapolatedXS[1]);
	EvolutionDataFromMLP.Setn2nXS(ExtrapolatedXS[2]);

	//cout<<"=====Evolution Data Built====="<<endl;
	return EvolutionDataFromMLP;
}
//________________________________________________________________________
EvolutionData XSM_MLP_PWR_MOX::GetCrossSections(IsotopicVector IV)
{
	EvolutionData EV;
	if(fIsStepTime)
		EV=GetCrossSectionsStep(IV);

	else
		EV=GetCrossSectionsTime(IV);

return EV;
}
//________________________________________________________________________
