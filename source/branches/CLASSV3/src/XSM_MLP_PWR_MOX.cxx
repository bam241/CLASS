
#include "XSModel.hxx"
#include "XSM_MLP_PWR_MOX.hxx"
#include "CLASSLogger.hxx"
#include "StringLine.hxx"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include <TGraph.h>
#include <TString.h>
#include "TFile.h"
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
XSM_MLP_PWR_MOX::XSM_MLP_PWR_MOX(CLASSLogger* Log,string TMVA_Weight_Directory,string InformationFile, bool IsTimeStep)
{

	SetLog(Log);

	fIsStepTime=IsTimeStep;
	fTMVAWeightFolder = TMVA_Weight_Directory;
	if(InformationFile=="")
		fMLPInformationFile = TMVA_Weight_Directory+"/Data_Base_Info.nfo";
	else
		fMLPInformationFile=fTMVAWeightFolder+InformationFile;

	GetMLPWeightFiles();
	GetDataBaseInformation();

	if(IsLog())
	{
		// Warning
		INFO << " A EvolutionData has been define : " << endl;
		INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
	}

}
//________________________________________________________________________
XSM_MLP_PWR_MOX::~XSM_MLP_PWR_MOX() 
{
	fMapOfTMVAVariableNames.clear();
}
//________________________________________________________________________
void XSM_MLP_PWR_MOX::GetDataBaseInformation()
{
	ifstream FILE(fMLPInformationFile.c_str());

	if(FILE.good())
	{
		while(!FILE.eof())
		{	
			string line;
			getline(FILE, line);
			size_t foundRType = line.find("ReactorType :");
			size_t foundFType = line.find("FuelType :");
			size_t foundHM    = line.find("Heavy Metal (t) :"); 
			size_t foundPower = line.find("Thermal Power (W) :");
			size_t foundTime  = line.find("Time (s) :");
			size_t foundZAI	  = line.find("Z A I Name (input MLP) :");
			int pos=0;
			if(foundRType != std::string::npos)
			{	StringLine::NextWord(line,pos,':'); 
				fDataBaseRType = atof( (StringLine::NextWord(line,pos,':')).c_str() );
			}
			if(foundFType != std::string::npos)
			{	StringLine::NextWord(line,pos,':'); 
				fDataBaseFType = atof( (StringLine::NextWord(line,pos,':')).c_str() );
			}
			if(foundHM != std::string::npos)
			{	StringLine::NextWord(line,pos,':'); 
				fDataBaseHMMass = atof( (StringLine::NextWord(line,pos,':')).c_str() );
			}
			if(foundPower !=std::string::npos)
			{	StringLine::NextWord(line,pos,':');
				fDataBasePower = atof( (StringLine::NextWord(line,pos,':') ).c_str() );
			}	
 			pos=0;
			if(foundTime!=std::string::npos)
			{
				StringLine::NextWord(line,pos,':');
				while( pos< (int)line.size() )
					fMLP_Time.push_back( atof( (StringLine::NextWord(line,pos,' ')).c_str() ));
			}	
				
			if(foundZAI != std::string::npos)
			{
				while(!FILE.eof())
				{
					int Z=-4;
					int A=-4;
					int I=-4;
					string Name;
					getline(FILE, line);
					stringstream ssline;
					ssline<<line;
					ssline>>Z>>A>>I>>Name;
				//	cout<<Z<<" "<<A<<" "<<I<<" "<<Name<<endl;
					fMapOfTMVAVariableNames.insert( pair<ZAI,string>(ZAI(Z,A,I),Name) );

				}
			}	

		}
	}
	else
	{
		ERROR << "Can't find/open file " << fMLPInformationFile << endl;
		exit(0);
	}

/********DEBUG************************************
	cout<<"Heavy Metal (t) :"<<fDataBaseHMMass<<endl;
	cout<<"Thermal Power (W) :"<<fDataBasePower<<endl;
	cout<<"Time (s) :"<<endl;
	for (int i = 0; i < (int)fMLP_Time.size(); ++i)
		cout<<fMLP_Time[i]<<endl;
	cout<<"Z A I Name (input MLP) :"<<endl;
		
	map<ZAI ,string >::iterator it;

	for (it= fMapOfTMVAVariableNames.begin();it!=fMapOfTMVAVariableNames.end();it++)
		cout<< it->first.Z()<<" "<<it->first.A()<<" "<<it->second<<endl;

*/
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
   	size_t foundext = sReaction.find(".weights.xml");
   	sReaction = sReaction.substr(0,foundext);

	if(sReaction=="fis")
		Reaction=0;
	if(sReaction=="cap")
		Reaction=1;
	if(sReaction=="n2n")
		Reaction=2;

	if(Z<=0 || A<=0 || I<0 || Reaction==-1)
	{
		ERROR << " wrong TMVA weight format " << endl;
		exit(0);
	}
}
//________________________________________________________________________
TTree* XSM_MLP_PWR_MOX::CreateTMVAInputTree(IsotopicVector isotopicvector,int TimeStep)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/
	TTree*   InputTree = new TTree("InTMP", "InTMP");

	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);

	float Time=0;

	IsotopicVector IVInputTMVA;
	map<ZAI ,string >::iterator it;
	int j=0;

	for( it = fMapOfTMVAVariableNames.begin()  ; it!=fMapOfTMVAVariableNames.end() ; it++)
	{	
		InputTree->Branch( ((*it).second).c_str()	,&InputTMVA[j], ((*it).second + "/F").c_str());
		IVInputTMVA+= ((*it).first)*1;
		j++;
	}
	
	if( !fIsStepTime)
		InputTree->Branch(	"Time"	,&Time	,"Time/F"	);

	IsotopicVector IVAccordingToUserInfoFile = isotopicvector.GetThisComposition(IVInputTMVA);

	double Ntot = IVAccordingToUserInfoFile.GetTotalMass()*1e6/IVAccordingToUserInfoFile.MeanMolar()*6.02214129e23;

	IVAccordingToUserInfoFile = IVAccordingToUserInfoFile/Ntot;

	j=0;
	map<ZAI ,string >::iterator it2;
	for( it2 = fMapOfTMVAVariableNames.begin() ; it2!=fMapOfTMVAVariableNames.end() ; it2++)
	{
		InputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it2).first ) ;
		j++;
	}

	Time=fMLP_Time[TimeStep];

	InputTree->Fill();
	return InputTree;
}
//________________________________________________________________________
double XSM_MLP_PWR_MOX::ExecuteTMVA(string WeightFile,TTree* InputTree)
{
	// --- Create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "V:Color" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);
	Float_t Time;

	map<ZAI ,string >::iterator it;
	int j=0;
	for( it = fMapOfTMVAVariableNames.begin()  ; it!=fMapOfTMVAVariableNames.end() ; it++)
	{	reader->AddVariable( ( (*it).second ).c_str(),&InputTMVA[j]);
		j++;
	}
	if(!fIsStepTime)	
		reader->AddVariable( "Time" ,&Time);

	// --- Book the MVA methods

	string dir    = fTMVAWeightFolder;
	if(dir[dir.size()-1]!='/')
   		dir+="/";

	// Book method MLP
	TString methodName = "MLP method";
	TString weightpath = dir + WeightFile ;
	reader->BookMVA( methodName, weightpath );

	map<ZAI ,string >::iterator it2;
	j=0;
	for( it2 = fMapOfTMVAVariableNames.begin()  ; it2!=fMapOfTMVAVariableNames.end() ; it2++)
	{
		InputTree->SetBranchAddress(( (*it2).second ).c_str(),&InputTMVA[j]);
		j++;
	}

	if(!fIsStepTime)	
		InputTree->SetBranchAddress( "Time" ,&Time  );

	InputTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;


	return (double)val;
}
//________________________________________________________________________
EvolutionData XSM_MLP_PWR_MOX::GetCrossSectionsTime(IsotopicVector IV)
{

	cout<<"=====Building Evolution Data From TMVA MLP====="<<endl;

	EvolutionData EvolutionDataFromMLP = EvolutionData();

	map<ZAI,TGraph*> ExtrapolatedXS[3];
	/*************DATA BASE INFO****************/
	EvolutionDataFromMLP.SetReactorType(fDataBaseRType);
	EvolutionDataFromMLP.SetFuelType(fDataBaseFType);
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

		for(int TimeStep=0;TimeStep<int(fMLP_Time.size());TimeStep++)
		{
			TTree* InputTree = CreateTMVAInputTree(IV,TimeStep);

			pair< map<ZAI, TGraph*>::iterator, bool> IResult;

			IResult = ExtrapolatedXS[Reaction].insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph() ) );

			double XSValue = ExecuteTMVA(fWeightFiles[i],InputTree );
			if(IResult.second )
			{
				(IResult.first)->second->SetPoint(0, (double)GetMLPTime()[TimeStep], XSValue );

			}
			else
			{
				(IResult.first)->second->SetPoint( (IResult.first)->second->GetN(), (double)GetMLPTime()[TimeStep], XSValue );
			}

			delete InputTree;
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


	cout<<"=====Evolution Data Built====="<<endl;
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

    if(Z==-1 || A==-1 || I==-1 || Reaction==-1 || TimeStep==-1)
	{
		ERROR << " wrong TMVA weight format " << endl;
		exit(0);
	}
}

//________________________________________________________________________
EvolutionData XSM_MLP_PWR_MOX::GetCrossSectionsStep(IsotopicVector IV)
{
	cout<<"Getting cross section step mode ..."<<endl;
	 TTree* InputTree=CreateTMVAInputTree(IV);
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
			(IResult.first)->second->SetPoint(0, (double)GetMLPTime()[TimeStep], ExecuteTMVA(fWeightFiles[i],InputTree) );
		}
		else
		{
			(IResult.first)->second->SetPoint( (IResult.first)->second->GetN(), (double)GetMLPTime()[TimeStep], ExecuteTMVA(fWeightFiles[i],InputTree) );
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

	cout<<"=====Evolution Data Built (step mode)====="<<endl;
	delete InputTree;
	return EvolutionDataFromMLP;
}
//________________________________________________________________________
EvolutionData XSM_MLP_PWR_MOX::GetCrossSections(IsotopicVector IV ,double t)
{	
	if(t!=0)
		WARNING << " argument t has non effect here " << endl;

	EvolutionData EV;
	if(fIsStepTime)
		EV=GetCrossSectionsStep(IV);

	else
		EV=GetCrossSectionsTime(IV);


return EV;
}
//________________________________________________________________________
