
#include "XSModel.hxx"
#include "XSM_MLP.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "StringLine.hxx"

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
//		XSM_MLP
//
//
//
//
//________________________________________________________________________
XSM_MLP::XSM_MLP(string TMVA_Weight_Directory,string InformationFile, bool IsTimeStep):XSModel(new CLASSLogger("XSM_MLP.log"))
{
	
	fIsStepTime = IsTimeStep;
	fTMVAWeightFolder = TMVA_Weight_Directory;
	
	fInformationFile = fTMVAWeightFolder+InformationFile;
	
	GetMLPWeightFiles();
	
	INFO << "__A cross section interpolator using" << endl;
	INFO << "Multi Layer Perceptron has been define__" << endl;
	INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
	
	LoadKeyword();
	ReadNFO();

    InitialiseTMVAReader();

}

//________________________________________________________________________
XSM_MLP::XSM_MLP(CLASSLogger* Log,string TMVA_Weight_Directory,string InformationFile, bool IsTimeStep):XSModel(Log)
{
	
	fIsStepTime = IsTimeStep;
	fTMVAWeightFolder = TMVA_Weight_Directory;

	fInformationFile = fTMVAWeightFolder+InformationFile;
	
	GetMLPWeightFiles();
	
	INFO << "__A cross section interpolator using" << endl;
	INFO << "Multi Layer Perceptron has been define__" << endl;
	INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
	
	LoadKeyword();
	ReadNFO();

    InitialiseTMVAReader();

}

//________________________________________________________________________
XSM_MLP::~XSM_MLP()
{
	DBGL
	fMapOfTMVAVariableNames.clear();
	fDKeyword.clear();
	DBGL
}


//________________________________________________________________________
void XSM_MLP::LoadKeyword()
{
	DBGL
	fDKeyword.insert( pair<string, XS_MLP_DMthPtr>( "k_timestep",	& XSM_MLP::ReadTimeSteps));
	fDKeyword.insert( pair<string, XS_MLP_DMthPtr>( "k_zainame",	& XSM_MLP::ReadZAIName)	 );
	DBGL
}


//________________________________________________________________________
void XSM_MLP::ReadZAIName(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_zainame" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_zainame\" not found !" << endl;
		exit(1);
	}
	
	int Z = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	
	string name = StringLine::NextWord(line, pos, ' ');
	
	fMapOfTMVAVariableNames.insert( pair<ZAI,string>( ZAI(Z, A, I), name ) );
	
	DBGL
}

//________________________________________________________________________
void XSM_MLP::ReadTimeSteps(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_timestep" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_timestep\" not found !" << endl;
		exit(1);
	}
	
	while( pos < (int)line.size() )
		fMLP_Time.push_back( atof( (StringLine::NextWord(line,pos,' ')).c_str() ));
	DBGL
}


//________________________________________________________________________
void XSM_MLP::ReadLine(string line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	
	map<string, XS_MLP_DMthPtr>::iterator it = fDKeyword.find(keyword);
	
	if(it != fDKeyword.end())
		(this->*(it->second))( line );
	
	DBGL
}


 
//________________________________________________________________________
void XSM_MLP::GetMLPWeightFiles()
{
	DBGL
	/**********Get All TMVA weight files*******************/
	//check if the folder containing weights exists
	DIR* rep = NULL;
	struct dirent* fichierLu = NULL;
	rep = opendir(fTMVAWeightFolder.c_str());
	
	if (rep ==  NULL)
	{
		ERROR << " Reading error for TMVA weight folder  " << fTMVAWeightFolder.c_str() << " : " << strerror(errno) << endl;
		exit(1);
	}
	
	/**Save file names of TMVA  weights*/
	fWeightFiles.resize(0);
	while ((fichierLu = readdir(rep)) != NULL)
	{
		string FileName = fichierLu->d_name ;
		if(FileName != "." && FileName != ".." )
		{
			if(FileName[FileName.size()-3] == 'x'  &&  FileName[FileName.size()-2] == 'm' && FileName[FileName.size()-1] == 'l' && FileName[0] != '.' )
				fWeightFiles.push_back(FileName);
			
		}
	}
	DBGL
}
//________________________________________________________________________
//________________________________________________________________________
//
//	Time  (MLP take time as parameter)
//________________________________________________________________________
//________________________________________________________________________
void XSM_MLP::ReadWeightFile(string Filename, int &Z, int &A, int &I, int &Reaction)
{
	Z = -1;
	A = -1;
	I = -1;
	Reaction = -1;
	
	size_t found = Filename.find("XS");
	
	string NameJOB;
	NameJOB = Filename.substr(found);
	int pos = 0;
	
	StringLine::NextWord(NameJOB, pos, '_');
	
	Z = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	
	A = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	
	I = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	
	string  sReaction = (StringLine::NextWord(NameJOB,pos,'_') ).c_str() ;
	size_t foundext = sReaction.find(".weights.xml");
	sReaction = sReaction.substr(0,foundext);
	
	if(sReaction == "fis")
		Reaction = 0;
	if(sReaction == "cap")
		Reaction = 1;
	if(sReaction == "n2n")
		Reaction = 2;
	
	if(Z<= 0 || A<= 0 || I<0 || Reaction == -1)
	{
		ERROR << " wrong TMVA weight format " << endl;
		exit(0);
	}
	
}
//________________________________________________________________________
void XSM_MLP::UpdateInputComposition(IsotopicVector TheFreshfuel,int TimeStep)
{
	
    IsotopicVector IVAccordingToUserInfoFile = TheFreshfuel.GetThisComposition(IVInputTMVA);

    double Ntot = IVAccordingToUserInfoFile.GetSumOfAll();

    IVAccordingToUserInfoFile *= 1/Ntot;
	

    DBGV("INPUT TMVA");

    int j = 0;
	map<ZAI ,string >::iterator it2;

	for( it2 = fMapOfTMVAVariableNames.begin() ; it2 != fMapOfTMVAVariableNames.end() ; it2++)
	{
		InputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it2).first );
		DBGV((*it2).first.Z() << " " << (*it2).first.A() << " " << InputTMVA[j]);
		j++;
	}
	
	Time = fMLP_Time[TimeStep];

}

void XSM_MLP::InitialiseTMVAReader()
{
    reader = new TMVA::Reader( "Silent" );

    for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
        InputTMVA.push_back(0);


    for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
        InputTMVA.push_back(0);

    map<ZAI ,string >::iterator it;
    int j = 0;
    for( it = fMapOfTMVAVariableNames.begin()  ; it != fMapOfTMVAVariableNames.end() ; it++)
        {	reader->AddVariable( ( (*it).second ).c_str(),&InputTMVA[j]);
            j++;
        }
    if(!fIsStepTime)
        reader->AddVariable( "Time" ,&Time);
}

//________________________________________________________________________
double XSM_MLP::ExecuteTMVA(string WeightFile, IsotopicVector isotopicvector, int TimeStep)
{

    UpdateInputComposition(isotopicvector, TimeStep);


    DBGV( "File :" << WeightFile);

	// --- Book the MVA methods
	
	string dir    = fTMVAWeightFolder;
	if(dir[dir.size()-1] != '/')
		dir+= "/";
	
	// Book method MLP
	TString methodName = "MLP method";
	TString weightpath = dir + WeightFile ;
	reader->BookMVA( methodName, weightpath );

	float val = (reader->EvaluateRegression( methodName ))[0];
	
	DBGL
	
	return (double)val;
}
//________________________________________________________________________
EvolutionData XSM_MLP::GetCrossSectionsTime(IsotopicVector IV)
{
	DBGL
	
	EvolutionData EvolutionDataFromMLP = EvolutionData();
	
	map<ZAI,TGraph*> ExtrapolatedXS[3];
	/*************DATA BASE INFO****************/
	EvolutionDataFromMLP.SetReactorType(fDBRType);
	EvolutionDataFromMLP.SetFuelType(fDBFType);
	EvolutionDataFromMLP.SetPower(fDBPower);
	EvolutionDataFromMLP.SetHeavyMetalMass(fDBHMMass);
	/************* The Cross sections***********/
	for(int i = 0;i<int(fWeightFiles.size());i++)
	{
		int Z = -2;
		int A = -2;
		int I = -2;
		int Reaction = -2;
		ReadWeightFile( fWeightFiles[i], Z, A, I, Reaction);

		if( Z >= GetZAIThreshold() )
		{

            TGraph* ZAIXSTGraph = new TGraph();

			for(int TimeStep = 0; TimeStep < int(fMLP_Time.size()) ; TimeStep++)
			{
				double XSValue = ExecuteTMVA(fWeightFiles[i], IV, TimeStep );
                ZAIXSTGraph->SetPoint( ZAIXSTGraph->GetN(), (double)fMLP_Time[TimeStep], XSValue );
            }

            pair< map<ZAI, TGraph*>::iterator, bool> IResult;

            IResult = ExtrapolatedXS[Reaction].insert( pair<ZAI ,TGraph* >(ZAI(Z,A,I), ZAIXSTGraph) );

		}
	}
	
	/**********Sorting TGraph*********/
	for(int x = 0;x<3;x++)
	{	map<ZAI,TGraph*>::iterator it;
		for(it = ExtrapolatedXS[x].begin(); it != ExtrapolatedXS[x].end(); it++)
			it->second->Sort();
		
	}
	/**********Filling Matrices*/
	EvolutionDataFromMLP.SetFissionXS(ExtrapolatedXS[0]);
	EvolutionDataFromMLP.SetCaptureXS(ExtrapolatedXS[1]);
	EvolutionDataFromMLP.Setn2nXS(ExtrapolatedXS[2]);
	
	
	DBGL
	return EvolutionDataFromMLP;
}
//________________________________________________________________________
//________________________________________________________________________
//
//	Time step (1 MLP per time step)
//________________________________________________________________________
//________________________________________________________________________
void XSM_MLP::ReadWeightFileStep(string Filename, int &Z, int &A, int &I, int &Reaction, int &TimeStep)
{
	Z = -1;
	A = -1;
	I = -1;
	Reaction = -1;
	
	size_t found = Filename.find("XS");
	
	string NameJOB;
	NameJOB = Filename.substr(found);
	int pos = 0;
	
	StringLine::NextWord(NameJOB, pos, '_');
	
	Z = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	A = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	I = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	
	string  sReaction = (StringLine::NextWord(NameJOB,pos,'_') ).c_str() ;
	if(sReaction == "fis")
		Reaction = 0;
	if(sReaction == "cap")
		Reaction = 1;
	if(sReaction == "n2n")
		Reaction = 2;
	
	TimeStep = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	
	if(Z == -1 || A == -1 || I == -1 || Reaction == -1 || TimeStep == -1)
	{
		ERROR << " wrong TMVA weight format " << endl;
		exit(0);
	}
}

//________________________________________________________________________
EvolutionData XSM_MLP::GetCrossSectionsStep(IsotopicVector IV)
{
	DBGL

	EvolutionData EvolutionDataFromMLP = EvolutionData();
	
	map<ZAI,TGraph*> ExtrapolatedXS[3];
	/*************DATA BASE INFO****************/
	EvolutionDataFromMLP.SetReactorType("PWR");
	EvolutionDataFromMLP.SetFuelType("MOX");
	EvolutionDataFromMLP.SetPower(fDBPower);
	EvolutionDataFromMLP.SetHeavyMetalMass(fDBHMMass);
	
	/************* The Cross sections***********/
	
	for(int i = 0;i<int(fWeightFiles.size());i++)
	{
		int Z = -2;
		int A = -2;
		int I = -2;
		int Reaction = -2;
		int TimeStep = -2;
		ReadWeightFileStep( fWeightFiles[i], Z, A, I, Reaction, TimeStep);
		if( Z >= GetZAIThreshold() )
		{
			ZAI zaitmp = ZAI(Z,A,I);
			
			pair< map<ZAI, TGraph*>::iterator, bool> IResult;
			
			IResult = ExtrapolatedXS[Reaction].insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph() ) );
			
			if( IResult.second )
			{
				(IResult.first)->second->SetPoint(0, (double)fMLP_Time[TimeStep], ExecuteTMVA(fWeightFiles[i],IV) );
			}
			else
			{
				(IResult.first)->second->SetPoint( (IResult.first)->second->GetN(), (double)fMLP_Time[TimeStep], ExecuteTMVA(fWeightFiles[i],IV) );
			}
		}
	}
	
	/**********Sorting TGraph*********/
	for(int x = 0;x<3;x++)
	{	map<ZAI,TGraph*>::iterator it;
		for(it = ExtrapolatedXS[x].begin(); it != ExtrapolatedXS[x].end(); it++)
			it->second->Sort();
	}
	/**********Filling Matrices*/
	EvolutionDataFromMLP.SetFissionXS(ExtrapolatedXS[0]);
	EvolutionDataFromMLP.SetCaptureXS(ExtrapolatedXS[1]);
	EvolutionDataFromMLP.Setn2nXS(ExtrapolatedXS[2]);
	
	DBGL
	return EvolutionDataFromMLP;
}
//________________________________________________________________________
EvolutionData XSM_MLP::GetCrossSections(IsotopicVector IV ,double t)
{
	DBGL
	if(t != 0)
		WARNING << " Argument t has non effect here " << endl;
	
	EvolutionData EV;
	if(fIsStepTime)
		EV = GetCrossSectionsStep(IV);
	
	else
		EV = GetCrossSectionsTime(IV);
	
	DBGL
	return EV;
}
//________________________________________________________________________
