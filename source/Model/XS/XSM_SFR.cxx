#include <dirent.h>
#include <errno.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include <TGraph.h>
#include <TString.h>
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "XSModel.hxx"
#include "XSM_SFR.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "CLASSReader.hxx"
#include "external/StringLine.hxx"

//________________________________________________________________________
//
//		XSM_SFR
//
//________________________________________________________________________
XSM_SFR::XSM_SFR(string TMVA_Weight_Directory,map<string,double> FixedParameters,string InformationFile, bool IsTimeStep):XSModel(new CLASSLogger("XSM_SFR.log"))
{
	
	fIsStepTime = IsTimeStep;
	fTMVAWeightFolder = TMVA_Weight_Directory;
	
	fInformationFile = fTMVAWeightFolder+InformationFile;

	GetMLPWeightFiles();
	
	INFO << "__A cross section interpolator using" << endl;
	INFO << "Neural Network has been define__" << endl;
	INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
	
	LoadKeyword();
	ReadNFO();

	int number_of_elements = fTMVAVariableNames.size();
	double default_value = -1;
	vector<double> v(number_of_elements, default_value);
	fTMVAFixedVariableValues = v;
	for(unsigned int i=0;i<fTMVAVariableNames.size();i++){
		if(FixedParameters.find(fTMVAVariableNames[i]) != FixedParameters.end() ){
			fTMVAFixedVariableValues[i]=FixedParameters[fTMVAVariableNames[i]];
			fTMVAFixedVariable[i]=true;
		}
	}
}
//________________________________________________________________________
XSM_SFR::XSM_SFR(CLASSLogger* Log,string TMVA_Weight_Directory,map<string,double> FixedParameters,string InformationFile, bool IsTimeStep):XSModel(Log)
{
	
	fIsStepTime = IsTimeStep;
	fTMVAWeightFolder = TMVA_Weight_Directory;

	fInformationFile = fTMVAWeightFolder+InformationFile;
	
	GetMLPWeightFiles();
	
	INFO << "__A cross section interpolator using" << endl;
	INFO << "Neural Network has been define__" << endl;
	INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
	
	LoadKeyword();
	ReadNFO();

	int number_of_elements = fTMVAVariableNames.size();
	double default_value = -1;
	vector<double> v(number_of_elements, default_value);
	fTMVAFixedVariableValues = v;
	for(unsigned int i=0;i<fTMVAVariableNames.size();i++){
		if(FixedParameters.find(fTMVAVariableNames[i]) != FixedParameters.end() ){
			fTMVAFixedVariableValues[i]=FixedParameters[fTMVAVariableNames[i]];
			fTMVAFixedVariable[i]=true;
		}
	}
}
//________________________________________________________________________
XSM_SFR::~XSM_SFR()
{
	DBGL
	fTMVAVariableNames.clear();
	fDKeyword.clear();
	DBGL
}
//________________________________________________________________________
void XSM_SFR::SetFixedVariablesValues(map<string,double> FixedParameters)
{
	int number_of_elements = fTMVAVariableNames.size();
	double default_value = -1;
	vector<double> v(number_of_elements, default_value);
	fTMVAFixedVariableValues=v;
	for(unsigned int i=0;i<fTMVAVariableNames.size();i++){
		if(FixedParameters.find(fTMVAVariableNames[i]) != FixedParameters.end() ){
			fTMVAFixedVariableValues[i]=FixedParameters[fTMVAVariableNames[i]];
			fTMVAFixedVariable[i]=true;
		}
	}
}
//________________________________________________________________________
void XSM_SFR::LoadKeyword()
{
	DBGL
	fDKeyword.insert( pair<string, XS_SFR_DMthPtr>( "k_timestep",	& XSM_SFR::ReadTimeSteps));
	fDKeyword.insert( pair<string, XS_SFR_DMthPtr>( "k_inputparameter",	& XSM_SFR::ReadInputParameter)	 );
	DBGL
}
//________________________________________________________________________
void XSM_SFR::ReadInputParameter(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_inputparameter" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_inputparameter\" not found !" << endl;
		exit(1);
	}
	
	string name = StringLine::NextWord(line, pos, ' ');
	
	fTMVAVariableNames.push_back( name );

	string IsFixed = StringLine::NextWord(line, pos, ' ');

	if(IsFixed=="F"){
		fTMVAFixedVariable.push_back(true);
	}
	else{
		fTMVAFixedVariable.push_back(false);
	}
	
	DBGL
}
//________________________________________________________________________
void XSM_SFR::ReadTimeSteps(const string &line)
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
void XSM_SFR::ReadLine(string line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	
	map<string, XS_SFR_DMthPtr>::iterator it = fDKeyword.find(keyword);
	if(it != fDKeyword.end())
		(this->*(it->second))( line );
	
	DBGL
}
//________________________________________________________________________
void XSM_SFR::GetMLPWeightFiles()
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
void XSM_SFR::ReadWeightFile(string Filename, int &Z, int &A, int &I, int &Reaction)
{
	int tempZ = -1;
	int tempA = -1;
	int tempI = -1;
	int tempReaction = -1;
	
	size_t found = Filename.find("XS");
	string NameJOB;
	NameJOB = Filename.substr(found);
	int pos = 0;
	
	StringLine::NextWord(NameJOB, pos, '_');
	tempZ = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	tempA = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	tempI = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
	string  sReaction = (StringLine::NextWord(NameJOB,pos,'_') ).c_str() ;
	size_t foundext = sReaction.find(".weights.xml");
	sReaction = sReaction.substr(0,foundext);
	
	if(sReaction == "fis")
		tempReaction = 0;
	if(sReaction == "cap")
		tempReaction = 1;
	if(sReaction == "n2n")
		tempReaction = 2;
	
	if(tempZ<= 0 || tempA<= 0 || tempI<0 || tempReaction == -1)
	{
		ERROR << " wrong TMVA weight format " << endl;
		exit(0);
	}
	else{
		Z = tempZ;
		A = tempA;
		I = tempI;
		Reaction = tempReaction;
	}

}
//________________________________________________________________________
TTree* XSM_SFR::CreateTMVAInputTree(IsotopicVector isotopicvector,int TimeStep)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/
	TTree*   InputTree = new TTree("InTMP", "InTMP");
	
	vector<float> 	InputTMVA;
	for(unsigned int i = 0 ; i< (int)fTMVAVariableNames.size() ; i++){
		InputTMVA.push_back(0);
	}
	float Time = 0;
	IsotopicVector IVInputTMVA;
	vector<ZAI*> myZAIs;
	
	for(unsigned int j = 0 ; j< fTMVAVariableNames.size() ; j++)
	{
		InputTree->Branch( fTMVAVariableNames[j].c_str() ,&InputTMVA[j], (fTMVAVariableNames[j] + "/F").c_str());

		string ZAIname=fTMVAVariableNames[j];
		// For this function to work, ZAI name in the xml file have to be AAAName*
		// AAA is A, the number of nucleon
		// Name is the symbol of the element (with letters : U, Pu, Am, Cm...) and give Z
		// * => I=1, no * => I=0
		int Z = -1;
		int A = -1;
		int I = 0;

		if(ZAIname.back()=='*'){//Check last charcter to see if it's *
			I = 1;
			ZAIname.pop_back();
		}
		if (isdigit(ZAIname.c_str()[0])){//read all numbers in the name to find A
			A = atoi (ZAIname.c_str());
		}
		//remove all numbers in the ZAIname
		ZAIname.erase(remove_if(ZAIname.begin(), ZAIname.end(),[](unsigned char c){ return isdigit(c); }),ZAIname.end());
		//remaing (i.e. the text in ZAIname) is the element symbol
		if(ZAIname=="Am"){
			Z = 95;
		}
		else if(ZAIname=="Pu"){
			Z = 94;
		}
		else if(ZAIname=="U"){
			Z = 92;
		}
		if((Z<= 0 || A<= 0 )&& fTMVAFixedVariable[j]==false)
		{
			ERROR << " wrong TMVA weight format " << endl;
			exit(0);
		}
		ZAI* myZAI= new ZAI(Z,A,I);
		myZAIs.push_back(myZAI);

		if(fTMVAFixedVariable[j]==false){
			IVInputTMVA.Add(*myZAI,1);
		}
	}
	
	if( !fIsStepTime){
		InputTree->Branch(	"Time"	,&Time	,"Time/F"	);
	}
	IsotopicVector IVAccordingToUserInfoFile = isotopicvector.GetThisComposition(IVInputTMVA);
	double Ntot = IVAccordingToUserInfoFile.GetSumOfAll();
	IVAccordingToUserInfoFile = IVAccordingToUserInfoFile/Ntot;
	
	DBGV("INPUT TMVA");
	for(unsigned int j = 0 ; j< fTMVAVariableNames.size() ; j++)
	{
		if(fTMVAFixedVariable[j]){
			InputTMVA[j] = fTMVAFixedVariableValues[j];
		}
		else{
			InputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( *myZAIs[j] ) ;			
		}
	}	
	Time = fMLP_Time[TimeStep];

	InputTree->Fill();
	
	return InputTree;
}
//________________________________________________________________________
double XSM_SFR::ExecuteTMVA(string WeightFile,TTree* InputTree)
{
	DBGV( "File :" << WeightFile);
	// --- Create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "Silent" );
	
	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	vector<float> 	InputTMVA;
	for(unsigned int i = 0 ; i< fTMVAVariableNames.size() ; i++){
		InputTMVA.push_back(0);
	}
	Float_t Time;
	
	for(unsigned int j = 0 ; j< fTMVAVariableNames.size() ; j++){	
		reader->AddVariable( fTMVAVariableNames[j].c_str(),&InputTMVA[j]);
	}
	if(!fIsStepTime){
		reader->AddVariable( "Time" ,&Time);
	}
	
	// --- Book the MVA methods
	string dir    = fTMVAWeightFolder;
	if(dir[dir.size()-1] != '/'){
		dir+= "/";
	}
	
	// Book method MLP
	TString methodName = "MLP method";
	TString weightpath = dir + WeightFile ;
	reader->BookMVA( methodName, weightpath );
	
	map<ZAI ,string >::iterator it2;
	for(unsigned int j = 0 ; j< fTMVAVariableNames.size() ; j++)
	{
		InputTree->SetBranchAddress(fTMVAVariableNames[j].c_str(),&InputTMVA[j]);
	}
	if(!fIsStepTime){
		InputTree->SetBranchAddress( "Time" ,&Time );
	}
	InputTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];
	
	delete reader;
	DBGL
	
	return (double)val;
}
//________________________________________________________________________
EvolutionData XSM_SFR::GetCrossSectionsTime(IsotopicVector IV)
{
	DBGL
	
	string dir = fTMVAWeightFolder;
	if(dir[dir.size()-1] != '/'){
		dir+= "/";
	}
	
	EvolutionData EvolutionDataFromMLP = EvolutionData();
	
	map<ZAI,TGraph*> ExtrapolatedXS[3];
	/*************DATA BASE INFO****************/
	EvolutionDataFromMLP.SetReactorType(fDBRType);
	EvolutionDataFromMLP.SetFuelType(fDBFType);
	EvolutionDataFromMLP.SetPower(fDBPower);
	EvolutionDataFromMLP.SetHeavyMetalMass(fDBHMMass);
	/************* The Cross sections***********/
	for(unsigned int i = 0;i<fWeightFiles.size();i++)
	{
		int Z = -2;
		int A = -2;
		int I = -2;
		int Reaction = -2;
		ReadWeightFile( fWeightFiles[i], Z, A, I, Reaction);
		if( Z >= GetZAIThreshold() )
		{
			CLASSReader * reader = new CLASSReader( fTMVAVariableNames );
			if(!fIsStepTime) { reader->AddVariable( "Time" ); }
			
			reader->BookMVA( "MLP method" , dir + fWeightFiles[i] );
		
			for(unsigned int TimeStep = 0;TimeStep<fMLP_Time.size();TimeStep++)
			{
				TTree* InputTree = CreateTMVAInputTree(IV,TimeStep);
				reader->SetInputData( InputTree );
				double XSValue = reader->EvaluateRegression( "MLP method" )[0];
				
				pair< map<ZAI, TGraph*>::iterator, bool> IResult = ExtrapolatedXS[Reaction].insert( pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph()) );
				
				if(IResult.second )
				{
					(IResult.first)->second->SetPoint(0, (double)fMLP_Time[TimeStep], XSValue );
					
				}
				else
				{
					(IResult.first)->second->SetPoint( (IResult.first)->second->GetN(), (double)fMLP_Time[TimeStep], XSValue );
				}
				
				delete InputTree;
			}
			
			delete reader;
		}
	}
	
	/**********Sorting TGraph*********/
	for(unsigned int x = 0;x<3;x++)
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
//	Time step (1 NN per time step)
//________________________________________________________________________
//________________________________________________________________________
void XSM_SFR::ReadWeightFileStep(string Filename, int &Z, int &A, int &I, int &Reaction, int &TimeStep)
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
	
	if(Z == -1 || A == -1 || I == -1 || Reaction == -1 || TimeStep == -1){
		ERROR << " wrong TMVA weight format " << endl;
		exit(0);
	}
}

//________________________________________________________________________
EvolutionData XSM_SFR::GetCrossSectionsStep(IsotopicVector IV)
{
	DBGL
	TTree* InputTree = CreateTMVAInputTree(IV);
	EvolutionData EvolutionDataFromMLP = EvolutionData();
	
	map<ZAI,TGraph*> ExtrapolatedXS[3];
	/*************DATA BASE INFO****************/
	EvolutionDataFromMLP.SetReactorType(fDBRType);
	EvolutionDataFromMLP.SetFuelType(fDBFType);
	EvolutionDataFromMLP.SetPower(fDBPower);
	EvolutionDataFromMLP.SetHeavyMetalMass(fDBHMMass);
	
	/************* The Cross sections***********/
	for(unsigned int i = 0;i<fWeightFiles.size();i++){ // JM2016 : need to  booker the reader at each step because the file is different each time so we don't use CLASSReader 
		int Z = -2;
		int A = -2;
		int I = -2;
		int Reaction = -2;
		int TimeStep = -2;
		ReadWeightFileStep( fWeightFiles[i], Z, A, I, Reaction, TimeStep);
		if( Z >= GetZAIThreshold() ){
			ZAI zaitmp = ZAI(Z,A,I);
			
			pair< map<ZAI, TGraph*>::iterator, bool> IResult;
			IResult = ExtrapolatedXS[Reaction].insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph() ) );
			
			if( IResult.second ){
				(IResult.first)->second->SetPoint(0, (double)fMLP_Time[TimeStep], ExecuteTMVA(fWeightFiles[i],InputTree) );
			}
			else{
				(IResult.first)->second->SetPoint( (IResult.first)->second->GetN(), (double)fMLP_Time[TimeStep], ExecuteTMVA(fWeightFiles[i],InputTree) );
			}
		}
	}
	
	/**********Sorting TGraph*********/
	for(unsigned int x = 0;x<3;x++){
		map<ZAI,TGraph*>::iterator it;
		for(it = ExtrapolatedXS[x].begin(); it != ExtrapolatedXS[x].end(); it++){
			it->second->Sort();
		}
	}
	/**********Filling Matrices*/
	EvolutionDataFromMLP.SetFissionXS(ExtrapolatedXS[0]);
	EvolutionDataFromMLP.SetCaptureXS(ExtrapolatedXS[1]);
	EvolutionDataFromMLP.Setn2nXS(ExtrapolatedXS[2]);
	
	delete InputTree;
	DBGL
	return EvolutionDataFromMLP;
}
//________________________________________________________________________
EvolutionData XSM_SFR::GetCrossSections(IsotopicVector IV ,double t)
{
	DBGL
	if(t != 0){
		WARNING << " Argument t has non effect here " << endl;
	}

	EvolutionData EV;
	if(fIsStepTime){
		EV = GetCrossSectionsStep(IV);
	}
	else{
		EV = GetCrossSectionsTime(IV);
	}

	DBGL
	return EV;
}
//________________________________________________________________________
