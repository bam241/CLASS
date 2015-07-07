#include "EQM_FBR_MLP_Keff_BOUND.hxx"
#include "CLASSMethod.hxx"
#include "CLASSLogger.hxx"
#include "StringLine.hxx"

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

//________________________________________________________________________
//
//		EQM_FBR_MLP_Keff_BOUND
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For FBR
//
//________________________________________________________________________



//________________________________________________________________________
//
//	Objects & Methods for content prediction with calculation of Keff bounded
//	over batches.
//
//________________________________________________________________________

//________________________________________________________________________
EQM_FBR_MLP_Keff_BOUND::EQM_FBR_MLP_Keff_BOUND(string TMVAWeightPath,  int NumOfBatch, double LowerKeffective, double UpperKeffective, string InformationFile):EquivalenceModel(new CLASSLogger("EQM_FBR_MLP_Keff_BOUND.log"))
{
	DBGL

	/** The tmva weight **/
	
	fTMVAWeightPath = TMVAWeightPath;
	
	
	/* INFORMATION FILE HANDLING */
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");
	
	fInformationFile = InformationFile;
	
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile
	
	
	/* OTHER MODEL PARAMETERS */
	
	fNumberOfBatch = NumOfBatch;
	fKmin = LowerKeffective ;
	fKmax = UpperKeffective ;
	

	/* MODEL PARAMETERS INITIALIZATION */
	
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;
	
	
	/* INFO */

	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of keff averaged over the number of batch" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fMLPInformationFile << " )" << endl;

	
	DBGL
}
//________________________________________________________________________
EQM_FBR_MLP_Keff_BOUND::EQM_FBR_MLP_Keff_BOUND(CLASSLogger* log, string TMVAWeightPath,  int NumOfBatch, double LowerKeffective, double UpperKeffective, string InformationFile):EquivalenceModel(log)
{
	DBGL
	
	/** The tmva weight **/
	
	fTMVAWeightPath = TMVAWeightPath;
	
	
	/* INFORMATION FILE HANDLING */
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");
	
	fInformationFile = InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile
	
	
	
	/* OTHER MODEL PARAMETERS */
	
	fNumberOfBatch = NumOfBatch;
	fKmin = LowerKeffective ;
	fKmax = UpperKeffective ;
	
	
	/* INFO */
	
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of keff averaged over the number of batch" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fMLPInformationFile << " )" << endl;
	
	
	DBGL
}


//________________________________________________________________________
TGraph* EQM_FBR_MLP_Keff_BOUND::BuildKeffGraph(IsotopicVector FreshFuel)
{
	DBGL
	
	TGraph * keffGraph = new TGraph();
	for(int i = 0 ; i < (int) fMLP_Time.size() ; i++)
	{
		double keff_t = ExecuteTMVA( CreateTMVAInputTree(FreshFuel,(float) fMLP_Time[i]), true );
		keffGraph->SetPoint(i, (double)fMLP_Time[i], keff_t );
	}
	
	DBGL
	return keffGraph;
}
//________________________________________________________________________
TGraph* EQM_FBR_MLP_Keff_BOUND::BuildAverageKeffGraph(TGraph* GRAPH_KEFF)
{
	DBGL
	
	TGraph * AveragekeffGraph = new TGraph();
	int NumberOfPoint = 50;
	
	int NumOfInputGraphPoint = GRAPH_KEFF->GetN()-1 ;
	double TimeFinal = 0;
	double KFinal = 0;
	
	GRAPH_KEFF->GetPoint(NumOfInputGraphPoint, TimeFinal, KFinal);
	
	double step = TimeFinal/NumberOfPoint;
	int p = 0;
	for(int n = 0 ;n<NumberOfPoint;n++)
	{
		double k_av = 0;
		for(int b = 0;b<fNumberOfBatch;b++)
		{
			if((step*p)> TimeFinal/(double)fNumberOfBatch)
				p = 0;
			k_av += GRAPH_KEFF->Eval( (step*p + b*TimeFinal/(double)fNumberOfBatch) , 0 , "S" );
			
		}
		p++;
		k_av/= (double)fNumberOfBatch;
		
		AveragekeffGraph->SetPoint(n, step*n, k_av);
	}
	
	DBGL
	return AveragekeffGraph;
}
//________________________________________________________________________
double EQM_FBR_MLP_Keff_BOUND::GetKeffAt(TGraph* GRAPH_KEFF, int Step)
{
	DBGL
	
	double Time = 0;
	double Keff = 0;
	GRAPH_KEFF->GetPoint(Step, Time,  Keff);
	
	DBGL
	return Keff;
}


//________________________________________________________________________
TTree* EQM_FBR_MLP_Keff_BOUND::CreateTMVAInputTree(IsotopicVector TheFreshfuel, double ThisTime)
{
	DBGL

	/******Create Input data tree to be interpreted by TMVA::Reader***/
	TTree*   InputTree = new TTree("InTMPKef", "InTMPKef");
	
	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);
	
	float Time = 0;
	
	IsotopicVector IVInputTMVA;
	map<ZAI ,string >::iterator it;
	int j = 0;
	
	for( it = fMapOfTMVAVariableNames.begin()  ; it != fMapOfTMVAVariableNames.end() ; it++)
	{
		InputTree->Branch( ((*it).second).c_str()	,&InputTMVA[j], ((*it).second + "/F").c_str());
		IVInputTMVA+=  ((*it).first)*1;
		j++;
	}
	
	if(ThisTime != -1)
		InputTree->Branch(	"Time"	,&Time	,"Time/F"	);
	
	IsotopicVector IVAccordingToUserInfoFile = TheFreshfuel.GetThisComposition(IVInputTMVA);
	
	double Ntot = IVAccordingToUserInfoFile.GetSumOfAll();
	
	IVAccordingToUserInfoFile = IVAccordingToUserInfoFile/Ntot;
	
	j = 0;
	map<ZAI ,string >::iterator it2;
	
	for( it2 = fMapOfTMVAVariableNames.begin() ; it2 != fMapOfTMVAVariableNames.end() ; it2++)
	{
		InputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it2).first ) ;
		j++;
	}
	
	Time = ThisTime;
	
	InputTree->Fill();
	
	DBGL
	return InputTree;
	
}
//________________________________________________________________________
double EQM_FBR_MLP_Keff_BOUND::ExecuteTMVA(TTree* InputTree, bool IsTimeDependent)
{
	DBGL

	// --- Create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "Silent" );
	
	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);
	Float_t Time;
	
	map<ZAI ,string >::iterator it;
	int j = 0;
	for( it = fMapOfTMVAVariableNames.begin()  ; it != fMapOfTMVAVariableNames.end() ; it++)
	{	reader->AddVariable( ( (*it).second ).c_str(),&InputTMVA[j]);
		j++;
	}
	
	if(IsTimeDependent)
		reader->AddVariable( "Time" ,&Time);
	
	// --- Book the MVA methods
	
	// Book method MLP
	TString methodName = "MLP method";
	reader->BookMVA( methodName, fTMVAWeightPath );
	
	map<ZAI ,string >::iterator it2;
	j = 0;
	for( it2 = fMapOfTMVAVariableNames.begin()  ; it2 != fMapOfTMVAVariableNames.end() ; it2++)
	{
		InputTree->SetBranchAddress(( (*it2).second ).c_str(),&InputTMVA[j]);
		j++;
	}
	
	if(IsTimeDependent)
		InputTree->SetBranchAddress( "Time" ,&Time );
	
	InputTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];
	
	delete reader;
	
	DBGL
	return (double)val;	//return k_{eff}(t = Time)
}

//________________________________________________________________________
void EQM_FBR_MLP_Keff_BOUND::LoadKeyword()
{
	DBGL

	fDKeyword.insert( pair<string, FBR_MLP_Keff_BOUND_DMthPtr>( "k_timestep",	& EQM_FBR_MLP_Keff_BOUND::ReadTimeSteps));
	fDKeyword.insert( pair<string, FBR_MLP_Keff_BOUND_DMthPtr>( "k_zainame",	& EQM_FBR_MLP_Keff_BOUND::ReadZAIName)	 );
	
	DBGL
}


//________________________________________________________________________
void EQM_FBR_MLP_Keff_BOUND::ReadZAIName(const string &line)
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
void EQM_FBR_MLP_Keff_BOUND::ReadTimeSteps(const string &line)
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
void EQM_FBR_MLP_Keff_BOUND::ReadLine(string line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	
	map<string, FBR_MLP_Keff_BOUND_DMthPtr>::iterator it = fDKeyword.find(keyword);
	
	if(it != fDKeyword.end())
		(this->*(it->second))( line );
	
	DBGL
}


//________________________________________________________________________
double EQM_FBR_MLP_Keff_BOUND::GetFissileMolarFraction(IsotopicVector Fissile, IsotopicVector Fertile, double TargetBU)
{
	DBGL
	if(TargetBU != 0)
		WARNING << "The third arguement : Burnup has no effect here.";
	
	/**Algorithm  not so clever ...**/
	/**need improvements to make it faster*/
	
	if(!fIsAverageKeff)
	{
		ERROR << " Can't be used with EQM_FBR_MLP_Keff_BOUND(string TMVAWeightPath, double keff_target, string InformationFile) constructor";
		exit(1);
	}
	
	double FissileContent = 0.01;
	double test_Keff_beg = fKmin - 0.01;
	double test_keff_end = fKmax + 0.01;
	
	double speedstep = 1;
	
	while(fKmin >= test_Keff_beg || fKmax <=  test_keff_end)
	{
		IsotopicVector FreshFuel = (1-FissileContent)*(Fertile/Fertile.GetSumOfAll()) + FissileContent*(Fissile/Fissile.GetSumOfAll());
		
		TGraph* KEFF  = BuildKeffGraph(FreshFuel);
		TGraph* KEFF_avg = BuildAverageKeffGraph(KEFF);
		
		test_Keff_beg = GetKeffAt(KEFF_avg , 0);
		test_keff_end = GetKeffAt(KEFF_avg , fMLP_Time.size()-1);
		delete KEFF;
		delete KEFF_avg;
		
		if(test_Keff_beg < 0.9) //why 0.9 ? exactly
			speedstep = 10;
		else
			speedstep = 1;
		
		FissileContent+= 0.001*speedstep;
		
		if( test_Keff_beg > 1.30 )
		{
			ERROR << "This plutonium can not satisfy the criticality condition imposed" << endl;
			FissileContent = -1;
			break;
		}
		
	}
	
	DBGL
	return FissileContent;
}
//________________________________________________________________________
