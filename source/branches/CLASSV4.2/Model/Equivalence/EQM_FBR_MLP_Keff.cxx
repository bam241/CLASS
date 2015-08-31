#include "EQM_FBR_MLP_Keff.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "StringLine.hxx"

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <map>

#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

//________________________________________________________________________
//
//		EQM_FBR_MLP_Keff
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For FBR 
//
//________________________________________________________________________


//________________________________________________________________________
//
//	Objects & Methods for content prediction with calculation of Keff at
//	one given time (often BOC or EOC)
//
//________________________________________________________________________

//________________________________________________________________________
EQM_FBR_MLP_Keff::EQM_FBR_MLP_Keff(string TMVAWeightPath, double keff_target, string InformationFile):EquivalenceModel(new CLASSLogger("EQM_FBR_MLP_Keff.log"))
{
	DBGL

	/**The tmva weight*/
	fTMVAWeightPath = TMVAWeightPath;


	/* INFORMATION FILE HANDLING */
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");

	fMaximalContent = 0;
	fInformationFile = InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile
	
	if(fMaximalContent == 0 )
	{
		ERROR<<"Can't find the k_maxfiscontent keyword in .nfo file\n this is mandatory"<<endl;
		exit(0);
	}
	
	fTargetKeff = keff_target;
	
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);	//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;

	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of keff at a specific time" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fInformationFile << " )" << endl;

	DBGL
}
//________________________________________________________________________
EQM_FBR_MLP_Keff::EQM_FBR_MLP_Keff(CLASSLogger* log, string TMVAWeightPath, double keff_target, string InformationFile):EquivalenceModel(log)
{
	DBGL

	/**The tmva weight*/
	fTMVAWeightPath = TMVAWeightPath;


	/* INFORMATION FILE HANDLING */
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");

	fMaximalContent = 0;
	fInformationFile = InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile
	
	if(fMaximalContent == 0 )
	{
		ERROR<<"Can't find the k_maxfiscontent keyword in .nfo file\n this is mandatory"<<endl;
		exit(0);
	}

	fTargetKeff = keff_target;
	
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);	//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;

	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of keff at a specific time" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | "  << fInformationFile << " )" << endl;

	DBGL
}


//________________________________________________________________________
TTree* EQM_FBR_MLP_Keff::CreateTMVAInputTree(IsotopicVector TheFreshfuel, double ThisTime)
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
double EQM_FBR_MLP_Keff::ExecuteTMVA(TTree* InputTree, bool IsTimeDependent)
{
	DBGL

	// --- Create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);
	
	Float_t Time = 0;

	
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
void EQM_FBR_MLP_Keff::LoadKeyword()
{
	DBGL

	fDKeyword.insert( pair<string, FBR_MLP_Keff_DMthPtr>( "k_zainame",	& EQM_FBR_MLP_Keff::ReadZAIName)	 );

	fDKeyword.insert( pair<string, FBR_MLP_Keff_DMthPtr>( "k_maxfiscontent", &EQM_FBR_MLP_Keff::ReadMaxFisContent) );

	DBGL
}


//________________________________________________________________________
void EQM_FBR_MLP_Keff::ReadZAIName(const string &line)
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
	
	fFissileList.Add(Z, A, I, 1.0);
	fStreamList.push_back(fFissileList);
	
	DBGL
}

//________________________________________________________________________
void EQM_FBR_MLP_Keff::ReadMaxFisContent(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_maxfiscontent" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_maxfiscontent\" not found !" << endl;
		exit(1);
	}
	
	fMaximalContent = atof(StringLine::NextWord(line, pos, ' ').c_str());
	
	DBGL
}

//________________________________________________________________________
void EQM_FBR_MLP_Keff::ReadLine(string line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	
	map<string, FBR_MLP_Keff_DMthPtr>::iterator it = fDKeyword.find(keyword);
	
	if(it != fDKeyword.end())
		(this->*(it->second))( line );
	
	DBGL
}




//________________________________________________________________________
map < string , double> EQM_FBR_MLP_Keff::GetMolarFraction(vector <IsotopicVector> IVStream,double TargetBU)
{
	DBGL
	
	if(TargetBU != 0)
		WARNING << "The third arguement : Burnup has no effect here.";

	IsotopicVector Fissile = IVStream["Fissile"];
	IsotopicVector Fertile = IVStream["Fertile"];

	//initialization
	double FissileContent = GetActualFissileContent();
	double OldFissileContentMinus = 0;
	double OldFissileContentPlus = fMaximalContent;
	double PredictedKeff = 0 ;
	IsotopicVector FreshFuel = (1-FissileContent)*(Fertile/Fertile.GetSumOfAll()) + FissileContent*(Fissile/Fissile.GetSumOfAll());
	double OldPredictedKeff = GetKeffAtFixedTime(FreshFuel);
	
	double Precision = fPCMprecision/1e5*fTargetKeff; //pcm to 1
	
	int count = 0;
	int MaximumLoopCount = 100;
	do
	{
		if(count > MaximumLoopCount )
		{
			ERROR << "CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on keff using :\nYourEQM_FBR_MLP_Keff->SetPCMPrecision(prop); with prop the precision  (default 0.5percent :  0.005) INCREASE IT \n If this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
			exit(1);
		}
		
		if( (OldPredictedKeff - fTargetKeff) < 0 ) //The Content can be increased
		{
			OldFissileContentMinus = FissileContent;
			FissileContent = FissileContent + fabs(OldFissileContentPlus-FissileContent)/2.;
		}
		else if( (OldPredictedKeff - fTargetKeff) > 0) //The Content is too high
		{
			OldFissileContentPlus = FissileContent;
			FissileContent = FissileContent - fabs(OldFissileContentMinus-FissileContent)/2.;
		}
		
		IsotopicVector FreshFuel = (1-FissileContent)*(Fertile/Fertile.GetSumOfAll()) + FissileContent*(Fissile/Fissile.GetSumOfAll());
		
		PredictedKeff = GetKeffAtFixedTime(FreshFuel);
		
		OldPredictedKeff = PredictedKeff;
		count ++;
		
	}while(fabs(fTargetKeff-PredictedKeff)>Precision);
	
	DBGV( "Predicted keff " << PredictedKeff << " FissileContent " << FissileContent << endl);
	
	map < string , double> MolarFraction;
	MolarFraction["Fissile"] = FissileContent;
	MolarFraction["Fertile"] = 1.- FissileContent;

	return MolarFraction; //return Molar content of each component in the fuel


}
//________________________________________________________________________
