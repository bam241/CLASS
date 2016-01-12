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

    InitialiseTMVAReader();

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
	INFO << "Maximal fissile content (molar proportion) : "<<fMaximalContent<<endl;
	EquivalenceModel::PrintInfo();

	if(fMapOfTMVAVariableNames.empty() || fFertileList.GetIsotopicQuantity().empty() || fFissileList.GetIsotopicQuantity().empty())	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}

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

    InitialiseTMVAReader();

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
	EquivalenceModel::PrintInfo();

	if(fMapOfTMVAVariableNames.empty() || fFertileList.GetIsotopicQuantity().empty() || fFissileList.GetIsotopicQuantity().empty())	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}

	DBGL
}

EQM_FBR_MLP_Keff::~EQM_FBR_MLP_Keff()
{
    delete freader;
}

void EQM_FBR_MLP_Keff::InitialiseTMVAReader()
{

    for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
        fInputTMVA.push_back(0);

    freader = new TMVA::Reader( "Silent" );


    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
        fInputTMVA.push_back(0);

    map<ZAI ,string >::iterator it;
    int j = 0;
    for( it = fMapOfTMVAVariableNames.begin()  ; it != fMapOfTMVAVariableNames.end() ; it++)
    {
        freader->AddVariable( ( (*it).second ).c_str(),&fInputTMVA[j]);
        fIVInputTMVA +=  ((*it).first)*1;
        j++;
    }

}

void EQM_FBR_MLP_Keff::UpdateInputComposition(IsotopicVector TheFreshfuel)
{


    IsotopicVector IVAccordingToUserInfoFile = TheFreshfuel.GetThisComposition(fIVInputTMVA);


    map<ZAI,string>::iterator it;
    int j = 0;

    for( it = fMapOfTMVAVariableNames.begin() ; it != fMapOfTMVAVariableNames.end() ; it++)
    {
        fInputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it).first ) ;
        j++;
    }

}


//________________________________________________________________________
double EQM_FBR_MLP_Keff::ExecuteTMVA(IsotopicVector TheFreshfuel)
{
	DBGL
	
	Float_t val = (freader->EvaluateRegression( "MLP Method" ))[0];

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
	
	string name = StringLine::NextWord(line, pos, ' ');
	
	fMapOfTMVAVariableNames.insert( pair<ZAI,string>( ZAI(Z, A, I), name ) );

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
double EQM_FBR_MLP_Keff::GetFissileMolarFraction(IsotopicVector Fissile,IsotopicVector Fertile,double TargetBU)
{
	DBGL
	
	if(TargetBU != 0)
		WARNING << "The third arguement : Burnup has no effect here.";


	//initialization
	double FissileContent = GetActualFissileContent();
	double OldFissileContentMinus = 0;
	double OldFissileContentPlus = fMaximalContent;
	double PredictedKeff = 0 ;
	IsotopicVector FreshFuel = (1-FissileContent)*(Fertile/Fertile.GetSumOfAll()) + FissileContent*(Fissile/Fissile.GetSumOfAll());
	double OldPredictedKeff = ExecuteTMVA(FreshFuel);
	
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
		
		PredictedKeff = ExecuteTMVA(FreshFuel);
		
		OldPredictedKeff = PredictedKeff;
		count ++;
		
	}while(fabs(fTargetKeff-PredictedKeff)>Precision);
	
	DBGV( "Predicted keff " << PredictedKeff << " FissileContent " << FissileContent << endl);
	return FissileContent;

}
//________________________________________________________________________