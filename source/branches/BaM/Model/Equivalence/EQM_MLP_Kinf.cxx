#include "EQM_MLP_Kinf.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "StringLine.hxx"

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "TSystem.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"


//________________________________________________________________________
//
//		EQM_MLP_Kinf
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For REP MOX use
//
//________________________________________________________________________
EQM_MLP_Kinf::EQM_MLP_Kinf(string WeightPathAlpha0, string WeightPathAlpha1, string WeightPathAlpha2, string InformationFile,  int NumOfBatch, double CriticalityThreshold):EquivalenceModel(new CLASSLogger("EQM_MLP_Kinf.log"))
{
	/**The information file and tmva weight*/
	fTMVAWeightPath.push_back(WeightPathAlpha0);
	fTMVAWeightPath.push_back(WeightPathAlpha1);
	fTMVAWeightPath.push_back(WeightPathAlpha2);
	
	fMaximalBU = 0;
	fMaximalContent = 0;
	fInformationFile = InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile
	
	fNumberOfBatch = NumOfBatch;
	fKThreshold = CriticalityThreshold ;
	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetBuildFuelFirstGuess(0.04);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA weight  files are :" << endl;
	INFO << "\t" << fTMVAWeightPath[0] << endl;
	INFO << "\t" << fTMVAWeightPath[1] << endl;
	INFO << "\t" << fTMVAWeightPath[2] << endl;
	INFO << "\tThe Information file is :" << endl;
	INFO << "\t" << fInformationFile << endl;
	INFO << "Maximal fissile content (molar proportion) : " << fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " <<  fMaximalBU << endl;
	EquivalenceModel::PrintInfo();

	if(fMapOfTMVAVariableNames.empty() || fFertileList.GetIsotopicQuantity().empty() || fFissileList.GetIsotopicQuantity().empty() || fMaximalBU == 0 || fMaximalContent == 0 )	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}

}
//________________________________________________________________________
EQM_MLP_Kinf::EQM_MLP_Kinf(CLASSLogger* log, string WeightPathAlpha0, string WeightPathAlpha1, string WeightPathAlpha2, string InformationFile,  int NumOfBatch, double CriticalityThreshold):EquivalenceModel(log)
{
	/**The information file and tmva weight*/
	fTMVAWeightPath.push_back(WeightPathAlpha0);
	fTMVAWeightPath.push_back(WeightPathAlpha1);
	fTMVAWeightPath.push_back(WeightPathAlpha2);
	
	fMaximalBU = 0;
	fMaximalContent = 0;
	fInformationFile = InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile
	
    InitialiseTMVAReader();

    fNumberOfBatch = NumOfBatch;
	fKThreshold = CriticalityThreshold ;
	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetBuildFuelFirstGuess(0.04);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA weight  files are :" << endl;
	INFO << "\t" << fTMVAWeightPath[0] << endl;
	INFO << "\t" << fTMVAWeightPath[1] << endl;
	INFO << "\t" << fTMVAWeightPath[2] << endl;
	INFO << "\tThe Information file is :" << endl;
	INFO << "\t" << fInformationFile << endl;
	INFO << "Maximal fissile content (molar proportion) : " << fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " <<  fMaximalBU << endl;
	EquivalenceModel::PrintInfo();

	if(fMapOfTMVAVariableNames.empty() || fFertileList.GetIsotopicQuantity().empty() || fFissileList.GetIsotopicQuantity().empty() || fMaximalBU == 0 || fMaximalContent == 0 )	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}

}
//________________________________________________________________________
EQM_MLP_Kinf::EQM_MLP_Kinf(string TMVAWeightPath,  int NumOfBatch, string InformationFile, double CriticalityThreshold):EquivalenceModel(new CLASSLogger("EQM_MLP_Kinf.log"))
{
	/**The information file and tmva weight*/
	fTMVAWeightPath.push_back(TMVAWeightPath);
	
	/* INFORMATION FILE HANDLING */
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");
	
	fMaximalBU = 0;
	fMaximalContent = 0;
	fInformationFile = InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile

    InitialiseTMVAReader();

	fNumberOfBatch = NumOfBatch;
	fKThreshold = CriticalityThreshold ;
	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.04);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fInformationFile << " )" << endl;
	INFO << "Maximal fissile content (molar proportion) : " << fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " <<  fMaximalBU << endl;
	EquivalenceModel::PrintInfo();

	if(fMapOfTMVAVariableNames.empty() || fFertileList.GetIsotopicQuantity().empty() || fFissileList.GetIsotopicQuantity().empty() || fMaximalBU == 0 || fMaximalContent == 0 )	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}

}
//________________________________________________________________________
EQM_MLP_Kinf::EQM_MLP_Kinf(CLASSLogger* log, string TMVAWeightPath,  int NumOfBatch, string InformationFile, double CriticalityThreshold):EquivalenceModel(log)
{
	
	/**The information file and tmva weight*/
	fTMVAWeightPath.push_back(TMVAWeightPath);
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");
	
	fMaximalBU = 0;
	fMaximalContent = 0;
	fInformationFile = InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile

    InitialiseTMVAReader();

	fNumberOfBatch = NumOfBatch;
	fKThreshold = CriticalityThreshold ;
	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.04);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fInformationFile << " )" << endl;
	INFO << "Maximal fissile content (molar proportion) : " << fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " <<  fMaximalBU << endl;


	EquivalenceModel::PrintInfo();

	if(fMapOfTMVAVariableNames.empty() || fFertileList.GetIsotopicQuantity().empty() || fFissileList.GetIsotopicQuantity().empty() || fMaximalBU == 0 || fMaximalContent == 0 )	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}

	
}
//________________________________________________________________________
void EQM_MLP_Kinf::LoadKeyword()
{
	DBGL
	
	fDKeyword.insert( pair<string, PWR_MLP_KINF_DMthPtr>( "k_zainame", &EQM_MLP_Kinf::ReadZAIName) );
	fDKeyword.insert( pair<string, PWR_MLP_KINF_DMthPtr>( "k_maxburnup", &EQM_MLP_Kinf::ReadMaxBurnUp) );
	fDKeyword.insert( pair<string, PWR_MLP_KINF_DMthPtr>( "k_maxfiscontent", &EQM_MLP_Kinf::ReadMaxFisContent) );
	DBGL
}
//________________________________________________________________________
void EQM_MLP_Kinf::ReadZAIName(const string &line)
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
void EQM_MLP_Kinf::ReadMaxBurnUp(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_maxburnup" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_maxburnup\" not found !" << endl;
		exit(1);
	}
	
	fMaximalBU = atof(StringLine::NextWord(line, pos, ' ').c_str());

	DBGL
}
//________________________________________________________________________
void EQM_MLP_Kinf::ReadMaxFisContent(const string &line)
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
void EQM_MLP_Kinf::ReadLine(string line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	
	map<string, PWR_MLP_KINF_DMthPtr>::iterator it = fDKeyword.find(keyword);
	
	if(it != fDKeyword.end())
		(this->*(it->second))( line );
	
	DBGL
}

void EQM_MLP_Kinf::InitialiseTMVAReader()
{

    for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
        InputTMVA.push_back(0);

    reader = new TMVA::Reader( "Silent" );


    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    vector<float> 	InputTMVA;
    for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
        InputTMVA.push_back(0);
    Float_t Time;

    map<ZAI ,string >::iterator it;
    int j = 0;
    for( it = fMapOfTMVAVariableNames.begin()  ; it != fMapOfTMVAVariableNames.end() ; it++)
        {
            reader->AddVariable( ( (*it).second ).c_str(),&InputTMVA[j]);
            IVInputTMVA +=  ((*it).first)*1;
            j++;
        }

    if(fTMVAWeightPath.size() == 1)
        reader->AddVariable( "Time" ,&Time);

}


void EQM_MLP_Kinf::UpdateInputComposition(IsotopicVector TheFreshfuel, double ThisTime)
{


    IsotopicVector IVAccordingToUserInfoFile = TheFreshfuel.GetThisComposition(IVInputTMVA);


    map<ZAI,string>::iterator it;
    int j = 0;

	for( it = fMapOfTMVAVariableNames.begin() ; it != fMapOfTMVAVariableNames.end() ; it++)
	{
		InputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it).first ) ;
		j++;
	}
	
	Time = ThisTime;

}
//________________________________________________________________________
double EQM_MLP_Kinf::ExecuteTMVA(IsotopicVector TheFreshfuel, double ThisTime, string TMVAWeightPath)
{
    UpdateInputComposition(TheFreshfuel,ThisTime);

    reader->BookMVA( "MLP Method", TMVAWeightPath );

    Float_t val = (reader->EvaluateRegression(  "MLP Method" ))[0];

	return (double)val;//retourn k_{inf}(t = Time)
}
//________________________________________________________________________
double EQM_MLP_Kinf::GetMaximumBurnUp_MLP(IsotopicVector TheFuel, double TargetBU)
{
	/**************************************************************************/
	//With a dichotomy, the maximal irradiation time (TheFinalTime) is calculated
	//When average Kinf is very close (according "Precision") to the threshold
	//then the corresponding irradiation time is convert in burnup and returned
	/**************************************************************************/
	//Algorithm initialization
	double TheFinalTime = BurnupToSecond(TargetBU);
	double OldFinalTimeMinus = 0;
	double MaximumBU = fMaximalBU;
	double OldFinalTimePlus = BurnupToSecond(MaximumBU);
	double k_av = 0; //average kinf
	double OldPredictedk_av = 0;
	for(int b = 0;b<fNumberOfBatch;b++)
	{
		float TheTime = (b+1)*TheFinalTime/fNumberOfBatch;
        OldPredictedk_av += ExecuteTMVA(TheFuel,TheTime, fTMVAWeightPath[0]);
	}
	OldPredictedk_av/= fNumberOfBatch;
	
	//Algorithm control
	int count = 0;
	int MaximumLoopCount = 500;
	do
	{
		if(count > MaximumLoopCount )
		{
			ERROR << "CRITICAL ! Can't manage to predict burnup\nHint : Try to increase the precision on k effective using :\n YourEQM_MLP_Kinf->SetPCMprecision(pcm); with pcm the precision in pcm (default 10) REDUCE IT\n If this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
			exit(1);
		}
		
		if( (OldPredictedk_av-fKThreshold)  > 0) //The burnup can be increased
		{
			OldFinalTimeMinus = TheFinalTime;
			TheFinalTime = TheFinalTime + fabs(OldFinalTimePlus - TheFinalTime)/2.;
			
			if(SecondToBurnup(TheFinalTime) >= (MaximumBU-MaximumBU*GetBurnUpPrecision() ) )
				return MaximumBU;
		}
		
		else if( (OldPredictedk_av-fKThreshold)  < 0)//The burnup is too high
		{
			OldFinalTimePlus = TheFinalTime;
			TheFinalTime = TheFinalTime - fabs(OldFinalTimeMinus-TheFinalTime)/2.;
			if( SecondToBurnup(TheFinalTime) < TargetBU*GetBurnUpPrecision() )
				return 0;
		}
		
		k_av = 0;
		for(int b = 0;b<fNumberOfBatch;b++)
		{
			float TheTime = (b+1)*TheFinalTime/fNumberOfBatch;
			k_av += ExecuteTMVA(TheFuel,TheTime, fTMVAWeightPath[0]);
		}
		k_av/= fNumberOfBatch;
		
		OldPredictedk_av = k_av;
		count++;
	}	while( fabs(OldPredictedk_av-fKThreshold) > GetPCMprecision() )  ;
	
	
	return SecondToBurnup(TheFinalTime);
}
//________________________________________________________________________
double EQM_MLP_Kinf::GetMaximumBurnUp_Pol2(IsotopicVector TheFuel,double TargetBU)
{
	
	double Alpha_0 = ExecuteTMVA(TheFuel, -1, fTMVAWeightPath[0]);
	double Alpha_1 = ExecuteTMVA(TheFuel, -1, fTMVAWeightPath[1]);
	double Alpha_2 = ExecuteTMVA(TheFuel, -1, fTMVAWeightPath[2]);

	if(Alpha_0 < fKThreshold) //not enought fissile for sure !!
		return 0;
	
	double Sum = 0;
	double SumSquare = 0;
	
	for(int i = 1 ; i<=   fNumberOfBatch; i++)
	{	Sum+= double(i)/double(fNumberOfBatch);
		SumSquare+= double(i*i)/double(fNumberOfBatch*fNumberOfBatch);
	}
	
	double C = (Alpha_0-fKThreshold)*fNumberOfBatch;
	double B = Alpha_1*Sum;
	double A = Alpha_2*SumSquare;
	
	double Delta = B*B-4*A*C;
	
	double T = 0;
	if( Delta > 0)
	{
		double T_1 = (-B + sqrt(Delta))/(2*A);
		double T_2 = (-B - sqrt(Delta))/(2*A);
		
		if(T_1 < 0 && T_2 > 0)
			T = T_2;
		else if(T_1 > 0 && T_2 < 0)
			T = T_1;
		
		else if( T_1 > 0 && T_2 > 0 )
		{
			if(T_2 < T_1)
				T = T_2;
			else
				T = T_1;
		}
		else
		{
			ERROR << "No positive solution" << endl;
			exit(1);
		}
	}
	else if(Delta == 0)
	{	T = -B/(2*A);
		if(T<0)
		{	ERROR << "No positive solution" << endl;
			exit(1);
		}
	}
	else
	{
		WARNING << "No real solution" << endl;
		double K_LongTime = Alpha_0+BurnupToSecond(10*TargetBU)*Alpha_1+Alpha_2*BurnupToSecond(10*TargetBU)*BurnupToSecond(10*TargetBU);
		DBGV("K_LongTime " << K_LongTime)
		
		if(K_LongTime > fKThreshold)
			return 10000;
		else
		{
			ERROR << " CRITICAL ! Can't find a physical solution ! \n Should not happening please contact BLG :" << endl;
			ERROR << "mail to baptiste.leniau@subatech.in2p2.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
			exit(1);
		}
		
	}
	
	return SecondToBurnup(T);
}
//________________________________________________________________________
double EQM_MLP_Kinf::GetMaximumBurnUp(IsotopicVector TheFuel, double TargetBU)
{
	double TheBurnUp = -1;
	if(fTMVAWeightPath.size() == 1)
		TheBurnUp  = 	GetMaximumBurnUp_MLP(TheFuel,TargetBU);
	else if(fTMVAWeightPath.size() == 3)
		TheBurnUp = GetMaximumBurnUp_Pol2(TheFuel,TargetBU);
	
	else
	{
		ERROR << "This method is not yet set up" << endl;
		exit(0);
	}
	
	return TheBurnUp;
}
//________________________________________________________________________
double EQM_MLP_Kinf::GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double TargetBU)
{
	
	//initialization
	double FissileContent = GetActualFissileContent();
	double OldFissileContentMinus = 0;
	double OldFissileContentPlus = fMaximalContent;
	double PredictedBU = 0 ;
	IsotopicVector FreshFuel = (1-FissileContent)*(Fertil/Fertil.GetSumOfAll()) + FissileContent*(Fissil/Fissil.GetSumOfAll());
	double OldPredictedBU = GetMaximumBurnUp(FreshFuel,TargetBU);
	
	double Precision = GetBurnUpPrecision()*TargetBU; //1 % of the targeted burnup
	int count = 0;
	int MaximumLoopCount = 500;
	do
	{
		if(count > MaximumLoopCount )
		{
			ERROR << "CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on burnup using :\nYourEQM_MLP_Kinf->SetBurnUpPrecision(prop); with prop the precision  (default 0.5percent :  0.005) INCREASE IT\nIf this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
			ERROR << "Targeted Burnup :" <<TargetBU<<endl;
			ERROR << "Last calculated Burnup :" <<OldPredictedBU<<endl;
			ERROR << "Last Fresh fuel composition :" <<endl;
			ERROR << FreshFuel.sPrint()<<endl;
			
			exit(1);
		}
		
		if( (OldPredictedBU - TargetBU) < 0 ) //The Content can be increased
		{
			OldFissileContentMinus = FissileContent;
			FissileContent = FissileContent + fabs(OldFissileContentPlus-FissileContent)/2.;
		}
		else if( (OldPredictedBU - TargetBU) > 0) //The Content is too high
		{
			OldFissileContentPlus = FissileContent;
			FissileContent = FissileContent - fabs(OldFissileContentMinus-FissileContent)/2.;
		}
		
		IsotopicVector FreshFuel = (1-FissileContent)*(Fertil/Fertil.GetSumOfAll()) + FissileContent*(Fissil/Fissil.GetSumOfAll());
		PredictedBU = GetMaximumBurnUp(FreshFuel,TargetBU);
		
		OldPredictedBU = PredictedBU;
		count ++;
		
	}while(fabs(TargetBU-PredictedBU)>Precision);
	
	
	
	DBGV("Predicted BU " << PredictedBU << " FissileContent " << FissileContent);
	return FissileContent;
}
//________________________________________________________________________