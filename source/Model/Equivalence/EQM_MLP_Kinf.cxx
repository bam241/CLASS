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
#include "TMVA/Reader.h"
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
	
	fMaximalBU 		= 0;
	fMaximalContent 	= 0;
	fNumberOfBatch 	= NumOfBatch;
	fKThreshold 		= CriticalityThreshold ;

	SetBurnUpPrecision(0.005);//1 % of the targeted burnup

	fInformationFile 	= InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile
	
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

	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	bool EmptyList = false;
	bool FirstContentNULL = false;

	for(  it_s_IV = fStreamList.begin();   it_s_IV != fStreamList.end();  it_s_IV++)
	{	
		if(fStreamList[(* it_s_IV).first].GetIsotopicQuantity().empty()){EmptyList=true;}
	}

	for(  it_s_D = fFirstGuessContent .begin();   it_s_D != fFirstGuessContent.end();  it_s_D++)
	{	
		if(fFirstGuessContent[(* it_s_D).first] ==0){FirstContentNULL=true;}
	}


	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || FirstContentNULL==true || fMaximalBU == 0 || fMaximalContent == 0 )	
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
	
	fMaximalBU 		= 0;
	fMaximalContent 	= 0;
	fNumberOfBatch 	= NumOfBatch;
	fKThreshold 		= CriticalityThreshold ;

	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	
	fInformationFile 	= InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile
	
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

	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	bool EmptyList = false;
	bool FirstContentNULL = false;

	for(  it_s_IV = fStreamList.begin();   it_s_IV != fStreamList.end();  it_s_IV++)
	{	
		if(fStreamList[(* it_s_IV).first].GetIsotopicQuantity().empty()){EmptyList=true;}
	}

	for(  it_s_D = fFirstGuessContent .begin();   it_s_D != fFirstGuessContent.end();  it_s_D++)
	{	
		if(fFirstGuessContent[(* it_s_D).first] ==0){FirstContentNULL=true;}
	}


	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || FirstContentNULL==true || fMaximalBU == 0 || fMaximalContent == 0 )	
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
	
	fMaximalBU 		= 0;
	fMaximalContent 	= 0;
	fNumberOfBatch 	= NumOfBatch;
	fKThreshold 		= CriticalityThreshold ;

	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetPCMPrecision(10);
	
	fInformationFile 	= InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fInformationFile << " )" << endl;
	INFO << "Maximal fissile content (molar proportion) : " << fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " <<  fMaximalBU << endl;
	EquivalenceModel::PrintInfo();

	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	bool EmptyList 		= false;
	bool FirstContentNULL 	= false;

	for(  it_s_IV = fStreamList.begin();   it_s_IV != fStreamList.end();  it_s_IV++)
	{	
		if(fStreamList[(* it_s_IV).first].GetIsotopicQuantity().empty()){EmptyList=true;}
	}

	for(  it_s_D = fFirstGuessContent .begin();   it_s_D != fFirstGuessContent.end();  it_s_D++)
	{	
		if(fFirstGuessContent[(* it_s_D).first] ==0){FirstContentNULL=true;}
	}


	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || FirstContentNULL==true || fMaximalBU == 0 || fMaximalContent == 0 )	
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
	
	fMaximalBU 		= 0;
	fMaximalContent 	= 0;
	fInformationFile 	= InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile
	
	fNumberOfBatch 	= NumOfBatch;
	fKThreshold 		= CriticalityThreshold ;
	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetPCMPrecision(10);
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fInformationFile << " )" << endl;
	INFO << "Maximal fissile content (molar proportion) : " << fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " <<  fMaximalBU << endl;

	EquivalenceModel::PrintInfo();

	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	bool EmptyList 		= false;
	bool FirstContentNULL 	= false;

	for(  it_s_IV = fStreamList.begin();   it_s_IV != fStreamList.end();  it_s_IV++)
	{	
		if(fStreamList[(* it_s_IV).first].GetIsotopicQuantity().empty()){EmptyList=true;}
	}

	for(  it_s_D = fFirstGuessContent .begin();   it_s_D != fFirstGuessContent.end();  it_s_D++)
	{	
		if(fFirstGuessContent[(* it_s_D).first] ==0){FirstContentNULL=true;}
	}


	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || FirstContentNULL==true || fMaximalBU == 0 || fMaximalContent == 0 )	
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
	int I  = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	
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
//________________________________________________________________________
TTree* EQM_MLP_Kinf::CreateTMVAInputTree(IsotopicVector TheFreshfuel, double ThisTime)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/
	TTree*   InputTree = new TTree("InTMPKinf", "InTMPKinf");
	
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
	
	return InputTree;
	
}
//________________________________________________________________________
double EQM_MLP_Kinf::ExecuteTMVA(TTree* InputTree,string WeightPath, bool IsTimeDependent)
{
	
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
	reader->BookMVA( methodName, WeightPath );
	
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
		TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
		OldPredictedk_av += ExecuteTMVA(InputTree,fTMVAWeightPath[0],true);
		delete InputTree;
	}
	OldPredictedk_av/= fNumberOfBatch;
	
	//Algorithm control
	int count = 0;
	int MaximumLoopCount = 500;
	do
	{
		if(count > MaximumLoopCount )
		{
			ERROR << "CRITICAL ! Can't manage to predict burnup\nHint : Try to increase the precision on k effective using :\n YourEQM_MLP_Kinf->SetPCMPrecision(pcm); with pcm the precision in pcm (default 10) REDUCE IT\n If this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
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
			TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
			k_av += ExecuteTMVA(InputTree,fTMVAWeightPath[0],true);
			delete InputTree;
		}
		k_av/= fNumberOfBatch;
		//cout<<SecondToBurnup(TheFinalTime)<<" ";
		OldPredictedk_av = k_av;
		count++;
	}	while( fabs(OldPredictedk_av-fKThreshold) > GetPCMPrecision() )  ;
	
	//cout<<endl;
	return SecondToBurnup(TheFinalTime);
}
//________________________________________________________________________
double EQM_MLP_Kinf::GetMaximumBurnUp_Pol2(IsotopicVector TheFuel,double TargetBU)
{
	
	TTree* InputTree = CreateTMVAInputTree(TheFuel,-1);
	double Alpha_0 = ExecuteTMVA(InputTree,fTMVAWeightPath[0],false);
	double Alpha_1 = ExecuteTMVA(InputTree,fTMVAWeightPath[1],false);
	double Alpha_2 = ExecuteTMVA(InputTree,fTMVAWeightPath[2],false);
	delete InputTree;
	
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
map < string , double>  EQM_MLP_Kinf::GetMolarFraction(map <string , IsotopicVector > IVStream, double TargetBU)
{
	//initialization
	IsotopicVector Fissil 	= IVStream["Fissile"];
	IsotopicVector Fertil 	= IVStream["Fertile"];

	map < string , double> MolarFraction ;

	double FissileContent 		= fActualMolarContentInFuel["Fissile"];
	double OldFissileContentMinus  	= 0;
	double OldFissileContentPlus 	= fMaximalContent;

	double PredictedBU 		= 0 ;
	IsotopicVector FreshFuel 	= (1-FissileContent)*(Fertil/Fertil.GetSumOfAll()) + FissileContent*(Fissil/Fissil.GetSumOfAll());
	double OldPredictedBU 		= GetMaximumBurnUp(FreshFuel,TargetBU);
	
	double Precision 		= GetBurnUpPrecision()*TargetBU; //1 % of the targeted burnup
	int count 			= 0;
	int MaximumLoopCount 	= 500;
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
	
	MolarFraction["Fissile"] 	= FissileContent;
	MolarFraction["Fertile"] 	= 1 - FissileContent;

	return MolarFraction;
}
//________________________________________________________________________