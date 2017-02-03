#include "EQM_MLP_Kinf.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "external/StringLine.hxx"

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

#include "CLASSReader.hxx"


//_________________________________________________________________________________
//
//		EQM_MLP_Kinf
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For REP MOX use
//
//_________________________________________________________________________________

//_________________________________________________________________________________
EQM_MLP_Kinf::EQM_MLP_Kinf(string TMVAWeightPath,  int NumOfBatch, string InformationFile, double CriticalityThreshold):EquivalenceModel(new CLASSLogger("EQM_MLP_Kinf.log"))
{
	/**The information file and tmva weight*/
	fTMVAWeightPath.push_back(TMVAWeightPath);
	
	/* INFORMATION FILE HANDLING */
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");
	
	fMaximalBU 		= 0;
	fNumberOfBatch 	= NumOfBatch;
	fKThreshold 		= CriticalityThreshold ;

	SetBurnUpPrecision(0.005);//0.8 % of the targeted burnup
	SetPCMPrecision(10);
	
	fInformationFile 	= InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fInformationFile << " )" << endl;
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

	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || fMaximalBU == 0 )	
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
	fInformationFile 	= InformationFile;
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile
	
	fNumberOfBatch 	= NumOfBatch;
	fKThreshold 		= CriticalityThreshold ;

	SetBurnUpPrecision(0.005);//0.8 % of the targeted burnup
	SetPCMPrecision(10);
	
	INFO << "__An equivalence model has been define__" << endl;
	INFO << "\tThis model is based on the prediction of kinf" << endl;
	INFO << "\tThe TMVA (weight | information) files are :" << endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " << fInformationFile << " )" << endl;
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


	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || fMaximalBU == 0 )	
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
		InputTree->Branch( ((*it).second).c_str(), &InputTMVA[j], ((*it).second + "/F").c_str());
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
double EQM_MLP_Kinf::CalculateTargetParameter(IsotopicVector TheFuel)
{
	/**************************************************************************/
	//With a dichotomy, the maximal irradiation time (TheFinalTime) is calculated
	//When average Kinf is very close (according "Precision") to the threshold
	//then the corresponding irradiation time is convert in burnup and returned
	/**************************************************************************/
	//Algorithm initialization
	double OldFinalTimeMinus 	= 0;
	double MaximumBU 		= fMaximalBU;
	double MinimumBU  		= 0 ;
	double TheFinalTime 		= BurnupToSecond((MaximumBU-MinimumBU)/2.);
	double OldFinalTimePlus 	= BurnupToSecond(MaximumBU);
	double k_av 			= 0; //average kinf
	double OldPredictedk_av 	= 0;
	
	CLASSReader * reader = new CLASSReader( fMapOfTMVAVariableNames );
	reader->AddVariable( "Time" );
	reader->BookMVA( "MLP method" , fTMVAWeightPath[0] );
	
	for(int b = 0;b<fNumberOfBatch;b++)
	{
		float TheTime = (b+1)*TheFinalTime/fNumberOfBatch;

		TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
		reader->SetInputData( InputTree );
		
		OldPredictedk_av += reader->EvaluateRegression( "MLP method" )[0];
		
		delete InputTree;
	}
	OldPredictedk_av /= fNumberOfBatch;
	
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
				{ delete reader; return MaximumBU; }
		}
		
		else if( (OldPredictedk_av-fKThreshold)  < 0)//The burnup is too high
		{
			OldFinalTimePlus = TheFinalTime;
			TheFinalTime = TheFinalTime - fabs(OldFinalTimeMinus-TheFinalTime)/2.;
			if( SecondToBurnup(TheFinalTime) < (MaximumBU-MinimumBU)/2.*GetBurnUpPrecision() )
				{ delete reader; return 0; }
		}
		
		k_av = 0;
		for(int b = 0;b<fNumberOfBatch;b++)
		{
			float TheTime = (b+1)*TheFinalTime/fNumberOfBatch;
			TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
			reader->SetInputData( InputTree );
			
			k_av += reader->EvaluateRegression("MLP method")[0];
			
			delete InputTree;
		}
		k_av/= fNumberOfBatch;
		//cout<<SecondToBurnup(TheFinalTime)<<" ";
		OldPredictedk_av = k_av;
		count++;
//std::clog << "-> " << k_av << "\t\t(" << count << ") \t [" << TheFinalTime << "]" << "\t" << OldPredictedk_av-fKThreshold << "\t" << GetPCMPrecision() << std::endl; 
	}	while( fabs(OldPredictedk_av-fKThreshold) > GetPCMPrecision() )  ;
	
	delete reader;
	//cout<<endl;
	return SecondToBurnup(TheFinalTime);
}
