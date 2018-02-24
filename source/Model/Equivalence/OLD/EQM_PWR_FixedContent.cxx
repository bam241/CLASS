#include "EQM_PWR_FixedContent.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "CLASSReader.hxx"
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



//________________________________________________________________________
EQM_PWR_FixedContent::EQM_PWR_FixedContent(string TMVAWeightPath,  int NumOfBatch, string InformationFile, double CriticalityThreshold):EquivalenceModel(new CLASSLogger("EQM_PWR_FixedContent.log"))
{
	/**The information file and tmva weight*/
	fTMVAWeightPath.push_back(TMVAWeightPath);
	
	/* INFORMATION FILE HANDLING */
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");
	
	fMaximalBU 		= 0;
	fMaximalContent 	= 0;
	fNumberOfBatch 		= NumOfBatch;
	fKThreshold 		= CriticalityThreshold ;

	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetPCMPrecision(10);
	
	fInformationFile 	= InformationFile;

	LoadKeyword();
	ReadNFO();//Getting information from fInformationFile
	
	INFO << "__An equivalence model has been define__" 		<< endl;
	INFO << "\tThis model is based on the prediction of kinf" 	<< endl;
	INFO << "\tThe TMVA (weight | information) files are :" 		<< endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " 		<< fInformationFile << " )" << endl;
	INFO << "Maximal fissile content (molar proportion) : " 		<< fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " 			<< fMaximalBU << endl;
	
	EquivalenceModel::PrintInfo();

	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	bool EmptyList 		= false;
	bool FirstContentNULL 	= false;
	bool FixedContentNULL = false;

	for(  it_s_IV = fStreamList.begin();   it_s_IV != fStreamList.end();  it_s_IV++)
	{	
		if(fStreamList[(* it_s_IV).first].GetIsotopicQuantity().empty()){EmptyList=true;}
	}

	for(  it_s_D = fFirstGuessContent .begin();   it_s_D != fFirstGuessContent.end();  it_s_D++)
	{	
		if(fFirstGuessContent[(* it_s_D).first] ==0){FirstContentNULL=true;}
	}

	for(  it_s_D = fFixedMassContent .begin();   it_s_D != fFixedMassContent.end();  it_s_D++)
	{	
		if(fFixedMassContent[(* it_s_D).first] ==0){FixedContentNULL=true;}
	}

	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || FirstContentNULL==true || FixedContentNULL==true || fMaximalBU == 0 || fMaximalContent == 0 )	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}

}
//________________________________________________________________________
EQM_PWR_FixedContent::EQM_PWR_FixedContent(CLASSLogger* log, string TMVAWeightPath,  int NumOfBatch, string InformationFile, double CriticalityThreshold):EquivalenceModel(log)
{
	
	/**The information file and tmva weight*/
	fTMVAWeightPath.push_back(TMVAWeightPath);
	
	if(InformationFile == "")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");
	
	fMaximalBU 		= 0;
	fMaximalContent 	= 0;
	fNumberOfBatch 	= NumOfBatch;
	fKThreshold 		= CriticalityThreshold;

	SetBurnUpPrecision(0.005);//1 % of the targeted burnup
	SetPCMPrecision(10);

	fInformationFile 	= InformationFile;
	
	LoadKeyword();
	ReadNFO();//Getting information from fMLPInformationFile
	
	INFO << "__An equivalence model has been define__" 		<< endl;
	INFO << "\tThis model is based on the prediction of kinf" 	<< endl;
	INFO << "\tThe TMVA (weight | information) files are :" 		<< endl;
	INFO << "\t" << "( " << fTMVAWeightPath[0] << " | " 		<< fInformationFile << " )" << endl;
	INFO << "Maximal fissile content (molar proportion) : " 		<< fMaximalContent << endl;
	INFO << "Maximal burnup (GWd/tHM) : " 			<<  fMaximalBU << endl;

	EquivalenceModel::PrintInfo();

	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	
	bool EmptyList 		= false;
	bool FirstContentNULL 	= false;
	bool FixedContentNULL = false;

	for(  it_s_IV = fStreamList.begin();   it_s_IV != fStreamList.end();  it_s_IV++)
	{	
		if(fStreamList[(* it_s_IV).first].GetIsotopicQuantity().empty()){EmptyList=true;}
	}

	for(  it_s_D = fFirstGuessContent .begin();   it_s_D != fFirstGuessContent.end();  it_s_D++)
	{	
		if(fFirstGuessContent[(* it_s_D).first] ==0){FirstContentNULL=true;}
	}

	for(  it_s_D = fFixedMassContent .begin();   it_s_D != fFixedMassContent.end();  it_s_D++)
	{	
		if(fFixedMassContent[(* it_s_D).first] ==0){FixedContentNULL=true;}
	}

	if(fMapOfTMVAVariableNames.empty() || EmptyList==true || FirstContentNULL==true || FixedContentNULL==true || fMaximalBU == 0 || fMaximalContent == 0 )	
	{
		ERROR<<"Missing information file in : "<<fInformationFile<<endl;
		exit(1);
	}
}
//________________________________________________________________________
void EQM_PWR_FixedContent::LoadKeyword()
{
	DBGL
	
	fDKeyword.insert( pair<string, PWR_Fixed_DMthPtr>( "k_zainame"		, &EQM_PWR_FixedContent::ReadZAIName) );
	fDKeyword.insert( pair<string, PWR_Fixed_DMthPtr>( "k_maxburnup"		, &EQM_PWR_FixedContent::ReadMaxBurnUp) );
	fDKeyword.insert( pair<string, PWR_Fixed_DMthPtr>( "k_maxfiscontent"	, &EQM_PWR_FixedContent::ReadMaxFisContent) );
	fDKeyword.insert( pair<string, PWR_Fixed_DMthPtr>( "k_fixedmasscontent"	, &EQM_PWR_FixedContent::ReadFixedMassContent) );
	DBGL
}
//________________________________________________________________________
void EQM_PWR_FixedContent::ReadZAIName(const string &line)
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
void EQM_PWR_FixedContent::ReadMaxBurnUp(const string &line)
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
void EQM_PWR_FixedContent::ReadMaxFisContent(const string &line)
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
void EQM_PWR_FixedContent::ReadLine(string line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	
	map<string, PWR_Fixed_DMthPtr>::iterator it = fDKeyword.find(keyword);
	
	if(it != fDKeyword.end())
		(this->*(it->second))( line );
	
	DBGL
}
//________________________________________________________________________
//________________________________________________________________________
void EQM_PWR_FixedContent::ReadFixedMassContent(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_fixedmasscontent" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_fixedmasscontent\" not found !" << endl;
		exit(1);
	}
	string ListName	= StringLine::NextWord(line, pos, ' ');
	double Q 	= atof(StringLine::NextWord(line, pos, ' ').c_str());
	fFixedMassContent[ListName] = Q;
	
	DBGL
}

//________________________________________________________________________
map <string , vector<double> > EQM_PWR_FixedContent::BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray)
{
	map <string , vector<double> > lambda ; // map containing name of the list and associated vector of proportions taken from stocks
	
	//Iterators declaration
	map < string , vector  <IsotopicVector> >::iterator it_s_vIV;
	map < string , vector  <double> >::iterator it_s_vD;
	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	map < string , bool >::iterator it_s_B;
	
	map <string , IsotopicVector > IVStream;
	map <string , double > MaterialMolarContent;
	map <string , double > MaterialMassContent;
	map <string , double > MaterialMassNeeded;
	map <string , double > LambdaNeeded;
	map <string , double > DeltaMass;	


	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{	
		for(size_t i=0; i<StreamArray[(*it_s_vIV).first].size(); i++)
		{
			lambda[(*it_s_vIV).first].push_back(0);
		}
	}	

	/*** Test if there is stocks **/
	bool BreakReturnLambda = false; 
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		if(StreamArray[(*it_s_vIV).first].size() == 0)
		{
			WARNING << " No stock available for stream : "<< (*it_s_vIV).first <<".  Fuel not build." << endl;
			SetLambdaToErrorCode(lambda[(*it_s_vIV).first]);
			BreakReturnLambda = true; 	
		}
	}
	if(BreakReturnLambda) { return lambda;}
	HMMass *=  1e6; //Unit onversion : tons to gram
	
	/**** Some initializations **/
	StocksTotalMassCalculation(StreamArray);

	MaterialMassContent 		= GetFixedMassContent();	
	bool FuelBuiltCorrectly	 	= false ; 
	
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		LambdaNeeded[(*it_s_vIV).first] = 0;

	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		MaterialMassNeeded[(*it_s_vIV).first] = HMMass*MaterialMassContent[(*it_s_vIV).first]; 
		DeltaMass[(*it_s_vIV).first] =  fTotalMassInStocks[(* it_s_vIV).first] - MaterialMassNeeded[(*it_s_vIV).first];
	}

	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		LambdaNeeded[(*it_s_vIV).first] = LambdaCalculation((*it_s_vIV).first, 0., MaterialMassNeeded[(*it_s_vIV).first], DeltaMass[(*it_s_vIV).first], StreamArray[(*it_s_vIV).first]);

	for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
		SetLambda(lambda[(*it_s_D).first], LambdaNeeded[(*it_s_D).first]); 

	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{	
		for(size_t i=0; i<StreamArray[(*it_s_vIV).first].size(); i++)
			IVStream[(*it_s_vIV).first]  +=  lambda[(*it_s_vIV).first][i] * StreamArray[(*it_s_vIV).first][i];	
	}

	IsotopicVector FreshFuel;

	for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
		FreshFuel +=  MaterialMassContent[(*it_s_IV).first]*(IVStream[(*it_s_IV).first]/IVStream[(*it_s_IV).first].GetSumOfAll());
	
	double PredictedBU =  GetMaximumBurnUp(FreshFuel,BurnUp); 
	
	if(PredictedBU<BurnUp)
	{
		WARNING << " Fixed contents are not adequate. Predicted BU is lower than target BU. "<< endl;
		WARNING << " Predicted BU :  "<<PredictedBU<<"  Target BU:  "<<BurnUp<< endl;
	}

	if(PredictedBU>BurnUp)
	{
		WARNING << " Fixed contents are not adequate. Predicted BU is greater than target BU. "<< endl;
		WARNING << " Predicted BU :  "<<PredictedBU<<"  Target BU:  "<<BurnUp<< endl;
	}

	return lambda;
}

TTree* EQM_PWR_FixedContent::CreateTMVAInputTree(IsotopicVector TheFreshfuel, double ThisTime)
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
double EQM_PWR_FixedContent::ExecuteTMVA(TTree* InputTree,string WeightPath)
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
	
	InputTree->SetBranchAddress( "Time" ,&Time );
	
	InputTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];
	
	delete reader;
	
	return (double)val;//retourn k_{inf}(t = Time)
}
//________________________________________________________________________
double EQM_PWR_FixedContent::GetMaximumBurnUp(IsotopicVector TheFuel, double TargetBU)
{
	/**************************************************************************/
	//With a dichotomy, the maximal irradiation time (TheFinalTime) is calculated
	//When average Kinf is very close (according "Precision") to the threshold
	//then the corresponding irradiation time is convert in burnup and returned
	/**************************************************************************/
	//Algorithm initialization
	double TheFinalTime      = BurnupToSecond(TargetBU);
	double OldFinalTimeMinus = 0;
	double MaximumBU         = fMaximalBU;
	double OldFinalTimePlus  = BurnupToSecond(MaximumBU); 
	double k_av              = 0; //average kinf
	double OldPredictedk_av  = 0;
	
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
	OldPredictedk_av/= fNumberOfBatch;

	//Algorithm control
	int count = 0;
	int MaximumLoopCount = 500;
	do
	{
		if(count > MaximumLoopCount )
		{
			ERROR << "CRITICAL ! Can't manage to predict burnup\nHint : Try to increase the precision on k effective using :\n YourEQM_PWR_FixedContent->SetPCMPrecision(pcm); with pcm the precision in pcm (default 10) REDUCE IT\n If this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
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
			if( SecondToBurnup(TheFinalTime) < TargetBU*GetBurnUpPrecision() )
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
		
		OldPredictedk_av = k_av;
		count++;
		
	}	while( fabs(OldPredictedk_av-fKThreshold) > GetPCMPrecision() )  ;
	
	delete reader;
	return SecondToBurnup(TheFinalTime);
}
//________________________________________________________________________
map < string , double>  EQM_PWR_FixedContent::GetMolarFraction(map <string , IsotopicVector > IVStream, double TargetBU)
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
			ERROR << "CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on burnup using :\nYourEQM_PWR_FixedContent->SetBurnUpPrecision(prop); with prop the precision  (default 0.5percent :  0.005) INCREASE IT\nIf this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
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


