#include "EquivalenceModel.hxx"
#include "EQM_MLP_PWR_MOxEUS.hxx"
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
//		EQM_MLP_PWR_MOxEUS
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For PWR MOxEUS use
//		-> MOxEUS building
//			- if PuMassContent calculated < Maximal Pu value => MOx fuel
//			- if PuMassContent calculated > Maximal Pu value => Uranium enrichment adjusted to reach BU
//
//________________________________________________________________________

EQM_MLP_PWR_MOxEUS::EQM_MLP_PWR_MOxEUS(string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold, double MaximalPuMassContent):EquivalenceModel(new CLASSLogger("EQM_MLP_PWR_MOxEUS.log"))
{
	fNumberOfBatch 			= NumOfBatch;
	fKThreshold 			= CriticalityThreshold ;
	fMaximalPuMassContent 	= MaximalPuMassContent;

	fTMVAWeightPath.push_back(TMVAWeightPath);

	SetMinimalU5Enrichment(0.003);

	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);

	fStreamList["PuList"] 			= Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1;
	fStreamList["EnrichmentList"] 		= U5*1;
	fStreamList["FertileList"] 		= U8*1;

	fFirstGuessContent["PuList"] 		=  0.02;
	fFirstGuessContent["EnrichmentList"]     =  (1-fFirstGuessContent["PuList"])*GetMinimalU5Enrichment();
	fFirstGuessContent["FertileList"] 	=  1 - fFirstGuessContent["PuList"] - fFirstGuessContent["EnrichmentList"];

	fSpecificPower 				= 34.24;
	fMaximalPuMolarContent 		= 0.50;
	fMaximalBU 				= 100;

	SetRelativMassPrecision(0.0008);
	SetMaximalU5Enrichment(0.10);
	SetBurnUpPrecision(0.008);//1 % of the targeted burnup
	SetPCMPrecision(10);

	INFO<<"__An equivalence model of PWR MOxEUS has been define__"<<endl;
	INFO<<"\tThis model is based on a multi layer perceptron"<<endl;
	INFO<<"\t\tThe TMVA weight file is :"<<endl;
	INFO<<"\t\t\t"<<fTMVAWeightPath[0]<<endl;

}

//________________________________________________________________________
EQM_MLP_PWR_MOxEUS::EQM_MLP_PWR_MOxEUS(CLASSLogger* log, string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold, double MaximalPuMassContent):EquivalenceModel(log)
{
	fNumberOfBatch 		= NumOfBatch;
	fKThreshold 			= CriticalityThreshold ;
	fMaximalPuMassContent 	= MaximalPuMassContent;

	fTMVAWeightPath.push_back(TMVAWeightPath);

	SetMinimalU5Enrichment(0.003);

	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);

	fStreamList["PuList"] 			= Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1;
	fStreamList["EnrichmentList"] 		= U5*1;
	fStreamList["FertileList"] 		= U8*1;

	fFirstGuessContent["PuList"] 		=  0.02;
	fFirstGuessContent["EnrichmentList"]     =  (1-fFirstGuessContent["PuList"])*GetMinimalU5Enrichment();
	fFirstGuessContent["FertileList"] 	=  1 - fFirstGuessContent["PuList"] - fFirstGuessContent["EnrichmentList"];

	fSpecificPower 				= 34.24;
	fMaximalPuMolarContent 		= 1.0;
	fMaximalBU 				= 100;

	SetRelativMassPrecision(0.0008);
	SetMaximalU5Enrichment(0.10);
	SetBurnUpPrecision(0.008);//1 % of the targeted burnup
	SetPCMPrecision(10);

	INFO<<"__An equivalence model of PWR MOxEUS has been define__"<<endl;
	INFO<<"\tThis model is based on a multi layer perceptron"<<endl;
	INFO<<"\t\tThe TMVA weight file is :"<<endl;
	INFO<<"\t\t\t"<<fTMVAWeightPath[0]<<endl;

}

//________________________________________________________________________
map <string , vector<double> > EQM_MLP_PWR_MOxEUS::BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray)
{

	map <string , vector<double> > lambda ; // map containing name of the list and associated vector of proportions taken from stocks
	
	//Iterators declaration
	map < string , vector  <IsotopicVector> >::iterator it_s_vIV;
	map < string , vector  <double> >::iterator it_s_vD;
	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	map < string , bool >::iterator it_s_B;
	
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{	
		for(int i=0; i<StreamArray[(*it_s_vIV).first].size(); i++)
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

	fActualMassContentInFuel 	= GetBuildFuelFirstGuess();	
	int loopCount 			= 0;
	bool FuelBuiltCorrectly	 	= false ; 

	map <string , IsotopicVector > IVStream;
	map <string , double > MaterialMolarContent;
	map <string , double > MaterialMassContent;
	map <string , double > MaterialMassNeeded;
	map <string , double > LambdaNeeded;
	map <string , double > DeltaMass;	
	
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		LambdaNeeded[(*it_s_vIV).first] = 0;
		DeltaMass[(*it_s_vIV).first] = 0;
	}

	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		MaterialMassNeeded[(*it_s_vIV).first] = HMMass*fActualMassContentInFuel[(*it_s_vIV).first]; 
		DeltaMass[(*it_s_vIV).first] =  - MaterialMassNeeded[(*it_s_vIV).first];
	}

	do
	{	
		map < string , double >  MaterialMeanMolar;
		map < string , double >  MaterialMassAvailable;
		map < string , bool > 	  CheckOnMass;

		map < string , double >  LambdaPreviousStep;

		double FuelMeanMolar   = 0;

		bool NotEnoughPuInStock 	= false;
		BreakReturnLambda 		= false;

		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
			LambdaPreviousStep[(*it_s_D ).first] =LambdaNeeded[(*it_s_D ).first]; 
		
		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
			LambdaNeeded[(*it_s_vIV).first] = LambdaCalculation((*it_s_vIV).first, LambdaPreviousStep[(*it_s_vIV).first], MaterialMassNeeded[(*it_s_vIV).first], DeltaMass[(*it_s_vIV).first], StreamArray[(*it_s_vIV).first]);


		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
		{		
			if( LambdaNeeded[(*it_s_D).first] == -1 )
			{
				SetLambdaToErrorCode(lambda[(*it_s_D).first]);
				WARNING << "Not enough : "<< (*it_s_D).first <<" to build fuel." << endl;	
				BreakReturnLambda = true; 
				if((*it_s_D).first=="PuList")	
					NotEnoughPuInStock = true;
			}
		}
		
		if(BreakReturnLambda && !NotEnoughPuInStock){ return lambda;}
		
		if(BreakReturnLambda && NotEnoughPuInStock)
		{ 
			MaterialMolarContent["PuList"] = -1;
			break;
		}
		
		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
			SetLambda(lambda[(*it_s_D).first], LambdaNeeded[(*it_s_D).first]); 
	
		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
			IVStream[(*it_s_vIV).first].Clear();

		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		{	
			for(int i=0; i<StreamArray[(*it_s_vIV).first].size(); i++)
				IVStream[(*it_s_vIV).first]  +=  lambda[(*it_s_vIV).first][i] * StreamArray[(*it_s_vIV).first][i];	
		}

		if (loopCount > fMaxInterration)
		{
			ERROR << "Too much iterration in BuildFuel Method !";
			ERROR << "Need improvement in fuel fabrication ! Ask for it or D.I.Y. !!" << endl;
			exit(1);
		}
		
		/* Calcul the quantity of this composition needed to reach the burnup */
		MaterialMolarContent = GetMolarFraction(IVStream, BurnUp);
		//cout<<MaterialMolarContent["PuList"]<<endl;
		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
		{	
			if( MaterialMolarContent[(*it_s_D).first] < 0 || MaterialMolarContent[(*it_s_D).first] > 1 )
			{
				SetLambdaToErrorCode(lambda[(*it_s_D).first]);
				WARNING << "GetMolarFraction return negative or greater than one value, at least for one material : "<<(*it_s_D).first;
				return lambda;
			}
		}
		for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
			MaterialMassAvailable[(*it_s_IV).first] = IVStream[(*it_s_IV).first].GetTotalMass()*1e06; 

		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
		{
			MaterialMeanMolar[(*it_s_D).first] = IVStream[(*it_s_D).first].GetMeanMolarMass();
			FuelMeanMolar +=  MaterialMolarContent[(*it_s_D).first] * MaterialMeanMolar[(*it_s_D).first];
		}

		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
			MaterialMassContent [(*it_s_D).first] =  MaterialMolarContent[(*it_s_D).first] * MaterialMeanMolar[(*it_s_D).first] / FuelMeanMolar; 

		fActualMolarContentInFuel  = MaterialMolarContent; //fActualMolarContentInFuel is used in GetMolarFraction
		
		for( it_s_D = MaterialMassContent.begin();  it_s_D != MaterialMassContent.end(); it_s_D++)
		{	
			MaterialMassNeeded[(*it_s_D).first] = HMMass*MaterialMassContent[(*it_s_D).first]; 
			DeltaMass[(*it_s_D).first] = MaterialMassAvailable[(*it_s_D).first] - MaterialMassNeeded[(*it_s_D).first];
		}

		for( it_s_D = MaterialMassNeeded.begin();  it_s_D != MaterialMassNeeded.end(); it_s_D++)
		{
			double DeltaM = fabs(MaterialMassNeeded[(*it_s_D).first] - MaterialMassAvailable[(*it_s_D).first]) / HMMass ; 
			if(DeltaM<fRelativMassPrecision) {CheckOnMass[(*it_s_D).first] = true;}
			else{CheckOnMass[(*it_s_D).first] = false;}
		}

		FuelBuiltCorrectly = true; 
		for( it_s_B = CheckOnMass.begin();  it_s_B != CheckOnMass.end(); it_s_B++)
		{
			if(!CheckOnMass[(*it_s_B).first]) {FuelBuiltCorrectly = false;}
		}
		
		loopCount++;

	}while(!FuelBuiltCorrectly);

	if(MaterialMassContent["PuList"]>fMaximalPuMassContent || MaterialMolarContent["PuList"]==-1)
	{
		for( it_s_vD = lambda.begin();  it_s_vD != lambda.end(); it_s_vD++)
		{
			for(int i=0; i<lambda[(*it_s_vD).first].size(); i++)
			{
				lambda[(*it_s_vD).first][i] = 0;
			}
		}

		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
			LambdaNeeded[(*it_s_D).first] = 0;

		for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
			IVStream[(*it_s_IV).first] = IsotopicVector();

		if(MaterialMassContent["PuList"]>fMaximalPuMassContent)
		{
			MaterialMolarContent["PuList"] = 0;
			MaterialMassNeeded["PuList"]  =  fMaximalPuMassContent*HMMass;
			LambdaNeeded["PuList"]          =  LambdaCalculation("PuList", 0., MaterialMassNeeded["PuList"], 1., StreamArray["PuList"]);
		}
		
		if(MaterialMolarContent["PuList"]==-1)
		{
			LambdaNeeded["PuList"]	 =  LambdaCalculation("PuList", 0., fTotalMassInStocks["PuList"], 1., StreamArray["PuList"]);
		}

		SetLambda(lambda["PuList"], LambdaNeeded["PuList"]);

		for( int i = 0 ; i < (int)lambda["PuList"].size() ; i++ )
			IVStream["PuList"] +=lambda["PuList"][i] * StreamArray["PuList"][i];
		
		fActualMassContentInFuel["PuList"] = (IVStream["PuList"].GetTotalMass()*1e6)/HMMass;

		LambdaNeeded["EnrichmentList"]	=  LambdaCalculation("EnrichmentList", 0., HMMass, 1., StreamArray["EnrichmentList"]);
		LambdaNeeded["FertileList"]	  	=  LambdaCalculation("FertileList", 0., HMMass, 1., StreamArray["FertileList"]);

		SetLambda(lambda["EnrichmentList"], LambdaNeeded["EnrichmentList"]);
		SetLambda(lambda["FertileList"], LambdaNeeded["FertileList"]);
		
		for( int i = 0 ; i < (int)lambda["EnrichmentList"].size() ; i++ )
			IVStream["EnrichmentList"] +=lambda["EnrichmentList"][i] * StreamArray["EnrichmentList"][i];

		for( int i = 0 ; i < (int)lambda["FertileList"].size() ; i++ )
			IVStream["FertileList"] +=lambda["FertileList"][i] * StreamArray["FertileList"][i];

		double U5Enrichment			= GetU5Enrichment(IVStream, BurnUp);
		double MeanMolarUdepleted		= U5Enrichment*IVStream["EnrichmentList"].GetMeanMolarMass() + (1 - U5Enrichment)*IVStream["FertileList"].GetMeanMolarMass();

		double UdepletedMassNeeded 		= (HMMass - IVStream["PuList"].GetTotalMass() * 1e6);
		double EnrichedSupportMassNeeded	= U5Enrichment*UdepletedMassNeeded*IVStream["EnrichmentList"].GetMeanMolarMass()/IVStream["FertileList"].GetMeanMolarMass();
		double FertilMassNeeded 		= UdepletedMassNeeded - EnrichedSupportMassNeeded;
		
		LambdaNeeded["EnrichmentList"] 	= LambdaCalculation("EnrichmentList", 0., EnrichedSupportMassNeeded, 1., StreamArray["EnrichmentList"]);
		LambdaNeeded["FertileList"] 		= LambdaCalculation("FertileList", 0., FertilMassNeeded, 1., StreamArray["FertileList"]);
		
		BreakReturnLambda = false; 
		
		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
		{		
			if( LambdaNeeded[(*it_s_D).first] == -1 )
			{
				SetLambdaToErrorCode(lambda[(*it_s_D).first]);
				WARNING << "Not enough : "<< (*it_s_D).first <<" to build fuel." << endl;	
				BreakReturnLambda = true; 				
			}
		}
		if(BreakReturnLambda) { return lambda;}

		SetLambda(lambda["EnrichmentList"], LambdaNeeded["EnrichmentList"]);
		SetLambda(lambda["FertileList"], LambdaNeeded["FertileList"]);
		
		IVStream["EnrichmentList"]	= IsotopicVector();
		IVStream["FertileList"]		= IsotopicVector();

		for( int i = 0 ; i < (int)lambda["EnrichmentList"].size() ; i++ )
			IVStream["EnrichmentList"] +=lambda["EnrichmentList"][i] * StreamArray["EnrichmentList"][i];

		for( int i = 0 ; i < (int)lambda["FertileList"].size() ; i++ )
			IVStream["FertileList"] +=lambda["FertileList"][i] * StreamArray["FertileList"][i];
	}
	
	for( it_s_D = MaterialMassNeeded.begin();  it_s_D != MaterialMassNeeded.end(); it_s_D++)
		DBGV( "Weight percent : "<<(*it_s_D).first<<" "<< MaterialMassNeeded[(*it_s_D).first]/HMMass);
	
	for( it_s_vD = lambda.begin();  it_s_vD != lambda.end(); it_s_vD++)
	{	
		DBGV( "Lambda vector : "<<(*it_s_vD).first );

		for(int i=0; i<lambda[(*it_s_vD).first].size(); i++)
		{
			DBGV(lambda[(*it_s_vD).first][i]); 
		}
	}		

	return lambda;
}

//________________________________________________________________________
TTree* EQM_MLP_PWR_MOxEUS::CreateTMVAInputTree(IsotopicVector TheFuel, double ThisTime)
{
	TTree*   InputTree = new TTree("EQTMP", "EQTMP");

	float Pu8   	= 0;
	float Pu9   	= 0;
	float Pu10  	= 0;
	float Pu11  	= 0;
	float Pu12  	= 0;
	float Am1   	= 0;
	float U5	= 0;
	float U8	= 0;
	float Time 	= 0;

	InputTree->Branch(	"U5"	 ,&U5	  ,"U5/F"	);
	InputTree->Branch(	"U8"	 ,&U8	  ,"U8/F"	);
	InputTree->Branch(	"Pu8"	 ,&Pu8	  ,"Pu8/F"	);
	InputTree->Branch(	"Pu9"	 ,&Pu9	  ,"Pu9/F"	);
	InputTree->Branch(	"Pu10"	 ,&Pu10  ,"Pu10/F"	);
	InputTree->Branch(	"Pu11"	 ,&Pu11  ,"Pu11/F"	);
	InputTree->Branch(	"Pu12"	 ,&Pu12  ,"Pu12/F"	);
	InputTree->Branch(	"Am1"	 ,&Am1	  ,"Am1/F"	);
	InputTree->Branch(	"Time"	 ,&Time  ,"Time/F"	);

	double Ntot = TheFuel.GetSumOfAll();

	U5	= TheFuel.GetZAIIsotopicQuantity(92,235,0)/Ntot;
	U8	= TheFuel.GetZAIIsotopicQuantity(92,238,0)/Ntot;
	Pu8    	= TheFuel.GetZAIIsotopicQuantity(94,238,0)/Ntot;
	Pu9    	= TheFuel.GetZAIIsotopicQuantity(94,239,0)/Ntot;
	Pu10   	= TheFuel.GetZAIIsotopicQuantity(94,240,0)/Ntot;
	Pu11   	= TheFuel.GetZAIIsotopicQuantity(94,241,0)/Ntot;
	Pu12   	= TheFuel.GetZAIIsotopicQuantity(94,242,0)/Ntot;
	Am1    	= TheFuel.GetZAIIsotopicQuantity(95,241,0)/Ntot;

	Time=ThisTime;

	InputTree->Fill();
	return InputTree;
}
//________________________________________________________________________
double EQM_MLP_PWR_MOxEUS::ExecuteTMVA(TTree* theTree, string WeightPath)
{
	// --- Create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "Silent" );
	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	Float_t U5, U8, Pu8,Pu9,Pu10,Pu11,Pu12,Am1, Time;

	reader->AddVariable( "U5"  		,&U5	);
	reader->AddVariable( "U8"  		,&U8	);
	reader->AddVariable( "Pu8"  		,&Pu8 	);
	reader->AddVariable( "Pu9"  		,&Pu9 	);
	reader->AddVariable( "Pu10" 		,&Pu10	);
	reader->AddVariable( "Pu11" 		,&Pu11	);
	reader->AddVariable( "Pu12" 		,&Pu12	);
	reader->AddVariable( "Am1"  		,&Am1	);
	reader->AddVariable( "Time"  		,&Time	);

	// --- Book the MVA methods

	// Book method MLP
	TString methodName = "MLP method";
	reader->BookMVA( methodName, fTMVAWeightPath[0] );
	theTree->SetBranchAddress( "U5"  		,&U5	 );
	theTree->SetBranchAddress( "U8"  		,&U8	 );
	theTree->SetBranchAddress( "Pu8"  		,&Pu8   	 );
	theTree->SetBranchAddress( "Pu9"  		,&Pu9   	 );
	theTree->SetBranchAddress( "Pu10" 		,&Pu10 	 );
	theTree->SetBranchAddress( "Pu11" 		,&Pu11 	 );
	theTree->SetBranchAddress( "Pu12" 		,&Pu12 	 );
	theTree->SetBranchAddress( "Am1"  		,&Am1   );
	theTree->SetBranchAddress( "Time"  		,&Time  );

	theTree->GetEntry(0);

	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;
	return (double)val; //return k_eff
}

//________________________________________________________________________
double EQM_MLP_PWR_MOxEUS::GetMaximumBurnUp(IsotopicVector TheFuel, double TargetBU)
{
	/**************************************************************************/
	//With a dichotomy, the maximal irradiation time (TheFinalTime) is calculated
	//When average Kinf is very close (according "Precision") to the threshold
	//then the corresponding irradiation time is convert in burnup and returned
	/**************************************************************************/
	//Algorithm initialization

	double TheFinalTime 		= BurnupToSecond(TargetBU); 
	double OldFinalTimeMinus 	= 0;
	double MaximumBU 		= fMaximalBU; 
	double OldFinalTimePlus 	= BurnupToSecond(MaximumBU);
	double k_av 			= 0; //average kinf
	double OldPredictedk_av 	= 0;
 
 	for(int b=0;b<fNumberOfBatch;b++)
 	{
 		float TheTime 		= (b+1)*TheFinalTime/fNumberOfBatch;
 		TTree* InputTree 	= CreateTMVAInputTree(TheFuel,TheTime);
 		OldPredictedk_av        += ExecuteTMVA(InputTree,fTMVAWeightPath[0]);
 		delete InputTree;
 	}
 	OldPredictedk_av/=fNumberOfBatch;	
 
 	//Algorithm control
	int count=0;
	int MaximumLoopCount = 500;
	do
	{	
		if(count > MaximumLoopCount )
		{
			ERROR<<"CRITICAL ! Can't manage to predict burnup\nHint : Try to increase the precision on k effective using :\n YourEQM_MLP_Kinf->SetPCMprecision(pcm); with pcm the precision in pcm (default 10) REDUCE IT\n If this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr "<<endl;
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
 			TheFinalTime = TheFinalTime - fabs(OldFinalTimeMinus - TheFinalTime)/2.;	
 			if( SecondToBurnup(TheFinalTime) < TargetBU*GetBurnUpPrecision() )
 				return 0;
 		}
 		
 		k_av = 0;
 		for(int b=0;b<fNumberOfBatch;b++)
 		{
 			float TheTime = (b+1)*TheFinalTime/fNumberOfBatch;
 			TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
 			k_av +=	ExecuteTMVA(InputTree,fTMVAWeightPath[0]);
 			delete InputTree;
 		}
 		k_av/=fNumberOfBatch;	

		OldPredictedk_av = k_av;
 		count++;
 	}	while( fabs(OldPredictedk_av-fKThreshold) > GetPCMPrecision() )  ;

 	return SecondToBurnup(TheFinalTime);
}

//________________________________________________________________________
map < string , double>  EQM_MLP_PWR_MOxEUS::GetMolarFraction(map <string , IsotopicVector > IVStream, double TargetBU)
{

	IsotopicVector Pu 			= IVStream["PuList"];
	IsotopicVector EnrichedSupport 	= IVStream["EnrichmentList"];
	IsotopicVector Fertil 			= IVStream["FertileList"];

	map < string , double> MolarFraction ;

	double OldPuContentMinus	= 0;
	double OldPuContentPlus	= fMaximalPuMolarContent;
	
	double PredictedBU 		= 0 ;
	double PuContent 		= fActualMolarContentInFuel["PuList"];
	double U5Enrichment	 	= fMinimalU5Enrichment; 

	IsotopicVector FreshFuel	= PuContent*(Pu/Pu.GetSumOfAll()) +  U5Enrichment*(1-PuContent) *(EnrichedSupport/EnrichedSupport.GetSumOfAll()) +  (1 - U5Enrichment)*(1-PuContent)*(Fertil/Fertil.GetSumOfAll());
	double OldPredictedBU		=  GetMaximumBurnUp(FreshFuel,TargetBU);

	double Precision 		= GetBurnUpPrecision()*TargetBU; //1 % of the targeted burnup
	int count 				= 0;
	int MaximumLoopCount 	= 500;
	do	
	{
		if(count > MaximumLoopCount )
		{
			ERROR<<"CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on burnup using :\nYourEQM_MLP_Kinf->SetBurnUpPrecision(prop); with prop the precision  (default 0.5percent :  0.005) INCREASE IT\nIf this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr "<<endl;
			exit(1);
		}
		if( (OldPredictedBU - TargetBU) < 0 ) //The Content can be increased
		{
			OldPuContentMinus = PuContent;
			PuContent = PuContent + fabs(OldPuContentPlus-PuContent)/2.;
		}
		else if( (OldPredictedBU - TargetBU) > 0) //The Content is too high
		{
			OldPuContentPlus = PuContent;
			PuContent = PuContent - fabs(OldPuContentMinus-PuContent)/2.;
		}

		IsotopicVector FreshFuel= PuContent*(Pu/Pu.GetSumOfAll()) +  U5Enrichment*(1-PuContent) *(EnrichedSupport/EnrichedSupport.GetSumOfAll()) +  (1 - U5Enrichment)*(1-PuContent)*(Fertil/Fertil.GetSumOfAll());
		PredictedBU = GetMaximumBurnUp(FreshFuel,TargetBU);
	
		OldPredictedBU = PredictedBU;
		count ++;
INFO<<"count "<<count<<" => Target BU "<<TargetBU<<" - Predicted BU "<<PredictedBU<<endl;

	}while(fabs(TargetBU-PredictedBU)>Precision);

	DBGV("Predicted BU "<<PredictedBU<<" PuContent "<<PuContent);
	
	MolarFraction["PuList"] 			= PuContent;
	MolarFraction["EnrichmentList"] 	= (1 - MolarFraction["PuList"])*fMinimalU5Enrichment;
	MolarFraction["FertileList"] 		= 1 - MolarFraction["PuList"] - MolarFraction["EnrichmentList"];

return MolarFraction;
}
//________________________________________________________________________

//________________________________________________________________________
double EQM_MLP_PWR_MOxEUS::GetU5Enrichment(map <string , IsotopicVector > IVStream, double TargetBU)
{
	IsotopicVector Pu 		 = IVStream["PuList"];
	IsotopicVector EnrichedSupport  = IVStream["EnrichmentList"];
	IsotopicVector Fertil 		 = IVStream["FertileList"];
	
	double OldU5EnrichmentMinus  = fMinimalU5Enrichment;
	double OldU5EnrichmentPlus	= fMaximalU5Enrichment;
	double U5Enrichment	 	= fMinimalU5Enrichment; 

	double PredictedBU 		= 0 ;
	double USupportMassNeeded	= (Pu.GetTotalMass()*1e6)/fActualMassContentInFuel["PuList"] - (Pu.GetTotalMass()*1e6);

	IsotopicVector USupport	= U5Enrichment*EnrichedSupport + (1 - U5Enrichment)*Fertil;
	double USupportMassContent	= USupportMassNeeded/(USupport.GetTotalMass()*1e6);	

	IsotopicVector FreshFuel	= Pu +  USupportMassContent*USupport;
	FreshFuel			= FreshFuel/FreshFuel.GetSumOfAll();

	double OldPredictedBU		= GetMaximumBurnUp(FreshFuel,TargetBU);
	
	double Precision 		= GetBurnUpPrecision()*TargetBU; //1 % of the targeted burnup
	int count 			= 0;
	int MaximumLoopCount 	= 500;

	do	
	{

		if(count > MaximumLoopCount )
		{
			ERROR<<"CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on burnup using :\nYourEQM_MLP_Kinf->SetBurnUpPrecision(prop); with prop the precision  (default 0.5percent :  0.005) INCREASE IT\nIf this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr "<<endl;
			exit(1);
		}
		if( (OldPredictedBU - TargetBU) < 0 ) //The Content can be increased
		{
			OldU5EnrichmentMinus = U5Enrichment;
			U5Enrichment = U5Enrichment + fabs(OldU5EnrichmentPlus-U5Enrichment)/2.;
		}
		else if( (OldPredictedBU - TargetBU) > 0) //The Content is too high
		{
			OldU5EnrichmentPlus = U5Enrichment;
			U5Enrichment = U5Enrichment - fabs(OldU5EnrichmentMinus-U5Enrichment)/2.;
		}

		IsotopicVector USupport	= U5Enrichment*EnrichedSupport + (1 - U5Enrichment)*Fertil;

		double USupportMassContent	= USupportMassNeeded/(USupport.GetTotalMass()*1e6);	
		IsotopicVector FreshFuel	= Pu +  USupportMassContent*USupport;
		FreshFuel			= FreshFuel/FreshFuel.GetSumOfAll();
		
		PredictedBU = GetMaximumBurnUp(FreshFuel,TargetBU);
		
		OldPredictedBU = PredictedBU;
		count ++;
	}while(fabs(TargetBU-PredictedBU)>Precision);

DBGV("Predicted BU "<<PredictedBU<<" U5Enrichment "<<U5Enrichment);
return U5Enrichment;

}
