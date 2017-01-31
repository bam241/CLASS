#include "EquivalenceModel.hxx"
#include "EQM_MLP_PWR_MOxEUS.hxx"
#include "CLASSLogger.hxx"
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

EQM_MLP_PWR_MOxEUS::EQM_MLP_PWR_MOxEUS(string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold):EquivalenceModel(new CLASSLogger("EQM_MLP_PWR_MOxEUS.log"))
{
	fNumberOfBatch 		= NumOfBatch;
	fKThreshold 			= CriticalityThreshold ;

	fTMVAWeightPath.push_back(TMVAWeightPath);

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

	fStreamListEqMMassFractionMin["PuList"]		= 0.0;				
	fStreamListEqMMassFractionMin["EnrichmentList"]	= 0.0025;					

	fStreamListEqMMassFractionMax["PuList"]		= 0.16; 
	fStreamListEqMMassFractionMax["EnrichmentList"]	= 0.05;

	fSpecificPower 				= 34.24;
	fMaximalBU 				= 75;

	SetBurnUpPrecision(0.008);//1 % of the targeted burnup
	SetPCMPrecision(10);

	INFO<<"__An equivalence model of PWR MOxEUS has been define__"<<endl;
	INFO<<"\tThis model is based on a multi layer perceptron"<<endl;
	INFO<<"\t\tThe TMVA weight file is :"<<endl;
	INFO<<"\t\t\t"<<fTMVAWeightPath[0]<<endl;

}

//________________________________________________________________________
EQM_MLP_PWR_MOxEUS::EQM_MLP_PWR_MOxEUS(CLASSLogger* log, string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold):EquivalenceModel(log)
{
	fNumberOfBatch 		= NumOfBatch;
	fKThreshold 			= CriticalityThreshold ;

	fTMVAWeightPath.push_back(TMVAWeightPath);

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

	fStreamListEqMMassFractionMin["PuList"]		= 0.0;				
	fStreamListEqMMassFractionMin["EnrichmentList"]	= 0.0025;					

	fStreamListEqMMassFractionMax["PuList"]		= 0.16; 
	fStreamListEqMMassFractionMax["EnrichmentList"]	= 0.05;

	fSpecificPower 				= 34.24;
	fMaximalBU 				= 75;

	SetBurnUpPrecision(0.008);//1 % of the targeted burnup
	SetPCMPrecision(10);

	INFO<<"__An equivalence model of PWR MOxEUS has been define__"<<endl;
	INFO<<"\tThis model is based on a multi layer perceptron"<<endl;
	INFO<<"\t\tThe TMVA weight file is :"<<endl;
	INFO<<"\t\t\t"<<fTMVAWeightPath[0]<<endl;
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
double EQM_MLP_PWR_MOxEUS::CalculateTargetParameter(IsotopicVector TheFuel)
{
	/**************************************************************************/
	//With a dichotomy, the maximal irradiation time (TheFinalTime) is calculated
	//When average Kinf is very close (according "Precision") to the threshold
	//then the corresponding irradiation time is convert in burnup and returned
	/**************************************************************************/
	//Algorithm initialization


	double OldFinalTimeMinus 	= 0;
	double MinimumBU 		= 0; 	
	double MaximumBU 		= fMaximalBU; 
	double TheFinalTime 		= BurnupToSecond((MaximumBU-MinimumBU)/2.); 	
	double OldFinalTimePlus 	= BurnupToSecond(MaximumBU);
	double k_av 			= 0; //average kinf
	double OldPredictedk_av 	= 0;
 
	CLASSReader * reader = new CLASSReader();
	reader->AddVariable( "U5"   );
	reader->AddVariable( "U8"   );
	reader->AddVariable( "Pu8"  );
	reader->AddVariable( "Pu9"  );
	reader->AddVariable( "Pu10" );
	reader->AddVariable( "Pu11" );
	reader->AddVariable( "Pu12" );
	reader->AddVariable( "Am1"  );
	reader->AddVariable( "Time" );
	reader->BookMVA( "MLP method" , fTMVAWeightPath[0] );
	
 	for(int b=0;b<fNumberOfBatch;b++)
 	{
 		float TheTime 		= (b+1)*TheFinalTime/fNumberOfBatch;
 		TTree* InputTree 	= CreateTMVAInputTree(TheFuel,TheTime);
 		reader->SetInputData( InputTree );
 		
 		OldPredictedk_av   += reader->EvaluateRegression( "MLP method" )[0];
 		
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
 				{ delete reader; return MaximumBU; }
 		}

 		else if( (OldPredictedk_av-fKThreshold)  < 0)//The burnup is too high
 		{	
 			OldFinalTimePlus = TheFinalTime;
 			TheFinalTime = TheFinalTime - fabs(OldFinalTimeMinus - TheFinalTime)/2.;	
 			if( SecondToBurnup(TheFinalTime) < (MaximumBU-MinimumBU)/2.*GetBurnUpPrecision() )
 				{ delete reader; return 0; }
 		}
 		
 		k_av = 0;
 		for(int b=0;b<fNumberOfBatch;b++)
 		{
 			float TheTime = (b+1)*TheFinalTime/fNumberOfBatch;
 			TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
 			reader->SetInputData( InputTree );
 			
 			k_av +=	reader->EvaluateRegression("MLP method")[0];
 			delete InputTree;
 		}
 		k_av/=fNumberOfBatch;	

		OldPredictedk_av = k_av;
 		count++;
 	}	while( fabs(OldPredictedk_av-fKThreshold) > GetPCMPrecision() )  ;


	delete reader;
 	return SecondToBurnup(TheFinalTime);
}

