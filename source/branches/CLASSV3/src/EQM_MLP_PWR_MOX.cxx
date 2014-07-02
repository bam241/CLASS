#include "EquivalenceModel.hxx"
#include "EQM_MLP_PWR_MOX.hxx"
#include "LogFile.hxx"
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


#include "CLASSHeaders.hxx"

//________________________________________________________________________
//
//		EQM_MLP_MOX
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For REP MOX use
//
//________________________________________________________________________

EQM_MLP_MOX::EQM_MLP_MOX(string TMVAWeightPath)
{
	fTMVAWeightPath =  TMVAWeightPath;
	
	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	double U5_enrich= 0.0025;

	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);

	fFissileList = Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1;
	fFertileList = U5*U5_enrich + U8*(1-U5_enrich);
	
}
//________________________________________________________________________
TTree* EQM_MLP_MOX::CreateTMVAInputTree(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/
	TTree*   InputTree = new TTree("EQTMP", "EQTMP");
	float Pu8   			 = 0;
	float Pu9   			 = 0;
	float Pu10  			 = 0;
	float Pu11  			 = 0;
	float Pu12  			 = 0;
	float Am1   			 = 0;
	float U5_enrichment 	 = 0;
	float BU  			 	 = 0;

	InputTree->Branch(	"Pu8"	,&Pu8	,"Pu8/F"	);
	InputTree->Branch(	"Pu9"	,&Pu9	,"Pu9/F"	);
	InputTree->Branch(	"Pu10"	,&Pu10	,"Pu10/F"	);
	InputTree->Branch(	"Pu11"	,&Pu11	,"Pu11/F"	);
	InputTree->Branch(	"Pu12"	,&Pu12	,"Pu12/F"	);
	InputTree->Branch(	"Am1"	,&Am1	,"Am1/F"	);
	InputTree->Branch(	"U5_enrichment"	,&U5_enrichment	,"U5_enrichment/F"	);
	InputTree->Branch(	"BU"	,&BU	,"BU/F"	);


	float U8     = Fertil.GetZAIIsotopicQuantity(92,238,0);
	float U5     = Fertil.GetZAIIsotopicQuantity(92,235,0);
	float U4     = Fertil.GetZAIIsotopicQuantity(92,234,0);

	float UTOT = U8 + U5 + U4;

	Pu8    	   = Fissil.GetZAIIsotopicQuantity(94,238,0);
	Pu9    	   = Fissil.GetZAIIsotopicQuantity(94,239,0);
	Pu10   	   = Fissil.GetZAIIsotopicQuantity(94,240,0);
	Pu11   	   = Fissil.GetZAIIsotopicQuantity(94,241,0);
	Pu12   	   = Fissil.GetZAIIsotopicQuantity(94,242,0);
	Am1        = Fissil.GetZAIIsotopicQuantity(95,241,0);

	double TOTPU=(Pu8+Pu9+Pu10+Pu11+Pu12+Am1);

	Pu8 = Pu8  / TOTPU;
	Pu9 = Pu9  / TOTPU;
	Pu10= Pu10 / TOTPU;
	Pu11= Pu11 / TOTPU;
	Pu12= Pu12 / TOTPU;
	Am1 = Am1  / TOTPU;

	U5_enrichment = U5 / UTOT;

	BU=BurnUp;

	if(Pu8 + Pu9 + Pu10 + Pu11 + Pu12 + Am1 > 1.00001 )//?????1.00001??? I don't know it! goes in condition if =1 !! may be float/double issue ...
	{
		cout<<"!!!!!!!!!!!ERRORR!!!!!!!!!!!!"<<endl;
		cout<<Pu8<<" "<<Pu9<<" "<<Pu10<<" "<<Pu11<<" "<<Pu12<<" "<<Am1<<endl;
		exit(0);
	}
	// All value are molar (!weight)

	InputTree->Fill();

	return InputTree;
}
//________________________________________________________________________
double EQM_MLP_MOX::ExecuteTMVA(TTree* theTree)
{
	// --- Create the Reader object

	TMVA::Reader *reader = new TMVA::Reader( "Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	Float_t Pu8,Pu9,Pu10,Pu11,Pu12,Am1,BU,U5_enrichment;
	reader->AddVariable( "Pu8"  		,&Pu8 );
	reader->AddVariable( "Pu9"  		,&Pu9 );
	reader->AddVariable( "Pu10" 		,&Pu10);
	reader->AddVariable( "Pu11" 		,&Pu11);
	reader->AddVariable( "Pu12" 		,&Pu12);
	reader->AddVariable( "Am1"  		,&Am1 );
	reader->AddVariable( "U5_enrichment",&U5_enrichment );
	reader->AddVariable( "BU"   		,&BU );


	// --- Book the MVA methods

	// Book method MLP
	TString methodName = "MLP method";
	reader->BookMVA( methodName, fTMVAWeightPath );

	//theTree->SetBranchAddress("teneur",&teneur);
	theTree->SetBranchAddress( "Pu8"  			,&Pu8   );
	theTree->SetBranchAddress( "Pu9"  			,&Pu9   );
	theTree->SetBranchAddress( "Pu10" 			,&Pu10  );
	theTree->SetBranchAddress( "Pu11" 			,&Pu11  );
	theTree->SetBranchAddress( "Pu12" 			,&Pu12  );
	theTree->SetBranchAddress( "Am1"  			,&Am1   );
	theTree->SetBranchAddress( "U5_enrichment"  ,&U5_enrichment   );
	theTree->SetBranchAddress( "BU"   			,&BU 	);

	theTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;
	delete theTree;
	//cout<<val<<endl;
	return (double)val; //retourne teneur
}
//________________________________________________________________________
double EQM_MLP_MOX::GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp)
{
	return	ExecuteTMVA(CreateTMVAInputTree(Fissil,Fertil,BurnUp));
}
