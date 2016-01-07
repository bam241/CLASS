#include "EquivalenceModel.hxx"
#include "EQM_PWR_MLP_MOX.hxx"
#include "CLASSLogger.hxx"
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
//		EQM_PWR_MLP_MOX
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For REP MOX use
//
//________________________________________________________________________

EQM_PWR_MLP_MOX::EQM_PWR_MLP_MOX(string TMVAWeightPath):EquivalenceModel(new CLASSLogger("EQM_PWR_MLP_MOX.log"))
{
	fTMVAWeightPath = TMVAWeightPath;

    reader = new TMVA::Reader( "Silent" );
    reader->BookMVA( "MLP Method", fTMVAWeightPath );
    reader->AddVariable( "BU"   		,&BU );
    reader->AddVariable( "U5_enrichment",&U5_enrichment );
    reader->AddVariable( "Pu8"  		,&Pu8 );
    reader->AddVariable( "Pu9"  		,&Pu9 );
    reader->AddVariable( "Pu10" 		,&Pu10);
    reader->AddVariable( "Pu11" 		,&Pu11);
    reader->AddVariable( "Pu12" 		,&Pu12);
    reader->AddVariable( "Am1"  		,&Am1 );


	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	double U5_enrich = 0.0025;

	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);

	fFissileList = Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1;
	fFertileList = U5*U5_enrich + U8*(1-U5_enrich);

	SetBuildFuelFirstGuess(0.04);

	INFO << "__An equivalence model of PWR MOX has been define__" << endl;
	INFO << "\tThis model is based on a multi layer perceptron" << endl;
	INFO << "\t\tThe TMVA weight file is :" << endl;
	INFO << "\t\t\t" << fTMVAWeightPath << endl;
	EquivalenceModel::PrintInfo();

}

//________________________________________________________________________
EQM_PWR_MLP_MOX::EQM_PWR_MLP_MOX(CLASSLogger* log, string TMVAWeightPath):EquivalenceModel(log)
{
	fTMVAWeightPath = TMVAWeightPath;

    reader = new TMVA::Reader( "Silent" );
    reader->BookMVA( "MLP Method", fTMVAWeightPath );
    reader->AddVariable( "BU"   		,&BU );
    reader->AddVariable( "U5_enrichment",&U5_enrichment );
    reader->AddVariable( "Pu8"  		,&Pu8 );
    reader->AddVariable( "Pu9"  		,&Pu9 );
    reader->AddVariable( "Pu10" 		,&Pu10);
    reader->AddVariable( "Pu11" 		,&Pu11);
    reader->AddVariable( "Pu12" 		,&Pu12);
    reader->AddVariable( "Am1"  		,&Am1 );



	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	double U5_enrich = 0.0025;

	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);

	fFissileList = Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1;
	fFertileList = U5*U5_enrich + U8*(1-U5_enrich);

	SetBuildFuelFirstGuess(0.04);

	INFO << "__An equivalence model of PWR MOX has been define__" << endl;
	INFO << "\tThis model is based on a multi layer perceptron" << endl;
	INFO << "\t\tThe TMVA weight file is :" << endl;
	INFO << "\t\t\t" << fTMVAWeightPath << endl;
	EquivalenceModel::PrintInfo();

}

//________________________________________________________________________
void EQM_PWR_MLP_MOX::UpdateInputComposition(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp)
{

    
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

    double TOTPU = (Pu8+Pu9+Pu10+Pu11+Pu12+Am1);

    Pu8 = Pu8  / TOTPU;
    Pu9 = Pu9  / TOTPU;
    Pu10 = Pu10 / TOTPU;
    Pu11 = Pu11 / TOTPU;
    Pu12 = Pu12 / TOTPU;
    Am1 = Am1  / TOTPU;

    U5_enrichment = U5 / UTOT;

	BU = BurnUp;
	if(Pu8 + Pu9 + Pu10 + Pu11 + Pu12 + Am1 > 1.00001 )//?????1.00001??? I don't know it! goes in condition if  = 1 !! may be float/double issue ...
	{
		ERROR << Pu8 << " " << Pu9 << " " << Pu10 << " " << Pu11 << " " << Pu12 << " " << Am1 << endl;
		exit(0);
	}
}
//________________________________________________________________________
double EQM_PWR_MLP_MOX::ExecuteTMVA(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp)
{
    UpdateInputComposition(Fissil, Fertil, BurnUp);
    Float_t val = (reader->EvaluateRegression( "MLP method" ))[0];

    return (double)val; //retourne teneur

}
//________________________________________________________________________
double EQM_PWR_MLP_MOX::GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp)
{DBGL
    return	ExecuteTMVA( Fissil, Fertil, BurnUp);
}
