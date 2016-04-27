#include "EquivalenceModel.hxx"
#include "EQM_ADS_MLP_RatioPuAM.hxx"
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
//		EQM_ADS_MLP_RatioPuAM
//
// The aim of these class is to constuct a fuel from a Multi layer perceptron (MLP).
//
//________________________________________________________________________

EQM_ADS_MLP_RatioPuAM::EQM_ADS_MLP_RatioPuAM(string TMVAWeightPath, double KeffAtBOC):EquivalenceModel(new CLASSLogger("EQM_ADS_MLP_RatioPuAM.log"))
{
	fKeffAtBOC 			= KeffAtBOC ;

	fTMVAWeightPath.push_back(TMVAWeightPath);

	ZAI N_93237(93,237,0);
	ZAI N_95241(95,241,0);
	ZAI N_95242(95,242,1);
	ZAI N_95243(95,243,0);
	ZAI N_96244(96,244,0);
	ZAI N_96245(96,245,0);
	ZAI N_96246(96,246,0);

	ZAI N_94238(94,238,0);
	ZAI N_94239(94,239,0);
	ZAI N_94240(94,240,0);
	ZAI N_94241(94,241,0);
	ZAI N_94242(94,242,0);

	fStreamList["Pu"] 			= N_94238*1+N_94239*1+N_94240*1+N_94241*1+N_94242*1;
	fStreamList["MA"] 			= N_93237*1+N_95241*1+N_95242*1+N_95243*1+N_96244*1+N_96245*1+N_96246*1;

	fFirstGuessContent["Pu"] 	=  0.50;
	fFirstGuessContent["MA"] 	=  1 - fFirstGuessContent["Pu"];
	
	fSpecificPower 				= 384e6 / 5.4e6; //Pth / M_HN = 71.1 W/g

	INFO<<"__An equivalence model of ADS has been define__"<<endl;
	INFO<<"\tThis model is based on a multi layer perceptron"<<endl;
	INFO<<"\t\tThe TMVA weight file is :"<<endl;
	INFO<<"\t\t\t"<<fTMVAWeightPath[0]<<endl;
}

//________________________________________________________________________
EQM_ADS_MLP_RatioPuAM::EQM_ADS_MLP_RatioPuAM(CLASSLogger* log, string TMVAWeightPath, double KeffAtBOC):EquivalenceModel(log)
{

	fKeffAtBOC 			= KeffAtBOC ;

	fTMVAWeightPath.push_back(TMVAWeightPath);

	ZAI N_93237(93,237,0);
	ZAI N_95241(95,241,0);
	ZAI N_95242(95,242,1);
	ZAI N_95243(95,243,0);
	ZAI N_96244(96,244,0);
	ZAI N_96245(96,245,0);
	ZAI N_96246(96,246,0);

	ZAI N_94238(94,238,0);
	ZAI N_94239(94,239,0);
	ZAI N_94240(94,240,0);
	ZAI N_94241(94,241,0);
	ZAI N_94242(94,242,0);

	fStreamList["Pu"] 			= N_94238*1+N_94239*1+N_94240*1+N_94241*1+N_94242*1;
	fStreamList["MA"] 			= N_93237*1+N_95241*1+N_95242*1+N_95243*1+N_96244*1+N_96245*1+N_96246*1;

	fFirstGuessContent["Pu"] 	=  0.50;
	fFirstGuessContent["MA"] 	=  1 - fFirstGuessContent["Pu"];
	
	fSpecificPower 				= 384e6 / 5.4e6; //Pth / M_HN = 71.1 W/g
/*
	SetRelativMassPrecision(0.0008);
	SetBurnUpPrecision(0.008);
	SetPCMPrecision(10);
*/
	INFO<<"__An equivalence model of ADS has been define__"<<endl;
	INFO<<"\tThis model is based on a multi layer perceptron"<<endl;
	INFO<<"\t\tThe TMVA weight file is :"<<endl;
	INFO<<"\t\t\t"<<fTMVAWeightPath[0]<<endl;

}

//________________________________________________________________________
map <string , vector<double> > EQM_ADS_MLP_RatioPuAM::BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray)
{
DBGL

	map <string , vector<double> > lambda ; // map containing name of the list and associated vector of proportions taken from stocks
cout<<1<<endl;	
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
cout<<2<<endl;	

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
cout<<3<<endl;	
	
	/**** Some initializations **/

	StocksTotalMassCalculation(StreamArray);

	fActualMassContentInFuel	= GetBuildFuelFirstGuess();
	fActualMolarContentInFuel	= GetBuildFuelFirstGuess();	
	int loopCount 			= 0;
	bool FuelBuiltCorrectly	 	= false ; 

	map <string , IsotopicVector > IVStream;
	map < string , double > MaterialMassNeeded ;
	map < string , double > LambdaNeeded;
	map < string , double > DeltaMass;
cout<<4<<endl;	
	
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		LambdaNeeded[(*it_s_vIV).first] = 0;

	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		MaterialMassNeeded[(*it_s_vIV).first] = HMMass*fActualMassContentInFuel[(*it_s_vIV).first]; 
		DeltaMass[(*it_s_vIV).first] =  - MaterialMassNeeded[(*it_s_vIV).first];
	}
cout<<5<<endl;	

	do
	{	
		map < string , double >  LambdaPreviousStep;
cout<<"5-1 "<<LambdaNeeded["Pu"]<<endl;	

		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
			LambdaPreviousStep[(*it_s_D).first] =LambdaNeeded[(*it_s_D).first]; 
cout<<"5-2 "<<LambdaNeeded["Pu"]<<endl;	
		
		BreakReturnLambda = false;
		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
			LambdaNeeded[(*it_s_vIV).first] = LambdaCalculation((*it_s_vIV).first, LambdaPreviousStep[(*it_s_vIV).first], MaterialMassNeeded[(*it_s_vIV).first], DeltaMass[(*it_s_vIV).first], StreamArray[(*it_s_vIV).first]);

for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		cout<<(*it_s_vIV).first<<endl;
		cout<<LambdaPreviousStep[(*it_s_vIV).first]<<endl;
		cout<<MaterialMassNeeded[(*it_s_vIV).first]<<endl;
		cout<<DeltaMass[(*it_s_vIV).first]<<endl;
		//cout<<StreamArray[(*it_s_vIV).first]<<endl;
	}

cout<<"5-3 "<<LambdaNeeded["Pu"]<<" "<<LambdaNeeded["MA"]<<endl;	

		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
		{		
			if( LambdaNeeded[(*it_s_D).first] == -1 )
			{
				SetLambdaToErrorCode(lambda[(*it_s_D).first]);
				WARNING << "Not enough : "<< (*it_s_D).first <<" to build fuel." << endl;	
				BreakReturnLambda = true; 				
			}
		}
cout<<"5-4 "<<LambdaNeeded["Pu"]<<endl;	
		
		if(BreakReturnLambda) { return lambda;}
		
		BreakReturnLambda = false; 
cout<<"5-5 "<<LambdaNeeded["Pu"]<<endl;	
		
		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
			SetLambda(lambda[(*it_s_D).first], LambdaNeeded[(*it_s_D).first]); 
cout<<"5-6 "<<LambdaNeeded["Pu"]<<endl;	
	
		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
			IVStream[(*it_s_vIV).first].Clear();
cout<<"5-7 "<<LambdaNeeded["Pu"]<<endl;	

		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		{	
			for(int i=0; i<StreamArray[(*it_s_vIV).first].size(); i++)
			{
				IVStream[(*it_s_vIV).first]  +=  lambda[(*it_s_vIV).first][i] * StreamArray[(*it_s_vIV).first][i];	
			}
		}
cout<<"5-8 "<<LambdaNeeded["Pu"]<<endl;	
		
		if (loopCount > fMaxInterration)
		{
			ERROR << "Too much iterration in BuildFuel Method !";
			ERROR << "Need improvement in fuel fabrication ! Ask for it or D.I.Y. !!" << endl;
			exit(1);
		}
cout<<"5-9 "<<LambdaNeeded["Pu"]<<endl;	
		
		/* Calcul the quantity of this composition needed to reach the burnup */

		map < string , double>  MaterialMolarContent = GetMolarFraction(IVStream, BurnUp);
cout<<"5-10 "<<LambdaNeeded["Pu"]<<endl;	
		
		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
		{	
			if( MaterialMolarContent[(*it_s_D).first] < 0 || MaterialMolarContent[(*it_s_D).first] > 1 )
			{
				SetLambdaToErrorCode(lambda[(*it_s_D).first]);
				WARNING << "GetFissileMolarFraction return negative or greater than one value, at least for one material : "<<(*it_s_D).first;
				BreakReturnLambda = true; 				
			}
		}
		if(BreakReturnLambda) { return lambda;}
cout<<"5-11 "<<LambdaNeeded["Pu"]<<endl;	
		
		map < string , double >  MaterialMassContent;
		map < string , double >  MaterialMeanMolar;
		map < string , double >  MaterialMassAvailable;
		map < string , bool > 	  CheckOnMass;

		double FuelMeanMolar   = 0;

		for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
			MaterialMassAvailable[(*it_s_IV).first] = IVStream[(*it_s_IV).first].GetTotalMass()*1e06; 
cout<<"5-12 "<<LambdaNeeded["Pu"]<<endl;	
		
		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
		{
			MaterialMeanMolar[(*it_s_D).first] = IVStream[(*it_s_D).first].GetMeanMolarMass();
			FuelMeanMolar +=  MaterialMolarContent[(*it_s_D).first] * MaterialMeanMolar[(*it_s_D).first];
		}
cout<<"5-13 "<<LambdaNeeded["Pu"]<<endl;	
		
		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
			MaterialMassContent [(*it_s_D).first] =  MaterialMolarContent[(*it_s_D).first] * MaterialMeanMolar[(*it_s_D).first] / FuelMeanMolar; 
cout<<"5-14 "<<LambdaNeeded["Pu"]<<endl;	

		fActualMolarContentInFuel = MaterialMolarContent; //fActualContent can be accessed by a derivated EquivalenModel to accelerate GetFissileMolarFraction function (exemple in EQM_MLP_Kinf)
cout<<"5-15 "<<LambdaNeeded["Pu"]<<endl;	
		
		for( it_s_D = MaterialMassContent.begin();  it_s_D != MaterialMassContent.end(); it_s_D++)
		{	
			MaterialMassNeeded[(*it_s_D).first] = HMMass*MaterialMassContent[(*it_s_D).first]; 
			DeltaMass[(*it_s_D).first] = MaterialMassAvailable[(*it_s_D).first] - MaterialMassNeeded[(*it_s_D).first];
		}
cout<<"5-16 "<<LambdaNeeded["Pu"]<<endl;	
		
		for( it_s_D = MaterialMassNeeded.begin();  it_s_D != MaterialMassNeeded.end(); it_s_D++)
		{
			double DeltaM = fabs(MaterialMassNeeded[(*it_s_D).first] - MaterialMassAvailable[(*it_s_D).first]) / HMMass ;
			
			if(DeltaM<fRelativMassPrecision) {CheckOnMass[(*it_s_D).first] = true;}
			else{CheckOnMass[(*it_s_D).first] = false;}
		}
cout<<"5-17 "<<LambdaNeeded["Pu"]<<endl;	
	
		FuelBuiltCorrectly = true; 
		for( it_s_B = CheckOnMass.begin();  it_s_B != CheckOnMass.end(); it_s_B++)
		{
			if(!CheckOnMass[(*it_s_B).first]) {FuelBuiltCorrectly = false;}
		}
cout<<"5-18 "<<LambdaNeeded["Pu"]<<endl;	

		loopCount++;
		
	}while(!FuelBuiltCorrectly);
cout<<6<<endl;	

	for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
		(*this).isIVInDomain(IVStream[(*it_s_IV).first]);
cout<<7<<endl;	
	
	for( it_s_D = MaterialMassNeeded.begin();  it_s_D != MaterialMassNeeded.end(); it_s_D++)
		DBGV( "Weight percent : "<<(*it_s_D).first<<" "<< MaterialMassNeeded[(*it_s_D).first]/HMMass);
cout<<8<<endl;	
	
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
TTree* EQM_ADS_MLP_RatioPuAM::CreateTMVAInputTree(IsotopicVector TheFuel)
{
	TTree*   InputTree = new TTree("EQTMP", "EQTMP");

	float MA_93_237	= 0;	InputTree->Branch(	"MA_93_237" ,&MA_93_237   ,"MA_93_237");
	float MA_95_241	= 0;	InputTree->Branch(	"MA_95_241" ,&MA_95_241   ,"MA_95_241");
	float MA_95_242	= 0;	InputTree->Branch(	"MA_95_242" ,&MA_95_242   ,"MA_95_242");
	float MA_95_243	= 0;	InputTree->Branch(	"MA_95_243" ,&MA_95_243   ,"MA_95_243");
	float MA_96_244	= 0;	InputTree->Branch(	"MA_96_244" ,&MA_96_244   ,"MA_96_244");
	float MA_96_245	= 0;	InputTree->Branch(	"MA_96_245" ,&MA_96_245   ,"MA_96_245");
	float MA_96_246	= 0;	InputTree->Branch(	"MA_96_246" ,&MA_96_246   ,"MA_96_246");
	float Pu_94_238	= 0;	InputTree->Branch(	"Pu_94_238" ,&Pu_94_238   ,"Pu_94_238");
	float Pu_94_239	= 0;	InputTree->Branch(	"Pu_94_239" ,&Pu_94_239   ,"Pu_94_239");
	float Pu_94_240	= 0;	InputTree->Branch(	"Pu_94_240" ,&Pu_94_240   ,"Pu_94_240");
	float Pu_94_241	= 0;	InputTree->Branch(	"Pu_94_241" ,&Pu_94_241   ,"Pu_94_241");
	float Pu_94_242	= 0;	InputTree->Branch(	"Pu_94_242" ,&Pu_94_242   ,"Pu_94_242");
	float KEFF 	= 0;		InputTree->Branch(	"KEFF" ,&KEFF ,"KEFF");

	MA_93_237	 = TheFuel.GetZAIIsotopicQuantity(93,237,0);
	MA_95_241	 = TheFuel.GetZAIIsotopicQuantity(95,241,0);
	MA_95_242	 = TheFuel.GetZAIIsotopicQuantity(95,242,1);
	MA_95_243	 = TheFuel.GetZAIIsotopicQuantity(95,243,0);
	MA_96_244	 = TheFuel.GetZAIIsotopicQuantity(96,244,0);
	MA_96_245	 = TheFuel.GetZAIIsotopicQuantity(96,245,0);
	MA_96_246	 = TheFuel.GetZAIIsotopicQuantity(96,246,0);
	Pu_94_238	 = TheFuel.GetZAIIsotopicQuantity(94,238,0);
	Pu_94_239	 = TheFuel.GetZAIIsotopicQuantity(94,239,0);
	Pu_94_240	 = TheFuel.GetZAIIsotopicQuantity(94,240,0);
	Pu_94_241	 = TheFuel.GetZAIIsotopicQuantity(94,241,0);
	Pu_94_242	 = TheFuel.GetZAIIsotopicQuantity(94,242,0);

	KEFF = fKeffAtBOC;

	InputTree->Fill();

	return InputTree;
}
//________________________________________________________________________
double EQM_ADS_MLP_RatioPuAM::ExecuteTMVA(TTree* TheTree, string WeightPath)
{
	// --- Create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "Silent" );
	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used

	Float_t MA_93_237=0; reader->AddVariable( "MA_93_237", &MA_93_237);
	Float_t MA_95_241=0; reader->AddVariable( "MA_95_241", &MA_95_241);
	Float_t MA_95_242=0; reader->AddVariable( "MA_95_242", &MA_95_242);
	Float_t MA_95_243=0; reader->AddVariable( "MA_95_243", &MA_95_243);
	Float_t MA_96_244=0; reader->AddVariable( "MA_96_244", &MA_96_244);
	Float_t MA_96_245=0; reader->AddVariable( "MA_96_245", &MA_96_245);
	Float_t MA_96_246=0; reader->AddVariable( "MA_96_246", &MA_96_246);
	Float_t Pu_94_238=0; reader->AddVariable( "Pu_94_238", &Pu_94_238);
	Float_t Pu_94_239=0; reader->AddVariable( "Pu_94_239", &Pu_94_239);
	Float_t Pu_94_240=0; reader->AddVariable( "Pu_94_240", &Pu_94_240);
	Float_t Pu_94_241=0; reader->AddVariable( "Pu_94_241", &Pu_94_241);
	Float_t Pu_94_242=0; reader->AddVariable( "Pu_94_242", &Pu_94_242); 
	Float_t KEFF = 0.00; reader->AddVariable( "KEFF", &KEFF);

	// --- Book the MVA methods

	// Book method MLP
	TString methodName = "MLP method";
	reader->BookMVA( methodName, fTMVAWeightPath[0] );

	TheTree->SetBranchAddress( "MA_93_237", &MA_93_237);
	TheTree->SetBranchAddress( "MA_95_241", &MA_95_241);
	TheTree->SetBranchAddress( "MA_95_242", &MA_95_242);
	TheTree->SetBranchAddress( "MA_95_243", &MA_95_243);
	TheTree->SetBranchAddress( "MA_96_244", &MA_96_244);
	TheTree->SetBranchAddress( "MA_96_245", &MA_96_245);
	TheTree->SetBranchAddress( "MA_96_246", &MA_96_246);
	TheTree->SetBranchAddress( "Pu_94_238", &Pu_94_238);
	TheTree->SetBranchAddress( "Pu_94_239", &Pu_94_239);
	TheTree->SetBranchAddress( "Pu_94_240", &Pu_94_240);
	TheTree->SetBranchAddress( "Pu_94_241", &Pu_94_241);
	TheTree->SetBranchAddress( "Pu_94_242", &Pu_94_242); 
	TheTree->SetBranchAddress( "KEFF", &KEFF);

	TheTree->GetEntry(0);

	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;


	return (double)val; //return Frac MA
}

//________________________________________________________________________
map < string , double>  EQM_ADS_MLP_RatioPuAM::GetMolarFraction(map <string , IsotopicVector > IVStream, double TargetBU)
{

cout<<"HERE..."<<endl;

	map < string , double> MolarFraction ;

	if (fTotalMassInStocks["Pu"] >= 5.4e6/2.)
	{

cout<<"1 - "<<fTotalMassInStocks["Pu"]<<" "<<fTotalMassInStocks["MA"]<<endl;

		MolarFraction["MA"] = 0.5;
		MolarFraction["Pu"] = 0.5;
	} 
	else
	{
		MolarFraction["Pu"] = fTotalMassInStocks["Pu"] / (5.4e6);
		MolarFraction["MA"] = 1-MolarFraction["Pu"];

cout<<"2 - "<<fTotalMassInStocks["Pu"]<<" "<<fTotalMassInStocks["MA"]<<" "<<MolarFraction["Pu"]<<endl;

	}

cout<<"NOT HERE..."<<endl;


// TO ADD
// 
// Model For Fraction / Keff = 0.96
// 
// Model where fractions are chosen
// 




return MolarFraction;
}
