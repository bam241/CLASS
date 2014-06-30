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
}
//________________________________________________________________________
void EQM_MLP_MOX::CreateTMVAInputTree(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/

	TFile*   InputFile = new TFile("./EQMTMP.root","RECREATE");
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

	InputFile->Write();
	delete InputTree;
	InputFile-> Close();
	delete InputFile;
}
//________________________________________________________________________
double EQM_MLP_MOX::ExecuteTMVA()
{
	// --- Create the Reader object

	TMVA::Reader *reader = new TMVA::Reader( "Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	Float_t Pu8,Pu9,Pu10,Pu11,Pu12,Am1,BU,teneur;
	reader->AddVariable( "Pu8"  ,&Pu8 );
	reader->AddVariable( "Pu9"  ,&Pu9 );
	reader->AddVariable( "Pu10" ,&Pu10);
	reader->AddVariable( "Pu11" ,&Pu11);
	reader->AddVariable( "Pu12" ,&Pu12);
	reader->AddVariable( "Am1"  ,&Am1 );
	reader->AddVariable( "BU"   ,&BU );


	// --- Book the MVA methods

	// Book method MLP
	TString methodName = "MLP method";
	reader->BookMVA( methodName, fTMVAWeightPath );

	// Prepare input tree
	TFile *input(0);
	if (!gSystem->AccessPathName( "./EQMTMP.root" )) {
		input = TFile::Open( "./EQMTMP.root" ); // check if file in local directory exists
	}
	if (!input) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}

	TTree* theTree = (TTree*)input->Get("EQTMP");

	theTree->SetBranchAddress("teneur",&teneur);
	theTree->SetBranchAddress( "Pu8"  ,&Pu8   );
	theTree->SetBranchAddress( "Pu9"  ,&Pu9   );
	theTree->SetBranchAddress( "Pu10" ,&Pu10  );
	theTree->SetBranchAddress( "Pu11" ,&Pu11  );
	theTree->SetBranchAddress( "Pu12" ,&Pu12  );
	theTree->SetBranchAddress( "Am1"  ,&Am1   );
	theTree->SetBranchAddress( "BU"   ,&BU 	  );

	theTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;
	delete theTree;
	input->Close();
	delete input;
	//cout<<"....done"<<endl;

	//cout<<val<<endl;
	return (double)val; //retourne teneur
}
//________________________________________________________________________
void EQM_MLP_MOX::GuessLambda(double& lambda, int& StockID,int FirstStockID, int LastStockID, double DeltaM,double StockHM)
{

	double Threshold = 50 ; //50 grams 


	/****Initialization***/
	if(lambda==0) 
		lambda=0.5;
	
	/********dichotomie**************/
	if(DeltaM>0)  //MassNeeded - AvailableMass
	{
		lambda+=lambda/2.;

		if(lambda >= 1) //this stock is not enought go to next one
		{
			lambda=1;

		if( !( StockID+1 == LastStockID) ) //if its the last stock don't try to look the next one (segfault otherwize)
				StockID++;
		}

	}

	else if(DeltaM<0) 
	{
		lambda-=lambda/2.;

		if(lambda*StockHM < Threshold) //if only 50 grams left in actual stock go to previous one
		{
			lambda=0;
			StockID--;
			if(StockID<FirstStockID)
			{
				cout<<"Critical error  EQM_MLP_MOX::GuessLambda"<<endl;
				cout<<"Contact BLG"<<endl;
				exit(1);
			}	

		}

	}

}
//________________________________________________________________________
vector<double> EQM_MLP_MOX::BuildFuel(double BurnUp, double HMMass,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray)
{

	HMMass*=1e6;//tons to gram
			
	vector<double> lambda ; //vector of portion of stocks taken (fissile & fertil)
	for(int i=0;FissilArray.size() + FertilArray.size();i++ );
		lambda.push_back(0);

		/*******************Depleted Uranium Vector**************************/
	/*	ZAI U5(92,235,0);
		ZAI U8(92,238,0);
		double U5_enrich= 0.0025;
		double MeanMolarDepletedU = U5.GetMass()*U5_enrich + (1-U5_enrich)*U8.GetMass();
		double AVOGADRO = 6.02214129e23;
		//Building a isotopic vector with a total mass of HMMass (we are going to take just a share of it but want to be sure we have enought)
		DepletedUranium.Add(U5,   U5_enrich   *HMMass/MeanMolarDepletedU*AVOGADRO);
		DepletedUranium.Add(U8,  (1-U5_enrich)*HMMass/MeanMolarDepletedU*AVOGADRO);
		*/
		IsotopicVector DepletedUranium;
		IsotopicVector PlutoniumVector;
		int StockID = 0 ;
		double AvailablePuMass=0;
		double PuMassNeeded=1000; //At present time I have no clue what is the requiered Pu mass, I assume at least 1kg is needed 
		double WeightPuContent=0;
		double MassPrecision=100; //Mass precision is 100 grams
		while(  abs(PuMassNeeded - AvailablePuMass) > MassPrecision )
		{
			//Increase the portion of the stock stokID taken, according to the followings variables
			double DeltaM=PuMassNeeded-AvailablePuMass;
			GuessLambda( lambda[StockID], StockID,0,FissilArray.size()-1, DeltaM, FissilArray[StockID].GetTotalMass() * 1e6 );

			//Build the Plutonium vector from stocks
			PlutoniumVector.Clear();
			for(int i=0;i<=StockID;i++) 
				PlutoniumVector += lambda[i] * FissilArray[i];

			AvailablePuMass = PlutoniumVector.GetTotalMass() * 1e6; //in grams

			//Build uranium vector from stocks
			int FertileStockID = FissilArray.size();
			double FertilMassNeeded = HMMass - AvailablePuMass;
			double AvailableFertilMass = 0;

			while( abs(FertilMassNeeded - AvailableFertilMass) > MassPrecision  )
			{	
				double DeltaM=FertilMassNeeded-AvailableFertilMass;
				GuessLambda( lambda[FertileStockID], FertileStockID,FissilArray.size(),FertilArray.size()-1, DeltaM, FertilArray[FertileStockID-FissilArray.size()].GetTotalMass() * 1e6 ) ;

				int j=-1;
				DepletedUranium.Clear();
				for(int i=int(FissilArray.size());i<=FertileStockID;i++)
				{	DepletedUranium += lambda[i] * FertilArray[j];
					j++;
				}	
				AvailableFertilMass=DepletedUranium.GetTotalMass() * 1e6; //in grams

				if( j+1 == int( FertilArray.size() ) && FertilMassNeeded  > AvailableFertilMass ) //if this is the last stock and it's not enought
				{
					cout<<"You requiere more depleted uranium "<<"("<<DeltaM/1e6<<" t needed) ! Reactor not fill"<<endl;
					for(int i=0 ; i<int(lambda.size()) ; i++)
						lambda[i]=-1;
					break;
				}

			}	
			/*Calcul the quantity of this composition needed to reach the burnup*/
			CreateTMVAInputTree(PlutoniumVector,DepletedUranium,BurnUp);
			double MolarPuContent = ExecuteTMVA();
			system("rm EQMTMP.root ");

			double MeanMolarPu = PlutoniumVector.MeanMolar();
			double MeanMolarDepletedU = DepletedUranium.MeanMolar();
			double MeanMolar   = MeanMolarPu*MolarPuContent + (1-MolarPuContent)*MeanMolarDepletedU;

			WeightPuContent = MolarPuContent * MeanMolarPu/MeanMolar ;

			PuMassNeeded = WeightPuContent  *  HMMass ;

			if( StockID+1 == int( FissilArray.size() ) && PuMassNeeded  > AvailablePuMass ) //if this is the last stock and it's not enought
			{
				cout<<"You requiere more (or better) plutonium !! Reactor not fill"<<endl;
				for(int i=0 ; i<int(lambda.size()) ; i++)
					lambda[i]=-1;
				break;
			}

		}

return lambda;
}
//________________________________________________________________________
