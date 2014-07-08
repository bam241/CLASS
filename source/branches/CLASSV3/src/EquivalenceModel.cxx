#include "EquivalenceModel.hxx"






EquivalenceModel::EquivalenceModel():CLASSObject()
{

}

EquivalenceModel::EquivalenceModel(CLASSLogger* log):CLASSObject(log)
{

}

EquivalenceModel::~EquivalenceModel()
{

}




//________________________________________________________________________
vector<double> EquivalenceModel::BuildFuel(double BurnUp, double HMMass,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray)
{	//cout<<"entrance in EquivalenceModel::BuildFuel"<<endl;

	HMMass*=1e6;//tons to gram

/****Some initializations**/
	vector<double> lambda ; //vector of portion of stocks taken (fissile & fertil)
	for(int i=0;i < (int)FissilArray.size() + (int)FertilArray.size();i++ )
		lambda.push_back(0);

	fOld_Lambda_Tot[0] = 0;
	
	FindLambdaMax(0,FissilArray.size()-1,FissilArray,HMMass);
	FindLambdaMax(FissilArray.size(),lambda.size()-1,FertilArray,HMMass);

	IsotopicVector Fertile;
	IsotopicVector Fissile;
	double AvailablePuMass=0;
	double PuMassNeeded=1000; //At present time I have no clue what is the requiered Pu mass, I assume at least 1kg is needed
	double WeightPuContent=0;
	double MassPrecision=100; //Mass precision is 100 grams

	while(  fabs(PuMassNeeded - AvailablePuMass) > MassPrecision )
	{	//cout<<"********FISSILE**********"<<endl;
		//Increase the portion of the stock stokID taken, according to the followings variables
		double DeltaM=PuMassNeeded-AvailablePuMass;
		//cout<<"DeltaM fis "<<DeltaM<<" PuMassNeeded "<<PuMassNeeded<<" Old Available Pu "<<AvailablePuMass;

		GuessLambda( lambda,0,FissilArray.size()-1, DeltaM, FissilArray ,HMMass);
		//Build the Plutonium vector from stocks
		Fissile.Clear();
		for(int i=0;i<(int)FissilArray.size();i++)
		{
			IsotopicVector IVTMP =	lambda[i] * FissilArray[i];
			Fissile += IVTMP;
		}	

		AvailablePuMass = Fissile.GetTotalMass() * 1e6; //in grams
		//cout<< "HMMass "<< HMMass <<" New Pu Mass ava : "<<AvailablePuMass<<endl;
		//Build uranium vector from stocks

	//cout<<"********CALCUL FERTILE**********"<<endl;
		/****Some initializations**/
		double FertilMassNeeded = HMMass - AvailablePuMass;
		double AvailableFertilMass = 0;
		fOld_Lambda_Tot[1] = 0;
		for(int i=(int)FissilArray.size() ;i < (int)FissilArray.size()+(int)FertilArray.size() ;i++ )
			lambda[i]=0;

		while( fabs(FertilMassNeeded - AvailableFertilMass) > MassPrecision  )
		{	
			double DeltaM=FertilMassNeeded-AvailableFertilMass;
			GuessLambda( lambda,  FissilArray.size(), lambda.size()-1, DeltaM, FertilArray,HMMass);

			int j=-1;
			Fertile.Clear();
			for(int i=(int)FissilArray.size();i<(int)FissilArray.size()+(int)FertilArray.size() ;i++)
			{	
				j++;
				IsotopicVector IVTMP = lambda[i] * FertilArray[j];
				Fertile += IVTMP ;
				//cout<<"lambda[i] "<< lambda[i]<<endl;
			}
			AvailableFertilMass=Fertile.GetTotalMass() * 1e6; //in grams

			if( j+1 == int( FertilArray.size() ) && FertilMassNeeded  > AvailableFertilMass ) //if this is the last stock and it's not enought
			{
				WARNING <<"You requiere more depleted uranium "<<"("<<DeltaM/1e6<<" t needed) ! Reactor not fill"<<endl;
				for(int i=0 ; i<int(lambda.size()) ; i++)
					lambda[i]=0;
				break;
			}

		}
	//cout<<"********FERTILE OK**********"<<endl;
		/*Calcul the quantity of this composition needed to reach the burnup*/
		double MolarPuContent = GetFissileMolarFraction(Fissile, Fertile, BurnUp);

		double MeanMolarPu = Fissile.MeanMolar();
		double MeanMolarDepletedU = Fertile.MeanMolar();
		double MeanMolar   = MeanMolarPu*MolarPuContent + (1-MolarPuContent)*MeanMolarDepletedU;

		WeightPuContent = MolarPuContent * MeanMolarPu/MeanMolar ;
		//cout<<"WeightPuContent  "<<WeightPuContent<<endl;
		PuMassNeeded = WeightPuContent  *  HMMass ;
		//cout<<"PuMassNeeded "<<PuMassNeeded<<endl;
//cout<<"********FIN  BOUCLE**********"<<endl;
	}

	//cout<<"exiting  EquivalenceModel::BuildFuel"<<endl;
	return lambda;
}
//________________________________________________________________________
void EquivalenceModel::SetLambda(vector<double>& lambda ,int FirstStockID, int LastStockID, double LAMBDA_TOT)
{
		if(LAMBDA_TOT > LastStockID - FirstStockID + 1 )
		{
			ERROR << " FATAL ERROR " <<endl;
			exit(0);
		}

		for(int i = FirstStockID ; i<=LastStockID ;i++) //set to 0 all non touched value (to be sure)
			lambda[i] = 0  ;

		int PartieEntier = floor(LAMBDA_TOT);
		double PartieDecimal = LAMBDA_TOT - PartieEntier;

		for(int i=FirstStockID  ; i< FirstStockID +PartieEntier ;i++ )
				lambda[i]=1;

		lambda[FirstStockID + PartieEntier] = PartieDecimal  ;

}
//________________________________________________________________________
void EquivalenceModel::GuessLambda(vector<double>& lambda,int FirstStockID, int LastStockID, double DeltaM, vector<IsotopicVector> Stocks, double  HMMass)
{


	int f=1;
	if(FirstStockID == 0)	//code :: 0 is fissile 1 fertile
		f=0;

	double	LAMBDA_TOT=0;
	for (int i = FirstStockID; i <=LastStockID ; i++)
		LAMBDA_TOT+=lambda[i]	;


	if( LAMBDA_TOT == 0 ) //Initialization je cherche les lambda telle que j'ai 5% de HM  !!!
	{	//cout<<"INITIALIZATION"<<endl;
		double	MASS = 0 ; 
		int ID_max = 0;
		while( MASS < 0.05*HMMass )
		{		
			double StockMass = Stocks[ID_max].GetTotalMass() * 1e6;

			if( StockMass > HMMass )
				LAMBDA_TOT += 0.01*HMMass/StockMass;
			else
				 LAMBDA_TOT += 0.01;

			if(LAMBDA_TOT > ID_max+1)	
			{	
				LAMBDA_TOT=ID_max+1;
				ID_max++ ;
				if( ID_max >=(int) Stocks.size())
				{
					cout<<"FATAL ERROR Not enought minerals"<<endl;
					exit(0);
				}
			}	

			SetLambda(lambda,FirstStockID,LastStockID,LAMBDA_TOT );	

			IsotopicVector test;
			int j=0;
			for(int i=FirstStockID;i<=LastStockID;i++)
			{
				IsotopicVector IVTMP =	lambda[i] * Stocks[j];
				test += IVTMP;
				j++;
			}	
			MASS = test.GetTotalMass() * 1e6; //in grams
			//cout<<"LAMBDA_TOT "<<LAMBDA_TOT<<" MASS "<<MASS<<endl;
		}
	}

	else  if( DeltaM > 0) 
	{
		fOld_Lambda_Tot[f] = LAMBDA_TOT;
		LAMBDA_TOT += (fLambda_max[f] - LAMBDA_TOT)/2.;
		//cout<<"DM>0 LAMBDA_TOT "<<LAMBDA_TOT<<endl;
		SetLambda(lambda,FirstStockID,LastStockID,LAMBDA_TOT );


	}

	else  if( DeltaM < 0) 
	{
		//cout <<"fOld_Lambda_Tot[f] "<<fOld_Lambda_Tot[f]<<" LAMBDA_TOT "<<LAMBDA_TOT<< endl;
		LAMBDA_TOT = LAMBDA_TOT - fabs(fOld_Lambda_Tot[f] - LAMBDA_TOT) / 2. ;
		//cout<<"DM<0 LAMBDA_TOT "<<LAMBDA_TOT<<endl;
		SetLambda(lambda,FirstStockID,LastStockID,LAMBDA_TOT );
	}


}
//________________________________________________________________________
void EquivalenceModel::FindLambdaMax(int FirstStockID, int LastStockID, vector<IsotopicVector> Stocks, double  HMMass)
{
	int f=1;
	if(FirstStockID == 0)	//code :: 0 is fissile 1 fertile
		f=0;

	//je cherche les lambda telle que j'ai 100% de HM  ou bien la mass total des stocks

	vector<double> lambda;
	double TotStockMass = 0;
	for (int i = 0 ; i <= LastStockID-FirstStockID  ; i++)
	{	TotStockMass += Stocks[i].GetTotalMass() * 1e6;
		lambda.push_back(0);
	}

	if(TotStockMass <= HMMass)
		fLambda_max[f]= LastStockID - FirstStockID + 1;

	else
	{
		double LAMBDA=0;
		double MASS = 0;
		int ID_max = 0;
		double MassPrecision = 0;
		while( fabs(MASS - HMMass) > MassPrecision )
		{		
			double StockMass = Stocks[ID_max].GetTotalMass() * 1e6;

			if( StockMass > HMMass )
			{	
				MassPrecision = 0.01*HMMass/StockMass * HMMass;
				LAMBDA += 0.005*HMMass/StockMass;

			}	
			else
				 LAMBDA += 1;

			if(LAMBDA > ID_max+1)	
			{	
				LAMBDA=ID_max+1;
				ID_max++ ;
				if( ID_max >=(int) Stocks.size())
				{
					ERROR << " FATAL ERROR BLG codes like shit " << endl;
					exit(0);
				}
			}	

			SetLambda(lambda,FirstStockID,LastStockID,LAMBDA );	

			IsotopicVector test;
			int j=0;
			for(int i=FirstStockID;i<=LastStockID;i++)
			{
				IsotopicVector IVTMP =	lambda[i] * Stocks[j];
				test += IVTMP;
				j++;
			}	

			MASS = test.GetTotalMass() * 1e6; //in grams
		}

		fLambda_max[f]=LAMBDA;
	//cout<<"MassPrecision "<<MassPrecision<<endl;

	}

}