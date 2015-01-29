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
{	
DBGL
	HMMass*=1e6;//tons to gram

	/****Some initializations**/
	vector<double> lambda ; //vector of portion of stocks taken (fissile & fertil)
	for(int i = 0 ; i < (int)FissilArray.size() + (int)FertilArray.size() ; i++ )
		lambda.push_back(0);


	if( (int)FissilArray.size()==0 )
	{	WARNING<<" No fissile stocks available ! Fuel not build"<<endl;
		lambda[0] = -1;
		return lambda;
	}


	fOld_Lambda_Tot = 0;
	fLambda_max = FindLambdaMax( FissilArray, HMMass );
	DBGV("fLambda_max "<<fLambda_max);

	IsotopicVector Fertile;
	IsotopicVector Fissile;
	double AvailablePuMass=0;
	double PuMassNeeded=1000; //At present time I have no clue what is the requiered Pu mass, I assume at least 1kg is needed
	double WeightPuContent=0;
	double MassPrecision=100; //Mass precision is 100 grams

	while(  fabs(PuMassNeeded - AvailablePuMass) > MassPrecision )
	{	//Increase the portion of the stock stokID taken, according to the followings variables
		double DeltaM=PuMassNeeded-AvailablePuMass;
		GuessLambda( lambda,0,FissilArray.size()-1, DeltaM, FissilArray ,HMMass);
		if(lambda[0]==-1)
		{	
			for(int i=0 ; i < (int)FissilArray.size() + (int)FertilArray.size();i++ )
				lambda[i ]= -1;
			WARNING<<"Not enought fissile material to build fuel"<<endl;
			return lambda;
		}
		//Build the Plutonium vector from stocks
		Fissile.Clear();
		
		for( int i = 0 ; i < (int)FissilArray.size() ; i++ )
			Fissile +=lambda[i] * FissilArray[i];

		AvailablePuMass = Fissile.GetTotalMass() * 1e6; //in grams
		//Building complementary Fertile from stocks
		double FertilMassNeeded = HMMass - AvailablePuMass;			
		double LAMBDA_FERTILE = FindLambdaMax(FertilArray, FertilMassNeeded);
		SetLambda(lambda,FissilArray.size(), lambda.size()-1,LAMBDA_FERTILE);
		int j=-1;
		Fertile.Clear();
		for(int i = (int)FissilArray.size() ; i < (int)FissilArray.size()+(int)FertilArray.size() ; i++)
		{	j++;
			Fertile +=  lambda[i] * FertilArray[j];
		}
	
		if(  fabs(Fertile.GetTotalMass()*1e6 - FertilMassNeeded) > FertilMassNeeded * 1e-6)//Not enought fertile in stocks
		{
			for( int i = 0 ; i < (int)FissilArray.size() + (int)FertilArray.size() ; i++ )
				lambda[i] = -1; //error code (read by FabricationPlant)
			WARNING<<"Not enought fertile material to build fuel"<<endl;
			return lambda;
		}
		/*Calcul the quantity of this composition needed to reach the burnup*/
		double MolarPuContent = GetFissileMolarFraction(Fissile, Fertile, BurnUp);
		double MeanMolarPu = Fissile.MeanMolar();
		double MeanMolarDepletedU = Fertile.MeanMolar();
		double MeanMolar   = MeanMolarPu * MolarPuContent + (1-MolarPuContent)*MeanMolarDepletedU;
		DBGV("MolarPuContent "<<MolarPuContent);
		WeightPuContent = MolarPuContent * MeanMolarPu / MeanMolar ;
		PuMassNeeded = WeightPuContent  *  HMMass ;
	}


	DBGV("Weight percent fissil : "<<PuMassNeeded/HMMass );
	DBGV("Lambda vector: ");
	for(int i = 0 ; i < (int)FissilArray.size() + (int)FertilArray.size() ; i++ )
			DBGV(lambda[i]);
	
DBGL
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

		for(int i = FirstStockID ; i <= LastStockID ; i++) //set to 0 all non touched value (to be sure)
			lambda[i] = 0  ;

		int PartieEntier = floor( LAMBDA_TOT );
		double PartieDecimal = LAMBDA_TOT - PartieEntier;

		for(int i=FirstStockID  ; i < FirstStockID +PartieEntier ;i++ )
				lambda[i]=1;

		lambda[FirstStockID + PartieEntier] = PartieDecimal  ;

}
//________________________________________________________________________
void EquivalenceModel::GuessLambda(vector<double>& lambda,int FirstStockID, int LastStockID, double DeltaM, vector<IsotopicVector> Stocks, double  HMMass)
{
	double	LAMBDA_TOT=0;
	for (int i = FirstStockID; i <=LastStockID ; i++)
		LAMBDA_TOT+=lambda[i]	;


	if( LAMBDA_TOT == 0 ) //Initialization Looking for lambda such as the fissile content is GetBuildFuelFirstGuess() of HM
	{	//cout<<"INITIALIZATION"<<endl;
		double	MASS = 0 ; 
		int ID_max = 0;
		while( MASS < GetBuildFuelFirstGuess()*HMMass )
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
				{	lambda[0]=-1;//error code;
					break;
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
		fOld_Lambda_Tot = LAMBDA_TOT;
		LAMBDA_TOT += (fLambda_max - LAMBDA_TOT)/2.;

		if(LAMBDA_TOT > 0.99*fLambda_max) //if we get close to the total of the stocks
		{	
			double MasseTot=0;
			for(int i=0;i<(int)Stocks.size();i++) 
				MasseTot+=Stocks[i].GetTotalMass()*1e6;

			if( (fLambda_max -LAMBDA_TOT )*MasseTot < 5 ) //If it remains less than 5gram in stocks => the stocks are not enought
				lambda[0]=-1;//error code;
		}
		//cout<<"DM>0 LAMBDA_TOT "<<LAMBDA_TOT<<endl;
		else
			SetLambda(lambda,FirstStockID,LastStockID,LAMBDA_TOT );
	}

	else  if( DeltaM < 0) 
	{
		//cout <<"fOld_Lambda_Tot[f] "<<fOld_Lambda_Tot[f]<<" LAMBDA_TOT "<<LAMBDA_TOT<< endl;
		LAMBDA_TOT = LAMBDA_TOT - fabs(fOld_Lambda_Tot - LAMBDA_TOT) / 2. ;
		//cout<<"DM<0 LAMBDA_TOT "<<LAMBDA_TOT<<endl;
		SetLambda(lambda,FirstStockID,LastStockID,LAMBDA_TOT );
	}
}
//________________________________________________________________________
double EquivalenceModel::FindLambdaMax(vector<IsotopicVector> Stocks, double  HMMass)
{	
	//je cherche les lambda telle que j'ai 100% de HMass  ou bien la mass total des stocks
	double LAMBDA = 0;
	double TotStockMass = 0;
	for (int i = 0 ; i < (int)Stocks.size() ; i++)
		TotStockMass += Stocks[i].GetTotalMass() * 1e6;
		

	if(TotStockMass <= HMMass)
		return Stocks.size() ;

	else
	{	
		double RemainHM=HMMass;
		for(int i=0 ;i<(int)Stocks.size();i++ )
		{
			if( RemainHM - Stocks[i].GetTotalMass() * 1e6 > 0)
			{
				RemainHM -= Stocks[i].GetTotalMass() * 1e6;
				LAMBDA++;
			}
			else
			{	LAMBDA+= RemainHM / (Stocks[i].GetTotalMass() * 1e6);
				return LAMBDA;
			}	
		}
	}

return -1;
}