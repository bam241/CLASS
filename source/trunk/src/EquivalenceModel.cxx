#include "EquivalenceModel.hxx"



EquivalenceModel::EquivalenceModel():CLASSObject()
{
	fRelativMassPrecision = 1/10000.; //Mass precision
	fMaxInterration = 1000; // Max iterration in build fueld algorythum
	fFirstGuessFissilContent = 0.04;

}

EquivalenceModel::EquivalenceModel(CLASSLogger* log):CLASSObject(log)
{
	fRelativMassPrecision = 1/10000.; //Mass precision
	fMaxInterration = 1000; // Max iterration in build fueld algorythm
	fFirstGuessFissilContent = 0.04;

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
	fActualFissileContent = GetBuildFuelFirstGuess(); //usefull for EQM_MLP_Kinf

	vector<double> lambda ; //vector of portion of stocks taken (fissile & fertil)
	for(int i = 0 ; i < (int)FissilArray.size() + (int)FertilArray.size() ; i++ )
		lambda.push_back(0);


	if( (int)FissilArray.size()==0 )
	{	WARNING<<" No fissile stocks available ! Fuel not build"<<endl;
		lambda[0] = -1;
		return lambda;
	}


	fLambda_max = FindLambdaMax( FissilArray, HMMass );
	fOld_Lambda_Tot_Minus = 0;
	fOld_Lambda_Tot_Plus  = fLambda_max;
	DBGV("fLambda_max "<<fLambda_max);

	IsotopicVector Fertile;
	IsotopicVector Fissile;
	double AvailablePuMass = 0;
	double PuMassNeeded = HMMass * 0.01; //At present time I have no clue what is the requiered Pu mass, I assume at least 1kg is needed
	double WeightPuContent = 0;

	int loopCount = 0;
	
	do
	{
		//Increase the portion of the stock stokID taken, according to the followings variables
		double DeltaM=PuMassNeeded-AvailablePuMass;
		GuessLambda( lambda, 0, (int)FissilArray.size()-1, DeltaM, FissilArray ,HMMass);
		
		if(lambda[0]==-1)	// Check if previous lambda was well calculated
		{	
			for(int i=0 ; i < (int)FissilArray.size() + (int)FertilArray.size();i++ )
				lambda[i]= -1;
			
			WARNING<<"Not enought fissile material to build fuel"<<endl;
			return lambda;
		}
		else if (loopCount > fMaxInterration)
		{
			ERROR << "Too much iterration in BuildFuel Method ! Need improvement in fuel fabrication ! Ask for it or D.I.Y. !!" << endl;
			exit(1);
		}
		
		
		//Build the Plutonium vector from stocks
		Fissile.Clear();
		
		for( int i = 0 ; i < (int)FissilArray.size() ; i++ )
			Fissile +=lambda[i] * FissilArray[i];

		AvailablePuMass = Fissile.GetTotalMass() * 1e6; //in grams
		//Building complementary Fertile from stocks
		double FertilMassNeeded = HMMass - AvailablePuMass;			
		double LAMBDA_FERTILE = FindLambdaMax(FertilArray, FertilMassNeeded);
		
		SetLambda(lambda, (int)FissilArray.size(), (int)lambda.size()-1,LAMBDA_FERTILE);
		
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
		fActualFissileContent = MolarPuContent; //usefull for EQM_MLP_Kinf
		double MeanMolarPu = Fissile.MeanMolar();
		double MeanMolarDepletedU = Fertile.MeanMolar();
		double MeanMolar   = MeanMolarPu * MolarPuContent + (1-MolarPuContent)*MeanMolarDepletedU;
		
		DBGV("MolarPuContent "<<MolarPuContent<<" DeltaM "<<PuMassNeeded-AvailablePuMass << " g");
		
		WeightPuContent = MolarPuContent * MeanMolarPu / MeanMolar ;
		PuMassNeeded = WeightPuContent  *  HMMass ;
		loopCount++;
	}while(  fabs(PuMassNeeded - AvailablePuMass)/HMMass > fRelativMassPrecision );


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
	{
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
		}
	}

	else  if( DeltaM > 0) 
	{
		fOld_Lambda_Tot_Minus = LAMBDA_TOT;
		LAMBDA_TOT += (fOld_Lambda_Tot_Plus - LAMBDA_TOT)/2.;

		if( ( (fLambda_max - LAMBDA_TOT ) < 1)
			&& ( fLambda_max - LAMBDA_TOT ) * Stocks.back().GetTotalMass() * 1e6 < HMMass *fRelativMassPrecision/2.) //if we get close to the total of the stocks
				lambda[0]=-1;//error code;
		else
			SetLambda(lambda,FirstStockID,LastStockID,LAMBDA_TOT );
	}

	else  if( DeltaM < 0) 
	{	fOld_Lambda_Tot_Plus =  LAMBDA_TOT;
		LAMBDA_TOT = LAMBDA_TOT - fabs(fOld_Lambda_Tot_Minus - LAMBDA_TOT) / 2. ;
		SetLambda(lambda,FirstStockID,LastStockID,LAMBDA_TOT );
	}
	DBGV("LAMBDA_TOT :"<<LAMBDA_TOT);
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