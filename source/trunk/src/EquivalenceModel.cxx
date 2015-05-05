#include "EquivalenceModel.hxx"


EquivalenceModel::EquivalenceModel():CLASSObject()
{
	fRelativMassPrecision = 1/10000.; // Mass precision
	fMaxInterration = 100; // Max iterration in build fueld algorythum
	fFirstGuessFissilContent = 0.02;
	
}

EquivalenceModel::EquivalenceModel(CLASSLogger* log):CLASSObject(log)
{
	fRelativMassPrecision = 1/10000.; // Mass precision
	fMaxInterration = 100; // Max iterration in build fueld algorythm
	fFirstGuessFissilContent = 0.02;
	
}

EquivalenceModel::~EquivalenceModel()
{
	
}
//________________________________________________________________________
double EquivalenceModel::LAMBDA_TOT_FOR(double MassNeeded, vector<IsotopicVector> Stocks, string FisOrFer)
{
	double Lambda_tot = 0;
	
	// Calculating total mass of stock once and for all
	if( fTotalFissileMassInStocks==0 || fTotalFertileMassInStocks==0 )
	{
		double TotalMassInStocks = 0;
		for( int i = 0 ; i<(int)Stocks.size() ; i++ )
			TotalMassInStocks += Stocks[i].GetTotalMass() ;
		
		if(FisOrFer == "Fis")
			fTotalFissileMassInStocks = TotalMassInStocks * 1e6; // in grams
		else
			fTotalFertileMassInStocks = TotalMassInStocks * 1e6; // in grams
	}
	
	double TotalMassInStocks = 0;
	
	if(FisOrFer == "Fis")
		TotalMassInStocks = fTotalFissileMassInStocks;
	else
		TotalMassInStocks = fTotalFertileMassInStocks;
	
	// If there is not enought matter in stocks construction fails
	if( MassNeeded > TotalMassInStocks )
	{
		WARNING<<"Not enought "<< FisOrFer <<" material to build fuel"<<endl;
		WARNING<<TotalMassInStocks<<endl;
		return -1;
	}
	
	for( int i = 0 ; i<(int)Stocks.size() ; i++ )
	{
		
		if( MassNeeded  >= (Stocks[i].GetTotalMass()*1e6) )
		{
			Lambda_tot+=1;
			MassNeeded -= (Stocks[i].GetTotalMass()*1e6);
		}
		else
		{
			Lambda_tot+=MassNeeded/(Stocks[i].GetTotalMass()*1e6);
			break;
		}
	}
	
	return Lambda_tot;
}
//________________________________________________________________________
bool EquivalenceModel::Build_Fuel_According_Lambda(vector<double> &lambda,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray, double HMMass, IsotopicVector &Fissile, IsotopicVector &Fertile)
{
	//Build the Plutonium vector from stocks
	Fissile.Clear();
	for( int i = 0; i < (int)FissilArray.size(); i++ )
		Fissile += lambda[i] * FissilArray[i];
	
	
	double AvailablePuMass = Fissile.GetTotalMass() * 1e6; // in grams
	
	// Building complementary Fertile from stocks
	double FertilMassNeeded = HMMass - AvailablePuMass;
	double LAMBDA_FERTILE = LAMBDA_TOT_FOR( FertilMassNeeded , FertilArray , "Fer" );
	
	SetLambda(lambda, (int)FissilArray.size(), (int)lambda.size()-1, LAMBDA_FERTILE);
	
	int j=-1;
	Fertile.Clear();
	for(int i = (int)FissilArray.size() ; i < (int)FissilArray.size()+(int)FertilArray.size() ; i++)
	{
		j++;
		Fertile +=  lambda[i] * FertilArray[j];
	}
	
	if(  fabs(Fertile.GetTotalMass()*1e6 - FertilMassNeeded) > FertilMassNeeded * 1e-6) // Not enought fertile in stocks
	{
		WARNING << "Not enought fertile material to build fuel" << endl;
		return false;
	}
	
	return true;
}
//________________________________________________________________________
vector<double> EquivalenceModel::BuildFuel(double BurnUp, double HMMass,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray)
{
	
	DBGL
	vector<double> lambda ; // vector of portion of stocks taken (fissile & fertil)
	
	
	/*** Test if there is a stock **/
	if( (int)FissilArray.size()==0 )
	{	WARNING << " No fissile stocks available ! Fuel not build" << endl;
		lambda[0] = -1;
		return lambda;
	}
	
	HMMass *= 1e6; // Unit onversion : tons to gram
	
	/**** Some initializations **/
	fTotalFissileMassInStocks = 0;
	fTotalFertileMassInStocks = 0;
	
	fActualFissileContent = GetBuildFuelFirstGuess();
	
	for(int i = 0 ; i < (int)FissilArray.size() + (int)FertilArray.size() ; i++ )
		lambda.push_back(0);
	
	IsotopicVector Fertile;
	IsotopicVector Fissile;
	
	double AvailablePuMass = 0;
	double PuMassNeeded = HMMass * fActualFissileContent;
	double WeightPuContent = 0;
	int loopCount = 0;
	
	do
	{
		double LAMBDA_NEEDED = LAMBDA_TOT_FOR(PuMassNeeded,FissilArray,"Fis");
		if( LAMBDA_NEEDED == -1 )	// Check if previous lambda was well calculated
		{
			for(int i=0 ; i < (int)FissilArray.size() + (int)FertilArray.size();i++ )
				lambda[i]= -1;
			
			WARNING << "Not enought fissile material to build fuel" << endl;
			
			return lambda;
		}
		
		SetLambda(lambda, 0, FissilArray.size()-1, LAMBDA_NEEDED );
		
		bool succeed = Build_Fuel_According_Lambda(lambda, FissilArray, FertilArray, HMMass, Fissile, Fertile);
		
		if(!succeed)
		{
			for( int i = 0; i < (int)FissilArray.size() + (int)FertilArray.size(); i++ )
				lambda[i] = -1; // Error code (readed by FabricationPlant)
			
			return lambda;
		}
		
		AvailablePuMass = Fissile.GetTotalMass() * 1e6; //in grams
		
		if (loopCount > fMaxInterration)
		{
			ERROR << "Too much iterration in BuildFuel Method !";
			ERROR << "Need improvement in fuel fabrication ! Ask for it or D.I.Y. !!" << endl;
			exit(1);
		}
		
		/* Calcul the quantity of this composition needed to reach the burnup */
		double MolarPuContent = GetFissileMolarFraction(Fissile, Fertile, BurnUp);
		
		double MeanMolarPu = Fissile.GetMeanMolarMass();
		double MeanMolarDepletedU = Fertile.GetMeanMolarMass();
		
		double MeanMolar = MeanMolarPu * MolarPuContent + (1-MolarPuContent)  *MeanMolarDepletedU;
		
		
		WeightPuContent = MolarPuContent * MeanMolarPu / MeanMolar;
		fActualFissileContent = MolarPuContent; //fActualFissileContent can be accessed by a derivated EquivalenModel to accelerate GetFissileMolarFraction function (exemple in EQM_MLP_Kinf)
		PuMassNeeded = WeightPuContent  *  HMMass ;
		
		DBGV( "MolarPuContent " << MolarPuContent << " DeltaM " << PuMassNeeded - AvailablePuMass << " g" );
		
		loopCount++;
		
	}while(  fabs( PuMassNeeded - AvailablePuMass )/HMMass > fRelativMassPrecision );
	
	
	DBGV( "Weight percent fissil : " << PuMassNeeded/HMMass );
	DBGV( "Lambda vector: " );
	for(int i = 0; i < (int)FissilArray.size() + (int)FertilArray.size(); i++ )
		DBGV(lambda[i]);
	
	return lambda;
}
//________________________________________________________________________
void EquivalenceModel::SetLambda(vector<double>& lambda ,int FirstStockID, int LastStockID, double LAMBDA_TOT)
{
	if( LAMBDA_TOT > LastStockID - FirstStockID + 1 )
	{
		ERROR << " FATAL ERROR " << endl;
		exit(0);
	}
	
	for( int i = FirstStockID; i <= LastStockID; i++ ) //set to 0 all non touched value (to be sure)
		lambda[i] = 0  ;
	
	int IntegerPart = floor( LAMBDA_TOT );
	double DecimalPart = LAMBDA_TOT - IntegerPart;
	
	for( int i=FirstStockID; i < FirstStockID +IntegerPart; i++ )
		lambda[i] = 1;
	
	lambda[FirstStockID + IntegerPart] = DecimalPart;
}
//________________________________________________________________________
