#include "EquivalenceModel.hxx"






EquivalenceModel::EquivalenceModel():CLASSObject()
{

}

EquivalenceModel::EquivalenceModel(LogFile* log):CLASSObject(log)
{

}

EquivalenceModel::~EquivalenceModel()
{

}




//________________________________________________________________________
vector<double> EquivalenceModel::BuildFuel(double BurnUp, double HMMass,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray)
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
	 Fertile.Add(U5,   U5_enrich   *HMMass/MeanMolarDepletedU*AVOGADRO);
	 Fertile.Add(U8,  (1-U5_enrich)*HMMass/MeanMolarDepletedU*AVOGADRO);
	 */
	IsotopicVector Fertile;
	IsotopicVector Fissile;
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
		Fissile.Clear();
		for(int i=0;i<=StockID;i++)
			Fissile += lambda[i] * FissilArray[i];

		AvailablePuMass = Fissile.GetTotalMass() * 1e6; //in grams

		//Build uranium vector from stocks
		int FertileStockID = FissilArray.size();
		double FertilMassNeeded = HMMass - AvailablePuMass;
		double AvailableFertilMass = 0;

		while( abs(FertilMassNeeded - AvailableFertilMass) > MassPrecision  )
		{
			double DeltaM=FertilMassNeeded-AvailableFertilMass;
			GuessLambda( lambda[FertileStockID], FertileStockID,FissilArray.size(),FertilArray.size()-1, DeltaM, FertilArray[FertileStockID-FissilArray.size()].GetTotalMass() * 1e6 ) ;

			int j=-1;
			Fertile.Clear();
			for(int i=int(FissilArray.size());i<=FertileStockID;i++)
			{	Fertile += lambda[i] * FertilArray[j];
				j++;
			}
			AvailableFertilMass=Fertile.GetTotalMass() * 1e6; //in grams

			if( j+1 == int( FertilArray.size() ) && FertilMassNeeded  > AvailableFertilMass ) //if this is the last stock and it's not enought
			{
				cout<<"You requiere more depleted uranium "<<"("<<DeltaM/1e6<<" t needed) ! Reactor not fill"<<endl;
				for(int i=0 ; i<int(lambda.size()) ; i++)
					lambda[i]=0;
				break;
			}

		}
		/*Calcul the quantity of this composition needed to reach the burnup*/
		double MolarPuContent = GetFissileMolarFraction(Fissile, Fertile, BurnUp);

		double MeanMolarPu = Fissile.MeanMolar();
		double MeanMolarDepletedU = Fertile.MeanMolar();
		double MeanMolar   = MeanMolarPu*MolarPuContent + (1-MolarPuContent)*MeanMolarDepletedU;

		WeightPuContent = MolarPuContent * MeanMolarPu/MeanMolar ;

		PuMassNeeded = WeightPuContent  *  HMMass ;

		if( StockID+1 == int( FissilArray.size() ) && PuMassNeeded  > AvailablePuMass ) //if this is the last stock and it's not enought
		{
			cout<<"You requiere more (or better) plutonium !! Reactor not fill"<<endl;
			for(int i=0 ; i<int(lambda.size()) ; i++)
				lambda[i]=0;
			break;
		}
		
	}
	
	return lambda;
}



//________________________________________________________________________
void EquivalenceModel::GuessLambda(double& lambda, int& StockID,int FirstStockID, int LastStockID, double DeltaM,double StockHM)
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
			lambda = 1;

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