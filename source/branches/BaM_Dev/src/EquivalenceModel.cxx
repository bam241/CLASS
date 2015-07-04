#include "EquivalenceModel.hxx"
#include "StringLine.hxx"
#include "CLASSMethod.hxx"

EquivalenceModel::EquivalenceModel():CLASSObject()
{
	fRelativMassPrecision = 1/10000.; // Mass precision
	fMaxInterration = 100; // Max iterration in build fueld algorythum
	fFirstGuessFissilContent = 0.02;
	freaded = false;
	
}

EquivalenceModel::EquivalenceModel(CLASSLogger* log):CLASSObject(log)
{
	fRelativMassPrecision = 1/10000.; // Mass precision
	fMaxInterration = 100; // Max iterration in build fueld algorythm
	fFirstGuessFissilContent = 0.02;
	freaded = false;
	
}

EquivalenceModel::~EquivalenceModel()
{
	
}



void EquivalenceModel::ReadNFO()
{
	DBGL
	ifstream NFO(fInformationFile.c_str());
	
	if(!NFO)
	{
		ERROR << "Can't find/open file " << fInformationFile << endl;
		exit(0);
	}
	
	do
	{
		string line;
		getline(NFO,line);
		
		EquivalenceModel::ReadLine(line);
		
	} while(!NFO.eof());
	
	DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadLine(string line)
{
	DBGL
	
	if (!freaded)
	{
		int pos = 0;
		string keyword = tlc(StringLine::NextWord(line, pos, ' '));
		
		map<string, EQM_MthPtr>::iterator it = fKeyword.find(keyword);
		
		if(it != fKeyword.end())
			(this->*(it->second))( line );
		
		freaded = true;
		ReadLine(line);
		
	}
	
	freaded = false;
	
	DBGL
}


void EquivalenceModel::LoadKeyword()
{
	DBGL
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_zail",	& EquivalenceModel::ReadZAIlimits));
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_reactor",	& EquivalenceModel::ReadType)	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fuel",	& EquivalenceModel::ReadType)	 );
	DBGL
}


void EquivalenceModel::ReadType(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_fuel" && keyword != "k_reactor" )	// Check the keyword
	{
		ERROR << " Bad keyword : " << keyword << " Not found !" << endl;
		exit(1);
	}
	if( keyword == "k_fuel" )
		fDBFType = StringLine::NextWord(line, pos, ' ');
	else if( keyword == "k_reactor" )
		fDBRType = StringLine::NextWord(line, pos, ' ');
	
	DBGL
}


void EquivalenceModel::ReadZAIlimits(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_zail" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_zail\" not found !" << endl;
		exit(1);
	}
	
	int Z = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	
	double downLimit = atof(StringLine::NextWord(line, pos, ' ').c_str());
	double upLimit = atof(StringLine::NextWord(line, pos, ' ').c_str());
	
	if (upLimit < downLimit)
	{
		double tmp = upLimit;
		upLimit = downLimit;
		downLimit = tmp;
	}
	fZAILimits.insert(pair<ZAI, pair<double, double> >(ZAI(Z,A,I), pair<double,double>(downLimit, upLimit)));
	DBGL
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
	
	int j = -1;
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
vector<double> EquivalenceModel::BuildFuel(double BurnUp, double HMMass, vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray)
{
	
	DBGL
	vector<double> lambda ; // vector of portion of stocks taken (fissile & fertil)
	for(int i = 0 ; i < (int)FissilArray.size() + (int)FertilArray.size() ; i++ )
		lambda.push_back(0);
	
	/*** Test if there is a stock **/
	if( (int)FissilArray.size()==0 )
	{
		WARNING << " No fissile stocks available ! Fuel not build" << endl;
		SetLambdaToErrorCode(lambda);
		return lambda;
	}
	
	HMMass *= 1e6; // Unit onversion : tons to gram
	
	/**** Some initializations **/
	fTotalFissileMassInStocks = 0;
	fTotalFertileMassInStocks = 0;
	
	fActualFissileContent = GetBuildFuelFirstGuess();
	
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
			SetLambdaToErrorCode(lambda);
			WARNING << "Not enought fissile material to build fuel" << endl;
			return lambda;
		}
		
		SetLambda(lambda, 0, FissilArray.size()-1, LAMBDA_NEEDED );
		
		bool succeed = Build_Fuel_According_Lambda(lambda, FissilArray, FertilArray, HMMass, Fissile, Fertile);
		
		if(!succeed)
		{
			SetLambdaToErrorCode(lambda);
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
		
		if( MolarPuContent < 0 || MolarPuContent > 1 )
		{	SetLambdaToErrorCode(lambda);
			WARNING<<"GetFissileMolarFraction return negative or greater than one value";
			return lambda;
		}
		
		double MeanMolarPu = Fissile.GetMeanMolarMass();
		double MeanMolarDepletedU = Fertile.GetMeanMolarMass();
		
		double MeanMolar = MeanMolarPu * MolarPuContent + (1-MolarPuContent)  * MeanMolarDepletedU;
		
		
		WeightPuContent = MolarPuContent * MeanMolarPu / MeanMolar;
		fActualFissileContent = MolarPuContent; //fActualFissileContent can be accessed by a derivated EquivalenModel to accelerate GetFissileMolarFraction function (exemple in EQM_PWR_MLP_Kinf)
		PuMassNeeded = WeightPuContent  *  HMMass ;
		
		DBGV( "MolarPuContent " << MolarPuContent << " DeltaM " << PuMassNeeded - AvailablePuMass << " g" );
		
		loopCount++;
		
	}while(  fabs( PuMassNeeded - AvailablePuMass )/HMMass > fRelativMassPrecision );
	
	(*this).isIVInDomain(Fissile);
	
	DBGV( "Weight percent fissile : " << PuMassNeeded/HMMass );
	DBGV( "Lambda vector : " );
	
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
void EquivalenceModel::SetLambdaToErrorCode(vector<double>& lambda)
{
	for(int i=0 ; i < (int)lambda.size() ;i++ )
		lambda[i]= -1;
}



//________________________________________________________________________
bool EquivalenceModel::isIVInDomain(IsotopicVector IV)
{
	DBGL
	bool IsInDomain=true;
	
	if(fZAILimits.empty())
	{
	 WARNING << "Fresh Fuel variation domain is not set" << endl;
	 WARNING << "CLASS has no clue if the computed evolution for this fresh fuel is correct" << endl;
	 WARNING << "Proceed finger crossed !!" << endl;
	 return true;
	}
	
	else
	{
		IsotopicVector IVNorm = IV /IV.GetSumOfAll();
		for (map< ZAI,pair<double,double> >::iterator Domain_it=fZAILimits.begin(); Domain_it!=fZAILimits.end(); Domain_it++)
		{
			double ThatZAIProp = IVNorm.GetIsotopicQuantity()[Domain_it->first]	;
			double ThatZAIMin  = Domain_it->second.first;
			double ThatZAIMax  = Domain_it->second.second;
			if( (ThatZAIProp > ThatZAIMax) || (ThatZAIProp <  ThatZAIMin) )
			{
				IsInDomain = false;
				
				WARNING << "Fresh fuel out of model range" << endl;
				WARNING << "\t AT LEAST this ZAI is accused to be outrange :" << endl;
				WARNING << "\t\t" << Domain_it->first.Z() << " "<<Domain_it->first.A() << " " << Domain_it->first.I() << endl;
				WARNING << "\t\t min=" << ThatZAIMin <<" value=" << ThatZAIProp << " max=" << ThatZAIMax << endl;
				WARNING << "\t IV accused :" << endl << endl;
				WARNING << IVNorm.sPrint() << endl;
				break;
			}
		}
	}
	DBGL
	return IsInDomain;
	
}

