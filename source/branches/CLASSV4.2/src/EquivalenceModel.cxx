#include "EquivalenceModel.hxx"
#include "StringLine.hxx"
#include "CLASSMethod.hxx"

EquivalenceModel::EquivalenceModel():CLASSObject()
{
	fRelativMassPrecision = 1/10000.; //Mass precision
	fMaxInterration = 10000; // Max iterration in build fueld algorythum
	freaded = false;

}

EquivalenceModel::EquivalenceModel(CLASSLogger* log):CLASSObject(log)
{
	fRelativMassPrecision = 1/10000.; //Mass precision
	fMaxInterration = 10000; // Max iterration in build fueld algorythm
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
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_reactor",	& EquivalenceModel::ReadType));
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fuel",	& EquivalenceModel::ReadType));
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fissil",	& EquivalenceModel::ReadFissil));
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fertil",	& EquivalenceModel::ReadFertil));
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_specpower",	& EquivalenceModel::ReadSpecificPower));
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fertil",	& EquivalenceModel::ReadMaximalContent));
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
	if( keyword ==  "k_fuel" )
		fDBFType = StringLine::NextWord(line, pos, ' ');
	else if( keyword ==  "k_reactor" )
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
void EquivalenceModel::ReadFissil(const string &line)
{
	DBGL
	int pos = 0;
	IsotopicVector FissileList;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_fissil" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_fissil\" not found !" << endl;
		exit(1);
	}
	
	int Z = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	double Q = atof(StringLine::NextWord(line, pos, ' ').c_str());
	
	FissileList.Add(Z, A, I, Q);
	
	DBGL
}


//________________________________________________________________________
void EquivalenceModel::ReadSpecificPower(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_specpower")	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_specpower\" Not found !" << endl;
		exit(1);
	}
	
	fSpecificPower = atof(StringLine::NextWord(line, pos, ' ').c_str());
	
	DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadMaximalContent(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' ').c_str());
	if( keyword != "k_contentmax")	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_contentmax\" Not found !" << endl;
		exit(1);
	}
	
	fMaximalContent = atof(StringLine::NextWord(line, pos, ' ').c_str());
	
	DBGL
}


//________________________________________________________________________
void EquivalenceModel::ReadFertil(const string &line)
{
	DBGL
	int pos = 0;
	IsotopicVector FertileList;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_fertil" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_fertil\" not found !" << endl;
		exit(1);
	}
	
	int Z 	= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A 	= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I 	= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	
	FertileList.Add(Z, A, I, 1.0);

	
	DBGL
}

//________________________________________________________________________
void EquivalenceModel::StocksTotalMassCalculation(map < string , vector <IsotopicVector> >  Stocks)
{
	// Calculating total mass of stock once and for all
	
	double TotalMassInStocks = 0;
	map < string , vector <IsotopicVector> >::iterator it_s_vIV;

	for(  it_s_vIV = Stocks.begin();   it_s_vIV != Stocks.end();  it_s_vIV++)
	{	
		TotalMassInStocks = 0;
		for(int i=0; i<Stocks[(* it_s_vIV).first].size(); i++)
		{
			TotalMassInStocks  +=  Stocks[(* it_s_vIV).first][i].GetTotalMass();
		}	
		fTotalMassInStocks[(* it_s_vIV).first] = TotalMassInStocks * 1e6; // in grams
	}
}

//________________________________________________________________________
double EquivalenceModel::LambdaCalculation(string MaterialDenomination, double MaterialMassNeeded, vector <IsotopicVector>  Stocks)
{
	double Lambda_tot = 0; 

	// If there is not enough matter in stocks construction fails
	if( MaterialMassNeeded > fTotalMassInStocks[MaterialDenomination] )
	{
		WARNING << "Not enough " << MaterialDenomination << " material to build fuel" << endl;
		WARNING << fTotalMassInStocks[MaterialDenomination] << endl;
		Lambda_tot = -1;
	}
	
	// Lambda calculation
	else
	{
		for( int i=0; i<Stocks.size(); i++)
		{	
			if( MaterialMassNeeded >= (Stocks[i].GetTotalMass()*1e6))
			{
				Lambda_tot +=  1;
				MaterialMassNeeded -=  (Stocks[i].GetTotalMass()*1e6);
			}
			else
			{
				Lambda_tot +=  MaterialMassNeeded/(Stocks[i].GetTotalMass()*1e6);
				break;
			}
		}
	}
	return Lambda_tot;
}
//________________________________________________________________________
map <string , vector<double> > EquivalenceModel::BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray)
{
DBGL

	map <string , vector<double> > lambda ; // map containing name of the list and associated vector of proportions taken from stocks
	
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
	
	/**** Some initializations **/
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		fTotalMassInStocks[(*it_s_vIV).first] = 0;

	StocksTotalMassCalculation(StreamArray);

	fActualMassContentInFuel	= GetBuildFuelFirstGuess();	
	int loopCount 			= 0;
	bool FuelBuiltCorrectly	 	= false ; 

	map <string , IsotopicVector > IVStream;
	map <string , double > MaterialMassNeeded ;
	map <string , double > LambdaNeeded;
	
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		LambdaNeeded[(*it_s_vIV).first] = 0;

	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		MaterialMassNeeded[(*it_s_vIV).first] = HMMass*fActualMassContentInFuel[(*it_s_vIV).first]; 

	do
	{	
		BreakReturnLambda = false;
		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
			LambdaNeeded[(*it_s_vIV).first] = LambdaCalculation((*it_s_vIV).first, MaterialMassNeeded[(*it_s_vIV).first], StreamArray[(*it_s_vIV).first]);

		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
		{		
			if( LambdaNeeded[(*it_s_D).first] == -1 )
			{
				SetLambdaToErrorCode(lambda[(*it_s_D).first]);
				WARNING << "Not enough : "<< (*it_s_D).first <<" to build fuel." << endl;	
				BreakReturnLambda = true; 				
			}
		}
		
		if(BreakReturnLambda) { return lambda;}
		
		BreakReturnLambda = false; 
		
		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
			SetLambda(lambda[(*it_s_D).first], LambdaNeeded[(*it_s_D).first]); 
	
		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
			IVStream[(*it_s_vIV).first].Clear();

		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		{	
			for(int i=0; i<StreamArray[(*it_s_vIV).first].size(); i++)
			{
				IVStream[(*it_s_vIV).first]  +=  lambda[(*it_s_vIV).first][i] * StreamArray[(*it_s_vIV).first][i];	
			}
		}

		if (loopCount > fMaxInterration)
		{
			ERROR << "Too much iterration in BuildFuel Method !";
			ERROR << "Need improvement in fuel fabrication ! Ask for it or D.I.Y. !!" << endl;
			exit(1);
		}
		
		/* Calcul the quantity of this composition needed to reach the burnup */
		map < string , double>  MaterialMolarContent = GetMolarFraction(IVStream, BurnUp);

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
		
		map < string , double >  MaterialMassContent;
		map < string , double >  MaterialMeanMolar;
		map < string , double >  MaterialMassAvailable;
		map < string , bool > 	   CheckOnMass;

		double FuelMeanMolar   = 0;

		for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
			MaterialMassAvailable[(*it_s_IV).first] = IVStream[(*it_s_IV).first].GetTotalMass()*1e06; 

		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
		{
			MaterialMeanMolar[(*it_s_D).first] = IVStream[(*it_s_D).first].GetMeanMolarMass();
			FuelMeanMolar +=  MaterialMolarContent[(*it_s_D).first] * MaterialMeanMolar[(*it_s_D).first];
		}

		for( it_s_D = MaterialMolarContent.begin();  it_s_D != MaterialMolarContent.end(); it_s_D++)
			MaterialMassContent [(*it_s_D).first] =  MaterialMolarContent[(*it_s_D).first] * MaterialMeanMolar[(*it_s_D).first] / FuelMeanMolar; 

		fActualMolarContentInFuel = MaterialMolarContent; //fActualContent can be accessed by a derivated EquivalenModel to accelerate GetFissileMolarFraction function (exemple in EQM_MLP_Kinf)
		
		for( it_s_D = MaterialMassContent.begin();  it_s_D != MaterialMassContent.end(); it_s_D++)
			MaterialMassNeeded[(*it_s_D).first] = HMMass*MaterialMassContent[(*it_s_D).first]; 

		for( it_s_D = MaterialMassNeeded.begin();  it_s_D != MaterialMassNeeded.end(); it_s_D++)
		{
			cout<<(*it_s_D).first<<"  "<<MaterialMassNeeded[(*it_s_D).first]<<"  "<<MaterialMassAvailable[(*it_s_D).first]<<" ";
			double DeltaM = fabs(MaterialMassNeeded[(*it_s_D).first] - MaterialMassAvailable[(*it_s_D).first]) / HMMass ; 
			if(DeltaM<fRelativMassPrecision) {CheckOnMass[(*it_s_D).first] = true;}
			else{CheckOnMass[(*it_s_D).first] = false;}
		}
		cout<<endl;
		FuelBuiltCorrectly = true; 
		for( it_s_B = CheckOnMass.begin();  it_s_B != CheckOnMass.end(); it_s_B++)
		{
			if(!CheckOnMass[(*it_s_B).first]) {FuelBuiltCorrectly = false;}
		}

		loopCount++;

	}while(!FuelBuiltCorrectly);
	
	for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
		(*this).isIVInDomain(IVStream[(*it_s_IV).first]);
	
	for( it_s_D = MaterialMassNeeded.begin();  it_s_D != MaterialMassNeeded.end(); it_s_D++)
		DBGV( "Weight percent : "<<(*it_s_D).first<<" "<< MaterialMassNeeded[(*it_s_D).first]/HMMass);
	
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
void EquivalenceModel::SetLambda(vector<double>& lambda, double Lambda_tot)
{
	if(Lambda_tot> lambda.size() )
	{
		ERROR << " FATAL ERROR " <<endl;
		exit(0);
	}

	for(int i = 0 ; i < lambda.size() ; i++) //set to 0 all non touched value (to be sure)
		lambda[i] = 0  ;

	int IntegerPart 		= floor( Lambda_tot );
	double DecimalPart 	= Lambda_tot - IntegerPart;

	for(int i=0  ; i < IntegerPart; i++ )
		lambda[i]=1;

	lambda[IntegerPart] = DecimalPart;
}

//________________________________________________________________________
void EquivalenceModel::SetLambdaToErrorCode(vector<double>& lambda)
{
	for( int i=0; i<lambda.size(); i++)
	{	
		lambda[i] = -1;	
	}
}

//________________________________________________________________________
bool EquivalenceModel::isIVInDomain(IsotopicVector IV)
{
	DBGL
	bool IsInDomain = true;
	
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
		for (map< ZAI,pair<double,double> >::iterator Domain_it = fZAILimits.begin(); Domain_it != fZAILimits.end(); Domain_it++)
		{
			double ThatZAIProp = IVNorm.GetIsotopicQuantity()[Domain_it->first];
			double ThatZAIMin  = Domain_it->second.first;
			double ThatZAIMax  = Domain_it->second.second;
			if( (ThatZAIProp > ThatZAIMax) || (ThatZAIProp <  ThatZAIMin) )
			{
				IsInDomain = false;
				
				WARNING << "Fresh fuel out of model range" << endl;
				WARNING << "\t AT LEAST this ZAI is accused to be outrange :" << endl;
				WARNING << "\t\t" << Domain_it->first.Z() << " " << Domain_it->first.A() << " " << Domain_it->first.I() << endl;
				WARNING << "\t\t min = " << ThatZAIMin  << " value = " << ThatZAIProp << " max = " << ThatZAIMax << endl;
				WARNING << "\t IV accused :" << endl << endl;
				WARNING << IVNorm.sPrint() << endl;
				break;
			}
		}
	}
	DBGL
	return IsInDomain;
	
}

