#include "EquivalenceModel.hxx"
#include "external/StringLine.hxx"
#include "CLASSMethod.hxx"

//________________________________________________________________________
EquivalenceModel::EquivalenceModel():CLASSObject()
{
	fRelativMassPrecision 	= 5/10000.; 	//Mass precision
	fMaxInterration 	= 10000; 	// Max iterration in build fueld algorythum
	freaded 		= false;
	
	EquivalenceModel::LoadKeyword();
}
//________________________________________________________________________
EquivalenceModel::EquivalenceModel(CLASSLogger* log):CLASSObject(log)
{
	fRelativMassPrecision 	= 5/10000.; 	//Mass precision
	fMaxInterration 	= 10000; 	// Max iterration in build fueld algorythm
	freaded 		= false;
	
	EquivalenceModel::LoadKeyword();
}
//________________________________________________________________________
EquivalenceModel::~EquivalenceModel()
{

}
//________________________________________________________________________
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
//________________________________________________________________________
void EquivalenceModel::LoadKeyword()
{
	DBGL
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_zail",			& EquivalenceModel::ReadZAIlimits)	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_reactor",			& EquivalenceModel::ReadType)	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fuel",			& EquivalenceModel::ReadType)	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_firstguesscontent",	& EquivalenceModel::ReadFirstGuessContent) 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_list",			& EquivalenceModel::ReadList) 	 		 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_specpower",		& EquivalenceModel::ReadSpecificPower)	 );
	DBGL
}
//________________________________________________________________________
void EquivalenceModel::PrintInfo()
{
	INFO << "Reactor Type : "<< fDBRType << endl;
	INFO << "Fuel Type : "<< fDBFType << endl;
	INFO << "Specific Power [W/g]: "<< fSpecificPower << endl;
	
	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;

	for(  it_s_IV = fStreamList.begin();   it_s_IV != fStreamList.end();  it_s_IV++)
	{	
		INFO <<(* it_s_IV).first<<"  (Z A I) :" << endl;
		map<ZAI ,double >::iterator it1;
		map<ZAI ,double > fMap1 = fStreamList[(* it_s_IV).first].GetIsotopicQuantity();
		for(it1 = fMap1.begin()  ; it1 != fMap1.end() ; it1++)
			INFO << (*it1).first.Z() <<" "<< (*it1).first.A() <<" "<< (*it1).first.I() << endl;
	}
	INFO<<"First guess content in fuel : "<<endl;
	for(  it_s_D = fFirstGuessContent.begin();   it_s_D != fFirstGuessContent.end();  it_s_D++)
	{
		INFO <<(* it_s_D).first<<" "<<fFirstGuessContent[(* it_s_D).first]<<endl;
	}	

	INFO<<"ZAI limits (validity domain)[prop in fresh fuel] (Z A I min max) :"<<endl;
	for (map< ZAI,pair<double,double> >::iterator Domain_it = fZAILimits.begin(); Domain_it != fZAILimits.end(); Domain_it++)
	{	
		double ThatZAIMin  = Domain_it->second.first;
		double ThatZAIMax  = Domain_it->second.second;
		int Z = Domain_it->first.Z();
		int A = Domain_it->first.A();
		int I = Domain_it->first.I();

		INFO <<ThatZAIMin<<" < ZAI ("<< Z<< " " << A << " " << I<<")"<<" < "<<ThatZAIMax<< endl;
	}

}
//________________________________________________________________________
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
//________________________________________________________________________
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
	
	int Z 	= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A 	= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I 	= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	
	double downLimit 	= atof(StringLine::NextWord(line, pos, ' ').c_str());
	double upLimit 		= atof(StringLine::NextWord(line, pos, ' ').c_str());
	
	if (upLimit < downLimit)
	{
		double tmp 	= upLimit;
		upLimit 	= downLimit;
		downLimit 	= tmp;
	}
	fZAILimits.insert(pair<ZAI, pair<double, double> >(ZAI(Z,A,I), pair<double,double>(downLimit, upLimit)));
	DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadList(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_list" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_list\" not found !" << endl;
		exit(1);
	}
	string ListName	= StringLine::NextWord(line, pos, ' ');
	int Z 		= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A 		= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I 		= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	double Q 	= atof(StringLine::NextWord(line, pos, ' ').c_str());
	fStreamList[ListName].Add(Z, A, I, Q);
	
	DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadFirstGuessContent(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_firstguesscontent" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_firstguesscontent\" not found !" << endl;
		exit(1);
	}
	string ListName	= StringLine::NextWord(line, pos, ' ');
	double Q 	= atof(StringLine::NextWord(line, pos, ' ').c_str());
	fFirstGuessContent[ListName] = Q;
	
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
void EquivalenceModel::StocksTotalMassCalculation(map < string , vector <IsotopicVector> >  Stocks)
{
	// Calculating total mass of stock once and for all
	
	
	double TotalMassInStocks = 0;
	map < string , vector <IsotopicVector> >::iterator it_s_vIV;

	for( it_s_vIV = Stocks.begin();  it_s_vIV != Stocks.end(); it_s_vIV++)
	{
		fTotalMassInStocks[(*it_s_vIV).first] = 0;
		fLambdaMax[(*it_s_vIV).first] = 0;

	}

	for(  it_s_vIV = Stocks.begin();   it_s_vIV != Stocks.end();  it_s_vIV++)
	{	
		TotalMassInStocks = 0;
		for(int i=0; i < (int)Stocks[(* it_s_vIV).first].size(); i++)
		{
			TotalMassInStocks  +=  Stocks[(* it_s_vIV).first][i].GetTotalMass();
		}
		fLambdaMax[(*it_s_vIV).first] = Stocks[(* it_s_vIV).first].size();	
		fTotalMassInStocks[(* it_s_vIV).first] = TotalMassInStocks * 1e6; // in grams
	}
}

//________________________________________________________________________
double EquivalenceModel::LambdaCalculation(string MaterialDenomination, double LambdaPreviousStep, double MaterialMassNeeded, double DeltaMass, vector <IsotopicVector>  Stocks)
{
	double Lambda_tot = 0; 

	// If there is not enough matter in stocks construction fails
	if( MaterialMassNeeded > fTotalMassInStocks[MaterialDenomination] )
	{
		if(DeltaMass > 0)
			Lambda_tot = LambdaPreviousStep - (fLambdaMax[MaterialDenomination] - LambdaPreviousStep)/2 ; 
		
		if(DeltaMass < 0)
			Lambda_tot = LambdaPreviousStep + (fLambdaMax[MaterialDenomination] - LambdaPreviousStep)/2 ; 
		
		if( (fLambdaMax[MaterialDenomination] - Lambda_tot)<1 && (fLambdaMax[MaterialDenomination]-Lambda_tot)*Stocks.back().GetTotalMass()*1e6<MaterialMassNeeded*fRelativMassPrecision/2.)
		{
			WARNING << "Not enough " << MaterialDenomination << " material to build fuel" << endl;
			WARNING << "Mass available "<<fTotalMassInStocks[MaterialDenomination] << endl;
			WARNING << "Mass needed "<<MaterialMassNeeded<< endl;
			Lambda_tot = -1;
		}
	}
	
	// Lambda calculation
	else
	{
		for( int i=0; i < (int)Stocks.size(); i++)
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
void EquivalenceModel::SetLambda(vector<double>& lambda, double Lambda_tot)
{
	if(Lambda_tot > (int)lambda.size() )
	{
		cout<<Lambda_tot<<"  "<<lambda.size()<<endl;
		ERROR << " FATAL ERROR " <<endl;
		exit(0);
	}

	for(int i = 0 ; i < (int)lambda.size() ; i++) //set to 0 all non touched value (to be sure)
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
	for( int i=0; i < (int)lambda.size(); i++)
	{	
		lambda[i] = -1;	
	}
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
		for(int i=0; i < (int)StreamArray[(*it_s_vIV).first].size(); i++)
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

	StocksTotalMassCalculation(StreamArray);

	fActualMassContentInFuel	= GetBuildFuelFirstGuess();
	fActualMolarContentInFuel	= GetBuildFuelFirstGuess();	
	int loopCount 			= 0;
	bool FuelBuiltCorrectly	 	= false ; 

	map <string , IsotopicVector > IVStream;
	map < string , double > MaterialMassNeeded ;
	map < string , double > LambdaNeeded;
	map < string , double > DeltaMass;
	
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
		LambdaNeeded[(*it_s_vIV).first] = 0;

	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		MaterialMassNeeded[(*it_s_vIV).first] = HMMass*fActualMassContentInFuel[(*it_s_vIV).first]; 
		DeltaMass[(*it_s_vIV).first] =  - MaterialMassNeeded[(*it_s_vIV).first];
	}
	do
	{	
		map < string , double >  LambdaPreviousStep;

		for( it_s_D = LambdaNeeded.begin();  it_s_D != LambdaNeeded.end(); it_s_D++)
			LambdaPreviousStep[(*it_s_D).first] =LambdaNeeded[(*it_s_D).first]; 
		
		BreakReturnLambda = false;
		for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
			LambdaNeeded[(*it_s_vIV).first] = LambdaCalculation((*it_s_vIV).first, LambdaPreviousStep[(*it_s_vIV).first], MaterialMassNeeded[(*it_s_vIV).first], DeltaMass[(*it_s_vIV).first], StreamArray[(*it_s_vIV).first]);

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
			for(int i=0; i < (int)StreamArray[(*it_s_vIV).first].size(); i++)
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
		map < string , bool > 	  CheckOnMass;

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
		{	
			MaterialMassNeeded[(*it_s_D).first] = HMMass*MaterialMassContent[(*it_s_D).first]; 
			DeltaMass[(*it_s_D).first] = MaterialMassAvailable[(*it_s_D).first] - MaterialMassNeeded[(*it_s_D).first];
		}
		
		for( it_s_D = MaterialMassNeeded.begin();  it_s_D != MaterialMassNeeded.end(); it_s_D++)
		{
			double DeltaM = fabs(MaterialMassNeeded[(*it_s_D).first] - MaterialMassAvailable[(*it_s_D).first]) / HMMass ;
			
			if(DeltaM<fRelativMassPrecision) {CheckOnMass[(*it_s_D).first] = true;}
			else{CheckOnMass[(*it_s_D).first] = false;}
		}
	
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

		for(int i=0; i < (int)lambda[(*it_s_vD).first].size(); i++)
		{
			DBGV(lambda[(*it_s_vD).first][i]); 
		}
	}		
	return lambda;
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

