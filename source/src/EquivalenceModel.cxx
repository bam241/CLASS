#include "EquivalenceModel.hxx"
#include "external/StringLine.hxx"
#include "CLASSMethod.hxx"

//________________________________________________________________________
EquivalenceModel::EquivalenceModel():CLASSObject()
{
	fMaxIterration 		= 500; 		// Max iterration in build fueld algorythum
	freaded 		= false;
	
	EquivalenceModel::LoadKeyword();
}
//________________________________________________________________________
EquivalenceModel::EquivalenceModel(CLASSLogger* log):CLASSObject(log)
{
	fMaxIterration 		= 500; 		// Max iterration in build fueld algorythm
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
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_reactor",		& EquivalenceModel::ReadType)	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fuel",			& EquivalenceModel::ReadType)	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_massfractionmin",	& EquivalenceModel::ReadEqMinFraction) 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_massfractionmax",	& EquivalenceModel::ReadEqMaxFraction) 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_list",			& EquivalenceModel::ReadList) 	 	 );
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
	INFO<<"Minimum fraction in the fuel for each material : "<<endl;
	for(  it_s_D = fStreamListEqMMassFractionMin.begin();   it_s_D != fStreamListEqMMassFractionMin.end();  it_s_D++)
	{
		INFO <<(* it_s_D).first<<" "<<fStreamListEqMMassFractionMin[(* it_s_D).first]<<endl;
	}	

	INFO<<"Maximum fraction in the fuel for each material : "<<endl;
	for(  it_s_D = fStreamListEqMMassFractionMax.begin();   it_s_D != fStreamListEqMMassFractionMax.end();  it_s_D++)
	{
		INFO <<(* it_s_D).first<<" "<<fStreamListEqMMassFractionMax[(* it_s_D).first]<<endl;
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
	string ListName= StringLine::NextWord(line, pos, ' ');
	int Z 		= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A 		= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I 		= atoi(StringLine::NextWord(line, pos, ' ').c_str());
	double Q 	= atof(StringLine::NextWord(line, pos, ' ').c_str());
	fStreamList[ListName].Add(Z, A, I, Q);
	
	DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadEqMinFraction(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_massfractionmin" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_massfractionmin\" not found !" << endl;
		exit(1);
	}
	string ListName= StringLine::NextWord(line, pos, ' ');
	double Q 	 = atof(StringLine::NextWord(line, pos, ' ').c_str());
	fStreamListEqMMassFractionMin[ListName] = Q;

	DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadEqMaxFraction(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_massfractionmax" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_massfractionmax\" not found !" << endl;
		exit(1);
	}
	string ListName= StringLine::NextWord(line, pos, ' ');
	double Q 	 = atof(StringLine::NextWord(line, pos, ' ').c_str());
	fStreamListEqMMassFractionMax[ListName] = Q;

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
void EquivalenceModel::SetLambdaToErrorCode(vector<double>& lambda)
{

	if(lambda.size() == 0) //then we have to add an element to send the error code to the fab (case for no storage in stream)
	{
		lambda.push_back(-1);
	}	

	else // other errors (no enough material or too many steps)
	{
		for( int i=0; i < (int)lambda.size(); i++)
		{	
			lambda[i] = -1;	
		}
	}	

}
//________________________________________________________________________
void EquivalenceModel::StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks)
{
	// Calculating total mass of stock once and for all
	
	
	double TotalMassInStocks = 0;
	map < string , vector <IsotopicVector> >::const_iterator it_s_vIV;

	for( it_s_vIV = Stocks.begin();  it_s_vIV != Stocks.end(); it_s_vIV++)
	{
		fTotalMassInStocks[ it_s_vIV->first ] = 0;
		fLambdaMax[ it_s_vIV->first ] = 0;
	}

	for(  it_s_vIV = Stocks.begin();   it_s_vIV != Stocks.end();  it_s_vIV++)
	{	
		TotalMassInStocks = 0;
		for(int i=0; i < (int)Stocks.at((* it_s_vIV).first).size(); i++)
		{
			TotalMassInStocks  +=  Stocks.at( it_s_vIV->first )[i].GetTotalMass();
		}
		fLambdaMax[(*it_s_vIV).first] = Stocks.at( it_s_vIV->first ).size();	
		fTotalMassInStocks[ it_s_vIV->first ] = TotalMassInStocks * 1e6; // in grams
	}
}


//________________________________________________________________________
map <string , vector<double> > EquivalenceModel::BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray,  map < string , double> StreamListFPMassFractionMin, map < string , double> StreamListFPMassFractionMax, map < int , string > StreamListPriority, map < string , bool> StreamListIsBuffer)
{
DBGL

	map <string , vector<double> > lambda ; // map containing name of the list and associated vector of proportions taken from stocks
	
	//Iterators declaration
	map < string , vector  <IsotopicVector> >::iterator it_s_vIV;
	map < string , vector  <double> >::iterator it_s_vD;
	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	map < int , string >::iterator it_i_s;
	
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{	
		for(int i=0; i < (int)StreamArray[(*it_s_vIV).first].size(); i++)
		{
			lambda[(*it_s_vIV).first].push_back(0);
		}
	}	

	// Test if there is at least one stock available in each list, otherwise fuel is not built //
	bool BreakReturnLambda = false; 
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{
		if(StreamArray[(*it_s_vIV).first].size() == 0)
		{
			WARNING << " No stock available for stream : "<< (*it_s_vIV).first <<".  Fuel not built." << endl;
			SetLambdaToErrorCode(lambda[(*it_s_vIV).first]);
			BreakReturnLambda = true; 	
		}
	}
	if(BreakReturnLambda) { return lambda;}
	
	/// Search for the minimum and maximum fraction of each material in fuel ///
	map < string, double >   StreamListMassFractionMin ; 
	map < string, double >   StreamListMassFractionMax ; 
	for( it_s_D = StreamListFPMassFractionMin.begin();  it_s_D != StreamListFPMassFractionMin.end(); it_s_D++)
	{	
		if(StreamListFPMassFractionMin[(*it_s_D).first] < fStreamListEqMMassFractionMin[(*it_s_D).first]) // if limits FP are lower than limits EqM => limits Eqm are applied
		{
			ERROR << " User mass fraction min requirement is lower than the model mass fraction min for list  : "<<(*it_s_D).first << endl;
			ERROR << " User mass fraction min requirement : "<<StreamListFPMassFractionMin[(*it_s_D).first]<<endl;
			ERROR << " Model mass fraction min requirement : "<<fStreamListEqMMassFractionMin[(*it_s_D).first]<<endl;			
			exit(1);
		}
		else
		{
			StreamListMassFractionMin[(*it_s_D).first] = StreamListFPMassFractionMin[(*it_s_D).first];
		}
	}	

	for( it_s_D = StreamListFPMassFractionMax.begin();  it_s_D != StreamListFPMassFractionMax.end(); it_s_D++)
	{	
		if(StreamListFPMassFractionMax[(*it_s_D).first] > fStreamListEqMMassFractionMax[(*it_s_D).first]) // if limits FP are higher than limits EqM => limits Eqm are applied
		{
			ERROR << " User mass fraction max requirement is higher than the model mass fraction max for list  : "<<(*it_s_D).first << endl;
			ERROR << " User mass fraction max requirement : "<<StreamListFPMassFractionMax[(*it_s_D).first]<<endl;
			ERROR << " Model mass fraction max requirement : "<<fStreamListEqMMassFractionMax[(*it_s_D).first]<<endl;			
			exit(1);
		}
		else
		{
			StreamListMassFractionMax[(*it_s_D).first] = StreamListFPMassFractionMax[(*it_s_D).first];
		}

	}

	// Check if there is enough material in stock to satisfy mass fraction min //
	BreakReturnLambda = false; 
	StocksTotalMassCalculation(StreamArray);
	for( it_s_D = StreamListMassFractionMin.begin();  it_s_D != StreamListMassFractionMin.end(); it_s_D++)
	{
		if(fTotalMassInStocks[(*it_s_D).first]< HMMass*StreamListMassFractionMin[(*it_s_D).first])
		{
			WARNING << " Not enough material  : "<< (*it_s_D).first << " in stocks to reach the build fuel lower limit of "<<StreamListMassFractionMin[(*it_s_D).first]<<" reactor mass.  Fuel not built." << endl;
			SetLambdaToErrorCode(lambda[(*it_s_D).first]);
			BreakReturnLambda = true; 	
		}
	}
	if(BreakReturnLambda) { return lambda;}

	//Check if Mass Fractions Min/Max lead to target BU//

	HMMass *=  1e6; //Unit onversion : tons to gram
	

	cout<<"Ca fonctionne !!! "<<endl;
	exit(1);

/*	//Search in the sorted stream array the point where calculated BU is higher than targeted BU//

	bool BurnUpExceeded 		= false ;
	int BurnUpExceededPosition 	= 0;
	string BurnUpExceededList	= "";
	double HigherLimitOnBurnUp	= 0;
	int StreamListNumber 		= 0;	

	int j=0;

	IsotopicVector FuelToTestWithoutBuffer  		= IsotopicVector();
	IsotopicVector PreviousFuelToTestWithoutBuffer 	= IsotopicVector();

	for( it_i_s = StreamListPriority.begin();  it_i_s != StreamArray.end(); it_i_s++)
	{	
		if(!BurnUpExceeded && !StreamListIsBuffer[(*it_i_s).second])
		{	
			while(FractionMassToReachMin[(*it_i_s).second][j] != 0.0)
			{
				FuelToTestWithoutBuffer += StreamArray[(*it_i_s).second][j] *  FractionMassToReachMin[(*it_i_s).second][j];
				j++;
			}
			IsotopicVector Buffer 		= BuildBuffer(FuelToTestWithoutBuffer , HMMass, StreamArray) ;
			IsotopicVector FuelToTest 	= FuelToTestWithoutBuffer + Buffer ; 
			FuelToTest 			= FuelToTest/FuelToTest.GetSumOfAll();
				
			double EqMMaximumBurnUp	= GetMaximumBurnUp (FuelToTest, BurnUp) ;
			
			if(EqMMaximumBurnUp=>BurnUp)
			{
				BurnUpExceeded 		= true;
				BurnUpExceededPosition	= j;
				HigherLimitOnBurnUp		= EqMMaximumBurnUp;
				BurnUpExceededList 		= (*it_i_s).second;
				break;
			}
			
			if(FractionMassToReachMin[(*it_i_s).second][j-1] < 1.0)
			{
				FuelToTestWithoutBuffer += StreamArray[(*it_i_s).second][j-1] *  (1-FractionMassToReachMin[(*it_i_s).second][j-1]);
				Buffer 				= BuildBuffer(FuelToTestWithoutBuffer , HMMass, StreamArray) ;
				FuelToTest 			= FuelToTestWithoutBuffer + Buffer ; 
				FuelToTest 			= FuelToTest/FuelToTest.GetSumOfAll();
				
				EqMMaximumBurnUp		= GetMaximumBurnUp (FuelToTest, BurnUp) ;

				if(EqMMaximumBurnUp=>BurnUp)
				{
					BurnUpExceeded 		= true;
					BurnUpExceededPosition	= i;
					HigherLimitOnBurnUp		= EqMMaximumBurnUp;
					BurnUpExceededList 		= (*it_i_s).second;
					break;
				}
			}

			for(int i=j; i < (int)StreamArray[(*it_i_s).second].size(); i++)
			{
				PreviousFuelToTestWithoutBuffer = FuelToTestWithoutBuffer ; //keep in memory fuel test during last step. When Burn-up is exceeded it will be the starting point.
				FuelToTestWithoutBuffer      +	= StreamArray[(*it_i_s).second][i] ;
				Buffer 				= BuildBuffer(FuelToTestWithoutBuffer , HMMass, StreamArray) ;
				FuelToTest 			= FuelToTestWithoutBuffer + Buffer ; 
				FuelToTest 			= FuelToTest/FuelToTest.GetSumOfAll();
				
				EqMMaximumBurnUp		= GetMaximumBurnUp (FuelToTest, BurnUp) ;

				if(EqMMaximumBurnUp=>BurnUp)
				{
					BurnUpExceeded 		= true;
					BurnUpExceededPosition	= i;
					HigherLimitOnBurnUp		= EqMMaximumBurnUp;
					BurnUpExceededList 		= (*it_i_s).second;
					break;
				}
			}
		}
		else
		{
			break;
		}
		StreamListNumber ++;
	}

	//No enough material to reach target burnup. Burnup never exceeded in loop.
	if (BurnUpExceeded==false)
	{
		WARNING << " Available material are not enough to reach the targeted BU.  Fuel not built." << endl;
		SetLambdaToErrorCode(lambda[(*it_s_vIV).first]);
		return lambda;	
	}

	//Lower limit leads to an higher burn-up than targeted burn-up 
	if (StreamListNumber == 0 && BurnUpExceededPosition==0)
	{
		WARNING << " Lower limit of first priority material is already to high for the target Burn-Up. Lower limit should be decreased.  Fuel not built." << endl;
		SetLambdaToErrorCode(lambda[(*it_s_vIV).first]);
		return lambda;	
	}

	double FractionOfLastIVToAdd 		= 1.0; //Start with 100% of the last IV in the fuel
	double LastFractionOfLastIVToAddMinus 	= 0.0; //Used in bissection method 
	double LastFractionOfLastIVToAddPlus	= 0.0; //Used in bissection method 	
	
	double CalculatedBurnUp 			= HigherLimitOnBurnUp ; 

	int count = 0;
	
	//Search in the fraction of last IV to add to the fuel to reach BU//
	do
	{
		if(count > fMaxIterration)
		{
			ERROR << "CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on burnup using :\nYourEquivalenceModel->SetBurnUpPrecision(Precision); " << endl;
			ERROR << "Targeted Burnup :" <<BurnUp<<endl;
			ERROR << "Last calculated Burnup :" <<CalculatedBurnUp<<endl;
			ERROR << "Last Fresh fuel composition without buffer:" <<endl;
			ERROR << FuelToTestWithoutBuffer .sPrint()<<endl;
			
			exit(1);
		}

		if( (CalculatedBurnUp - BurnUp) < 0 ) //Need to add more of last IV in fuel
		{
			LastFractionOfLastIVToAddMinus 	= FractionOfLastIVToAdd;
			FractionOfLastIVToAdd 		= FractionOfLastIVToAdd + fabs(LastFractionOfLastIVToAddPlus - FractionOfLastIVToAdd)/2.;
		}
		else if( (CalculatedBurnUp - BurnUp) > 0) //Need to add less of last IV in fuel
		{
			LastFractionOfLastIVToAddPlus 	= FractionOfLastIVToAdd;
			FractionOfLastIVToAdd 		= FractionOfLastIVToAdd - fabs(LastFractionOfLastIVToAddMinus - FractionOfLastIVToAdd)/2.;
		}
		FuelToTestWithoutBuffer  	= PreviousFuelToTestWithoutBuffer + FractionOfLastIVToAdd*StreamArray[BurnUpExceededList][i]; // Fuel tested at the step before Burn-up targeted is exceeded +  a fraction of last IV.
		IsotopicVector Buffer 		= BuildBuffer(FuelToTestWithoutBuffer , HMMass, StreamArray) ;
		IsotopicVector FuelToTest 	= FuelToTestWithoutBuffer + Buffer ; 
		CalculatedBurnUp 		= GetMaximumBurnUp(FuelToTest, BurnUp);
		count ++;

	}while(fabs(BurnUp - CalculatedBurnUp) > GetBurnUpPrecision()*BurnUp);

	for( it_s_IV = IVStream.begin();  it_s_IV != IVStream.end(); it_s_IV++)
		(*this).isIVInDomain(IVStream[(*it_s_IV).first]);
		
	for( it_s_vD = lambda.begin();  it_s_vD != lambda.end(); it_s_vD++)
	{	
		DBGV( "Lambda vector : "<<(*it_s_vD).first );

		for(int i=0; i < (int)lambda[(*it_s_vD).first].size(); i++)
		{
			DBGV(lambda[(*it_s_vD).first][i]); 
		}
	}
*/
	return lambda;
}
//________________________________________________________________________
IsotopicVector EquivalenceModel::BuildBuffer(IsotopicVector FuelToTestWithoutBuffer, double HMMass, map < string, vector < IsotopicVector > > StreamArray)
{

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

