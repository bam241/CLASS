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
	DBGL
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
	DBGL
}
//________________________________________________________________________
void EquivalenceModel::StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks)
{
	DBGL	
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
	DBGL
}

//________________________________________________________________________
void EquivalenceModel::ConvertMassToLambdaVector(string MaterialDenomination, vector<double>& lambda, double MaterialMassNeeded, vector <IsotopicVector>  Stocks)
{
	DBGL
	double Lambda_tot = 0; 

	// Calculation of Lambda tot associated to the required mass MaterialMassNeeded
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
	
	// Calculate lambda vector associated to the lambda tot 
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
	DBGL	
}

//________________________________________________________________________
IsotopicVector EquivalenceModel::BuildFuelToTest(map < string, vector<double> >& lambda, map < string , vector <IsotopicVector> > const& StreamArray, double HMMass, map <string, bool> StreamListIsBuffer)
{
	DBGL
	//Iterators declaration
	map < string , vector  <IsotopicVector> >::const_iterator it_s_vIV;
	map < string , bool >::iterator it_s_B;

	//Find the buffer and set its lambda to 0
	string BufferDenomination ="";
	for( it_s_B = StreamListIsBuffer.begin();  it_s_B != StreamListIsBuffer.end(); it_s_B++)
	{
		if((*it_s_B ).second==true){ BufferDenomination = (*it_s_B).first; }	
	}

	for(int i = 0; i< lambda[BufferDenomination].size(); i++)
	{
		lambda[BufferDenomination][i]=0;
	}
	
	//Build an IV with all materials besides buffer to get the total mass of others materials
	IsotopicVector IV;
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{	
		for(int i=0; i < (int)StreamArray.at( it_s_vIV->first ).size(); i++)
		{
			IV  +=  lambda[(*it_s_vIV).first][i] * StreamArray.at( it_s_vIV->first )[i];	
		}
	}

	//Calculate MassBuffer
	double MassBuffer = HMMass - IV.GetTotalMass()*1e06;

	//Set buffer lambda according to MassBuffer
	ConvertMassToLambdaVector(BufferDenomination, lambda[BufferDenomination], MassBuffer, StreamArray.at(BufferDenomination));

	IV.Clear();

	//Build fuel with all materials
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{	
		for(int i=0; i < (int)StreamArray.at( it_s_vIV->first ).size(); i++)
		{
			IV  +=  lambda[(*it_s_vIV).first][i] * StreamArray.at( it_s_vIV->first )[i];	
		}
	}
	DBGL
	return IV; 

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

	// Initialize lambda to 0 //
	for( it_s_vIV = StreamArray.begin();  it_s_vIV != StreamArray.end(); it_s_vIV++)
	{	
		for(int i=0; i < (int)StreamArray[(*it_s_vIV).first].size(); i++)
		{
			lambda[(*it_s_vIV).first].push_back(0);
		}
	}	
	// Check if the targeted burn-up is lower than maximum burn-up of model //
	if(BurnUp > fMaximalBU)
	{
		ERROR << " Targeted burn-up is higher than maximum burn-up handles by EqM."<< endl;
		ERROR << " Targeted burn-up : "<<BurnUp<<" GWd/t"<<endl;
		ERROR << " Maximum burn-up : "<<GetEqMHigherLimitOnBU()<<" GWd/t"<<endl;			
		exit(1);	
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

	StocksTotalMassCalculation(StreamArray);
	
	// Check if there is enough material in stock to satisfy mass fraction min //
	BreakReturnLambda = false; 
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

	// Check if there is enough material in stock to satisfy mass fraction max, if not mass fraction max is set to  MassINStock/MassReactor//
	for( it_s_D = StreamListMassFractionMax.begin();  it_s_D != StreamListMassFractionMax.end(); it_s_D++)
	{
		if(fTotalMassInStocks[(*it_s_D).first]< HMMass*StreamListMassFractionMax[(*it_s_D).first])
		{			
			StreamListMassFractionMax[(*it_s_D).first] = fTotalMassInStocks[(*it_s_D).first]/HMMass;
			WARNING << " Not enough material  : "<< (*it_s_D).first << " in stocks to reach the build fuel higher limit of "<<StreamListMassFractionMax[(*it_s_D).first]<<" reactor mass. " << endl;
			WARNING << " Mass fraction max of material :  "<< (*it_s_D).first << " is set to MassInStock/HMMassReactor : "<< StreamListMassFractionMax[(*it_s_D).first]<< endl;
		}
	}

	//Check targeted BU is inside [BUMin, BUMax] associated to fraction Min et Max//
	HMMass *=  1e6; //Unit conversion : tons to gram

	map < string , double > MassMin; 	
	map < string , double > MassMax;	 

	map < string , double > BUMin; 
	map < string , double > BUMax; 

	IsotopicVector FuelToTest;

	bool BurnUpMaxFound = false;

	for( it_i_s = StreamListPriority.begin();  it_i_s != StreamListPriority.end(); it_i_s++)
	{	
		//Calculate BU min for each possibility : min1 ; max1 + min2 ;  max1 + max2 + min3 ....
		MassMin[(*it_i_s ).second] 	=  HMMass * StreamListMassFractionMin[(*it_i_s).second];
		ConvertMassToLambdaVector((*it_i_s ).second, lambda[(*it_i_s ).second], MassMin[(*it_i_s ).second], StreamArray[(*it_i_s ).second]);
		FuelToTest 			= BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
		FuelToTest 			= FuelToTest/FuelToTest.GetSumOfAll();
		BUMin[(*it_i_s ).second] 	=  CalculateTargetParameter(FuelToTest);

		//Check is BUMin < BUTarget
		if(BUMin[(*it_i_s ).second]>BurnUp)
		{
			if((*it_i_s).first ==1) //Minimum of first material is too high
			{
				ERROR << "CRITICAL ! Minimum burn-up associated to the first priority material ( "<<(*it_i_s ).second <<" ) is higher than targeted burn-up."<< endl;
				ERROR << "Targeted Burn-up : " <<BurnUp<<endl;
				ERROR << "Minimum Burn-up : " <<BUMin[(*it_i_s ).second]<<endl;
				ERROR << "Try to increase targeted burn-up" <<endl;
				exit(1);
			}
			else if ((*it_i_s).first >1) //BU target is located between max n-1 and min n
			{
				ERROR << "CRITICAL ! Targeted burn-up is located between 2 materials. "<<endl;
				it_i_s --;
				ERROR << "Maximum Burn-up of : "<< (*it_i_s).second<<" ---> "<<BUMax[(*it_i_s ).second]<<endl;
				it_i_s ++;
				ERROR << "Minimum Burn-up of : "<< (*it_i_s ).second<<" ---> "<<BUMin[(*it_i_s ).second]<<endl;
				ERROR << "Targeted Burn-up : " <<BurnUp<<endl;					
				ERROR << "Try to decrease mimim fraction of : "<< (*it_i_s ).second<<endl;
				exit(1);
			}
		}
		FuelToTest.Clear();

		//Calculate BU max for each possibility : max1 ; max1 + max2 ;  max1 + max2 + max3 ....
		MassMax[(*it_i_s ).second] 	=  HMMass * StreamListMassFractionMax[(*it_i_s).second];	
		ConvertMassToLambdaVector((*it_i_s ).second, lambda[(*it_i_s ).second], MassMax[(*it_i_s ).second], StreamArray[(*it_i_s ).second]);
		FuelToTest 			= BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
		FuelToTest 			= FuelToTest/FuelToTest.GetSumOfAll();
		BUMax[(*it_i_s ).second] 	=  CalculateTargetParameter(FuelToTest);
		
		if(BUMax[(*it_i_s ).second]>BurnUp)
		{
			BurnUpMaxFound = true ; 
			break;
		}
	}

	if(!BurnUpMaxFound) 
	{
		ERROR << "CRITICAL ! Maximum reachable burn-up is lower than target BU. "<< endl;
		ERROR << "Targeted Burn-up : " <<BurnUp<<endl;
		ERROR << "Maximum reachable burn-up : " <<BUMax[(*--StreamListPriority.end()).second]<<endl;
		ERROR << "Try to increase maximum fraction of materials, or decrease burn-up. " <<endl;
		exit(1);
	}

	//Search the BU max //	
	string MaterialToSearch 	= (*it_i_s ).second;
	double CalculatedBurnUp 	= BUMax[MaterialToSearch] ;   //Algo start with maximum point
	double MassToAdd 		= MassMax[MaterialToSearch]; //Algo start with maximum point
	
	double LastMassMinus 	= MassMin[MaterialToSearch]; //Used in bissection method 
	double LastMassPlus		= MassMax[MaterialToSearch]; //Used in bissection method 	

	int count = 0;
	
	FuelToTest.Clear();

	do
	{
		if(count > fMaxIterration)
		{
			ERROR << "CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on burnup using :\nYourEquivalenceModel->SetBurnUpPrecision(Precision); " << endl;
			ERROR << "Targeted Burnup : "  <<BurnUp<<endl;
			ERROR << "Last calculated Burnup : " <<CalculatedBurnUp<<endl;
			ERROR << "Last Fresh fuel normalized composition : " <<endl;
			ERROR << FuelToTest.sPrint()<<endl;	
			exit(1);
		}

		if( (CalculatedBurnUp - BurnUp) < 0 ) //Need to add more fissile material in fuel
		{
			LastMassMinus = MassToAdd;
			MassToAdd 	= MassToAdd + fabs(LastMassPlus - MassToAdd)/2.;
		}
		else if( (CalculatedBurnUp - BurnUp) > 0) //Need to add less fissile material in fuel
		{
			LastMassPlus 	= MassToAdd;
			MassToAdd 	= MassToAdd - fabs(LastMassMinus - MassToAdd)/2.;
		}
		ConvertMassToLambdaVector(MaterialToSearch, lambda[MaterialToSearch], MassToAdd, StreamArray[MaterialToSearch]);
		
		FuelToTest 		= BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
		FuelToTest 		= FuelToTest/FuelToTest.GetSumOfAll();
		CalculatedBurnUp 	= CalculateTargetParameter(FuelToTest);
		
		count ++;

	}while(fabs(BurnUp - CalculatedBurnUp) > GetBurnUpPrecision()*BurnUp);
	
	//Final builded fuel 
	IsotopicVector IVStream;
	for( it_s_vD = lambda.begin();  it_s_vD != lambda.end(); it_s_vD++)
	{
		for(int i=0; i<(int)lambda[(*it_s_vD).first].size(); i++) 
		{
			IVStream +=lambda[(*it_s_vD).first][i] * StreamArray[(*it_s_vD).first][i];
		}
	}	
	
	//Check if BuildedFuel is in Model isotopic bounds 
	(*this).isIVInDomain(IVStream);
		
	for( it_s_vD = lambda.begin();  it_s_vD != lambda.end(); it_s_vD++)
	{	
		DBGV( "Lambda vector : "<<(*it_s_vD).first );

		for(int i=0; i < (int)lambda[(*it_s_vD).first].size(); i++)
		{
			DBGV(lambda[(*it_s_vD).first][i]); 
		}
	}
	DBGL
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

