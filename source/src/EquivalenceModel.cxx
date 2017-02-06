#include "EquivalenceModel.hxx"
#include "external/StringLine.hxx"
#include "CLASSMethod.hxx"

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "CLASSReader.hxx"

//________________________________________________________________________
//________________________________________________________________________
EquivalenceModel::EquivalenceModel(string TMVAXMLFilePath, string TMVANFOFilePath):CLASSObject()
{
	fMaxIterration 	= 500;
	freaded 	= false;		
	fPCMprecision  = 10;

    fTMVAXMLFilePath = TMVAXMLFilePath;
    fTMVANFOFilePath = TMVANFOFilePath;

    LoadKeyword();  // Load Key words defineds in NFO file
    ReadNFO();      //Getting information from file NFO

    //Check if k_targetparameter is defined in file NFO
    if(fTargetParameter.empty())  
    {
        ERROR<<"Missing information for k_targetparameter in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_modelparameter is defined in file NFO
    if(fModelParameter.empty())  
    {
        ERROR<<"Missing information for k_modelparameter in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_predictortype is defined in file NFO
    if(fPredictorType.empty())  
    {
        ERROR<<"Missing information for k_predictortype in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_output is defined in file NFO
    if(fOutput.empty())  
    {
        ERROR<<"Missing information for k_output in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_buffer is defined in file NFO
    if(fBuffer.empty())  
    {
        ERROR<<"Missing information for k_buffer in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_massfractionmin and k_massfractionmax are defined in file NFO
    if(fStreamListEqMMassFractionMin.empty() || fStreamListEqMMassFractionMax.empty())  
    {
        ERROR<<"Missing information for k_massfractionmin and/or k_massfractionmax in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_zainame are defined in file NFO
    if(fMapOfTMVAVariableNames.empty())  
    {
        ERROR<<"Missing information for k_zainame in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_list are defined in file NFO
    if(fStreamList.empty())  
    {
        ERROR<<"Missing information for k_list in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }

    INFO << "__An equivalence model has been define__" << endl;
    INFO << "\tThe TMVA weights file is :" << fTMVAXMLFilePath << endl;
    INFO << "\tThe TMVA NFO file is :" << fTMVANFOFilePath << endl;
    PrintInfo();

}
//________________________________________________________________________
EquivalenceModel::EquivalenceModel(CLASSLogger* log, string TMVAXMLFilePath, string TMVANFOFilePath):CLASSObject(log)
{
	fMaxIterration 	= 500;
	freaded 	= false;		
	fPCMprecision  = 10;

    fTMVAXMLFilePath = TMVAXMLFilePath;
    fTMVANFOFilePath = TMVANFOFilePath;

    LoadKeyword();  // Load Key words defineds in NFO file
    ReadNFO();      //Getting information from file NFO

    //Check if k_targetparameter is defined in file NFO
    if(fTargetParameter.empty())  
    {
        ERROR<<"Missing information for k_targetparameter in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_modelparameter is defined in file NFO
    if(fModelParameter.empty())  
    {
        ERROR<<"Missing information for k_modelparameter in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_predictortype is defined in file NFO
    if(fPredictorType.empty())  
    {
        ERROR<<"Missing information for k_predictortype in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_output is defined in file NFO
    if(fOutput.empty())  
    {
        ERROR<<"Missing information for k_output in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_buffer is defined in file NFO
    if(fBuffer.empty())  
    {
        ERROR<<"Missing information for k_buffer in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_massfractionmin and k_massfractionmax are defined in file NFO
    if(fStreamListEqMMassFractionMin.empty() || fStreamListEqMMassFractionMax.empty())  
    {
        ERROR<<"Missing information for k_massfractionmin and/or k_massfractionmax in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_zainame are defined in file NFO
    if(fMapOfTMVAVariableNames.empty())  
    {
        ERROR<<"Missing information for k_zainame in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }
    //Check if k_list are defined in file NFO
    if(fStreamList.empty())  
    {
        ERROR<<"Missing information for k_list in : "<<fTMVANFOFilePath<<endl;
        exit(1);
    }

    INFO << "__An equivalence model has been define__" << endl;
    INFO << "\tThe TMVA weights file is :" << fTMVAXMLFilePath << endl;
    INFO << "\tThe TMVA NFO file is :" << fTMVANFOFilePath << endl;
    PrintInfo();}
//________________________________________________________________________
EquivalenceModel::~EquivalenceModel()
{
							
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

    // Check if EquivalenceModel->SetTMVAXMLFilePath() and/or EquivalenceModel->SetTMVANFOFilePath() have been defined
    if (fTMVAXMLFilePath.empty() || fTMVANFOFilePath.empty())
    {
        ERROR << " TMVA XML and/or NFO File path are not defined..."<< endl;
        ERROR << " You have to use EquivalenceModel->SetTMVAXMLFilePath() and/or EquivalenceModel->SetTMVANFOFilePath() methods."<<endl;
        exit(1);
    }

	//Iterators declaration
	map < string , vector  <IsotopicVector> >::iterator it_s_vIV;
	map < string , vector  <double> >::iterator it_s_vD;
	map < string , IsotopicVector >::iterator it_s_IV;
	map < string , double >::iterator it_s_D;
	map < int , string >::iterator it_i_s;

	double TargetParameterValue = 0;

	for( it_s_D = fModelParameter.begin();  it_s_D != fModelParameter.end(); it_s_D++)
	{	
		if(fModelParameter[(*it_s_D).first] == -1)
		{
			ERROR<< "Model parameter ( "<<fModelParameter[(*it_s_D).first] << " ) value is not defined in the input." <<endl;
			ERROR<< "Use EqM->SetModelParameter( \" "<<(*it_s_D).first<<" \", value) to define it." <<endl;
			exit(1);			
		}
	}

	if(fTargetParameter=="BurnUpMax") {TargetParameterValue  = BurnUp;}
	else if (fTargetParameter=="keffBOC") {TargetParameterValue = fModelParameter["keffBOC"];}
	else
	{
		ERROR<< "Target parameter defined in InformationFile ( "<<fTargetParameter<<" ) doesn't exist." <<endl;
		ERROR<< "Possible target parameters for the moment are : BurnUpMax and keffBOC." <<endl;
		exit(1);			
	}

	map <string , vector<double> > lambda ; // map containing name of the list and associated vector of proportions taken from stocks
	
	// Initialize lambda to 0 //
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
	
	// Check if the targeted burn-up is lower than maximum burn-up of model //
	if(BurnUp > fMaximalBU)
	{
		ERROR << " Targeted burn-up is higher than maximum burn-up handles by EqM."<< endl;
		ERROR << " Targeted burn-up : "<<BurnUp<<" GWd/t"<<endl;
		ERROR << " Maximum burn-up : "<<fMaximalBU<<" GWd/t"<<endl;			
		exit(1);	
	}

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

	//Calculate Total mass in stock for each stream and fill fTotalMassInStocks
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

	// Check if there is enough material in stock to satisfy mass fraction max, if not mass fraction max is set to MassINStock/MassReactor//
	for( it_s_D = StreamListMassFractionMax.begin();  it_s_D != StreamListMassFractionMax.end(); it_s_D++)
	{
		if(fTotalMassInStocks[(*it_s_D).first]< HMMass*StreamListMassFractionMax[(*it_s_D).first])
		{			
			StreamListMassFractionMax[(*it_s_D).first] = fTotalMassInStocks[(*it_s_D).first]/HMMass;
			WARNING << " Not enough material  : "<< (*it_s_D).first << " in stocks to reach the build fuel higher limit of "<<StreamListMassFractionMax[(*it_s_D).first]<<" reactor mass. " << endl;
			WARNING << " Mass fraction max of material :  "<< (*it_s_D).first << " is set to MassInStock/HMMassReactor : "<< StreamListMassFractionMax[(*it_s_D).first]<< endl;
		}
	}
	
	//Check if TargetParameter is inside [TargetParameterMin, TargetParameterMax] associated to fraction Min et Max//
	HMMass *=  1e6; //Unit conversion : tons to gram

	map < string , double > MassMin; 	
	map < string , double > MassMax;	 

	map < string , double > TargetParameterMin; 
	map < string , double > TargetParameterMax; 

	IsotopicVector FuelToTest;

	bool TargetParameterIncluded = false;

	for( it_i_s = StreamListPriority.begin();  it_i_s != StreamListPriority.end(); it_i_s++)
	{	
		//Calculate TargetParameterMin for each possibility : min1 ; max1 + min2 ;  max1 + max2 + min3 ....
		MassMin[(*it_i_s ).second] 		=  HMMass * StreamListMassFractionMin[(*it_i_s).second];
		ConvertMassToLambdaVector((*it_i_s ).second, lambda[(*it_i_s ).second], MassMin[(*it_i_s ).second], StreamArray[(*it_i_s ).second]);
		FuelToTest 				= BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
		FuelToTest 				= FuelToTest/FuelToTest.GetSumOfAll();
		TargetParameterMin[(*it_i_s ).second] =  CalculateTargetParameter(FuelToTest, fTargetParameter);

		//Check is TargetParameterMin < TargetParameter
		if(TargetParameterMin[(*it_i_s ).second]>TargetParameterValue)
		{
			if((*it_i_s).first ==1) //Minimum of first material is too high
			{
				ERROR << "CRITICAL ! Minimum parameter value associated to the first priority material ( "<<(*it_i_s ).second <<" ) is higher than targeted parameter."<< endl;
				ERROR << "Targeted parameter : "<<fTargetParameter<<" = "<<TargetParameterValue<<endl;
				ERROR << "Minimum parameter value : " <<TargetParameterMin[(*it_i_s ).second]<<endl;
				ERROR << "Try to increase targeted burn-up" <<endl;
				exit(1);
			}
			else if ((*it_i_s).first >1) //TargetParameter is located between max n-1 and min n
			{
				ERROR << "CRITICAL ! Targeted parameter value ( "<<fTargetParameter<<" ) is located between 2 materials. "<<endl;
				it_i_s --;
				ERROR << fTargetParameter <<" of max fraction of material : "<< (*it_i_s).second<<" ---> "<<TargetParameterMax[(*it_i_s ).second]<<endl;
				it_i_s ++;
				ERROR << fTargetParameter<<  " of min fraction of material : "<< (*it_i_s ).second<<" ---> "<<TargetParameterMin[(*it_i_s ).second]<<endl;
				ERROR << "Targeted "<<fTargetParameter<<" : " <<TargetParameterValue<<endl;					
				ERROR << "Try to decrease mimimum fraction of : "<< (*it_i_s ).second<<endl;
				exit(1);
			}
		}
		FuelToTest.Clear();

		//Calculate TargetParameter max for each possibility : max1 ; max1 + max2 ;  max1 + max2 + max3 ....
		MassMax[(*it_i_s ).second] 	=  HMMass * StreamListMassFractionMax[(*it_i_s).second];	
		ConvertMassToLambdaVector((*it_i_s ).second, lambda[(*it_i_s ).second], MassMax[(*it_i_s ).second], StreamArray[(*it_i_s ).second]);
		FuelToTest 			= BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
		FuelToTest 			= FuelToTest/FuelToTest.GetSumOfAll();
		TargetParameterMax[(*it_i_s ).second] 	=  CalculateTargetParameter(FuelToTest, fTargetParameter);
		
		if(TargetParameterMax[(*it_i_s ).second]>BurnUp)
		{
			TargetParameterIncluded = true ; 
			break;
		}
	}

	//Check if target parameter increases monotously with the material mass
	CheckTargetParameterConsistency(StreamListPriority, TargetParameterMin, TargetParameterMax);

	if(!TargetParameterIncluded) 
	{
		ERROR << "CRITICAL ! Maximum reachable "<<fTargetParameter<<" is lower than targeted "<< fTargetParameter<<". "<< endl;
		ERROR << "Targeted "<<fTargetParameter<<" = "<<TargetParameterValue<<endl;
		ERROR << "Maximum reachable "<<fTargetParameter<<" : "<<TargetParameterMax[(*--StreamListPriority.end()).second]<<endl;
		ERROR << "Try to increase maximum fraction of materials, or decrease "<< fTargetParameter<<" ." <<endl;
		exit(1);
	}

	//Search the TargetParameterValue location in the mass damain //	
	string MaterialToSearch 		= (*it_i_s ).second;
	double CalculatedTargetParameter 	= TargetParameterMax[MaterialToSearch] ;   //Algo start with maximum point
	double MassToAdd 			= MassMax[MaterialToSearch]; //Algo start with maximum point
	
	double LastMassMinus 		= MassMin[MaterialToSearch]; //Used in bissection method 
	double LastMassPlus			= MassMax[MaterialToSearch]; //Used in bissection method 	

	int count = 0;
	
	FuelToTest.Clear();

	do
	{
                      if(count > fMaxIterration-10)
                      {
                          cout << "CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on the target parameter using :\nYourEquivalenceModel->SetTargetParameterStDev(Precision); " << endl;
                          cout << "Targeted "<<fTargetParameter<<" : "<<TargetParameterValue<<endl;
                          cout << "Last calculated "<<fTargetParameter<<" : "<<CalculatedTargetParameter<<endl;
                          cout << "Last Fresh fuel normalized composition : " <<endl;
                          cout << FuelToTest.sPrint()<<endl; 
                      }

		if(count > fMaxIterration)
		{
			ERROR << "CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on the target parameter using :\nYourEquivalenceModel->SetTargetParameterStDev(Precision); " << endl;
			ERROR << "Targeted "<<fTargetParameter<<" : "<<TargetParameterValue<<endl;
			ERROR << "Last calculated "<<fTargetParameter<<" : "<<CalculatedTargetParameter<<endl;
			ERROR << "Last Fresh fuel normalized composition : " <<endl;
			ERROR << FuelToTest.sPrint()<<endl;	
			exit(1);
		}

		if( (CalculatedTargetParameter - TargetParameterValue) < 0 ) //Need to add more fissile material in fuel
		{
			LastMassMinus = MassToAdd;
			MassToAdd 	= MassToAdd + fabs(LastMassPlus - MassToAdd)/2.;
		}
		else if( (CalculatedTargetParameter - TargetParameterValue) > 0) //Need to add less fissile material in fuel
		{
			LastMassPlus 	= MassToAdd;
			MassToAdd 	= MassToAdd - fabs(LastMassMinus - MassToAdd)/2.;
		}
		ConvertMassToLambdaVector(MaterialToSearch, lambda[MaterialToSearch], MassToAdd, StreamArray[MaterialToSearch]);
		
		FuelToTest 			= BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
		FuelToTest 			= FuelToTest/FuelToTest.GetSumOfAll();
		CalculatedTargetParameter 	= CalculateTargetParameter(FuelToTest, fTargetParameter);
		
		count ++;

	}while(fabs(TargetParameterValue - CalculatedTargetParameter) > GetTargetParameterStDev()*TargetParameterValue);
	
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
TTree* EquivalenceModel::CreateTMVAInputTree(IsotopicVector TheFreshfuel, double ThisTime)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/
	TTree*   InputTree = new TTree(fOutput.c_str(), fOutput.c_str());
	
	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);
	
	float Time = 0;
	
	IsotopicVector IVInputTMVA;
	map<ZAI ,string >::iterator it_ZAI_s;
	int j = 0;
	
	for( it_ZAI_s = fMapOfTMVAVariableNames.begin()  ; it_ZAI_s != fMapOfTMVAVariableNames.end() ; it_ZAI_s++)
	{
		InputTree->Branch( ((*it_ZAI_s).second).c_str(), &InputTMVA[j], ((*it_ZAI_s).second + "/F").c_str());
		IVInputTMVA+=  ((*it_ZAI_s).first)*1;
		j++;
	}
	
	if(ThisTime != -1)
		InputTree->Branch("Time" ,&Time ,"Time/F");
	
	IsotopicVector IVAccordingToUserInfoFile 	= TheFreshfuel.GetThisComposition(IVInputTMVA);
	double Ntot 					= IVAccordingToUserInfoFile.GetSumOfAll();
	IVAccordingToUserInfoFile 			= IVAccordingToUserInfoFile/Ntot;
	
	j = 0;
	
	for( it_ZAI_s = fMapOfTMVAVariableNames.begin() ; it_ZAI_s != fMapOfTMVAVariableNames.end() ; it_ZAI_s++)
	{
		InputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it_ZAI_s).first ) ;
		j++;
	}
	
	Time = ThisTime;
	InputTree->Fill();
	
	return InputTree;
	
}
//________________________________________________________________________
void EquivalenceModel::CheckTargetParameterConsistency(map < int , string >  StreamListPriority, map < string , double >  TargetParameterMin, map < string , double > TargetParameterMax)
{
	map < int , string >::iterator it_i_s;
	
	double TargetParameterUp 		= 0;
	double TargetParameterDown 		= 0;

	//Loop on priority order to check if target parameter increases monotously with the material mass
	for( it_i_s = StreamListPriority.begin();  it_i_s != StreamListPriority.end(); it_i_s++)
	{	
		TargetParameterDown = TargetParameterMin[(*it_i_s ).second];
		if (TargetParameterDown < TargetParameterUp )
		{
			ERROR<< "Target parameter evolution as a function of material mass is not monotonous." <<endl;
			ERROR<< "Check the its evolution." <<endl;
			exit(1);			

		}
		
		TargetParameterUp 	= TargetParameterMax[(*it_i_s ).second];
		if (TargetParameterDown > TargetParameterUp )
		{			
			ERROR<< "Target parameter evolution as a function of material mass is not monotonous." <<endl;
			ERROR<< "Check the its evolution." <<endl;
			exit(1);			
		}
	}
}
//________________________________________________________________________
double EquivalenceModel::CalculateTargetParameter(IsotopicVector TheFuel, string TargetParameterName)
{
	double ParameterToCalculate = 0; 
	if(TargetParameterName=="BurnUpMax") ParameterToCalculate 	= CalculateBurnUpMax(TheFuel, fModelParameter);
	else if(TargetParameterName=="keffBOC") ParameterToCalculate 	= CalculateKeffAtBOC(TheFuel);
	else
	{
		ERROR<< "Target parameter defined in InformationFile ( "<<TargetParameterName<<" ) doesn't exist" <<endl;
		ERROR<< "Possible target parameters for the moment are : BurnUpMax and keffBOC." <<endl;
		exit(1);			
	}

	return ParameterToCalculate ;
}
//________________________________________________________________________
double EquivalenceModel::CalculateBurnUpMax(IsotopicVector TheFuel, map<string, double> ModelParameter)
{
	/**************************************************************************/
	//With a dichotomy, the maximal irradiation time (TheFinalTime) is calculated
	//When average Kinf is very close (according "Precision") to the threshold
	//then the corresponding irradiation time is convert in burnup and returned
	/**************************************************************************/
	//Algorithm initialization
	double KThreshold  		= fModelParameter["kThreshold"];
	int NumberOfBatch	 	= (int)fModelParameter["NumberOfBatch"];
	double OldFinalTimeMinus 	= 0;
	double MaximumBU 		= fMaximalBU;
	double MinimumBU  		= 0 ;
	double TheFinalTime 		= BurnupToSecond((MaximumBU-MinimumBU)/2.);
	double OldFinalTimePlus 	= BurnupToSecond(MaximumBU);
	double k_av 			= 0; //average kinf
	double OldPredictedk_av 	= 0;
	
	CLASSReader * reader = new CLASSReader( fMapOfTMVAVariableNames );
	reader->AddVariable( "Time" );
	reader->BookMVA( "MLP method" , fTMVAXMLFilePath );
	
	for(int b = 0;b<NumberOfBatch;b++)
	{
		float TheTime = (b+1)*TheFinalTime/NumberOfBatch;

		TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
		reader->SetInputData( InputTree );
		
		OldPredictedk_av += reader->EvaluateRegression( "MLP method" )[0];
		
		delete InputTree;
	}
	OldPredictedk_av /= NumberOfBatch;
	
	//Algorithm control
	int count = 0;
	int MaximumLoopCount = 500;
	do
	{
		if(count > MaximumLoopCount )
		{
			ERROR << "CRITICAL ! Can't manage to predict burnup\nHint : Try to increase the precision on k effective using :\n YourEQM_MLP_Kinf->SetPCMPrecision(pcm); with pcm the precision in pcm (default 10) REDUCE IT\n If this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr " << endl;
			exit(1);
		}
		
		if( (OldPredictedk_av-KThreshold)  > 0) //The burnup can be increased
		{
			OldFinalTimeMinus = TheFinalTime;
			TheFinalTime = TheFinalTime + fabs(OldFinalTimePlus - TheFinalTime)/2.;
			
			if(SecondToBurnup(TheFinalTime) >= (MaximumBU-MaximumBU*GetTargetParameterStDev() ) )
				{ delete reader; return MaximumBU; }
		}
		
		else if( (OldPredictedk_av-KThreshold)  < 0)//The burnup is too high
		{
			OldFinalTimePlus = TheFinalTime;
			TheFinalTime = TheFinalTime - fabs(OldFinalTimeMinus-TheFinalTime)/2.;
			if( SecondToBurnup(TheFinalTime) < (MaximumBU-MinimumBU)/2.*GetTargetParameterStDev() )
				{ delete reader; return 0; }
		}
		
		k_av = 0;
		for(int b = 0;b<NumberOfBatch;b++)
		{
			float TheTime = (b+1)*TheFinalTime/NumberOfBatch;
			TTree* InputTree = CreateTMVAInputTree(TheFuel,TheTime);
			reader->SetInputData( InputTree );
			
			k_av += reader->EvaluateRegression("MLP method")[0];
			
			delete InputTree;
		}
		k_av/= NumberOfBatch;
		//cout<<SecondToBurnup(TheFinalTime)<<" ";
		OldPredictedk_av = k_av;
		count++;
//std::clog << "-> " << k_av << "\t\t(" << count << ") \t [" << TheFinalTime << "]" << "\t" << OldPredictedk_av-KThreshold << "\t" << GetPCMPrecision() << std::endl; 
	}	while( fabs(OldPredictedk_av-KThreshold) > GetPCMPrecision() )  ;
	
	delete reader;
	//cout<<endl;
	return SecondToBurnup(TheFinalTime);
}

//________________________________________________________________________
double 	EquivalenceModel::CalculateKeffAtBOC(IsotopicVector FreshFuel)
{ 
	CLASSReader * reader = new CLASSReader( fMapOfTMVAVariableNames );
	reader->BookMVA( "MLP method" , fTMVAXMLFilePath );

	TTree* InputTree = CreateTMVAInputTree(FreshFuel,-1) ; 
	reader->SetInputData( InputTree );

	double keff =  reader->EvaluateRegression( "MLP method" )[0];

	delete InputTree;

	return keff;
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



//________________________________________________________________________
void EquivalenceModel::ReadNFO()
{
	DBGL
	ifstream NFO(fTMVANFOFilePath.c_str());
	
	if(!NFO)
	{
		ERROR << "Can't find/open file " << fTMVANFOFilePath << endl;
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
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_zail",				& EquivalenceModel::ReadZAIlimits)	 	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_reactor",			& EquivalenceModel::ReadType)	 	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_fuel",				& EquivalenceModel::ReadType)	 	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_massfractionmin",		& EquivalenceModel::ReadEqMinFraction) 	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_massfractionmax",		& EquivalenceModel::ReadEqMaxFraction) 	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_list",				& EquivalenceModel::ReadList) 	 	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_specpower",			& EquivalenceModel::ReadSpecificPower)	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_zainame", 			& EquivalenceModel::ReadZAIName) 		 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_maxburnup", 			& EquivalenceModel::ReadMaxBurnUp) 	 	 );	
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_targetparameter",		& EquivalenceModel::ReadTargetParameter) 	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_predictortype",		& EquivalenceModel::ReadPredictorType)	 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_output", 			& EquivalenceModel::ReadOutput) 		 	 );
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_buffer", 			& EquivalenceModel::ReadBuffer)	 	 	 );	
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_modelparameter", 		& EquivalenceModel::ReadModelParameter) 	 	 );	
	fKeyword.insert( pair<string, EQM_MthPtr>( "k_targetparameterstdev", 	& EquivalenceModel::ReadTargetParameterStDev) 	 );	

	DBGL
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
	double upLimit 	= atof(StringLine::NextWord(line, pos, ' ').c_str());
	
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
void EquivalenceModel::ReadZAIName(const string &line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_zainame" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_zainame\" not found !" << endl;
		exit(1);
	}
	
	int Z = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	int I  = atoi(StringLine::NextWord(line, pos, ' ').c_str());
	
	string name = StringLine::NextWord(line, pos, ' ');
	
	fMapOfTMVAVariableNames.insert( pair<ZAI,string>( ZAI(Z, A, I), name ) );

	DBGL	
}
//________________________________________________________________________
void EquivalenceModel::ReadMaxBurnUp(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_maxburnup" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_maxburnup\" not found !" << endl;
		exit(1);
	}
	
	fMaximalBU = atof(StringLine::NextWord(line, pos, ' ').c_str());

	DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadTargetParameter(const string &line)
{
	DBGL
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_targetparameter" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_targetparameter\" not found !" << endl;
		exit(1);
	}
	
	fTargetParameter = StringLine::NextWord(line, pos, ' ');

	DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadPredictorType(const string &line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_predictortype" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_predictortype\" not found !" << endl;
		exit(1);
	}
		
	fPredictorType = StringLine::NextWord(line, pos, ' ');
	
	DBGL	
}
//________________________________________________________________________
void EquivalenceModel::ReadOutput(const string &line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_output" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_output\" not found !" << endl;
		exit(1);
	}
		
	fOutput = StringLine::NextWord(line, pos, ' ');
	
	DBGL	
}
//________________________________________________________________________
void EquivalenceModel::ReadBuffer(const string &line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_buffer" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_buffer\" not found !" << endl;
		exit(1);
	}
		
	fBuffer = StringLine::NextWord(line, pos, ' ');
	
	DBGL	
}
//________________________________________________________________________
void EquivalenceModel::ReadModelParameter(const string &line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_modelparameter" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_modelparameter\" not found !" << endl;
		exit(1);
	}
		
	keyword = StringLine::NextWord(line, pos, ' ');

	fModelParameter[keyword] = -1;
	
	DBGL	
}

//________________________________________________________________________
void EquivalenceModel::ReadTargetParameterStDev(const string &line)
{
	DBGL
	
	int pos = 0;
	string keyword = tlc(StringLine::NextWord(line, pos, ' '));
	if( keyword != "k_targetparameterstdev" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_targetparameterstdev\" not found !" << endl;
		exit(1);
	}
		
	fTargetParameterStDev = atof(StringLine::NextWord(line, pos, ' ').c_str());;
	
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
