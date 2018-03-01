/*****************************************************************/
//              DESCRIPTION
// Illustrative scenario to show how to use the
// SeparationPlant and BackEnd Facilities (Pool/Storage/Separation)
// This park is constituted by a simple PWR UOX Reactor, a pool
// 4 Storage and 2 FabricationPlant.
//   _______     _______________  
//  |       |   |    	     	 | =>       Pool      => Storage1
//  |Reactor| =>|SeparationPlant1| => Storage2  
//  |_______|   |_______________ | =>|SeparationPlant2|=>Storage3
//                                   |                |=>Storage4
//
//
//@author BaL
/****************************************************************/
#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
	//seconds in one year
	cSecond year = 3600*24.*365.25; 
	/******LOG MANAGEMENT**********************************/
	//Definition of the Log file : CLASS messages output 
	int Std_output_level 	= 0; //Only error are shown in terminal
	int File_output_level 	= 2; // Error + Warning + Info are shown in the file CLASS_OUTPUT.log
	CLASSLogger *Logger 	= new CLASSLogger("CLASS_OUTPUT.log",Std_output_level,File_output_level);

	/******SCENARIO**********************************/
	// The scenario start at year 1977
	Scenario *gCLASS=new Scenario(1977*year,Logger);
	gCLASS->SetStockManagement(true);					//If false all the IsotopicVector in stocks are mixed together.
	gCLASS->SetTimeStep(year/4.);						//the scenario calculation is updated every 3 months
	gCLASS->SetOutputFileName("Separation.root");		//Set the name of the output file

	/******DATA BASES**********************************/
	//Geting CLASS to path
	string CLASS_PATH = getenv("CLASS_PATH");
	if (CLASS_PATH=="")
   	{
		cout<<" Please setenv CLASS_PATH to your CLASS installation folder in your .bashs or .tcshrc"<<endl;
   	 	exit(0);
   	}
   	string PATH_TO_DATA = CLASS_PATH + "/DATA_BASES/";

	/*===Decay data base===*/
	//The decay data base is taken from the file Decay.idx
	DecayDataBank* DecayDB = new DecayDataBank(gCLASS->GetLog(), PATH_TO_DATA + "DECAY/ALL/Decay.idx");
	gCLASS->SetDecayDataBase(DecayDB);//This decay data base will be used for all the decay calculations in this Scenario

	/*===Reactor data base===*/
	//The file STD900.dat correspond to a fuel evolution of a UOX PWR (see manual for details)
	EvolutionData *STD900 = new EvolutionData(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/FixedFuel/STD900.dat");
	/******FACILITIES*********************************/
	/*===The 4 Stocks ===*/
	Storage* MyStorage1 = new Storage(Logger);
	MyStorage1->SetName("Storage1");
	gCLASS->Add(MyStorage1); //Adding the stock to the Scenario 

	Storage* MyStorage2 = new Storage(Logger);
	MyStorage2->SetName("Storage2");
	gCLASS->Add(MyStorage2); //Adding the stock to the Scenario 

	Storage* MyStorage3 = new Storage(Logger);
	MyStorage3->SetName("Storage3");
	gCLASS->Add(MyStorage3); //Adding the stock to the Scenario 

	Storage* MyStorage4 = new Storage(Logger);	
	MyStorage4->SetName("Storage4");
	gCLASS->Add(MyStorage4); //Adding the stock to the Scenario 


	/*===A Pool===*/
	Pool* MyPool1 = new Pool ( Logger , MyStorage1 , 5*year) ; //Defined the Pool MyPool1 connected to the Storage MyStorage1. 5year of cooling
															   //After 5 years of cooling, the pool sends its content to "MyStorage1"
	MyPool1->SetName("MyPool1");
	gCLASS->Add(MyPool1);


	//*=======2 Separation Plants===*/
	/***Separation Efficencies Definition***/
	IsotopicVector UraniumSeparationEfficiency;
	UraniumSeparationEfficiency.Add(92,232,0,0.99);//Efficiency of 99%
	UraniumSeparationEfficiency.Add(92,233,0,0.99);
	UraniumSeparationEfficiency.Add(92,234,0,0.99);
	UraniumSeparationEfficiency.Add(92,235,0,0.99);
	UraniumSeparationEfficiency.Add(92,236,0,0.99);
	UraniumSeparationEfficiency.Add(92,238,0,0.99);

	IsotopicVector PlutoniumSeparationEfficiency;
	PlutoniumSeparationEfficiency.Add(94,236,0,0.98);//Efficiency of 98%
	PlutoniumSeparationEfficiency.Add(94,238,0,0.98);//Efficiency of 98%
	PlutoniumSeparationEfficiency.Add(94,239,0,0.98);
	PlutoniumSeparationEfficiency.Add(94,240,0,0.98);
	PlutoniumSeparationEfficiency.Add(94,241,0,0.98);
	PlutoniumSeparationEfficiency.Add(94,242,0,0.98);

	IsotopicVector AmericiumSeparationEfficiency;
	AmericiumSeparationEfficiency.Add(95,241,0,0.98);//Efficiency of 98%
	AmericiumSeparationEfficiency.Add(95,243,0,0.98);//Efficiency of 98%

	//Yes we assume that may be one day, such isotopic separation will be doable !!
	IsotopicVector Am241SeparationEfficiency;
	Am241SeparationEfficiency.Add(95,241,0,0.98);//Efficiency of 98%
	IsotopicVector Am243SeparationEfficiency;
	Am243SeparationEfficiency.Add(95,243,0,0.98);//Efficiency of 98%


	//Separation 2
	SeparationPlant* MySeparation2 = new SeparationPlant ( Logger ) ;
	MySeparation2->SetName("Separation2");

	double AvoidReprocessingUntilYear = 1990*year;
	MySeparation2->SetBackEndDestination( MyStorage3	, Am241SeparationEfficiency	,	AvoidReprocessingUntilYear	);
		//Extract 241Am from incoming material and send it to MyStorage3,  wait year 1990 to begin extraction
	MySeparation2->SetBackEndDestination( MyStorage4	, Am243SeparationEfficiency	,	AvoidReprocessingUntilYear	);
		//Extract 243Am from incoming material and send it to MyStorage4,  wait year 1990 to begin extraction
	//The material flux wich not goes to a BackeEnDestination is sent to the WASTE
	gCLASS->Add(MySeparation2);

	//Separation 1
	SeparationPlant* MySeparation1 = new SeparationPlant ( Logger ) ;
	MySeparation1->SetName("Separation1");


	MySeparation1->SetBackEndDestination( MyPool1	, UraniumSeparationEfficiency	,	0	);
	//Extract uranium from incoming material and send it to MyPool1, Don't wait before extraction 
	MySeparation1->SetBackEndDestination( MyStorage2 , PlutoniumSeparationEfficiency , 0 ) ;
	//Extract Plutonium from incoming material and send it to MyStorage2 , Don't wait before extraction 
	MySeparation1->SetBackEndDestination( MySeparation2 , AmericiumSeparationEfficiency , AvoidReprocessingUntilYear ) ;
	//Extract Americium from incoming material and send it to MySeparation2 , wait year 1990 to begin extraction 
	gCLASS->Add(MySeparation1);


	/*===A Reactor : PWR_UOX===*/
	double  HMMass = 72.5;//heavy metal mass (in tons)
	double	Power_CP0 = 2660e6;//Thermal power (in W)
	double  BurnUp = 33; //33 GWd/tHM

	cSecond StartingTime =  1978*year;
	cSecond LifeTime     =  40*year;
					
	Reactor* PWR_UOX = new Reactor(gCLASS->GetLog(),	//Log
							   STD900,				// Data base
							   MySeparation1,			// Connected Backend facility : The reactor discharge its fuel into the SeparationPlant1 "MySeparation1"
							   StartingTime,		// Starting time
							   LifeTime,			// time of reactor life time
							   Power_CP0,			// Power
							   HMMass,// HM mass
							   BurnUp,				// BurnUp
							   0.8);			// Load Factor
					

	PWR_UOX->SetName("a_PWR_Reactor");// name of the reactor (as it will show up in the CLASSGui)
	gCLASS->AddReactor(PWR_UOX);//Add this reactor to the scenario
					
	gCLASS->Evolution((double)year*2018);//Perform the calculation from year 1977(defined in Scenario declaration) to year 2018

	delete gCLASS;

}


//==========================================================================================
// Compilation
//==========================================================================================
/*
 
 \rm CLASS* ; g++ -o CLASS_Exec Separation.cxx -I $CLASS_CFLAG
 
 */
