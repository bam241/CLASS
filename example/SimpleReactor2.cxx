/*************************************************/
//              DESCRIPTION
// Simple scenario :
// This park is constituted by a simple PWR UOX
// Reactor, a pool and a Storage.
//   _______     ____    _______
//  |       |   |    |  |       | 
//  |Reactor| =>|Pool|=>|Storage| 
//  |_______|   |____|  |_______| 
//
// The spent fuel goes to the pool for 5 y
// then it goes to the Storage
//The scenario is run for 40 years
/*************************************************/
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
	int Std_output_level=0; //Only error are shown in terminal
	int File_output_level=2; // Error + Warning + Info are shown in the file CLASS_OUTPUT.log
	CLASSLogger *Logger		 = new CLASSLogger("CLASS_OUTPUT.log",Std_output_level,File_output_level);

	/******SCENARIO**********************************/
	// The scenario start at year 1977
	Scenario *gCLASS=new Scenario(1977*year,Logger);
	gCLASS->SetStockManagement(true);//If false all the IsotopicVector in stocks are mixed together.
	gCLASS->SetTimeStep(year/4.);//the scenario calculation is updated every 3 months

	/******DATA BASES**********************************/
	/*===Decay data base===*/
	//The decay data base is taken from the file Decay.idx
	DecayDataBank* DecayDB = new DecayDataBank(gCLASS->GetLog(), "../DATA_BASES/DECAY/ALL/Decay.idx");
	gCLASS->SetDecayDataBase(DecayDB);//This decay data base will be used for all the decay calculations in this Scenario

	/*===Reactor data base===*/
	//The file STD900.dat correspond to a fuel evolution of a UOX PWR (see manual for details)
	EvolutionData *STD900 = new EvolutionData(gCLASS->GetLog(), "../DATA_BASES/PWR/UOX/STD900/STD900.dat");
	/******FACILITIES*********************************/
	/*===A Stock===*/
	Storage *Stock = new Storage(gCLASS->GetLog()); //Definition of the stock
	Stock->SetName("Stock_UOX"); //Its name
	Stock->AddToStock(ZAI(92,238,0) * 10);// Just to illustrate, here we add ten nuclei of 238U in this stock
	gCLASS->Add(Stock); //Adding the stock to the Scenario 

	/*===A Pool===*/
	Pool *Cooling_UOX = new Pool(gCLASS->GetLog(),Stock, 5*year); //After 5 years of cooling, the pool sends its content to "Stock"
	Cooling_UOX->SetName("Pool_UOX");
	gCLASS->Add(Cooling_UOX);

	/*===A Reactor : PWR_UOX===*/
	double  HMMass = 72.5;//heavy metal mass (in tons)
	double	Power_CP0 = 2660e6;//Thermal power (in W)
	double  BurnUp = 33; //33 GWd/tHM

	cSecond StartingTime =  1978*year;
	cSecond LifeTime     =  40*year;
					
	Reactor* PWR_UOX = new Reactor(gCLASS->GetLog(),	//Log
							   STD900,				// Data base
							   Cooling_UOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
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
 
 \rm CLASS* ; g++ -o CLASS_Exec SimpleReactor2.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
 
 
 */
