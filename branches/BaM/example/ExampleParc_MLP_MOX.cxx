/************************************************************/
//              DESCRIPTION
// 
// 
// Simple Scenario with 8 PWR UOX and one PWR MOX
// 
// 
//   _______     ____    ______    	_______________
//  | 10 x  |   |    |  |       |      |		|
//  |Reactor| =>|Pool|=>|Storage|=====>|FabricationPlant| 
//  | UOx   |	|UOX |	| UOX 	|      |________________|
//  |_______|   |____|  |_______|   		||
//  						||
//						\/
//		 ______  	 ____ 	  _______
//		|       |	|    |	 |  1 x  |
//		|Storage|<===   |Pool|<==|Reactor|
//		| MOX 	|	|MOX |	 | MOX	 |
//		|_______|	|____|	 |_______|
//
//
//@author BaL
/***********************************************************/
#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include "XS/XSM_MLP.hxx"			//Load the include for Neural network cross section predictor
#include "Irradiation/IM_RK4.hxx"		//Load the include for Runge Kutta 4 resolution
#include "Equivalence/EQM_PWR_MLP_MOX.hxx"	//Load the include for Neural Network Equivalence Model (PWRMOX)
using namespace std;

int main(int argc, char** argv)
{
	//seconds in one year
	cSecond year = 3600*24.*365.25; 
	/******LOG MANAGEMENT**********************************/
	//Definition of the Log file : CLASS messages output 
	int Std_output_level 	= 0;  // Only error are shown in terminal
	int File_output_level 	= 2; // Error + Warning + Info are shown in the file CLASS_OUTPUT.log
	CLASSLogger *Logger 	= new CLASSLogger("CLASS_OUTPUT.log",Std_output_level,File_output_level);

	/******SCENARIO**********************************/
	// The scenario start at year 1977
	Scenario *gCLASS=new Scenario(1977*year,Logger);
	gCLASS->SetStockManagement(true);					//If false all the IsotopicVector in stocks are mixed together.
	gCLASS->SetTimeStep(year/4.);	 					//the scenario calculation is updated every 3 months
	cSecond EndOfScenarioTime=2040*year;				//Scenario ends in year 2040
	gCLASS->SetOutputFileName("ExampleParc_MLP.root");	//Set the name of the output file

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
	DecayDataBank* DecayDB = new DecayDataBank(gCLASS->GetLog(), PATH_TO_DATA + "DECAY/ALL/Decay.idx"); //you may have to open this file and do the proper changes according your path
	gCLASS->SetDecayDataBase(DecayDB);//This decay data base will be used for all the decay calculations in this Scenario

	/*===Reactor data base===*/

		// Reprocessed fuel PWR MOX
	XSM_MLP* XSMOX = new XSM_MLP(gCLASS->GetLog(), PATH_TO_DATA + "PWR/MOX/XSModel/30Wg_FullMOX");//Defining the XS Predictor
	IM_RK4 *IMRK4 = new IM_RK4(gCLASS->GetLog());	//Bateman's equation solver method (RungeKutta4)
	EQM_PWR_MLP_MOX* EQMMLPPWRMOX = new EQM_PWR_MLP_MOX(gCLASS->GetLog(),PATH_TO_DATA + "PWR/MOX/EQModel/MLP/EQM_MLP_PWR_MOX_3batch.xml");//Defining the EquivalenceModel
	PhysicsModels* PHYMOD = new PhysicsModels(XSMOX, EQMMLPPWRMOX, IMRK4); 							//The PhysicsModels containing the 3 object previously defined

		//Fixed fuel : PWR UOX
	EvolutionData *CYCLADE = new EvolutionData(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/FixedFuel/CYCLADES.dat");
	EvolutionData *GARANCE = new EvolutionData(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/FixedFuel/GARANCE.dat");
	EvolutionData *STD900 = new EvolutionData(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/FixedFuel/STD900.dat");

	/******FACILITIES*********************************/
	/*=== Stock===*/
	//Storage for UOX
	Storage *StockUOX = new Storage(gCLASS->GetLog());	// Definition of the stock
	StockUOX->SetName("StockUOX");				// Its name
	gCLASS->Add(StockUOX);					//Adding the stock to the Scenario
	//Storage for MOX
	Storage *StockMOX = new Storage(gCLASS->GetLog());	// Definition of the stock
	StockMOX->SetName("StockMOX"); 				// Its name
	gCLASS->Add(StockMOX);					//Adding the stock to the Scenario

	/*===Pool===*/
	//Pool for UOX
	Pool *Cooling_UOX = new Pool(gCLASS->GetLog(),StockUOX, 5*year); //After 5 years of cooling, the pool sends its content to "StockUOX"
	Cooling_UOX->SetName("Pool_UOX");
	gCLASS->Add(Cooling_UOX);
	//Pool for MOX
	Pool *Cooling_MOX = new Pool(gCLASS->GetLog(),StockMOX, 5*year); //After 5 years of cooling, the pool sends its content to "StockMOX"
	Cooling_MOX->SetName("Pool_MOX");
	gCLASS->Add(Cooling_MOX);

	/*===A FabricationPlant===*/
	FabricationPlant *FP_MOX = new FabricationPlant(gCLASS->GetLog(), 3*year); //Declare a FabricationPlant. After the build of the fuel, it decays during 3years before to be loaded in Reactor
	FP_MOX->SetFiFo(false); //The latest isotopicVector to enter in "Stock" will be used to build the fuel (Opposite of First In First Out)
	FP_MOX->SetName("Fab_MOX");
	FP_MOX->AddFissileStorage(StockUOX);	//Tell the FP to look in StockUOX for fissionable material 
	//FP_MOX->AddFertileStorage(Stock2);//Tell the FP to look in Stock2 for fertile material 
	//If fertile stock is not defined (like here), CLASS get fertile from nature (OUTCOMING vector)
	//FP_MOX->SetReUsableStorage(wastestock);//By default the fabricationplant get the list of nuclei defined in the EquivalenceModel (here EQM_MLP_MOX) from stock and send the others nuclei in WASTE. If user want these nuclei to go in another stock  he can use the SetReUsableStorage function
	gCLASS->AddFabricationPlant(FP_MOX);


	/*=== Reactors : Pallier CP0===*/
	double	Power_CP0 = 2.66e9;	    //Thermal power (in W)

		//Combustibles type : CYCLADE
		double  BurnUp_Cyclade    = 47; 		// GWd/tHM
		double  HMMass_Cyclade    = 72.3;		//heavy metal mass (in tons)
	
			//Fessenheim power plant
				//reactor n°1
				cSecond StartingTime =  1978*year;
				cSecond LifeTime     =  (EndOfScenarioTime - StartingTime);
								
				Reactor* Fessenheim_1 = new Reactor(gCLASS->GetLog(),// Log
										   CYCLADE,			// The DataBase used
										   Cooling_UOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor life time
										   Power_CP0,			// Power
										   HMMass_Cyclade,		// HM mass
										   BurnUp_Cyclade,		// BurnUp
										   0.8);			// Load Factor
								
			
				Fessenheim_1->SetName("Fessenheim_1");// name of the reactor (as it will show up in the CLASSGui)
				gCLASS->AddReactor(Fessenheim_1);//Add this reactor to the scenario
								
				//reactor n°2
				StartingTime =  1978*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
								
				Reactor* Fessenheim_2 = new Reactor(gCLASS->GetLog(),// Log
										   CYCLADE,			// The DataBase used
										   Cooling_UOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor life time
										   Power_CP0,			// Power
										   HMMass_Cyclade,		// HM mass
										   BurnUp_Cyclade,		// BurnUp
										   0.8);			// Load Factor
								
			
				Fessenheim_2->SetName("Fessenheim_2");// name of the reactor (as it will show up in the CLASSGui)
				gCLASS->AddReactor(Fessenheim_2);//Add this reactor to the scenario

			//Bugey power plant
				//reactor n°2
				StartingTime =  1979*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Bugey_2 = new Reactor(gCLASS->GetLog(),// Log
										   CYCLADE,			// The DataBase used
										   Cooling_UOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor life time
										   Power_CP0,			// Power
										   HMMass_Cyclade,		// HM mass
										   BurnUp_Cyclade,		// BurnUp
										   0.8);			// Load Factor
								
			
				Bugey_2->SetName("Bugey_2");// name of the reactor (as it will show up in the CLASSGui)
				gCLASS->AddReactor(Bugey_2);//Add this reactor to the scenario
				//reactor n° 3
				StartingTime =  1979*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Bugey_3 = new Reactor(gCLASS->GetLog(),// Log
										   CYCLADE,			// The DataBase used
										   Cooling_UOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor life time
										   Power_CP0,			// Power
										   HMMass_Cyclade,		// HM mass
										   BurnUp_Cyclade,		// BurnUp
										   0.8);			// Load Factor
								
			
				Bugey_3->SetName("Bugey_3");// name of the reactor (as it will show up in the CLASSGui)
				gCLASS->AddReactor(Bugey_3);//Add this reactor to the scenario

				//reactor n° 4
				StartingTime =  1979*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Bugey_4 = new Reactor(gCLASS->GetLog(),// Log
										   CYCLADE,			// The DataBase used
										   Cooling_UOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor life time
										   Power_CP0,			// Power
										   HMMass_Cyclade,		// HM mass
										   BurnUp_Cyclade,		// BurnUp
										   0.8);			// Load Factor
								
			
				Bugey_4->SetName("Bugey_4");// name of the reactor (as it will show up in the CLASSGui)
				gCLASS->AddReactor(Bugey_4);//Add this reactor to the scenario

				//reactor n° 5
				StartingTime =  1980*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Bugey_5 = new Reactor(gCLASS->GetLog(),// Log
										   CYCLADE,			// The DataBase used
										   Cooling_UOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor life time
										   Power_CP0,			// Power
										   HMMass_Cyclade,		// HM mass
										   BurnUp_Cyclade,		// BurnUp
										   0.8);			// Load Factor
								
			
				Bugey_5->SetName("Bugey_5");// name of the reactor (as it will show up in the CLASSGui)
				gCLASS->AddReactor(Bugey_5);//Add this reactor to the scenario

	/*=== Reactors : Pallier CPY===*/
	double	Power_CPY = 2.785e9;	    //Thermal power (in W)
	 //Combustibles type : GARANCE
		double  BurnUp_GARANCE    = 42; 		// GWd/tHM
		double  HMMass_GARANCE    = 72.3;		//heavy metal mass (in tons)
	

			// Gravelines power plant
				//reactor n° 1
				StartingTime =  1980*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Gravelines_1 = new Reactor(gCLASS->GetLog(),// Log
										   GARANCE,			// The DataBase used
										   Cooling_UOX,			// Connected Backend
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor l
										   Power_CPY,			// Power
										   HMMass_GARANCE,		// HM mass
										   BurnUp_GARANCE,		// BurnUp
										   0.8);			// Load Factor
								
			
				Gravelines_1->SetName("Gravelines_1");// name of the reactor (as it will show 
				gCLASS->AddReactor(Gravelines_1);//Add this reactor to the scenario

		//Combustibles type : STANDARD 900
		double  BurnUp_STD900    = 33; 		// GWd/tHM
		double  HMMass_STD900    = 72.3;		//heavy metal mass (in tons)

				//reactor n° 2
				StartingTime =  1980*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Gravelines_2 = new Reactor(gCLASS->GetLog(),// Log
										   STD900,			// The DataBase used
										   Cooling_UOX,			// Connected Backend
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor l
										   Power_CPY,			// Power
										   HMMass_STD900,		// HM mass
										   BurnUp_STD900,		// BurnUp
										   0.8);			// Load Factor
								
			
				Gravelines_2->SetName("Gravelines_2");// name of the reactor (as it will show
				gCLASS->AddReactor(Gravelines_2);//Add this reactor to the scenario

			// Tricastin power plant

				//reactor n° 1
				StartingTime =  1980*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Tricastin_1 = new Reactor(gCLASS->GetLog(),// Log
										   STD900,			// The DataBase used
										   Cooling_UOX,			// Connected Backend
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor l
										   Power_CPY,			// Power
										   HMMass_STD900,		// HM mass
										   BurnUp_STD900,		// BurnUp
										   0.8);			// Load Factor
								
			
				Tricastin_1->SetName("Tricastin_1");// name of the reactor (as it will ll sh
				gCLASS->AddReactor(Tricastin_1);//Add this reactor to the scenario

				//reactor n° 2
				StartingTime =  1980*year;
				LifeTime     =  (EndOfScenarioTime - StartingTime);
				Reactor* Tricastin_2 = new Reactor(gCLASS->GetLog(),// Log
										   STD900,			// The DataBase used
										   Cooling_UOX,			// Connected Backend
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor l
										   Power_CPY,			// Power
										   HMMass_STD900,		// HM mass
										   BurnUp_STD900,		// BurnUp
										   0.8);			// Load Factor
								
			
				Tricastin_2->SetName("Tricastin_2");// name of the reactor (as it will l
				gCLASS->AddReactor(Tricastin_2);//Add this reactor to the scenario




				// Dampierre Power Plant UOX puis MOX
				//UOX
				cSecond Dampierre_MOX_Time		 =  1991*year;

				StartingTime =  1980*year;
				LifeTime     =  Dampierre_MOX_Time - StartingTime;
				double BunrUpMOX = 35;
				Reactor* Dampierre_UOX = new Reactor(gCLASS->GetLog(),// Log
										   STD900,			// The DataBase used
										   Cooling_UOX,			// Connected Backend
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor l
										   Power_CPY,			// Power
										   HMMass_STD900,		// HM mass
										   BurnUp_STD900,		// BurnUp
										   0.8);			// Load Factor


				Dampierre_UOX->SetName("Dampierre_UOX");// name of the reactor (as it will l
				gCLASS->AddReactor(Dampierre_UOX);//Add this reactor to the scenario

				//the PWR MOX
				StartingTime =  Dampierre_MOX_Time;
				LifeTime     =  EndOfScenarioTime - StartingTime;
				Reactor* Dampierre_MOX = new Reactor(gCLASS->GetLog(),// Log
							 			   PHYMOD,			// The models used to build the fuel & to calculate its evolution
							 			   FP_MOX,			// The FabricationPlant
										   Cooling_MOX,			// Connected Backend
										   StartingTime,		// Starting time
										   LifeTime,			// time of reactor l
										   Power_CPY,			// Power
										   HMMass_STD900,		// HM mass
										   BunrUpMOX,			// BurnUp
										   0.8);			// Load Factor

				Dampierre_MOX->SetName("Dampierre_MOX");// name of the reactor (as it will l
				gCLASS->AddReactor(Dampierre_MOX);//Add this reactor to the scenario



	gCLASS->Evolution((double)EndOfScenarioTime);//Perform the calculation from year 1977(defined in Scenario declaration) to year 2018

	delete gCLASS;

}


//==========================================================================================
// Compilation
//==========================================================================================
/*
 
 \rm CLASS* ; g++ -o CLASS_Exec ExampleParc_MLP_MOX.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
 
 
 */
