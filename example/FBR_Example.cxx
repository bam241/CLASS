/************************************************************/
//              DESCRIPTION
// Close Fuel cycle scenario :
// This park is constituted by a FBR-Na MOX which
// multi-recycle its own fuel.
// The Storage is initially filled with Pu in order
// to this scenario to be doable
//         _______________	_______     ____    _______
//        |		   |   | FBR   |   |    |  |       |
//  ||===>|FabricationPlant| =>|Reactor| =>|Pool|=>|Storage|===||
//  ||    |________________|   |_______|   |____|  |_______|   ||
//  ||=========================================================||
//
// The spent fuel goes to the pool for 5 y
// then it goes to the Storage
//
//@author BaL
/***********************************************************/
#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include "XS/XSM_MLP.hxx"				//Load the include for Neural network cross section predictor
#include "Irradiation/IM_RK4.hxx"			//Load the include for Runge Kutta 4 resolution
#include "Equivalence/EQM_BakerRoss_FBR_MOX.hxx"	//Load the include for Neural Network Equivalence Model (PWRMOX)
using namespace std;

int main(int argc, char** argv)
{
	//seconds in one year
	cSecond year = 3600*24.*365.25; 
	/******LOG MANAGEMENT**********************************/
	//Definition of the Log file : CLASS messages output 
	int Std_output_level=0;  // Only error are shown in terminal
	int File_output_level=3; // Error + Warning + Info are shown in the file CLASS_OUTPUT.log
	CLASSLogger *Logger	 = new CLASSLogger("CLASS_OUTPUT.log",Std_output_level,File_output_level);

	/******SCENARIO**********************************/
	// The scenario start at year 1977
	Scenario *gCLASS=new Scenario(1977*year,Logger);
	gCLASS->SetStockManagement(true);//If false all the IsotopicVector in stocks are mixed together.
	gCLASS->SetTimeStep(year/4.);	 //the scenario calculation is updated every 3 months
	gCLASS->SetOutputFileName("FBR_Example.root");	//Set the name of the output file

	/******DATA BASES**********************************/
	/*===Decay data base===*/
	//The decay data base is taken from the file Decay.idx
	DecayDataBank* DecayDB = new DecayDataBank(gCLASS->GetLog(), "../DATA_BASES/DECAY/ALL/Decay.idx");
	gCLASS->SetDecayDataBase(DecayDB);//This decay data base will be used for all the decay calculations in this Scenario

	/*===Reactor data base===*/

	XSM_MLP* XS_FBRMOX = new XSM_MLP(gCLASS->GetLog(), "../DATA_BASES/FBR_Na/MOX/XSModel/ESFR_48Wg");//Defining the XS Predictor
	IM_RK4 *IMRK4 = new IM_RK4(gCLASS->GetLog());							 //Bateman's equation solver method (RungeKutta4)
	IMRK4->SetSpectrumType("fast");									 //Set the spectrum to fast for reactions isomeric branching ratios (can be fast or thermal)
	IMRK4->LoadFPYield("","../data/FPyield_Fast_JEFF3.1.dat");						 //Add the handling of fission procuct and gets fission yields from this file (the first argument is for spontaneousfission yield : here is not handle)
	
	EQM_BakerRoss_FBR_MOX* EQM_FBRMOX = new EQM_BakerRoss_FBR_MOX(gCLASS->GetLog());//Defining the EquivalenceModel
	EQM_FBRMOX->SetBuildFuelFirstGuess(0.12);					//Set the first guess of fissile content to 12 per cent of plutonium (default : 5%)
	PhysicsModels* PHYMOD = new PhysicsModels(XS_FBRMOX, EQM_FBRMOX, IMRK4);	//The PhysicsModels containing the 3 object previously defined


	/******FACILITIES*********************************/
	/*===A Stock===*/
	Storage *Stock = new Storage(gCLASS->GetLog()); // Definition of the stock
	Stock->SetName("Stock"); 			// Its name
	//Fill the stock with an initial amount of Plutonium
	// In order to allow the FBR_MOX to work
	IsotopicVector InitialIV;
	InitialIV.Add(94,238,0,4e27  );
	InitialIV.Add(94,239,0,6.4e28);
	InitialIV.Add(94,240,0,1.6e28);
	InitialIV.Add(94,241,0,9.0e27);
	InitialIV.Add(94,242,0,6e27  );
	InitialIV.Add(95,241,0,1e27  );
	gCLASS->Add(Stock); 	//Adding the stock to the Scenario 
	Stock->AddToStock(InitialIV);

	/*===A Pool===*/
	Pool *Cooling_MOX = new Pool(gCLASS->GetLog(),Stock, 5*year); //After 5 years of cooling, the pool sends its content to "Stock"
	Cooling_MOX->SetName("Pool_MOX");
	gCLASS->Add(Cooling_MOX);

	/*===A FabricationPlant===*/
	FabricationPlant *FP_MOX = new FabricationPlant(gCLASS->GetLog(), 3*year); //Declare a FabricationPlant. After the build of the fuel, it decays during 3years before to be loaded in Reactor
	FP_MOX->SetFiFo(false); //The latest isotopicVector to enter in "Stock" will be used to build the fuel (Opposite of First In First Out)
	FP_MOX->SetName("Fab_MOX");
	FP_MOX->AddFissileStorage(Stock);	//Tell the FP to look in Stock for fissionable material 
	//FP_MOX->AddFertileStorage(Stock2);//Tell the FP to look in Stock2 for fertile material 
	//If fertile stock is not defined (like here), CLASS get fertile from nature (OUTCOMING vector)
	//FP_MOX->SetReUsableStorage(wastestock);//By default the fabricationplant get the list of nuclei defined in the EquivalenceModel (here EQM_MLP_MOX) from stock and send the others nuclei in WASTE. If user want these nuclei to go in another stock  he can use the SetReUsableStorage function
	gCLASS->AddFabricationPlant(FP_MOX);


	/*===A Reactor : FBR_MOX===*/
	double  HMMass    = 7.48336500000000058e+01;	//heavy metal mass (in tons)
	double	Power_CP0 = 3.6e9;			//Thermal power (in W)
	double  BurnUp    = 100;			//100 GWd/tHM

	cSecond StartingTime =  1985*year;
	cSecond LifeTime     =  40*year;
					
	Reactor* FBR_MOX = new Reactor(gCLASS->GetLog(),				// Log
							   PHYMOD,			// The models used to build the fuel & to calculate its evolution
							   FP_MOX,			// The FabricationPlant
							   Cooling_MOX,			// Connected Backend facility : The reactor discharge its fuel into the Pool "Cooling_UOX"
							   StartingTime,		// Starting time
							   LifeTime,			// time of reactor life time
							   Power_CP0,			// Power
							   HMMass,			// HM mass
							   BurnUp,			// BurnUp
							   0.8);			// Reactor efficiency (% of time it is working @ full power)
					

	FBR_MOX->SetName("a_FBR_MOX");	// name of the reactor (as it will show up in the CLASSGui)
	gCLASS->AddReactor(FBR_MOX);	//Add this reactor to the scenario
					
	gCLASS->Evolution((double)year*2018);//Perform the calculation from year 1977(defined in Scenario declaration) to year 2018

	delete gCLASS;

}


//==========================================================================================
// Compilation
//==========================================================================================
/*
 
 \rm CLASS* ; g++ -o CLASS_Exec FBR_Example.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result
 
 
 */
