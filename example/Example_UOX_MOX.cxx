/*************************************************/
//              DESCRIPTION
// Simple scenario :
// This park is constituted by a PWR UOX, a PWR MOx for PERMIS exo 1
// Pool, FP, and storage
//   _______     _______      _______      _______      _______     _______
//  |       |   |       |    |       |    |       |    |       |   |       |
//  |Reactor| =>|  POOL | => | Stock | => |Reactor| => | POOL  | =>| Stock | 
//  | UOX   |   |  UOX  |    |  UOX  | FP | MOX   |    |  MOX  |   |  MOX  |
//  |_______|   |_______|    |_______|    |_______|    |_______|   |_______|
//
//
// Some reactor and fleet parameters are in argument :
// Burn-Up UOX
// Burn-Up MOX
// Fraction of Reactor MOX in the fleet
// Cooling time in the pool
// Stock Management (1 = LIFO ; 2 FIFO ; 3 MIX ; 4 Random)
//
//
//@author Nico
/*************************************************/
#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>

#include "XSM_MLP.hxx"			
#include "IM_RK4.hxx"	
#include "EQ_OneParameter.hxx"

using namespace std;

int main(int argc, char** argv)
{

	cSecond year = 3600*24.*365.25; 

	//##########################################################################################
	//####### USE ##############################################################################
	//##########################################################################################

	if (argc != 6)
	{
		cout<<"#############################################"<<endl;
		cout<<"#############################################"<<endl<<endl;
		cout<<"USE : "<<endl<<endl;
		cout<<"CLASS_Exec BU_UOx BU_MOX FractionMOX TCooling StockManagment"<<endl<<endl;;
		cout<<"EXAMPLE : "<<endl<<endl;
		cout<<"CLASS_Exec 40 45 0.10 5 1"<<endl<<endl;;
		cout<<"#############################################"<<endl<<endl;;
		cout<<"#############################################"<<endl;
		exit(1);
	}

	//##########################################################################################
	//####### SCENARIO DATA ####################################################################
	//##########################################################################################

	cSecond d_TimeStep 		= year / 12.;
	cSecond d_TimeScenario	= 100.*year  ;
	string s_FileName 		= string("OUT.root");

	//##########################################################################################
	//####### PARAMETERS #######################################################################
	//##########################################################################################

	double BU_UOX 			= atof(argv[1]);
	double BU_MOX 			= atof(argv[2]);
	double FractionMOX 		= atof(argv[3]);
	cSecond TCooling 		= (cSecond)(atof(argv[4])*year);

	int StockManagment 		= atoi(argv[5]); // 1 = kpLiFo, 2 = kpFiFo, 3 = kpMix
	StorageManagement SM;
	if (StockManagment == 1)	SM = kpLiFo;
	else if (StockManagment == 2)	SM = kpFiFo;
	else if (StockManagment == 3)	SM = kpMix;
	else if (StockManagment == 4)	SM = kpRand;
	else {cout<<endl<<"StockManagement should be 1, 2, 3 or 4... EXIT"<<endl<<endl; exit(1);}

	double LoadFactorUOx 	= 0.75;
	double LoadFactorMOx 	= 0.75;
	double KseuilUOx		= 1.03;
	double KseuilMOx		= 1.03;
	int NumberOfBatchesUOx	= 3;
	int NumberOfBatchesMOx	= 3;
	cSecond TFabMOX			= (cSecond)(2*year);

	double N_UOx 		= 1;
	double HMMassUOx	= 72.3	 * N_UOx;
	double PowerUOx 	= 2785e6 * N_UOx;

	cSecond StartingTimeUOx	=  0*year;
	cSecond LifeTimeUOx		=  d_TimeScenario;

	double N_MOx 		= FractionMOX;
	double HMMassMOx	= 72.3	 * N_MOx;
	double PowerMOx 	= 2785e6 * N_MOx;

	cSecond TimeStartMox 	= 20.*year;
	cSecond LifeTimeMOx		=  d_TimeScenario;

	//##########################################################################################
	//####### LOG MANAGEMENT ###################################################################
	//##########################################################################################

	int Std_output_level 	= 0;
	int File_output_level 	= 0;
	CLASSLogger *Logger 	= new CLASSLogger("CLASS_OUTPUT.log",Std_output_level,File_output_level);

	//##########################################################################################
	//####### SCENARIO #########################################################################
	//##########################################################################################

	Scenario *gCLASS=new Scenario(0*year,Logger);
	gCLASS->SetStockManagement(true);
	gCLASS->SetTimeStep(d_TimeStep);
	gCLASS->SetOutputFileName(s_FileName);
	gCLASS->SetZAIThreshold(82);

	//##########################################################################################
	//####### DATABASE #########################################################################
	//##########################################################################################

    string CLASS_PATH = getenv("CLASS_PATH");
    string PATH_TO_DATA = CLASS_PATH + "/DATA_BASES/";

   	// DECAY
	DecayDataBank* DecayDB = new DecayDataBank(gCLASS->GetLog(), PATH_TO_DATA + "DECAY/ALL/Decay.idx");
	gCLASS->SetDecayDataBase(DecayDB);

	// Bateman Solver
	IM_RK4*			IMRK4 		= new IM_RK4(gCLASS->GetLog());

	// REP UOX
	XSM_MLP* 		XSMUOX 			= new XSM_MLP(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/XSModel/30Wg_FullUOX");
	EQ_OneParameter*	EQMPWRUOX 	= new EQ_OneParameter(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/EQModel/XML/PWR_UOX.xml", PATH_TO_DATA + "PWR/UOX/EQModel/NFO/PWR_UOX.nfo");
	EQMPWRUOX->SetModelParameter("kThreshold",KseuilUOx);
	EQMPWRUOX->SetModelParameter("NumberOfBatch",NumberOfBatchesUOx);
	PhysicsModels*		PMUOX		= new PhysicsModels(XSMUOX, EQMPWRUOX, IMRK4);

	// REP MOx 
	XSM_MLP* 			XSMMOX 		= new XSM_MLP(gCLASS->GetLog(), PATH_TO_DATA + "PWR/MOX/XSModel/30Wg_FullMOX");
	EQ_OneParameter* 	EQMPWRMOX 	= new EQ_OneParameter(gCLASS->GetLog(), PATH_TO_DATA + "PWR/MOX/EQModel/XML/PWR_MOX.xml", PATH_TO_DATA + "PWR/MOX/EQModel/NFO/PWR_MOX.nfo");
	EQMPWRMOX->SetModelParameter("kThreshold",KseuilMOx);
	EQMPWRMOX->SetModelParameter("NumberOfBatch",NumberOfBatchesMOx);
	PhysicsModels*		PMMOX		= new PhysicsModels(XSMMOX, EQMPWRMOX, IMRK4);

	//##########################################################################################
	//####### IV ###############################################################################
	//##########################################################################################

	//##########################################################################################
	//####### STORAGE ##########################################################################
	//##########################################################################################
 
	Storage *StockUOx = new Storage(gCLASS->GetLog());
	StockUOx->SetName("StockUOx");
	gCLASS->Add(StockUOx);

	Storage *StockMOx = new Storage(gCLASS->GetLog());
	StockMOx->SetName("StockMOx");
	gCLASS->Add(StockMOx);

	//##########################################################################################
	//####### SEPARATION PLANT #################################################################
	//##########################################################################################

	//##########################################################################################
	//####### POOL #############################################################################
	//##########################################################################################

	Pool *PoolUOx = new Pool(gCLASS->GetLog(),StockUOx, TCooling);
	PoolUOx->SetName("PoolUOx");
	gCLASS->Add(PoolUOx);

	Pool *PoolMOx = new Pool(gCLASS->GetLog(),StockMOx, 5*year);
	PoolMOx->SetName("PoolMOx");
	gCLASS->Add(PoolMOx);

	//##########################################################################################
	//####### FABRICATION PLANT #################################################################
	//##########################################################################################

	FabricationPlant *FP_UOX = new FabricationPlant(gCLASS->GetLog(), 0*year);
	FP_UOX->SetName("Fab_UOX");
	FP_UOX->AddInfiniteStorage("Fissile",0.02,0.06,1);
	FP_UOX->AddFuelBuffer("Fertile");
	gCLASS->AddFabricationPlant(FP_UOX);

	FabricationPlant *FP_MOX = new FabricationPlant(gCLASS->GetLog(), 2.*year);
	FP_MOX->SetName("FP_MOX");
	FP_MOX->SetSeparationManagement(true);
	FP_MOX->SetStorageManagement(SM);
	FP_MOX->AddStorage("Fissile", StockUOx, 0.04, 0.16,1);
	FP_MOX->AddFuelBuffer("Fertile");
	gCLASS->AddFabricationPlant(FP_MOX);

	//##########################################################################################
	//####### REP UOX ##########################################################################
	//##########################################################################################

	double HMMass 		= HMMassUOx;
	double Power 		= PowerUOx;
	double BurnUp 		= BU_UOX;
	double LoadFactor 	= LoadFactorUOx;

	cSecond StartingTime =  StartingTimeUOx;
	cSecond LifeTime     =  LifeTimeUOx;
	
	Reactor* PWR_UOX = new Reactor(gCLASS->GetLog(),PMUOX,FP_UOX,PoolUOx,StartingTime,LifeTime,Power,HMMass,BurnUp,LoadFactor);
	PWR_UOX->AddScheduleEntry(TimeStartMox, PMUOX, BurnUp, Power*(1 - FractionMOX), HMMass*(1 - FractionMOX));

	PWR_UOX->SetName("PWR_UOx");
	gCLASS->AddReactor(PWR_UOX);

	//##########################################################################################
	//####### REP MOX ##########################################################################
	//##########################################################################################

	HMMass 		= HMMassMOx;
	Power 		= PowerMOx;
	BurnUp 		= BU_MOX;
	LoadFactor 	= LoadFactorMOx;

	StartingTime =  TimeStartMox;
	LifeTime     =  LifeTimeMOx;
	
	Reactor* PWR_MOX = new Reactor(gCLASS->GetLog(),PMMOX,FP_MOX,PoolMOx,StartingTime,LifeTime,Power,HMMass,BurnUp,LoadFactor);

	PWR_MOX->SetName("PWR_MOX");
	gCLASS->AddReactor(PWR_MOX);

	//##########################################################################################
	//####### EVOLUTION ########################################################################
	//##########################################################################################

	gCLASS->Evolution((double)d_TimeScenario);

	//##########################################################################################
	//####### OUTPUT ###########################################################################
	//##########################################################################################

	delete gCLASS;
}

//==========================================================================================
// Compilation
//==========================================================================================
/*
 
 g++ -o CLASS_Exec Example_UOX_MOX.cxx $CLASS_CFLAG
 
 
 */
