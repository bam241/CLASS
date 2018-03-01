#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>

#include "XSM_MLP.hxx"
#include "IM_RK4.hxx"
#include "EQ_OneParameter.hxx"

TEST ( Fleet, FullCalculation ) {

	cSecond year = 3600 * 24.*365.25;

//##########################################################################################
//####### WRITE TEST CHARACTERISTICS #######################################################
//##########################################################################################

	cout << endl;
	cout << "####################################################################" << endl;
	cout << "### GOOGLE TEST : test_PWR_UOX_MOX..." << endl;
	cout << "####################################################################" << endl << endl;

//##########################################################################################
//####### SCENARIO DATA ####################################################################
//##########################################################################################

	cSecond d_TimeStep 		= year / 12.;
	cSecond d_TimeScenario	= 50.*year  ;
	string s_FileName 		= string("OUT.root");

//##########################################################################################
//####### PARAMETERS #######################################################################
//##########################################################################################

	double BU_UOX 			= 40.;
	double BU_MOX 			= 45;
	double FractionMOX 		= 0.10;
	cSecond TCooling 		= (cSecond)(5 * year);

	StorageManagement SM = kpLiFo;

	double LoadFactorUOx 	= 0.75;
	double LoadFactorMOx 	= 0.75;
	double KseuilUOx		= 1.03;
	double KseuilMOx		= 1.03;
	int NumberOfBatchesUOx	= 3;
	int NumberOfBatchesMOx	= 3;
	cSecond TFabMOX			= (cSecond)(2 * year);

	double N_UOx 		= 1;
	double HMMassUOx	= 72.3	 * N_UOx;
	double PowerUOx 	= 2785e6 * N_UOx;

	cSecond StartingTimeUOx	=  0 * year;
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
	CLASSLogger *Logger 	= new CLASSLogger("CLASS_OUTPUT.log", Std_output_level, File_output_level);

//##########################################################################################
//####### SCENARIO #########################################################################
//##########################################################################################

	Scenario *gCLASS = new Scenario(0 * year, Logger);
	gCLASS->SetStockManagement(true);
	gCLASS->SetTimeStep(d_TimeStep);
	gCLASS->SetOutputFileName(s_FileName);
	gCLASS->SetZAIThreshold(82);

//##########################################################################################
//####### DATABASE #########################################################################
//##########################################################################################

// Global Path
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
	EQMPWRUOX->SetModelParameter("kThreshold", KseuilUOx);
	EQMPWRUOX->SetModelParameter("NumberOfBatch", NumberOfBatchesUOx);
	PhysicsModels*		PMUOX		= new PhysicsModels(XSMUOX, EQMPWRUOX, IMRK4);

// REP MOx
	XSM_MLP* 			XSMMOX 		= new XSM_MLP(gCLASS->GetLog(), PATH_TO_DATA + "PWR/MOX/XSModel/30Wg_FullMOX");
	EQ_OneParameter* 	EQMPWRMOX 	= new EQ_OneParameter(gCLASS->GetLog(), PATH_TO_DATA + "PWR/MOX/EQModel/XML/PWR_MOX.xml", PATH_TO_DATA + "PWR/MOX/EQModel/NFO/PWR_MOX.nfo");
	EQMPWRMOX->SetModelParameter("kThreshold", KseuilMOx);
	EQMPWRMOX->SetModelParameter("NumberOfBatch", NumberOfBatchesMOx);
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

	Pool *PoolUOx = new Pool(gCLASS->GetLog(), StockUOx, TCooling);
	PoolUOx->SetName("PoolUOx");
	gCLASS->Add(PoolUOx);

	Pool *PoolMOx = new Pool(gCLASS->GetLog(), StockMOx, 5 * year);
	PoolMOx->SetName("PoolMOx");
	gCLASS->Add(PoolMOx);

//##########################################################################################
//####### FABRICATION PLANT #################################################################
//##########################################################################################

	FabricationPlant *FP_UOX = new FabricationPlant(gCLASS->GetLog(), 0 * year);
	FP_UOX->SetName("Fab_UOX");
	FP_UOX->AddInfiniteStorage("Fissile", 0.02, 0.06, 1);
	FP_UOX->AddFuelBuffer("Fertile");
	gCLASS->AddFabricationPlant(FP_UOX);

	FabricationPlant *FP_MOX = new FabricationPlant(gCLASS->GetLog(), 2.*year);
	FP_MOX->SetName("FP_MOX");
	FP_MOX->SetSeparationManagement(true);
	FP_MOX->SetStorageManagement(SM);
	FP_MOX->AddStorage("Fissile", StockUOx, 0.04, 0.16, 1);
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

	Reactor* PWR_UOX = new Reactor(gCLASS->GetLog(), PMUOX, FP_UOX, PoolUOx, StartingTime, LifeTime, Power, HMMass, BurnUp, LoadFactor);
	PWR_UOX->AddScheduleEntry(TimeStartMox, PMUOX, BurnUp, Power * (1 - FractionMOX), HMMass * (1 - FractionMOX));

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

	Reactor* PWR_MOX = new Reactor(gCLASS->GetLog(), PMMOX, FP_MOX, PoolMOx, StartingTime, LifeTime, Power, HMMass, BurnUp, LoadFactor);

	PWR_MOX->SetName("PWR_MOX");
	gCLASS->AddReactor(PWR_MOX);

//##########################################################################################
//####### EVOLUTION ########################################################################
//##########################################################################################

	gCLASS->Evolution((double)d_TimeScenario);
	delete gCLASS;

//##########################################################################################
//####### TEST #############################################################################
//##########################################################################################

	// Actinides number of atoms at End Of Scenario In Reactor PWR MOX
	double N_Pu9_RMOX_EOS_Expected = 8.21547e+26;
	double N_Np7_RMOX_EOS_Expected = 1.19262e+23;
	double N_Am1_RMOX_EOS_Expected = 1.52076e+25;
	double N_Cm4_RMOX_EOS_Expected = 3.04428e+22;

	double N_Pu9_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetQuantity(94, 239, 0);
	double N_Np7_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetQuantity(93, 237, 0);
	double N_Am1_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetQuantity(95, 241, 0);
	double N_Cm4_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetQuantity(96, 244, 0);

	// Actinides element masses at End Of Scenario In Reactor PWR MOX
	double M_Pu_RMOX_EOS_Expected = 0.561599;
	double M_Np_RMOX_EOS_Expected = 0.000401171;
	double M_Am_RMOX_EOS_Expected = 0.00666378;
	double M_Cm_RMOX_EOS_Expected = 0.000191061;

	double M_Pu_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetSpeciesComposition(94).GetTotalMass();
	double M_Np_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetSpeciesComposition(93).GetTotalMass();
	double M_Am_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetSpeciesComposition(95).GetTotalMass();
	double M_Cm_RMOX_EOS_Calculated = PWR_MOX->GetIVReactor().GetSpeciesComposition(96).GetTotalMass();

	// Actinides element masses at End Of Scenario In Spent UOX Storage
	double M_Pu_SUOX_EOS_Expected = 2.84185;
	double M_Np_SUOX_EOS_Expected = 0.17761;
	double M_Am_SUOX_EOS_Expected = 0.400192;
	double M_Cm_SUOX_EOS_Expected = 0.00503383;

	double M_Pu_SUOX_EOS_Calculated = StockUOx->GetInsideIV().GetSpeciesComposition(94).GetTotalMass();
	double M_Np_SUOX_EOS_Calculated = StockUOx->GetInsideIV().GetSpeciesComposition(93).GetTotalMass();
	double M_Am_SUOX_EOS_Calculated = StockUOx->GetInsideIV().GetSpeciesComposition(95).GetTotalMass();
	double M_Cm_SUOX_EOS_Calculated = StockUOx->GetInsideIV().GetSpeciesComposition(96).GetTotalMass();

	// Actinides element masses at End Of Scenario In Spent MOX Storage
	double M_Pu_SMOX_EOS_Expected = 1.95429;
	double M_Np_SMOX_EOS_Expected = 0.0104991;
	double M_Am_SMOX_EOS_Expected = 0.289419;
	double M_Cm_SMOX_EOS_Expected = 0.020052;

	double M_Pu_SMOX_EOS_Calculated = StockMOx->GetInsideIV().GetSpeciesComposition(94).GetTotalMass();
	double M_Np_SMOX_EOS_Calculated = StockMOx->GetInsideIV().GetSpeciesComposition(93).GetTotalMass();
	double M_Am_SMOX_EOS_Calculated = StockMOx->GetInsideIV().GetSpeciesComposition(95).GetTotalMass();
	double M_Cm_SMOX_EOS_Calculated = StockMOx->GetInsideIV().GetSpeciesComposition(96).GetTotalMass();

	// TEST
	ASSERT_NEAR(N_Pu9_RMOX_EOS_Expected, N_Pu9_RMOX_EOS_Calculated, N_Pu9_RMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(N_Np7_RMOX_EOS_Expected, N_Np7_RMOX_EOS_Calculated, N_Np7_RMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(N_Am1_RMOX_EOS_Expected, N_Am1_RMOX_EOS_Calculated, N_Am1_RMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(N_Cm4_RMOX_EOS_Expected, N_Cm4_RMOX_EOS_Calculated, N_Cm4_RMOX_EOS_Expected * 1e-5);

	ASSERT_NEAR(M_Pu_RMOX_EOS_Expected, M_Pu_RMOX_EOS_Calculated, M_Pu_RMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Np_RMOX_EOS_Expected, M_Np_RMOX_EOS_Calculated, M_Np_RMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Am_RMOX_EOS_Expected, M_Am_RMOX_EOS_Calculated, M_Am_RMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Cm_RMOX_EOS_Expected, M_Cm_RMOX_EOS_Calculated, M_Cm_RMOX_EOS_Expected * 1e-5);

	ASSERT_NEAR(M_Pu_SUOX_EOS_Expected, M_Pu_SUOX_EOS_Calculated, M_Pu_SUOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Np_SUOX_EOS_Expected, M_Np_SUOX_EOS_Calculated, M_Np_SUOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Am_SUOX_EOS_Expected, M_Am_SUOX_EOS_Calculated, M_Am_SUOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Cm_SUOX_EOS_Expected, M_Cm_SUOX_EOS_Calculated, M_Cm_SUOX_EOS_Expected * 1e-5);

	ASSERT_NEAR(M_Pu_SMOX_EOS_Expected, M_Pu_SMOX_EOS_Calculated, M_Pu_SMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Np_SMOX_EOS_Expected, M_Np_SMOX_EOS_Calculated, M_Np_SMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Am_SMOX_EOS_Expected, M_Am_SMOX_EOS_Calculated, M_Am_SMOX_EOS_Expected * 1e-5);
	ASSERT_NEAR(M_Cm_SMOX_EOS_Expected, M_Cm_SMOX_EOS_Calculated, M_Cm_SMOX_EOS_Expected * 1e-5);

//##########################################################################################
//####### WRITE OUTPUT #####################################################################
//##########################################################################################

	cout << endl;
	cout << "####################################################################" << endl;
	cout << "### GOOGLE TEST : END OF TEST" << endl;
	cout << "####################################################################" << endl << endl;

//##########################################################################################
//####### CLEANING #########################################################################
//##########################################################################################

	system("rm -rf CLASS_OUTPUT.log CLASS_TimeStep DecayDataBank.log OUT.root");
}
