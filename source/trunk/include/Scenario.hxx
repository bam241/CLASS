#ifndef _SCENARIO_HXX_
#define _SCENARIO_HXX_
/*!
 \file
 \brief Header file for CLASS classes. 
 */

#include "CLASSObject.hxx"
#include "IsotopicVector.hxx"

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <map>
#include <iostream>



using namespace std;
typedef long long int cSecond;

class DecayDataBank;
class FabricationPlant;
class SeparationPlant;
class Reactor;
class Pool;
class Storage;

//-----------------------------------------------------------------------------//
/*!
Define a Scenario .
The aim of these Scenario is to manage the park and its evolution and to lead all Storage, FabricationPlant, Reactor, Pool.


 @author BaM
 @author Marc
 @version 2.0
 */
//________________________________________________________________________


class Scenario : public CLASSObject
{
public :

//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	Scenario();	///< Normal Constructor.

	Scenario(CLASSLogger* Log, cSecond abstime = 0);	///< Log Constructor.
 	/*!
	 Use to load a CLASSLogger
	 \param CLASSLogger: CLASSLogger used for the log...
	 */
	Scenario(cSecond abstime);	///< Time Constructor.


	Scenario(cSecond abstime, CLASSLogger* log);	///< Time Constructor.
 	/*!
	 Use to set the starting time of the Parc
	 \param abstime: Starting time of the Parc in second
	 */

	~Scenario();	///< Normal Destructor.
	//@}


//********* Get Method *********//
	/*!
	 \name Get Function
	 */
	//@{
	cSecond				GetAbsoluteTime()	{ return fAbsoluteTime; }	///< Return the Absolute Clock
	map<cSecond, int>		GetTimeStep()		{ return fTimeStep; }		///< Return the Time Step Vector
	vector<Reactor*>		GetReactor()		{ return fReactor; }		///< Return the Reactor Vector
	vector<Storage*>		GetStorage()		{ return fStorage; }		///< Return the Storage Vector
	vector<Pool*>			GetPool()		{ return fPool; }		///< Return the Pool Vector
	vector<FabricationPlant*>	GetFabricationPlant()	{ return fFabricationPlant; }	///< Return the FabricationPlant Vector
	DecayDataBank*			GetDecayDataBase() 	{ return fDecayDataBase; }	//!< Return the Pointer to the Decay DataBase

	cSecond				GetPrintSet()		{ return fPrintStep; }		///< Return the Print Step Periodicity
	bool				GetStockManagement()	{ return fStockManagement; }	///< Return the StockManagement method (True or False)
	string				GetOutputFileName()	{ return fOutputFileName; }	///< Return the Output File name
	string				GetOutputTreeName()	{ return fOutputTreeName; }	///< Return the Output ROOT TTree name

	IsotopicVector			GetWaste()		{ return fWaste;}		///< Return the waste IsotopicVcetor

	//@}




//********* Set Method *********//
	/*!
	 \name Set Function
	 */
	//@{

	//{
	/// Set the Printing Step periodicity
	/*!
	 Use to set the periodicity of the output
	 \param timestep: periodicity of outpout in second
	 */
	void	SetTimeStep(cSecond timestep) 				{ fPrintStep = timestep; }
	//}

	//{
	/// Set the StockManagement method
	/*!
	 Use to define the stock managment method : true all fuel are stored individualy and false all fuel are mixed in a stock, and one can separate each isotope as needed
	 \param val: true or false depending on the stock management method used
	 */
	void	SetStockManagement(bool val)				{ fStockManagement = val; }
	//}

	//{
	/// Set the Decay DataBank
	/*!
	 Use to define Decay DataBank to be used
	 \param decaydatabase: a DecayDataBank which should contain the evolution of each nuclei of the chart
	 */
	void	SetDecayDataBase(DecayDataBank* decaydatabase) { fDecayDataBase = decaydatabase; }
	//}

	//{
	/// Set the Output File Name
	/*!
	 Use to define name of the output file
	 \param name: a string which correspond to the output file name
	 */
	void	SetOutputFileName(string name)	{ fOutputFileName = name; }
	//}


	//{
	/// Set the Output TTree Name
	/*!
	 Use to define name of the output ROOT TTree
	 \param name: a string which correspond to the output ROOT TTree name
	 */
	void	SetOutputTreeName(string name)	{ fOutputTreeName = name; }
	//}
	//@}

	void SetLogTimeStep(bool val = true)	{fLogTimeStep = true;}

//********* Add Method *********//
	/*!
	 \name Adding Facilities
	 */
	//@{

	void	AddPool(Pool* Pool);						///< Add a Pool to the Park
	void	AddReactor(Reactor* reactor);					///< Add a Reactor to the Park 
	void 	AddStorage(Storage* storage);					///< Add a Storage to the Park
	void 	AddFabricationPlant(FabricationPlant* fabricationplant);	///< Add a Storage to the Park
	void	AddSeparationPlant(SeparationPlant* separationplant);

	void	Add(Pool* Pool)					{AddPool(Pool);}	///< Add a Pool to the Park
	void	Add(Reactor* reactor)				{AddReactor(reactor);}	///< Add a Reactor to the Park
	void 	Add(Storage* storage)				{AddStorage(storage);}	///< Add a Storage to the Park
	void 	Add(FabricationPlant* fabricationplant)		{AddFabricationPlant(fabricationplant);}///< Add a Storage to the Park
	void 	Add(SeparationPlant* separationplant)		{AddSeparationPlant(separationplant);}///< Add a Storage to the Park

	//@}

	
 	
//********* Evolution Method *********//
	/*!
	 \name Evolution Method
	 */
	//@{

	void	BuildTimeVector(cSecond t);		///< Build the Time Evolution Vector where :
							/// \li 1 printing,
							/// \li 2 reactor Studown
							/// \li 4 start/End of reactor cycle,
							/// \li 8 end of Cooling,
							/// \li 16 fuel Fabrication

	void	Evolution(cSecond t);			///< Perform the Evolution
	void	BackEndEvolution();			///< Perform BackEnd Evolution
	void	PoolEvolution();			///< Perform Pool Evolution
	void	PoolDump();

	void	ReactorEvolution();			///< Perform the Reactor Evolution
	void	FabricationPlantEvolution();		///< Perform the FabricationPlant Evolution
	void	StorageEvolution();			///< Perform the Storage Evolution

	//@}



//-------- IsotopicVector --------//

	/*!
	 \name  IsotopicVector Sum
	 */
	//@{


	IsotopicVector		GetOutIncome() const		{ return fOutIncome; }		//!< Return the OutIncome Providings IsotopicVector
	void AddOutIncome(ZAI zai, double quantity)		{ AddOutIncome(zai*quantity); }	//!< Add a ZAI*quantity to OutIncomeIncome
	void AddOutIncome(IsotopicVector isotopicvector)	{ fOutIncome.Add(isotopicvector); }	//!< Add a isotopicVector to OutIncomeIncome
	void AddWaste(ZAI zai, double quantity)		{ AddWaste(zai*quantity); }	//!< Add a ZAI*quantity to Waste
	void AddWaste(IsotopicVector isotopicvector)	{ fWaste.Add(isotopicvector); }	//!< Add a isotopicVector to Waste
	void AddToPower(double power)			{ fParcPower += power;}		//!< Add power to the installed power in the Parc

	//@}



//********* In/Out related Method *********//

	/*!
	 \name  In/Out Method
	 */
	//@{


	void	ProgressPrintout(cSecond t);		//!< Update the prompt output to the time t
	
	void	Print();				//!< Print some information about the Parc
	void	Write();				//!< Write information in a file
	void	UpdateParc();				//!< Update the Global IsotopicVector
	void	OpenOutputTree();			//!< Open and define the Ouput ROOT TTree
	void	CloseOutputTree();			//!< Close and delete the Ouput ROOT TTree
	void	OutAttach();				//!< Attach the Branch to the Ouput ROOT TTree
	void	ResetQuantity();			//!< Reset the values of the GLobal IsotopicVector
	
	//@}

	
	
protected :
	bool			fNewTtree;		//!< Tru if we want to define a new TTree in the output File
	bool			fStockManagement;	///< True if real StockManagement false unstead (Default = true)
	bool			fLogTimeStep;
	
	cSecond			fPrintStep;		///< Time interval between two output update
	cSecond			fAbsoluteTime;		///< Absolute Clock
	cSecond			fStartingTime;		///< Starting Time
	map<cSecond, int>	fTimeStep;		///< Time Step Vector for the evolution :
							/// \li 1 printing,
							/// \li 2 reactor Studown
							/// \li 4 start/End of reactor cycle,
							/// \li 8 end of Cooling,
							/// \li 16 fuel Fabrication

	
	vector<Storage*>		fStorage;		///< Vector of Storages
	vector<Pool*>			fPool;			///< Vector of Pool
	vector<Reactor*>		fReactor;		///< Vector of Reactor
	vector<FabricationPlant*>	fFabricationPlant;	///< Vector of FabricationPlant
	vector<SeparationPlant*>	fSeparationPlant;	///< Vector of FabricationPlant
	DecayDataBank*			fDecayDataBase;		//!< Pointer to the Decay DataBase

	
	TFile*		fOutFile;		///< Pointer to the Root Output File
	string		fOutputFileName;	//! Name of the Output File
	TTree*		fOutT;			///< Pointer to the Root Output TTr3ee
	string		fOutputTreeName;	//! Name of the Output TTree
	string		fOutputLogName;		///< Name of the Ouput log File
	
	IsotopicVector	fWaste;			///< Waste IV
	IsotopicVector	fTotalStorage;		///< Sum of all IV in Storage IV
	IsotopicVector	fOutIncome;			///< OutIncomeIncome IV
	IsotopicVector	fTotalCooling;		///< Sum of all IV in Cooling IV
	IsotopicVector	fFuelFabrication;	///< Sum of all IV in Fabrication IV
	IsotopicVector	fTotalInReactor;	///< Sum of all IV in Reactor IV


	IsotopicVector	fIVInCycleTotal;	///< Sum of all IV in the cycle (without Waste) IV
	IsotopicVector	fIVTotal;		///< Sum of all IV in the parc (including Waste) IV
	double		fParcPower;		///< Sum of the Power of all reactor in the parc

};


#endif
