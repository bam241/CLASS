#ifndef _CLASS_HXX_
#define _CLASS_HXX_

/*!
 \file
 \brief Header file for CLASS classes. 
 Define a CLASS Parc.
 The aim of thes Class is to manage the parc, store reactor, TreatmentFactory, process the evolution, build Isotopic vector

 
 @author BaM
 @version 0.
 */
#include "IsotopicVector.hxx"

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <map>
#include <iostream>


using namespace std;

template <class T> class EvolutionDataBase;
class FabricationPlant;
class Reactor;
class TreatmentFactory;
class LogFile;
class Storage;

class CLASS 
{
public :
	///< Normal Constructor.
	CLASS();
	CLASS(double abstime);
	CLASS(string name);
	CLASS(string name, double abstime);

	///< Normal Destructor.
	~CLASS();


//********* Get Method *********//
	double				GetAbsoluteTime()	{ return fAbsoluteTime; }	///< Return the Absolute Clock
	map<double, int>		GetTimeStep()		{ return fTimeStep; }		///< Return the Time Step Vector
	vector<TreatmentFactory*>	GetTreatmentFactory()	{ return fTreatmentFactory; }	///< Return the TF Vector
	vector<Reactor*>		GetReactor()		{ return fReactor; }		///< Return the Reactor Vector
	vector<Storage*>		GetStorage()		{ return fStorage; }		///< Return the Reactor Vector
	EvolutionDataBase<ZAI>*		GetDecayDataBase() 	{ return fDecayDataBase; }	//!< Return the Pointer to the Decay DataBase

	double				GetPrintSet()		{ return fPrintStep; }		///< Return the Print Step Periodicity
	bool				GetStockManagement()	{ return fStockManagement; }
//	bool				GetOpenCycle()		{ return fOpenCycle; }
	string				GetOutputName()		{ return fOutputName; }
	LogFile*			GetLog()		{ return fLog; }


//********* Set Method *********//
	void	SetTimeStep(double timestep) 				{ fPrintStep = timestep; }		///< Set the Printing Step periodicity
	void	SetStockManagement(bool val)				{ fStockManagement = val; }
	void	SetDecayDataBase(EvolutionDataBase<ZAI>* decaydatabase) { fDecayDataBase = decaydatabase; }	//!< Set the Pointer to the Decay DataBase
//	void	SetOpenCycle(bool val)		{ fOpenCycle = val; }
	
	void	SetBuildingMethod(int val)	{ fStockManagement = true; fBuildingMethod = val; }
	void	SetOutpurName(string name)	{ fOutputName =name; }


//********* Add Method *********//
	void	AddTreatmentFactory(TreatmentFactory* treatmentfactory);	///< Add A TF to the Park
	void	AddReactor(Reactor* reactor);					///< Add a Reactor to the Park 
	void 	AddStorage(Storage* storage);					///< Add a Storage to the Park
	void 	AddFabricationPlant(FabricationPlant* fabricationplant);	///< Add a Storage to the Park

//********* IV Creation Method *********//
//	IsotopicVector	BuildIsotopicVector(IsotopicVector isotopicvector );	///< Build The needed Isotopic Vector from the Storage
//	IsotopicVector	MonteCarloBuild(IsotopicVector isotopicvector );	///< Build the needed IV with the MonteCarlo Method
//	IsotopicVector	MinimizationBuild(IsotopicVector isotopicvector );	///< Build the needed IV with the Minimization Method
//	void	UpdateParcStorage();
//	void	DumpParcStorage();
	
	
	
	
	
 	
//********* Evolution Method *********//
	void	BuildTimeVector(double t);				///< Build the Time Evolution Vector
	void	Evolution(double t);					///< Do the Evolution
	void	TreatmentEvolution();					///< Do TF Evolution
	void	ReactorEvolution();					///< Do the Reactor Evolution
	void	FabricationPlantEvolution();					///< Do the FabricationPlant Evolution
	void	StorageEvolution();					///< Do the Storage Evolution


	void	AddTotalStorage(IsotopicVector IV)		{ fIVTotal +=IV; fIVInCycleTotal +=IV; fTotalStorage +=IV; }
	void	RemoveTotalStorage(IsotopicVector IV)		{ fIVTotal -=IV; fIVInCycleTotal -=IV; fTotalStorage -=IV; }
	void	AddTotalCooling(IsotopicVector IV)		{ fIVTotal +=IV; fIVInCycleTotal +=IV; fTotalCooling +=IV; }
	void	RemoveTotalCooling(IsotopicVector IV)		{ fIVTotal -=IV; fIVInCycleTotal -=IV; fTotalCooling -=IV; }
	void	AddTotalInReactor(IsotopicVector IV)		{ fIVTotal +=IV; fIVInCycleTotal +=IV; fTotalInReactor +=IV; }
	void	RemoveTotalInReactor(IsotopicVector IV)		{ fIVTotal -=IV; fIVInCycleTotal -=IV; fTotalInReactor -=IV; }
	void	AddIVInCycleTotal(IsotopicVector IV)		{ fIVTotal +=IV; fIVInCycleTotal +=IV; }
	void	AddIVTotal(IsotopicVector IV)			{ fIVTotal +=IV; }

//-------- GodIncome --------//

	IsotopicVector		GetGodIncome() const		{ return fGodIncome; }	//!< Return the God Providings IsotopicVector
	void AddGodIncome(ZAI zai, double quantity)		{ AddGodIncome(zai*quantity); }	//!< Add a ZAI*quantity to GodIncome
	void AddGodIncome(IsotopicVector isotopicvector)	{ fGodIncome.Add(isotopicvector); }	//!< Add a isotopicVector to GodIncome
	void AddWaste(ZAI zai, double quantity)			{ AddWaste(zai*quantity); }	//!< Add a ZAI*quantity to Waste
	void AddWaste(IsotopicVector isotopicvector)		{ fWaste.Add(isotopicvector); }	//!< Add a isotopicVector to Waste



//********* In/Out related Method *********//
	void	ProgressPrintout( double t);
	
	void	Print();
	void	Write();
	void	UpdateParc();
	void	OpenOutputTree();
	void	CloseOutputTree();
	void	OutAttach();
	void	ResetQuantity();
	
	
	
	
protected :
	double			fPrintStep;
	double			fAbsoluteTime;		///< Absolute Clock
	double 			fStartingTime;
	map<double, int>	fTimeStep;		///< Time Step Vector for the evolution : 
							///< 1 Printing, 
							///< 2 Reactor Studown
							///< 4 Start/End of reactor cycle,
							///< 8 End of Cooling,
							///< 16 Fuel Fabrication 

	vector<Storage*>		fStorage;			///< Vector of Storages
	vector<TreatmentFactory*>	fTreatmentFactory;		///< Vector of Treament Factory
	vector<Reactor*>		fReactor;			///< Vector of Reactor
	vector<FabricationPlant*>	fFabricationPlant;		///< Vector of FabricationPlant

	EvolutionDataBase<ZAI>*		fDecayDataBase;			//!< Pointer to the Decay DataBase

	bool				fStockManagement;		///< ...
//	bool 				fOpenCycle;
	int				fBuildingMethod;
	
	
	LogFile*			fLog;
	
	string				fOutputName;
	TFile*				fOutTree;
	TTree*				fOutT;

	IsotopicVector			fWaste;
	IsotopicVector			fTotalStorage;
	IsotopicVector			fGodIncome;
	IsotopicVector			fTotalCooling;
	IsotopicVector			fFuelFabrication;
	
	IsotopicVector			fTotalInReactor;
	
	IsotopicVector			fIVInCycleTotal;
	IsotopicVector			fIVTotal;
	
	map< pair<int,int>, IsotopicVector >	fParcStorage;


};


#endif
