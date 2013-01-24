#ifndef __FabricationPlant_HXX__
#define __FabricationPlant_HXX__

/*!
 \file 
 \brief Header file for FabricationPlant class.

 The aim of the Class is to manage evolution of FabricationPlant

 
 @author BaM
 @version 0.
 */



#include "IsotopicVector.hxx"
#include "EvolutiveProduct.hxx"

#include <vector>
#include <map>
#include "TObject.h"
using namespace std;

class Storage;
class CLASS;
class LogFile;
class ZAI;
class IsotopicVector;		
template <class T> 
class EvolutionDataBase;





class FabricationPlant : public TObject
{

public :
	///< Normal constructor
	FabricationPlant();
	
	FabricationPlant(Storage* storage, Storage* reusable, double fabricationtime = 365.25*24*3600*2);
	///< Normal Destructor.
	~FabricationPlant();




//********* Set Method *********//
	void	SetId(int id)				{ fId = id; }				//!< Set The FB Parc'Id
	void	SetParc(CLASS* parc)			{ fParc = parc; }			//!< Set the Pointer to the Parc
	void	SetLog(LogFile* Log)			{ fLog = Log; }				//!< Set the Pointer to the Log
	void	SetStorage(Storage* storage)		{ fStorage = storage; }			//!< Set the Pointer to the Storage
	void	SetDecayDataBase(EvolutionDataBase<ZAI>* ddb)	{ fDecayDataBase = ddb; }	//!< Set the pointer to the Decay DataBase
	

	
	void	AddReactor(int reactorid, double creationtime)
			{ fReactorNextStep.insert( pair<int,double> (reactorid, creationtime-fFabricationTime) ); }	//!< Add a new reactor

//********* Get Method *********//

	LogFile*	GetLog()		{ return fLog; }		//!< Return the Pointer to the Log
	CLASS*		GetParc()		{ return fParc; }		//!< Return the Pointer to the Parc
	Storage*	GetStorage()		{ return fStorage; }		//!< Return the Pointer to the Storage

	double	GetInternalTime() const		{ return fInternalTime; }	//!< Return Creation Time
	double	GetFabricationTime() const	{ return fFabricationTime; }	//!< Return the Fabrication Time
	
	
	map<int, IsotopicVector >	GetReactorFuturIncome() const
						{ return fReactorFuturIV;}	//!< Return the List of the Futur Fuel IV
	EvolutionDataBase<ZAI>* 	GeDecayDataBase() const
						{ return fDecayDataBase; }	//!< Return the pointer to the DecayDB

	
	EvolutiveProduct GetReactorEvolutionDB(int ReactorId);			//!< Return the EvolutiveProduct of Reactor ReactorId

	

//---------- FabricationPlant ----------//

	void	AddValorisableIV(ZAI zai, double factor);						///< Add Valorisable Element
	void	Evolution(double t);									//!< Perform the Evolution
	void	BuildFuelForReactor(int ReactorId);							//!< Build a Fuel for the reactor ReactorId
	void	RecycleStock(double fraction);								//!< Take a franction of the current stock
	IsotopicVector	GetStockToRecycle();								//!< Get the next stock to recycle
	void 	DumpStock();										//!< Update the storage
	EvolutiveProduct	BuildEvolutiveDB(int ReactorId, IsotopicVector isotopicvector);		//!< Build the Evolution Database for the Reactir ReactorId Fuel
	void	TakeReactorFuel(int Id) ;								//!< Remove the Fuel of reactor ReactorId
	


//********* Other Method *********//






protected :
	int		fId;			//!< Identity of the FabricationPlant inside the Parc
	double 		fInternalTime;		///< Internal Clock


//********* Internal Parameter *********//
	CLASS* 		fParc;				//!< Pointer to the Parc
	LogFile*	fLog;				//!< Pointer to the Log

	map<ZAI ,double>	fValorisableIV;		///< The Valorisable Table
	map<int, double >	fReactorNextStep;	//!< Next Time Step to Build a New Fuel

	map< int,EvolutiveProduct >	fReactorFuturDB; //!< List of the Futur EvolutiveProduct use in the reactor
	map< int,IsotopicVector >	fReactorFuturIV; //!< List of the Futur Fuel Isotopic Vector used in the reactor

	EvolutionDataBase<ZAI>*		fDecayDataBase;	//!< Pointer to the Decay DataBase


	Storage*	fStorage;			//!< Pointer to the Storage to recycle
	Storage*	fReUsable;			//!< Pointer to the Storage using for recycling unused Product

	vector< pair<int, double> >	fFractionToTake;	//!< The Temporary Storage IsotopicVector

	double		fFabricationTime;		///< Fabrication Duration Time
	bool		fChronologicalTimePriority;	//!< Set the Chronological Priotity (for the Stock Management) or the anti-chronological one


//********* Private Method *********//
	IsotopicVector GetDecay(IsotopicVector isotopicvector, double t);	//!< Get IsotopicVector Decay at the t time
	void	FabricationPlantEvolution(double t);				//!< Deal the FabricationPlant Evolution
	pair<IsotopicVector, IsotopicVector> Separation(IsotopicVector isotopicvector);	//!< Make the Separation 
						//!< return IV[0] -> To Stock / IV[1] -> To Waste

	ClassDef(FabricationPlant,0);

};

#endif
