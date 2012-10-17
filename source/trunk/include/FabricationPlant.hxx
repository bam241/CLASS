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
	
	void	SetDecayDataBase(EvolutionDataBase<ZAI>* ddb)		{ fDecayDataBase = ddb; }		//!< Set the pointer to the Decay DataBase
	void	AddReactor(int reactorid, double creationtime)		{ fReactorNextStep.insert( pair<int,double> (reactorid, creationtime-fFabricationTime) ); }	//!< ...

//********* Get Method *********//

	LogFile*	GetLog()		{ return fLog; }		//!< Return the Pointer to the Log
	CLASS*		GetParc()		{ return fParc; }		//!< Return the Pointer to the Parc
	Storage*	GetStorage()		{ return fStorage; }		//!< Return the Pointer to the Storage

	double	GetInternalTime()	const	{ return fInternalTime; }	//!< Return Creation Time
	double	GetFabricationTime()	const	{ return fFabricationTime; }	//!< Return the Fabrication Time
	
	
	map<int, IsotopicVector >	GetReactorFuturIncome()	const	{return fReactorFuturIncome;}	//!< ...
	EvolutionDataBase<ZAI>* 	GeDecayDataBase()	const	{ return fDecayDataBase; }	//!< Return the pointer to the Decay DataBase

	EvolutiveProduct		GetReactorEvolutionDB(int ReactorId);				//!<



//---------- FabricationPlant ----------//

	void		AddValorisableIV(ZAI zai, double factor);						///< Add Valorisable Element
	void		Evolution(double t);									//!< ...
	void		BuildFuelForReactor(int ReactorId);							//!< ...
	void		BuildMOXFuelForReactor(int ReactorId, EvolutionDataBase<IsotopicVector>* FuelType);	//!< ...
	void		BuildADSFuelForReactor(int ReactorId, EvolutionDataBase<IsotopicVector>* FuelType);	//!< ...
	void		RecycleStock(double fraction);								//!< ...
	IsotopicVector	GetStockToRecycle();									//!< ...
	void 		DumpStock();										//!< ...
	EvolutiveProduct	BuildEvolutiveDB(int ReactorId, IsotopicVector isotopicvector);			//!< ...
	void TakeReactorFuel(int Id) ;										//!< ...
	


//********* Other Method *********//






protected :
	int		fId;			//!< Identity of the FabricationPlant inside the Parc
	double 		fInternalTime;		///< Internal Clock


//********* Internal Parameter *********//
	CLASS* 		fParc;				//!< Pointer to the Parc
	LogFile*	fLog;				//!< Pointer to the Log

	map<ZAI ,double>	fValorisableIV;		///< The Valorisable Table
	map<int, double >	fReactorNextStep;	//!< ...

	map< int,EvolutiveProduct >	fReactorIncome; 	//!<
	map< int,IsotopicVector >	fReactorFuturIncome; 	//!<

	EvolutionDataBase<ZAI>*		fDecayDataBase;		//!< Pointer to the Decay DataBase


	Storage*	fStorage;				//!< Pointer to the Storage to recycle
	Storage*	fReUsable;				//!< 

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
