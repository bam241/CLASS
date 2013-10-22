#ifndef __FabricationPlant_HXX__
#define __FabricationPlant_HXX__

/*!
 \file 
 \brief Header file for FabricationPlant class.

 The aim of the Class is to manage evolution of FabricationPlant

 
 @author BaM, Marc
 @version 2.0
 */



#include <vector>
#include <map>

#include "CLSSFacility.hxx"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"
#include "CLASS.hxx"
#include "Storage.hxx"
#include "Reactor.hxx"
#include "DataBank.hxx"
#include "LogFile.hxx"
#include "ZAI.hxx"

using namespace std;
typedef long long int cSecond;





class FabricationPlant : public CLSSFacility
{

public :
	///< Normal constructor
	FabricationPlant();
	FabricationPlant(LogFile* log);
	
	FabricationPlant(LogFile* log, Storage* storage, Storage* reusable, double fabricationtime = 365.25*24*3600*2);
	///< Normal Destructor.
	~FabricationPlant();




//********* Set Method *********//
	void	SetUpdateReferenceDBatEachStep(bool val){ fUpdateReferenceDBatEachStep = val;}
	void	SetStorage(Storage* storage)		{ fStorage = storage; }			//!< Set the Pointer to the Storage
	void	SetDecayDataBase(DataBank<ZAI>* ddb)	{ fDecayDataBase = ddb; }	//!< Set the pointer to the Decay DataBase
	
	void	SetChronologicalTimePriority(bool bval)	{ fChronologicalTimePriority = bval;}	//!< Set the chronological priority (true for chronological, false unstead)
	
	void	SetSubstitutionFuel(EvolutionData fuel);
	
	void	AddReactor(int reactorid, double creationtime)
			{ fReactorNextStep.insert( pair<int,cSecond> (reactorid, (cSecond)creationtime-GetCycleTime() ) ); }	//!< Add a new reactor

//********* Get Method *********//

	Storage*	GetStorage()		{ return fStorage; }		//!< Return the Pointer to the Storage	
	
	map<int, IsotopicVector >	GetReactorFuturIncome() const
						{ return fReactorFuturIV;}	//!< Return the List of the Futur Fuel IV
	DataBank<ZAI>* 	GeDecayDataBase() const
						{ return fDecayDataBase; }	//!< Return the pointer to the DecayDB

	IsotopicVector GetFullFabrication();

	EvolutionData GetReactorEvolutionDB(int ReactorId);			//!< Return the EvolutionData of Reactor ReactorId

	

//---------- FabricationPlant ----------//

	void	AddValorisableIV(ZAI zai, double factor);						///< Add Valorisable Element
	void	Evolution(cSecond t);									//!< Perform the Evolution
	void	BuildFuelForReactor(int ReactorId);							//!< Build a Fuel for the reactor ReactorId
	void	RecycleStock(double fraction);								//!< Take a franction of the current stock
	IsotopicVector	GetStockToRecycle();								//!< Get the next stock to recycle
	void 	DumpStock();										//!< Update the storage
	EvolutionData	BuildEvolutiveDB(int ReactorId, IsotopicVector isotopicvector);		//!< Build the Evolution Database for the Reactir ReactorId Fuel
	void	TakeReactorFuel(int ReactorId) ;								//!< Remove the Fuel of reactor ReactorId
	


//********* Other Method *********//






protected :
	bool		fUpdateReferenceDBatEachStep;

//********* Internal Parameter *********//
	map<ZAI ,double>	fValorisableIV;		///< The Valorisable Table
	map<int, cSecond >	fReactorNextStep;	///< Next Time Step to Build a New Fuel

	map< int,EvolutionData >	fReactorFuturDB; ///< List of the Futur EvolutionData use in the reactor
	map< int,IsotopicVector >	fReactorFuturIV; ///< List of the Futur Fuel Isotopic Vector used in the reactor

	DataBank<ZAI>*		fDecayDataBase;	//!< Pointer to the Decay DataBase


	Storage*	fStorage;			//!< Pointer to the Storage to recycle
	Storage*	fReUsable;			//!< Pointer to the Storage using for recycling unused Product

	vector< pair<int, double> >	fFractionToTake;	//!< The Temporary Storage IsotopicVector

		//	double		fFabricationTime;		///< Fabrication Duration Time
	bool		fChronologicalTimePriority;	//!< Set the Chronological Priotity (for the Stock Management) or the anti-chronological one

	bool		fSubstitutionFuel;
	EvolutionData	fSubstitutionEvolutionData;
	

//********* Private Method *********//
	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time
	void	FabricationPlantEvolution(cSecond t);				//!< Deal the FabricationPlant Evolution
	pair<IsotopicVector, IsotopicVector> Separation(IsotopicVector isotopicvector);	//!< Make the Separation 
						//!< return IV[0] -> To Stock / IV[1] -> To Waste

	ClassDef(FabricationPlant,2);

};

#endif
