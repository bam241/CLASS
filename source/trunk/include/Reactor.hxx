#ifndef __Reactor_HXX__
#define __Reactor_HXX__

/*!
 \file
 \brief Header file for reactor classes. 
  Define a reactor.
 
 
 @author BaM
 @version 0.
 */

#include "TObject.h"
#include <string>
#include <map>
#include "IsotopicVector.hxx"
#include "EvolutiveProduct.hxx"


using namespace std;
typedef long long int cSecond;


class CLASS;
class TreatmentFactory;
class EvolutiveProduct;
template <class T> 
class EvolutionDataBase;
class FabricationPlant;
class Storage;
class LogFile;

class Reactor : public TObject
{
public :
	///< Normal Constructor.
	Reactor();
	///< Advbanced Constructor.
	Reactor(EvolutionDataBase<IsotopicVector>* 	fueltypeDB,
		FabricationPlant* fabricationplant, TreatmentFactory* treatmentfactory,
		double creationtime , double lifetime);						//!
	Reactor(EvolutionDataBase<IsotopicVector>* 	fueltypeDB,
		FabricationPlant* fabricationplant, TreatmentFactory* treatmentfactory,
		double creationtime , double lifetime, double cycletime,
		double HMMass, double BurnUp);						//!<
	Reactor(double Power, EvolutionDataBase<IsotopicVector>* 	fueltypeDB,
		FabricationPlant* fabricationplant, TreatmentFactory* treatmentfactory,
		double creationtime , double lifetime,
		double HMMass, double BurnUp);						//!<

	
	Reactor(EvolutiveProduct evolutivedb, TreatmentFactory* treatmentfactory,	//!<
		double creationtime, double lifetime, double cycletime);		//!<


	///< Normal Destructor
	~Reactor();
	

//********* Get Method *********//
	int 			GetId()		const		{ return fId; }			//!< Return the Reactor Parc'Is
	IsotopicVector 		GetIVReactor()	const		{ return fIVReactor; } 		//!< Return the IV contain in the Reactor
	IsotopicVector		GetIVBeginCycle() const		{ return fIVBeginCycle; }	//!< Return the Starting Cycle IV (Note : IVBegin != IVIn, only if using charging plan)
	IsotopicVector		GetIVOutCycle()	const		{ return fIVOutCycle; }		//!< Return the Out Cycle IV 
	IsotopicVector		GetIVInCycle()	const		{ return fIVInCycle; }		//!< Return the In Cycle IV (Note : IVIn != IVBegin, only if using charging plan)
	cSecond 		GetCycleTime()	const		{ return fCycleTime; } 		//!< Return the cycle time of the Reactor
	cSecond 		GetCreationTime() const		{ return fCreationTime; }	//!< Return the creation time of the Reactor
	cSecond 		GetLifeTime()	const		{ return fLifeTime; }		//!< Return the creation time of the Reactor

	EvolutiveProduct	GetEvolutionDB()		const	{ return fEvolutionDB; }		//!< Return the Evolution database of the Fuel
	TreatmentFactory*	GetAssociedTreatmentFactory()	const	{ return fAssociedTreatmentFactory; }	//!< Return the pointer to Associeted TF
	LogFile*		GetLog()			const	{ return fLog; }			//!< Return the Pointer to Log

	bool 			IsFuelFixed()				{ return fFixedFuel; }			//!< True if using fixed Fuel, False otherwise
	FabricationPlant*	GetFabricationPlant()		const	{ return fFabricationPlant; }		//!< Return the Pointer to the FabricationPlant
	EvolutionDataBase<IsotopicVector>* GetFuelType()	const	{ return fFuelTypeDB; }			//!< Return the Fuel Type DB of the reactor
	double			GetHeavyMetalMass()		const	{ return fHeavyMetalMass; }		//!< Return the HeavyMetal Mass in the Core at the begining of the cycle
	double			GetBurnUp()			const	{ return fBurnUp; }			//!< Return the Burn Up of the Fuel at the end of the cycle


//********* Set Method *********//
	void SetId(int id)					{ fId = id; }					//!< Set The Reactor Parc'Id
	void SetParc(CLASS* parc)				{ fParc = parc; }				//!< Set the Pointer to the Parc
	void SetStorage(Storage* storage)			{ fStorage = storage; fIsStorage = true;}	//!< Set the Pointer to the Storage
	void SetLog(LogFile* LOG)				{ fLog = LOG; }				//!< Set the Pointer to the Log
	void SetIVReactor(IsotopicVector isotopicvector)	{ fIVReactor = isotopicvector; }	//!< Set the IV inside the Reactor Core
	void SetIVBeginCycle(IsotopicVector isotopicvector)	{ fIVBeginCycle = isotopicvector; }	//!< Set the IV at the Beginging of the Reactor Cycle
	void SetIVOutCycle(IsotopicVector isotopicvector)	{ fIVOutCycle = isotopicvector; }	//!< Set the IV Going Out at the End of the Cycle
	void SetIVInCycle(IsotopicVector isotopicvector)	{ fIVInCycle = isotopicvector; }	//!< Set the IV Coming In at the Beginning of the Cycle 
	void SetCycleTime(double cycletime);								//!< Set the Power time (Cycle of the loading Plan)
	void SetPower(double Power);									//!< Set the Power
	void SetHMMass(double Mass)		{fHeavyMetalMass = Mass;}				//!< Set the Mass
	void SetBurnUp(double BU)		{fBurnUp = BU;}						//!< Set the Mass

	void SetEvolutionDB(EvolutiveProduct evolutionDB);						//!< Set the Pointer to the DB Evolution of the Reactor
	
//********* Modification Method *********//
	void Evolution(cSecond t);									//!< Performe the Evolution until the Time t
	void Dump();											//!< Write Modification (IV In/Out, filling the TF...)
	void SetNewFuel(EvolutiveProduct ivdb);								//!< Change the Evolutive DB of the Reactor
//********* Other Method *********//
	
	
protected :
	int		fId;			//!< Identity of the Reactor inside the Parc
	cSecond		fInternalTime;		///< Internal Clock
	cSecond		fInCycleTime;		///< Time spend since the beginning of the last Cycle
	bool		fIsStarted;		///< True if Running, False Otherwise
	bool		fShutDown;		///< True if ShutDown
	bool		fEndOfCycle;		///< True if Reaching the End of a Reactor Cycle
	
	bool		fFixedFuel;
	bool		fIsStorage;
	
//********* Internal Parameter *********//
 	LogFile*		fLog;				//!< Pointer to the Log
	CLASS*			fParc;				//!< Pointer to the main Parc
	TreatmentFactory*	fAssociedTreatmentFactory;	//!< Pointer to the TF which collect the spend fuel
	Storage*		fStorage;			//!< Pointer to the Stock
								//!<
	EvolutiveProduct	fEvolutionDB;			//!< Pointer to the Evolution DataBase
	EvolutionDataBase<IsotopicVector>* 	fFuelTypeDB;	//! Pointer to a Fuel Type Database
	
	cSecond		fCreationTime;		///< CLASS Universal Time of Creation
	cSecond		fLifeTime;		///< LifeTime Of the Reactor
	cSecond		fCycleTime;		///< Cycle Time
	double 		fPower;			///< Power
	
	IsotopicVector	fIVReactor;		///< Fuel evoluated IV in the reactor
	IsotopicVector	fIVBeginCycle;		///< Fuel IV at the Beginning of a Cycle
	IsotopicVector	fIVInCycle;		///< IVBegin add at the Beginning of the Cycle
	IsotopicVector	fIVOutCycle;		///< IV wich get out at the End of a Cycle
	

//********* Unfixed Fuel Parameter *********//


	FabricationPlant*	fFabricationPlant;		//!< Poitner to the FabricationPlant
	double			fHeavyMetalMass;		///< In tons
	double			fBurnUp;			///< In GWd/tHM

 	ClassDef(Reactor,1);
 };


#endif
