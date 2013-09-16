#ifndef __Reactor_HXX__
#define __Reactor_HXX__

/*!
 \file
 \brief Header file for reactor classes. 
  Define a reactor.
 
 
 @author BaM
 @version 2.0
 */

#include <string>
#include <map>

#include "CLSSFacility.hxx"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"

using namespace std;
typedef long long int cSecond;


class CLASS;
class Pool;
class EvolutionData;
template <class T> 
class DataBank;
class FabricationPlant;
class Storage;
class LogFile;

class Reactor : public CLSSFacility
{
public :
	///< Normal Constructor.
	Reactor();
	Reactor(LogFile* log);
	///< Advbanced Constructor.
	Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB,
		FabricationPlant* fabricationplant, Pool* Pool,
		double creationtime , double lifetime);				//!<
	
	Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB,
		FabricationPlant* fabricationplant, Pool* Pool,
		double creationtime , double lifetime, double cycletime,
		double HMMass, double BurnUp);					//!<
	
	Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB,
		FabricationPlant* fabricationplant, Pool* Pool,
		double creationtime , double lifetime,
		double Power, double HMMass, double BurnUp, double ChargeFactor);	//!<

	Reactor(LogFile* log, EvolutionData evolutivedb, Pool* Pool,
		double creationtime, double lifetime,
		double power, double HMMass, double BurnUp, double ChargeFactor);

	///< Normal Destructor
	~Reactor();
	

//********* Get Method *********//
	IsotopicVector 	GetIVReactor()		const	{ return GetInsideIV(); } 	//!< Return the IV contain in the Reactor
	IsotopicVector	GetIVBeginCycle()	const	{ return fIVBeginCycle; }	//!< Return the Starting Cycle IV
											//!< (Note : IVBegin != IVIn, only if using charging plan)
	
	IsotopicVector	GetIVOutCycle()		const	{ return fIVOutCycle; }		//!< Return the Out Cycle IV
	IsotopicVector	GetIVInCycle()		const	{ return fIVInCycle; }		//!< Return the In Cycle IV
											//!< (Note : IVIn != IVBegin, only if using charging plan)
	

	EvolutionData	GetEvolutionDB()		const	{ return fEvolutionDB; }	//!< Return the Evolution database of the Fuel
	DataBank<IsotopicVector>*	GetFuelType()	const	{ return fFuelTypeDB; }		//!< Return the Fuel Type DB of the reactor

	Pool*			GetAssociedPool()	const	{ return fAssociedPool; }	//!< Return the pointer to Associeted TF
	FabricationPlant*	GetFabricationPlant()	const	{ return fFabricationPlant; }	//!< Return the Pointer to the FabricationPlant

	bool	IsFuelFixed()		const	{ return fFixedFuel; }		//!< True if using fixed Fuel, False otherwise
	double	GetHeavyMetalMass()	const	{ return fHeavyMetalMass; }	//!< Return the HeavyMetal Mass in the Core at the begining of the cycle
	double	GetBurnUp()		const	{ return fBurnUp; }		//!< Return the Burn Up of the Fuel at the end of the cycle
	double	GetPower()		const	{ return fPower; } 		//!< Return the cycle time of the Reactor


//********* Set Method *********//
	void SetStorage(Storage* storage)	{ fStorage = storage; fIsStorage = true;}	//!< Set the Pointer to the Storage

	void SetIVReactor(IsotopicVector isotopicvector)	{ fInsideIV = isotopicvector; }	//!< Set the IV inside the Reactor Core
	void SetIVBeginCycle(IsotopicVector isotopicvector)	{ fIVBeginCycle = isotopicvector; }	//!< Set the IV at the Beginging of the Reactor Cycle
	void SetIVOutCycle(IsotopicVector isotopicvector)	{ fIVOutCycle = isotopicvector; }	//!< Set the IV Going Out at the End of the Cycle
	void SetIVInCycle(IsotopicVector isotopicvector)	{ fIVInCycle = isotopicvector; }	//!< Set the IV Coming In at the Beginning of the Cycle
	
	void SetCycleTime(double cycletime);					//!< Set the Power time (Cycle of the loading Plan)
	
	void SetEvolutionDB(EvolutionData evolutionDB);				//!< Set the Pointer to the DB Evolution of the Reactor
	void SetPower(double Power);						//!< Set the Power
	void SetHMMass(double Mass)		{fHeavyMetalMass = Mass;}	//!< Set the HeavyMetal Mass in the Core at the begining of the cycle
	void SetBurnUp(double BU)		{fBurnUp = BU;}			//!< Set the the Burn Up of the Fuel at the end of the cycle

//********* Modification Method *********//
	void Evolution(cSecond t);						//!< Performe the Evolution until the Time t
	void Dump();								//!< Write Modification (IV In/Out, filling the TF...)
	void SetNewFuel(EvolutionData ivdb);					//!< Change the Evolutive DB of the Reactor
	
	
//********* Other Method *********//
	
	
protected :
	
	bool		fFixedFuel;
	bool		fIsStorage;
	
//********* Internal Parameter *********//
	Pool*		fAssociedPool;		//!< Pointer to the TF which collect the spend fuel
	Storage*	fStorage;		//!< Pointer to the Stock
						
	EvolutionData	fEvolutionDB;			//!< Pointer to the Evolution DataBase
	DataBank<IsotopicVector>* 	fFuelTypeDB;	//! Pointer to a Fuel Type Database
	
	double 		fPower;			///< Power (in Watt)
	
	IsotopicVector	fIVBeginCycle;		///< Fuel IV at the Beginning of a Cycle
	IsotopicVector	fIVInCycle;		///< IVBegin add at the Beginning of the Cycle
	IsotopicVector	fIVOutCycle;		///< IV wich get out at the End of a Cycle
	map<ZAI, double> fZAImass;


//********* Unfixed Fuel Parameter *********//


	FabricationPlant*	fFabricationPlant;		//!< Poitner to the FabricationPlant
	double			fHeavyMetalMass;		///< In tons
	double			fBurnUp;			///< In GWd/tHM

 	ClassDef(Reactor,2);
 };


#endif
