#ifndef _Reactor_
#define _Reactor_

/*!
 \file
 \brief Header file for reactor classes.
 */

#include <string>
#include <map>

#include "CLASSFacility.hxx"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"
#include "PhysicsModels.hxx"
#include "ReactorScheduler.hxx"

using namespace std;

class CLASSBackEnd;
class FabricationPlant;
class Storage;
class CLASSLogger;

//-----------------------------------------------------------------------------//
//! Defines the Reactor

/*!
 The aim of this class is to deal the evolution of the fuel inside a reactor.
 The fuel state in the reactor is describe in the IsotopicVector. Its evolution is contained in an EvolutionData
 
 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class Reactor : public CLASSFacility
{
	public :
	
	
	//********* Constructor/Destructor Method *********//
	
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	Reactor();		///< Normal Constructor.
	
	//{
	/// CLASSLogger Constructor.
	/*!
	 Use create an empty Reactor loading a CLASSLogger
	 \param log: used for the log.
	 */
	Reactor(CLASSLogger* log);
	//}
#ifndef __ROOTCLING__	
	//{
	/// Special Constructor for reprocessed fuel using cycletime and Burn-Up.
	/*!
	 Make a new reactor
	 \param log: used for the log,
	 \param backend: CLASSBackend which get the fuel after the cooling,
	 \param creationtime: creation time in [s],
	 \param lifetime: working time duration in [s],
	 \param Power: Thermal power of the reactor in [W],
	 \param HMMass: Mass of Heavy Metal in the Reactor in [t] of heavy metal,
	 \param CapacityFactor effective charge of the reactor, fraction between 0 & 1.
	 */
	Reactor(CLASSLogger* log, CLASSBackEnd* backend,
		cSecond creationtime , cSecond lifetime, double Power,
		double HMMass, double CapacityFactor = 1);
	//}
	
	//{
	/// Special Constructor for reprocessed fuel using cycletime and Burn-Up.
	/*!
	 Make a new reactor
	 \param log: used for the log.
	 \param backend: CLASSBackend which get the fuel after the cooling,
	 \param creationtime: creation time in [s],
	 \param lifetime: working time duration in [s],
	 \param cycletime: duration of a cycle in [s],
	 \param Power: Thermal power of the reactor in [W],
	 \param HMMass: Mass of Heavy Metal in the Reactor in [t] of heavy metal,
	 \param CapacityFactor effective charge of the reactor, fraction between 0 & 1.
	 */
	Reactor(CLASSLogger* log,
		FabricationPlant* fabricationplant, CLASSBackEnd* backend,
		cSecond creationtime , cSecond lifetime, double Power,
		double HMMass, double CapacityFactor = 1);
	//}
	
	//{
	/// Special Constructor for reprocessed fuel using cycletime and Burn-Up.
	/*!
	 Make a new reactor
	 \param log: used for the log,
	 \param fueltypeDB Databank describing the evolution of the fuel;
	 \param backend: CLASSBackend which get the fuel after the cooling,
	 \param creationtime: creation time in [s],
	 \param lifetime: working time duration in [s],
	 \param cycletime: duration of a cycle in [s],
	 \param HMMass: Mass of Heavy Metal in the Reactor in [t] of heavy metal,
	 \param BurnUp: Burnup reach by the fuel at the end of the cycle in [GWd/t].
	 */
	Reactor(CLASSLogger* log, PhysicsModels* fueltypeDB,
		FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		cSecond creationtime , cSecond lifetime, cSecond cycletime,
		double HMMass, double BurnUp);
	//}
	
	//{
	/// Special Constructor for reprocessed fuel using Power and Burn-Up.
	/*!
	 Make a new reactor
	 \param log: used for the log,
	 \param fueltypeDB Databank describing the evolution of the fuel,
	 \param backend: CLASSBackend which get the fuel after the cooling,
	 \param creationtime: creation time in [s],
	 \param lifetime: working time duration in [s],
	 \param Power: Thermal power of the reactor in [W],
	 \param HMMass: Mass of Heavy Metal in the Reactor in [t] of heavy metal,
	 \param BurnUp: Burnup reach by the fuel at the end of the cycle in [GWd/t],
	 \param CapacityFactor effective charge of the reactor, fraction between 0 & 1.
	 */
	Reactor(CLASSLogger* log, PhysicsModels* fueltypeDB,
		FabricationPlant* fabricationplant, CLASSBackEnd* backend,
		cSecond creationtime , cSecond lifetime,
		double Power, double HMMass, double BurnUp, double CapacityFactor);
	//}
	
	//{
	/// Special Constructor for fixed fuel using Power and Burn-Up.
	/*!
	 Make a new reactor
	 \param log: used for the log,
	 \param evolutivedb: EvolutionData describing the evolution of the fuel,
	 \param backend: CLASSBackend which get the fuel after the cooling,
	 \param creationtime: creation time in [s],
	 \param lifetime: working time duration in [s],
	 \param Power: Thermal power of the reactor in [W],
	 \param HMMass: Mass of Heavy Metal in the Reactor in [t] of heavy metal,
	 \param BurnUp: Burnup reach by the fuel at the end of the cycle in [GWd/t],
	 \param CapacityFactor: effective charge of the reactor, fraction between 0 & 1.
	 */
	Reactor(CLASSLogger* log, EvolutionData* evolutivedb, CLASSBackEnd* backend,
		cSecond creationtime, cSecond lifetime,
		double power, double HMMass, double BurnUp, double CapacityFactor);
	//}
	
	//{
	/// Special Constructor for fixed fuel using Power and Burn-Up.
	/*!
	 Make a new reactor
	 \param log: used for the log,
	 \param evolutivedb: EvolutionData describing the evolution of the fuel,
	 \param backend: CLASSBackend which get the fuel after the cooling,
	 \param creationtime: creation time in [s],
	 \param lifetime: working time duration in [s]
	 \param Power: Thermal power of the reactor in [W],
	 \param HMMass: Mass of Heavy Metal in the Reactor in [t] of heavy metal,
	 \param BurnUp: Burnup reach by the fuel at the end of the cycle in [GWd/t],
	 */
	Reactor(CLASSLogger* log, EvolutionData* evolutivedb, CLASSBackEnd* backend,
		cSecond creationtime, cSecond lifetime,
		cSecond cycletime, double HMMass, double BurnUp);
	//}
#endif
	~Reactor();	///< Normal Destructor
	
	//@}
	
	
	
	
	//********* Get Method *********//
	
	/*!
	 \name Get Method
	 */
	//@{
	
	IsotopicVector 	GetIVReactor()		const	{ return GetInsideIV(); } 	//!< Return the IV contain in the Reactor
	IsotopicVector	GetIVBeginCycle()	const	{ return fIVBeginCycle; }	//!< Return the Starting Cycle IV
	//!< (Note : IVBegin != IVIn, only if using ReactorScheduler)
	IsotopicVector	GetIVOutCycle()		const	{ return fIVOutCycle; }		//!< Return the Out Cycle IV
	IsotopicVector	GetIVInCycle()		const	{ return fIVInCycle; }		//!< Return the In Cycle IV
	//!< (Note : IVIn != IVBegin, only if using ReactorScheduler)
	
	cSecond GetNextCycleTime(cSecond time);
	
	bool	IsFuelFixed()	const	{ return fFixedFuel; }		//!< True if using fixed fuel, False otherwise
	double	GetHeavyMetalMass() const { return fHeavyMetalMass; }	//!< Return the HeavyMetal Mass in the Core at the begining of the cycle
    double  GetBurnUp() const   { return fBurnUp; }    			 //!< Return the Burn Up of the Fuel at the end of the cycle
	double	GetCapacityFactor()	const	{ return fCapacityFactor; }	//!< Return the capacity factor of the reactor
	double	GetPower()	const	{ return fPower; } 			//!< Return the cycle time of the Reactor
	
	map<cSecond, double> GetPowerEvolution(){return fPowerEvolution;} //!<Get Reactor Power Evolution
	map<cSecond, double> GetHMMassEvolution(){return fHMMassEvolution;} //!<Get Reactor HMMassEvolution
	map<cSecond, double> GetSchedulePowerEvolution(){return fSchedulePowerEvolution;} //!<Get Reactor theoritical Power evolution following Scheduler
	map<cSecond, double> GetScheduleHMMassEvolution(){return fScheduleHMMassEvolution;} //!<Get Reactor theoritical HMMass Evolution following Scheduler


#ifndef __ROOTCLING__
	
	EvolutionData	GetEvolutionDB()	const	{ return fEvolutionDB; }	//!< Return the Evolution database of the fuel
	CLASSBackEnd*	GetOutBackEndFacility()	const	{ return fOutBackEndFacility; }	//!< Return the pointer to Associeted BackEnd Facility
	FabricationPlant*	GetFabricationPlant()	const	{ return fFabricationPlant; }			//!< Return the pointer to the FabricationPlant
	ReactorScheduler*	GetScheduler()		const	{ return fReactorScheduler; }			///< return the ReactorScheduler
	
#endif
	//@}
	
	
	
	
	//********* Set Method *********//
	
	/*!
	 \name Set Method
	 */
	//@{
	
	void SetHMMass(double Mass)		{fHeavyMetalMass = Mass;}							//!< Set the heavy metal mass in the core at the begining of the cycle
	void SetCycleTime(double cycletime);											//!< Set the cycle time (Power fixed)
  	void SetPower(double Power);                   						 				//!< Set the power (burnup cte)
	void SetBurnUp(double BU);					   							 //!< Set the burnUp reach at end of cycle (Power cte)
	
	void SetIVReactor(IsotopicVector isotopicvector) { fInsideIV = isotopicvector; }	//!< Set the IV inside the Reactor core
	void SetIVBeginCycle(IsotopicVector isotopicvector) { fIVBeginCycle = isotopicvector;}	//!< Set the IV at the beginging of the Reactor cycle
	void SetIVOutCycle(IsotopicVector isotopicvector){ fIVOutCycle = isotopicvector;}	//!< Set the IV Going out at the end of the cycle
	void SetIVInCycle(IsotopicVector isotopicvector) { fIVInCycle = isotopicvector;}	//!< Set the IV coming In at the beginning of the cycle
	
	
	
	
#ifndef __ROOTCLING__
	
	void SetReactorScheduler(ReactorScheduler* reactorscheduler)	{ fReactorScheduler = reactorscheduler; }	//!< Set the ReactorScheduler
	void	SetOutBackEndFacility(CLASSBackEnd* pool)	{ fOutBackEndFacility = pool; }				//!< Return the pointer to OutBackEnd Facility
	void	SetStorage(Storage* storage)			{ fStorage = storage; fIsStorage = true;}			//!< Set the pointer to the Storage
	void	SetFabricationPlant(FabricationPlant* FP)	{ fFabricationPlant = FP;}					//!< Set the Pointer to the FabricationPlant
	void	SetEvolutionDB(EvolutionData evolutionDB);									//!< Set the pointer to the DB evolution of the Reactor

#endif
	
	using CLASSFacility::SetName;
	using CLASSFacility::GetName;
	
	//@}
	
	
	
	
	//********* Evolution & Modification Method *********//
	
	/*!
	 \name Evolution & Modification Method
	 */
	//@{
	
	void Evolution(cSecond t);						//!< Performs the Evolution until time t
	void Dump();								//!< Write modification (IV In/Out, filling the TF...)
	void SetNewFuel(EvolutionData ivdb);					//!< Change the Evolutive DB of the Reactor
	void CheckListConsistency(PhysicsModels* fueltypeDB, FabricationPlant* fabricationplant); //!< Check if FP/EqM lists are identical
#ifndef __ROOTCLING__

	void AddScheduleEntry(cSecond time,  ReactorModel* Model, double BurnUp, double Power, double HMMass)	//!< Add A new schedule entry 
	{ 	if(Model->GetPhysicsModels() && !fFabricationPlant )
			{ERROR<<"You must add a FabicationPlant to Reactor " <<GetName()
		   <<"\n"<< "Use YourReactor->AddFabricationPlant(YourFP);"<<endl;exit(0);}
		fReactorScheduler->AddEntry( time, Model, BurnUp, Power*fCapacityFactor, HMMass); 
	}			//!< Change the Model and/or BurnUp and/or Power and/or HMMass from time time
	void AddScheduleEntry(cSecond time,  EvolutionData* Model, double BurnUp, double Power, double HMMass)
	{  fReactorScheduler->AddEntry( time, new ReactorModel(Model), BurnUp, Power*fCapacityFactor, HMMass); }			//!< Change the Model and/or BurnUp and/or Power and/or HMMass from time time
	void AddScheduleEntry(cSecond time,  PhysicsModels* Model, double BurnUp, double Power, double HMMass)
	{	if(!fFabricationPlant) 
			{ERROR<<"You must add a FabicationPlant to Reactor " <<GetName()
		   <<"\n"<< "Use YourReactor->AddFabricationPlant(YourFP);"<<endl;exit(0);}
		fReactorScheduler->AddEntry( time, new ReactorModel(Model), BurnUp,  Power*fCapacityFactor, HMMass); 
	}			//!< Change the Model and/or BurnUp and/or Power and/or HMMass from time time
#endif

	//@}
	
	
	
	protected :
	
	bool		fFixedFuel;		//!< true if the fuel is fixed (not reprocessed)
	bool		fIsStorage;		//!< true if a storage has been define (to approximate the reprocessing using fixed fuel)
	
	//********* Internal Parameter *********//
	
	double 		fPower;			///< Power (in Watt)
	double 		fElectricPower;		///< ElectrocPower (in Watt)
	double 		fEfficiencyFactor;	///< ElectrocPower (in Watt)
	double		fCapacityFactor;    ///< Capacity factor  [0-1]

	IsotopicVector	fIVBeginCycle;		///< Fuel IV at the beginning of a cycle
	IsotopicVector	fIVInCycle;		///< IVBegin add at the beginning of the cycle
	IsotopicVector	fIVOutCycle;		///< IV wich get out at the end of a cycle
	
#ifndef __ROOTCLING__
	EvolutionData	fEvolutionDB;			//!< Pointer to the actual evolution DataBase
	
	CLASSBackEnd*	fOutBackEndFacility;	//!< Pointer to the BackEnd Facility which collect the spend fuel
		
	FabricationPlant*	fFabricationPlant;		//!< Pointer to the FabricationPlant
	
	Storage*	fStorage;		//!< Pointer to the Stock (only for reprocessing fuel in fixed base...)

    ReactorScheduler*   fReactorScheduler;      //!< Pointer to the ReactorScheduler
	
#endif
	//********* Unfixed Fuel Parameter *********//
	
	
	double			fHeavyMetalMass;		///< In tons
	double			fBurnUp;			///< In GWd/tHM
	
	ClassDef(Reactor,3);
};


#endif
