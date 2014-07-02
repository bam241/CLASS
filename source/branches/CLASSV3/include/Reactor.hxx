#ifndef __Reactor_HXX__
#define __Reactor_HXX__

/*!
 \file
 \brief Header file for reactor classes. 
 */

#include <string>
#include <map>

#include "CLASSFacility.hxx"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"

using namespace std;
typedef long long int cSecond;


class CLASSBackEnd;
//class Pool;
class EvolutionData;
class PhysicModels;
class FabricationPlant;
class Storage;
class LogFile;

//-----------------------------------------------------------------------------//
/*!
 Define a reactor.
 The aim of this class is to deal the evolution of the fuel inside a reactor.
 The fuel state of the reactor is describe in the IsotopicVector. Its evolution is contain in the EvolutionData

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
	/// LogFile Constructor.
	/*!
	 Use create an empty Reactor loading a LogFile
	 \param LogFile LogFile used for the log...
	 */
	Reactor(LogFile* log);
	//}

	//{
	/// Special Constructor for reprocessed fuel.
	/*!
	 Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param fueltypeDB Databank describing the evolution of the fuel
	 \param CLASSBAckEnd Pool used facility wich get the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 */
	Reactor(LogFile* log, PhysicModels* 	fueltypeDB,
		FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		cSecond creationtime , cSecond lifetime);
	//}

	//{
	 /// Special Constructor for reprocessed fuel using cycletime and Burn-Up.
	 /*!
	  Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param fueltypeDB Databank describing the evolution of the fuel
	 \param CLASSBAckEnd Pool used facility wich get the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 \param cycletime duration of a cycle
	 \param HMMass Mass of Heavy Metal in the Reactor
	 \param BurnUp Burnup reach by the fuel at the end of the cycle
	 */
	Reactor(LogFile* log, PhysicModels* 	fueltypeDB,
		FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		cSecond creationtime , cSecond lifetime, cSecond cycletime,
		double HMMass, double BurnUp);
	//}

	//{
	/// Special Constructor for reprocessed fuel using Power and Burn-Up.
	/*!
	 Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param fueltypeDB Databank describing the evolution of the fuel
	 \param CLASSBAckEnd Pool used facility wich get the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 \param Power Thermal power of the reactor
	 \param HMMass Mass of Heavy Metal in the Reactor
	 \param BurnUp Burnup reach by the fuel at the end of the cycle
	 \param ChargeFactor effective charge of the reactor.
	 */
	Reactor(LogFile* log, PhysicModels* 	fueltypeDB,
		FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		cSecond creationtime , cSecond lifetime,
		double Power, double HMMass, double BurnUp, double ChargeFactor);
	//}

	//{
	/// Special Constructor for fixed fuel using Power and Burn-Up.
	/*!
	 Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param evolutivedb EvolutionData describing the evolution of the fuel
	 \param CLASSBAckEnd Pool used facility wich get the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 \param Power Thermal power of the reactor
	 \param HMMass Mass of Heavy Metal in the Reactor
	 \param BurnUp Burnup reach by the fuel at the end of the cycle
	 \param ChargeFactor effective charge of the reactor.
	 */
	Reactor(LogFile* log, EvolutionData evolutivedb, CLASSBackEnd* Pool,
		cSecond creationtime, cSecond lifetime,
		double power, double HMMass, double BurnUp, double ChargeFactor = 1);
	//}

	~Reactor();	///< Normal Destructor

	//@}




//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{
	
	IsotopicVector 	GetIVReactor()		const	{ return GetInsideIV(); } 	//!< Return the IV contain in the Reactor
	IsotopicVector	GetIVBeginCycle()	const	{ return fIVBeginCycle; }	//!< Return the Starting Cycle IV
											//!< (Note : IVBegin != IVIn, only if using charging plan)
	
	IsotopicVector	GetIVOutCycle()		const	{ return fIVOutCycle; }		//!< Return the Out Cycle IV
	IsotopicVector	GetIVInCycle()		const	{ return fIVInCycle; }		//!< Return the In Cycle IV
											//!< (Note : IVIn != IVBegin, only if using charging plan)
	

	EvolutionData	GetEvolutionDB()	const	{ return fEvolutionDB; }	//!< Return the Evolution database of the Fuel
	PhysicModels*	GetFuelType()		const	{ return fFuelTypeDB; }		//!< Return the Fuel Type DB of the reactor

	CLASSBackEnd*		GetOutBackEndFacility()	const	{ return fOutBackEndFacility; }	//!< Return the pointer to Associeted BackEnd Facility
	FabricationPlant*	GetFabricationPlant()	const	{ return fFabricationPlant; }	//!< Return the Pointer to the FabricationPlant

	bool	IsFuelFixed()		const	{ return fFixedFuel; }		//!< True if using fixed Fuel, False otherwise
	double	GetHeavyMetalMass()	const	{ return fHeavyMetalMass; }	//!< Return the HeavyMetal Mass in the Core at the begining of the cycle
	double	GetBurnUp()		const	{ return fBurnUp; }		//!< Return the Burn Up of the Fuel at the end of the cycle
	double	GetPower()		const	{ return fPower; } 		//!< Return the cycle time of the Reactor

#ifndef __CINT__
	map<cSecond, pair<EvolutionData, double> >	GetLoadingPlan()		const
						{ return fLoadingPlan; }	//!< return the LoadingPlan
	map<cSecond, pair<EvolutionData, double> >::iterator	GetNextPlan()	const
						{ return fNextPlan; }	//!< return the next fuel in the Plan

#endif
	//@}




//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{
	void	SetOutBackEndFacility(CLASSBackEnd* pool)
					{ fOutBackEndFacility = pool; }	//!< Return the pointer to OutBackEnd Facility

	void	SetStorage(Storage* storage)
					{ fStorage = storage; fIsStorage = true;}	//!< Set the Pointer to the Storage

	void	SetHMMass(double Mass)		{fHeavyMetalMass = Mass;}	//!< Set the HeavyMetal Mass in the Core at the begining of the cycle

	void	SetIVReactor(IsotopicVector isotopicvector)
					{ fInsideIV = isotopicvector; }		//!< Set the IV inside the Reactor Core
	void	SetIVBeginCycle(IsotopicVector isotopicvector)
					{ fIVBeginCycle = isotopicvector; }	//!< Set the IV at the Beginging of the Reactor Cycle
	void	SetIVOutCycle(IsotopicVector isotopicvector)
					{ fIVOutCycle = isotopicvector; }	//!< Set the IV Going Out at the End of the Cycle
	void	SetIVInCycle(IsotopicVector isotopicvector)
					{ fIVInCycle = isotopicvector; }	//!< Set the IV Coming In at the Beginning of the Cycle

	void	SetEvolutionDB(EvolutionData evolutionDB);			//!< Set the Pointer to the DB Evolution of the Reactor

	void	SetCycleTime(double cycletime);					//!< Set the Cycle time (Power fixed)
	void	SetPower(double Power);						//!< Set the Power (BurnUp cte)
	void	SetBurnUp(double BU);						//!< Set the BurnUp reach at end of cycle (Power cte)




	void	SetLoadingPlan(map<cSecond, pair<EvolutionData, double> > loadingplan)
					{ fLoadingPlan = loadingplan; fNextPlan = fLoadingPlan.begin(); }
										//!< Set a LaodingPlan to change the Fuel after some cycle
	//@}




//********* Evolution & Modification Method *********//

	/*!
	 \name Evolution & Modification Method
	 */
	//@{
	
	void Evolution(cSecond t);						//!< Performe the Evolution until the Time t
	void Dump();								//!< Write Modification (IV In/Out, filling the TF...)
	void SetNewFuel(EvolutionData ivdb);					//!< Change the Evolutive DB of the Reactor
	
	//@}



protected :
	
	bool		fFixedFuel;		//!< true if the fuel is fixed (not reprocessed)
	bool		fIsStorage;		//!< true if a storage has been define (to approximate the reprocessing using fixed fuel)
	
//********* Internal Parameter *********//
	CLASSBackEnd*	fOutBackEndFacility;	//!< Pointer to the BackEnd Facility which collect the spend fuel
	Storage*	fStorage;		//!< Pointer to the Stock (only for reprocessing fuel in fixed base...)
						
	EvolutionData	fEvolutionDB;		//!< Pointer to the Evolution DataBase
	PhysicModels* 	fFuelTypeDB;		//! Pointer to a Fuel Type Database
	
	double 		fPower;			///< Power (in Watt)
	
	IsotopicVector	fIVBeginCycle;		///< Fuel IV at the Beginning of a Cycle
	IsotopicVector	fIVInCycle;		///< IVBegin add at the Beginning of the Cycle
	IsotopicVector	fIVOutCycle;		///< IV wich get out at the End of a Cycle

#ifndef __CINT__
	map<cSecond, pair<EvolutionData, double> >	fLoadingPlan;	///< Loading PLan to change the EvolutionData (and the associetedBurnup) according to the Plan
	map<cSecond, pair<EvolutionData, double> >::iterator	fNextPlan;	///< Next EvolutionData, and time until it should be load (at the end of the last cycle)
#endif
//********* Unfixed Fuel Parameter *********//


	FabricationPlant*	fFabricationPlant;		//!< Poitner to the FabricationPlant
	double			fHeavyMetalMass;		///< In tons
	double			fBurnUp;			///< In GWd/tHM

 	ClassDef(Reactor,3);
 };


#endif
