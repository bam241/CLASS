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


class CLASS;
class CLASSBackEnd;
//class Pool;
class EvolutionData;
class FuelDataBank;
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
	 \param Pool Pool used for the cooling of the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 */
	Reactor(LogFile* log, FuelDataBank* 	fueltypeDB,
		FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
//		FabricationPlant* fabricationplant, Pool* Pool,
		double creationtime , double lifetime);
	//}

	//{
	 /// Special Constructor for reprocessed fuel using cycletime and Burn-Up.
	 /*!
	  Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param fueltypeDB Databank describing the evolution of the fuel
	 \param Pool Pool used for the cooling of the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 \param cycletime duration of a cycle
	 \param HMMass Mass of Heavy Metal in the Reactor
	 \param BurnUp Burnup reach by the fuel at the end of the cycle
	 */
	Reactor(LogFile* log, FuelDataBank* 	fueltypeDB,
//		FabricationPlant* fabricationplant, Pool* Pool,
		FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		double creationtime , double lifetime, double cycletime,
		double HMMass, double BurnUp);
	//}

	//{
	/// Special Constructor for reprocessed fuel using Power and Burn-Up.
	/*!
	 Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param fueltypeDB Databank describing the evolution of the fuel
	 \param Pool Pool used for the cooling of the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 \param Power Thermal power of the reactor
	 \param HMMass Mass of Heavy Metal in the Reactor
	 \param BurnUp Burnup reach by the fuel at the end of the cycle
	 \param ChargeFactor effective charge of the reactor.
	 */
	Reactor(LogFile* log, FuelDataBank* 	fueltypeDB,
//		FabricationPlant* fabricationplant, Pool* Pool,
		FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		double creationtime , double lifetime,
		double Power, double HMMass, double BurnUp, double ChargeFactor);
	//}

	//{
	/// Special Constructor for fixed fuel using Power and Burn-Up.
	/*!
	 Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param evolutivedb EvolutionData describing the evolution of the fuel
	 \param Pool Pool used for the cooling of the fuel after iradiation
	 \param creationtime creation time
	 \param lifetime working time duration.
	 \param Power Thermal power of the reactor
	 \param HMMass Mass of Heavy Metal in the Reactor
	 \param BurnUp Burnup reach by the fuel at the end of the cycle
	 \param ChargeFactor effective charge of the reactor.
	 */
	Reactor(LogFile* log, EvolutionData evolutivedb, CLASSBackEnd* Pool,
//	Reactor(LogFile* log, EvolutionData evolutivedb, Pool* Pool,
		double creationtime, double lifetime,
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
	

	EvolutionData	GetEvolutionDB()		const	{ return fEvolutionDB; }	//!< Return the Evolution database of the Fuel
	FuelDataBank*	GetFuelType()	const	{ return fFuelTypeDB; }		//!< Return the Fuel Type DB of the reactor

//	Pool*		GetAssociedPool()	const	{ return fAssociedPool; }	//!< Return the pointer to Associeted Pool
	CLASSBackEnd*		GetOutBackEndFacility()	const	{ return fOutBackEndFacility; }	//!< Return the pointer to Associeted BackEnd Facility
	FabricationPlant*	GetFabricationPlant()	const	{ return fFabricationPlant; }	//!< Return the Pointer to the FabricationPlant

	bool	IsFuelFixed()		const	{ return fFixedFuel; }		//!< True if using fixed Fuel, False otherwise
	double	GetHeavyMetalMass()	const	{ return fHeavyMetalMass; }	//!< Return the HeavyMetal Mass in the Core at the begining of the cycle
	double	GetBurnUp()		const	{ return fBurnUp; }		//!< Return the Burn Up of the Fuel at the end of the cycle
	double	GetPower()		const	{ return fPower; } 		//!< Return the cycle time of the Reactor

	//@}




//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{
	void	SetOutBackEndFacility(CLASSBackEnd* pool)	{ fOutBackEndFacility = pool; }	//!< Return the pointer to Associeted TF
//	void	SetAssociedPool(Pool* pool)	{ fAssociedPool = pool; }	//!< Return the pointer to Associeted TF

	void	SetStorage(Storage* storage)		{ fStorage = storage; fIsStorage = true;}	//!< Set the Pointer to the Storage

	void	SetIVReactor(IsotopicVector isotopicvector)	{ fInsideIV = isotopicvector; }	//!< Set the IV inside the Reactor Core
	void	SetIVBeginCycle(IsotopicVector isotopicvector)	{ fIVBeginCycle = isotopicvector; }	//!< Set the IV at the Beginging of the Reactor Cycle
	void	SetIVOutCycle(IsotopicVector isotopicvector)	{ fIVOutCycle = isotopicvector; }	//!< Set the IV Going Out at the End of the Cycle
	void	SetIVInCycle(IsotopicVector isotopicvector)	{ fIVInCycle = isotopicvector; }	//!< Set the IV Coming In at the Beginning of the Cycle
	
	void	SetCycleTime(double cycletime);					//!< Set the Power time (Cycle of the loading Plan)
	
	void	SetEvolutionDB(EvolutionData evolutionDB);				//!< Set the Pointer to the DB Evolution of the Reactor
	void	SetPower(double Power);						//!< Set the Power
	void	SetHMMass(double Mass)		{fHeavyMetalMass = Mass;}	//!< Set the HeavyMetal Mass in the Core at the begining of the cycle
	void	SetBurnUp(double BU)		{fBurnUp = BU;}			//!< Set the the Burn Up of the Fuel at the end of the cycle
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
//	Pool*		fAssociedPool;		//!< Pointer to the TF which collect the spend fuel
	CLASSBackEnd*	fOutBackEndFacility;		//!< Pointer to the TF which collect the spend fuel
	Storage*	fStorage;		//!< Pointer to the Stock
						
	EvolutionData	fEvolutionDB;			//!< Pointer to the Evolution DataBase
	FuelDataBank* 	fFuelTypeDB;	//! Pointer to a Fuel Type Database
	
	double 		fPower;			///< Power (in Watt)
	
	IsotopicVector	fIVBeginCycle;		///< Fuel IV at the Beginning of a Cycle
	IsotopicVector	fIVInCycle;		///< IVBegin add at the Beginning of the Cycle
	IsotopicVector	fIVOutCycle;		///< IV wich get out at the End of a Cycle


//********* Unfixed Fuel Parameter *********//


	FabricationPlant*	fFabricationPlant;		//!< Poitner to the FabricationPlant
	double			fHeavyMetalMass;		///< In tons
	double			fBurnUp;			///< In GWd/tHM

 	ClassDef(Reactor,2);
 };


#endif
