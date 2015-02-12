
#ifndef _CLASSFACILITY_
#define _CLASSFACILITY_


/*!
 \file
 \brief Header file for CLASSFacility class.

 */

#include <string>
#include <fstream>

#include "CLASSObject.hxx"
#include "IsotopicVector.hxx"
#include "DecayDataBank.hxx"

#include "TNamed.h"

using namespace std;
typedef long long int cSecond;

class Scenario;
//-----------------------------------------------------------------------------//
//! Defines the common properties of all facilities

/*!
 Define a CLASS Facility.
 The aim of these class is to gather all the commom properties of the facilities.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________




class CLASSFacility : public CLASSObject
{
public :
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	//{
	/// Normal Constructor.
	/*!
	Make a new Facility
	\param type identification type of the facility :
	\li 4 Reactor,
	\li 8 Pool,
	\li 16 FabricationPlant.
	 */
	
	CLASSFacility(int type = 0);
	//}

	//{
	/// Special Constructor.
	/*!
	 Make a new Facility
	 \param log : used for the log.
	 \param type identification type of the facility :
	 \li 4 Reactor,
	 \li 8 Pool,
	 \li 16 FabricationPlant.
	 
	 */
	CLASSFacility(CLASSLogger* log, int type = 0);
	//}
	
	//{
	/// Special Constructor.
	/*!
	 Make a new Facility
	 \param log : used for the log.
	 \param cycletime duration of the cycle [s],
	 \param type identification type of the facility :
	 \li 4 Reactor,
	 \li 8 Pool,
	 \li 16 FabricationPlant.
	 
	 */
	CLASSFacility(CLASSLogger* log, cSecond cycletime, int type = 0);
	//}
	
	//{
	/// Special Constructor.
	/*!
	 Make a new Facility
	 \param log : used for the log.
	 \param creationtime creation date of the Facility [s],
	 \param lifetime operating duration [s],
	 \param type identification type of the facility :
	 \li 4 Reactor,
	 \li 8 Pool,
	 \li 16 FabricationPlant.
	 
	 */
	CLASSFacility(CLASSLogger* log, cSecond creationtime, cSecond lifetime, int type = 0);
	//}
	
	//{
	/// Special Constructor.
	/*!
	 Make a new Facility
	 \param log : used for the log.
	 \param creationtime creation date of the Facility [s],
	 \param lifetime operating duration [s],
	 \param cycletime duration of the cycle [s],
	 \param type identification type of the facility :
	 \li 4 Reactor,
	 \li 8 Pool,
	 \li 16 FabricationPlant.
	 
	 */
	CLASSFacility(CLASSLogger* log, cSecond startingtime, cSecond lifetime, cSecond cycletime, int type = 0);
	//}
	
	//********* Get Method *********//
	/*!
	 \name Get Function
	 */
	//@{

	int 		GetId()			const	{ return fId; }			//!< Return the Facility Parc'Is
	IsotopicVector 	GetInsideIV()		const	{ return fInsideIV; } 		//!< Return the IV contained in the Facility

	int		GetFacilityType()	const	{ return fFacilityType; }	//!< Return the Facility Type id

	cSecond		GetInternalTime()	const	{ return fInternalTime; }	//!< Return Creation Time

	cSecond		GetCycleTime()		const	{ return fCycleTime; } 		//!< Return the cycle time of the Facility
	cSecond 	GetCreationTime()	const	{ return fCreationTime; }	//!< Return the creation time of the Facility
	cSecond 	GetLifeTime()		const	{ return fLifeTime; }		//!< Return the life time of the Facility
	Scenario*	GetParc()			{ return fParc; }		//!< return the pointer to the Park


	IsotopicVector GetCumulativeIVIn() { return fCumulativeIVIn;}			//!< return the cumulative sum of all incoming IV
	IsotopicVector GetCumulativeIVOut() { return fCumulativeIVOut;}			//!< return the cumulative sum of all outcoming IV
	//@}
	
		//********* Set Method *********//
	/*!
	 \name Set Function
	 */
	//@{
	void	SetId(int id)			{ fId = id; }				//!< Set The Facility Parc'Id
	void	SetParc(Scenario* parc)		{ fParc = parc; }			//!< Set the Pointer to the Parc
	void	SetFacilityType(int type)	{ fFacilityType = type; }		//!< Set the facility type :
											/// \li 2 reactor Studown
											/// \li 4 start/End of reactor cycle,
											/// \li 8 end of Cooling,
											/// \li 16 fuel Fabrication
	using CLASSObject::SetName;
	using CLASSObject::GetName;


	void SetInsideIV(IsotopicVector isotopicvector)	{ fInsideIV = isotopicvector; }	//!< Set the IV inside the Facility
	void SetCreationTime(double creationtime)	{ fCreationTime = (cSecond)creationtime;} //!< Set the creation Time
	void SetLifeTime(double lifetime)		{ fLifeTime = (cSecond)lifetime; }	//!< Set the life time of the facility
	virtual void SetCycleTime(double cycletime)	{ fCycleTime = (cSecond)cycletime; }	//!< Set the cycle time (Cycle of the loading Plan)
	void SetInCycleTime(double incycletime)		{ fInCycleTime = (cSecond)incycletime; fIsStarted = true; }	//!< Set the cycle time (Cycle of the loading Plan)
	void SetInternalTime(double internaltime)	{ fInternalTime = (cSecond)internaltime; }	//!< Set the cycle time (Cycle of the loading Plan)

	//@}


	/*!
	 \name Evolution Method
	 */
	//@{

	void AddCumulativeIVIn(IsotopicVector IV) { fCumulativeIVIn += IV;}		//!< Add the Input IsotopicVector in the cumulative IV IN
	void AddCumulativeIVOut(IsotopicVector IV) { fCumulativeIVOut += IV;}		//!< Add the Input IsotopicVector in the cumulative IV OUT
	virtual void Evolution(cSecond t)	= 0;	//!< Performs the Evolution to time t
	virtual void Dump()			{ }	//!< Write Modification (IV In/Out, filling the TF...)
		
	//@}
protected :
	bool		fIsStarted;		///< True if Running, False Otherwise
	bool		fIsShutDown;		///< True if the facility is stoped, False Otherwise
	bool		fIsAtEndOfCycle;	///< True if Reaching the end of a Facility cycle

		
	cSecond		fInternalTime;		///< Internal Clock [s]
	cSecond		fInCycleTime;		///< Time spent since the beginning of the last Cycle [s]
	cSecond		fCycleTime;		///< Cycle duration Time [s]

	IsotopicVector	fInsideIV;		///< All IV in the Facility (fuel for reactor, total for all others...)
	IsotopicVector	fCumulativeIVIn;	///< All IV in the Facility (fuel for reactor, total for all others...)
	IsotopicVector	fCumulativeIVOut;	///< All IV in the Facility (fuel for reactor, total for all others...)

		//********* Internal Parameter *********//
private :
	int		fId;			//!< Identity of the Facility inside the Scenario
	int		fFacilityType;		///< Type of facility :
						/// \li 4 reactor,
						/// \li 8 Pool,
						/// \li 16 FabricationPlant.


	Scenario*	fParc;			//!< Pointer to the main Scenario

	cSecond		fCreationTime;		///< CLASS Universal Time of Creation [s]
	cSecond		fLifeTime;		///< Time of life Of the Reactor (operating's duration) [s]

	ClassDef(CLASSFacility,1);
};

#endif

