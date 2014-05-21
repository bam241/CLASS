
#ifndef _CLASSFACILITY_HXX
#define _CLASSFACILITY_HXX


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

class CLASS;
//-----------------------------------------------------------------------------//
/*!
 Define a CLASS Facility.
 The aim of these class is synthetyse all the commum properties of the nuclear facilities.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________




class CLASSFacility : public CLASSObject
{
public :
		///< Normal Constructor.
	CLASSFacility();
	
		//********* Get Method *********//
	/*!
	 \name Get Function
	 */
	//@{

	int 		GetId()			const	{ return fId; }			//!< Return the Facility Parc'Is
	IsotopicVector 	GetInsideIV()		const	{ return fInsideIV; } 		//!< Return the IV contain in the Facility

	int		GetFacilityType()	const	{ return fFacilityType; }

	cSecond		GetInternalTime()	const	{ return fInternalTime; }	//!< Return Creation Time

	cSecond		GetCycleTime()		const	{ return fCycleTime; } 		//!< Return the cycle time of the Facility
	cSecond 	GetCreationTime()	const	{ return fCreationTime; }	//!< Return the creation time of the Facility
	cSecond 	GetLifeTime()		const	{ return fLifeTime; }		//!< Return the life time of the Facility
	CLASS*		GetParc()			{ return fParc; }		//!< return the pointer to the Park


	IsotopicVector GetCumulativeIVIn() { return fCumulativeIVIn;}			//!< return the culative sum of all incoming IV
	IsotopicVector GetCumulativeIVOut() { return fCumulativeIVOut;}			//!< return the culative sum of all outcoming IV
	//@}
	
		//********* Set Method *********//
	/*!
	 \name Set Function
	 */
	//@{
	void	SetId(int id)			{ fId = id; }				//!< Set The Facility Parc'Id
	void	SetParc(CLASS* parc)		{ fParc = parc; }			//!< Set the Pointer to the Parc
	void	SetFacilityType(int type)		{ fFacilityType = type; }	//!< Set the facility type :
											/// \li 2 reactor Studown
											/// \li 4 start/End of reactor cycle,
											/// \li 8 end of Cooling,
											/// \li 16 fuel Fabrication


	void SetInsideIV(IsotopicVector isotopicvector)	{ fInsideIV = isotopicvector; }	//!< Set the IV inside the Facility Core
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

	void AddCumulativeIVIn(IsotopicVector IV) { fCumulativeIVIn += IV;}		//!< Add the Input IsotopicVector the The cumulative IV IN
	void AddCumulativeIVOut(IsotopicVector IV) { fCumulativeIVOut += IV;}		//!< Add the Input IsotopicVector the The cumulative IV OUT
	virtual void Evolution(cSecond t)	{ }	//!< Performe the Evolution to the Time t
	virtual void Dump()			{ }	//!< Write Modification (IV In/Out, filling the TF...)
		
	//@}
protected :
	bool		fIsStarted;		///< True if Running, False Otherwise
	bool		fIsShutDown;		///< True if the facility is stoped, False Otherwise
	bool		fIsAtEndOfCycle;	///< True if Reaching the End of a Facility Cycle

		
	cSecond		fInternalTime;		///< Internal Clock
	cSecond		fInCycleTime;		///< Time spend since the beginning of the last Cycle
	cSecond		fCycleTime;		///< Cycle duration Time

	IsotopicVector	fInsideIV;		///< All IV in the Facility (fuel for reactor, total for all others...)
	IsotopicVector	fCumulativeIVIn;	///< All IV in the Facility (fuel for reactor, total for all others...)
	IsotopicVector	fCumulativeIVOut;	///< All IV in the Facility (fuel for reactor, total for all others...)

		//********* Internal Parameter *********//
private :
	int		fId;			//!< Identity of the Facility inside the Parc
	int		fFacilityType;		///< Type of facility :
						/// \li 2 reactor Studown
						/// \li 4 start/End of reactor cycle,
						/// \li 8 end of Cooling,
						/// \li 16 fuel Fabrication

	CLASS*		fParc;			//!< Pointer to the main Parc

	cSecond		fCreationTime;		///< CLASS Universal Time of Creation
	cSecond		fLifeTime;		///< LifeTime Of the Reactor (Operating's Duration)

	ClassDef(CLASSFacility,1);
};

#endif

