
#ifndef _CLASSFACILITY_HXX
#define _CLASSFACILITY_HXX


/*!
 \file
 \brief Header file for CLSSFacility class.

 */

#include <string>
#include <fstream>

#include "CLSSObject.hxx"
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




class CLSSFacility : public CLSSObject
{
public :
		///< Normal Constructor.
	CLSSFacility();
	
		//********* Get Method *********//
	/*!
	 \name Get Function
	 */
	//@{

	int 		GetId()			const	{ return fId; }			//!< Return the Facility Parc'Is
	IsotopicVector 	GetInsideIV()		const	{ return fInsideIV; } 		//!< Return the IV contain in the Facility

	cSecond	GetInternalTime() const		{ return fInternalTime; }	//!< Return Creation Time

	cSecond		GetCycleTime()		const	{ return fCycleTime; } 		//!< Return the cycle time of the Facility
	cSecond 	GetCreationTime()	const	{ return fCreationTime; }	//!< Return the creation time of the Facility
	cSecond 	GetLifeTime()		const	{ return fLifeTime; }		//!< Return the life time of the Facility
	CLASS*		GetParc()			{ return fParc; }		//! retiurn the pointer to the Park
	DecayDataBank*	GetDecayDataBank()		{ return fDecayDataBase;}	//!< Return the pointer to the decay DataBank


	IsotopicVector GetCumulativeIVIn() { return fCumulativeIVIn;}
	IsotopicVector GetCumulativeIVOut() { return fCumulativeIVOut;}
	//@}
	
		//********* Set Method *********//
	/*!
	 \name Set Function
	 */
	//@{
	void SetId(int id)			{ fId = id; }				//!< Set The Facility Parc'Id
	void SetParc(CLASS* parc)		{ fParc = parc; }			//!< Set the Pointer to the Parc
	
	void SetInsideIV(IsotopicVector isotopicvector)	{ fInsideIV = isotopicvector; }	//!< Set the IV inside the Facility Core
	void SetCreationTime(double creationtime)	{ fCreationTime = (cSecond)creationtime;}
	void SetLifeTime(double lifetime)		{ fLifeTime = (cSecond)lifetime; }	//!< Set the life time of the facility
	virtual void SetCycleTime(double cycletime)	{ fCycleTime = (cSecond)cycletime; }	//!< Set the cycle time (Cycle of the loading Plan)
	void SetInCycleTime(double incycletime)		{ fInCycleTime = (cSecond)incycletime; fIsStarted = true; }	//!< Set the cycle time (Cycle of the loading Plan)
	void SetInternalTime(double internaltime)	{ fInternalTime = (cSecond)internaltime; }	//!< Set the cycle time (Cycle of the loading Plan)

	void SetDecayDataBank(DecayDataBank* decayDB) {fDecayDataBase = decayDB;} //! Set the Decay DataBank
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
	bool		fShutDown;		///< True if the facility is stoped, False Otherwise
	bool		fEndOfCycle;		///< True if Reaching the End of a Facility Cycle

		
	cSecond		fInternalTime;		///< Internal Clock
	cSecond		fInCycleTime;		///< Time spend since the beginning of the last Cycle
	cSecond		fCycleTime;		///< Cycle duration Time

	IsotopicVector	fInsideIV;		///< All IV in the Facility (fuel for reactor, total for all others...)
	IsotopicVector	fCumulativeIVIn;	///< All IV in the Facility (fuel for reactor, total for all others...)
	IsotopicVector	fCumulativeIVOut;	///< All IV in the Facility (fuel for reactor, total for all others...)


	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time

		//********* Internal Parameter *********//
private :
	int		fId;			//!< Identity of the Facility inside the Parc
		
	CLASS*		fParc;			//!< Pointer to the main Parc
	DecayDataBank*	fDecayDataBase;		//!< Pointer to the Decay DataBase

	cSecond		fCreationTime;		///< CLASS Universal Time of Creation
	cSecond		fLifeTime;		///< LifeTime Of the Reactor (Operating's Duration)

	ClassDef(CLSSFacility,1);
};

#endif

