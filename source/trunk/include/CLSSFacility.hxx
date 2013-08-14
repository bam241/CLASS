
#ifndef _CLASSFACILITY_HXX
#define _CLASSFACILITY_HXX


/*!
 \file
 \brief Header file for CLSSFacility class.
 
 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>

#include "CLSSObject.hxx"
#include "IsotopicVector.hxx"

#include "TNamed.h"

using namespace std;
typedef long long int cSecond;

class CLASS;

class CLSSFacility : public CLSSObject
{
public :
		///< Normal Constructor.
	CLSSFacility();
	
		//********* Get Method *********//
	int 		GetId()			const	{ return fId; }			//!< Return the Facility Parc'Is
	virtual IsotopicVector 	GetInsideIV()	const	{ return fInsideIV; } 		//!< Return the IV contain in the Facility

	cSecond	GetInternalTime() const		{ return fInternalTime; }	//!< Return Creation Time

	cSecond		GetCycleTime()		const	{ return fCycleTime; } 		//!< Return the cycle time of the Facility
	cSecond 	GetCreationTime()	const	{ return fCreationTime; }	//!< Return the creation time of the Facility
	cSecond 	GetLifeTime()		const	{ return fLifeTime; }		//!< Return the life time of the Facility
	CLASS*		GetParc()			{ return fParc; }
		
	
		//********* Set Method *********//
	void SetId(int id)			{ fId = id; }				//!< Set The Facility Parc'Id
	void SetParc(CLASS* parc)		{ fParc = parc; }			//!< Set the Pointer to the Parc
	
	void SetInsideIV(IsotopicVector isotopicvector)	{ fInsideIV = isotopicvector; }	//!< Set the IV inside the Facility Core
	void SetCreationTime(double creationtime)	{ fCreationTime = (cSecond)creationtime;}
	void SetLifeTime(double lifetime)		{ fLifeTime = (cSecond)lifetime; }	//!< Set the life time of the facility
	virtual void SetCycleTime(double cycletime)	{ fCycleTime = (cSecond)cycletime; }	//!< Set the cycle time (Cycle of the loading Plan)
	
	
		//********* Modification Method *********//
	virtual void Evolution(cSecond t)	{ }	//!< Performe the Evolution to the Time t
	virtual void Dump()			{ }			//!< Write Modification (IV In/Out, filling the TF...)
		
		
protected :
	bool		fIsStarted;		///< True if Running, False Otherwise
	bool		fShutDown;		///< True if ShutDown
	bool		fEndOfCycle;		///< True if Reaching the End of a Facility Cycle

		
	cSecond		fInternalTime;		///< Internal Clock
	cSecond		fInCycleTime;		///< Time spend since the beginning of the last Cycle
	cSecond		fCycleTime;		///< Cycle Time

	IsotopicVector	fInsideIV;		///< All IV in the Facility (fuel for reactor, total for all others...)

		//********* Internal Parameter *********//
private :
	int		fId;			//!< Identity of the Facility inside the Parc
		
	CLASS*		fParc;			//!< Pointer to the main Parc
	
	cSecond		fCreationTime;		///< CLASS Universal Time of Creation
	cSecond		fLifeTime;		///< LifeTime Of the Reactor (Operating's Duration)

	ClassDef(CLSSFacility,1);
};

#endif

