
#ifndef _CLASSOBJECT_HXX
#define _CLASSOBJECT_HXX


/*!
 \file
 \brief Header file for CLASSObject class.
 
 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>

#include "LogFile.hxx"

#include "TNamed.h"

using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a CLASS Object.
 The aim of these class is synthetyse all the commum properties to all CLASS Element.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class CLASSObject : public TNamed
{
public :
	///< Normal Constructor.
	CLASSObject();
	CLASSObject(LogFile* log);
	
	virtual CLASSObject* Clone()	{ return new CLASSObject(*this); } //!< Correct way to copy a CLASSObject in case of derivation


	void		SetLog(LogFile* log)	{ fLog = log; fIsLog = true; }		//!< Set the LogFile

	LogFile*	GetLog()		{ return fLog; }		//!< Return the Pointer to the Log
	bool		IsLog()			{ return fIsLog; }		//!< reutrn true if a LogFile is defined
	
private :
 	LogFile*	fLog;			//!< Pointer to the Log
	bool		fIsLog;			//!< Set at true if a LogFile are define
	
	ClassDef(CLASSObject,0);
};

#endif

