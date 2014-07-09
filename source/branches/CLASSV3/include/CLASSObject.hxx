
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

#include "CLASSLogger.hxx"

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
	CLASSObject(CLASSLogger* log);

	
	virtual CLASSObject* Clone()	{ return new CLASSObject(*this); } //!< Correct way to copy a CLASSObject in case of derivation


	void		SetLog(CLASSLogger* log)	{ fLog = log;}		//!< Set the CLASSLogger

	CLASSLogger*	GetLog()		{ return fLog; }		//!< Return the Pointer to the Log

	using TNamed::SetName;
protected :
	CLASSLogger*	fLog;			//!< Pointer to the Log


private :

	ClassDef(CLASSObject,0);
};

#endif

