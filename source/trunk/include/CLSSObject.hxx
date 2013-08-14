
#ifndef _CLASSOBJECT_HXX
#define _CLASSOBJECT_HXX


/*!
 \file
 \brief Header file for CLSSObject class.
 
 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>

#include "LogFile.hxx"

#include "TNamed.h"

using namespace std;


class CLSSObject : public TNamed
{
public :
	///< Normal Constructor.
	CLSSObject();
	CLSSObject(LogFile* log);
	virtual CLSSObject* Clone()	{ return new CLSSObject(*this); } //!< Correct way to copy a CLSSObject in case of derivation

	virtual ~CLSSObject(){}			//!< destructor

	void		SetLog(LogFile* log)	{ fLog = log; }
	
	LogFile*	GetLog()		{ return fLog; }		//!< Return the Pointer to the Log
	
private :
 	LogFile*	fLog;			//!< Pointer to the Log
	
	ClassDef(CLSSObject,0);
};

#endif

